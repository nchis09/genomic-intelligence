"""
db_sync.py — API-to-Database bridge for the Genomic Epidemic Intelligence System.

Orchestrates data syncing into the PostgreSQL database. By default it runs the
Ebola protein variant calling pipeline (incremental, skip-existing). Optional
sources such as WHO DON outbreak alerts can be enabled once API access is
configured.

Design principles:
  - CURATED data is NEVER overwritten by API data (data_source = 'curated' is protected)
  - API records are upserted: insert new, update existing API-sourced rows with latest values
  - Unknown pathogens are logged and skipped — never silently inserted
  - Extra columns from the API that don't exist in the DB are ignored
  - Missing columns from the API are filled with NULL
  - Switching API providers is transparent — all providers implement the same interface
  - Epidemiological context (outbreaks, country indicators) is fetched on demand
    by the intelligence engine, not stored in the reference database

Usage:
    python scripts/db_sync.py
    python scripts/db_sync.py --days-back 180
    python scripts/db_sync.py --dry-run

Environment variable alternative (set DB URL once, never type it again):
    export PGIRL_DB_URL=postgresql://localhost:5432/pgirl
    python scripts/db_sync.py

DB URL is loaded from config.py automatically — no hardcoded paths.
"""

import argparse
import logging
import os
import subprocess
import sys
from pathlib import Path
from typing import Optional

try:
    import psycopg2
    import psycopg2.extras
except ImportError:
    print("ERROR: psycopg2 not installed. Run: pip install psycopg2-binary")
    sys.exit(1)

# Resolve project root and load config — works regardless of where project lives
PROJECT_ROOT = Path(__file__).parent.parent.resolve()
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "database" / "api_sources"))

try:
    from config import DB_URL as CONFIG_DB_URL, API_SOURCES_DIR
except ImportError:
    CONFIG_DB_URL = "postgresql://localhost:5432/pgirl"
    API_SOURCES_DIR = PROJECT_ROOT / "database" / "api_sources"

# Optional additional sources — each is independent; if a module is missing the
# corresponding sync is simply skipped (never crashes the whole run).
try:
    from ncbi_api import NCBISource
except ImportError:
    NCBISource = None
try:
    from nextstrain_api import NextstrainSource
except ImportError:
    NextstrainSource = None

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)


# =============================================================================
# PATHOGEN ROUTING TABLE
# Maps API disease name strings → (pathogen_id, species_id) in the DB.
# Add new pathogens here as you expand beyond Ebola.
# =============================================================================

PATHOGEN_ROUTING = {
    # Ebola species
    "zaire ebolavirus":      ("ebola", "EBOV"),
    "ebola virus disease":   ("ebola", "EBOV"),
    "ebola":                 ("ebola", "EBOV"),
    "ebov":                  ("ebola", "EBOV"),
    "sudan ebolavirus":      ("ebola", "SUDV"),
    "sudan virus":           ("ebola", "SUDV"),
    "sudv":                  ("ebola", "SUDV"),
    "bundibugyo ebolavirus": ("ebola", "BDBV"),
    "bdbv":                  ("ebola", "BDBV"),
    "reston ebolavirus":     ("ebola", "RESTV"),
    "restv":                 ("ebola", "RESTV"),
    "tai forest ebolavirus": ("ebola", "TAFV"),
    "tafv":                  ("ebola", "TAFV"),
    "bombali ebolavirus":    ("ebola", "BOMV"),
    "bomv":                  ("ebola", "BOMV"),
    # Future pathogens — uncomment and add species_id when you add them to the DB
    # "dengue":              ("dengue", "DENV"),
    # "dengue virus":        ("dengue", "DENV"),
    # "influenza a":         ("influenza", "IAV"),
    # "influenza b":         ("influenza", "IBV"),
    # "mpox":                ("mpox", "MPXV"),
    # "monkeypox":           ("mpox", "MPXV"),
    # "rift valley fever":   ("rvf", "RVFV"),
}


# =============================================================================
# PATHOGEN ROUTING
# =============================================================================

def route_pathogen(species_str: str) -> Optional[tuple[str, str]]:
    """
    Map an API species string to (pathogen_id, species_id).
    Returns None if the species is unrecognised — caller should skip the record.
    """
    if not species_str:
        return None
    key = species_str.strip().lower()
    result = PATHOGEN_ROUTING.get(key)
    if result is None:
        log.warning(
            "Unknown species '%s' — no pathogen routing found. "
            "Add it to PATHOGEN_ROUTING in db_sync.py to enable this pathogen.",
            species_str,
        )
    return result


# =============================================================================
# GENOMES  (NCBI -> reference_genomes)
# =============================================================================

def sync_genomes(conn, pathogen: str = "ebola", dry_run: bool = False, retmax: int = 30):
    """Fetch RefSeq + GenBank genomes from NCBI and upsert into reference_genomes."""
    if NCBISource is None:
        log.warning("ncbi_api not available — skipping genome sync")
        return
    src = NCBISource()
    cur = conn.cursor()

    # Which species belong to this pathogen (drive queries off the DB backbone)
    cur.execute("SELECT species_id, species_name FROM species WHERE pathogen_id = %s", (pathogen,))
    species_rows = cur.fetchall()
    if not species_rows:
        log.warning("No species found for pathogen '%s' — skipping genome sync", pathogen)
        cur.close()
        return

    inserted = updated = skipped = 0
    for species_id, species_name in species_rows:
        genomes = src.fetch_refseq_genomes(species_name, retmax=5)
        genomes += src.fetch_genbank_genomes(species_name, retmax=retmax)
        for g in genomes:
            row = {
                "accession":          g.accession,
                "pathogen_id":        pathogen,
                "species_id":         species_id,
                "genome_role":        "canonical_reference" if g.is_refseq else "sequence",
                "genome_length":      g.length,
                "collection_year":    g.collection_year,
                "collection_country": g.collection_country or None,
                "source_database":    g.source_database,
            }
            r = _upsert_generic(cur, "reference_genomes", "accession", row, dry_run)
            inserted += r == "inserted"; skipped += r.startswith("skipped")

    if not dry_run:
        conn.commit()
    cur.close()
    log.info("Genomes: inserted=%d skipped=%d", inserted, skipped)


# =============================================================================
# LINEAGES  (Nextstrain -> lineages)
# =============================================================================

def sync_lineages(conn, pathogen: str = "ebola", dry_run: bool = False):
    """Fetch clade/lineage aggregates from Nextstrain and upsert into lineages."""
    if NextstrainSource is None:
        log.warning("nextstrain_api not available — skipping lineage sync")
        return
    src = NextstrainSource()
    cur = conn.cursor()

    lineages = src.fetch_lineages(pathogen)
    inserted = updated = skipped = 0
    for l in lineages:
        routing = route_pathogen(l.species)
        if routing is None:
            continue
        pathogen_id, species_id = routing
        row = {
            "lineage_id":               l.lineage_id,
            "pathogen_id":              pathogen_id,
            "species_id":               species_id,
            "lineage_name":             l.lineage_name,
            "clade":                    l.clade or None,
            "first_country_detected":   l.first_country_detected or None,
            "countries_reported":       l.countries_reported or None,
            "first_detected":           l.first_detected,
            "last_detected":            l.last_detected,
            "number_genomes_available": l.number_genomes_available,
            "primary_host":             l.primary_host or None,
            "last_updated":             l.last_updated,
        }
        r = _upsert_generic(cur, "lineages", "lineage_id", row, dry_run)
        inserted += r == "inserted"; skipped += r.startswith("skipped")

    if not dry_run:
        conn.commit()
    cur.close()
    log.info("Lineages: inserted=%d skipped=%d", inserted, skipped)


# =============================================================================
# GENERIC UPSERT HELPERS
# =============================================================================

def _upsert_generic(cur, table: str, pk: str, row: dict, dry_run: bool) -> str:
    """
    Insert a row keyed on a single PK column, or skip if it already exists.
    Existing rows are never overwritten — reference data is curated and immutable.
    """
    pk_val = row.get(pk)
    if not pk_val:
        return "skipped_no_id"

    cur.execute(f"SELECT 1 FROM {table} WHERE {pk} = %s", (pk_val,))
    if cur.fetchone():
        return "skipped_existing"

    if not dry_run:
        cols = list(row.keys())
        cur.execute(
            f"INSERT INTO {table} ({', '.join(cols)}) VALUES ({', '.join(['%s'] * len(cols))})",
            [row[c] for c in cols])
    return "inserted"


# =============================================================================
# CLI
# =================================================
# VARIANTS  (NCBI -> protein_variants / genome_metadata)
# =================================================

def _run_subprocess(cmd: list, label: str, dry_run: bool = False):
    """Run a Python helper and log its output."""
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=False)
        if result.stdout:
            log.info(result.stdout)
        if result.returncode != 0:
            log.error("%s failed:\n%s", label, result.stderr)
            return False
        log.info("%s completed", label)
        return True
    except Exception as e:
        log.error("Could not run %s: %s", label, e)
        return False


def _flag_pubtator_only(cur, dry_run: bool = False):
    """Mark rows that are not corroborated by UniProt genotype evidence.

    Note: data_source column was removed during schema cleanup. This function
    now uses record_flagged + flag_reason to identify PubTator-origin rows that
    still need corroboration. Rows previously seeded from UniProt will have
    record_flagged = false; PubTator-origin rows start with record_flagged = true.
    """
    if dry_run:
        return
    # Corroborated rows: clear manual-review flag (same protein/position/alt_aa exists with flag cleared)
    cur.execute(
        """
        UPDATE genotype_phenotype gp
        SET record_flagged = false,
            flag_reason = 'Corroborated by high-confidence seed'
        WHERE gp.record_flagged = true
          AND EXISTS (
              SELECT 1 FROM genotype_phenotype up
              WHERE up.record_flagged = false
                AND up.pathogen_id = gp.pathogen_id
                AND up.species_id = gp.species_id
                AND up.protein = gp.protein
                AND up.position = gp.position
                AND up.alt_aa = gp.alt_aa
          )
        """
    )
    # Uncorroborated rows: keep flagged for manual review
    cur.execute(
        """
        UPDATE genotype_phenotype gp
        SET flag_reason = 'No corroboration found (manual review needed)'
        WHERE gp.record_flagged = true
          AND NOT EXISTS (
              SELECT 1 FROM genotype_phenotype up
              WHERE up.record_flagged = false
                AND up.pathogen_id = gp.pathogen_id
                AND up.species_id = gp.species_id
                AND up.protein = gp.protein
                AND up.position = gp.position
                AND up.alt_aa = gp.alt_aa
          )
        """
    )


def sync_genotype_phenotype(conn, pathogen: str = "ebola", dry_run: bool = False):
    """Run UniProt seeding + PubTator-based extraction + verification.

    For Ebola, this first inserts curated UniProt mutagenesis/natural-variant
    associations (high-confidence seeds), then runs PubTator as a candidate
    generator. PubTator rows are verified and any that are not corroborated by
    a UniProt row on protein/position/alt_aa are flagged for manual review.
    """
    if pathogen != "ebola":
        log.warning("Genotype-phenotype extraction is currently only implemented for Ebola")
        return

    python = Path("/Users/christianndekezi/anaconda3/bin/python3")
    gp_dir = PROJECT_ROOT / "database" / "ebola" / "genotype_phenotype"
    seed = gp_dir / "seed_from_uniprot.py"
    extractor = gp_dir / "extract_from_pubtator.py"
    verifier = gp_dir / "verify_extracted_pmids.py"

    # 1. Seed curated UniProt associations (idempotent)
    log.info("=== SOURCE: genotype_phenotype (UniProt high-confidence seed) ===")
    seed_cmd = [str(python), str(seed)]
    if dry_run:
        seed_cmd.append("--dry-run")
    if not _run_subprocess(seed_cmd, "UniProt seed", dry_run):
        return

    # 2. PubTator candidate generation
    log.info("=== SOURCE: genotype_phenotype (PubTator literature extraction) ===")
    extract_cmd = [str(python), str(extractor), "--max-results", "50"]
    if dry_run:
        extract_cmd.append("--dry-run")
    if not _run_subprocess(extract_cmd, "PubTator extraction", dry_run):
        return

    if not dry_run:
        # 3. Validate literature PMIDs
        if not _run_subprocess([str(python), str(verifier), "--batch-size", "100"], "PMID verification"):
            log.warning("Continuing despite PMID verification failure")

        # 4. Flag PubTator-only rows
        log.info("Flagging uncorroborated PubTator rows...")
        cur = conn.cursor()
        _flag_pubtator_only(cur)
        conn.commit()
        cur.close()

        # 5. Refresh dependent materialized views
        log.info("Refreshing phenotype surveillance materialized views...")
        cur = conn.cursor()
        for pheno_view in (
            "mv_mutation_with_phenotype",
            "mv_phenotype_surveillance",
            "mv_phenotype_geo_temporal",
        ):
            try:
                cur.execute(f"REFRESH MATERIALIZED VIEW CONCURRENTLY {pheno_view};")
                conn.commit()
                log.info("%s refreshed", pheno_view)
            except Exception as e:
                conn.rollback()
                try:
                    cur.execute(f"REFRESH MATERIALIZED VIEW {pheno_view};")
                    conn.commit()
                    log.info("%s refreshed", pheno_view)
                except Exception as e2:
                    conn.rollback()
                    log.error("Could not refresh %s: %s", pheno_view, e2)
        cur.close()


def sync_gene_function(conn, pathogen: str = "ebola", dry_run: bool = False):
    """Seed gene/protein functional context from UniProt.

    For each ebolavirus species, species-specific UniProt entries are fetched
    so functional annotation (protein name, function, domains, motifs, PDB IDs,
    literature) matches the correct organism. Coordinates and protein lengths
    are taken from the reference_proteomes table. Domains/sites beyond the
    reference protein length are filtered out.
    """
    if pathogen != "ebola":
        log.warning("Gene-function seeding is currently only implemented for Ebola")
        return

    python = Path("/Users/christianndekezi/anaconda3/bin/python3")
    script = PROJECT_ROOT / "database" / "ebola" / "gene_function" / "seed_gene_function_from_uniprot.py"
    cmd = [str(python), str(script)]
    if dry_run:
        cmd.append("--dry-run")

    log.info("=== SOURCE: gene_function (UniProt functional annotation) ===")
    log.info("Running: %s", " ".join(cmd))
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=False)
        if result.stdout:
            log.info(result.stdout)
        if result.returncode != 0:
            log.error("Gene function seeding failed: %s", result.stderr)
        else:
            log.info("Gene function seeding completed")
    except Exception as e:
        log.error("Could not run gene function seeding: %s", e)


def sync_variants(conn, pathogen: str = "ebola", dry_run: bool = False):
    """Run the protein variant calling pipeline for the pathogen.

    For Ebola, this invokes database/ebola/protein_variants/call_variants.py
    which fetches genomes from NCBI, extracts proteins, calls amino acid
    variants, and stores only metadata + variant calls in the DB.
    """
    if pathogen != "ebola":
        log.warning("Variant calling is currently only implemented for Ebola")
        return

    script = PROJECT_ROOT / "database" / "ebola" / "protein_variants" / "call_variants.py"
    python = Path("/Users/christianndekezi/anaconda3/bin/python3")

    cmd = [str(python), str(script), "--batch-size", "50"]
    if dry_run:
        cmd.append("--dry-run")
    else:
        # Check if protein_variants is empty — if so, use --rebuild to
        # reprocess all genomes (fetch from NCBI, call variants, aggregate)
        cur = conn.cursor()
        cur.execute("SELECT count(*) FROM protein_variants WHERE pathogen_id = %s", (pathogen,))
        variant_count = cur.fetchone()[0]
        cur.close()
        if variant_count == 0:
            log.info("protein_variants is empty — using REBUILD mode to repopulate aggregated variants")
            cmd.append("--rebuild")
        else:
            cmd.append("--skip-existing")

    log.info("=== SOURCE: variants (NCBI variant calling) ===")
    log.info("Running: %s", " ".join(cmd))
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=False)
        if result.stdout:
            log.info(result.stdout)
        if result.returncode != 0:
            log.error("Variant calling pipeline failed: %s", result.stderr)
            return
        log.info("Variant calling pipeline completed")
    except Exception as e:
        log.error("Could not run variant calling pipeline: %s", e)
        return

    # Run Nextclade lineage assignment on any new genomes missing lineage_id.
    # This is a separate subprocess so Nextclade errors don't break variant sync.
    if not dry_run:
        log.info("=== SOURCE: lineages (Nextclade clade assignment) ===")
        lineage_script = PROJECT_ROOT / "database" / "ebola" / "protein_variants" / "assign_lineages.py"
        lineage_cmd = [str(python), str(lineage_script), "--batch-size", "100"]
        log.info("Running: %s", " ".join(lineage_cmd))
        try:
            lineage_result = subprocess.run(lineage_cmd, capture_output=True, text=True, check=False)
            if lineage_result.stdout:
                log.info(lineage_result.stdout)
            if lineage_result.returncode != 0:
                log.error("Lineage assignment failed or Nextclade not installed:\n%s", lineage_result.stderr)
            else:
                log.info("Lineage assignment completed")
        except Exception as e:
            log.error("Could not run lineage assignment: %s", e)

    # No materialized views to refresh — protein_variants is now pre-aggregated
    # and lightweight views (v_variant_summary, v_variant_with_phenotype) query
    # the base tables directly at query time.
    log.info("Variant sync complete — data is pre-aggregated, no view refresh needed.")


# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Sync biological reference data into the PostgreSQL database."
    )
    parser.add_argument(
        "--db-url",
        default=os.environ.get("PGIRL_DB_URL", CONFIG_DB_URL),
        help="PostgreSQL connection URL (default loaded from config.py)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Fetch and validate records but do NOT write to the database",
    )
    parser.add_argument(
        "--pathogens",
        nargs="+",
        default=None,
        metavar="PATHOGEN_ID",
        help="Pathogen(s) to sync (default: all supported). Example: --pathogens ebola",
    )
    parser.add_argument(
        "--sources",
        nargs="+",
        default=["variants", "genotype_phenotype", "gene_function"],
        choices=["variants", "genomes", "lineages", "genotype_phenotype", "gene_function", "all"],
        help="Which sources to sync (default: variants genotype_phenotype gene_function). "
             "'all' runs every source. Example: --sources variants genotype_phenotype gene_function.",
    )
    args = parser.parse_args()

    log.info("Connecting to: %s", args.db_url[:40] + "..." if len(args.db_url) > 40 else args.db_url)

    try:
        conn = psycopg2.connect(args.db_url)
        conn.autocommit = False
    except psycopg2.OperationalError as e:
        log.error("Cannot connect to database: %s", e)
        sys.exit(1)

    sources = set(args.sources)
    run_all = "all" in sources
    pathogen = (args.pathogens or ["ebola"])[0]

    try:
        if run_all or "genomes" in sources:
            log.info("=== SOURCE: genomes (NCBI) ===")
            sync_genomes(conn, pathogen=pathogen, dry_run=args.dry_run)
        if run_all or "lineages" in sources:
            log.info("=== SOURCE: lineages (Nextstrain) ===")
            sync_lineages(conn, pathogen=pathogen, dry_run=args.dry_run)
        if run_all or "genotype_phenotype" in sources:
            sync_genotype_phenotype(conn, pathogen=pathogen, dry_run=args.dry_run)
        if run_all or "gene_function" in sources:
            sync_gene_function(conn, pathogen=pathogen, dry_run=args.dry_run)
        if run_all or "variants" in sources:
            sync_variants(conn, pathogen=pathogen, dry_run=args.dry_run)

        # No materialized views to refresh — all data is pre-aggregated
        # and lightweight views query base tables directly.
        if not args.dry_run:
            log.info("Sync complete — data is pre-aggregated, no view refresh needed.")
    except KeyboardInterrupt:
        log.warning("Interrupted by user")
        conn.rollback()
    finally:
        conn.close()


if __name__ == "__main__":
    main()
