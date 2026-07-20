"""
bioinformatics_query.py — Local PostgreSQL query engine for bioinformatics output.

Reads each per-sample bioinformatics JSON under a bioinformatics output directory,
extracts the normalised sample summary, and runs deterministic SQL queries against
the PGIRL reference database. Results are written as one JSON file per sample
under the data_query output directory, structured for downstream use by the
epidemiological query engine and the intelligence pipeline.

Usage:
    python -m intelligence_engine.data_engine.sql_querying.bioinformatics_query \
        --bioinformatics-dir output/bioinformatics \
        --output-dir output/data_query

Output:
    output/data_query/<sample_id>/db_query_results.json
"""

import argparse
import json
import logging
import os
import re
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional

# Resolve project root so config.py can be imported regardless of cwd
PROJECT_ROOT = Path(__file__).parents[3].resolve()
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

try:
    from config import DB_URL as CONFIG_DB_URL
except Exception:
    CONFIG_DB_URL = None

try:
    import psycopg2
    import psycopg2.extras
except Exception as exc:  # pragma: no cover - deployment may install psycopg2 later
    psycopg2 = None
    _psycopg2_import_error = exc

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _get_db_url(args_url: Optional[str]) -> str:
    """Return the DB URL from CLI, environment, or config, in that order."""
    return args_url or os.environ.get("PGIRL_DB_URL") or CONFIG_DB_URL or "postgresql://localhost:5432/pgirl"


def _iso_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")


def _parse_hgvs(hgvs: str) -> Optional[dict[str, Any]]:
    """Parse a simple protein HGVS string like 'GP:A82V' into components."""
    if not hgvs:
        return None
    m = re.match(r"^([A-Za-z0-9_]+):([A-Z*])(\d+)([A-Z*])$", hgvs.strip())
    if not m:
        return None
    return {
        "gene": m.group(1),
        "ref_aa": m.group(2),
        "position": int(m.group(3)),
        "alt_aa": m.group(4),
    }


def _extract_stage9(bio_output: dict[str, Any]) -> dict[str, Any]:
    """Return the normalised stage9 section, or the whole dict if missing."""
    if not isinstance(bio_output, dict):
        return {}
    if "stage9_normalised_output" in bio_output:
        return bio_output["stage9_normalised_output"] or {}
    return bio_output


def _extract_variants(stage9: dict[str, Any]) -> list[dict[str, Any]]:
    """Return amino-acid variants from stage9."""
    if not isinstance(stage9, dict):
        return []
    mutations = stage9.get("mutations") or []
    variants = []
    for v in mutations:
        if not isinstance(v, dict):
            continue
        parsed = _parse_hgvs(v.get("hgvs_p", ""))
        if parsed:
            variants.append({**parsed, **v})
        else:
            # Accept explicit fields if HGVS parsing fails
            gene = v.get("gene")
            pos = v.get("position")
            ref = v.get("ref_aa")
            alt = v.get("alt_aa")
            if gene and pos is not None and ref and alt:
                variants.append(v)
    return variants


def _extract_negative_findings(stage9: dict[str, Any]) -> list[dict[str, Any]]:
    """Return variants explicitly marked NOT_detected in stage9."""
    if not isinstance(stage9, dict):
        return []
    findings = []
    for nf in stage9.get("negative_findings", []):
        if not isinstance(nf, dict):
            continue
        if nf.get("status") == "NOT_detected":
            parsed = _parse_hgvs(nf.get("hgvs_p", ""))
            if parsed:
                findings.append({**parsed, **nf})
            else:
                gene = nf.get("gene")
                pos = nf.get("position")
                ref = nf.get("ref_aa") or ""
                alt = nf.get("alt_aa") or ""
                if gene and pos is not None:
                    findings.append({
                        "gene": gene,
                        "position": pos,
                        "ref_aa": ref,
                        "alt_aa": alt,
                        **nf,
                    })
    return findings


def _rows_to_dicts(rows: list[Any]) -> list[dict[str, Any]]:
    """Convert psycopg2 RealDictRow / tuples to plain JSON-serialisable dicts."""
    out = []
    for row in rows:
        if isinstance(row, dict):
            out.append({k: v for k, v in row.items()})
        else:
            out.append(str(row))
    return out


# ---------------------------------------------------------------------------
# DB client
# ---------------------------------------------------------------------------


class BioinformaticsQuery:
    """Run the deterministic SQL layer for one or more bioinformatics outputs."""

    def __init__(self, db_url: str, dry_run: bool = False):
        self.db_url = db_url
        self.dry_run = dry_run
        self.conn = None

    def __enter__(self):
        if psycopg2 is None:
            raise RuntimeError("psycopg2 is not installed; cannot query the local DB")
        if not self.dry_run:
            try:
                self.conn = psycopg2.connect(self.db_url)
            except psycopg2.OperationalError as exc:
                raise RuntimeError(f"Could not connect to PGIRL DB at {self.db_url}: {exc}")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.conn is not None:
            self.conn.close()
            self.conn = None

    def _execute(self, sql: str, params: tuple = ()) -> list[dict[str, Any]]:
        """Execute a read-only query and return dict rows."""
        if self.dry_run or self.conn is None:
            return []
        with self.conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
            try:
                cur.execute(sql, params)
                return _rows_to_dicts(cur.fetchall())
            except psycopg2.Error as exc:
                log.warning("SQL failed: %s | params=%s | error=%s", sql, params, exc)
                return []

    # ------------------------------------------------------------------
    # Layer queries
    # ------------------------------------------------------------------

    def _variant_summary(
        self, pathogen_id: str, species_id: str, variants: list[dict[str, Any]]
    ) -> list[dict[str, Any]]:
        if not variants:
            return []
        rows = []
        for v in variants:
            gene = v.get("gene")
            pos = v.get("position")
            alt = v.get("alt_aa")
            if not gene or pos is None or not alt:
                continue
            # Try the post-aggregation view first, then fall back to raw protein_variants
            for sql in (
                """
                SELECT * FROM v_variant_summary
                WHERE pathogen_id = %s AND species_id = %s
                  AND gene = %s AND position = %s AND alt_aa = %s
                """,
                """
                SELECT pathogen_id, species_id, gene, position, ref_aa, alt_aa,
                       variant_type, hgvs_p, genome_count, species_total_genomes,
                       round(genome_count::numeric / NULLIF(species_total_genomes, 0), 4) AS global_frequency,
                       first_seen_date, last_seen_date, first_seen_year, last_seen_year,
                       countries_seen, country_codes, lineage_ids
                FROM protein_variants
                WHERE pathogen_id = %s AND species_id = %s
                  AND gene = %s AND position = %s AND alt_aa = %s
                """,
            ):
                result = self._execute(sql, (pathogen_id, species_id, gene, pos, alt))
                if result:
                    rows.extend(result)
                    break
        return rows

    def _variant_phenotypes(
        self, pathogen_id: str, species_id: str, variants: list[dict[str, Any]]
    ) -> list[dict[str, Any]]:
        if not variants:
            return []
        rows = []
        for v in variants:
            gene = v.get("gene")
            pos = v.get("position")
            alt = v.get("alt_aa")
            if not gene or pos is None or not alt:
                continue
            for sql in (
                """
                SELECT * FROM v_variant_with_phenotype
                WHERE pathogen_id = %s AND species_id = %s
                  AND gene = %s AND position = %s AND alt_aa = %s
                """,
                """
                SELECT pv.pathogen_id, pv.species_id, pv.gene, pv.position, pv.ref_aa, pv.alt_aa,
                       pv.hgvs_p, gp.association_id, gp.phenotype_category, gp.phenotype_specific,
                       gp.evidence_strength, gp.effect_size, gp.genotype_description, gp.literature_refs
                FROM protein_variants pv
                LEFT JOIN genotype_phenotype gp
                    ON pv.pathogen_id = gp.pathogen_id
                    AND pv.species_id = gp.species_id
                    AND pv.gene = gp.protein
                    AND pv.position = gp.position
                    AND pv.alt_aa = gp.alt_aa
                WHERE pv.pathogen_id = %s AND pv.species_id = %s
                  AND pv.gene = %s AND pv.position = %s AND pv.alt_aa = %s
                """,
            ):
                result = self._execute(sql, (pathogen_id, species_id, gene, pos, alt))
                if result:
                    rows.extend(result)
                    break
        return rows

    def _negative_finding_lookup(
        self, pathogen_id: str, species_id: str, negative_findings: list[dict[str, Any]]
    ) -> list[dict[str, Any]]:
        if not negative_findings:
            return []
        rows = []
        for nf in negative_findings:
            gene = nf.get("gene")
            pos = nf.get("position")
            alt = nf.get("alt_aa")
            if not gene or pos is None or not alt:
                continue
            for sql in (
                """
                SELECT * FROM v_variant_with_phenotype
                WHERE pathogen_id = %s AND species_id = %s
                  AND gene = %s AND position = %s AND alt_aa = %s
                """,
                """
                SELECT pv.pathogen_id, pv.species_id, pv.gene, pv.position, pv.ref_aa, pv.alt_aa,
                       pv.hgvs_p, pv.genome_count, pv.species_total_genomes,
                       round(pv.genome_count::numeric / NULLIF(pv.species_total_genomes, 0), 4) AS global_frequency,
                       gp.phenotype_category, gp.phenotype_specific, gp.evidence_strength
                FROM protein_variants pv
                LEFT JOIN genotype_phenotype gp
                    ON pv.pathogen_id = gp.pathogen_id
                    AND pv.species_id = gp.species_id
                    AND pv.gene = gp.protein
                    AND pv.position = gp.position
                    AND pv.alt_aa = gp.alt_aa
                WHERE pv.pathogen_id = %s AND pv.species_id = %s
                  AND pv.gene = %s AND pv.position = %s AND pv.alt_aa = %s
                """,
            ):
                result = self._execute(sql, (pathogen_id, species_id, gene, pos, alt))
                if result:
                    rows.extend(result)
                    break
        return rows

    def _lineage_metadata(
        self, pathogen_id: str, species_id: str, lineage_id: str
    ) -> Optional[dict[str, Any]]:
        if not lineage_id:
            return None
        sql = """
            SELECT lineage_id, lineage_name, clade, countries_reported, regions_reported,
                   first_detected, last_detected, primary_host, reservoir,
                   human_to_human, animal_to_human, number_genomes_available,
                   known_aliases, endemic_regions
            FROM lineages
            WHERE pathogen_id = %s AND species_id = %s AND lineage_id = %s
        """
        rows = self._execute(sql, (pathogen_id, species_id, lineage_id))
        return rows[0] if rows else None

    def _lineage_genomes_summary(
        self, pathogen_id: str, species_id: str, lineage_id: str
    ) -> dict[str, Any]:
        """Return aggregated collection_date/country stats for a lineage."""
        out = {
            "lineage_last_seen_date": None,
            "lineage_countries_in_db": [],
            "lineage_total_genomes": 0,
        }
        if not lineage_id:
            return out
        sql = """
            SELECT
                COUNT(*) AS total_genomes,
                MAX(collection_date) AS last_seen_date,
                ARRAY_AGG(DISTINCT collection_country ORDER BY collection_country) AS countries
            FROM genome_metadata
            WHERE pathogen_id = %s AND species_id = %s AND lineage_id = %s
        """
        for row in self._execute(sql, (pathogen_id, species_id, lineage_id)):
            out["lineage_last_seen_date"] = row.get("last_seen_date")
            out["lineage_countries_in_db"] = row.get("countries") or []
            out["lineage_total_genomes"] = row.get("total_genomes") or 0
        return out

    def _lineage_defining_variants(
        self, pathogen_id: str, species_id: str, lineage_id: str, top_n: int = 20
    ) -> list[dict[str, Any]]:
        if not lineage_id:
            return []
        for sql in (
            """
            SELECT * FROM v_variant_summary
            WHERE pathogen_id = %s AND species_id = %s AND %s = ANY(lineage_ids)
            ORDER BY genome_count DESC NULLS LAST
            LIMIT %s
            """,
            """
            SELECT pathogen_id, species_id, gene, position, ref_aa, alt_aa,
                   variant_type, hgvs_p, genome_count, species_total_genomes,
                   round(genome_count::numeric / NULLIF(species_total_genomes, 0), 4) AS global_frequency,
                   first_seen_date, last_seen_date, countries_seen
            FROM protein_variants
            WHERE pathogen_id = %s AND species_id = %s AND %s = ANY(lineage_ids)
            ORDER BY genome_count DESC NULLS LAST
            LIMIT %s
            """,
        ):
            rows = self._execute(sql, (pathogen_id, species_id, lineage_id, top_n))
            if rows:
                return rows
        return []

    def _lineage_in_country(
        self, pathogen_id: str, species_id: str, lineage_id: str, country: str
    ) -> list[dict[str, Any]]:
        if not lineage_id or not country:
            return []
        sql = """
            SELECT collection_country, collection_year,
                   COUNT(*) AS genome_count,
                   MIN(collection_date) AS earliest_date,
                   MAX(collection_date) AS latest_date
            FROM genome_metadata
            WHERE pathogen_id = %s AND species_id = %s AND lineage_id = %s
              AND collection_country = %s
            GROUP BY collection_country, collection_year
            ORDER BY collection_year
        """
        return self._execute(sql, (pathogen_id, species_id, lineage_id, country))

    def _lineages_in_country(
        self, pathogen_id: str, species_id: str, country: str
    ) -> list[dict[str, Any]]:
        if not country:
            return []
        sql = """
            SELECT lineage_id, COUNT(*) AS genome_count,
                   MIN(collection_date) AS earliest_date,
                   MAX(collection_date) AS latest_date
            FROM genome_metadata
            WHERE pathogen_id = %s AND species_id = %s
              AND collection_country = %s
              AND lineage_id IS NOT NULL
            GROUP BY lineage_id
            ORDER BY genome_count DESC
        """
        return self._execute(sql, (pathogen_id, species_id, country))

    def _species_country_totals(
        self, pathogen_id: str, species_id: str, top_n: int = 20
    ) -> list[dict[str, Any]]:
        sql = """
            SELECT collection_country, collection_country_code,
                   COUNT(*) AS total_genomes,
                   MIN(collection_date) AS earliest_date,
                   MAX(collection_date) AS latest_date
            FROM genome_metadata
            WHERE pathogen_id = %s AND species_id = %s
            GROUP BY collection_country, collection_country_code
            ORDER BY total_genomes DESC
            LIMIT %s
        """
        return self._execute(sql, (pathogen_id, species_id, top_n))

    def _species_year_totals(
        self, pathogen_id: str, species_id: str
    ) -> list[dict[str, Any]]:
        sql = """
            SELECT collection_year, COUNT(*) AS total_genomes
            FROM genome_metadata
            WHERE pathogen_id = %s AND species_id = %s AND collection_year IS NOT NULL
            GROUP BY collection_year
            ORDER BY collection_year
        """
        return self._execute(sql, (pathogen_id, species_id))

    def _gene_function(
        self, pathogen_id: str, species_id: str, genes: list[str]
    ) -> list[dict[str, Any]]:
        if not genes:
            return []
        # dedupe while preserving order
        seen = set()
        unique_genes = [g for g in genes if g and not (g in seen or seen.add(g))]
        sql = """
            SELECT gene, protein_name, protein_function, protein_length_aa,
                   key_domains, functional_sites, known_hotspots, pdb_ids
            FROM gene_function
            WHERE pathogen_id = %s AND species_id = %s AND gene = ANY(%s)
        """
        return self._execute(sql, (pathogen_id, species_id, unique_genes))

    # ------------------------------------------------------------------
    # Main per-sample runner
    # ------------------------------------------------------------------

    def query_sample(self, bio_output: dict[str, Any]) -> dict[str, Any]:
        """Run all deterministic DB queries for a single bioinformatics output."""
        stage9 = _extract_stage9(bio_output)
        if not stage9:
            raise ValueError("Bioinformatics output contains no stage9_normalised_output")

        # Prefer the top-level sample block (which matches the folder) over stage9
        sample_id = bio_output.get("sample", {}).get("sample_id") or stage9.get("sample_id", "unknown")
        pathogen_id = stage9.get("pathogen", stage9.get("pathogen_id"))
        species_id = stage9.get("species_id")
        species_name = stage9.get("species", stage9.get("species_name"))
        lineage_id = stage9.get("lineage")
        # The bioinformatics stage now emits `clade` (e.g. "Ebov-1976") while
        # leaving `lineage` empty. The curated DB keys lineages by the full
        # "{species_id}-{clade}" id (e.g. "EBOV-Ebov-1976"), matched exactly in
        # SQL, so derive that id from the clade when no explicit lineage exists.
        if not str(lineage_id or "").strip():
            clade = str(stage9.get("clade") or "").strip()
            if clade and species_id:
                sid = str(species_id).strip()
                # DB ids use the uppercase species prefix (e.g. "EBOV-Ebov-1976")
                # while Nextclade emits the title-cased short clade ("Ebov-1976").
                # Use a case-sensitive check so the short clade is prefixed and an
                # already-full id is left untouched (no double-prefixing).
                lineage_id = clade if clade.startswith(f"{sid}-") else f"{sid}-{clade}"
            elif clade:
                lineage_id = clade
        metadata = stage9.get("metadata", {}) if isinstance(stage9.get("metadata"), dict) else {}
        country = metadata.get("country", stage9.get("collection_country"))
        collection_date = metadata.get("collection_date", stage9.get("collection_date"))

        if not pathogen_id or not species_id:
            raise ValueError("stage9 is missing pathogen_id and/or species_id")

        variants = _extract_variants(stage9)
        negative_findings = _extract_negative_findings(stage9)
        genes = [v.get("gene") for v in variants if v.get("gene")]

        # Layer 1: variant-level lookup
        variant_summaries = self._variant_summary(pathogen_id, species_id, variants)
        variant_phenotypes = self._variant_phenotypes(pathogen_id, species_id, variants)
        negative_lookup = self._negative_finding_lookup(pathogen_id, species_id, negative_findings)

        # Layer 2: lineage context
        lineage_meta = self._lineage_metadata(pathogen_id, species_id, lineage_id) or {}
        lineage_genomes = self._lineage_genomes_summary(pathogen_id, species_id, lineage_id)
        lineage_defining = self._lineage_defining_variants(pathogen_id, species_id, lineage_id)

        # Layer 3: geographic/temporal context
        lineage_in_country = self._lineage_in_country(pathogen_id, species_id, lineage_id, country)
        lineages_in_country = self._lineages_in_country(pathogen_id, species_id, country)

        # Layer 4: gene function
        gene_function = self._gene_function(pathogen_id, species_id, genes)

        # Layer 5: species-level surveillance
        country_totals = self._species_country_totals(pathogen_id, species_id)
        year_totals = self._species_year_totals(pathogen_id, species_id)

        # Build the summary used by the epidemiological query engine
        lineage_countries = lineage_genomes.get("lineage_countries_in_db") or []
        lineage_not_in_db_for = []
        if country and country not in lineage_countries:
            lineage_not_in_db_for = [country]

        variant_frequencies = {}
        for row in variant_summaries:
            hgvs = row.get("hgvs_p")
            if hgvs:
                variant_frequencies[hgvs] = {
                    "global_frequency": row.get("global_frequency"),
                    "genome_count": row.get("genome_count"),
                    "species_total_genomes": row.get("species_total_genomes"),
                    "first_seen_date": row.get("first_seen_date"),
                    "last_seen_date": row.get("last_seen_date"),
                    "countries_seen": row.get("countries_seen") or [],
                }

        phenotype_associations = {}
        for row in variant_phenotypes:
            hgvs = row.get("hgvs_p")
            if not hgvs:
                continue
            pheno = {
                "phenotype_category": row.get("phenotype_category"),
                "phenotype_specific": row.get("phenotype_specific"),
                "evidence_strength": row.get("evidence_strength"),
                "effect_size": row.get("effect_size"),
                "genotype_description": row.get("genotype_description"),
                "association_id": row.get("association_id"),
            }
            phenotype_associations.setdefault(hgvs, []).append(pheno)

        local_db_results = {
            "lineage_last_seen_date": lineage_genomes.get("lineage_last_seen_date"),
            "lineage_countries_in_db": lineage_countries,
            "lineage_not_in_db_for": lineage_not_in_db_for,
            "variant_frequencies": variant_frequencies,
            "phenotype_associations": phenotype_associations,
        }

        result = {
            "sample_id": sample_id,
            "pathogen_id": pathogen_id,
            "species_id": species_id,
            "species_name": species_name,
            "lineage_id": lineage_id,
            "collection_country": country,
            "collection_date": collection_date,
            "query_timestamp": _iso_now(),
            "db_url": None,
            "local_db_results": local_db_results,
            "layer1_variant_lookup": {
                "variant_summaries": variant_summaries,
                "variant_phenotypes": variant_phenotypes,
                "negative_findings": negative_lookup,
            },
            "layer2_lineage_context": {
                "lineage_metadata": lineage_meta,
                "lineage_genomes_summary": lineage_genomes,
                "lineage_defining_variants": lineage_defining,
            },
            "layer3_geographic_temporal": {
                "lineage_in_country": lineage_in_country,
                "lineages_in_country": lineages_in_country,
                "variant_temporal": {
                    "first_seen_summary": [
                        {
                            "hgvs_p": row.get("hgvs_p"),
                            "first_seen_date": row.get("first_seen_date"),
                            "last_seen_date": row.get("last_seen_date"),
                            "first_seen_year": row.get("first_seen_year"),
                            "last_seen_year": row.get("last_seen_year"),
                            "genome_count": row.get("genome_count"),
                            "countries_seen": row.get("countries_seen") or [],
                        }
                        for row in variant_summaries
                    ],
                },
            },
            "layer4_gene_function_context": {
                "gene_function": gene_function,
            },
            "layer5_species_surveillance": {
                "country_totals": country_totals,
                "year_totals": year_totals,
            },
        }

        return result

    def run_batch(
        self,
        bioinformatics_dir: Path,
        output_dir: Path,
        sample_filter: Optional[str] = None,
    ) -> list[Path]:
        """Process all samples under bioinformatics_dir and write results to output_dir."""
        written = []
        if not bioinformatics_dir.exists():
            raise FileNotFoundError(f"Bioinformatics directory not found: {bioinformatics_dir}")

        sample_dirs = [d for d in bioinformatics_dir.iterdir() if d.is_dir()]
        if not sample_dirs:
            log.warning("No sample directories found under %s", bioinformatics_dir)
            return []

        for sample_dir in sorted(sample_dirs):
            sample_id = sample_dir.name
            if sample_filter and sample_id != sample_filter:
                continue

            bio_json = sample_dir / "bio_output.json"
            if not bio_json.exists():
                log.warning("Skipping %s — no bio_output.json", sample_dir)
                continue

            log.info("Processing %s", sample_id)
            try:
                with open(bio_json) as f:
                    bio_output = json.load(f)
            except json.JSONDecodeError as exc:
                log.error("Could not parse %s: %s", bio_json, exc)
                continue

            try:
                result = self.query_sample(bio_output)
            except Exception as exc:
                log.error("Query failed for %s: %s", sample_id, exc)
                continue

            out_sample_dir = output_dir / sample_id
            out_sample_dir.mkdir(parents=True, exist_ok=True)
            out_path = out_sample_dir / "db_query_results.json"
            with open(out_path, "w") as f:
                json.dump(result, f, indent=2, default=str, ensure_ascii=False)
            log.info("Wrote %s", out_path)
            written.append(out_path)

        return written


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(
        description="Run deterministic local DB queries for bioinformatics pipeline outputs."
    )
    parser.add_argument(
        "--bioinformatics-dir",
        default="output/bioinformatics",
        help="Directory containing per-sample bioinformatics output folders (default: output/bioinformatics)",
    )
    parser.add_argument(
        "--output-dir",
        default="output/data_query",
        help="Directory to write per-sample db_query_results.json files (default: output/data_query)",
    )
    parser.add_argument(
        "--sample",
        default=None,
        help="Process only a single sample directory (optional)",
    )
    parser.add_argument(
        "--db-url",
        default=None,
        help="PostgreSQL connection URL (default: PGIRL_DB_URL env var, then config.DB_URL)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Parse inputs and build queries without connecting to the database",
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Enable debug logging",
    )
    args = parser.parse_args()

    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)

    if psycopg2 is None and not args.dry_run:
        log.error("psycopg2 is not installed. Use the anaconda python for DB scripts: /Users/christianndekezi/anaconda3/bin/python3")
        sys.exit(1)

    bioinformatics_dir = Path(args.bioinformatics_dir)
    output_dir = Path(args.output_dir)
    db_url = _get_db_url(args.db_url)

    try:
        with BioinformaticsQuery(db_url=db_url, dry_run=args.dry_run) as engine:
            written = engine.run_batch(
                bioinformatics_dir=bioinformatics_dir,
                output_dir=output_dir,
                sample_filter=args.sample,
            )
    except Exception as exc:
        log.error("Batch failed: %s", exc)
        sys.exit(1)

    if written:
        log.info("Complete. Wrote %d result file(s) to %s", len(written), output_dir)
    else:
        log.info("Complete. No result files were written.")


if __name__ == "__main__":
    main()
