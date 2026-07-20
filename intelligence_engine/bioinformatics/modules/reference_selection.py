#!/usr/bin/env python3
"""
Stage 2: Reference selection & context gathering.

Selects the appropriate reference genome and proteome for a newly classified
sample from the PGIRL reference database, and gathers contextual genomes for
downstream phylogenetic analysis.

All heavy lifting is delegated to validated external tools / public databases:
  - Reference metadata comes from the PGIRL PostgreSQL database.
  - Reference nucleotide sequences are fetched from NCBI Entrez on demand.
  - Reference protein sequences are fetched from the PGIRL reference_proteomes
    table (populated by database/ebola/protein_variants/fetch_reference_proteomes.py).

Output is a plain JSON / FASTA bundle consumed by the variant-calling and
phylogenetics modules.
"""

from __future__ import annotations

import json
import logging
import os
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Optional

from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

PROJECT_ROOT = Path(__file__).resolve().parents[4]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from config import DB_URL as CONFIG_DB_URL  # noqa: E402

try:
    import psycopg2
    import psycopg2.extras
except Exception:  # pragma: no cover
    psycopg2 = None

Entrez.email = os.environ.get("PGIRL_ENTREZ_EMAIL", "pgirl_pipeline@local")
Entrez.tool = "pgirl_bioinformatics_pipeline"

log = logging.getLogger(__name__)


@dataclass
class ReferenceGenome:
    accession: str
    species_id: str
    pathogen_id: str
    genome_role: Optional[str] = None
    genome_length: Optional[int] = None
    collection_year: Optional[int] = None
    collection_country: Optional[str] = None
    source_database: str = "NCBI"
    gene_coordinates: Optional[dict[str, Any]] = None
    sequence: Optional[str] = None

    def to_dict(self) -> dict[str, Any]:
        d = asdict(self)
        # sequence is kept only if explicitly requested / cached
        return d


@dataclass
class ReferenceProtein:
    gene: str
    protein_name: Optional[str]
    sequence: str
    genome_start: int
    genome_end: int
    strand: int = 1

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


class ReferenceSelector:
    """DB-backed reference genome and context-genome selector."""

    def __init__(self, db_url: Optional[str] = None, dry_run: bool = False):
        self.db_url = db_url or os.environ.get("PGIRL_DB_URL") or CONFIG_DB_URL
        self.dry_run = dry_run
        self.conn = None
        if not dry_run and psycopg2 is None:
            raise RuntimeError(
                "psycopg2 is not installed. Use anaconda python for DB scripts: "
                "/Users/christianndekezi/anaconda3/bin/python3"
            )

    def __enter__(self):
        if not self.dry_run:
            try:
                self.conn = psycopg2.connect(self.db_url)
            except psycopg2.OperationalError as exc:
                raise RuntimeError(f"Could not connect to PGIRL DB at {self.db_url}: {exc}")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.conn is not None:
            self.conn.close()

    def _query(self, sql: str, params: tuple = ()) -> list[dict[str, Any]]:
        if self.dry_run or self.conn is None:
            return []
        with self.conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
            cur.execute(sql, params)
            rows = cur.fetchall()
            return [dict(r) for r in rows]

    def select_reference(self, species_id: str) -> ReferenceGenome:
        """Return the curated reference genome for a species.

        Prefers rows marked as 'reference', then falls back to any reference
        genome row for the species.
        """
        sql = """
            SELECT accession, pathogen_id, species_id, genome_role,
                   genome_length, collection_year, collection_country,
                   source_database, gene_coordinates
            FROM reference_genomes
            WHERE species_id = %s
            ORDER BY CASE genome_role WHEN 'reference' THEN 0 ELSE 1 END,
                     accession
            LIMIT 1;
        """
        rows = self._query(sql, (species_id,))
        if not rows:
            raise ValueError(f"No reference genome found in PGIRL for species_id={species_id}")

        row = rows[0]
        return ReferenceGenome(
            accession=row["accession"],
            species_id=row["species_id"],
            pathogen_id=row["pathogen_id"],
            genome_role=row.get("genome_role"),
            genome_length=row.get("genome_length"),
            collection_year=row.get("collection_year"),
            collection_country=row.get("collection_country"),
            source_database=row.get("source_database") or "NCBI",
            gene_coordinates=row.get("gene_coordinates"),
        )

    def fetch_reference_fasta(self, accession: str, cache_dir: Optional[Path] = None) -> SeqRecord:
        """Fetch a reference nucleotide sequence from NCBI, with local caching."""
        if cache_dir:
            cache_path = cache_dir / f"{accession}.fasta"
            if cache_path.exists():
                log.info(f"Using cached reference sequence: {cache_path}")
                return SeqIO.read(cache_path, "fasta")

        log.info(f"Fetching reference sequence {accession} from NCBI")
        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
        record = SeqIO.read(handle, "fasta")
        handle.close()
        record.id = accession
        record.description = ""

        if cache_dir:
            cache_dir.mkdir(parents=True, exist_ok=True)
            SeqIO.write(record, cache_path, "fasta")

        return record

    def fetch_reference_proteome(
        self, reference_accession: str
    ) -> dict[str, ReferenceProtein]:
        """Return the curated reference proteome for a reference accession."""
        sql = """
            SELECT gene, protein_name, protein_sequence, genome_start, genome_end, strand
            FROM reference_proteomes
            WHERE reference_accession = %s
            ORDER BY genome_start;
        """
        rows = self._query(sql, (reference_accession,))
        if not rows:
            raise ValueError(
                f"No reference proteome found in PGIRL for {reference_accession}. "
                "Run database/ebola/protein_variants/fetch_reference_proteomes.py first."
            )

        return {
            row["gene"]: ReferenceProtein(
                gene=row["gene"],
                protein_name=row.get("protein_name"),
                sequence=row["protein_sequence"],
                genome_start=row["genome_start"],
                genome_end=row["genome_end"],
                strand=row.get("strand", 1),
            )
            for row in rows
        }

    def gather_context_genomes(
        self,
        species_id: str,
        exclude_sample_id: Optional[str] = None,
        limit: int = 50,
        min_length: int = 15000,
    ) -> list[dict[str, Any]]:
        """Return curated context genomes for phylogenetic analysis.

        Currently selects the most recently collected, near-complete genomes of
        the same species from the PGIRL genome_metadata table. In future this
        can be expanded to prefer outbreak-linked or geographically close
        genomes.
        """
        sql = """
            SELECT genome_accession, strain, isolate, collection_date,
                   collection_country, collection_region, host,
                   genome_length, lineage_id, source_db
            FROM genome_metadata
            WHERE species_id = %s
              AND genome_length >= %s
            ORDER BY collection_date DESC NULLS LAST
            LIMIT %s;
        """
        rows = self._query(sql, (species_id, min_length, limit))
        if exclude_sample_id:
            rows = [r for r in rows if r.get("genome_accession") != exclude_sample_id]
        return rows


def run_reference_selection(
    species_id: str,
    output_dir: Path,
    db_url: Optional[str] = None,
    cache_dir: Optional[Path] = None,
) -> dict[str, Any]:
    """Run stage 2 for a single sample and write outputs to disk.

    Returns a dict describing the selected reference genome, paths to saved
    FASTA/proteome files, and context genomes.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    cache_dir = Path(cache_dir) if cache_dir else output_dir / "reference_cache"

    with ReferenceSelector(db_url=db_url) as selector:
        ref = selector.select_reference(species_id)
        record = selector.fetch_reference_fasta(ref.accession, cache_dir=cache_dir)
        ref.sequence = str(record.seq)

        proteome = selector.fetch_reference_proteome(ref.accession)
        context_genomes = selector.gather_context_genomes(species_id)

    ref_fasta = output_dir / "reference.fasta"
    SeqIO.write(record, ref_fasta, "fasta")

    proteome_json = output_dir / "reference_proteome.json"
    proteome_json.write_text(
        json.dumps(
            {gene: prot.to_dict() for gene, prot in proteome.items()},
            indent=2,
        )
    )

    context_json = output_dir / "context_genomes.json"
    context_json.write_text(json.dumps(context_genomes, indent=2, default=str))

    summary = {
        "reference_genome": {
            "accession": ref.accession,
            "species_id": ref.species_id,
            "pathogen_id": ref.pathogen_id,
            "genome_role": ref.genome_role,
            "genome_length": ref.genome_length,
            "collection_year": ref.collection_year,
            "collection_country": ref.collection_country,
            "source_database": ref.source_database,
            "gene_coordinates": ref.gene_coordinates,
        },
        "reference_fasta": str(ref_fasta),
        "reference_proteome": {gene: prot.to_dict() for gene, prot in proteome.items()},
        "reference_proteome_json": str(proteome_json),
        "context_genomes": context_genomes,
        "context_genomes_json": str(context_json),
    }

    summary_json = output_dir / "reference_selection_summary.json"
    summary_json.write_text(json.dumps(summary, indent=2, default=str))

    log.info(f"Stage 2 complete: reference={ref.accession}, n_context={len(context_genomes)}")
    return summary


def main():
    parser = argparse.ArgumentParser(description="Stage 2: Reference selection & context gathering")
    parser.add_argument("--species-id", required=True, help="Species identifier (e.g. EBOV, SUDV)")
    parser.add_argument("--output-dir", required=True, help="Directory to write stage 2 outputs")
    parser.add_argument("--db-url", default=None, help="PostgreSQL connection URL")
    parser.add_argument("--cache-dir", default=None, help="Directory to cache NCBI reference sequences")
    parser.add_argument("--dry-run", action="store_true", help="Parse arguments without querying the DB")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    summary = run_reference_selection(
        species_id=args.species_id,
        output_dir=Path(args.output_dir),
        db_url=args.db_url,
        cache_dir=Path(args.cache_dir) if args.cache_dir else None,
    )

    print(f"Reference genome: {summary['reference_genome']['accession']}")
    print(f"Reference FASTA:  {summary['reference_fasta']}")
    print(f"Context genomes:  {len(summary['context_genomes'])}")


if __name__ == "__main__":
    import argparse
    main()
