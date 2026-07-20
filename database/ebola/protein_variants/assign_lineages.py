#!/usr/bin/env python3
"""
Assign Nextclade clades/lineages to Ebola genomes stored in the database.

This script:
  1. Finds genomes in genome_metadata missing lineage_id (or all if --reassign)
  2. Fetches their FASTA sequences from NCBI
  3. Runs Nextclade per species
  4. Inserts newly observed clades into the lineages table
  5. Updates genome_metadata.lineage_id with the assigned clade

Requires Nextclade CLI:
    /Users/christianndekezi/anaconda3/bin/conda install -c bioconda nextclade -y

Run with anaconda python:
    /Users/christianndekezi/anaconda3/bin/python3 database/ebola/protein_variants/assign_lineages.py

Options:
    --species EBOV          Only assign for one species
    --batch-size 100        Number of genomes to process per Nextclade run
    --reassign              Re-assign lineages even if already present
    --dry-run               Don't write to DB, just report assignments
"""

import argparse
import logging
import sys
from pathlib import Path

import psycopg2

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))
from config import DB_URL

from intelligence_engine.bioinformatics.nextclade_runner import (
    assign_lineages_from_accessions,
    check_nextclade,
    SPECIES_DATASET,
)

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s", datefmt="H:%M:%S")
log = logging.getLogger(__name__)


EBOLA_SPECIES = ["EBOV", "SUDV", "BDBV", "RESTV", "TAFV", "BOMV"]


def get_genomes_to_assign(cur, species_id: str, reassign: bool = False):
    """Fetch accessions for genomes needing lineage assignment."""
    sql = """
        SELECT genome_accession
        FROM genome_metadata
        WHERE pathogen_id = 'ebola' AND species_id = %s
    """
    if not reassign:
        sql += " AND lineage_id IS NULL"
    cur.execute(sql, (species_id,))
    return [row[0] for row in cur.fetchall()]


def ensure_lineage(cur, lineage_id: str, pathogen_id: str, species_id: str):
    """Insert a clade into lineages table if it does not exist."""
    cur.execute("SELECT 1 FROM lineages WHERE lineage_id = %s", (lineage_id,))
    if cur.fetchone():
        return

    # Use clade as lineage_name
    cur.execute(
        """
        INSERT INTO lineages (
            lineage_id, pathogen_id, species_id, lineage_name, clade,
            data_source, verification_status
        ) VALUES (%s, %s, %s, %s, %s, %s, %s)
        ON CONFLICT (lineage_id) DO NOTHING
        """,
        (lineage_id, pathogen_id, species_id, lineage_id, lineage_id, "nextclade", "unverified"),
    )
    log.info(f"Registered new lineage: {lineage_id}")


def update_genome_lineage(cur, accession: str, lineage_id: str):
    """Update genome_metadata.lineage_id for one accession."""
    cur.execute(
        "UPDATE genome_metadata SET lineage_id = %s WHERE genome_accession = %s",
        (lineage_id, accession),
    )


def run_assignment(conn, species_filter: str = None, batch_size: int = 100,
                   reassign: bool = False, dry_run: bool = False):
    """Run lineage assignment for Ebola genomes."""
    if not check_nextclade():
        log.error(
            "Nextclade CLI not found. Install with:\n"
            "  /Users/christianndekezi/anaconda3/bin/conda install -c bioconda nextclade -y\n"
            "Then run: nextclade dataset get --name ebola --output-dir ~/.nextclade/datasets/ebola"
        )
        return

    species_list = [species_filter] if species_filter else EBOLA_SPECIES
    cur = conn.cursor()
    total_assigned = 0
    total_skipped_no_dataset = 0

    for species_id in species_list:
        if species_id not in SPECIES_DATASET:
            log.warning(f"No Nextclade dataset for {species_id}; skipping")
            total_skipped_no_dataset += 1
            continue

        accessions = get_genomes_to_assign(cur, species_id, reassign=reassign)
        log.info(f"{species_id}: {len(accessions)} genomes to assign")
        if not accessions:
            continue

        assigned = 0
        for i in range(0, len(accessions), batch_size):
            batch = accessions[i:i+batch_size]
            try:
                assignments = assign_lineages_from_accessions(batch, species_id)
            except Exception as e:
                log.error(f"Nextclade failed for {species_id} batch {i//batch_size + 1}: {e}")
                continue

            for acc, info in assignments.items():
                clade = info.get("clade")
                if not clade:
                    log.warning(f"No clade returned for {acc}")
                    continue

                lineage_id = f"{species_id}-{clade}"
                if not dry_run:
                    ensure_lineage(cur, lineage_id, "ebola", species_id)
                    update_genome_lineage(cur, acc, lineage_id)
                log.info(f"  {acc} -> {lineage_id} (qc={info.get('qc_status')})")
                assigned += 1

            if not dry_run:
                conn.commit()
                log.info(f"  Committed batch {i//batch_size + 1}: {len(batch)} genomes")

        total_assigned += assigned
        log.info(f"{species_id}: {assigned} lineages assigned")

    cur.close()
    log.info(f"=== Lineage assignment complete: {total_assigned} genomes assigned ===")
    if dry_run:
        log.info("(dry-run: nothing was written to DB)")


def main():
    parser = argparse.ArgumentParser(description="Assign Nextclade lineages to Ebola DB genomes")
    parser.add_argument("--species", type=str, default=None,
                        help="Only process one species (e.g. EBOV)")
    parser.add_argument("--batch-size", type=int, default=100,
                        help="Number of genomes per Nextclade batch")
    parser.add_argument("--reassign", action="store_true",
                        help="Re-assign lineages even if already present")
    parser.add_argument("--dry-run", action="store_true",
                        help="Don't write to DB")
    args = parser.parse_args()

    conn = psycopg2.connect(DB_URL)
    conn.autocommit = False
    try:
        run_assignment(
            conn,
            species_filter=args.species,
            batch_size=args.batch_size,
            reassign=args.reassign,
            dry_run=args.dry_run,
        )
    finally:
        conn.close()


if __name__ == "__main__":
    main()
