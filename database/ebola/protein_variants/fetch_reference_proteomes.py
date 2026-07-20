#!/usr/bin/env python3
"""
Fetch the six canonical Ebola species reference genomes and their ORFs from NCBI,
then populate the reference_genomes and reference_proteomes tables.

Run with anaconda python (which has psycopg2 and biopython):
    /Users/christianndekezi/anaconda3/bin/python3 database/ebola/protein_variants/fetch_reference_proteomes.py
"""

import json
import logging
import sys
from pathlib import Path

from Bio import Entrez, SeqIO
import psycopg2

# Use project config for DB URL
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))
from config import DB_URL

Entrez.email = "pgirl_pipeline@local"  # NCBI requires an email

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s", datefmt="%H:%M:%S")
log = logging.getLogger(__name__)

EBOLA_SPECIES = [
    {"species_id": "EBOV", "species_name": "Zaire ebolavirus", "common_name": "Ebola virus", "accession": "NC_002549.1"},
    {"species_id": "SUDV", "species_name": "Sudan ebolavirus", "common_name": "Sudan virus", "accession": "NC_006432.1"},
    {"species_id": "BDBV", "species_name": "Bundibugyo ebolavirus", "common_name": "Bundibugyo virus", "accession": "NC_014373.1"},
    {"species_id": "RESTV", "species_name": "Reston ebolavirus", "common_name": "Reston virus", "accession": "NC_004161.1"},
    {"species_id": "TAFV", "species_name": "Tai Forest ebolavirus", "common_name": "Tai Forest virus", "accession": "NC_014372.1"},
    {"species_id": "BOMV", "species_name": "Bombali ebolavirus", "common_name": "Bombali virus", "accession": "NC_039345.1"},
]


def fetch_refseq(accession: str):
    """Fetch a GenBank record from NCBI by accession."""
    log.info(f"Fetching {accession} from NCBI...")
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()
    return record


def parse_cds_features(record):
    """Return a list of gene/protein info dicts from GenBank CDS features."""
    proteins = []
    for feature in record.features:
        if feature.type != "CDS":
            continue
        gene = feature.qualifiers.get("gene", [""])[0]
        protein_id = feature.qualifiers.get("protein_id", [""])[0]
        product = feature.qualifiers.get("product", [""])[0]
        translation = feature.qualifiers.get("translation", [""])[0]
        location = feature.location
        if not gene or not translation:
            continue
        proteins.append({
            "gene": gene,
            "protein_id": protein_id,
            "protein_name": product,
            "protein_sequence": translation,
            "genome_start": int(location.start) + 1,  # 0-based to 1-based
            "genome_end": int(location.end),
            "strand": int(location.strand) if location.strand else 1,
        })
    return proteins


def ensure_pathogen(cur):
    cur.execute("SELECT 1 FROM pathogens WHERE pathogen_id = 'ebola'")
    if cur.fetchone() is None:
        cur.execute(
            "INSERT INTO pathogens (pathogen_id, family, genus, ncbi_taxonomy_id, notes, last_curated) VALUES (%s, %s, %s, %s, %s, CURRENT_DATE)",
            ("ebola", "Filoviridae", "Orthoebolavirus", 186536, "6 species: EBOV, SUDV, BDBV, RESTV, TAFV, BOMV"),
        )
        log.info("Inserted pathogens row for ebola")


def ensure_species(cur, sp):
    cur.execute("SELECT 1 FROM species WHERE species_id = %s", (sp["species_id"],))
    if cur.fetchone() is None:
        cur.execute(
            """
            INSERT INTO species (species_id, pathogen_id, species_name, common_name, abbreviation, ncbi_refseq_accession, biosafety_level, human_pathogen, last_curated)
            VALUES (%s, %s, %s, %s, %s, %s, %s, %s, CURRENT_DATE)
            """,
            (sp["species_id"], "ebola", sp["species_name"], sp["common_name"], sp["species_id"], sp["accession"], 4, True),
        )
        log.info(f"Inserted species {sp['species_id']}")


def load_reference(cur, sp, record, proteins):
    # Reference genome record
    gene_coords = {p["gene"]: {"start": p["genome_start"], "end": p["genome_end"], "strand": p["strand"]} for p in proteins}
    cur.execute(
        """
        INSERT INTO reference_genomes (accession, pathogen_id, species_id, genome_role, genome_length, segmented, source_database, gene_coordinates, notes, last_curated)
        VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, CURRENT_DATE)
        ON CONFLICT (accession) DO UPDATE SET
            genome_length = EXCLUDED.genome_length,
            gene_coordinates = EXCLUDED.gene_coordinates,
            last_curated = EXCLUDED.last_curated
        """,
        (
            sp["accession"],
            "ebola",
            sp["species_id"],
            "canonical_reference",
            len(record.seq),
            False,
            "NCBI",
            json.dumps(gene_coords),
            f"Canonical reference for {sp['species_name']}",
        ),
    )

    # Proteome records
    for p in proteins:
        cur.execute(
            """
            INSERT INTO reference_proteomes (reference_accession, species_id, pathogen_id, gene, protein_name, protein_sequence, genome_start, genome_end, strand, last_curated)
            VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, CURRENT_DATE)
            ON CONFLICT (reference_accession, gene) DO UPDATE SET
                protein_sequence = EXCLUDED.protein_sequence,
                protein_name = EXCLUDED.protein_name,
                genome_start = EXCLUDED.genome_start,
                genome_end = EXCLUDED.genome_end,
                last_curated = EXCLUDED.last_curated
            """,
            (
                sp["accession"],
                sp["species_id"],
                "ebola",
                p["gene"],
                p["protein_name"],
                p["protein_sequence"],
                p["genome_start"],
                p["genome_end"],
                p["strand"],
            ),
        )
    log.info(f"Loaded {len(proteins)} reference proteins for {sp['species_id']}")


def main():
    conn = psycopg2.connect(DB_URL)
    conn.autocommit = False
    cur = conn.cursor()

    try:
        ensure_pathogen(cur)
        for sp in EBOLA_SPECIES:
            ensure_species(cur, sp)
            record = fetch_refseq(sp["accession"])
            proteins = parse_cds_features(record)
            load_reference(cur, sp, record, proteins)
            conn.commit()
        log.info("Reference proteomes loaded successfully.")
    except Exception as e:
        conn.rollback()
        log.error(f"Failed to load reference proteomes: {e}")
        raise
    finally:
        cur.close()
        conn.close()


if __name__ == "__main__":
    main()
