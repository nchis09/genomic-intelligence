#!/usr/bin/env python3
"""
Seed genotype_phenotype with curated UniProt mutagenesis / natural-variant
features for Ebola virus proteins. This is a deterministic, LLM-free source of
experimentally validated genotype-phenotype associations.

Each feature is checked against the local reference_proteomes table so that
only positions that match the reference sequence are inserted.
"""
from __future__ import annotations

import logging
import re
import sys
from pathlib import Path

import psycopg2
import requests

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from config import DB_URL

sys.path.insert(0, str(Path(__file__).parent))
import extract_from_pubtator as ext

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s", datefmt="%H:%M:%S")
log = logging.getLogger(__name__)

# Canonical Zaire ebolavirus (strain Mayinga-76) UniProt accessions
EBOV_UNIPROT_ACCESSIONS = {
    "NP": "P18272",
    "VP35": "Q05127",
    "VP40": "Q05128",
    "GP": "Q05320",
    "VP30": "Q05323",
    "VP24": "Q05322",
    "L": "A9QPM4",
}

UNIPROT_API = "https://rest.uniprot.org/uniprotkb"

PHENOTYPE_KEYWORDS = [
    ("vaccine_effectiveness", ["vaccine", "vaccination", "immunization"]),
    ("immune_escape", ["antibody", "neutraliz", "immunoglobulin", "escape", "IFN", "antagonize", "antagonist"]),
    ("disease_severity", ["severity", "severe", "mortality", "lethal", "fatal", "virulence", "virulent", "pathogenicity", "pathogen"]),
    ("increased_transmission", ["transmission", "transmissible"]),
    ("drug_resistance", ["resistance", "resistant"]),
    ("drug_susceptibility", ["susceptib", "sensitive", "sensitivity"]),
    ("host_adaptation", ["host adaptation", "adapted to", "mouse-adapted", "guinea pig-adapted", "isolate"]),
    ("host_range", ["host range", "host tropism", "species specificity"]),
    ("virulence", ["replication", "replicative", "release", "entry", "enter", "infectivity", "oligomerization", "trimerization", "processing", "localization", "rna-binding", "binding", "fusion", "synthesis", "interaction", "shedding", "secretion", "counteract", "budding"]),
]


def fetch_uniprot_entry(accession: str) -> dict:
    """Fetch a UniProt entry in JSON; return empty dict on failure."""
    url = f"{UNIPROT_API}/{accession}"
    try:
        resp = requests.get(url, headers={"Accept": "application/json"}, timeout=30)
        resp.raise_for_status()
        return resp.json()
    except Exception as e:
        log.warning("Could not fetch UniProt %s: %s", accession, e)
        return {}


def description_to_phenotype(description: str) -> str:
    """Map a UniProt feature description to a phenotype_category."""
    desc_lower = description.lower()
    for phen, keywords in PHENOTYPE_KEYWORDS:
        if any(kw in desc_lower for kw in keywords):
            return phen
    return "unknown_significance"


def parse_literature_refs(feature: dict) -> list[str]:
    """Extract PMID refs from feature evidences."""
    refs = []
    for ev in feature.get("evidences", []):
        if ev.get("source") == "PubMed" and ev.get("id"):
            refs.append(f"PMID:{ev['id']}")
    return refs


def seed_uniprot(pathogen_id: str = "ebola", species_id: str = "EBOV", dry_run: bool = False) -> dict:
    ref_proteomes = ext.load_reference_proteomes(pathogen_id)
    stats = {"mutagenesis": 0, "natural_variant": 0, "skipped": 0}

    conn = psycopg2.connect(DB_URL)
    cur = conn.cursor()

    for protein, accession in EBOV_UNIPROT_ACCESSIONS.items():
        entry = fetch_uniprot_entry(accession)
        if not entry:
            continue
        sequence = entry.get("sequence", {}).get("value", "")
        ref_seq = ref_proteomes.get((species_id, protein), "")
        if not ref_seq:
            log.warning("No reference sequence for %s/%s; skipping validation", species_id, protein)

        for feature in entry.get("features", []):
            ftype = feature.get("type")
            if ftype not in ("Mutagenesis", "Natural variant"):
                continue

            loc = feature.get("location", {}).get("start", {})
            pos = loc.get("value")
            if not pos or not isinstance(pos, int):
                continue

            alt_seq = feature.get("alternativeSequence", {})
            ref_aa = alt_seq.get("originalSequence")
            alts = alt_seq.get("alternativeSequences", [])
            description = feature.get("description") or ""

            # Skip entries that explicitly report no phenotype
            if re.search(r"\bno effect\b|\bnot affect\b|\bdoes not affect\b|\bno loss\b", description, re.IGNORECASE):
                continue

            # Validate against local reference proteome if available
            if ref_seq:
                if pos < 1 or pos > len(ref_seq):
                    log.debug("Skipping %s %s%d: outside reference length", protein, ref_aa, pos)
                    stats["skipped"] += 1
                    continue
                if ref_aa and ref_seq[pos - 1].upper() != ref_aa.upper():
                    log.warning(
                        "Reference mismatch for %s %s%d (UniProt says %s, reference has %s); skipping",
                        protein, ref_aa, pos, ref_aa, ref_seq[pos - 1],
                    )
                    stats["skipped"] += 1
                    continue

            category = description_to_phenotype(description)
            refs = parse_literature_refs(feature)

            for alt_aa in alts:
                if not alt_aa or alt_aa == "*":
                    continue
                genotype_desc = f"{protein} {ref_aa}{pos}{alt_aa}"
                aid = ext.association_id(
                    pathogen_id,
                    species_id,
                    protein,
                    pos,
                    ref_aa,
                    alt_aa,
                    f"uniprot:{genotype_desc}",
                    category,
                )

                if dry_run:
                    log.info("Would insert %s -> %s", genotype_desc, category)
                    stats[ftype.lower().replace(" ", "_")] += 1
                    continue

                # Skip if already seeded
                cur.execute("SELECT 1 FROM genotype_phenotype WHERE association_id = %s LIMIT 1", (aid,))
                if cur.fetchone():
                    continue

                cur.execute(
                    """
                    INSERT INTO genotype_phenotype (
                        association_id, pathogen_id, species_id, lineage_id,
                        protein, position, ref_aa, alt_aa,
                        genotype_description, phenotype_category, phenotype_specific,
                        effect_size, evidence_strength,
                        literature_refs, source_url, notes,
                        record_flagged, flag_reason,
                        verification_status, data_source, last_updated, ingested_at
                    ) VALUES (
                        %s, %s, %s, %s,
                        %s, %s, %s, %s,
                        %s, %s, %s,
                        %s, %s, %s,
                        %s, %s, %s,
                        %s, %s,
                        %s, %s, %s, %s
                    )
                    """,
                    (
                        aid,
                        pathogen_id,
                        species_id,
                        None,
                        protein,
                        pos,
                        ref_aa,
                        alt_aa,
                        genotype_desc,
                        category,
                        description,
                        "not quantified",
                        "strong" if ftype == "Mutagenesis" else "moderate",
                        refs,
                        f"https://www.uniprot.org/uniprotkb/{accession}",
                        f"Imported from UniProt {accession}. {description}",
                        False,
                        "",
                        "verified",
                        "uniprot",
                        ext.date.today(),
                        ext.date.today(),
                    ),
                )
                stats[ftype.lower().replace(" ", "_")] += 1

    if not dry_run:
        conn.commit()
    cur.close()
    conn.close()
    return stats


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Seed genotype_phenotype from UniProt mutagenesis features")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    stats = seed_uniprot(dry_run=args.dry_run)
    log.info("UniProt seed complete: %s", stats)
