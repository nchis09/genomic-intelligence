#!/usr/bin/env python3
"""
Seed gene_function for Ebola virus proteins from UniProt.

For each ebolavirus species, species-specific UniProt entries are fetched
so that functional annotation, domains, PDB IDs, and literature references
match the correct organism.  When no species-specific entry exists, the
Zaire ebolavirus (Mayinga-76) canonical entry is used as a fallback.

Coordinates and protein length are taken from the local reference_proteomes
table so each species keeps its own genome positions.  Domains and functional
sites that fall outside the reference protein length are filtered out, and
any length mismatch between UniProt and reference_proteomes is noted in the
protein_function text.
"""
from __future__ import annotations

import json
import logging
import sys
from datetime import date
from pathlib import Path
from typing import Any

import psycopg2
import requests

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))
from config import DB_URL  # noqa: E402

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s", datefmt="%H:%M:%S")
log = logging.getLogger(__name__)

# Canonical Zaire ebolavirus (strain Mayinga-76) UniProt accessions –
# used as fallback when no species-specific entry is found.
EBOV_UNIPROT_ACCESSIONS = {
    "NP": "P18272",
    "VP35": "Q05127",
    "VP40": "Q05128",
    "GP": "Q05320",
    "VP30": "Q05323",
    "VP24": "Q05322",
    "L": "A9QPM4",
}

# Map species_id → UniProt organism search name.
# UniProt uses "common name" style for some species.
SPECIES_UNIPROT_NAMES = {
    "EBOV":  "Zaire ebolavirus",
    "BDBV":  "Bundibugyo virus",
    "SUDV":  "Sudan ebolavirus",
    "RESTV": "Reston ebolavirus",
    "TAFV":  "Tai Forest ebolavirus",
    "BOMV":  "Bombali virus",
}

# Protein name keywords used to match the correct gene when searching.
GENE_PROTEIN_KEYWORDS = {
    "NP":   "Nucleoprotein",
    "VP35": "Polymerase cofactor VP35",
    "VP40": "Matrix protein VP40",
    "GP":   "Envelope glycoprotein",
    "VP30": "Transcriptional activator VP30",
    "VP24": "Membrane-associated protein VP24",
    "L":    "RNA-directed RNA polymerase L",
}

UNIPROT_API = "https://rest.uniprot.org/uniprotkb"

DOMAIN_FEATURE_TYPES = {"Domain", "Region", "Motif", "Compositional bias", "Topological domain", "Transmembrane", "Signal", "Coiled coil"}
SITE_FEATURE_TYPES = {"Site", "Active site", "Binding site", "Glycosylation", "Lipidation", "Disulfide bond"}


def fetch_uniprot_entry(accession: str) -> dict[str, Any]:
    """Fetch a UniProt entry in JSON; return empty dict on failure."""
    try:
        resp = requests.get(f"{UNIPROT_API}/{accession}", headers={"Accept": "application/json"}, timeout=30)
        resp.raise_for_status()
        return resp.json()
    except Exception as e:
        log.warning("Could not fetch UniProt %s: %s", accession, e)
        return {}


def _extract_text(block: dict[str, Any]) -> str:
    """Concatenate text blocks from a UniProt comment."""
    parts = []
    for t in block.get("texts", []):
        value = t.get("value")
        if value:
            parts.append(value)
    return " ".join(parts)


def _feature_location(feature: dict[str, Any]) -> dict[str, Any]:
    """Return start/end position if available."""
    loc = feature.get("location") or {}
    start = (loc.get("start") or {}).get("value")
    end = (loc.get("end") or {}).get("value")
    return {"start": start, "end": end}


def parse_protein_name(entry: dict[str, Any]) -> str:
    """Return the recommended fullName."""
    try:
        return entry["proteinDescription"]["recommendedName"]["fullName"]["value"]
    except (KeyError, TypeError):
        return entry.get("uniProtkbId", "")


def parse_comment_text(entry: dict[str, Any], comment_type: str) -> str:
    """Join all text for a given comment type."""
    texts = []
    for comment in entry.get("comments", []):
        if comment.get("commentType") == comment_type:
            text = _extract_text(comment)
            if text:
                texts.append(text)
    return " ".join(texts)


def parse_features(entry: dict[str, Any], max_length: int | None = None) -> tuple[list[dict], list[dict]]:
    """Return (key_domains, functional_sites) extracted from feature list.

    If *max_length* is given, features whose end position exceeds it are
    filtered out so that domains/sites never reference positions beyond the
    actual protein stored in reference_proteomes.
    """
    key_domains: list[dict] = []
    functional_sites: list[dict] = []
    for feature in entry.get("features", []):
        ftype = feature.get("type")
        if ftype not in DOMAIN_FEATURE_TYPES and ftype not in SITE_FEATURE_TYPES:
            continue
        loc = _feature_location(feature)
        start = loc["start"]
        end = loc["end"]
        if max_length is not None and end is not None and end > max_length:
            continue
        item = {
            "type": ftype,
            "description": feature.get("description") or "",
            "start": start,
            "end": end,
        }
        if ftype in DOMAIN_FEATURE_TYPES:
            key_domains.append(item)
        else:
            functional_sites.append(item)
    return key_domains, functional_sites


def parse_pdb_ids(entry: dict[str, Any]) -> list[str]:
    """Return PDB identifiers from cross-references."""
    pdb_ids = []
    for ref in entry.get("uniProtKBCrossReferences", []):
        if ref.get("database") == "PDB":
            pdb_id = ref.get("id")
            if pdb_id and pdb_id not in pdb_ids:
                pdb_ids.append(pdb_id)
    return pdb_ids


def parse_literature_refs(entry: dict[str, Any]) -> list[str]:
    """Collect PubMed IDs from references."""
    pmids: set[str] = set()
    for ref in entry.get("references", []):
        citation = ref.get("citation", {})
        for cross in citation.get("citationCrossReferences", []):
            if cross.get("database") == "PubMed" and cross.get("id"):
                pmids.add(cross["id"])
    # Also pick up PMIDs attached to feature evidences
    for feature in entry.get("features", []):
        for ev in feature.get("evidences", []):
            if ev.get("source") == "PubMed" and ev.get("id"):
                pmids.add(ev["id"])
    return sorted(pmids)


def load_species_coordinates(pathogen_id: str = "ebola") -> dict[tuple[str, str], dict[str, Any]]:
    """Return {(species_id, gene): {protein_length, genome_start, genome_end}}."""
    coords: dict[tuple[str, str], dict[str, Any]] = {}
    try:
        conn = psycopg2.connect(DB_URL)
        cur = conn.cursor()
        cur.execute(
            """
            SELECT species_id, gene, protein_length, genome_start, genome_end, protein_name
            FROM reference_proteomes
            WHERE pathogen_id = %s
            """,
            (pathogen_id,),
        )
        for species_id, gene, protein_length, genome_start, genome_end, protein_name in cur.fetchall():
            coords[(species_id, gene)] = {
                "protein_length": protein_length,
                "genome_start": genome_start,
                "genome_end": genome_end,
                "protein_name": protein_name,
            }
        cur.close()
        conn.close()
    except Exception as e:
        log.error("Could not load reference proteome coordinates: %s", e)
    return coords


def search_uniprot_for_gene(organism: str, gene: str, keyword: str) -> dict[str, Any] | None:
    """Search UniProt for a species-specific entry matching *gene*.

    Returns the first reviewed (Swiss-Prot) entry whose protein name
    contains *keyword*, or the first entry if no name match is found.
    Falls back to None if nothing is found.
    """
    import urllib.parse
    query = urllib.parse.quote(f'organism_name:"{organism}" AND {keyword}')
    try:
        resp = requests.get(
            f"{UNIPROT_API}/search?query={query}&format=json&size=10",
            timeout=30,
        )
        resp.raise_for_status()
        data = resp.json()
        entries = data.get("results", [])
        if not entries:
            return None
        # Prefer reviewed (Swiss-Prot) entries
        reviewed = [e for e in entries if "reviewed" in e.get("entryType", "").lower()]
        pool = reviewed or entries
        # Prefer entries whose protein name contains the keyword
        for e in pool:
            pname = parse_protein_name(e).lower()
            if keyword.lower() in pname:
                return e
        return pool[0]
    except Exception as e:
        log.warning("UniProt search failed for %s/%s: %s", organism, gene, e)
        return None


def seed_gene_function(pathogen_id: str = "ebola", dry_run: bool = False) -> dict[str, int]:
    """Fetch species-specific UniProt annotations and upsert gene_function rows."""
    stats = {"fetched": 0, "inserted": 0, "updated": 0, "skipped": 0, "fallback": 0}

    coords = load_species_coordinates(pathogen_id)
    if not coords:
        log.error("No reference proteome coordinates found; aborting.")
        return stats

    # Group coordinates by species_id → {gene: coord}
    species_genes: dict[str, dict[str, dict[str, Any]]] = {}
    for (species_id, gene), coord in coords.items():
        species_genes.setdefault(species_id, {})[gene] = coord

    # Cache: (species_id, gene) → UniProt entry dict
    entry_cache: dict[tuple[str, str], dict[str, Any]] = {}

    for species_id, gene_map in species_genes.items():
        organism = SPECIES_UNIPROT_NAMES.get(species_id, "")
        for gene in gene_map:
            entry: dict[str, Any] | None = None
            if organism:
                keyword = GENE_PROTEIN_KEYWORDS.get(gene, gene)
                entry = search_uniprot_for_gene(organism, gene, keyword)
            if entry:
                stats["fetched"] += 1
                log.info("%s/%s: species-specific entry %s (len=%s)",
                         species_id, gene, entry.get("primaryAccession"),
                         entry.get("sequence", {}).get("length"))
            else:
                # Fallback to EBOV canonical entry
                accession = EBOV_UNIPROT_ACCESSIONS.get(gene)
                if accession:
                    entry = fetch_uniprot_entry(accession)
                if entry:
                    stats["fallback"] += 1
                    log.info("%s/%s: fallback to EBOV entry %s", species_id, gene, accession)
                else:
                    stats["skipped"] += 1
                    log.warning("%s/%s: no UniProt entry found, skipping", species_id, gene)
                    continue
            entry_cache[(species_id, gene)] = entry

    if dry_run:
        for (species_id, gene), entry in entry_cache.items():
            coord = coords[(species_id, gene)]
            ref_len = coord["protein_length"]
            uniprot_len = entry.get("sequence", {}).get("length", 0)
            key_domains, functional_sites = parse_features(entry, max_length=ref_len)
            pdb_ids = parse_pdb_ids(entry)
            mismatch = "" if uniprot_len == ref_len else f" [LENGTH MISMATCH: UniProt={uniprot_len}, ref={ref_len}]"
            log.info("Would seed %s/%s: %s (%d domains, %d sites, %d PDBs)%s",
                     species_id, gene, parse_protein_name(entry),
                     len(key_domains), len(functional_sites), len(pdb_ids), mismatch)
        return stats

    conn = psycopg2.connect(DB_URL)
    cur = conn.cursor()
    try:
        for (species_id, gene), entry in entry_cache.items():
            coord = coords[(species_id, gene)]
            ref_len = coord["protein_length"]
            uniprot_len = entry.get("sequence", {}).get("length", 0)

            protein_name = parse_protein_name(entry) or coord["protein_name"] or gene
            protein_function = parse_comment_text(entry, "FUNCTION")

            # Filter domains/sites to fit within reference protein length
            key_domains, functional_sites = parse_features(entry, max_length=ref_len)
            pdb_ids = parse_pdb_ids(entry)
            literature_refs = parse_literature_refs(entry)

            # Flag length mismatch in protein_function text
            if uniprot_len and uniprot_len != ref_len:
                protein_function = (
                    f"{protein_function}\n\n"
                    f"[NOTE: UniProt entry {entry.get('primaryAccession', '?')} "
                    f"({entry.get('organism', {}).get('scientificName', '?')}) "
                    f"has sequence length {uniprot_len}, but reference_proteomes "
                    f"stores {ref_len} aa. Domains/sites beyond position {ref_len} "
                    f"have been filtered out.]"
                ) if protein_function else (
                    f"[NOTE: UniProt entry {entry.get('primaryAccession', '?')} "
                    f"has sequence length {uniprot_len}, but reference_proteomes "
                    f"stores {ref_len} aa. Domains/sites beyond position {ref_len} "
                    f"have been filtered out.]"
                )

            row = {
                "pathogen_id": pathogen_id,
                "species_id": species_id,
                "gene": gene,
                "protein_name": protein_name,
                "protein_function": protein_function,
                "genome_start": coord["genome_start"],
                "genome_end": coord["genome_end"],
                "protein_length_aa": ref_len,
                "key_domains": json.dumps(key_domains) if key_domains else None,
                "functional_sites": json.dumps(functional_sites) if functional_sites else None,
                "known_hotspots": None,
                "conserved_regions": None,
                "pdb_ids": pdb_ids if pdb_ids else None,
                "literature_refs": literature_refs if literature_refs else None,
                "last_curated": date.today(),
                "curator": "UniProt",
            }

            cols = list(row.keys())
            cur.execute(
                f"""
                INSERT INTO gene_function ({', '.join(cols)})
                VALUES ({', '.join(['%s'] * len(cols))})
                ON CONFLICT (species_id, gene) DO UPDATE SET
                {', '.join(f"{c} = EXCLUDED.{c}" for c in cols if c not in ("species_id", "gene"))}
                """,
                [row[c] for c in cols],
            )
            if cur.rowcount == 1:
                stats["inserted"] += 1
            else:
                stats["updated"] += 1

        conn.commit()
    except Exception as e:
        conn.rollback()
        log.error("Gene function seeding failed: %s", e)
        raise
    finally:
        cur.close()
        conn.close()

    return stats


if __name__ == "__main__":
    import argparse
    import sys

    parser = argparse.ArgumentParser(description="Seed gene_function from UniProt for Ebola proteins")
    parser.add_argument("--dry-run", action="store_true", help="Show what would be inserted without writing to DB")
    args = parser.parse_args()

    stats = seed_gene_function(dry_run=args.dry_run)
    log.info("Gene function seeding complete: %s", stats)
