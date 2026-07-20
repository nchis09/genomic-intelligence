#!/usr/bin/env python3
"""
Extract candidate genotype-phenotype associations from PubTator-annotated papers.

Pipeline:
  1. Search PubMed for Ebola papers that mention genotype and phenotype terms.
  2. Collect PMIDs and fetch their PubTator3 annotations.
  3. Extract variants, genes, species, and clade/strain/motif signals from text.
  4. Detect phenotype keywords from title + abstract.
  5. Insert candidate associations into genotype_phenotype as `unverified` and
     `record_flagged=True` for curator review.

PubTator handles the mutation/gene entity recognition, so the downstream
extraction is deterministic and stays aligned with the no-AI/ML engine rule.
"""

from __future__ import annotations

import argparse
import logging
import re
import sys
import time
import uuid
from collections import defaultdict
from datetime import date
from pathlib import Path
from typing import Optional

import psycopg2
import requests

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from config import DB_URL

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s", datefmt="%H:%M:%S")
log = logging.getLogger(__name__)

PUBTATOR_BIOCJSON_URL = "https://www.ncbi.nlm.nih.gov/research/pubtator-api/publications/export/biocjson"
PUBMED_ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"

# Pathogen / species mapping for the DB
EBOLA_SPECIES = {
    "Zaire ebolavirus": "EBOV",
    "Sudan ebolavirus": "SUDV",
    "Bundibugyo ebolavirus": "BDBV",
    "Reston ebolavirus": "RESTV",
    "Tai Forest ebolavirus": "TAFV",
    "Bombali ebolavirus": "BOMV",
}

# Genotype keywords used to enrich the PubMed search
GENOTYPE_SEARCH_TERMS = [
    "A82V", "R141G", "D531N", "I242T", "V95I", "T142A",
    "glycan cap", "mucin-like domain", "receptor-binding site",
    "fusion loop", "internal fusion loop", "NP-VP35",
    "Makona", "Mayinga", "Kikwit", "Gulu", "Mubende",
    "EBOV-2014", "EBOV-2018",
    "GP mutation", "NP mutation", "VP24 mutation", "VP35 mutation",
    "monoclonal antibody", "antibody resistance", "neutralization",
    "vaccine efficacy", "vaccine response", "drug resistance",
]

# Phenotype keywords for PubMed search
PHENOTYPE_SEARCH_TERMS = [
    "transmission", "virulence", "pathogenicity", "severity", "mortality",
    "immune escape", "antibody resistance", "neutralization", "vaccine efficacy",
    "drug resistance", "antiviral resistance", "host range", "replication fitness",
    "glycan cap", "monoclonal antibody", "diagnostic",
]

# Phenotype detection patterns used on title+abstract text
PHENOTYPE_PATTERNS = {
    "increased_transmission": [r"transmiss(ion|ibility)", r"\bspread\b", r"\bR0\b"],
    "decreased_transmission": [r"reduced\s+transmission", r"lower\s+transmission"],
    "virulence": [r"virulence", r"pathogenicity", r"virulent", r"more\s+virulent", r"increased\s+virulence", r"decreased\s+virulence"],
    "disease_severity": [r"severity", r"severe", r"mortality", r"case\s+fatality", r"\bCFR\b", r"\bfatal\w*\b"],
    "host_adaptation": [r"human\s+adaptation", r"host\s+adaptation", r"adaptation\s+to\s+human"],
    "host_range": [r"host\s+range", r"host\s+tropism", r"species\s+tropism"],
    "reservoir_association": [r"reservoir", r"\bbat\b", r"frugivorous"],
    "vector_competence": [r"vector\s+competence", r"vector\s+transmission"],
    "immune_escape": [r"immune\s+(escape|evasion)", r"antibody\s+(resist|resistance|escape)", r"neutrali[sz](ation|ing|e)"],
    "vaccine_escape": [r"vaccine\s+(escape|resistance)", r"\bErvebo\b", r"rVSV-ZEBOV"],
    "vaccine_effectiveness": [r"vaccine\s+(efficacy|effectiveness|response|induced)", r"vaccination"],
    "drug_resistance": [r"drug\s+(resist|resistance)", r"antiviral\s+resistance", r"\bremdesivir\b", r"\bfavipiravir\b"],
    "drug_susceptibility": [r"drug\s+susceptib", r"susceptib\w*\s+to", r"antiviral\s+activity", r"\binhibited\b"],
    "diagnostic_escape": [r"diagnostic\s+(escape|failure)", r"assay\s+failure", r"false\s+negative"],
    "diagnostic_performance": [r"diagnostic\s+assay", r"\bPCR\b", r"RT-PCR"],
    "unknown_significance": [r"unknown\s+significance", r"no\s+known\s+significance"],
}

# Clade / strain / motif detection
CLADE_STRAIN_NAMES = [
    "Makona", "Mayinga", "Kikwit", "Gulu", "Mubende",
    "EBOV-2014", "EBOV-2018", "EBOV-2013", "EBOV-1976", "EBOV-1995",
    "B.1", "A.3", "A.4", "C15", "C07", "C14",
]

MOTIF_PATTERNS = [
    "glycan cap", "mucin-like domain", "mucin-like", "receptor-binding site",
    "receptor binding site", "RBS", "fusion loop", "internal fusion loop",
    "IFL", "base", "head", "chalice", "GP1", "GP2", "GP1/GP2",
    "NPC1 binding", "cathepsin cleavage", "furin cleavage", "cleavage site",
    "membrane fusion", "matrix layer", "ribonucleoprotein", "NP-VP35",
    "VP24 interferon", "VP35 IFN", "VP35 interferon", "L polymerase",
    "RNA polymerase", "editing site", "transcriptional editing",
]



# ---------------------------------------------------------------------------
# Protein/species inference helpers
# ---------------------------------------------------------------------------

# Map a motif phrase (from MOTIF_PATTERNS) to its canonical protein.
# Ambiguous motifs (NP-VP35, ribonucleoprotein) may be overridden by context.
MOTIF_TO_PROTEIN = {
    "glycan cap": "GP",
    "mucin-like domain": "GP",
    "mucin-like": "GP",
    "receptor-binding site": "GP",
    "receptor binding site": "GP",
    "RBS": "GP",
    "fusion loop": "GP",
    "internal fusion loop": "GP",
    "IFL": "GP",
    "base": "GP",
    "head": "GP",
    "chalice": "GP",
    "GP1": "GP",
    "GP2": "GP",
    "GP1/GP2": "GP",
    "NPC1 binding": "GP",
    "NPC1": "GP",
    "cathepsin cleavage": "GP",
    "furin cleavage": "GP",
    "cleavage site": "GP",
    "membrane fusion": "GP",
    "editing site": "GP",
    "transcriptional editing": "GP",
    "matrix layer": "VP40",
    "matrix": "VP40",
    "ribonucleoprotein": "NP",
    "NP-VP35": "NP",
    "VP24 interferon": "VP24",
    "VP35 IFN": "VP35",
    "VP35 interferon": "VP35",
    "L polymerase": "L",
    "RNA polymerase": "L",
}

# Context patterns that indicate the variant is in an antibody, not a viral protein.
ANTIBODY_CONTEXT_RE = re.compile(
    r"\b(light[- ]chain|heavy[- ]chain|CDR|framework|complementarity[- ]determining|framework[- ]region|FR[- ]grafting|ADI-\d+)\b",
    re.IGNORECASE,
)
ANTIBODY_VARIANT_RE = re.compile(r"[A-Z]\d+[A-Z]\s*[-–]\s*(LC|HC|VL|VH)\b", re.IGNORECASE)

# Canonical Ebola proteins and text aliases used to find nearby mentions.
PROTEIN_ALIASES = {
    "GP": ["GP", "glycoprotein", "envelope glycoprotein", "sGP", "GP1", "GP2",
           "glycan cap", "mucin-like", "mucin-like domain", "mucin domain",
           "receptor-binding site", "receptor binding site", "RBS",
           "fusion loop", "internal fusion loop", "IFL",
           "NPC1 binding", "NPC1", "cathepsin cleavage", "furin cleavage",
           "cleavage site", "membrane fusion", "editing site",
           "transcriptional editing", "base", "head", "chalice"],
    "NP": ["NP", "nucleoprotein", "ribonucleoprotein", "RNP", "ribonucleoprotein complex"],
    "VP35": ["VP35", "VP35 IFN", "VP35 interferon", "NP-VP35 interface", "NP-VP35 interaction"],
    "VP40": ["VP40", "matrix layer", "matrix protein"],
    "VP30": ["VP30"],
    "VP24": ["VP24", "VP24 interferon", "VP24 interferon inhibitory domain"],
    "L": ["L polymerase", "RNA polymerase", "RNA-dependent RNA polymerase",
          "RdRp", "polymerase", "L protein", "L gene"],
}

SPECIES_ALIASES = {
    "EBOV": ["Zaire ebolavirus", "Ebola virus", "EBOV", "ZEBOV", "Zaire"],
    "SUDV": ["Sudan ebolavirus", "Sudan virus", "SUDV", "SEBOV", "Sudan"],
    "BDBV": ["Bundibugyo ebolavirus", "Bundibugyo virus", "BDBV"],
    "RESTV": ["Reston ebolavirus", "Reston virus", "RESTV", "REBOV", "Reston"],
    "TAFV": ["Tai Forest ebolavirus", "Tai Forest", "TAFV"],
    "BOMV": ["Bombali ebolavirus", "Bombali", "BOMV"],
}

CONTEXT_WINDOW = 120

def build_pubmed_query(species: str, genotype_terms: list[str], phenotype_terms: list[str]) -> str:
    """Build a PubMed esearch query that is broad but phenotype-enriched."""
    species_terms = [species, "Ebolavirus", "EBOV", "Ebola virus"]
    species_query = " OR ".join(f'"{t}"[Title/Abstract]' for t in species_terms)
    genotype_query = " OR ".join(f'"{t}"[Title/Abstract]' for t in genotype_terms)
    phenotype_query = " OR ".join(f'"{t}"[Title/Abstract]' for t in phenotype_terms)
    return f"({species_query}) AND ({genotype_query}) AND ({phenotype_query})"


def search_pubmed_pmids(
    query: str,
    max_results: int = 50,
    api_key: Optional[str] = None,
) -> list[str]:
    """Run PubMed esearch and return a list of PMIDs."""
    params = {
        "db": "pubmed",
        "term": query,
        "retmax": max_results,
        "sort": "pub_date",
        "retmode": "json",
    }
    if api_key:
        params["api_key"] = api_key

    try:
        resp = requests.get(PUBMED_ESEARCH_URL, params=params, timeout=30)
        resp.raise_for_status()
        data = resp.json()
        ids = data.get("esearchresult", {}).get("idlist", [])
        return ids
    except Exception as e:
        log.error("PubMed esearch failed: %s", e)
        return []


def fetch_pubtator_annotations(pmids: list[str]) -> list[dict]:
    """Fetch PubTator3 BiocJSON annotations for a list of PMIDs."""
    if not pmids:
        return []
    chunk_size = 100
    all_pubs = []
    for i in range(0, len(pmids), chunk_size):
        chunk = pmids[i:i + chunk_size]
        ids = ",".join(chunk)
        url = f"{PUBTATOR_BIOCJSON_URL}?pmids={ids}"
        try:
            resp = requests.get(url, timeout=60)
            resp.raise_for_status()
            data = resp.json()
            pubs = data.get("PubTator3", [])
            all_pubs.extend(pubs)
        except Exception as e:
            log.error("PubTator request failed for %s: %s", ids[:100], e)
        time.sleep(0.5)
    return all_pubs


def parse_variants_from_pubtator(pub: dict) -> list[dict]:
    """Extract variant annotations from a PubTator publication.
    Each variant includes the global offset of its first location so callers
    can resolve the local context window.
    """
    variants = []
    seen = set()
    for passage in pub.get("passages", []):
        for ann in passage.get("annotations", []):
            infons = ann.get("infons", {})
            if infons.get("type") != "Variant":
                continue
            text = ann.get("text", "")
            locations = ann.get("locations", [])
            if not locations:
                continue
            offset = locations[0].get("offset", 0)
            hgvs = infons.get("hgvs") or infons.get("HGVS") or ""
            parsed = parse_hgvs(hgvs)
            if parsed:
                key = (parsed["ref_aa"], parsed["position"], parsed["alt_aa"])
                if key in seen:
                    continue
                seen.add(key)
                parsed["text"] = text
                parsed["offset"] = offset
                variants.append(parsed)
    return variants
def parse_hgvs(hgvs: str) -> Optional[dict]:
    """Parse HGVS p. strings like p.A82V or p.Gly82Val."""
    if not hgvs.startswith("p."):
        return None
    m = re.match(r"p\.([A-Z][a-z]{2}|[A-Z])(\d+)([A-Z][a-z]{2}|[A-Z])", hgvs)
    if not m:
        return None
    ref = m.group(1)
    pos = int(m.group(2))
    alt = m.group(3)
    return {
        "ref_aa": one_letter(ref),
        "position": pos,
        "alt_aa": one_letter(alt),
    }


def parse_tmvar(tmvar: str) -> Optional[dict]:
    """Parse tmVar identifier like tmVar:p|SUB|A|82|V."""
    parts = tmvar.split("|")
    if len(parts) < 5 or parts[0] != "tmVar:p":
        return None
    if parts[1] == "SUB":
        return {"ref_aa": parts[2], "position": int(parts[3]), "alt_aa": parts[4]}
    return None


AA3 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
}


def one_letter(aa: str) -> str:
    if len(aa) == 1:
        return aa
    return AA3.get(aa[0].upper() + aa[1:].lower(), aa)


def parse_genes_from_pubtator(pub: dict) -> list[dict]:
    """Extract gene/protein annotations with their global offsets."""
    genes = []
    for passage in pub.get("passages", []):
        for ann in passage.get("annotations", []):
            infons = ann.get("infons", {})
            if infons.get("type") != "Gene":
                continue
            locations = ann.get("locations", [])
            offset = locations[0].get("offset", 0) if locations else None
            genes.append({
                "text": ann.get("text", ""),
                "gene_id": infons.get("identifier", ""),
                "name": infons.get("name", "") or ann.get("text", ""),
                "offset": offset,
            })
    return genes
def get_publication_text(pub: dict) -> str:
    """Concatenate title and abstract."""
    parts = []
    for passage in pub.get("passages", []):
        text = passage.get("text", "")
        if text:
            parts.append(text)
    return " ".join(parts)


def detect_phenotypes(text: str) -> list[tuple[str, str]]:
    """Detect phenotype categories from text."""
    text_lower = text.lower()
    found = []
    for term, patterns in PHENOTYPE_PATTERNS.items():
        for pat in patterns:
            if re.search(pat, text_lower, re.IGNORECASE):
                found.append((term, pat))
                break
    return found


def detect_clades(text: str) -> list[str]:
    found = []
    for name in CLADE_STRAIN_NAMES:
        if re.search(rf"\b{re.escape(name)}\b", text, re.IGNORECASE):
            found.append(name)
    return found


def detect_motifs(text: str) -> list[str]:
    found = []
    for motif in MOTIF_PATTERNS:
        if re.search(rf"\b{re.escape(motif)}\b", text, re.IGNORECASE):
            found.append(motif)
    return found


def load_lineage_aliases() -> dict[str, str]:
    """Load lineage_id by common/alias names from the lineages table."""
    aliases: dict[str, str] = {}
    try:
        conn = psycopg2.connect(DB_URL)
        cur = conn.cursor()
        cur.execute(
            "SELECT lineage_id, known_aliases FROM lineages WHERE known_aliases IS NOT NULL"
        )
        for lineage_id, known_aliases in cur.fetchall():
            for alias in known_aliases or []:
                aliases[alias.lower().strip()] = lineage_id
            # Also allow matching by the technical lineage_id itself
            aliases[lineage_id.lower()] = lineage_id
        cur.close()
        conn.close()
    except Exception as e:
        log.warning("Could not load lineage aliases: %s", e)
    return aliases


def resolve_lineage(clade: str, aliases: dict[str, str]) -> Optional[str]:
    """Map a clade/strain name to a lineage_id. Returns None if unknown."""
    if not clade:
        return None
    return aliases.get(clade.lower().strip())


def load_reference_proteomes(pathogen_id: str = "ebola") -> dict[tuple[str, str], str]:
    """Load reference proteome sequences keyed by (species_id, gene)."""
    ref: dict[tuple[str, str], str] = {}
    try:
        conn = psycopg2.connect(DB_URL)
        cur = conn.cursor()
        cur.execute(
            "SELECT species_id, gene, protein_sequence FROM reference_proteomes WHERE pathogen_id = %s",
            (pathogen_id,),
        )
        for species_id, gene, sequence in cur.fetchall():
            ref[(species_id, gene)] = sequence
        cur.close()
        conn.close()
    except Exception as e:
        log.warning("Could not load reference proteomes: %s", e)
    return ref


def get_context_window(text: str, offset: int, target: str, window: int = CONTEXT_WINDOW) -> str:
    """Return a slice of text around the given global offset, padded safely."""
    start = max(0, offset - window)
    end = min(len(text), offset + len(target) + window)
    return text[start:end]


def is_antibody_context(window: str, variant_text: str = "") -> bool:
    """Return True if the variant is in an antibody-engineering context."""
    if not window:
        return False
    if ANTIBODY_VARIANT_RE.search(window):
        return True
    if ANTIBODY_CONTEXT_RE.search(window):
        return True
    if variant_text and re.search(
        rf"{re.escape(variant_text)}\s*[-–]\s*(LC|HC|VL|VH)\b", window, re.IGNORECASE
    ):
        return True
    return False


def _closest_protein_in_window(window: str, target_center: float) -> Optional[str]:
    """Return the canonical protein whose alias appears closest to target_center."""
    best = None
    best_dist = None
    for protein, aliases in PROTEIN_ALIASES.items():
        for alias in aliases:
            pattern = rf"\b{re.escape(alias)}\b" if len(alias) >= 2 else re.escape(alias)
            for m in re.finditer(pattern, window, re.IGNORECASE):
                dist = abs(m.start() - target_center)
                if best_dist is None or dist < best_dist:
                    best_dist = dist
                    best = protein
    return best


def resolve_protein_from_context(text: str, offset: int, target: str) -> Optional[str]:
    """Find the nearest protein alias in the context window and return canonical protein."""
    if not text:
        return None
    window = get_context_window(text, offset, target, window=CONTEXT_WINDOW)
    window_start = max(0, offset - CONTEXT_WINDOW)
    target_center = (offset - window_start) + len(target) / 2.0
    return _closest_protein_in_window(window, target_center)


def resolve_species_from_context(text: str, offset: int, target: str, default_species_id: str) -> str:
    """Use species mentions in text to choose species_id."""
    species, _ = _closest_species_distance(text, offset, target)
    return species or default_species_id


def _closest_species_distance(text: str, offset: int, target: str) -> tuple[Optional[str], Optional[float]]:
    """Return the closest species mention to the target and the distance."""
    if not text:
        return None, None
    window = get_context_window(text, offset, target, window=CONTEXT_WINDOW)
    window_start = max(0, offset - CONTEXT_WINDOW)
    target_center = (offset - window_start) + len(target) / 2.0
    best = None
    best_dist = None
    for species_id, aliases in SPECIES_ALIASES.items():
        for alias in aliases:
            for m in re.finditer(rf"\b{re.escape(alias)}\b", text, re.IGNORECASE):
                dist = abs(m.start() - target_center)
                if best_dist is None or dist < best_dist:
                    best_dist = dist
                    best = species_id
    return best, best_dist


def find_reference_matches(
    ref_aa: str, position: int, ref_proteomes: dict[tuple[str, str], str]
) -> list[tuple[str, str]]:
    """Return all (species_id, gene) where the reference has ref_aa at position (1-indexed)."""
    matches = []
    if not ref_proteomes or not ref_aa or not position:
        return matches
    for (species_id, gene), seq in ref_proteomes.items():
        if position < 1 or position > len(seq):
            continue
        if seq[position - 1].upper() == ref_aa.upper():
            matches.append((species_id, gene))
    return matches


def resolve_variant_protein_and_species(
    var: dict,
    text: str,
    ref_proteomes: dict[tuple[str, str], str],
    default_species_id: str,
) -> tuple[Optional[str], str, str]:
    """Resolve protein, species_id, and a reason for a variant."""
    variant_text = var.get("text") or f"{var['ref_aa']}{var['position']}{var['alt_aa']}"
    offset = var.get("offset", 0)
    window = get_context_window(text, offset, variant_text, window=CONTEXT_WINDOW)
    if is_antibody_context(window, variant_text):
        return None, default_species_id, "antibody variant"
    if not ref_proteomes:
        protein = resolve_protein_from_context(text, offset, variant_text)
        return protein, default_species_id, "no reference data"
    matches = find_reference_matches(var["ref_aa"], var["position"], ref_proteomes)
    if not matches:
        protein = resolve_protein_from_context(text, offset, variant_text)
        return (
            protein,
            default_species_id,
            f"ref_aa {var['ref_aa']}@{var['position']} not in reference",
        )
    # Determine species from context and matches
    species_ids = sorted({s for s, _ in matches})
    if len(species_ids) == 1:
        species_id = species_ids[0]
    else:
        species_id = resolve_species_from_context(text, offset, variant_text, default_species_id)
        if species_id not in species_ids:
            # Fall back to the matching species with the closest context mention
            best = None
            best_dist = None
            for s in species_ids:
                for alias in SPECIES_ALIASES.get(s, []):
                    for m in re.finditer(rf"\b{re.escape(alias)}\b", text, re.IGNORECASE):
                        dist = abs(m.start() - offset)
                        if best_dist is None or dist < best_dist:
                            best_dist = dist
                            best = s
            species_id = best or species_ids[0]
    # Candidate proteins for the selected species
    candidates = sorted({g for s, g in matches if s == species_id})
    protein = resolve_protein_from_context(text, offset, variant_text)
    if protein and protein not in candidates:
        # Context protein doesn't validate against reference; choose closest candidate
        best = None
        best_dist = None
        for c in candidates:
            for alias in PROTEIN_ALIASES.get(c, []):
                for m in re.finditer(rf"\b{re.escape(alias)}\b", text, re.IGNORECASE):
                    dist = abs(m.start() - offset)
                    if best_dist is None or dist < best_dist:
                        best_dist = dist
                        best = c
        protein = best or candidates[0]
    elif not protein:
        protein = candidates[0] if len(candidates) == 1 else resolve_protein_from_context(text, offset, variant_text)
        if not protein:
            protein = candidates[0]
    return protein, species_id, "resolved from reference and context"


def resolve_motif_protein(
    motif: str, text: str, ref_proteomes: dict[tuple[str, str], str]
) -> Optional[str]:
    """Resolve protein for a motif, using curated map plus context for ambiguous motifs."""
    protein = MOTIF_TO_PROTEIN.get(motif)
    if not text:
        return protein
    for m in re.finditer(re.escape(motif), text, re.IGNORECASE):
        offset = m.start()
        target = m.group(0)
        if motif in ("NP-VP35", "ribonucleoprotein"):
            ctx_protein = resolve_protein_from_context(text, offset, target)
            if ctx_protein:
                return ctx_protein
        break
    return protein
def association_id(pathogen_id: str, species_id: str, protein: Optional[str], position: Optional[int], ref: Optional[str], alt: Optional[str], description: str, phenotype: str) -> str:
    tokens = [
        pathogen_id, species_id,
        protein or "", str(position or ""), ref or "", alt or "",
        description, phenotype,
    ]
    return "gpa_" + uuid.uuid5(uuid.NAMESPACE_URL, "|".join(tokens)).hex[:12]


def build_candidates(
    pub: dict,
    species_id: str = "EBOV",
    pathogen_id: str = "ebola",
    lineage_aliases: Optional[dict[str, str]] = None,
    ref_proteomes: Optional[dict[tuple[str, str], str]] = None,
) -> list[dict]:
    """Build candidate genotype-phenotype associations from a PubTator publication.
    Variants and motifs are resolved to a protein (when evidence is strong) using
    local context + reference proteome validation; antibody variants are skipped.
    """
    pmid = pub.get("id", "")
    text = get_publication_text(pub)
    phenotypes = detect_phenotypes(text)
    if not phenotypes:
        return []

    variants = parse_variants_from_pubtator(pub)
    motifs = detect_motifs(text)
    clades = detect_clades(text)

    refs = [f"PMID:{pmid}"]
    source_url = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid else ""
    candidates = []
    seen = set()

    if ref_proteomes is None:
        ref_proteomes = load_reference_proteomes(pathogen_id)

    for var in variants:
        variant_text = var.get("text") or f"{var['ref_aa']}{var['position']}{var['alt_aa']}"
        protein, resolved_species, reason = resolve_variant_protein_and_species(
            var, text, ref_proteomes, species_id
        )
        if reason == "antibody variant":
            log.info("Skipping antibody variant %s in PMID:%s", variant_text, pmid)
            continue
        for term, _ in phenotypes:
            aid = association_id(
                pathogen_id,
                resolved_species,
                protein,
                var["position"],
                var["ref_aa"],
                var["alt_aa"],
                variant_text,
                term,
            )
            if aid in seen:
                continue
            seen.add(aid)
            notes = f"Extracted via PubTator from PMID:{pmid}. {reason}. Flagged for review."
            candidates.append(
                {
                    "association_id": aid,
                    "pathogen_id": pathogen_id,
                    "species_id": resolved_species,
                    "lineage_id": None,
                    "protein": protein,
                    "position": var["position"],
                    "ref_aa": var["ref_aa"],
                    "alt_aa": var["alt_aa"],
                    "genotype_description": variant_text,
                    "phenotype_category": term,
                    "phenotype_specific": ", ".join(k for _, k in phenotypes),
                    "evidence_strength": "preliminary",
                    "literature_refs": refs,
                    "source_url": source_url,
                    "notes": notes,
                }
            )

    for motif in motifs:
        motif_protein = resolve_motif_protein(motif, text, ref_proteomes)
        motif_species = species_id
        best_species_dist = None
        for m in re.finditer(re.escape(motif), text, re.IGNORECASE):
            s, d = _closest_species_distance(text, m.start(), m.group(0))
            if d is not None and (best_species_dist is None or d < best_species_dist):
                best_species_dist = d
                motif_species = s or motif_species
        if best_species_dist is None:
            # No ebolavirus species mention near any motif occurrence; skip cross-pathogen hits
            log.info("Skipping motif '%s' in PMID:%s (no ebolavirus species in context)", motif, pmid)
            continue
        for term, _ in phenotypes:
            aid = association_id(
                pathogen_id, motif_species, motif_protein, None, None, None, f"motif:{motif}", term
            )
            if aid in seen:
                continue
            seen.add(aid)
            candidates.append(
                {
                    "association_id": aid,
                    "pathogen_id": pathogen_id,
                    "species_id": motif_species,
                    "lineage_id": None,
                    "protein": motif_protein,
                    "position": None,
                    "ref_aa": None,
                    "alt_aa": None,
                    "genotype_description": f"motif: {motif}",
                    "phenotype_category": term,
                    "phenotype_specific": ", ".join(k for _, k in phenotypes),
                    "evidence_strength": "preliminary",
                    "literature_refs": refs,
                    "source_url": source_url,
                    "notes": f"Extracted via PubTator from PMID:{pmid}. Flagged for review.",
                }
            )

    for clade in clades:
        lineage_id = resolve_lineage(clade, lineage_aliases or {})
        for term, _ in phenotypes:
            aid = association_id(
                pathogen_id, species_id, None, None, None, None, f"clade:{clade}", term
            )
            if aid in seen:
                continue
            seen.add(aid)
            candidates.append(
                {
                    "association_id": aid,
                    "pathogen_id": pathogen_id,
                    "species_id": species_id,
                    "lineage_id": lineage_id,
                    "protein": None,
                    "position": None,
                    "ref_aa": None,
                    "alt_aa": None,
                    "genotype_description": f"clade/strain: {clade}",
                    "phenotype_category": term,
                    "phenotype_specific": ", ".join(k for _, k in phenotypes),
                    "evidence_strength": "preliminary",
                    "literature_refs": refs,
                    "source_url": source_url,
                    "notes": (
                        f"Extracted via PubTator from PMID:{pmid}. Flagged for review."
                        + ("" if lineage_id else " No matching lineage_id in DB; kept at species level.")
                    ),
                }
            )

    return candidates
def save_candidates(candidates: list[dict], dry_run: bool = False) -> tuple[int, int]:
    if not candidates:
        return 0, 0
    if dry_run:
        return 0, len(candidates)

    conn = psycopg2.connect(DB_URL)
    cur = conn.cursor()
    inserted = 0
    skipped = 0
    try:
        for cand in candidates:
            cur.execute("SELECT 1 FROM genotype_phenotype WHERE association_id = %s LIMIT 1", (cand["association_id"],))
            if cur.fetchone():
                skipped += 1
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
                    cand["association_id"], cand["pathogen_id"], cand["species_id"], cand["lineage_id"],
                    cand["protein"], cand["position"], cand["ref_aa"], cand["alt_aa"],
                    cand["genotype_description"], cand["phenotype_category"], cand["phenotype_specific"],
                    "not quantified", cand["evidence_strength"],
                    cand["literature_refs"], cand["source_url"], cand["notes"],
                    True, "Auto-extracted via PubTator; requires human review.",
                    "unverified", "pubtator", date.today(), date.today(),
                ),
            )
            inserted += 1
        conn.commit()
        return inserted, skipped
    except Exception as e:
        conn.rollback()
        raise
    finally:
        cur.close()
        conn.close()


def run(
    species: str = "Zaire ebolavirus",
    species_id: str = "EBOV",
    pathogen_id: str = "ebola",
    max_results: int = 50,
    dry_run: bool = False,
    api_key: Optional[str] = None,
) -> tuple[int, int, int]:
    query = build_pubmed_query(species, GENOTYPE_SEARCH_TERMS, PHENOTYPE_SEARCH_TERMS)
    log.info("PubMed query: %s", query)
    pmids = search_pubmed_pmids(query, max_results=max_results, api_key=api_key)
    log.info("Found %d PMIDs", len(pmids))

    pubs = fetch_pubtator_annotations(pmids)
    log.info("Fetched PubTator annotations for %d publications", len(pubs))

    lineage_aliases = load_lineage_aliases()
    ref_proteomes = load_reference_proteomes(pathogen_id)
    all_candidates = []
    for pub in pubs:
        candidates = build_candidates(pub, species_id, pathogen_id, lineage_aliases, ref_proteomes)
        all_candidates.extend(candidates)

    if dry_run:
        log.info("DRY RUN: would insert %d candidates from %d publications", len(all_candidates), len(pubs))
        return len(all_candidates), 0, 0

    inserted, skipped = save_candidates(all_candidates, dry_run=False)
    log.info("Saved %d candidates, skipped %d duplicates", inserted, skipped)
    return len(all_candidates), inserted, skipped
def main():
    parser = argparse.ArgumentParser(description="Extract genotype-phenotype candidates from PubTator")
    parser.add_argument("--species", default="Zaire ebolavirus")
    parser.add_argument("--species-id", default="EBOV")
    parser.add_argument("--pathogen-id", default="ebola")
    parser.add_argument("--max-results", type=int, default=50)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--api-key", default=None, help="NCBI API key (optional)")
    args = parser.parse_args()

    total, inserted, skipped = run(
        species=args.species,
        species_id=args.species_id,
        pathogen_id=args.pathogen_id,
        max_results=args.max_results,
        dry_run=args.dry_run,
        api_key=args.api_key,
    )
    if args.dry_run:
        print(f"Dry run: {total} candidate associations would be inserted")
    else:
        print(f"Done: {total} candidates, {inserted} inserted, {skipped} skipped")


if __name__ == "__main__":
    main()
