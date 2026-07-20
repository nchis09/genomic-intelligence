#!/usr/bin/env python3
"""
Ebola Variant Calling Pipeline
==============================

Streaming pipeline that:
  1. Queries NCBI for all available Ebola genomes per species
  2. For each genome (in batches, sequences NOT stored permanently):
     a. Fetches the GenBank record from NCBI
     b. Extracts CDS features and translates to protein sequences
     c. Keeps only Ebola proteins (NP, VP35, VP40, GP, VP30, VP24, L)
     d. Aligns each protein to the matching reference protein (same gene, same species)
     e. Calls amino acid variants (substitutions, insertions, deletions, stops)
     f. Stores genome metadata + variant calls in the database
     g. Discards the sequence (no sequence storage)

Only variant calls and genome metadata are persisted — not sequences.

Run with anaconda python:
    /Users/christianndekezi/anaconda3/bin/python3 database/ebola/protein_variants/call_variants.py

Options:
    --species EBOV          Only process one species
    --batch-size 50         Number of genomes to fetch per NCBI batch
    --dry-run               Don't write to DB, just report counts
    --skip-existing         Skip genomes already in genome_metadata
"""

import argparse
import logging
import re
import sys
import time
from datetime import date
from pathlib import Path
from typing import Optional

from Bio import Entrez, SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
import psycopg2
import psycopg2.extras

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))
from config import DB_URL

Entrez.email = "pgirl_pipeline@local"
Entrez.tool = "pgirl_variant_pipeline"

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s", datefmt="%H:%M:%S")
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

EBOLA_SPECIES = [
    {"species_id": "EBOV", "species_name": "Zaire ebolavirus", "taxid": 186536,  "accession": "NC_002549.1"},
    {"species_id": "SUDV", "species_name": "Sudan ebolavirus", "taxid": 186538,  "accession": "NC_006432.1"},
    {"species_id": "BDBV", "species_name": "Bundibugyo ebolavirus", "taxid": 186545, "accession": "NC_014373.1"},
    {"species_id": "RESTV", "species_name": "Reston ebolavirus", "taxid": 186539,  "accession": "NC_004161.1"},
    {"species_id": "TAFV", "species_name": "Tai Forest ebolavirus", "taxid": 186544, "accession": "NC_014372.1"},
    {"species_id": "BOMV", "species_name": "Bombali ebolavirus", "taxid": 186549,  "accession": "NC_039345.1"},
]

# Only keep these Ebola proteins — ignore host or unrelated CDS
EBOLA_GENES = {"NP", "VP35", "VP40", "GP", "VP30", "VP24", "L"}

# Map common /product values to gene names when /gene qualifier is missing
PRODUCT_TO_GENE = {
    "nucleoprotein": "NP",
    "polymerase complex protein": "VP35",
    "matrix protein": "VP40",
    "spike glycoprotein": "GP",
    "sGP": "GP",
    "secreted glycoprotein": "GP",
    "minor nucleoprotein": "VP30",
    "membrane-associated protein": "VP24",
    "rna-dependent rna polymerase": "L",
    "polymerase": "L",
}

# Amino acid 3-letter to 1-letter mapping for HGVS notation
AA3_TO_AA1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Sec": "U", "Pyl": "O", "Ter": "*",
}
AA1_TO_AA3 = {v: k for k, v in AA3_TO_AA1.items()}

# Month name → number mapping for NCBI date parsing (DD-Mon-YYYY format)
MONTHS = {
    "jan": 1, "feb": 2, "mar": 3, "apr": 4, "may": 5, "jun": 6,
    "jul": 7, "aug": 8, "sep": 9, "oct": 10, "nov": 11, "dec": 12,
}

# ISO 3166 alpha-3 country codes → country names (Ebola-relevant subset)
COUNTRY_CODES = {
    "COD": "Democratic Republic of the Congo", "COG": "Republic of the Congo",
    "GAB": "Gabon", "SLE": "Sierra Leone", "LBR": "Liberia",
    "GIN": "Guinea", "MLI": "Mali", "SEN": "Senegal",
    "UGA": "Uganda", "SDN": "Sudan", "SSD": "South Sudan",
    "CIV": "Cote d'Ivoire", "CMR": "Cameroon", "CAF": "Central African Republic",
    "TZA": "Tanzania", "KEN": "Kenya", "NGA": "Nigeria",
    "USA": "United States", "PHL": "Philippines", "CHN": "China",
    "AGO": "Angola", "ZMB": "Zambia", "ZWE": "Zimbabwe",
    "BWA": "Botswana", "NAM": "Namibia", "ZAF": "South Africa",
}

# Rate limiting for NCBI E-utilities (3 req/s without API key)
NCBI_SLEEP = 0.35


# ---------------------------------------------------------------------------
# Step 1: Query NCBI for all genome accessions per species
# ---------------------------------------------------------------------------

def search_genomes(species_name: str, taxid: int, batch_size: int = 500) -> list[str]:
    """Search NCBI nucleotide DB for all Ebola genomes for this species.

    Uses the species scientific name (e.g. 'Zaire ebolavirus') rather than raw
    taxid to avoid taxonomy contamination from other species. Searches in two
    passes:
      1. Full/near-full genomes (>15kb) — extract ORFs using reference coordinates
      2. Partial sequences (500bp–15kb) — use NCBI CDS annotations when available
    Excludes the canonical RefSeq reference (loaded separately) by accession.
    """
    all_uids = []
    ref_acc = next(s["accession"] for s in EBOLA_SPECIES if s["taxid"] == taxid)

    # Pass 1: Full genomes >=15kb
    full_term = f'{species_name}[Organism] AND biomol_genomic[PROP] AND 15000:100000[SLEN] NOT {ref_acc}[Accession]'
    full_uids = _esearch_all(full_term, batch_size)
    log.info(f"Found {len(full_uids)} full genomes (>=15kb) for {species_name}")
    all_uids.extend(full_uids)

    # Pass 2: Partial sequences 500bp–15kb (may contain individual gene sequences)
    partial_term = f'{species_name}[Organism] AND biomol_genomic[PROP] AND 500:14999[SLEN] NOT {ref_acc}[Accession]'
    partial_uids = _esearch_all(partial_term, batch_size)
    log.info(f"Found {len(partial_uids)} partial sequences (500bp–15kb) for {species_name}")
    all_uids.extend(partial_uids)

    # Deduplicate (a UID could appear in both searches if at boundary)
    seen = set()
    unique = []
    for uid in all_uids:
        if uid not in seen:
            seen.add(uid)
            unique.append(uid)

    return unique


def _esearch_all(search_term: str, batch_size: int = 500) -> list[str]:
    """Run an NCBI esearch and paginate through all results."""
    accessions = []
    handle = Entrez.esearch(
        db="nucleotide",
        term=search_term,
        retmax=batch_size,
        usehistory="y",
    )
    results = Entrez.read(handle)
    handle.close()
    time.sleep(NCBI_SLEEP)

    count = int(results["Count"])
    webenv = results["WebEnv"]
    query_key = results["QueryKey"]
    accessions.extend(results["IdList"])

    fetched = len(results["IdList"])
    while fetched < count:
        handle = Entrez.esearch(
            db="nucleotide",
            term=search_term,
            retmax=batch_size,
            retstart=fetched,
            usehistory="y",
            WebEnv=webenv,
            query_key=query_key,
        )
        batch = Entrez.read(handle)
        handle.close()
        time.sleep(NCBI_SLEEP)
        accessions.extend(batch["IdList"])
        fetched += len(batch["IdList"])
        if not batch["IdList"]:
            break

    return accessions


def fetch_genome_batch(uids: list[str]) -> list:
    """Fetch GenBank records for a batch of NCBI UIDs."""
    records = []
    # Fetch in sub-batches of 20 to keep response sizes manageable
    for i in range(0, len(uids), 20):
        batch = uids[i:i+20]
        handle = Entrez.efetch(
            db="nucleotide",
            id=",".join(batch),
            rettype="gb",
            retmode="text",
        )
        for rec in SeqIO.parse(handle, "genbank"):
            records.append(rec)
        handle.close()
        time.sleep(NCBI_SLEEP)
    return records


# ---------------------------------------------------------------------------
# Step 2: Extract ORFs / CDS from a genome
# ---------------------------------------------------------------------------

def _gene_from_feature(feature) -> Optional[str]:
    """Resolve Ebola gene name from /gene or /product qualifiers."""
    gene = feature.qualifiers.get("gene", [""])[0]
    if gene in EBOLA_GENES:
        return gene
    product = feature.qualifiers.get("product", [""])[0].strip().lower()
    if product in PRODUCT_TO_GENE:
        return PRODUCT_TO_GENE[product]
    return None


def extract_cds(record) -> list[dict]:
    """Extract CDS features from a GenBank record, keeping only Ebola genes.

    This is used for partial sequences where NCBI has pre-annotated CDS.
    For full genomes without CDS annotations, use extract_orfs_from_reference().
    """
    proteins = []
    seen_genes = set()
    for feature in record.features:
        if feature.type != "CDS":
            continue
        gene = _gene_from_feature(feature)
        if gene is None:
            continue
        translation = feature.qualifiers.get("translation", [""])[0]
        if not translation:
            continue
        # Prefer the first/main CDS per gene (e.g., spike glycoprotein over sGP)
        if gene in seen_genes:
            continue
        seen_genes.add(gene)
        protein_id = feature.qualifiers.get("protein_id", [""])[0]
        product = feature.qualifiers.get("product", [""])[0]
        location = feature.location
        proteins.append({
            "gene": gene,
            "protein_id": protein_id,
            "protein_name": product,
            "protein_sequence": translation,
            "genome_start": int(location.start) + 1,
            "genome_end": int(location.end),
            "strand": int(location.strand) if location.strand else 1,
        })
    return proteins


def extract_orfs_from_reference(record, ref_proteins: dict[str, dict]) -> list[dict]:
    """Extract Ebola proteins from a genome without NCBI CDS annotations.

    Translates the genome in all 6 reading frames, finds continuous ORFs, and
    matches them to the reference proteins by identity. This automatically
    handles both forward and reverse-complement orientations (some submitters
    deposit the negative-sense strand).
    """
    from Bio.Seq import Seq

    genome_seq = record.seq
    genome_len = len(genome_seq)
    if genome_len < 1000:
        return []

    # Build reference info: min length, identity matrix
    ref_info_list = []
    for gene, info in ref_proteins.items():
        ref_seq = info["sequence"].rstrip("*")
        if len(ref_seq) < 50:
            continue
        ref_info_list.append((gene, ref_seq, len(ref_seq)))

    # Collect ORFs from both strands
    orf_candidates = []  # (gene, protein, start_1based, end_1based, strand, score)

    for strand_name, strand_seq in [("+", genome_seq), ("-", genome_seq.reverse_complement())]:
        strand = 1 if strand_name == "+" else -1
        for frame in range(3):
            translated = str(strand_seq[frame:].translate())
            # Split into ORFs at stop codons
            start_aa = 0
            for stop_aa in range(len(translated)):
                if translated[stop_aa] == "*":
                    orf = translated[start_aa:stop_aa]
                    # Approximate nucleotide coordinates
                    nt_start = frame + start_aa * 3 + 1
                    nt_end = frame + stop_aa * 3
                    if strand == -1:
                        # Reverse map coordinates to original genome
                        nt_start_orig = genome_len - nt_end + 1
                        nt_end_orig = genome_len - nt_start + 1
                        nt_start, nt_end = nt_start_orig, nt_end_orig
                    # Clean trailing X's
                    orf = orf.rstrip("X")
                    if len(orf) >= 150:
                        orf_candidates.append((orf, nt_start, nt_end, strand))
                    start_aa = stop_aa + 1

    # Score each ORF against each reference protein and assign best per gene
    best_per_gene = {}
    for orf, nt_start, nt_end, strand in orf_candidates:
        # Skip if too many ambiguous residues
        if orf.count("X") / len(orf) > 0.3:
            continue
        for gene, ref_seq, ref_len in ref_info_list:
            if len(orf) < ref_len * 0.5:
                continue
            # Global alignment-ish score: count matches at best offset
            best_match = 0
            min_len = min(len(orf), ref_len)
            # Try a few offsets to account for indels/frames
            for offset in range(-30, 31, 10):
                matches = 0
                for i in range(max(0, -offset), min(min_len, len(orf) - offset, ref_len)):
                    if orf[i + offset] == ref_seq[i] and orf[i + offset] != "X":
                        matches += 1
                if matches > best_match:
                    best_match = matches
            score = best_match / ref_len if ref_len > 0 else 0
            if score > 0.25:  # 25% identity threshold
                if gene not in best_per_gene or score > best_per_gene[gene]["score"]:
                    best_per_gene[gene] = {
                        "protein": orf,
                        "start": nt_start,
                        "end": nt_end,
                        "strand": strand,
                        "score": score,
                    }

    proteins = []
    for gene, data in best_per_gene.items():
        proteins.append({
            "gene": gene,
            "protein_id": "",
            "protein_name": "",
            "protein_sequence": data["protein"],
            "genome_start": min(data["start"], data["end"]),
            "genome_end": max(data["start"], data["end"]),
            "strand": data["strand"],
        })

    return proteins


# ---------------------------------------------------------------------------
# Step 3: Extract metadata from GenBank record
# ---------------------------------------------------------------------------

def extract_metadata(record, species: dict) -> dict:
    """Extract collection metadata from GenBank record qualifiers."""
    meta = {
        "genome_accession": record.name,  # e.g. "KC242791"
        "pathogen_id": "ebola",
        "species_id": species["species_id"],
        "reference_accession": species["accession"],
        "ncbi_taxonomy_id": species["taxid"],
        "strain": "",
        "isolate": "",
        "collection_date": None,
        "collection_year": None,
        "collection_country": "",
        "collection_country_code": "",
        "collection_region": "",
        "host": "",
        "isolation_source": "",
        "genome_length": len(record.seq),
        "completeness": "complete",
        "source_db": "NCBI",
        "release_date": None,
    }

    # Parse source feature for metadata
    for feature in record.features:
        if feature.type != "source":
            continue
        q = feature.qualifiers
        meta["strain"] = q.get("strain", [""])[0]
        meta["isolate"] = q.get("isolate", [""])[0]
        meta["host"] = q.get("host", [""])[0]
        meta["isolation_source"] = q.get("isolation_source", [""])[0]
        meta["collection_country"] = q.get("country", [""])[0]
        meta["collection_region"] = q.get("region", [""])[0]

        # Parse collection date
        coll_date = q.get("collection_date", [""])[0]
        if coll_date:
            parsed = parse_date(coll_date)
            if parsed:
                meta["collection_date"] = parsed
                meta["collection_year"] = parsed.year

        # Country code (sometimes in /country qualifier as "Gabon: Booue")
        country_val = q.get("country", [""])[0]
        if ":" in country_val:
            parts = country_val.split(":", 1)
            meta["collection_country"] = parts[0].strip()
            if not meta["collection_region"]:
                meta["collection_region"] = parts[1].strip()

    # Parse release date from record annotations
    if record.annotations.get("date"):
        parsed = parse_date(record.annotations["date"])
        if parsed:
            meta["release_date"] = parsed

    # If country is still missing, try to extract from isolate/strain name
    # NCBI isolate format: EBOV/H.sapiens-tc/COD/1977/Bonduni
    #                       Ebola virus/H.sapiens-wt/SLE/2014/Makona-G3770.2
    if not meta["collection_country"]:
        for name_field in [meta["isolate"], meta["strain"]]:
            if not name_field:
                continue
            # Look for 3-letter country code in the name
            parts = name_field.split("/")
            for part in parts:
                code = part.strip().upper()
                if code in COUNTRY_CODES:
                    meta["collection_country"] = COUNTRY_CODES[code]
                    meta["collection_country_code"] = code
                    break
            if meta["collection_country"]:
                break

    # If date is still missing but we have it in isolate name, extract year
    if not meta["collection_date"]:
        for name_field in [meta["isolate"], meta["strain"]]:
            if not name_field:
                continue
            # Look for 4-digit year in the name
            match = re.search(r"\b(19\d{2}|20\d{2})\b", name_field)
            if match:
                year = int(match.group(1))
                meta["collection_year"] = year
                meta["collection_date"] = date(year, 1, 1)
                break

    # Check completeness
    desc = (record.description or "").lower()
    if "partial" in desc:
        meta["completeness"] = "partial"
    elif "nearly complete" in desc or "near complete" in desc:
        meta["completeness"] = "nearly complete"

    # Truncate string fields to prevent DB truncation errors
    MAX_LEN = {
        "strain": 255, "isolate": 255, "collection_country": 255,
        "collection_region": 255, "host": 255, "isolation_source": 255,
        "completeness": 30,
    }
    for k, maxlen in MAX_LEN.items():
        if meta.get(k) and len(str(meta[k])) > maxlen:
            meta[k] = str(meta[k])[:maxlen]

    # Remove empty strings → None for DB
    for k, v in meta.items():
        if v == "":
            meta[k] = None

    return meta


def parse_date(s: str) -> Optional[date]:
    """Parse NCBI date strings in various formats:
      - YYYY-MM-DD  (e.g. 2014-06-14)
      - DD-Mon-YYYY (e.g. 14-Jun-2014)
      - YYYY        (e.g. 1977)
    """
    s = s.strip()
    if not s:
        return None
    parts = s.split("-")
    try:
        if len(parts) == 3:
            # Could be YYYY-MM-DD or DD-Mon-YYYY
            if parts[0].isdigit() and len(parts[0]) == 4:
                # YYYY-MM-DD
                return date(int(parts[0]), int(parts[1]), int(parts[2]))
            else:
                # DD-Mon-YYYY (e.g. 14-Jun-2014)
                day = int(parts[0])
                month = MONTHS.get(parts[1].lower())
                year = int(parts[2])
                if month:
                    return date(year, month, day)
        elif len(parts) == 1 and parts[0].isdigit():
            return date(int(parts[0]), 1, 1)
    except (ValueError, IndexError):
        pass
    return None


# ---------------------------------------------------------------------------
# Step 4: Align protein to reference and call variants
# ---------------------------------------------------------------------------

aligner = PairwiseAligner()
aligner.mode = "global"
aligner.open_gap_score = -10
aligner.extend_gap_score = -0.5
aligner.substitution_matrix = None  # default BLOSUM62


def load_reference_proteins(cur, species_id: str) -> dict[str, dict]:
    """Load reference proteins for a species from the database."""
    cur.execute(
        "SELECT gene, protein_sequence, reference_accession FROM reference_proteomes WHERE species_id = %s",
        (species_id,),
    )
    refs = {}
    for row in cur.fetchall():
        refs[row[0]] = {
            "sequence": row[1],
            "reference_accession": row[2],
        }
    return refs


def call_variants(ref_seq: str, query_seq: str) -> list[dict]:
    """
    Align query protein to reference protein and call amino acid variants.

    Returns a list of variant dicts with position, ref_aa, alt_aa, variant_type, hgvs_p.
    Skips positions with ambiguous amino acids (X).
    """
    # Clean sequences — remove trailing stop codons for alignment
    ref_clean = ref_seq.rstrip("*")
    query_clean = query_seq.rstrip("*")

    if not query_clean or len(query_clean) < 10:
        return []

    # Skip if query is mostly ambiguous
    non_x = sum(1 for aa in query_clean if aa != "X")
    if non_x < len(query_clean) * 0.5:
        log.warning(f"Skipping protein: >50% ambiguous residues")
        return []

    alignments = aligner.align(ref_clean, query_clean)
    if not alignments:
        log.warning("No alignment produced")
        return []

    best = alignments[0]
    ref_aligned, query_aligned = best[0], best[1]

    variants = []
    ref_pos = 0  # 1-based position in reference protein

    for i in range(len(ref_aligned)):
        ref_aa = ref_aligned[i]
        qry_aa = query_aligned[i]

        if ref_aa == "-":
            # Insertion in query (relative to reference)
            ref_pos += 0  # don't advance ref position
            if qry_aa != "-" and qry_aa != "X":
                variants.append({
                    "position": ref_pos,  # position after which insertion occurs
                    "ref_aa": "-",
                    "alt_aa": qry_aa,
                    "variant_type": "insertion",
                    "hgvs_p": f"p.{ref_pos}_{ref_pos+1}ins{qry_aa}" if ref_pos > 0 else f"p.1ins{qry_aa}",
                })
        elif qry_aa == "-":
            # Deletion in query (relative to reference)
            ref_pos += 1
            if ref_aa != "X":
                variants.append({
                    "position": ref_pos,
                    "ref_aa": ref_aa,
                    "alt_aa": "-",
                    "variant_type": "deletion",
                    "hgvs_p": f"p.{AA1_TO_AA3.get(ref_aa, ref_aa)}{ref_pos}del",
                })
        elif ref_aa == qry_aa:
            ref_pos += 1
            # Match — no variant
        else:
            ref_pos += 1
            # Skip ambiguous positions
            if ref_aa == "X" or qry_aa == "X":
                continue
            # Determine variant type
            vtype = "substitution"
            if qry_aa == "*":
                vtype = "stop_gained"
            elif ref_aa == "*":
                vtype = "stop_lost"

            hgvs = f"p.{AA1_TO_AA3.get(ref_aa, ref_aa)}{ref_pos}{AA1_TO_AA3.get(qry_aa, qry_aa)}"
            variants.append({
                "position": ref_pos,
                "ref_aa": ref_aa,
                "alt_aa": qry_aa,
                "variant_type": vtype,
                "hgvs_p": hgvs,
            })

    return variants


# ---------------------------------------------------------------------------
# Step 5: Store results in database
# ---------------------------------------------------------------------------

def store_genome_metadata(cur, meta: dict, dry_run: bool = False):
    """Insert or update genome metadata."""
    if dry_run:
        return
    cur.execute(
        """
        INSERT INTO genome_metadata (
            genome_accession, pathogen_id, species_id, reference_accession,
            ncbi_taxonomy_id, strain, isolate,
            collection_date, collection_year, collection_country,
            collection_country_code, collection_region, host, isolation_source,
            genome_length, completeness, source_db, release_date, last_updated
        ) VALUES (
            %(genome_accession)s, %(pathogen_id)s, %(species_id)s, %(reference_accession)s,
            %(ncbi_taxonomy_id)s, %(strain)s, %(isolate)s,
            %(collection_date)s, %(collection_year)s, %(collection_country)s,
            %(collection_country_code)s, %(collection_region)s, %(host)s, %(isolation_source)s,
            %(genome_length)s, %(completeness)s, %(source_db)s, %(release_date)s, CURRENT_DATE
        )
        ON CONFLICT (genome_accession) DO UPDATE SET
            collection_date = EXCLUDED.collection_date,
            collection_country = EXCLUDED.collection_country,
            host = EXCLUDED.host,
            genome_length = EXCLUDED.genome_length,
            last_updated = CURRENT_DATE
        """,
        meta,
    )


def upsert_aggregated_variants(cur, pathogen_id: str, species_id: str,
                                  reference_accession: str, gene: str,
                                  variants: list[dict], genome_meta: dict,
                                  dry_run: bool = False):
    """Upsert aggregated variant counts for one genome's variants.

    For each unique variant (gene+position+ref_aa+alt_aa+variant_type), either:
    - INSERT a new row with genome_count=1 and this genome's metadata, or
    - UPDATE the existing row: increment count, merge countries/dates/lineages
    """
    if dry_run or not variants:
        return

    coll_date = genome_meta.get("collection_date")
    coll_year = genome_meta.get("collection_year")
    coll_country = genome_meta.get("collection_country")
    coll_country_code = genome_meta.get("collection_country_code")
    lineage_id = genome_meta.get("lineage_id")

    for v in variants:
        is_stop = v["variant_type"] in ("stop_gained", "stop_lost")
        cur.execute(
            """
            INSERT INTO protein_variants (
                pathogen_id, species_id, gene, reference_accession,
                position, ref_aa, alt_aa, variant_type, hgvs_p, is_stop,
                genome_count, first_seen_date, last_seen_date,
                first_seen_year, last_seen_year,
                countries_seen, country_codes, lineage_ids, last_updated
            ) VALUES (
                %s, %s, %s, %s, %s, %s, %s, %s, %s, %s,
                1, %s, %s, %s, %s,
                %s, %s, %s, CURRENT_DATE
            )
            ON CONFLICT (pathogen_id, species_id, gene, position, ref_aa, alt_aa, variant_type)
            DO UPDATE SET
                genome_count = protein_variants.genome_count + 1,
                first_seen_date = LEAST(protein_variants.first_seen_date, EXCLUDED.first_seen_date),
                last_seen_date = GREATEST(protein_variants.last_seen_date, EXCLUDED.last_seen_date),
                first_seen_year = LEAST(protein_variants.first_seen_year, EXCLUDED.first_seen_year),
                last_seen_year = GREATEST(protein_variants.last_seen_year, EXCLUDED.last_seen_year),
                countries_seen =
                    CASE
                        WHEN EXCLUDED.countries_seen[1] IS NULL THEN protein_variants.countries_seen
                        WHEN protein_variants.countries_seen = '{}' THEN EXCLUDED.countries_seen
                        ELSE array(
                            SELECT DISTINCT unnest(
                                protein_variants.countries_seen || EXCLUDED.countries_seen
                            ) ORDER BY 1
                        )
                    END,
                country_codes =
                    CASE
                        WHEN EXCLUDED.country_codes[1] IS NULL THEN protein_variants.country_codes
                        WHEN protein_variants.country_codes = '{}' THEN EXCLUDED.country_codes
                        ELSE array(
                            SELECT DISTINCT unnest(
                                protein_variants.country_codes || EXCLUDED.country_codes
                            ) ORDER BY 1
                        )
                    END,
                lineage_ids =
                    CASE
                        WHEN EXCLUDED.lineage_ids[1] IS NULL THEN protein_variants.lineage_ids
                        WHEN protein_variants.lineage_ids = '{}' THEN EXCLUDED.lineage_ids
                        ELSE array(
                            SELECT DISTINCT unnest(
                                protein_variants.lineage_ids || EXCLUDED.lineage_ids
                            ) ORDER BY 1
                        )
                    END,
                last_updated = CURRENT_DATE
            """,
            (
                pathogen_id, species_id, gene, reference_accession,
                v["position"], v["ref_aa"], v["alt_aa"],
                v["variant_type"], v["hgvs_p"], is_stop,
                coll_date, coll_date, coll_year, coll_year,
                [coll_country] if coll_country else [],
                [coll_country_code] if coll_country_code else [],
                [lineage_id] if lineage_id else [],
            ),
        )


def genome_exists(cur, accession: str) -> bool:
    """Check if a genome is already in genome_metadata."""
    cur.execute("SELECT 1 FROM genome_metadata WHERE genome_accession = %s", (accession,))
    return cur.fetchone() is not None


# ---------------------------------------------------------------------------
# Step 6: Update species_total_genomes after processing
# ---------------------------------------------------------------------------

def update_species_totals(cur, pathogen_id: str, dry_run: bool = False):
    """Update species_total_genomes on all protein_variants rows for a pathogen."""
    if dry_run:
        return
    cur.execute(
        """
        UPDATE protein_variants pv
        SET species_total_genomes = sub.total
        FROM (
            SELECT species_id, count(*) AS total
            FROM genome_metadata
            WHERE pathogen_id = %s
            GROUP BY species_id
        ) sub
        WHERE pv.pathogen_id = %s AND pv.species_id = sub.species_id
        """,
        (pathogen_id, pathogen_id),
    )


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def run_pipeline(species_filter: Optional[str] = None, batch_size: int = 50,
                 dry_run: bool = False, skip_existing: bool = True,
                 rebuild: bool = False):
    """Run the full variant calling pipeline.

    If rebuild=True, reprocess ALL genomes (fetch from NCBI, call variants,
    aggregate) even if they already exist in genome_metadata. This is used
    when protein_variants has been truncated and needs to be repopulated.
    """
    conn = psycopg2.connect(DB_URL)
    conn.autocommit = False
    cur = conn.cursor()

    if rebuild:
        skip_existing = False
        log.info("REBUILD mode: will reprocess all genomes to repopulate aggregated variants")

    species_list = EBOLA_SPECIES
    if species_filter:
        species_list = [s for s in EBOLA_SPECIES if s["species_id"] == species_filter]
        if not species_list:
            log.error(f"Unknown species: {species_filter}")
            return

    total_genomes = 0
    total_variants = 0

    for sp in species_list:
        log.info(f"=== Processing {sp['species_id']} ({sp['species_name']}) ===")

        # Load reference proteins for this species
        ref_proteins = load_reference_proteins(cur, sp["species_id"])
        if not ref_proteins:
            log.warning(f"No reference proteins found for {sp['species_id']}. Run fetch_reference_proteomes.py first.")
            continue

        log.info(f"Loaded {len(ref_proteins)} reference proteins: {list(ref_proteins.keys())}")

        # Search NCBI for all genomes
        uids = search_genomes(sp["species_name"], sp["taxid"], batch_size)
        log.info(f"Found {len(uids)} genomes to process for {sp['species_id']}")

        # Process in batches
        for batch_start in range(0, len(uids), batch_size):
            batch_uids = uids[batch_start:batch_start + batch_size]
            log.info(f"  Batch {batch_start//batch_size + 1}: fetching {len(batch_uids)} genomes...")

            try:
                records = fetch_genome_batch(batch_uids)
            except Exception as e:
                log.error(f"  Failed to fetch batch: {e}")
                continue

            for rec in records:
                acc = rec.name
                meta_in_db = genome_exists(cur, acc)
                if skip_existing and meta_in_db:
                    log.debug(f"  Skipping {acc} (already in DB)")
                    continue

                # Extract metadata
                meta = extract_metadata(rec, sp)

                # Set completeness based on genome length
                if meta["genome_length"] >= 18000:
                    meta["completeness"] = "complete"
                elif meta["genome_length"] >= 12000:
                    meta["completeness"] = "nearly complete"
                else:
                    meta["completeness"] = "partial"

                # If rebuild and metadata already exists, skip metadata insert
                # but still process variants
                if rebuild and meta_in_db:
                    # Load lineage_id from existing metadata for aggregation
                    cur.execute(
                        "SELECT lineage_id FROM genome_metadata WHERE genome_accession = %s",
                        (acc,),
                    )
                    row = cur.fetchone()
                    if row:
                        meta["lineage_id"] = row[0]

                # Extract proteins: try NCBI CDS annotations first, then ORF extraction
                proteins = extract_cds(rec)
                extraction_method = "CDS"
                if not proteins and meta["genome_length"] >= 1000:
                    # Full or near-full genome without CDS annotations — extract ORFs ourselves
                    proteins = extract_orfs_from_reference(rec, ref_proteins)
                    extraction_method = "ORF"
                if not proteins:
                    log.warning(f"  {acc}: no proteins extracted (CDS or ORF), skipping")
                    continue

                # Store metadata (skip if rebuild and already exists)
                if not (rebuild and meta_in_db):
                    store_genome_metadata(cur, meta, dry_run)

                # For each protein, align to reference and call variants
                genome_variants = 0
                for prot in proteins:
                    gene = prot["gene"]
                    if gene not in ref_proteins:
                        log.debug(f"  {acc}: no reference for gene {gene}, skipping")
                        continue

                    ref_seq = ref_proteins[gene]["sequence"]
                    ref_acc = ref_proteins[gene]["reference_accession"]
                    query_seq = prot["protein_sequence"]

                    variants = call_variants(ref_seq, query_seq)
                    upsert_aggregated_variants(cur, "ebola", sp["species_id"], ref_acc, gene, variants, meta, dry_run)
                    genome_variants += len(variants)

                total_genomes += 1
                total_variants += genome_variants
                log.info(f"  {acc}: {len(proteins)} proteins ({extraction_method}), {genome_variants} variants")

            # Commit after each batch
            if not dry_run:
                conn.commit()
                log.info(f"  Committed batch {batch_start//batch_size + 1}")

    # Update species_total_genomes on all variants
    if not dry_run:
        log.info("Updating species_total_genomes...")
        update_species_totals(cur, "ebola", dry_run)
        conn.commit()

    cur.close()
    conn.close()

    log.info(f"=== Pipeline complete: {total_genomes} genomes, {total_variants} variants ===")
    if dry_run:
        log.info("(dry-run: nothing was written to DB)")


def main():
    parser = argparse.ArgumentParser(description="Ebola variant calling pipeline")
    parser.add_argument("--species", type=str, default=None,
                        help="Only process one species (e.g. EBOV)")
    parser.add_argument("--batch-size", type=int, default=50,
                        help="Number of genomes per NCBI fetch batch")
    parser.add_argument("--dry-run", action="store_true",
                        help="Don't write to DB, just report counts")
    parser.add_argument("--skip-existing", action="store_true", default=True,
                        help="Skip genomes already in genome_metadata")
    parser.add_argument("--rebuild", action="store_true", default=False,
                        help="Reprocess ALL genomes to repopulate aggregated variants (ignores --skip-existing)")
    args = parser.parse_args()

    run_pipeline(
        species_filter=args.species,
        batch_size=args.batch_size,
        dry_run=args.dry_run,
        skip_existing=args.skip_existing,
        rebuild=args.rebuild,
    )


if __name__ == "__main__":
    main()
