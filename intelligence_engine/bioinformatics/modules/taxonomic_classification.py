#!/usr/bin/env python3
"""
Stage 1: Taxonomic classification.

Determines the pathogen species from a consensus FASTA using Kraken2 (local
k-mer classifier) and/or BLASTn (sequence similarity). The module follows the
same pattern as the other bioinformatics modules: a `run_taxonomic_classification()`
function returns a dict of per-sample classifications, and a CLI `main()` is
provided for standalone use.
"""

from __future__ import annotations

import argparse
import json
import logging
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import yaml
from Bio import SeqIO
from Bio.Blast import NCBIXML, NCBIWWW

log = logging.getLogger(__name__)

DEFAULT_BLAST_EVALUE = 0.001
DEFAULT_BLAST_HITLIST_SIZE = 10
DEFAULT_BLAST_RETRIES = 2
DEFAULT_BLAST_RETRY_DELAY = 5
DEFAULT_BLAST_INTER_QUERY_DELAY = 10
SUPPORTED_EBOLA_SPECIES = {"EBOV", "SUDV", "BDBV", "RESTV", "TAFV", "BOMV"}


def parse_fasta(fasta_path: Path) -> Dict[str, str]:
    """Parse FASTA file into {header: sequence}."""
    return {record.id: str(record.seq) for record in SeqIO.parse(fasta_path, "fasta")}


def load_taxonomy(taxonomy_path: Optional[Path] = None) -> Dict[str, Any]:
    """Load taxonomy reference library.

    If no path is provided, load the bundled Ebola taxonomy YAML.
    """
    if taxonomy_path is None:
        taxonomy_path = Path(__file__).with_suffix("").parent / "schemas" / "ebola_taxonomy.yaml"
    with open(taxonomy_path) as f:
        return yaml.safe_load(f)


def _species_id_from_name(species_name: str) -> str:
    """Map common/scientific Ebola species name to species_id."""
    mapping = {
        "Zaire ebolavirus": "EBOV",
        "Sudan ebolavirus": "SUDV",
        "Bundibugyo ebolavirus": "BDBV",
        "Reston ebolavirus": "RESTV",
        "Tai Forest ebolavirus": "TAFV",
        "Bombali ebolavirus": "BOMV",
    }
    return mapping.get(species_name, "")


def run_kraken2(fasta_path: Path, kraken_db: Optional[Path] = None) -> Optional[Dict[str, str]]:
    """Run Kraken2 classification if a local database is available."""
    if kraken_db is None or not Path(kraken_db).exists():
        log.info("Kraken2 database not found, skipping")
        return None

    log.info(f"Running Kraken2 with database: {kraken_db}")
    try:
        result = subprocess.run(
            ["kraken2", "--db", str(kraken_db), "--output", "-", "--report", "-", str(fasta_path)],
            capture_output=True,
            text=True,
            check=True,
        )
        classifications: Dict[str, str] = {}
        for line in result.stdout.strip().split("\n"):
            if line.strip():
                parts = line.split("\t")
                if len(parts) >= 3:
                    sample_id = parts[1]
                    classification = parts[2]
                    classifications[sample_id] = classification
        log.info(f"Kraken2 classification complete ({len(classifications)} samples)")
        return classifications
    except FileNotFoundError:
        log.warning("Kraken2 not found. Install with: conda install -c bioconda kraken2")
    except subprocess.CalledProcessError as e:
        log.warning(f"Kraken2 failed: {e.stderr}")
    return None


def run_blastn_online(
    sequence: str,
    sample_id: str,
    output_dir: Optional[Path] = None,
    retries: int = DEFAULT_BLAST_RETRIES,
    retry_delay: int = DEFAULT_BLAST_RETRY_DELAY,
) -> List[Dict[str, Any]]:
    """Run a remote BLASTn search against NCBI nt for one sequence.

    Retries on failure (NCBI can throttle or time out consecutive queries)
    and caches the raw XML response when *output_dir* is provided so failures
    can be inspected.
    """
    last_exc: Optional[Exception] = None
    for attempt in range(1, retries + 1):
        log.info(f"Running NCBI online BLAST for {sample_id} (attempt {attempt}/{retries})")
        try:
            blast_result = NCBIWWW.qblast(
                program="blastn",
                database="nt",
                sequence=sequence,
                expect=DEFAULT_BLAST_EVALUE,
                hitlist_size=DEFAULT_BLAST_HITLIST_SIZE,
                format_type="XML",
            )
            xml_bytes = blast_result.read() if hasattr(blast_result, "read") else blast_result
            if isinstance(xml_bytes, str):
                xml_bytes = xml_bytes.encode("utf-8")

            if output_dir is not None:
                output_dir.mkdir(parents=True, exist_ok=True)
                xml_path = output_dir / f"{sample_id}.blast.xml"
                xml_path.write_bytes(xml_bytes)
                log.info(f"Saved raw BLAST XML for {sample_id}: {xml_path}")

            import io
            blast_records = NCBIXML.parse(io.StringIO(xml_bytes.decode("utf-8")))
            hits: List[Dict[str, Any]] = []
            for record in blast_records:
                query_length = record.query_length
                for alignment in record.alignments:
                    for hsp in alignment.hsps:
                        coverage = (hsp.align_length / query_length * 100) if query_length else 0
                        mismatches = hsp.align_length - hsp.identities
                        hits.append({
                            "hit_id": alignment.hit_id,
                            "hit_def": alignment.hit_def,
                            "accession": alignment.accession,
                            "hit_length": alignment.length,
                            "percent_identity": round((hsp.identities / hsp.align_length) * 100, 2),
                            "alignment_length": hsp.align_length,
                            "mismatches": mismatches,
                            "gap_opens": getattr(hsp, "gaps", 0),
                            "query_start": hsp.query_start,
                            "query_end": hsp.query_end,
                            "hit_start": hsp.sbjct_start,
                            "hit_end": hsp.sbjct_end,
                            "evalue": hsp.expect,
                            "bit_score": hsp.bits,
                            "coverage": round(coverage, 2),
                            "query_length": query_length,
                            "positives": getattr(hsp, "positives", None),
                            "gaps": getattr(hsp, "gaps", 0),
                            "strand": getattr(hsp, "strand", None),
                        })
                        break  # best HSP per alignment
                    if len(hits) >= DEFAULT_BLAST_HITLIST_SIZE:
                        break
            log.info(f"Online BLAST complete for {sample_id} ({len(hits)} hits)")
            return hits
        except Exception as exc:
            last_exc = exc
            log.warning(f"Online BLAST attempt {attempt} failed for {sample_id}: {exc}")
            if attempt < retries:
                sleep_seconds = retry_delay * attempt
                log.info(f"Retrying NCBI BLAST for {sample_id} in {sleep_seconds}s...")
                time.sleep(sleep_seconds)

    log.error(f"Online BLAST failed for {sample_id} after {retries} attempts: {last_exc}")
    return []


def determine_ebolavirus_species(blast_hits: List[Dict[str, Any]], taxonomy: Dict[str, Any]) -> Tuple[str, str, str]:
    """Return (species_name, confidence, method) from BLAST hit text."""
    if not blast_hits:
        return "Unknown", "LOW", "no_hits"

    best_hit = blast_hits[0]
    hit_def = best_hit.get("hit_def", "")
    hit_id = best_hit.get("hit_id", "")
    percent_identity = best_hit.get("percent_identity", 0)

    species_mapping = {
        "Zaire ebolavirus": ["EBOV", "Zaire ebolavirus", "Zaire", "Mayinga", "Ebolavirus"],
        "Sudan ebolavirus": ["SUDV", "Sudan ebolavirus", "Sudan", "Boniface"],
        "Bundibugyo ebolavirus": ["BDBV", "Bundibugyo ebolavirus", "Bundibugyo"],
        "Taï Forest ebolavirus": ["TAFV", "Taï Forest ebolavirus", "Tai Forest", "Ivory Coast"],
        "Reston ebolavirus": ["RESTV", "Reston ebolavirus", "Reston"],
        "Bombali ebolavirus": ["BOMV", "Bombali ebolavirus", "Bombali"],
    }

    detected_species = "Unknown"
    search_text = f"{hit_def} {hit_id}".lower()
    for species, keywords in species_mapping.items():
        for keyword in keywords:
            if keyword.lower() in search_text:
                detected_species = species
                break
        if detected_species != "Unknown":
            break

    if percent_identity >= 98:
        confidence = "HIGH"
    elif percent_identity >= 95:
        confidence = "MODERATE"
    else:
        confidence = "LOW"

    return detected_species, confidence, "blast_ncbi_online"


def classify_sample(
    sample_id: str,
    kraken_result: Optional[str] = None,
    blast_hits: Optional[List[Dict[str, Any]]] = None,
    taxonomy: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Classify a single sample from Kraken2 and/or BLAST results."""
    result: Dict[str, Any] = {
        "sample_id": sample_id,
        "species": "Unknown",
        "species_id": "",
        "pathogen_id": "ebola",
        "pathogen_family": "Filoviridae",
        "pathogen_genus": "Orthoebolavirus",
        "confidence": "LOW",
        "method": "none",
        "kraken2_species": None,
        "kraken2_confidence": None,
        "blast_species": None,
        "blast_confidence": None,
        "blast_hits": [],
        "agreement": False,
    }

    if kraken_result:
        result["kraken2_species"] = kraken_result
        result["kraken2_confidence"] = "HIGH"

    if blast_hits:
        result["blast_hits"] = blast_hits
        blast_species, blast_conf, method = determine_ebolavirus_species(blast_hits, taxonomy)
        result["blast_species"] = blast_species
        result["blast_confidence"] = blast_conf
        result["confidence"] = blast_conf
        result["method"] = method
        result["species"] = blast_species
        result["species_id"] = _species_id_from_name(blast_species)

        if blast_hits:
            best_hit = blast_hits[0]
            result["best_hit"] = {
                "accession": best_hit.get("accession"),
                "percent_identity": best_hit.get("percent_identity"),
                "alignment_length": best_hit.get("alignment_length"),
                "evalue": best_hit.get("evalue"),
                "bit_score": best_hit.get("bit_score"),
                "coverage": best_hit.get("coverage"),
                "mismatches": best_hit.get("mismatches"),
                "gaps": best_hit.get("gaps"),
            }

    if kraken_result and result["blast_species"]:
        if kraken_result == result["blast_species"]:
            result["species"] = result["blast_species"]
            result["confidence"] = result["blast_confidence"]
            result["agreement"] = True
            result["method"] = "kraken2_blast_agreement"
        else:
            result["species"] = f"Conflict: Kraken2={kraken_result}, BLAST={result['blast_species']}"
            result["confidence"] = "LOW"
            result["agreement"] = False
            result["method"] = "kraken2_blast_conflict"
    elif kraken_result:
        result["species"] = kraken_result
        result["confidence"] = "HIGH"
        result["method"] = "kraken2_only"

    # Enrich taxonomy-derived fields if we resolved a clean species
    if result["species"] not in {"Unknown", ""} and not result["species"].startswith("Conflict:"):
        result["species_id"] = _species_id_from_name(result["species"])
        if taxonomy and "species" in taxonomy:
            for entry in taxonomy["species"]:
                if entry.get("scientific_name") == result["species"]:
                    result["pathogen_id"] = entry.get("pathogen_id", "ebola")
                    result["pathogen_family"] = entry.get("family", "Filoviridae")
                    result["pathogen_genus"] = entry.get("genus", "Orthoebolavirus")
                    result["is_supported"] = entry.get("is_supported", False)
                    break

    return result


def run_taxonomic_classification(
    fasta_path: Path,
    output_dir: Path,
    taxonomy_path: Optional[Path] = None,
    kraken_db: Optional[Path] = None,
    skip_kraken: bool = False,
    skip_blast: bool = False,
) -> Dict[str, Dict[str, Any]]:
    """Run Stage 1 taxonomic classification.

    Returns a dict mapping sample_id -> classification record.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    taxonomy = load_taxonomy(taxonomy_path)
    sequences = parse_fasta(fasta_path)
    log.info(f"Loaded {len(sequences)} sequences")

    # Kraken2 (optional)
    kraken_results: Optional[Dict[str, str]] = None
    if not skip_kraken:
        kraken_results = run_kraken2(fasta_path, kraken_db)

    # BLAST (optional)
    blast_results: Dict[str, List[Dict[str, Any]]] = {}
    blast_cache_dir = output_dir / "blast_xml_cache"
    if not skip_blast:
        for idx, (sample_id, sequence) in enumerate(sequences.items()):
            if idx > 0:
                # NCBI can throttle rapid sequential full-genome queries.
                log.info(
                    f"Waiting {DEFAULT_BLAST_INTER_QUERY_DELAY}s before BLAST query for {sample_id}..."
                )
                time.sleep(DEFAULT_BLAST_INTER_QUERY_DELAY)
            hits = run_blastn_online(sequence, sample_id, output_dir=blast_cache_dir)
            if hits:
                blast_results[sample_id] = hits

    # Classify each sample
    classifications: Dict[str, Dict[str, Any]] = {}
    for sample_id in sequences:
        kraken_result = kraken_results.get(sample_id) if kraken_results else None
        blast_hits = blast_results.get(sample_id, [])
        classifications[sample_id] = classify_sample(
            sample_id,
            kraken_result=kraken_result,
            blast_hits=blast_hits,
            taxonomy=taxonomy,
        )
        log.info(
            f"{sample_id}: species={classifications[sample_id]['species']}, "
            f"confidence={classifications[sample_id]['confidence']}, "
            f"method={classifications[sample_id]['method']}"
        )

    # Write outputs
    species_output = output_dir / "species_id.json"
    species_output.write_text(json.dumps(classifications, indent=2))

    report_output = output_dir / "classification_report.txt"
    lines = ["=" * 70, "PATHOGEN CLASSIFICATION REPORT", "=" * 70, ""]
    for sample_id, result in classifications.items():
        lines.extend([
            f"Sample: {sample_id}",
            f"  Species: {result['species']}",
            f"  Confidence: {result['confidence']}",
            f"  Method: {result['method']}",
        ])
        if result["kraken2_species"]:
            lines.append(f"  Kraken2: {result['kraken2_species']}")
        if result["blast_species"]:
            lines.append(f"  BLAST: {result['blast_species']} ({result['blast_confidence']})")
        lines.append("")
    report_output.write_text("\n".join(lines))

    log.info(f"Stage 1 complete: {species_output}")
    return classifications


def main():
    parser = argparse.ArgumentParser(description="Stage 1: Pathogen Classification")
    parser.add_argument("--fasta", required=True, help="Path to FASTA file")
    parser.add_argument("--taxonomy", help="Path to taxonomy.yaml (defaults to bundled Ebola taxonomy)")
    parser.add_argument("--kraken-db", help="Path to Kraken2 database (optional)")
    parser.add_argument("--output", required=True, help="Output directory")
    parser.add_argument("--skip-kraken", action="store_true", help="Skip Kraken2 classification")
    parser.add_argument("--skip-blast", action="store_true", help="Skip BLAST classification")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    run_taxonomic_classification(
        fasta_path=Path(args.fasta),
        output_dir=Path(args.output),
        taxonomy_path=Path(args.taxonomy) if args.taxonomy else None,
        kraken_db=Path(args.kraken_db) if args.kraken_db else None,
        skip_kraken=args.skip_kraken,
        skip_blast=args.skip_blast,
    )


if __name__ == "__main__":
    main()
