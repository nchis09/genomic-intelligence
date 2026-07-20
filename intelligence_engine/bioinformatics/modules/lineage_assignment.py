#!/usr/bin/env python3
"""
Stage 5: Lineage & clade assignment.

Wraps Nextclade as the validated tool for Ebola clade assignment, QC and
mutation calling. For each sample Nextclade produces:
  - clade / lineage assignment
  - QC metrics (status, missing data, mixed sites, etc.)
  - amino-acid and nucleotide substitutions/deletions/insertions
  - aligned genome
  - a small placement tree

This module runs Nextclade once and exposes the parsed results so downstream
stages (variant calling, phylogenetics, normalisation) can reuse them without
rerunning the tool.
"""

from __future__ import annotations

import csv
import json
import logging
import re
import shutil
import sys
import tempfile
from pathlib import Path
from typing import Any, Optional

PROJECT_ROOT = Path(__file__).resolve().parents[4]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from intelligence_engine.bioinformatics.nextclade_runner import (  # noqa: E402
    check_nextclade,
    ensure_dataset,
)

log = logging.getLogger(__name__)

# Nextclade dataset name per species_id. Keep in sync with nextclade_runner.py.
SPECIES_DATASET = {
    "EBOV": "ebola",
    "SUDV": "ebola-sudan",
    "BDBV": "ebola-bundibugyo",
}


def _ensure_nextclade() -> str:
    nxt = check_nextclade()
    if not nxt:
        raise RuntimeError(
            "Nextclade CLI not found. Install with:\n"
            "  /Users/christianndekezi/anaconda3/bin/conda install -c bioconda nextclade -y"
        )
    return nxt


def _parse_mutation_list(value: str) -> list[dict[str, Any]]:
    """Parse Nextclade comma-separated mutation strings like 'GP:A82V,NP:R123K'.

    Ignores empty/invalid tokens.
    """
    if not value or value.strip() in {"", "nan", "NA"}:
        return []

    results = []
    for token in value.split(","):
        token = token.strip()
        if not token:
            continue
        # Match gene:refAAposaltAA or gene:refposalt (nucleotide)
        m = re.match(r"^([A-Za-z0-9_]+):([A-Z*]?)(\d+)([A-Z*]?)$", token)
        if not m:
            # Keep raw token for diagnostics but do not attempt further parsing
            results.append({"raw": token})
            continue
        gene, ref, pos, alt = m.groups()
        results.append(
            {
                "gene": gene,
                "position": int(pos),
                "ref": ref or None,
                "alt": alt or None,
                "hgvs_p": f"{gene}:{ref}{pos}{alt}" if ref and alt else None,
            }
        )
    return results


def _parse_numeric(value: str, default: Any = None) -> Any:
    if value is None or value == "":
        return default
    try:
        return int(value)
    except ValueError:
        try:
            return float(value)
        except ValueError:
            return value


def run_nextclade_analysis(
    fasta_path: Path,
    species_id: str,
    output_dir: Path,
    nextclade_bin: Optional[str] = None,
) -> dict[str, Any]:
    """Run Nextclade on a FASTA file and parse the outputs into a plain dict.

    Args:
        fasta_path: Path to the sample consensus FASTA.
        species_id: Species identifier (EBOV, SUDV, BDBV, ...).
        output_dir: Directory where Nextclade outputs will be written.
        nextclade_bin: Optional path to the Nextclade executable.

    Returns:
        Dict keyed by sample_id containing parsed Nextclade results.
    """
    fasta_path = Path(fasta_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    nextclade = nextclade_bin or _ensure_nextclade()
    dataset_name = SPECIES_DATASET.get(species_id)
    if not dataset_name:
        raise ValueError(f"No Nextclade dataset configured for species_id={species_id}")

    # Ensure dataset is present (downloads if missing)
    ensure_dataset(dataset_name, nextclade_bin=nextclade)

    # Run Nextclade; the runner writes output-all files into output_dir
    from intelligence_engine.bioinformatics.nextclade_runner import run_nextclade_on_fasta

    tsv_path = run_nextclade_on_fasta(
        input_fasta=fasta_path,
        output_dir=output_dir,
        dataset_name=dataset_name,
        nextclade_bin=nextclade,
    )

    results: dict[str, dict[str, Any]] = {}

    # Parse the TSV summary (one row per sequence)
    with open(tsv_path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            sample_id = row.get("seqName") or row.get("seq_name", "unknown")
            qc_overall = row.get("qc.overallStatus", row.get("qc_overallStatus", ""))

            # Mutation strings are kept in their original Nextclade format so that
            # the variant-calling module (Stage 4) can parse them consistently.
            aa_subs = row.get("aaSubstitutions", "")
            aa_dels = row.get("aaDeletions", "")
            aa_ins = row.get("aaInsertions", "")
            nt_subs = row.get("substitutions", "")
            nt_dels = row.get("deletions", "")
            nt_ins = row.get("insertions", "")

            results[sample_id] = {
                "sample_id": sample_id,
                "species_id": species_id,
                "clade": row.get("clade", ""),
                "lineage": row.get("lineage", ""),
                "qc": {
                    "overall_status": qc_overall,
                    "alignment_score": _parse_numeric(row.get("alignmentScore")),
                    "total_substitutions": _parse_numeric(row.get("totalSubstitutions"), 0),
                    "total_deletions": _parse_numeric(row.get("totalDeletions"), 0),
                    "total_insertions": _parse_numeric(row.get("totalInsertions"), 0),
                    "total_missing": _parse_numeric(row.get("totalMissing"), 0),
                    "total_non_acgt": _parse_numeric(row.get("totalNonACGTNs"), 0),
                    "nearest_node_id": row.get("nearestNodeId", ""),
                    "missing_positions": [int(x) for x in row.get("missing", "").split(",") if x],
                    "qc_missing_data": row.get("qc.missingData.status", ""),
                    "qc_mixed_sites": row.get("qc.mixedSites.status", ""),
                    "qc_private_mutations": row.get("qc.privateMutations.status", ""),
                    "qc_snp_clusters": row.get("qc.snpClusters.status", ""),
                    "qc_frame_shifts": row.get("qc.frameShifts.status", ""),
                    "qc_stop_codons": row.get("qc.stopCodons.status", ""),
                },
                "mutations": {
                    "aa_substitutions": aa_subs,
                    "aa_deletions": aa_dels,
                    "aa_insertions": aa_ins,
                    "nt_substitutions": nt_subs,
                    "nt_deletions": nt_dels,
                    "nt_insertions": nt_ins,
                },
                "nextclade": {
                    "tsv": str(tsv_path),
                    "aligned_fasta": str(output_dir / "nextclade.aligned.fasta"),
                    "tree": str(output_dir / "nextclade.nwk"),
                    "json": str(output_dir / "nextclade.json"),
                    "auspice_json": str(output_dir / "nextclade.auspice.json"),
                    "output_dir": str(output_dir),
                },
            }

    # Save a stable JSON copy of the parsed results for traceability
    parsed_path = output_dir / "nextclade_parsed.json"
    parsed_path.write_text(json.dumps(results, indent=2, default=str))
    log.info(f"Parsed Nextclade results written to {parsed_path}")

    return results


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Stage 5: Lineage & clade assignment via Nextclade")
    parser.add_argument("--fasta", required=True, help="Input consensus FASTA")
    parser.add_argument("--species-id", required=True, help="Species identifier (e.g. EBOV, SUDV)")
    parser.add_argument("--output-dir", required=True, help="Directory for Nextclade outputs")
    parser.add_argument("--nextclade-bin", default=None, help="Path to nextclade executable")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    results = run_nextclade_analysis(
        fasta_path=Path(args.fasta),
        species_id=args.species_id,
        output_dir=Path(args.output_dir),
        nextclade_bin=args.nextclade_bin,
    )

    for sample_id, res in results.items():
        print(f"{sample_id}: clade={res['clade']}, qc={res['qc']['overall_status']}")


if __name__ == "__main__":
    main()
