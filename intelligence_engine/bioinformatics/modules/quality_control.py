#!/usr/bin/env python3
"""
Stage 0: Input & Quality Control.

Modular implementation with the same architecture as other bioinformatics
modules: a `run_quality_control()` function that takes file paths and returns
a dict of results, plus a CLI `main()` entry point.

Input:
  - CSV/TSV metadata file
  - FASTA file (one or more sequences)

Output:
  - quality_metrics.json
  - validated_metadata.csv
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import yaml

log = logging.getLogger(__name__)

DEFAULT_EXPECTED_LENGTH = 18900
DEFAULT_TOLERANCE = 0.1
MISSING_REGION_MIN_RUN = 50


def parse_fasta(fasta_path: Path) -> Dict[str, str]:
    """Parse FASTA file and return {header: sequence}."""
    sequences: Dict[str, str] = {}
    current_header: Optional[str] = None
    current_seq: List[str] = []

    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_header is not None:
                    sequences[current_header] = "".join(current_seq)
                current_header = line[1:]
                current_seq = []
            else:
                current_seq.append(line)

        if current_header is not None:
            sequences[current_header] = "".join(current_seq)

    return sequences


def parse_metadata(csv_path: Path) -> List[Dict[str, str]]:
    """Parse metadata CSV/TSV file and return list of row dicts."""
    text = csv_path.read_text()
    dialect = csv.Sniffer().sniff(text[:2048], delimiters=",\t")
    return list(csv.DictReader(text.splitlines(), dialect=dialect))


def load_metadata_schema(schema_path: Optional[Path] = None) -> Dict[str, Any]:
    """Load the metadata schema YAML and return the list of required fields."""
    if schema_path is None:
        schema_path = Path(__file__).with_suffix("").parent / "schemas" / "metadata_schema.yaml"
    with open(schema_path) as f:
        return yaml.safe_load(f)


def get_required_fields(schema: Dict[str, Any]) -> List[str]:
    return [f["name"] for f in schema.get("fields", []) if f.get("tier") == "required"]


def validate_sample_id_matching(
    sequences: Dict[str, str], metadata: List[Dict[str, str]]
) -> Tuple[bool, List[str]]:
    errors = []
    metadata_ids = {row["sample_id"] for row in metadata if "sample_id" in row}
    fasta_ids = set(sequences.keys())

    missing_in_fasta = metadata_ids - fasta_ids
    if missing_in_fasta:
        errors.append(f"Sample IDs in metadata but not in FASTA: {missing_in_fasta}")

    extra_in_fasta = fasta_ids - metadata_ids
    if extra_in_fasta:
        errors.append(f"Sample IDs in FASTA but not in metadata: {extra_in_fasta}")

    return len(errors) == 0, errors


def validate_metadata_schema(
    metadata: List[Dict[str, str]], required_fields: List[str]
) -> Tuple[bool, List[str]]:
    errors = []
    if not metadata:
        errors.append("Metadata CSV is empty")
        return False, errors

    header = set(metadata[0].keys())
    for field in required_fields:
        if field not in header:
            errors.append(f"Required field '{field}' missing from metadata")

    for row in metadata:
        sample_id = row.get("sample_id", "UNKNOWN")
        for field in required_fields:
            if field in row and not row[field].strip():
                errors.append(f"Sample {sample_id}: required field '{field}' is empty")

    return len(errors) == 0, errors


def run_fastqc(fastq_path: Path, output_dir: Path) -> Optional[str]:
    """Run FastQC on FASTQ reads, if available."""
    log.info(f"Running FastQC on {fastq_path}")
    try:
        subprocess.run(
            ["fastqc", "-o", str(output_dir), str(fastq_path)],
            capture_output=True,
            text=True,
            check=True,
        )
        report_path = output_dir / f"{fastq_path.stem}_fastqc.html"
        log.info(f"FastQC report: {report_path}")
        return str(report_path)
    except FileNotFoundError:
        log.warning("FastQC not found. Install with: conda install -c bioconda fastqc")
        return None
    except subprocess.CalledProcessError as e:
        log.warning(f"FastQC failed: {e.stderr}")
        return None


def run_seqkit_stats(fasta_path: Path) -> Optional[Dict[str, Dict[str, Any]]]:
    """Run seqkit fx2tab and return per-sequence stats."""
    log.info(f"Running seqkit stats on {fasta_path}")
    try:
        result = subprocess.run(
            ["seqkit", "fx2tab", "-n", "-l", "-g", str(fasta_path)],
            capture_output=True,
            text=True,
            check=True,
        )
        stats: Dict[str, Dict[str, Any]] = {}
        for line in result.stdout.strip().split("\n"):
            parts = line.split("\t")
            if len(parts) >= 3:
                seq_name = parts[0]
                seq_length = int(parts[1])
                gc_content = float(parts[2])
                n_count = 0
                if len(parts) > 3:
                    n_count = parts[3].upper().count("N")
                n_content = (n_count / seq_length * 100) if seq_length > 0 else 0
                stats[seq_name] = {
                    "genome_length": seq_length,
                    "gc_content": gc_content,
                    "percent_ambiguous_bases": round(n_content, 2),
                }
        log.info(f"seqkit stats complete ({len(stats)} sequences)")
        return stats
    except FileNotFoundError:
        log.warning("seqkit not found. Install with: conda install -c bioconda seqkit")
        return None
    except subprocess.CalledProcessError as e:
        log.warning(f"seqkit failed: {e.stderr}")
        return None


def identify_missing_regions(sequence: str, min_run_length: int = MISSING_REGION_MIN_RUN) -> List[str]:
    """Identify contiguous N runs."""
    regions: List[str] = []
    in_missing = False
    start_pos: Optional[int] = None

    for i, base in enumerate(sequence.upper()):
        if base == "N":
            if not in_missing:
                in_missing = True
                start_pos = i + 1  # 1-indexed
        else:
            if in_missing and start_pos is not None:
                end_pos = i
                run_length = end_pos - start_pos + 1
                if run_length >= min_run_length:
                    regions.append(f"Position {start_pos}-{end_pos} ({run_length} bp)")
                in_missing = False

    if in_missing and start_pos is not None:
        end_pos = len(sequence)
        run_length = end_pos - start_pos + 1
        if run_length >= min_run_length:
            regions.append(f"Position {start_pos}-{end_pos} ({run_length} bp)")

    return regions


def compute_quality_metrics(
    sequence: str,
    expected_length: int = DEFAULT_EXPECTED_LENGTH,
    tolerance: float = DEFAULT_TOLERANCE,
    seqkit_stats: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Compute per-sample quality metrics."""
    if seqkit_stats:
        seq_length = seqkit_stats["genome_length"]
        n_content = seqkit_stats["percent_ambiguous_bases"]
        gc_content = seqkit_stats.get("gc_content")
    else:
        seq_length = len(sequence)
        n_count = sequence.upper().count("N")
        n_content = (n_count / seq_length * 100) if seq_length > 0 else 0
        gc_content = None

    if gc_content is None:
        seq_upper = sequence.upper()
        gc_count = seq_upper.count("G") + seq_upper.count("C")
        gc_content = (gc_count / seq_length * 100) if seq_length > 0 else 0.0

    completeness = (seq_length / expected_length * 100) if expected_length > 0 else 0
    min_length = expected_length * (1 - tolerance)
    max_length = expected_length * (1 + tolerance)

    if completeness >= 90 and n_content <= 5:
        quality_flag = "HIGH"
    elif completeness >= 70 and n_content <= 10:
        quality_flag = "MODERATE"
    else:
        quality_flag = "LOW"

    missing_regions = identify_missing_regions(sequence)

    return {
        "genome_length": seq_length,
        "percent_ambiguous_bases": round(n_content, 2),
        "gc_content_pct": round(gc_content, 2),
        "genome_completeness": round(completeness, 2),
        "quality_flag": quality_flag,
        "missing_regions": missing_regions,
        "within_expected_range": min_length <= seq_length <= max_length,
        "expected_length": expected_length,
        "length_tolerance": tolerance,
    }


def run_quality_control(
    fasta_path: Path,
    metadata_csv: Path,
    output_dir: Path,
    fastq_path: Optional[Path] = None,
    expected_length: int = DEFAULT_EXPECTED_LENGTH,
    tolerance: float = DEFAULT_TOLERANCE,
    schema_path: Optional[Path] = None,
) -> Dict[str, Any]:
    """Run Stage 0 input & quality control.

    Returns a dict with:
      - sample_ids
      - quality_metrics: {sample_id: metrics}
      - metadata_rows
      - output_files
      - summary
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    schema = load_metadata_schema(schema_path)
    required_fields = get_required_fields(schema)

    sequences = parse_fasta(fasta_path)
    metadata = parse_metadata(metadata_csv)

    log.info(f"Found {len(sequences)} sequences and {len(metadata)} metadata rows")

    # Validate sample IDs
    is_valid, errors = validate_sample_id_matching(sequences, metadata)
    if not is_valid:
        for error in errors:
            log.error(error)
        raise ValueError("Sample ID mismatch between FASTA and metadata")
    log.info("All sample IDs match between FASTA and metadata")

    # Validate metadata schema
    is_valid, errors = validate_metadata_schema(metadata, required_fields)
    if not is_valid:
        for error in errors[:5]:
            log.warning(error)
        if len(errors) > 5:
            log.warning(f"... and {len(errors) - 5} more metadata errors")
    else:
        log.info("All required metadata fields present")

    # Optional FastQC
    fastqc_report: Optional[str] = None
    if fastq_path:
        fastqc_report = run_fastqc(fastq_path, output_dir)

    # seqkit stats
    seqkit_stats = run_seqkit_stats(fasta_path)

    # Compute metrics per sample
    quality_metrics: Dict[str, Any] = {}
    for sample_id, sequence in sequences.items():
        sample_stats = seqkit_stats.get(sample_id) if seqkit_stats else None
        metrics = compute_quality_metrics(
            sequence,
            expected_length=expected_length,
            tolerance=tolerance,
            seqkit_stats=sample_stats,
        )
        quality_metrics[sample_id] = metrics
        log.info(
            f"{sample_id}: length={metrics['genome_length']} bp, "
            f"completeness={metrics['genome_completeness']}%, flag={metrics['quality_flag']}"
        )

    high = sum(1 for m in quality_metrics.values() if m["quality_flag"] == "HIGH")
    moderate = sum(1 for m in quality_metrics.values() if m["quality_flag"] == "MODERATE")
    low = sum(1 for m in quality_metrics.values() if m["quality_flag"] == "LOW")

    summary = {
        "total_samples": len(quality_metrics),
        "high_quality": high,
        "moderate_quality": moderate,
        "low_quality": low,
        "fastqc_report": fastqc_report,
    }
    log.info(f"Stage 0 complete: total={summary['total_samples']}, HIGH={high}, MODERATE={moderate}, LOW={low}")

    # Write quality_metrics.json (includes aggregate summary under _summary)
    quality_output = output_dir / "quality_metrics.json"
    quality_metrics_with_summary = dict(quality_metrics)
    quality_metrics_with_summary["_summary"] = summary
    with open(quality_output, "w") as f:
        json.dump(quality_metrics_with_summary, f, indent=2)

    # Write validated_metadata.csv (enriched with QC metrics)
    metadata_output = output_dir / "validated_metadata.csv"
    if metadata:
        fieldnames = list(metadata[0].keys()) + [
            "genome_length",
            "gc_content_pct",
            "percent_ambiguous_bases",
            "genome_completeness",
            "quality_flag",
        ]
        with open(metadata_output, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for row in metadata:
                sample_id = row["sample_id"]
                if sample_id in quality_metrics:
                    row.update({
                        "genome_length": quality_metrics[sample_id]["genome_length"],
                        "gc_content_pct": quality_metrics[sample_id]["gc_content_pct"],
                        "percent_ambiguous_bases": quality_metrics[sample_id]["percent_ambiguous_bases"],
                        "genome_completeness": quality_metrics[sample_id]["genome_completeness"],
                        "quality_flag": quality_metrics[sample_id]["quality_flag"],
                    })
                writer.writerow(row)

    return {
        "sample_ids": list(quality_metrics.keys()),
        "quality_metrics": quality_metrics,
        "metadata_rows": metadata,
        "output_files": {
            "quality_metrics_json": str(quality_output),
            "validated_metadata_csv": str(metadata_output),
        },
        "summary": summary,
    }


def main():
    parser = argparse.ArgumentParser(description="Stage 0: Input & Quality Control")
    parser.add_argument("--metadata", required=True, help="Path to metadata CSV/TSV")
    parser.add_argument("--fasta", required=True, help="Path to FASTA file")
    parser.add_argument("--fastq", help="Path to FASTQ file (optional, future assembly)")
    parser.add_argument("--output", required=True, help="Output directory")
    parser.add_argument("--expected-length", type=int, default=DEFAULT_EXPECTED_LENGTH, help="Expected genome length")
    parser.add_argument("--tolerance", type=float, default=DEFAULT_TOLERANCE, help="Allowed length deviation")
    parser.add_argument("--schema", help="Path to metadata_schema.yaml")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    run_quality_control(
        fasta_path=Path(args.fasta),
        metadata_csv=Path(args.metadata),
        output_dir=Path(args.output),
        fastq_path=Path(args.fastq) if args.fastq else None,
        expected_length=args.expected_length,
        tolerance=args.tolerance,
        schema_path=Path(args.schema) if args.schema else None,
    )


if __name__ == "__main__":
    main()
