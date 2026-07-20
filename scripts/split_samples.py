#!/usr/bin/env python3
"""
split_samples.py — Split a multi-sample FASTA + metadata CSV into per-sample
FASTA + metadata files for Nextflow per-sample processing.

Usage:
    python scripts/split_samples.py \
        --fasta input/input_FASTA.fasta \
        --metadata input/input_metadata.csv \
        --output-dir output/split_input

Output:
    output/split_input/<sample_id>/sample.fasta
    output/split_input/<sample_id>/metadata.csv
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Split multi-sample inputs into per-sample files")
    parser.add_argument("--fasta", required=True, help="Multi-sample FASTA")
    parser.add_argument("--metadata", required=True, help="Sample metadata CSV")
    parser.add_argument("--output-dir", required=True, help="Directory to write per-sample files")
    return parser.parse_args()


def read_fasta(path: Path) -> dict[str, str]:
    """Return dict {sample_id: sequence} from a multi-FASTA."""
    records: dict[str, str] = {}
    current_id: str | None = None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current_id = line[1:].split()[0].strip()
                records[current_id] = []
            elif current_id is not None:
                records[current_id].append(line)
    return {sid: "".join(seq) for sid, seq in records.items()}


def main() -> int:
    args = parse_args()
    fasta_path = Path(args.fasta)
    metadata_path = Path(args.metadata)
    output_dir = Path(args.output_dir)

    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA not found: {fasta_path}")
    if not metadata_path.exists():
        raise FileNotFoundError(f"Metadata not found: {metadata_path}")

    sequences = read_fasta(fasta_path)

    with open(metadata_path, newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    fieldnames = reader.fieldnames or []
    output_dir.mkdir(parents=True, exist_ok=True)

    for row in rows:
        sample_id = row.get("sample_id", "").strip()
        if not sample_id:
            raise ValueError("Metadata row missing sample_id")
        if sample_id not in sequences:
            raise KeyError(f"FASTA has no record for sample_id '{sample_id}'")

        sample_dir = output_dir / sample_id
        sample_dir.mkdir(parents=True, exist_ok=True)

        (sample_dir / "sample.fasta").write_text(f">{sample_id}\n{sequences[sample_id]}\n")

        sample_meta = sample_dir / "metadata.csv"
        with open(sample_meta, "w", newline="") as out:
            writer = csv.DictWriter(out, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerow(row)

    print(f"Split {len(rows)} sample(s) into {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
