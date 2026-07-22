#!/usr/bin/env python3
"""
Copy a sample's genomic intelligence and bioinformatics outputs into a
consolidated report folder.

Creates:
  <output_dir>/<sample_id>/final_report.txt        (brief + report concatenated)
  <output_dir>/<sample_id>/figures/                (genomic-intelligence + bioinformatics figures)
  <output_dir>/<sample_id>/data/                   (bio_output.json + module JSON files)
  <output_dir>/<sample_id>/context_used.json       (grounding context for traceability)
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path


def _copy_bioinformatics_outputs(bioinformatics_dir: Path, report_dir: Path) -> dict[str, list[str]]:
    """Copy figures and JSON files from the bioinformatics sample directory."""
    copied: dict[str, list[str]] = {"figures": [], "data": []}
    if not bioinformatics_dir.is_dir():
        return copied

    figures_dst = report_dir / "figures" / "bioinformatics"
    figures_dst.mkdir(parents=True, exist_ok=True)
    data_dst = report_dir / "data" / "bioinformatics"
    data_dst.mkdir(parents=True, exist_ok=True)

    # Copy the consolidated bio_output.json at the top of data/
    bio_output = bioinformatics_dir / "bio_output.json"
    if bio_output.exists():
        shutil.copy2(bio_output, report_dir / "data" / "bio_output.json")
        copied["data"].append(str(report_dir / "data" / "bio_output.json"))

    # Copy module figures and JSON files, preserving subdirectory structure.
    for item in sorted(bioinformatics_dir.rglob("*")):
        if not item.is_file():
            continue
        rel = item.relative_to(bioinformatics_dir)
        if item.suffix == ".png":
            dst = figures_dst / rel.name
            counter = 1
            while dst.exists():
                dst = figures_dst / f"{rel.stem}_{counter}{rel.suffix}"
                counter += 1
            shutil.copy2(item, dst)
            copied["figures"].append(str(dst))
        elif item.suffix == ".json":
            dst = data_dst / rel.name
            counter = 1
            while dst.exists():
                dst = data_dst / f"{rel.stem}_{counter}{rel.suffix}"
                counter += 1
            shutil.copy2(item, dst)
            copied["data"].append(str(dst))

    return copied


def publish(genomic_dir: Path, report_dir: Path, bioinformatics_dir: Path | None = None) -> int:
    if not genomic_dir.is_dir():
        print(f"ERROR: genomic intelligence directory not found: {genomic_dir}", file=sys.stderr)
        return 1

    report_dir.mkdir(parents=True, exist_ok=True)

    brief = genomic_dir / "genomic_intelligence_brief.txt"
    report = genomic_dir / "genomic_intelligence_report.txt"
    final = report_dir / "final_report.txt"
    with open(final, "w") as out:
        if brief.exists():
            out.write(brief.read_text())
            out.write("\n\n")
        if report.exists():
            out.write(report.read_text())

    if not final.exists() or final.stat().st_size == 0:
        final.write_text("No narrative report was generated.\n")

    figures_src = genomic_dir / "figures"
    figures_dst = report_dir / "figures"
    if figures_src.exists():
        if figures_dst.exists():
            shutil.rmtree(figures_dst)
        shutil.copytree(figures_src, figures_dst)

    context_src = genomic_dir / "context_used.json"
    if context_src.exists():
        shutil.copy2(context_src, report_dir / "context_used.json")

    bio_copied: dict[str, list[str]] = {"figures": [], "data": []}
    if bioinformatics_dir:
        bio_copied = _copy_bioinformatics_outputs(bioinformatics_dir, report_dir)

    manifest = report_dir / "bioinformatics_manifest.json"
    try:
        manifest.write_text(json.dumps(bio_copied, indent=2))
    except Exception as exc:
        print(f"WARNING: could not write bioinformatics manifest: {exc}", file=sys.stderr)

    print(f"Published final report for {genomic_dir.name} -> {report_dir}")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description="Publish consolidated per-sample report")
    parser.add_argument("--genomic-dir", required=True, help="Path to genomic_intelligence/<sample_id> output")
    parser.add_argument("--output-dir", required=True, help="Target report directory")
    parser.add_argument("--bioinformatics-dir", default=None, help="Path to bioinformatics/<sample_id> output (optional)")
    args = parser.parse_args()
    return publish(
        Path(args.genomic_dir),
        Path(args.output_dir),
        Path(args.bioinformatics_dir) if args.bioinformatics_dir else None,
    )


if __name__ == "__main__":
    raise SystemExit(main())
