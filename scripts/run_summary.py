#!/usr/bin/env python3
"""
Generate a JSON run summary from the pipeline output directory.

The summary lists each sample, key genomic findings, pipeline stage status,
and paths to the final consolidated report files.
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path


def _load_json(path: Path) -> dict | None:
    try:
        with open(path) as f:
            return json.load(f)
    except Exception:
        return None


def _sample_summary(outdir: Path, sample_id: str) -> dict:
    io_path = outdir / "evidence_integration" / sample_id / "intelligence_object.json"
    io = _load_json(io_path) or {}
    sample = io.get("sample") or {}

    risk = io.get("risk_assessment") or {}
    risk_tier = risk.get("risk_tier", "unknown")

    report_dir = outdir / "reports" / sample_id
    report_path = report_dir / "final_report.txt" if report_dir.exists() else None

    return {
        "sample_id": sample_id,
        "status": "completed" if io_path.exists() else "failed",
        "pathogen": sample.get("species") or sample.get("pathogen") or "unknown",
        "lineage": sample.get("lineage") or "unknown",
        "clade": sample.get("clade") or "unknown",
        "country": sample.get("country") or "unknown",
        "collection_date": sample.get("collection_date") or "unknown",
        "quality_flag": sample.get("quality_flag") or "unknown",
        "genome_completeness_pct": sample.get("genome_completeness_pct"),
        "mean_depth": sample.get("mean_depth"),
        "total_aa_variants": sample.get("total_aa_variants"),
        "matched_phenotypes_count": len(io.get("matched_phenotypes", [])),
        "risk_tier": risk_tier,
        "data_quality_warnings_count": len(io.get("data_quality_warnings", [])),
        "report_file": str(report_path) if report_path and report_path.exists() else None,
        "evidence_integration_file": str(io_path) if io_path.exists() else None,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate PGIRL run summary")
    parser.add_argument("--outdir", required=True, help="Pipeline output directory to scan for results")
    parser.add_argument("--output", default="run_summary.json", help="Where to write the summary JSON")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    if not outdir.is_dir():
        print(f"ERROR: output directory not found: {outdir}", file=sys.stderr)
        return 1

    evidence_dir = outdir / "evidence_integration"
    sample_ids = sorted([d.name for d in evidence_dir.iterdir() if d.is_dir()]) if evidence_dir.exists() else []

    samples = [_sample_summary(outdir, sid) for sid in sample_ids]
    completed = [s for s in samples if s["status"] == "completed"]
    failed = [s for s in samples if s["status"] != "completed"]

    summary = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "output_directory": str(outdir.resolve()),
        "total_samples": len(samples),
        "completed": len(completed),
        "failed": len(failed),
        "samples": samples,
    }

    summary_path = Path(args.output)
    summary_path.write_text(json.dumps(summary, indent=2))
    print(f"Run summary written to {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
