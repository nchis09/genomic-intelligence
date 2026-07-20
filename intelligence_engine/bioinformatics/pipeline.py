#!/usr/bin/env python3
"""
Top-level bioinformatics pipeline for the PGIRL system.

Usage:
    /Users/christianndekezi/anaconda3/bin/python3 -m intelligence_engine.bioinformatics.pipeline \
        --fasta input/input_FASTA.fasta \
        --metadata input/input_metadata.csv \
        --output-dir output/bioinformatics

    # Optional explicit hints (normally auto-detected):
    #   --pathogen ebola --species-id EBOV

The pipeline performs:
  1. Input & quality control (stage 0)
  2. Taxonomic classification / species confirmation (stage 1)
  3. Reference selection & context gathering (stage 2)
  4. Pathogen-specific workflow dispatch
  5. Normalised `bio_output.json` per sample

Currently only the Ebola workflow is fully wired. Other pathogens can be added
by implementing a module under pathogen_workflows/.
"""

from __future__ import annotations

import argparse
import csv
import importlib
import json
import logging
import sys
from pathlib import Path
from typing import Any, Optional

PROJECT_ROOT = Path(__file__).resolve().parents[3]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

log = logging.getLogger(__name__)

SUPPORTED_PATHOGENS = ["ebola", "dengue", "influenza", "rvf", "mpox"]
SUPPORTED_EBOLA_SPECIES = ["EBOV", "SUDV", "BDBV", "RESTV", "TAFV", "BOMV"]


def _infer_pathogen_from_metadata(metadata_path: Path) -> Optional[str]:
    """Return pathogen from metadata if a supported value is present.

    Looks for 'pathogen' or 'pathogen_id' columns. Falls back to None if the
    column is absent or contains unsupported values.
    """
    if not metadata_path.exists():
        return None
    try:
        with open(metadata_path, newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                for key in ("pathogen", "pathogen_id"):
                    value = (row.get(key) or "").strip().lower()
                    if value and value in SUPPORTED_PATHOGENS:
                        return value
        return None
    except Exception as exc:
        log.warning(f"Could not read metadata for pathogen inference: {exc}")
        return None


def load_pathogen_workflow(pathogen: str):
    """Dynamically import the workflow module for the requested pathogen.

    Each workflow module must expose:
      - PATHOGEN_ID (str)
      - run_workflow(**kwargs) -> dict[str, Path]
      - REQUIRED_ARGS (list[str], optional)
    """
    module_name = f"intelligence_engine.bioinformatics.pathogen_workflows.{pathogen}"
    try:
        module = importlib.import_module(module_name)
    except ImportError as exc:
        raise ValueError(
            f"Pathogen workflow '{pathogen}' is not implemented. "
            f"Supported: {', '.join(SUPPORTED_PATHOGENS)}"
        ) from exc
    return module


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="PGIRL bioinformatics pipeline (consensus FASTA -> bio_output.json)"
    )
    parser.add_argument("--fasta", required=True, help="Input multi-FASTA of consensus genomes")
    parser.add_argument("--metadata", required=True, help="Sample metadata CSV/TSV")
    parser.add_argument(
        "--pathogen",
        choices=SUPPORTED_PATHOGENS,
        default=None,
        help="Pathogen-specific workflow to run. If omitted, the pipeline infers it "
             "from a 'pathogen' or 'pathogen_id' column in the metadata, or falls back to ebola.",
    )
    parser.add_argument(
        "--species-id",
        default=None,
        help="Optional species hint (e.g. EBOV, SUDV). When omitted, stage1 classification "
             "(BLAST) detects the species and downstream stages use the detected species."
    )
    parser.add_argument(
        "--output-dir",
        default="output/bioinformatics",
        help="Root output directory (default: output/bioinformatics)",
    )
    parser.add_argument(
        "--db-url",
        default=None,
        help="PostgreSQL URL (default: config.DB_URL or PGIRL_DB_URL env var)",
    )
    parser.add_argument("--skip-qc", action="store_true", help="Skip stage 0 QC if outputs exist")
    parser.add_argument("--skip-classification", action="store_true", help="Skip stage1 classification and trust --species-id")
    parser.add_argument("--skip-blast", action="store_true", help="Skip NCBI BLASTn classification")
    parser.add_argument("--skip-kraken", action="store_true", help="Skip Kraken2 classification")
    parser.add_argument("--debug", action="store_true", help="Enable debug logging")
    return parser.parse_args()


def _validate_workflow_args(workflow: Any, args: argparse.Namespace) -> None:
    """Raise if a workflow's required arguments are missing."""
    required = getattr(workflow, "REQUIRED_ARGS", [])
    for arg in required:
        value = getattr(args, arg, None)
        if not value:
            raise ValueError(
                f"Workflow '{workflow.PATHOGEN_ID}' requires --{arg.replace('_', '-')}"
            )


def main() -> int:
    args = parse_args()
    logging.basicConfig(
        level=logging.DEBUG if args.debug else logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
    )

    fasta = Path(args.fasta)
    metadata = Path(args.metadata)
    output_dir = Path(args.output_dir)

    if not fasta.exists():
        log.error(f"FASTA file not found: {fasta}")
        return 1
    if not metadata.exists():
        log.error(f"Metadata file not found: {metadata}")
        return 1

    output_dir.mkdir(parents=True, exist_ok=True)

    # Resolve pathogen: explicit flag > metadata column > default Ebola
    pathogen = args.pathogen or _infer_pathogen_from_metadata(metadata) or "ebola"
    log.info(f"Using pathogen workflow: {pathogen}")

    # Write a manifest for reproducibility
    manifest = {
        "fasta": str(fasta),
        "metadata": str(metadata),
        "pathogen": pathogen,
        "species_id": args.species_id,
        "output_dir": str(output_dir),
        "db_url": args.db_url,
    }
    manifest_path = output_dir / "pipeline_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2))
    log.info(f"Pipeline manifest written to {manifest_path}")

    workflow = load_pathogen_workflow(pathogen)
    try:
        _validate_workflow_args(workflow, args)
    except ValueError as exc:
        log.error(str(exc))
        return 1

    # Classification is now run by default (BLASTn + optional Kraken2).
    # Use --skip-* flags to fall back to the user-provided --species-id.
    results = workflow.run_workflow(
        fasta_path=fasta,
        metadata_csv=metadata,
        output_dir=output_dir,
        db_url=args.db_url,
        skip_qc=args.skip_qc,
        skip_classification=args.skip_classification,
        skip_blast=args.skip_blast,
        skip_kraken=args.skip_kraken,
        species_id=args.species_id,
    )

    # Summary
    print("=" * 60)
    print("BIOINFORMATICS PIPELINE COMPLETE")
    print("=" * 60)
    for sample_id, bio_output_path in results.items():
        print(f"  {sample_id}: {bio_output_path}")

    summary_path = output_dir / "pipeline_summary.json"
    summary_path.write_text(
        json.dumps({sid: str(p) for sid, p in results.items()}, indent=2)
    )
    log.info(f"Summary written to {summary_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
