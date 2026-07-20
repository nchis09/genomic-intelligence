"""
run_per_sample_batch.py — Run the single-sample intelligence pipeline for every
bioinformatics sample independently.

Each sample is treated as its own entity: it gets its own bio_output,
db_query_results, epi_output, tree, and final intelligence_object.json + outputs.
No grouping or aggregation across samples is performed.

Usage:
    /Users/christianndekezi/anaconda3/bin/python3 \
        -m intelligence_engine.evidence_integration.pipeline.run_per_sample_batch \
        --bioinformatics-dir output/bioinformatics \
        --data-query-dir output/data_query \
        --output-dir output/genomic_intelligence
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import Optional

from intelligence_engine.evidence_integration.pipeline.intelligence_pipeline import (
    IntelligencePipeline,
)

log = logging.getLogger(__name__)


def run_batch(
    bioinformatics_dir: Path,
    data_query_dir: Path,
    output_dir: Path,
    associations: Path,
    variants: Path,
    lineages: Path,
    genome_metadata: Path,
    db_url: Optional[str] = None,
) -> list[Path]:
    """Run the single-sample intelligence pipeline for every sample directory."""
    output_dir.mkdir(parents=True, exist_ok=True)
    completed: list[Path] = []

    if not bioinformatics_dir.is_dir():
        log.error("Bioinformatics directory not found: %s", bioinformatics_dir)
        sys.exit(1)

    sample_dirs = sorted(
        d for d in bioinformatics_dir.iterdir() if d.is_dir() and not d.name.startswith("_")
    )
    if not sample_dirs:
        log.warning("No sample directories found under %s", bioinformatics_dir)
        return completed

    for sample_dir in sample_dirs:
        sample_id = sample_dir.name
        bio_path = sample_dir / "bio_output.json"
        epi_path = data_query_dir / sample_id / "epi_output.json"
        tree_path = sample_dir / "tree.nwk"
        sample_output_dir = output_dir / sample_id

        if not bio_path.exists():
            log.warning("[%s] Skipping — no bio_output.json", sample_id)
            continue
        if not epi_path.exists():
            log.warning("[%s] Skipping — no epi_output.json in %s", sample_id, data_query_dir)
            continue

        log.info("[%s] Running intelligence pipeline...", sample_id)
        try:
            pipeline = IntelligencePipeline(
                epi_output_path=str(epi_path),
                bio_output_path=str(bio_path),
                associations_csv_path=str(associations),
                protein_variants_csv_path=str(variants),
                lineages_csv_path=str(lineages),
                genome_metadata_csv_path=str(genome_metadata),
                tree_file_path=str(tree_path) if tree_path.exists() else None,
                db_url=db_url,
            )
            result = pipeline.run(output_dir=str(sample_output_dir))
            log.info("[%s] Wrote intelligence object: %s", sample_id, result["json_path"])
            completed.append(Path(result["json_path"]))
        except Exception as exc:
            log.error("[%s] Intelligence pipeline failed: %s", sample_id, exc)
            continue

    log.info("Batch complete. Processed %d sample(s).", len(completed))
    return completed


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run the single-sample intelligence pipeline for every sample independently."
    )
    parser.add_argument(
        "--bioinformatics-dir",
        default="output/bioinformatics",
        help="Directory containing per-sample bioinformatics output folders.",
    )
    parser.add_argument(
        "--data-query-dir",
        default="output/data_query",
        help="Directory containing per-sample db_query_results.json and epi_output.json.",
    )
    parser.add_argument(
        "--output-dir",
        default="output/evidence_integration",
        help="Directory to write per-sample evidence integration outputs (intelligence object, evidence package, figures, analysis outputs).",
    )
    parser.add_argument(
        "--associations",
        default="database/exports/genotype_phenotype.csv",
        help="Path to genotype-phenotype associations CSV.",
    )
    parser.add_argument(
        "--variants",
        default="database/exports/protein_variants.csv",
        help="Path to protein variant frequency CSV.",
    )
    parser.add_argument(
        "--lineages",
        default="database/exports/lineages.csv",
        help="Path to lineage metadata CSV.",
    )
    parser.add_argument(
        "--genome-metadata",
        default="database/exports/genome_metadata.csv",
        help="Path to genome metadata CSV.",
    )
    parser.add_argument(
        "--db-url",
        default=None,
        help="PostgreSQL connection URL (defaults to PGIRL_DB_URL env var, then config.DB_URL).",
    )
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    run_batch(
        bioinformatics_dir=Path(args.bioinformatics_dir),
        data_query_dir=Path(args.data_query_dir),
        output_dir=Path(args.output_dir),
        associations=Path(args.associations),
        variants=Path(args.variants),
        lineages=Path(args.lineages),
        genome_metadata=Path(args.genome_metadata),
        db_url=args.db_url,
    )


if __name__ == "__main__":
    main()
