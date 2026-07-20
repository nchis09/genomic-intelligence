#!/usr/bin/env python3
"""
Ebola-specific bioinformatics workflow.

Chains the reusable modules for quality control, reference selection, lineage
assignment (Nextclade), phylogenetics, comparative genomics and normalisation
into a single per-sample pipeline. All outputs are written under
output_dir/<sample_id>/.
"""

from __future__ import annotations

import csv
import json
import logging
import math
import sys
from pathlib import Path
from typing import Any, Optional

PROJECT_ROOT = Path(__file__).resolve().parents[4]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from intelligence_engine.bioinformatics.modules import (  # noqa: E402
    comparative_genomics,
    lineage_assignment,
    normalisation,
    phylogenetics,
    quality_control,
    reference_selection,
    taxonomic_classification,
    variant_calling,
)

log = logging.getLogger(__name__)

# Pathogen workflow metadata used by the top-level dispatcher.
PATHOGEN_ID = "ebola"
SUPPORTED_SPECIES = ["EBOV", "SUDV", "BDBV", "RESTV", "TAFV", "BOMV"]
REQUIRED_ARGS: list[str] = []


def run_ebola_workflow(
    fasta_path: Path,
    metadata_csv: Path,
    output_dir: Path,
    species_id: Optional[str] = None,
    db_url: Optional[str] = None,
    skip_qc: bool = False,
    skip_classification: bool = False,
    skip_blast: bool = False,
    skip_kraken: bool = False,
) -> dict[str, Path]:
    """Run the Ebola bioinformatics pipeline for all samples in the FASTA.

    Args:
        fasta_path: Multi-FASTA of consensus genomes.
        metadata_csv: Metadata CSV/TSV with one row per sample.
        output_dir: Root output directory (e.g. output/bioinformatics).
        species_id: Optional species hint (EBOV, SUDV, BDBV, ...). If omitted or
            BLAST succeeds, the species detected by classification is used for
            reference selection and Nextclade lineage assignment.
        db_url: Optional PostgreSQL URL.
        skip_qc: If True, reuse existing stage0 outputs if present.
        skip_classification: If True, skip stage1 (species_id hint is trusted).
        skip_blast: If True, skip NCBI BLASTn classification.
        skip_kraken: If True, skip Kraken2 classification.

    Returns:
        Dict mapping sample_id -> bio_output.json path.
    """
    from Bio import SeqIO

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # Stage 0: Quality control
    # ------------------------------------------------------------------
    stage0_dir = output_dir / "_stage0_qc"
    stage0_dir.mkdir(parents=True, exist_ok=True)

    if not skip_qc or not (stage0_dir / "quality_metrics.json").exists():
        quality_control.run_quality_control(
            fasta_path=fasta_path,
            metadata_csv=metadata_csv,
            output_dir=stage0_dir,
            expected_length=18900,
        )

    quality_metrics_json = stage0_dir / "quality_metrics.json"
    validated_metadata_csv = stage0_dir / "validated_metadata.csv"

    # ------------------------------------------------------------------
    # Stage 1: Taxonomic classification (BLAST default, Kraken2 optional)
    # ------------------------------------------------------------------
    stage1_dir = output_dir / "_stage1_classification"
    stage1_dir.mkdir(parents=True, exist_ok=True)
    species_id_json = stage1_dir / "species_id.json"

    if not skip_classification or not species_id_json.exists():
        if not skip_classification:
            # Run actual taxonomic classification (Kraken2 + NCBI BLAST).
            # BLAST is remote and can be slow; Kraken2 runs if a local DB exists.
            taxonomic_classification.run_taxonomic_classification(
                fasta_path=fasta_path,
                output_dir=stage1_dir,
                skip_kraken=skip_kraken,
                skip_blast=skip_blast,
            )

        if not species_id_json.exists():
            if species_id:
                # Fallback: build classification from user-provided species_id, but
                # derive a meaningful confidence score from QC metrics instead of
                # hard-coding "HIGH".
                qc_metrics: dict[str, Any] = {}
                if quality_metrics_json.exists():
                    try:
                        qc_metrics = json.loads(quality_metrics_json.read_text())
                    except json.JSONDecodeError:
                        qc_metrics = {}

                classifications = {}
                for record in SeqIO.parse(fasta_path, "fasta"):
                    sample_id = record.id
                    sample_qc = qc_metrics.get(sample_id, {})
                    completeness = sample_qc.get("genome_completeness", 0.0)
                    n_content = sample_qc.get("percent_ambiguous_bases", 100.0)
                    within_range = sample_qc.get("within_expected_range", False)

                    # Weighted quality score (0-1)
                    score = (
                        (completeness / 100.0) * 0.50
                        + max(0.0, 1.0 - (n_content / 5.0)) * 0.30
                        + (1.0 if within_range else 0.0) * 0.20
                    )
                    score = round(min(1.0, max(0.0, score)), 4)

                    if score >= 0.90 and completeness >= 98.0 and n_content <= 1.0 and within_range:
                        confidence = "HIGH"
                    elif score >= 0.70:
                        confidence = "MODERATE"
                    else:
                        confidence = "LOW"

                    classifications[record.id] = {
                        "sample_id": record.id,
                        "species": _species_name(species_id),
                        "species_id": species_id,
                        "pathogen_id": "ebola",
                        "pathogen_family": "Filoviridae",
                        "pathogen_genus": "Orthoebolavirus",
                        "confidence": confidence,
                        "confidence_score": score,
                        "confidence_metrics": {
                            "genome_completeness_pct": completeness,
                            "n_content_pct": n_content,
                            "genome_length_bp": sample_qc.get("genome_length"),
                            "expected_length_bp": sample_qc.get("expected_length"),
                            "within_expected_range": within_range,
                        },
                        "method": "user_specified_with_qc_validation",
                    }
                species_id_json.write_text(json.dumps(classifications, indent=2, default=str))
            else:
                raise FileNotFoundError(
                    f"No classification output found at {species_id_json} and no --species-id hint provided. "
                    "Run classification or supply a species hint."
                )

    # Ensure every classification record has the structured fields expected downstream.
    classifications = json.loads(species_id_json.read_text())
    for sample_id, c in classifications.items():
        c.setdefault("species_id", _species_id_from_name(c.get("species", "")))
        c.setdefault("pathogen_id", "ebola")
        c.setdefault("pathogen_family", "Filoviridae")
        c.setdefault("pathogen_genus", "Orthoebolavirus")
    species_id_json.write_text(json.dumps(classifications, indent=2))

    # Resolve the species_id to use for each sample. BLAST-detected species takes
    # precedence; the user's --species-id hint is only a fallback.
    sample_species_ids = {
        sid: _effective_species_id(c, species_id)
        for sid, c in classifications.items()
    }
    unique_species_ids = sorted(set(sample_species_ids.values()))

    # ------------------------------------------------------------------
    # Stage 2: Reference selection & context gathering (per detected species)
    # ------------------------------------------------------------------
    stage2_dir = output_dir / "_stage2_reference_selection"
    stage2_dir.mkdir(parents=True, exist_ok=True)

    ref_summaries: dict[str, dict[str, Any]] = {}
    ref_summary_json_by_species: dict[str, Path] = {}
    for sid in unique_species_ids:
        species_ref_dir = stage2_dir / sid
        species_ref_dir.mkdir(parents=True, exist_ok=True)
        ref_summary_json = species_ref_dir / "reference_selection_summary.json"
        ref_summary_json_by_species[sid] = ref_summary_json
        if not ref_summary_json.exists():
            log.info(f"Running reference selection for detected species {sid}")
            reference_selection.run_reference_selection(
                species_id=sid,
                output_dir=species_ref_dir,
                db_url=db_url,
            )
        ref_summaries[sid] = json.loads(ref_summary_json.read_text())

    # ------------------------------------------------------------------
    # Stages 3-6, 8: per-sample modules
    # ------------------------------------------------------------------
    results: dict[str, Path] = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        sample_id = record.id
        sample_dir = output_dir / sample_id
        sample_dir.mkdir(parents=True, exist_ok=True)

        sample_fasta = sample_dir / "sample.fasta"
        if not sample_fasta.exists() or sample_fasta.stat().st_mtime < fasta_path.stat().st_mtime:
            SeqIO.write([record], sample_fasta, "fasta")

        effective_species_id = sample_species_ids.get(sample_id, species_id or "EBOV")
        ref_summary = ref_summaries[effective_species_id]
        ref_summary_json = ref_summary_json_by_species[effective_species_id]

        # Stage 5: Lineage / clade assignment via Nextclade
        stage5_dir = sample_dir / "_stage5_nextclade"
        stage5_dir.mkdir(parents=True, exist_ok=True)
        lineage_result_json = stage5_dir / "nextclade_parsed.json"

        if not lineage_result_json.exists():
            lineage_assignment.run_nextclade_analysis(
                fasta_path=sample_fasta,
                species_id=effective_species_id,
                output_dir=stage5_dir,
            )

        lineage_result = json.loads(lineage_result_json.read_text()).get(sample_id, {})

        # Stage 6: Phylogenetics (currently reuses Nextclade placement tree)
        stage6_dir = sample_dir / "_stage6_phylogenetics"
        stage6_dir.mkdir(parents=True, exist_ok=True)
        tree_file = sample_dir / "tree.nwk"

        if not tree_file.exists():
            phylo_result = phylogenetics.run_phylogenetics(
                sample_fasta=sample_fasta,
                reference_fasta=Path(ref_summary["reference_fasta"]),
                context_genomes=ref_summary.get("context_genomes", []),
                output_dir=stage6_dir,
                lineage_result=lineage_result,
            )
            if phylo_result.get("tree_file"):
                tree_file.symlink_to(Path(phylo_result["tree_file"]).resolve())
        else:
            # Reconstruct minimal stage6 result from the existing tree file.
            phylo_result = {
                "method": "nextclade_placement_stub",
                "tree_file": str(tree_file.resolve()),
                "sequences_in_tree": phylogenetics._count_tips(tree_file),
                "time_scaled_tree": None,
                "note": "Reused existing Nextclade placement tree.",
            }

        # Stage 4: Variant calling (amino-acid-centric, with NT support)
        stage4_dir = sample_dir / "_stage4_variant_calling"
        stage4_dir.mkdir(parents=True, exist_ok=True)
        variant_call_result = variant_calling.run_variant_calling(
            lineage_result=lineage_result,
            reference_fasta=Path(ref_summary["reference_fasta"]),
            reference_proteome=ref_summary.get("reference_proteome", {}),
        )
        variant_call_result["sample_id"] = sample_id
        variant_call_result["raw_output_file"] = str(stage4_dir / "variant_calling_result.json")
        stage4_json = stage4_dir / "variant_calling_result.json"
        stage4_json.write_text(json.dumps(variant_call_result, indent=2, default=str))

        # Stage 8: Comparative genomics
        comp_result = comparative_genomics.run_comparative_genomics(
            sample_fasta=sample_fasta,
            reference_fasta=Path(ref_summary["reference_fasta"]),
            gene_coordinates=ref_summary.get("reference_genome", {}).get("gene_coordinates", {}),
        )

        # Stage 9: Normalisation
        bio_output_path = normalisation.write_bio_output(
            sample_id=sample_id,
            metadata_csv=validated_metadata_csv,
            quality_metrics_json=quality_metrics_json,
            classification_json=species_id_json,
            reference_summary_json=ref_summary_json,
            lineage_result_json=lineage_result_json,
            output_dir=sample_dir,
            db_url=db_url,
            phylogenetic_tree_result=phylo_result,
        )

        # Enrich bio_output with stage 4 and stage 8 outputs
        bio_output = json.loads(bio_output_path.read_text())
        bio_output["stage4_variant_calling"] = variant_call_result
        bio_output["stage8_comparative_genomics"] = comp_result
        bio_output["stage9_normalised_output"]["comparative"] = {
            "gene_content": comp_result.get("gene_content"),
            "gc_content": comp_result.get("gc_content"),
            "dn_ds": comp_result.get("dn_ds"),
            "selection_pressure": comp_result.get("selection_pressure"),
        }
        bio_output_path.write_text(json.dumps(bio_output, indent=2, default=str))

        results[sample_id] = bio_output_path
        log.info(f"Sample {sample_id}: {bio_output_path}")

    return results


def _species_name(species_id: str) -> str:
    mapping = {
        "EBOV": "Zaire ebolavirus",
        "SUDV": "Sudan ebolavirus",
        "BDBV": "Bundibugyo ebolavirus",
        "RESTV": "Reston ebolavirus",
        "TAFV": "Tai Forest ebolavirus",
        "BOMV": "Bombali ebolavirus",
    }
    return mapping.get(species_id, species_id)


def _species_id_from_name(name: str) -> str:
    mapping = {v: k for k, v in {
        "EBOV": "Zaire ebolavirus",
        "SUDV": "Sudan ebolavirus",
        "BDBV": "Bundibugyo ebolavirus",
        "RESTV": "Reston ebolavirus",
        "TAFV": "Tai Forest ebolavirus",
        "BOMV": "Bombali ebolavirus",
    }.items()}
    return mapping.get(name, "")


def _effective_species_id(
    classification: dict[str, Any],
    default_species_id: Optional[str] = None,
) -> str:
    """Return the species_id to use for a sample.

    Prefer the species resolved by stage1 classification, then the user's
    species hint, then a conservative default.
    """
    species_id = _safe_str(classification.get("species_id"))
    if species_id and species_id in SUPPORTED_SPECIES:
        return species_id
    species_name = _safe_str(classification.get("species"))
    species_id = _species_id_from_name(species_name)
    if species_id and species_id in SUPPORTED_SPECIES:
        return species_id
    if default_species_id and default_species_id in SUPPORTED_SPECIES:
        return default_species_id
    return "EBOV"


def _safe_str(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return "" if math.isnan(value) else str(value).strip()
    return str(value).strip()


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Ebola bioinformatics workflow")
    parser.add_argument("--fasta", required=True, help="Input consensus FASTA")
    parser.add_argument("--metadata", required=True, help="Sample metadata CSV/TSV")
    parser.add_argument(
        "--species-id",
        default=None,
        help="Optional species hint (EBOV, SUDV, ...). If omitted, classification will detect it.",
    )
    parser.add_argument("--output-dir", required=True, help="Root output directory")
    parser.add_argument("--db-url", default=None, help="PostgreSQL URL")
    parser.add_argument("--skip-qc", action="store_true", help="Reuse existing QC outputs")
    parser.add_argument("--skip-classification", action="store_true", help="Skip stage1 classification and trust --species-id")
    parser.add_argument("--skip-blast", action="store_true", help="Skip NCBI BLASTn classification")
    parser.add_argument("--skip-kraken", action="store_true", help="Skip Kraken2 classification")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    results = run_ebola_workflow(
        fasta_path=Path(args.fasta),
        metadata_csv=Path(args.metadata),
        output_dir=Path(args.output_dir),
        species_id=args.species_id,
        db_url=args.db_url,
        skip_qc=args.skip_qc,
        skip_classification=args.skip_classification,
        skip_blast=args.skip_blast,
        skip_kraken=args.skip_kraken,
    )

    print("=" * 60)
    print("EBOLA BIOINFORMATICS WORKFLOW COMPLETE")
    print("=" * 60)
    for sample_id, path in results.items():
        print(f"  {sample_id}: {path}")


# Generic dispatcher alias so pipeline.py can call module.run_workflow().
run_workflow = run_ebola_workflow


if __name__ == "__main__":
    main()
