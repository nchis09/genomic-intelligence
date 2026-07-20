#!/usr/bin/env python3
"""
Stage 9: Normalised output assembler.

Takes the outputs from all upstream bioinformatics modules and assembles the
single `bio_output.json` consumed by the data engine and intelligence engine.

This module is intentionally thin: it does not perform new analysis; it only
harmonises and formats results produced by earlier stages.
"""

from __future__ import annotations

import csv
import datetime
import json
import logging
import math
import os
import sys
from pathlib import Path
from typing import Any, Optional

PROJECT_ROOT = Path(__file__).resolve().parents[4]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from config import DB_URL as CONFIG_DB_URL  # noqa: E402

from intelligence_engine.bioinformatics.modules import variant_calling  # noqa: E402

try:
    import psycopg2
    import psycopg2.extras
except Exception:  # pragma: no cover
    psycopg2 = None

log = logging.getLogger(__name__)

CONCERNING_PHENOTYPE_CATEGORIES = {
    "vaccine_escape",
    "immune_escape",
    "drug_resistance",
    "virulence",
    "disease_severity",
    "increased_transmission",
    "host_adaptation",
}


def _load_json(path: Path) -> dict[str, Any]:
    with open(path) as f:
        return json.load(f)


def _load_csv_metadata(path: Path) -> list[dict[str, Any]]:
    text = path.read_text()
    dialect = csv.Sniffer().sniff(text[:2048], delimiters=",\t")
    return list(csv.DictReader(text.splitlines(), dialect=dialect))


def _to_iso_date(value: Optional[str]) -> Optional[str]:
    """Convert common date strings to ISO-8601 YYYY-MM-DD."""
    if not value:
        return None
    value = value.strip()
    for fmt in ("%Y-%m-%d", "%d/%m/%Y", "%m/%d/%Y", "%d-%m-%Y"):
        try:
            from datetime import datetime
            return datetime.strptime(value, fmt).strftime("%Y-%m-%d")
        except ValueError:
            continue
    return value


def _redact_sequences(reference_summary: dict[str, Any]) -> dict[str, Any]:
    """Return a copy of the reference summary with protein sequences removed
    from the embedded proteome to keep bio_output.json compact."""
    summary = dict(reference_summary)
    if "reference_proteome" in summary:
        summary["reference_proteome"] = {
            gene: {k: v for k, v in prot.items() if k != "sequence"}
            for gene, prot in summary["reference_proteome"].items()
        }
    return summary


def _find_metadata(metadata_rows: list[dict[str, Any]], sample_id: str) -> dict[str, Any]:
    for row in metadata_rows:
        if row.get("sample_id") == sample_id:
            return row
    return {}


def _pct_identity_from_snps(snps: int, ref_length: int) -> float:
    if ref_length <= 0:
        return 0.0
    return round(100.0 * (ref_length - snps) / ref_length, 2)


def _get_concerning_mutations(
    species_id: str, db_url: Optional[str] = None
) -> list[dict[str, Any]]:
    """Return curated mutations with concerning phenotype associations."""
    if psycopg2 is None:
        return []

    db_url = db_url or os.environ.get("PGIRL_DB_URL") or CONFIG_DB_URL
    try:
        with psycopg2.connect(db_url) as conn:
            with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
                cur.execute(
                    """
                    SELECT DISTINCT gp.protein AS gene, gp.position, gp.ref_aa, gp.alt_aa
                    FROM genotype_phenotype gp
                    WHERE gp.species_id = %s
                      AND gp.phenotype_category = ANY(%s)
                      AND gp.verification_status != 'rejected';
                    """,
                    (species_id, list(CONCERNING_PHENOTYPE_CATEGORIES)),
                )
                return [dict(r) for r in cur.fetchall()]
    except Exception as exc:
        log.warning(f"Could not query concerning mutations: {exc}")
        return []


def _annotate_mutations_with_domains(
    mutations: list[dict[str, Any]],
    species_id: str,
    db_url: Optional[str] = None,
) -> list[dict[str, Any]]:
    """Add domain and hotspot annotations from the gene_function table."""
    if psycopg2 is None or not mutations:
        return mutations

    db_url = db_url or os.environ.get("PGIRL_DB_URL") or CONFIG_DB_URL
    genes = list({m["gene"] for m in mutations if m.get("gene")})
    try:
        with psycopg2.connect(db_url) as conn:
            with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
                cur.execute(
                    """
                    SELECT gene, key_domains, known_hotspots
                    FROM gene_function
                    WHERE species_id = %s AND gene = ANY(%s);
                    """,
                    (species_id, genes),
                )
                gene_info = {row["gene"]: row for row in cur.fetchall()}
    except Exception as exc:
        log.warning(f"Could not query gene_function: {exc}")
        return mutations

    annotated = []
    for m in mutations:
        gene = m.get("gene", "")
        pos = m.get("position")
        info = gene_info.get(gene, {})
        key_domains = info.get("key_domains") or []
        hotspots = info.get("known_hotspots") or []

        in_domain = False
        domain_name = None
        is_hotspot = False

        for domain in key_domains:
            if not isinstance(domain, dict):
                continue
            start = domain.get("start")
            end = domain.get("end")
            name = domain.get("name") or domain.get("domain")
            if start is not None and end is not None and pos is not None:
                if start <= pos <= end:
                    in_domain = True
                    domain_name = name

        for hotspot in hotspots:
            if not isinstance(hotspot, dict):
                continue
            if hotspot.get("position") == pos:
                is_hotspot = True
                if not domain_name:
                    domain_name = hotspot.get("domain")

        annotated_m = dict(m)
        annotated_m["domain"] = domain_name
        annotated_m["in_domain"] = in_domain
        annotated_m["known_hotspot"] = is_hotspot
        annotated.append(annotated_m)

    return annotated


def _build_negative_findings(
    sample_mutations: list[dict[str, Any]],
    concerning_mutations: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """List concerning curated mutations that are NOT present in the sample."""
    sample_keys = {
        (m.get("gene"), m.get("position"), m.get("ref_aa"), m.get("alt_aa"))
        for m in sample_mutations
    }
    findings = []
    for cm in concerning_mutations:
        key = (cm.get("gene"), cm.get("position"), cm.get("ref_aa"), cm.get("alt_aa"))
        if key not in sample_keys:
            hgvs = f"{cm['gene']}:{cm['ref_aa']}{cm['position']}{cm['alt_aa']}"
            findings.append(
                {
                    "gene": cm["gene"],
                    "position": cm["position"],
                    "ref_aa": cm["ref_aa"],
                    "alt_aa": cm["alt_aa"],
                    "hgvs_p": hgvs,
                    "status": "NOT_detected",
                    "note": "Concerning curated variant not observed in this sample",
                }
            )
    return findings


def assemble_bio_output(
    sample_id: str,
    metadata_row: dict[str, Any],
    quality_metrics: dict[str, Any],
    classification: dict[str, Any],
    reference_summary: dict[str, Any],
    lineage_result: dict[str, Any],
    output_dir: Path,
    db_url: Optional[str] = None,
    quality_metrics_json: Optional[Path] = None,
    validated_metadata_csv: Optional[Path] = None,
    phylogenetic_tree_result: Optional[dict[str, Any]] = None,
) -> dict[str, Any]:
    """Assemble the full `bio_output.json` for one sample."""
    ref_genome = reference_summary.get("reference_genome", {})
    ref_length = ref_genome.get("genome_length") or quality_metrics.get("expected_length", 0)
    sample_qc = quality_metrics.get(sample_id, quality_metrics)

    # Extract aggregate QC summary written by quality_control module
    qc_summary = quality_metrics.pop("_summary", {}) if isinstance(quality_metrics, dict) else {}

    # Variant calling: amino-acid variants are primary, nucleotide changes are kept
    # as supporting evidence.
    mutations = variant_calling.extract_amino_acid_variants(lineage_result)

    species_id = classification.get("species_id") or lineage_result.get("species_id", "")
    mutations = _annotate_mutations_with_domains(mutations, species_id, db_url=db_url)

    concerning = _get_concerning_mutations(species_id, db_url=db_url)
    negative_findings = _build_negative_findings(mutations, concerning)

    total_subs = lineage_result.get("qc", {}).get("total_substitutions", 0)
    identity = _pct_identity_from_snps(total_subs, ref_length)

    genome_quality = {
        "completeness_pct": sample_qc.get("genome_completeness", 0.0),
        "gc_content_pct": sample_qc.get("gc_content_pct"),
        "depth": sample_qc.get("mean_depth") if "mean_depth" in sample_qc else None,
        "flag": sample_qc.get("quality_flag", "UNKNOWN"),
        "missing": sample_qc.get("missing_regions", []),
    }

    collection_date = _to_iso_date(metadata_row.get("collection_date"))
    symptom_onset_date = _to_iso_date(metadata_row.get("symptom_onset_date"))

    metadata = {
        "country": metadata_row.get("country"),
        "admin1": metadata_row.get("admin1"),
        "admin2": metadata_row.get("admin2"),
        "collection_date": collection_date,
        "host": metadata_row.get("host"),
        "host_species": metadata_row.get("host_species"),
        "sample_type": metadata_row.get("sample_type"),
        "vaccination_status": metadata_row.get("vaccination_status"),
        "travel_history": metadata_row.get("travel_history"),
        "travel_locations": metadata_row.get("travel_locations"),
        "outcome": metadata_row.get("outcome"),
        "epi_link_id": metadata_row.get("epi_link_id"),
        "suspected_exposure": metadata_row.get("suspected_exposure"),
        "symptom_onset_date": symptom_onset_date,
    }

    stage9 = {
        "sample_id": sample_id,
        "pathogen": classification.get("pathogen_id", "ebola"),
        "species": classification.get("species", lineage_result.get("species", "")),
        "species_id": species_id,
        "lineage": lineage_result.get("lineage", ""),
        "clade": lineage_result.get("clade", ""),
        "closest_reference": {
            "accession": ref_genome.get("accession"),
            "name": ref_genome.get("accession"),
            "identity_pct": identity,
            "snps": total_subs,
        },
        "closest_outbreak_genome": None,  # populated by phylogenetics module later
        "mutations": mutations,
        "negative_findings": negative_findings,
        "genome_quality": genome_quality,
        "metadata": metadata,
        "tree_file": str(output_dir / "tree.nwk"),
        "recombination": {"detected": False, "events": []},
        "comparative": {
            "gene_content": "annotation pending",
            "gc_content": None,
            "dn_ds": None,
            "selection_pressure": None,
        },
        "collection_country": metadata_row.get("country"),
        "collection_date": collection_date,
    }

    bio_output = {
        "pipeline_version": "1.0.0",
        "pipeline_run_timestamp": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "sample": {
            "sample_id": sample_id,
            "pathogen_id": classification.get("pathogen_id", "ebola"),
            "species_id": species_id,
            "species_name": classification.get("species", lineage_result.get("species", "")),
            "pathogen_family": classification.get("pathogen_family", ""),
            "pathogen_genus": classification.get("pathogen_genus", ""),
        },
        "stage0_quality_control": {
            "input_format": "consensus_fasta",
            "genome_length_bp": sample_qc.get("genome_length"),
            "expected_length_bp": sample_qc.get("expected_length"),
            "genome_completeness_pct": sample_qc.get("genome_completeness"),
            "gc_content_pct": sample_qc.get("gc_content_pct"),
            "n_content_pct": sample_qc.get("percent_ambiguous_bases"),
            "mean_depth": sample_qc.get("mean_depth"),
            "quality_flag": sample_qc.get("quality_flag"),
            "missing_regions": sample_qc.get("missing_regions", []),
            "within_expected_range": sample_qc.get("within_expected_range"),
            "length_tolerance": sample_qc.get("length_tolerance"),
            "aggregate_summary": {
                "total_samples": qc_summary.get("total_samples"),
                "high_quality": qc_summary.get("high_quality"),
                "moderate_quality": qc_summary.get("moderate_quality"),
                "low_quality": qc_summary.get("low_quality"),
                "fastqc_report": qc_summary.get("fastqc_report"),
            },
            "output_files": {
                "quality_metrics_json": str(quality_metrics_json) if quality_metrics_json else None,
                "validated_metadata_csv": str(validated_metadata_csv) if validated_metadata_csv else None,
            },
        },
        "stage1_classification": classification,
        "stage2_reference_context": _redact_sequences(reference_summary),
        "stage5_lineage_clade": lineage_result,
        "stage6_phylogenetic_tree": phylogenetic_tree_result or {
            "method": None,
            "tree_file": None,
            "sequences_in_tree": None,
            "time_scaled_tree": None,
            "note": "No phylogenetic tree result available.",
        },
        "stage9_normalised_output": stage9,
    }

    return bio_output


def write_bio_output(
    sample_id: str,
    metadata_csv: Path,
    quality_metrics_json: Path,
    classification_json: Path,
    reference_summary_json: Path,
    lineage_result_json: Path,
    output_dir: Path,
    db_url: Optional[str] = None,
    phylogenetic_tree_result: Optional[dict[str, Any]] = None,
) -> Path:
    """Convenience wrapper that loads upstream JSON/CSV outputs and writes
    `bio_output.json` to output_dir."""
    metadata_rows = _load_csv_metadata(metadata_csv)
    metadata_row = _find_metadata(metadata_rows, sample_id)

    quality_metrics = _load_json(quality_metrics_json)
    classification = _load_json(classification_json).get(sample_id, {})
    reference_summary = _load_json(reference_summary_json)
    lineage_result = _load_json(lineage_result_json).get(sample_id, {})

    bio_output = assemble_bio_output(
        sample_id=sample_id,
        metadata_row=metadata_row,
        quality_metrics=quality_metrics,
        classification=classification,
        reference_summary=reference_summary,
        lineage_result=lineage_result,
        output_dir=output_dir,
        db_url=db_url,
        quality_metrics_json=quality_metrics_json,
        validated_metadata_csv=metadata_csv,
        phylogenetic_tree_result=phylogenetic_tree_result,
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    out_path = output_dir / "bio_output.json"
    out_path.write_text(json.dumps(bio_output, indent=2, default=str))
    log.info(f"Wrote bio_output.json to {out_path}")
    return out_path


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Stage 9: Assemble bio_output.json")
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--metadata", required=True, help="Metadata CSV/TSV")
    parser.add_argument("--quality-metrics", required=True, help="stage0 quality_metrics.json")
    parser.add_argument("--classification", required=True, help="stage1 species_id.json")
    parser.add_argument("--reference-summary", required=True, help="stage2 reference_selection_summary.json")
    parser.add_argument("--lineage-result", required=True, help="stage5 parsed Nextclade result JSON")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--db-url", default=None)
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    path = write_bio_output(
        sample_id=args.sample_id,
        metadata_csv=Path(args.metadata),
        quality_metrics_json=Path(args.quality_metrics),
        classification_json=Path(args.classification),
        reference_summary_json=Path(args.reference_summary),
        lineage_result_json=Path(args.lineage_result),
        output_dir=Path(args.output_dir),
        db_url=args.db_url,
    )
    print(path)


if __name__ == "__main__":
    main()
