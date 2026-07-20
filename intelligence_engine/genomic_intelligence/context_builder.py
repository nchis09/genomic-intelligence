"""
genomic_intelligence/context_builder.py — Assemble the grounding context for
LLM-based evidence synthesis.

Curates the intelligence_object.json (from evidence_integration) and
evidence_package.json (harmonized evidence + cross-evidence statistics) into
a single, size-bounded JSON bundle that is fed verbatim to the LLM as the
*only* source of truth. Deliberately excludes risk tiers, evidence-weighted
threat scores, and public-health-implication content: this stage synthesizes
evidence, it does not recommend actions.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any, Optional

# Keys that carry risk-tier / recommendation-flavored content and must not
# reach the synthesis prompt.
_EXCLUDED_ANALYSIS_KEYS = {"evidence_weighted_threat"}
_EXCLUDED_TOP_LEVEL_KEYS = {"risk_assessment"}

# Hard caps on how many findings / rows make it into the prompt.  Evidence
# packages for well-studied pathogens (e.g. EBOV) can contain hundreds of
# thousands of cross-evidence pairings; feeding them verbatim stalls synthesis.
_MAX_CONTEXT_FINDINGS = 20
_MAX_CSV_ROWS = 50
_MAX_PHENOTYPE_ANCHORS = 40
_MAX_VARIANT_ANCHORS = 25


def _compact_finding(f: dict) -> dict:
    return {
        "assertion": f.get("assertion"),
        "source": f.get("source"),
        "confidence": f.get("confidence"),
        "biological_relevance": f.get("biological_relevance") or None,
        "epidemiological_relevance": f.get("epidemiological_relevance") or None,
        "finding_type": f.get("finding_type"),
        "supporting_refs": f.get("supporting_refs") or None,
        "record_flagged": f.get("record_flagged", False),
    }


def _compact_analysis(analysis: dict) -> dict:
    findings = analysis.get("findings", [])[:_MAX_CONTEXT_FINDINGS]
    return {
        "title": analysis.get("title"),
        "summary": analysis.get("summary") or None,
        "findings": [_compact_finding(f) for f in findings],
        "gaps": analysis.get("gaps") or None,
        "metrics": analysis.get("metrics") or None,
    }


def _compact_statistical_finding(f: dict) -> dict:
    return {
        "metric": f.get("metric"),
        "value": f.get("value"),
        "method": f.get("method") or None,
        "sample_size": f.get("sample_size"),
        "p_value": f.get("p_value"),
        "supporting_evidence": f.get("supporting_evidence") or None,
        "notes": f.get("notes") or None,
    }


def _compact_evidence_object(obj: dict) -> dict:
    return {
        "key": obj.get("key"),
        "level": obj.get("level"),
        "variant": obj.get("variant") or None,
        "lineage": obj.get("lineage") or None,
        "phenotype_associations": obj.get("phenotype_associations") or None,
        "historical_outbreaks": obj.get("historical_outbreaks") or None,
        "molecular_epidemiology": obj.get("molecular_epidemiology") or None,
        "epidemiological_context": obj.get("epidemiological_context") or None,
        "interventions": obj.get("interventions") or None,
        "sources": obj.get("sources") or None,
    }


def _collect_figures(evidence_integration_dir: Path) -> list:
    """Scan the evidence_integration output directory for generated figure
    images (evidence_package/*.png and figures/*.png) and number them
    consistently so the LLM can reference them as "Ref: Figure N – Title"
    without inventing figure numbers or filenames."""
    candidates = []
    figures_subdir = evidence_integration_dir / "figures"
    if figures_subdir.is_dir():
        candidates.extend(sorted(figures_subdir.glob("*.png")))
    ep_subdir = evidence_integration_dir / "evidence_package"
    if ep_subdir.is_dir():
        candidates.extend(sorted(ep_subdir.glob("*.png")))

    figures = []
    for i, path in enumerate(candidates, start=1):
        stem = path.stem
        # Strip a leading numeric prefix like "01_" before title-casing.
        parts = stem.split("_")
        if parts and parts[0].isdigit():
            parts = parts[1:]
        title = " ".join(parts).replace("-", " ").title() or stem
        figures.append({
            "figure_number": i,
            "title": title,
            "filename": path.name,
            "source_path": str(path),
        })
    return figures


def build_context(
    intelligence_object: dict,
    evidence_package: Optional[dict] = None,
    bio_output: Optional[dict] = None,
    data_query_results: Optional[dict] = None,
    epi_output: Optional[dict] = None,
    analysis_outputs_dir: Optional[Path] = None,
    evidence_integration_dir: Optional[Path] = None,
    max_chars: Optional[int] = 50000,
) -> dict:
    """Build the curated, evidence-only context dict for LLM synthesis.

    Args:
        intelligence_object: Parsed intelligence_object.json from
            evidence_integration.pipeline.intelligence_pipeline.
        evidence_package: Parsed evidence_package.json from the same run
            (harmonization + cross-evidence analysis output), if available.
        bio_output: Optional bio_output.json from the bioinformatics pipeline.
        data_query_results: Optional db_query_results.json from the data engine.
        epi_output: Optional epi_output.json from online epidemiological querying.
        analysis_outputs_dir: Optional directory containing evidence_integration
            per-analysis CSV outputs.
        max_chars: Soft cap on the serialized context size; if exceeded,
            the largest lists are trimmed evenly so the prompt stays within
            the LLM's context window.

    Returns:
        A dict ready to be ``json.dumps``'d directly into the LLM prompt.
    """
    io = intelligence_object or {}

    analyses = {
        key: _compact_analysis(value)
        for key, value in (io.get("analyses") or {}).items()
        if key not in _EXCLUDED_ANALYSIS_KEYS
    }

    # The anchors already contain the essential sample, variant and phenotype
    # numbers. Keep only a compact summary of the sample QC details and the
    # epidemiological summary at the top level to avoid duplication.
    sample = io.get("sample") or {}
    compact_sample = {
        "sample_id": sample.get("sample_id"),
        "genome_qc": sample.get("genome_qc") or None,
    } if sample.get("genome_qc") else None

    context: dict[str, Any] = {
        "quantitative_anchors": _extract_quantitative_anchors(io, evidence_package),
        "bioinformatics_summary": _compact_bio_output(bio_output),
        "data_query_summary": _compact_data_query(data_query_results),
        "epidemiological_query_summary": _compact_epi_output(epi_output),
        "epidemiological_summary": _compact_epi_summary(io.get("epi_summary")),
        "genomic_analyses": analyses,
        "negative_findings": io.get("negative_findings"),
        "data_quality_warnings": io.get("data_quality_warnings"),
        "tree_available": io.get("tree_available"),
        "references": io.get("references"),
    }
    if compact_sample:
        context["sample_qc_details"] = compact_sample

    if evidence_package:
        context["cross_evidence_statistics"] = {
            k: {
                "title": v.get("title"),
                "findings": [
                    _compact_statistical_finding(f)
                    for f in (v.get("findings", [])[:_MAX_CONTEXT_FINDINGS])
                ],
                "data_gaps": v.get("data_gaps") or None,
                "metrics": v.get("metrics") or None,
            }
            for k, v in (evidence_package.get("cross_evidence_analysis") or {}).items()
        }

    if analysis_outputs_dir and Path(analysis_outputs_dir).is_dir():
        context["analysis_output_tables"] = _compact_analysis_outputs_csvs(Path(analysis_outputs_dir))

    if evidence_integration_dir and Path(evidence_integration_dir).is_dir():
        figures = _collect_figures(Path(evidence_integration_dir))
        if figures:
            context["figures"] = figures

    context = _drop_empty(context)

    # Add a catalog of real source/metric keys the model is allowed to cite.
    context["valid_citation_keys"] = _extract_valid_citation_keys(io, evidence_package, context)

    serialized = json.dumps(context, default=str)
    if max_chars and len(serialized) > max_chars:
        context = _trim_context(context, max_chars)

    return context


def _extract_quantitative_anchors(io: dict, evidence_package: Optional[dict]) -> dict:
    """Surface the most important numeric anchors at the top of the prompt.

    This gives the LLM an explicit, compact table of the key quantitative
    values it must use for cross-domain synthesis: sample quality, reference
    distances, variant frequencies, phenotype evidence, lineage metrics,
    cross-evidence statistics, and data quality flags. Do NOT invent values
    here; only copy numbers that already exist in the evidence.
    """
    sample = io.get("sample") or {}
    qc = sample.get("genome_qc") or {}
    closest = sample.get("closest_reference") or {}
    lineage_meta = sample.get("lineage_metadata") or {}

    anchors: dict[str, Any] = {
        "sample_id": sample.get("sample_id"),
        "pathogen_species": sample.get("species"),
        "species_id": sample.get("species_id"),
        "lineage": sample.get("lineage"),
        "clade": sample.get("clade"),
        "collection_date": sample.get("collection_date"),
        "country": sample.get("country"),
        "admin1": sample.get("admin1") or None,
        "admin2": sample.get("admin2") or None,
        "host": sample.get("host"),
        "genome_completeness_pct": qc.get("genome_completeness_pct") or sample.get("genome_completeness_pct"),
        "mean_depth": qc.get("mean_depth") or sample.get("mean_depth"),
        "quality_flag": qc.get("quality_flag") or sample.get("quality_flag"),
        "total_aa_variants": sample.get("total_aa_variants"),
        "closest_reference_name": closest.get("name"),
        "closest_reference_identity_pct": closest.get("identity_pct"),
        "closest_reference_snps": closest.get("snps"),
        "lineage_total_genomes": lineage_meta.get("total_genomes"),
        "lineage_first_detected": lineage_meta.get("first_detected"),
        "lineage_last_detected": lineage_meta.get("last_detected"),
        "lineage_countries_reported": lineage_meta.get("countries_reported"),
        "travel_history": sample.get("travel_history"),
        "travel_locations": sample.get("travel_locations"),
        "suspected_exposure": sample.get("suspected_exposure"),
    }

    # Variant-level anchors: cap to a high-signal subset so the synthesis prompt
    # does not drown out the output-template instruction. Prioritise variants
    # with curated database context (frequency/countries) or known-hotspot flags.
    variant_anchors = []
    for v in (io.get("variants") or []):
        ctx = v.get("_curated_context") or {}
        variant_anchors.append({
            "hgvs_p": v.get("hgvs_p"),
            "gene": v.get("gene"),
            "ref_aa": v.get("ref_aa"),
            "position": v.get("position"),
            "alt_aa": v.get("alt_aa"),
            "known_hotspot": v.get("known_hotspot", False),
            "domain": v.get("domain") or None,
            "frequency": ctx.get("frequency") or None,
            "genome_count": ctx.get("genome_count") or None,
            "total_genomes": ctx.get("total_genomes") or None,
            "first_seen_date": ctx.get("first_seen_date") or None,
            "last_seen_date": ctx.get("last_seen_date") or None,
            "countries_seen": ctx.get("countries_seen") or None,
            "lineage_ids": ctx.get("lineage_ids") or None,
        })
    total_variants = len(variant_anchors)
    if total_variants > _MAX_VARIANT_ANCHORS:
        scored = [
            (v, (1 if v.get("known_hotspot") else 0) +
                (1 if v.get("frequency") else 0) +
                (1 if v.get("countries_seen") else 0))
            for v in variant_anchors
        ]
        scored.sort(key=lambda x: x[1], reverse=True)
        variant_anchors = [v for v, _ in scored[:_MAX_VARIANT_ANCHORS]]
    anchors["total_aa_variants"] = total_variants or anchors.get("total_aa_variants")
    anchors["variants"] = variant_anchors
    if total_variants > _MAX_VARIANT_ANCHORS:
        anchors["variants_note"] = (
            f"Showing {_MAX_VARIANT_ANCHORS} of {total_variants} variants; "
            "prioritised curated/hotspot variants."
        )

    # Phenotype anchors grouped by variant
    phenotype_anchors = []
    seen = set()
    for p in (io.get("matched_phenotypes") or []):
        key = (p.get("hgvs_p"), p.get("phenotype_category"), p.get("phenotype_specific"))
        if key in seen:
            continue
        seen.add(key)
        ctx = p.get("curated_context") or {}
        phenotype_anchors.append({
            "hgvs_p": p.get("hgvs_p"),
            "gene": p.get("gene"),
            "phenotype_category": p.get("phenotype_category"),
            "phenotype_specific": p.get("phenotype_specific"),
            "evidence_strength": p.get("evidence_strength"),
            "literature_refs": p.get("literature_refs") or None,
            "frequency": ctx.get("frequency") or None,
            "genome_count": ctx.get("genome_count") or None,
            "total_genomes": ctx.get("total_genomes") or None,
        })
    total_phenotypes = len(phenotype_anchors)
    if total_phenotypes > _MAX_PHENOTYPE_ANCHORS:
        # Prioritise associations that carry a literature reference and a
        # stronger evidence tier so the highest-signal ones survive.
        strength_rank = {"strong": 3, "moderate": 2, "preliminary": 1}
        phenotype_anchors.sort(
            key=lambda a: (
                1 if a.get("literature_refs") else 0,
                strength_rank.get(str(a.get("evidence_strength")).lower(), 0),
            ),
            reverse=True,
        )
        phenotype_anchors = phenotype_anchors[:_MAX_PHENOTYPE_ANCHORS]
    anchors["phenotype_associations"] = phenotype_anchors
    if total_phenotypes > _MAX_PHENOTYPE_ANCHORS:
        anchors["phenotype_associations_note"] = (
            f"Showing {_MAX_PHENOTYPE_ANCHORS} of {total_phenotypes} phenotype associations; "
            "prioritised those with literature references and stronger evidence."
        )

    # Analysis metrics
    analysis_metrics: dict[str, Any] = {}
    for key, analysis in (io.get("analyses") or {}).items():
        if key in _EXCLUDED_ANALYSIS_KEYS:
            continue
        metrics = analysis.get("metrics") or {}
        if metrics:
            # Keep only scalar / small metric values, not huge nested dicts.
            compact_metrics = {}
            for mkey, mval in metrics.items():
                if isinstance(mval, (str, int, float, bool)) or mval is None:
                    compact_metrics[mkey] = mval
                elif isinstance(mval, list) and len(mval) <= 5:
                    compact_metrics[mkey] = mval
                elif isinstance(mval, dict):
                    # Keep small dicts (country -> year/count, etc.) but cap size
                    compact_metrics[mkey] = dict(list(mval.items())[:10])
            if compact_metrics:
                analysis_metrics[key] = compact_metrics
    anchors["analysis_metrics"] = analysis_metrics

    # Cross-evidence statistics
    cross_stats: dict[str, Any] = {}
    if evidence_package:
        for key, block in (evidence_package.get("cross_evidence_analysis") or {}).items():
            metrics = block.get("metrics") or {}
            findings_summary = []
            for f in (block.get("findings") or [])[:3]:
                summary = {
                    "metric": f.get("metric"),
                    "value": f.get("value"),
                    "method": f.get("method") or None,
                    "p_value": f.get("p_value") or None,
                    "sample_size": f.get("sample_size") or None,
                }
                findings_summary.append(_drop_empty(summary))
            cross_stats[key] = {
                "metrics": _drop_empty(metrics) if metrics else None,
                "top_findings": findings_summary,
            }
    anchors["cross_evidence_statistics"] = _drop_empty(cross_stats)

    anchors["data_quality_warnings"] = io.get("data_quality_warnings") or None
    anchors["negative_findings"] = io.get("negative_findings") or None

    return _drop_empty(anchors)


def _extract_valid_citation_keys(
    io: dict,
    evidence_package: Optional[dict],
    context: dict,
) -> dict:
    """Build a catalog of real keys the model can use for citations.

    This discourages invented metric/source names. The catalog includes top-level
    evidence sections, analysis names, source fields from findings, actual metric
    names from cross-evidence statistics, and PMIDs/literature references.
    """
    keys = {
        "evidence_sections": sorted([k for k in context.keys() if not k.startswith("_")]),
        "anchor_keys": sorted([k for k in (context.get("quantitative_anchors") or {}).keys()]),
        "analysis_names": sorted([k for k in (io.get("analyses") or {}).keys() if k not in _EXCLUDED_ANALYSIS_KEYS]),
        "cross_evidence_metric_names": [],
        "source_fields": set(),
        "literature_refs": set(),
    }

    for analysis in (io.get("analyses") or {}).values():
        for f in (analysis.get("findings") or []):
            src = f.get("source")
            if src:
                keys["source_fields"].add(src)
            ftype = f.get("finding_type")
            if ftype:
                keys["source_fields"].add(ftype)
            refs = f.get("supporting_refs") or []
            if isinstance(refs, list):
                keys["literature_refs"].update(str(r) for r in refs if r)

    if evidence_package:
        for block in (evidence_package.get("cross_evidence_analysis") or {}).values():
            for f in (block.get("findings", [])[:_MAX_CONTEXT_FINDINGS]):
                metric = f.get("metric")
                if metric:
                    keys["cross_evidence_metric_names"].append(metric)
                src = f.get("source") or block.get("title")
                if src:
                    keys["source_fields"].add(src)

    # Also collect PMIDs from matched phenotypes and anchors
    for p in (io.get("matched_phenotypes") or []):
        refs = p.get("literature_refs") or []
        if isinstance(refs, list):
            keys["literature_refs"].update(str(r) for r in refs if r)
    anchors = context.get("quantitative_anchors") or {}
    for p in (anchors.get("phenotype_associations") or []):
        refs = p.get("literature_refs") or []
        if isinstance(refs, list):
            keys["literature_refs"].update(str(r) for r in refs if r)

    keys["source_fields"] = sorted(keys["source_fields"])
    keys["literature_refs"] = sorted(keys["literature_refs"])
    return _drop_empty(keys)


def _drop_empty(obj: Any) -> Any:
    """Recursively drop None/empty values so the prompt isn't padded with noise."""
    if isinstance(obj, dict):
        cleaned = {k: _drop_empty(v) for k, v in obj.items()}
        return {k: v for k, v in cleaned.items() if v not in (None, {}, [], "")}
    if isinstance(obj, list):
        return [_drop_empty(v) for v in obj if v not in (None, {}, [], "")]
    return obj


def _trim_context(context: dict, max_chars: int) -> dict:
    """Trim the largest list found anywhere in the context until the serialized
    context fits within max_chars, noting the truncation explicitly so the LLM
    (and any human reviewer) knows evidence was omitted, not fabricated.

    The previous implementation only trimmed ``findings`` lists inside
    ``genomic_analyses``/``cross_evidence_statistics``. When the raw
    bioinformatics/data-query dumps dominated the payload, that loop would strip
    every analytical finding while leaving the raw dumps intact -- leaving the
    LLM with almost no interpretable evidence. This version walks the whole
    structure and pops from whichever list is currently the longest, so bulky
    raw sections are reduced first and the analytical findings survive.
    """
    trimmed_note = "Some lower-priority list entries were omitted to fit the context window."

    def _size() -> int:
        return len(json.dumps(context, default=str))

    # Never trim these small, high-value scalar sections.
    def _find_largest_list(node: Any, protect: bool = False):
        """Return (list_obj, length) for the longest list in the tree."""
        best = (None, 0)
        if isinstance(node, dict):
            for k, v in node.items():
                # Protect the top-level scalar anchors' identity fields, but
                # still allow trimming their large sub-lists (variants, etc.).
                child = _find_largest_list(v)
                if child[1] > best[1]:
                    best = child
        elif isinstance(node, list):
            if len(node) > best[1]:
                best = (node, len(node))
            for item in node:
                child = _find_largest_list(item)
                if child[1] > best[1]:
                    best = child
        return best

    guard = 0
    while _size() > max_chars and guard < 10000:
        guard += 1
        largest, length = _find_largest_list(context)
        if largest is None or length <= 1:
            break
        largest.pop()
        context["_truncated"] = trimmed_note

    return context


def _compact_bio_output(bio_output: Optional[dict]) -> Optional[dict]:
    """Extract the bioinformatics evidence most relevant for synthesis."""
    if not bio_output:
        return None

    qc = bio_output.get("stage0_quality_control") or {}
    classification = bio_output.get("stage1_classification") or {}
    lineage = bio_output.get("stage5_lineage_clade") or {}
    tree = bio_output.get("stage6_phylogenetic_tree") or {}
    norm = bio_output.get("stage9_normalised_output") or {}

    # The full per-mutation list can contain 1000+ entries with long
    # nucleotide-change arrays; the detected variants are already surfaced in
    # quantitative_anchors, so keep only a small preview here.
    norm_mutations = norm.get("mutations")
    if isinstance(norm_mutations, list) and len(norm_mutations) > _MAX_VARIANT_ANCHORS:
        norm_mutations = norm_mutations[:_MAX_VARIANT_ANCHORS]

    return _drop_empty({
        "quality_control": {
            "genome_length_bp": qc.get("genome_length_bp"),
            "expected_length_bp": qc.get("expected_length_bp"),
            "genome_completeness_pct": qc.get("genome_completeness_pct"),
            "gc_content_pct": qc.get("gc_content_pct"),
            "n_content_pct": qc.get("n_content_pct"),
            "mean_depth": qc.get("mean_depth"),
            "quality_flag": qc.get("quality_flag"),
            "missing_regions": qc.get("missing_regions"),
            "assembly_method": qc.get("assembly_method"),
        },
        "classification": {
            "method": classification.get("method"),
            "species": classification.get("species"),
            "species_id": classification.get("species_id"),
            "family": classification.get("pathogen_family"),
            "genus": classification.get("pathogen_genus"),
            "confidence": classification.get("confidence"),
            "best_hit": classification.get("best_hit") if isinstance(classification.get("best_hit"), dict) else None,
            "kraken2_species": classification.get("kraken2_species"),
            "kraken2_confidence": classification.get("kraken2_confidence"),
            "blast_species": classification.get("blast_species"),
            "blast_confidence": classification.get("blast_confidence"),
            "agreement": classification.get("agreement"),
        },
        "lineage_and_clade": {
            "lineage": lineage.get("lineage"),
            "clade": lineage.get("clade"),
        },
        "phylogenetic_tree": {
            "sequences_in_tree": tree.get("sequences_in_tree"),
            "tree_method": tree.get("tree_method"),
            "model_selected": tree.get("model_selected"),
            "bootstrap_replicates": tree.get("bootstrap_replicates"),
            "time_scaled_tree": tree.get("time_scaled_tree"),
        },
        "normalised_output": {
            "sample_id": norm.get("sample_id"),
            "pathogen": norm.get("pathogen"),
            "species": norm.get("species"),
            "lineage": norm.get("lineage"),
            "clade": norm.get("clade"),
            "collection_country": norm.get("collection_country"),
            "collection_date": norm.get("collection_date"),
            "metadata": norm.get("metadata"),
            "genome_quality": norm.get("genome_quality"),
            "closest_reference": norm.get("closest_reference"),
            "closest_outbreak_genome": norm.get("closest_outbreak_genome"),
            "mutations": norm_mutations,
            "comparative": norm.get("comparative"),
            "recombination": norm.get("recombination"),
        },
    })


def _compact_data_query(data: Optional[dict]) -> Optional[dict]:
    """Summarize the local database query results for synthesis.

    Keeps only counts and top-level summaries. The full variant/gene records are
    already represented in the intelligence_object/evidence_package.
    """
    if not data:
        return None

    ldb = data.get("local_db_results") or {}
    l1 = data.get("layer1_variant_lookup") or {}
    l2 = data.get("layer2_lineage_context") or {}
    l3 = data.get("layer3_geographic_temporal") or {}
    l4 = data.get("layer4_gene_function_context") or {}
    l5 = data.get("layer5_species_surveillance") or {}

    summaries = l1.get("variant_summaries") or []
    phenotypes = l1.get("variant_phenotypes") or []
    lineage_meta = (l2.get("lineage_metadata") or {}).get("lineage_metadata") or {}
    genomes_summary = l2.get("lineage_genomes_summary") or {}

    # Reduce a full variant summary record to a few key fields.
    def _mini_variant(rec: dict) -> dict:
        return {
            "hgvs_p": rec.get("hgvs_p"),
            "gene": rec.get("gene"),
            "global_frequency": rec.get("global_frequency"),
            "genome_count": rec.get("genome_count"),
            "species_total_genomes": rec.get("species_total_genomes"),
            "first_seen_date": rec.get("first_seen_date"),
            "last_seen_date": rec.get("last_seen_date"),
            "countries_seen": rec.get("countries_seen"),
            "lineage_ids": rec.get("lineage_ids"),
            "is_stop": rec.get("is_stop"),
            "variant_type": rec.get("variant_type"),
        }

    # Keep a small preview of variant summaries; the anchors already cover
    # detected variants, but this helps the model see the DB context.
    variant_preview = [_mini_variant(v) for v in summaries[:3]]

    # Several DB layers can return very large per-variant maps/lists (hundreds of
    # KB) that would otherwise dominate the prompt and starve the analytical
    # findings. Cap them to a small preview here.
    def _cap(value: Any, n: int = 10) -> Any:
        if isinstance(value, list):
            return value[:n]
        if isinstance(value, dict):
            return dict(list(value.items())[:n])
        return value

    # Limit lineage-defining variants preview.
    ldv = l2.get("lineage_defining_variants") or []
    ldv_preview = []
    for v in (ldv if isinstance(ldv, list) else [])[:10]:
        if isinstance(v, dict):
            ldv_preview.append({
                "hgvs_p": v.get("hgvs_p"),
                "gene": v.get("gene"),
                "frequency_in_lineage": v.get("frequency_in_lineage"),
                "is_defining": v.get("is_defining"),
            })

    return _drop_empty({
        "lineage_surveillance": {
            "lineage_last_seen_date": ldb.get("lineage_last_seen_date"),
            "lineage_countries_in_db": ldb.get("lineage_countries_in_db"),
            "lineage_not_in_db_for": ldb.get("lineage_not_in_db_for"),
            "variant_frequencies": _cap(ldb.get("variant_frequencies")),
        },
        "variant_lookup": {
            "variant_summaries_count": len(summaries),
            "variant_summaries_preview": variant_preview,
            "variant_phenotypes_count": len(phenotypes),
            "negative_findings": _cap(l1.get("negative_findings")) or None,
        },
        "lineage_context": {
            "lineage_metadata": {
                "lineage_id": lineage_meta.get("lineage_id"),
                "lineage_name": lineage_meta.get("lineage_name"),
                "parent_lineage": lineage_meta.get("parent_lineage"),
                "first_detected": lineage_meta.get("first_detected"),
                "last_detected": lineage_meta.get("last_detected"),
                "total_genomes": lineage_meta.get("total_genomes"),
                "countries_reported": lineage_meta.get("countries_reported"),
                "primary_host": lineage_meta.get("primary_host"),
            },
            "lineage_genomes_summary": genomes_summary,
            "lineage_defining_variants_preview": ldv_preview,
        },
        "geographic_temporal": {
            "lineage_in_country": l3.get("lineage_in_country"),
            "lineages_in_country": _cap(l3.get("lineages_in_country")),
            "variant_temporal": _cap(l3.get("variant_temporal")),
        },
        "gene_function": _compact_gene_function(l4.get("gene_function")),
        "species_surveillance": {
            "country_totals": (l5.get("country_totals") or [])[:10],
            "year_totals": (l5.get("year_totals") or [])[:10],
        },
    })


def _compact_epi_summary(epi_summary: Optional[dict]) -> Optional[dict]:
    """Compact the epidemiological summary to a preview of key rows.

    The detailed epidemiological information is already available in the
    epidemiological_query_summary and genomic_analyses; this just provides a
    lightweight index so the model can cross-reference outbreak/surveillance rows.
    """
    if not epi_summary:
        return None

    blocks = epi_summary.get("_raw_text_blocks") or []
    preview = []
    for block in blocks[:3]:
        rows = block.get("rows") or []
        preview.append({
            "question": block.get("question"),
            "preview_rows": rows[:3],
        })
    return _drop_empty({
        "questions_answered": epi_summary.get("questions_answered"),
        "total_epi_rows": epi_summary.get("total_epi_rows"),
        "preview_blocks": preview,
        "references": (epi_summary.get("references") or [])[:5],
    })


def _compact_gene_function(gf: Any) -> Any:
    """Compact the gene_function dict/list to a small summary."""
    if not gf:
        return None
    if isinstance(gf, dict):
        return {k: (v[:200] if isinstance(v, str) and len(v) > 200 else v) for k, v in gf.items()}
    if isinstance(gf, list):
        return gf[:3]
    return gf


def _compact_epi_output(epi: Optional[dict]) -> Optional[dict]:
    """Summarize online epidemiological query results for synthesis."""
    if not epi:
        return None

    def _top(items, n=2):
        if isinstance(items, list):
            return items[:n]
        return items

    return _drop_empty({
        "pathogen_profile": epi.get("pathogen_profile"),
        "molecular_epidemiology": _top(epi.get("molecular_epidemiology") or []),
        "outbreaks": _top(epi.get("outbreaks") or []),
        "transmission": epi.get("transmission"),
        "demographics": _top(epi.get("demographics") or [], n=10),
        "clinical": _top(epi.get("clinical") or [], n=5),
        "genomic_links": _top(epi.get("genomic_links") or [], n=10),
        "surveillance": _top(epi.get("surveillance") or []),
        "interventions": _top(epi.get("interventions") or []),
        "vaccines": _top(epi.get("vaccines") or []),
        "knowledge_assertions": epi.get("knowledge_assertions"),
        "references": epi.get("references"),
    })


def _compact_analysis_outputs_csvs(analysis_dir: Path) -> dict:
    """Read evidence_integration CSV outputs and convert them to compact tables.

    Each table is limited to a header plus a small number of rows so the prompt
    stays bounded. Risk/tier and combined-dump outputs are deliberately excluded.
    """
    tables: dict[str, Any] = {}
    excluded_files = {
        "risk_assessment.csv",
        "data_quality_warnings.csv",
        "combined_analysis_output.csv",
        "sample_metadata.csv",
    }

    csv_files = sorted([p for p in analysis_dir.glob("*.csv") if p.name not in excluded_files])[:10]
    for csv_path in csv_files:
        try:
            with open(csv_path, newline="", encoding="utf-8") as f:
                reader = csv.DictReader(f)
                rows = list(reader)
        except Exception:
            continue

        if not rows:
            continue
        # Keep only the most informative row for each table.
        limited = []
        for row in rows[:1]:
            limited.append({k: (v[:120] if isinstance(v, str) and len(v) > 120 else v) for k, v in row.items()})
        key = csv_path.stem
        tables[key] = {
            "columns": list(rows[0].keys()),
            "rows": limited,
            "total_rows": len(rows),
        }
    return _drop_empty(tables)


def build_context_str(
    intelligence_object: dict,
    evidence_package: Optional[dict] = None,
    bio_output: Optional[dict] = None,
    data_query_results: Optional[dict] = None,
    epi_output: Optional[dict] = None,
    analysis_outputs_dir: Optional[Path] = None,
    evidence_integration_dir: Optional[Path] = None,
    max_chars: Optional[int] = 50000,
) -> str:
    """Convenience wrapper returning the pretty-printed JSON string for the prompt."""
    return json.dumps(
        build_context(
            intelligence_object,
            evidence_package,
            bio_output=bio_output,
            data_query_results=data_query_results,
            epi_output=epi_output,
            analysis_outputs_dir=analysis_outputs_dir,
            evidence_integration_dir=evidence_integration_dir,
            max_chars=max_chars,
        ),
        indent=2,
        default=str,
    )
