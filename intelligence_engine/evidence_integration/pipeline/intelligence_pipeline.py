#!/usr/bin/env python3
"""Assessment-logic pipeline that builds the structured intelligence object.

This pipeline loads the epi + bioinformatics inputs, runs all assessment analyzers,
writes `intelligence_object.json`, exports plain CSV analysis outputs, and
(optionally) generates decision-oriented R figures. Report generation (LLM/template)
consumes the intelligence object downstream.
"""

from __future__ import annotations

import argparse
import json
import logging
import re
import shutil
import subprocess
import sys
import tempfile
from dataclasses import asdict
from pathlib import Path
from typing import Any, Optional

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from intelligence_engine.evidence_integration.engine import (  # noqa: E402
    GenomicIntelligenceEngine,
    _match_lineage,
    _parse_pg_array,
)
from intelligence_engine.evidence_integration.harmonization import (  # noqa: E402
    build_evidence_objects,
)
from intelligence_engine.evidence_integration.cross_evidence import (  # noqa: E402
    ConfidenceScoringAnalyzer,
    GeographicDistributionAnalyzer,
    InterventionAssociationAnalyzer,
    LineagePhenotypeAssociationAnalyzer,
    MutationCooccurrenceAnalyzer,
    MutationPersistenceAnalyzer,
    TemporalTrendAnalyzer,
)
from intelligence_engine.evidence_integration.visualization import (  # noqa: E402
    save_evidence_network,
    save_geographic_distribution_chart,
    save_temporal_trend_chart,
)

log = logging.getLogger(__name__)

# Vocabulary kept in sync with database/_vocabularies/risk_tiers.yaml
RISK_TIERS = {
    "routine": 1,
    "monitor": 2,
    "investigate": 3,
    "high_priority": 4,
    "emergency": 5,
}

CONCERNING_PHENOTYPES = {
    "vaccine_escape",
    "vaccine_effectiveness",
    "drug_resistance",
    "immune_escape",
    "virulence",
    "disease_severity",
    "host_adaptation",
}

# Single-letter -> three-letter amino-acid code for matching HGVS p. strings.
AA_3LETTER = {
    "A": "Ala",
    "C": "Cys",
    "D": "Asp",
    "E": "Glu",
    "F": "Phe",
    "G": "Gly",
    "H": "His",
    "I": "Ile",
    "K": "Lys",
    "L": "Leu",
    "M": "Met",
    "N": "Asn",
    "P": "Pro",
    "Q": "Gln",
    "R": "Arg",
    "S": "Ser",
    "T": "Thr",
    "V": "Val",
    "W": "Trp",
    "Y": "Tyr",
}


def _load_json(path: Path) -> dict:
    with open(path) as f:
        return json.load(f)


def _clean_for_json(obj: Any) -> Any:
    """Recursively make an object JSON-serialisable (NaN/NaT/inf -> None)."""
    import math

    if isinstance(obj, dict):
        return {k: _clean_for_json(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_clean_for_json(v) for v in obj]
    if isinstance(obj, float):
        if math.isnan(obj) or math.isinf(obj):
            return None
        return obj
    if obj is pd.NaT:
        return None
    if isinstance(obj, pd.Timestamp):
        return obj.isoformat()
    return obj


def _extract_variants(stage9: dict, bio: Optional[dict] = None) -> list[dict]:
    """Return amino-acid variant records from the bioinformatics output.

    Prefer the normalised stage9 mutation list, but fall back to the richer
    stage4 variant calling table when stage9 is empty or incomplete.
    """
    variants = []
    for v in stage9.get("mutations", []):
        ref_aa = str(v.get("ref_aa", "")).strip().upper()
        alt_aa = str(v.get("alt_aa", "")).strip().upper()
        if ref_aa == alt_aa:
            # Skip synonymous/no-change positions; they carry no protein-level
            # signal for phenotype matching or risk assessment.
            continue
        variants.append(
            {
                "gene": str(v.get("gene", "")).strip(),
                "position": v.get("position"),
                "ref_aa": ref_aa,
                "alt_aa": alt_aa,
                "hgvs_p": v.get("hgvs_p"),
                "domain": v.get("domain"),
                "known_hotspot": bool(v.get("known_hotspot")),
            }
        )

    # Fall back to stage4 amino-acid variants if stage9 supplied no usable variants.
    if not variants and bio is not None:
        stage4 = bio.get("stage4_variant_calling", {})
        for v in stage4.get("amino_acid_variants", []):
            dc = v.get("domain_context", {})
            ref_aa = str(v.get("ref_aa", "")).strip().upper()
            alt_aa = str(v.get("alt_aa", "")).strip().upper()
            if ref_aa == alt_aa:
                continue
            variants.append(
                {
                    "gene": str(v.get("gene", "")).strip(),
                    "position": v.get("position"),
                    "ref_aa": ref_aa,
                    "alt_aa": alt_aa,
                    "hgvs_p": v.get("hgvs_p"),
                    "domain": dc.get("domain") if isinstance(dc, dict) else None,
                    "known_hotspot": dc.get("known_hotspot") if isinstance(dc, dict) else False,
                }
            )
    return variants


def _enrich_variants(variants: list[dict], protein_variants_df: pd.DataFrame) -> list[dict]:
    """Add curated frequency/country context to each variant."""
    df = protein_variants_df
    if df.empty:
        return variants

    for v in variants:
        gene = str(v.get("gene", "")).strip()
        pos = v.get("position")
        alt = str(v.get("alt_aa", "")).strip()
        if not gene or pd.isna(pos) or not alt:
            continue
        mask = (
            (df["gene"].str.upper() == gene.upper())
            & (df["position"] == pos)
            & (df["alt_aa"].str.upper() == alt.upper())
        )
        rows = df[mask]
        if not rows.empty:
            top = rows.iloc[0]
            v["_curated_context"] = {
                "genome_count": int(top.get("genome_count", 0) or 0),
                "total_genomes": int(top.get("species_total_genomes", 0) or 0),
                "first_seen_date": top.get("first_seen_date"),
                "last_seen_date": top.get("last_seen_date"),
                "countries_seen": top.get("countries_seen"),
                "lineage_ids": top.get("lineage_ids"),
                "frequency": _calc_frequency(
                    top.get("genome_count"), top.get("species_total_genomes")
                ),
            }
    return variants


def _calc_frequency(count, total) -> str:
    try:
        c = int(count or 0)
        t = int(total or 0)
        if t == 0:
            return "unknown"
        pct = 100 * c / t
        return f"{c}/{t} ({pct:.1f}%)"
    except (ValueError, TypeError):
        return "unknown"


def _to_hgvs_three_letter(hgvs_p: str) -> Optional[str]:
    """Convert e.g. 'GP:A82V' to 'p.Ala82Val'."""
    if not hgvs_p:
        return None
    m = re.match(r"([A-Za-z0-9_]+:)?([A-Z])(\d+)([A-Z])", str(hgvs_p))
    if not m:
        return None
    ref = AA_3LETTER.get(m.group(2))
    alt = AA_3LETTER.get(m.group(4))
    if not ref or not alt:
        return None
    return f"p.{ref}{m.group(3)}{alt}"


def _match_phenotypes(
    variants: list[dict],
    associations: pd.DataFrame,
    epi_summary: dict,
) -> list[dict]:
    """Join detected variants to curated genotype-phenotype associations.

    Three sources are considered:
      1. Exact amino-acid change in the associations CSV.
      2. Motif/domain-level associations when the variant lies in that domain.
      3. Phenotype associations already computed by the DB query engine.
    """
    if associations.empty and not epi_summary:
        return []

    variant_index = {}
    for v in variants:
        gene = str(v.get("gene", "")).strip()
        pos = v.get("position")
        alt = str(v.get("alt_aa", "")).strip()
        hgvs = str(v.get("hgvs_p", "")).strip()
        if not gene or pd.isna(pos) or not alt:
            continue
        variant_index[(gene.upper(), int(pos), alt.upper())] = v
        variant_index[hgvs.upper()] = v
        three_letter = _to_hgvs_three_letter(hgvs)
        if three_letter:
            variant_index[three_letter] = v

    matched: list[dict] = []

    # 1. CSV associations: exact protein + position + alt_aa.
    if not associations.empty:
        for _, row in associations.iterrows():
            protein = str(row.get("protein", "")).strip()
            pos = row.get("position")
            alt = str(row.get("alt_aa", "")).strip()
            if protein and not pd.isna(pos) and alt:
                key = (protein.upper(), int(pos), alt.upper())
                if key in variant_index:
                    v = variant_index[key]
                    matched.append(_phenotype_row(v, row))
                    continue

            # 2. CSV associations: motif/domain-level when position is missing.
            desc = str(row.get("genotype_description", "")).strip().lower()
            if desc and (pd.isna(pos) or not alt):
                for v in variants:
                    gene = str(v.get("gene", "")).strip()
                    domain = str(v.get("domain", "")).strip().lower()
                    if gene.upper() != protein.upper():
                        continue
                    # Extract a keyword after "motif:" or "domain:".
                    keyword = None
                    if "motif:" in desc:
                        keyword = desc.split("motif:")[1].split(",")[0].strip()
                    elif "domain:" in desc:
                        keyword = desc.split("domain:")[1].split(",")[0].strip()
                    if keyword and (keyword in domain or domain in keyword):
                        matched.append(_phenotype_row(v, row))
                        break

    # 3. DB query results: the data engine already mapped variants to phenotypes.
    epi = epi_summary or {}
    if isinstance(epi, dict):
        local_db = epi.get("local_db_results", {})
        db_phenotypes = local_db.get("phenotype_associations", {})
        # Build a literature-ref lookup from the richer variant_phenotypes records
        # (nested under layer1_variant_lookup in the DB query output).
        refs_by_id: dict[str, Any] = {}
        layer1 = epi.get("layer1_variant_lookup", {})
        for vp in layer1.get("variant_phenotypes", []):
            if isinstance(vp, dict) and vp.get("association_id"):
                refs_by_id[str(vp["association_id"])] = vp.get("literature_refs")
        if isinstance(db_phenotypes, dict):
            for hgvs_key, records in db_phenotypes.items():
                v = variant_index.get(hgvs_key) or variant_index.get(hgvs_key.upper())
                if not v:
                    # Try the three-letter form stored in the DB query.
                    for candidate in list(variant_index.values()):
                        if _to_hgvs_three_letter(candidate.get("hgvs_p")) == hgvs_key:
                            v = candidate
                            break
                if not v:
                    continue
                for rec in records:
                    if not isinstance(rec, dict):
                        continue
                    assoc_id = str(rec.get("association_id", ""))
                    refs = refs_by_id.get(assoc_id)
                    matched.append(
                        {
                            "gene": v.get("gene"),
                            "position": v.get("position"),
                            "ref_aa": v.get("ref_aa"),
                            "alt_aa": v.get("alt_aa"),
                            "hgvs_p": v.get("hgvs_p"),
                            "curated_context": v.get("_curated_context"),
                            "genotype_description": rec.get("genotype_description"),
                            "phenotype_category": rec.get("phenotype_category"),
                            "phenotype_specific": rec.get("phenotype_specific"),
                            "evidence_strength": rec.get("evidence_strength"),
                            "literature_refs": refs,
                            "record_flagged": False,
                            "source": "db_query:phenotype_associations",
                        }
                    )

    # Deduplicate by gene/position/alt/phenotype_category.
    seen = set()
    deduped = []
    for p in matched:
        key = (
            p.get("gene", ""),
            p.get("position"),
            p.get("alt_aa", ""),
            p.get("phenotype_category", ""),
            p.get("genotype_description", ""),
        )
        if key not in seen:
            seen.add(key)
            deduped.append(p)
    return deduped


def _phenotype_row(variant: dict, row: pd.Series) -> dict:
    return {
        "gene": variant.get("gene"),
        "position": variant.get("position"),
        "ref_aa": variant.get("ref_aa"),
        "alt_aa": variant.get("alt_aa"),
        "hgvs_p": variant.get("hgvs_p"),
        "curated_context": variant.get("_curated_context"),
        "genotype_description": row.get("genotype_description"),
        "phenotype_category": row.get("phenotype_category"),
        "phenotype_specific": row.get("phenotype_specific"),
        "evidence_strength": row.get("evidence_strength"),
        "literature_refs": row.get("literature_refs"),
        "record_flagged": row.get("record_flagged") == "t" or row.get("record_flagged") is True,
        "source": "PostgreSQL PGIRL database (genotype_phenotype)",
    }


def _summarise_epi(epi: dict) -> dict:
    """Flatten the epi JSON into a compact summary for the intelligence object."""
    blocks = []
    total_rows = 0
    for key, value in epi.items():
        if key.startswith("_") or key == "references":
            continue
        if isinstance(value, list):
            question = key.replace("Q", "").replace("_", " ").strip().capitalize()
            blocks.append({"question": question, "rows": value})
            total_rows += len(value)
    return {
        "questions_answered": len(blocks),
        "total_epi_rows": total_rows,
        "_raw_text_blocks": blocks,
        "references": epi.get("references", []),
    }


def _assess_risk(
    stage9: dict,
    variants: list[dict],
    matched_phenotypes: list[dict],
    epi_summary: dict,
) -> dict:
    """Apply deterministic rules to derive a risk tier and confidence."""
    rank = RISK_TIERS["routine"]
    reasons = []
    confidence = "high"

    # Concerning phenotypes observed in the sample
    concerning_records = []
    concerning_variant_keys = set()
    for p in matched_phenotypes:
        cat = (p.get("phenotype_category") or "").lower().strip()
        strength = (p.get("evidence_strength") or "").lower().strip()
        if cat in CONCERNING_PHENOTYPES and strength in ("strong", "moderate", "preliminary"):
            concerning_records.append(p)
            concerning_variant_keys.add(
                (p.get("gene"), p.get("position"), p.get("alt_aa"))
            )
            if cat in ("vaccine_escape", "drug_resistance", "virulence"):
                rank = max(rank, RISK_TIERS["high_priority"])
            elif cat in ("immune_escape", "disease_severity", "host_adaptation", "vaccine_effectiveness"):
                rank = max(rank, RISK_TIERS["investigate"])

    if concerning_records:
        n_variants = len(concerning_variant_keys)
        reasons.append(
            f"{n_variants} detected variant(s) have {len(concerning_records)} known concerning phenotype association record(s)."
        )

    # Known hotspot mutations
    hotspots = [v for v in variants if v.get("known_hotspot")]
    if hotspots:
        rank = max(rank, RISK_TIERS["monitor"])
        reasons.append(f"{len(hotspots)} variant(s) fall in known functional hotspots.")

    # Lineage not seen in the collection country before
    country_history = epi_summary.get("country_history") or []
    if country_history:
        for row in country_history[:1]:
            text = json.dumps(row).lower()
            if any(
                phrase in text
                for phrase in ["not seen", "no prior", "novel introduction", "first time"]
            ):
                rank = max(rank, RISK_TIERS["investigate"])
                reasons.append("Lineage/variant combination not previously reported in this country.")

    # Re-emergence of an old lineage
    lineage_meta = stage9.get("lineage_metadata", {})
    last_detected = lineage_meta.get("last_detected")
    collection_date = stage9.get("metadata", {}).get("collection_date")
    if last_detected and collection_date:
        try:
            last_year = int(str(last_detected)[:4])
            coll_year = int(str(collection_date)[:4])
            if coll_year - last_year > 5:
                rank = max(rank, RISK_TIERS["investigate"])
                reasons.append(
                    f"Lineage last detected in {last_detected}; sample collected {coll_year}. "
                    "Consider re-emergence, persistent infection, or undetected transmission."
                )
            elif coll_year - last_year > 0:
                rank = max(rank, RISK_TIERS["monitor"])
                reasons.append(
                    f"Lineage last detected in {last_detected}; sample collected {coll_year}. "
                    "Recurrence should be confirmed against current circulation."
                )
        except (ValueError, TypeError):
            pass

    # Variant context: if any variant is known in other countries/lineages but not in
    # this lineage/country, it may be a novel introduction.
    sample_country = stage9.get("metadata", {}).get("country", "")
    sample_lineage = stage9.get("lineage", "")
    for v in variants:
        ctx = v.get("_curated_context") or {}
        if not ctx:
            continue
        countries_seen = ctx.get("countries_seen") or []
        lineage_ids = ctx.get("lineage_ids") or []
        if isinstance(countries_seen, str):
            countries_seen = _parse_pg_array(countries_seen)
        if isinstance(lineage_ids, str):
            lineage_ids = _parse_pg_array(lineage_ids)
        if sample_country and sample_country not in countries_seen:
            rank = max(rank, RISK_TIERS["investigate"])
            reasons.append(
                f"{v.get('hgvs_p')} not previously reported in {sample_country}; possible novel introduction."
            )
        elif sample_lineage and sample_lineage not in lineage_ids:
            rank = max(rank, RISK_TIERS["monitor"])
            reasons.append(
                f"{v.get('hgvs_p')} not in curated lineage {sample_lineage} context; warrants verification."
            )

    # Travel history and epidemiological exposure
    sample_meta = stage9.get("metadata", {})
    travel = str(sample_meta.get("travel_history", "")).lower()
    travel_locations = str(sample_meta.get("travel_locations", "")).lower()
    exposure = str(sample_meta.get("suspected_exposure", "")).lower()
    if travel in ("yes", "true", "1") and any(
        loc in travel_locations for loc in ["drc", "congo", "democratic republic", "guinea", "sierra leone", "liberia"]
    ):
        rank = max(rank, RISK_TIERS["investigate"])
        reasons.append(
            f"Travel history to an EVD-affected region ({sample_meta.get('travel_locations')}) supports importation risk."
        )
    if "funeral" in exposure or "contact" in exposure:
        rank = max(rank, RISK_TIERS["investigate"])
        reasons.append(
            f"Suspected exposure ({sample_meta.get('suspected_exposure')}) indicates potential human-to-human transmission chain."
        )

    # Lower confidence when phenotypic data are sparse
    if not matched_phenotypes and not hotspots:
        confidence = "low" if len(variants) > 0 else "medium"
        reasons.append("No curated phenotype associations found for detected variants.")
    elif matched_phenotypes and not any(
        p.get("evidence_strength") == "strong" for p in matched_phenotypes
    ):
        confidence = "medium"

    tier_name = [k for k, v in RISK_TIERS.items() if v == rank][0]

    limitations = (
        "The risk assessment is contingent on the completeness and accuracy of curated reference data. "
        "In vivo or epidemiological verification may be necessary for confirming phenotype associations based on this genomic data."
    )

    return {
        "risk_tier": tier_name,
        "rank": rank,
        "rationale": " ".join(reasons) if reasons else "No concerning features identified.",
        "confidence": confidence,
        "limitations": limitations,
        "concerning_phenotype_count": len(concerning_variant_keys),
    }


def _build_intelligence_object(
    stage9: dict,
    variants: list[dict],
    matched_phenotypes: list[dict],
    epi_summary: dict,
    risk: dict,
    all_results: dict[str, Any],
    threat: Any,
    tree_input: Any,
    data_quality: list[str],
    bio: Optional[dict] = None,
) -> dict:
    """Assemble the structured object that feeds report generation and R figures."""
    sample_meta = stage9.get("metadata", {})
    genome_quality = stage9.get("genome_quality", {})
    closest_ref = stage9.get("closest_reference", {})
    lineage_meta = stage9.get("lineage_metadata", {})
    qc = bio.get("stage0_quality_control", {}) if bio else {}
    stage2 = bio.get("stage2_reference_context", {}) if bio else {}

    sample = {
        "sample_id": stage9.get("sample_id"),
        "pathogen": stage9.get("pathogen"),
        "species": stage9.get("species"),
        "species_id": stage9.get("species_id"),
        "lineage": stage9.get("lineage"),
        "clade": stage9.get("clade"),
        "collection_date": sample_meta.get("collection_date"),
        "country": sample_meta.get("country"),
        "admin1": sample_meta.get("admin1"),
        "admin2": sample_meta.get("admin2"),
        "host": sample_meta.get("host"),
        "host_species": sample_meta.get("host_species"),
        "sample_type": sample_meta.get("sample_type"),
        "vaccination_status": sample_meta.get("vaccination_status"),
        "travel_history": sample_meta.get("travel_history"),
        "travel_locations": sample_meta.get("travel_locations"),
        "outcome": sample_meta.get("outcome"),
        "suspected_exposure": sample_meta.get("suspected_exposure"),
        "epi_link_id": sample_meta.get("epi_link_id"),
        "genome_completeness_pct": genome_quality.get("completeness_pct"),
        "mean_depth": genome_quality.get("depth"),
        "quality_flag": genome_quality.get("flag"),
        "total_aa_variants": len(variants),
        "closest_reference": closest_ref,
        "blast_top_hits": stage2.get("blast_top_hits", []),
        "genome_qc": {
            "genome_length_bp": qc.get("genome_length_bp"),
            "expected_length_bp": qc.get("expected_length_bp"),
            "genome_completeness_pct": qc.get("genome_completeness_pct"),
            "n_content_pct": qc.get("n_content_pct"),
            "mean_depth": qc.get("mean_depth"),
            "quality_flag": qc.get("quality_flag"),
            "missing_regions": qc.get("missing_regions", []),
        },
        "lineage_metadata": lineage_meta,
        "summary": (
            f"Sample {stage9.get('sample_id')} is {stage9.get('species')} lineage "
            f"{stage9.get('lineage')} collected from {sample_meta.get('host')} in "
            f"{sample_meta.get('country')} on {sample_meta.get('collection_date')}."
        ),
    }

    analyses = {key: asdict(result) for key, result in all_results.items()}
    analyses["evidence_weighted_threat"] = asdict(threat)

    # Extract negative findings (absence of known concerning variants/motifs) as a
    # concise surveillance-oriented list in the top-level intelligence object.
    negative_findings = []
    gs = all_results.get("genomic_significance")
    if gs and hasattr(gs, "findings"):
        for f in gs.findings:
            fd = f if isinstance(f, dict) else asdict(f)
            if fd.get("finding_type") == "negative_finding":
                negative_findings.append(fd)

    tree_tips = []
    if tree_input:
        tree_tips = [asdict(t) for t in tree_input.tips]

    return {
        "sample": sample,
        "variants": variants,
        "matched_phenotypes": matched_phenotypes,
        "epi_summary": epi_summary,
        "risk_assessment": risk,
        "analyses": analyses,
        "negative_findings": negative_findings,
        "tree_available": tree_input.has_tree if tree_input else False,
        "tree_file": str(tree_input.file_path) if tree_input and tree_input.file_path else None,
        "tree_tips": tree_tips,
        "references": epi_summary.get("references", []),
        "data_quality_warnings": data_quality,
    }


def _build_epi_text(epi_summary: dict) -> str:
    """Create a concise text block for the brief template from epi findings."""
    parts = []
    for block in epi_summary.get("_raw_text_blocks", []):
        q = block["question"]
        rows = block["rows"]
        if not rows:
            continue
        parts.append(f"### {q}")
        for row in rows:
            parts.append(
                "; ".join(f"{k}: {v}" for k, v in row.items() if not k.startswith("_"))
            )
    return "\n".join(parts)


def _check_data_quality(
    stage9: dict,
    variants: list[dict],
    matched_phenotypes: list[dict],
    tree_input: Any,
) -> list[str]:
    """Return human-readable warnings about data gaps that affect confidence."""
    warnings: list[str] = []
    sample_id = stage9.get("sample_id")
    if not sample_id:
        warnings.append("stage9_normalised_output is missing sample_id.")

    lineage = stage9.get("lineage")
    lineage_meta = stage9.get("lineage_metadata", {})
    if not lineage_meta.get("last_detected"):
        warnings.append(f"Lineage metadata for {lineage} has no last_detected date; re-emergence inference uses genome metadata only.")
    else:
        collection_date = stage9.get("metadata", {}).get("collection_date")
        if collection_date:
            try:
                last_year = int(str(lineage_meta.get("last_detected"))[:4])
                coll_year = int(str(collection_date)[:4])
                gap = coll_year - last_year
                if gap > 5:
                    warnings.append(
                        f"Lineage {lineage} was last detected in {lineage_meta.get('last_detected')} "
                        f"and the sample was collected in {collection_date} (gap: {gap} years). "
                        "Consider laboratory contamination, archival/persistent infection, or cryptic transmission."
                    )
            except (ValueError, TypeError):
                pass
    if not lineage_meta.get("countries_reported"):
        warnings.append(f"Lineage metadata for {lineage} has no reported countries; geographic novelty is inferred from genome metadata only.")

    if not variants:
        warnings.append("No amino-acid variants were extracted from the bioinformatics output.")

    if not matched_phenotypes:
        warnings.append("No phenotype associations were matched; risk assessment relies on lineage/history alone.")

    if tree_input and tree_input.has_tree:
        annotated = sum(1 for t in tree_input.tips if t.country)
        total = len(tree_input.tips)
        if annotated < total / 2:
            warnings.append(
                f"Only {annotated}/{total} tree tips could be annotated with country/date metadata; "
                "phylogeographic and transmission-network figures will be sparse."
            )
    else:
        warnings.append("No phylogenetic tree available; tree-aware analyses use metadata fallback only.")

    return warnings


def _flatten_metrics(metrics: dict, prefix: str = "") -> list[dict]:
    """Flatten a nested metrics dict into a list of {metric, value} rows."""
    rows = []
    for key, value in metrics.items():
        full_key = f"{prefix}.{key}" if prefix else key
        if isinstance(value, dict):
            rows.extend(_flatten_metrics(value, full_key))
        elif isinstance(value, list):
            if value and isinstance(value[0], dict):
                list_str = "; ".join(
                    ", ".join(f"{k}={v}" for k, v in sorted(item.items()) if not k.startswith("_"))
                    for item in value
                )
            else:
                list_str = "; ".join(str(v) for v in value)
            rows.append({"metric": full_key, "value": list_str})
        else:
            rows.append({"metric": full_key, "value": value if value is not None else ""})
    return rows


def _build_evidence_package(
    variant_findings: list[dict],
    matched_phenotypes: list[dict],
    lineage_row,
    associations_df: pd.DataFrame,
    genome_metadata_df: pd.DataFrame,
    epi: dict,
    output_dir: Path,
) -> dict:
    """Run the evidence-integration layer: harmonize evidence into unified
    objects, run all cross-evidence statistical analyzers, render
    visualizations, and write the resulting evidence package to disk.

    This is purely a scientific evidence synthesis step -- no risk tiers or
    public-health conclusions are produced here.
    """
    lineage_id = str(lineage_row.get("lineage_id")) if lineage_row is not None else None

    evidence_objects = build_evidence_objects(
        variants=variant_findings,
        matched_phenotypes=matched_phenotypes,
        lineage_row=lineage_row,
        epi_output=epi or {},
    )

    outbreaks = (epi or {}).get("outbreaks") or []
    surveillance = (epi or {}).get("surveillance") or []

    cross_evidence_results = {
        "mutation_cooccurrence": MutationCooccurrenceAnalyzer().analyze(variant_findings).to_dict(),
        "lineage_phenotype_association": LineagePhenotypeAssociationAnalyzer(associations_df).analyze(lineage_id).to_dict(),
        "temporal_trend": TemporalTrendAnalyzer().analyze(genome_metadata_df, outbreaks).to_dict(),
        "geographic_distribution": GeographicDistributionAnalyzer().analyze(genome_metadata_df, outbreaks, surveillance).to_dict(),
        "mutation_persistence": MutationPersistenceAnalyzer().analyze(variant_findings).to_dict(),
        "intervention_association": InterventionAssociationAnalyzer().analyze(matched_phenotypes).to_dict(),
        "confidence_scoring": ConfidenceScoringAnalyzer().analyze(evidence_objects).to_dict(),
    }

    evidence_dir = output_dir / "evidence_package"
    evidence_dir.mkdir(parents=True, exist_ok=True)

    network_path = save_evidence_network(evidence_objects, str(evidence_dir))

    geo_findings = cross_evidence_results["geographic_distribution"]["findings"]
    genome_country_counts = {
        f["value"]["country"]: f["value"]["count"]
        for f in geo_findings if f["metric"] == "genome_country_distribution"
    }
    geo_chart_path = save_geographic_distribution_chart(genome_country_counts, str(evidence_dir))

    trend_findings = cross_evidence_results["temporal_trend"]["findings"]
    trend_pairs = []
    for f in trend_findings:
        if f["metric"] == "genome_collection_trend":
            trend_pairs = [
                (int(e.split(":")[0].split("=")[1]), int(e.split(":")[1].split("=")[1]))
                for e in f.get("supporting_evidence", [])
            ]
            break
    trend_chart_path = save_temporal_trend_chart(trend_pairs, str(evidence_dir), ylabel="genomes/year")

    package = {
        "evidence_objects": [obj.to_dict() for obj in evidence_objects],
        "cross_evidence_analysis": cross_evidence_results,
        "visualizations": {
            "evidence_network": network_path,
            "geographic_distribution_chart": geo_chart_path,
            "temporal_trend_chart": trend_chart_path,
        },
    }

    package_path = evidence_dir / "evidence_package.json"
    with open(package_path, "w") as f:
        json.dump(_clean_for_json(package), f, indent=2, default=str)
    log.info("Wrote evidence package to %s", package_path)

    package["_path"] = str(package_path)
    return package


def _write_analysis_outputs(
    output_dir: Path,
    all_results: dict[str, Any],
    threat: Any,
    risk: dict,
    sample: dict,
    data_quality: list[str],
) -> Path:
    """Write each analytical result as plain CSV outputs."""
    outputs_dir = output_dir / "analysis_outputs"
    outputs_dir.mkdir(parents=True, exist_ok=True)

    findings_cols = [
        "analysis",
        "assertion",
        "source",
        "confidence",
        "biological_relevance",
        "epidemiological_relevance",
        "finding_type",
        "record_flagged",
        "supporting_refs",
    ]
    all_findings_rows: list[dict] = []

    results_to_export = {**all_results, "evidence_weighted_threat": threat}

    for idx, (name, result) in enumerate(results_to_export.items(), start=1):
        result_dict = asdict(result) if not isinstance(result, dict) else result
        title = result_dict.get("title", name)
        findings = result_dict.get("findings", [])
        metrics = result_dict.get("metrics", {})

        finding_rows = []
        for f in findings:
            fdict = asdict(f) if not isinstance(f, dict) else f
            row = {
                "analysis": title,
                "assertion": fdict.get("assertion", ""),
                "source": fdict.get("source", ""),
                "confidence": fdict.get("confidence", ""),
                "biological_relevance": fdict.get("biological_relevance", ""),
                "epidemiological_relevance": fdict.get("epidemiological_relevance", ""),
                "finding_type": fdict.get("finding_type", ""),
                "record_flagged": fdict.get("record_flagged", False),
                "supporting_refs": "; ".join(str(x) for x in fdict.get("supporting_refs", [])),
            }
            finding_rows.append(row)
            all_findings_rows.append(row)

        prefix = f"{idx:02d}_{name}"
        pd.DataFrame(finding_rows if finding_rows else [], columns=findings_cols).to_csv(
            outputs_dir / f"{prefix}_output.csv", index=False
        )

        metric_rows = _flatten_metrics(metrics)
        if metric_rows:
            pd.DataFrame(metric_rows).to_csv(
                outputs_dir / f"{prefix}_metrics.csv", index=False
            )
        else:
            (outputs_dir / f"{prefix}_metrics.csv").write_text("metric,value\n")

    if all_findings_rows:
        pd.DataFrame(all_findings_rows, columns=findings_cols).to_csv(
            outputs_dir / "combined_analysis_output.csv", index=False
        )
    else:
        pd.DataFrame(columns=findings_cols).to_csv(
            outputs_dir / "combined_analysis_output.csv", index=False
        )

    risk_rows = [{"key": k, "value": v} for k, v in risk.items()]
    pd.DataFrame(risk_rows).to_csv(outputs_dir / "risk_assessment.csv", index=False)

    sample_rows = [{"key": k, "value": v} for k, v in sample.items()]
    pd.DataFrame(sample_rows).to_csv(outputs_dir / "sample_metadata.csv", index=False)

    # Always write data quality file so stale warnings from previous runs are removed.
    pd.DataFrame({"warning": data_quality}).to_csv(
        outputs_dir / "data_quality_warnings.csv", index=False
    )

    log.info("Analysis outputs written to %s", outputs_dir)
    return outputs_dir


def _run_decision_figures(
    output_dir: Path,
    json_path: Path,
) -> Optional[Path]:
    """Run the R decision-figures script to generate visual outputs."""
    if shutil.which("Rscript") is None:
        log.warning("Rscript not found; skipping decision figures")
        return None

    r_script = (
        Path(__file__).resolve().parent.parent
        / "figures"
        / "decision_figures.R"
    )
    if not r_script.exists():
        log.warning("R script not found: %s; skipping figures", r_script)
        return None

    figures_dir = output_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "Rscript",
        str(r_script),
        str(json_path),
        str(figures_dir),
    ]
    log.info("Running decision figures: %s", " ".join(cmd))
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        log.info("Figures written to %s", figures_dir)
        return figures_dir
    except subprocess.CalledProcessError as e:
        log.warning("R decision figures failed: %s\n%s", e, e.stderr)
        return None


class IntelligencePipeline:
    """Run all assessment analyzers and produce the intelligence object and outputs."""

    def __init__(
        self,
        epi_output_path: str,
        bio_output_path: str,
        associations_csv_path: str = "database/exports/genotype_phenotype.csv",
        protein_variants_csv_path: str = "database/exports/protein_variants.csv",
        lineages_csv_path: str = "database/exports/lineages.csv",
        genome_metadata_csv_path: str = "database/exports/genome_metadata.csv",
        tree_file_path: Optional[str] = None,
        db_url: Optional[str] = None,
    ):
        self.epi_output_path = Path(epi_output_path)
        self.bio_output_path = Path(bio_output_path)
        self.tree_file_path = tree_file_path
        self.associations_csv_path = Path(associations_csv_path)
        self.protein_variants_csv_path = Path(protein_variants_csv_path)
        self.lineages_csv_path = Path(lineages_csv_path)
        self.genome_metadata_csv_path = Path(genome_metadata_csv_path)

        self.engine = GenomicIntelligenceEngine(
            associations_csv_path=str(self.associations_csv_path),
            protein_variants_csv_path=str(self.protein_variants_csv_path),
            lineages_csv_path=str(self.lineages_csv_path),
            genome_metadata_csv_path=str(self.genome_metadata_csv_path),
            epi_output_path=str(self.epi_output_path),
            db_url=db_url,
        )

    def run(self, output_dir: str = "output") -> dict:
        """Execute the pipeline and write all outputs."""
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        epi = _load_json(self.epi_output_path)
        bio = _load_json(self.bio_output_path)
        stage9 = bio.get("stage9_normalised_output", bio)

        # The bioinformatics stage now populates `clade` (e.g. "Ebov-1976") from
        # Nextclade but may leave `lineage` empty. The curated PGIRL lineage table
        # keys its records by clade-equivalent ids (matched via _match_lineage
        # aliases), so treat the clade as the lineage identifier when no explicit
        # lineage was assigned. Without this, every lineage-based analyzer
        # (molecular epidemiology, surveillance, lineage behaviour, etc.) is
        # starved and returns empty findings.
        if not str(stage9.get("lineage") or "").strip() and str(stage9.get("clade") or "").strip():
            stage9["lineage"] = stage9["clade"]

        # Set pathogen/lineage context so the engine can query the PostgreSQL database directly.
        self.engine.set_pathogen_context(
            pathogen_id=stage9.get("pathogen"),
            species_id=stage9.get("species_id"),
        )
        self.engine.set_lineage_context(stage9.get("lineage"))

        # Expose stage6 tree metadata and phylogenetic placement to tree-aware analyzers
        stage6 = bio.get("stage6_phylogenetic_tree")
        if stage6:
            stage9.setdefault("stage6_phylogenetic_tree", stage6)
        stage5 = bio.get("stage5_lineage_clade", {})
        placement = stage5.get("phylogenetic_placement", {})
        if placement:
            stage9.setdefault("phylogenetic_placement", placement)

        # Variant extraction happens first so the engine can issue narrow DB queries.
        variant_findings = _extract_variants(stage9, bio)
        self.engine.set_detected_variants(variant_findings)

        # Curated lineage metadata (now loaded with targeted DB context)
        curated_lineage_meta = self.engine.get_lineage_metadata(stage9.get("lineage"))
        if curated_lineage_meta:
            stage9["lineage_metadata"] = curated_lineage_meta
        else:
            stage9.setdefault(
                "lineage_metadata",
                bio.get("stage5_lineage_clade", {})
                .get("lineage_determination", {})
                .get("lineage_metadata", {}),
            )

        # Load the remaining targeted reference tables.
        associations, protein_variants_df, _lineages_df, _genome_meta_df = self.engine.get_reference_data()
        variant_findings = _enrich_variants(variant_findings, protein_variants_df)
        epi_summary = _summarise_epi(epi)
        matched_phenotypes = _match_phenotypes(variant_findings, associations, epi)

        # Evidence integration: harmonize + cross-evidence statistical analysis.
        # Purely evidentiary/statistical output -- no risk tiers or conclusions.
        lineage_row = _match_lineage(stage9.get("lineage", ""), _lineages_df)
        evidence_package = _build_evidence_package(
            variant_findings=variant_findings,
            matched_phenotypes=matched_phenotypes,
            lineage_row=lineage_row,
            associations_df=associations,
            genome_metadata_df=_genome_meta_df,
            epi=epi,
            output_dir=output_dir,
        )

        risk = _assess_risk(stage9, variant_findings, matched_phenotypes, epi_summary)

        # Core analyses
        genomic_significance = self.engine.analyze_genomic_significance(
            stage9, variant_findings, matched_phenotypes
        )
        molecular_epi = self.engine.analyze_molecular_epidemiology(stage9)

        # Tree-aware analyses
        tree_input = self.engine.load_tree_input(stage9, tree_path=self.tree_file_path)
        phylogeo = self.engine.analyze_phylogeography(stage9, tree_input)
        genetic_relatedness = self.engine.analyze_genetic_relatedness(stage9, tree_input)
        molecular_clock = self.engine.analyze_molecular_clock(stage9, tree_input)

        # Extended analyses
        epidemic_dynamics = self.engine.analyze_epidemic_dynamics(stage9)
        genomic_signal = self.engine.analyze_genomic_signal(stage9, variant_findings, matched_phenotypes)
        surveillance = self.engine.analyze_surveillance(stage9)
        lineage_behavior = self.engine.analyze_lineage_behavior(stage9)

        all_results = {
            "genomic_significance": genomic_significance,
            "molecular_epidemiology": molecular_epi,
            "phylogeographic_analysis": phylogeo,
            "genetic_relatedness_analysis": genetic_relatedness,
            "molecular_clock_analysis": molecular_clock,
            "epidemic_dynamics": epidemic_dynamics,
            "genomic_signal": genomic_signal,
            "surveillance_metrics": surveillance,
            "lineage_behavior": lineage_behavior,
        }

        # Integrated inference depends on all prior results
        integrated = self.engine.analyze_integrated_inference(stage9, tree_input, all_results)
        all_results["integrated_inference"] = integrated

        # Cross-cutting analyses that consume all prior results
        evidence_consistency = self.engine.analyze_evidence_consistency(stage9, all_results)
        all_results["evidence_consistency"] = evidence_consistency

        comparative_outbreak = self.engine.analyze_comparative_outbreak(stage9)
        all_results["comparative_outbreak_analysis"] = comparative_outbreak

        knowledge_gaps = self.engine.analyze_knowledge_gaps(stage9, all_results)
        all_results["knowledge_gaps"] = knowledge_gaps

        # Final evidence-weighted threat
        threat = self.engine.analyze_evidence_weighted_threat(stage9, all_results, risk)
        all_results["evidence_weighted_threat"] = threat

        data_quality = _check_data_quality(stage9, variant_findings, matched_phenotypes, tree_input)
        if data_quality:
            for w in data_quality:
                log.warning("Data quality: %s", w)

        intelligence_object = _build_intelligence_object(
            stage9,
            variant_findings,
            matched_phenotypes,
            epi_summary,
            risk,
            all_results,
            threat,
            tree_input,
            data_quality,
            bio,
        )

        json_path = output_dir / "intelligence_object.json"
        with open(json_path, "w") as f:
            json.dump(_clean_for_json(intelligence_object), f, indent=2, default=str)
        log.info("Wrote intelligence object to %s", json_path)

        # Generate decision-oriented figures.
        figures_dir = _run_decision_figures(output_dir, json_path)

        # Export the analytical results as CSV outputs.
        outputs_dir = _write_analysis_outputs(
            output_dir,
            all_results,
            threat,
            risk,
            intelligence_object.get("sample", {}),
            data_quality,
        )

        return {
            "intelligence_object": intelligence_object,
            "json_path": str(json_path),
            "outputs_dir": str(outputs_dir),
            "figures_dir": str(figures_dir) if figures_dir else None,
            "evidence_package_path": evidence_package.get("_path"),
        }


def main():
    parser = argparse.ArgumentParser(
        description="Run the assessment-logic pipeline and write intelligence_object.json, CSVs and figures."
    )
    parser.add_argument(
        "--epi-output",
        default="output/data_query/epi_output_inspection.json",
        help="Path to the epidemiological output JSON.",
    )
    parser.add_argument(
        "--bio-output",
        default="output/bioinformatics/EBOV-UGA-2027-001/bio_output.json",
        help="Path to the bioinformatics pipeline output JSON.",
    )
    parser.add_argument(
        "--associations",
        default="database/exports/genotype_phenotype.csv",
        help="Path to the genotype-phenotype associations CSV.",
    )
    parser.add_argument(
        "--variants",
        default="database/exports/protein_variants.csv",
        help="Path to the protein variant frequency CSV.",
    )
    parser.add_argument(
        "--lineages",
        default="database/exports/lineages.csv",
        help="Path to the lineage metadata CSV.",
    )
    parser.add_argument(
        "--genome-metadata",
        default="database/exports/genome_metadata.csv",
        help="Path to the genome metadata CSV.",
    )
    parser.add_argument(
        "--tree-file",
        default="output/bioinformatics/EBOV-UGA-2027-001/tree.nwk",
        help="Path to the Newick tree file (4th engine input).",
    )
    parser.add_argument(
        "--output-dir",
        default="output/evidence_integration/EBOV-UGA-2027-001",
        help="Directory to write the evidence integration outputs (intelligence object, evidence package, figures, analysis outputs).",
    )
    parser.add_argument(
        "--db-url",
        default=None,
        help="PostgreSQL connection URL (defaults to PGIRL_DB_URL env var, then config.DB_URL).",
    )
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    pipeline = IntelligencePipeline(
        epi_output_path=args.epi_output,
        bio_output_path=args.bio_output,
        associations_csv_path=args.associations,
        protein_variants_csv_path=args.variants,
        lineages_csv_path=args.lineages,
        genome_metadata_csv_path=args.genome_metadata,
        tree_file_path=args.tree_file,
        db_url=args.db_url,
    )
    result = pipeline.run(output_dir=args.output_dir)
    print(f"Intelligence object written to: {result['json_path']}")
    print(f"Analysis outputs written to: {result['outputs_dir']}")
    if result["figures_dir"]:
        print(f"Figures written to: {result['figures_dir']}")
    if result.get("evidence_package_path"):
        print(f"Evidence package written to: {result['evidence_package_path']}")


if __name__ == "__main__":
    main()
