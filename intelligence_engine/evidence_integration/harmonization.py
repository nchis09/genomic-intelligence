"""
evidence_integration/harmonization.py — Evidence harmonization.

Links genomic variants, lineage information, protein annotations, phenotype
associations, historical outbreaks, and epidemiological context into unified
``EvidenceObject`` records, one per detected variant (plus one lineage-level
object for the sample as a whole).

This module deliberately produces *evidence*, not conclusions: no risk tiers,
no public-health implications. It is the input layer for
``evidence_integration.cross_evidence``.
"""

from __future__ import annotations

import logging
from dataclasses import asdict, dataclass, field
from typing import Any, Optional

import pandas as pd

log = logging.getLogger(__name__)


@dataclass
class EvidenceObject:
    """Unified evidence record for one detected variant (or the lineage as a
    whole, when ``variant`` is None), linking every evidence source that
    mentions it."""

    key: str  # e.g. "GP:A82V" or "lineage:EBOV-Makona"
    level: str  # "variant" | "lineage"
    variant: Optional[dict] = None
    lineage: Optional[dict] = None
    phenotype_associations: list[dict] = field(default_factory=list)
    historical_outbreaks: list[dict] = field(default_factory=list)
    molecular_epidemiology: list[dict] = field(default_factory=list)
    epidemiological_context: dict = field(default_factory=dict)
    interventions: list[dict] = field(default_factory=list)
    genomic_links: list[dict] = field(default_factory=list)
    sources: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _lineage_matches(lineage_field: Any, lineage_id: Optional[str], lineage_name: Optional[str]) -> bool:
    """Return True if a free-text lineage/pathogen field plausibly refers to
    the sample's lineage."""
    if not lineage_field or not (lineage_id or lineage_name):
        return False
    text = str(lineage_field).lower()
    return (lineage_id and lineage_id.lower() in text) or (lineage_name and lineage_name.lower() in text)


def build_evidence_objects(
    variants: list[dict],
    matched_phenotypes: list[dict],
    lineage_row: Optional[pd.Series],
    epi_output: dict,
) -> list[EvidenceObject]:
    """Build one ``EvidenceObject`` per detected variant, plus one lineage-level
    object, harmonizing all available evidence sources.

    Args:
        variants: Detected variants (as enriched by ``_enrich_variants``),
                   each with ``_curated_context`` when available.
        matched_phenotypes: Output of ``_match_phenotypes``.
        lineage_row: The matched row from the curated lineages table, if any.
        epi_output: The normalized entity-based epi object (outbreaks,
                    molecular_epidemiology, demographics, transmission,
                    interventions, surveillance, genomic_links, ...).
    """
    epi_output = epi_output or {}
    lineage_id = str(lineage_row.get("lineage_id")) if lineage_row is not None else None
    lineage_name = str(lineage_row.get("lineage_name")) if lineage_row is not None else None

    outbreaks = epi_output.get("outbreaks") or []
    molecular_epi = epi_output.get("molecular_epidemiology") or []
    interventions = (epi_output.get("interventions") or []) + \
        (epi_output.get("vaccines") or []) + (epi_output.get("therapeutics") or [])
    genomic_links = epi_output.get("genomic_links") or []
    transmission = epi_output.get("transmission") or {}
    demographics = epi_output.get("demographics") or []
    surveillance = epi_output.get("surveillance") or []

    objects: list[EvidenceObject] = []

    # ── Variant-level evidence objects ──
    phenotype_by_variant: dict[tuple, list[dict]] = {}
    for p in matched_phenotypes:
        key = (str(p.get("gene", "")).upper(), p.get("position"), str(p.get("alt_aa", "")).upper())
        phenotype_by_variant.setdefault(key, []).append(p)

    for v in variants:
        gene = str(v.get("gene", "")).strip()
        pos = v.get("position")
        alt = str(v.get("alt_aa", "")).strip()
        hgvs = v.get("hgvs_p") or f"{gene}:{v.get('ref_aa', '')}{pos}{alt}"
        key = (gene.upper(), pos, alt.upper())
        phenotypes = phenotype_by_variant.get(key, [])

        matching_links = [
            gl for gl in genomic_links
            if gl.get("mutations") and (hgvs in str(gl.get("mutations")) or gene.upper() in str(gl.get("mutations")).upper())
        ]

        sources = ["bioinformatics_call"]
        if v.get("_curated_context"):
            sources.append("curated_protein_variants")
        if phenotypes:
            sources.append("genotype_phenotype")
        if matching_links:
            sources.append("epi_genomic_links")

        objects.append(
            EvidenceObject(
                key=str(hgvs),
                level="variant",
                variant={
                    "gene": gene, "position": pos, "ref_aa": v.get("ref_aa"),
                    "alt_aa": alt, "hgvs_p": hgvs, "domain": v.get("domain"),
                    "curated_context": v.get("_curated_context"),
                },
                lineage={"lineage_id": lineage_id, "lineage_name": lineage_name} if lineage_id else None,
                phenotype_associations=phenotypes,
                genomic_links=matching_links,
                sources=sources,
            )
        )

    # ── One lineage-level evidence object aggregating epi context ──
    if lineage_id or lineage_name:
        matching_outbreaks = [
            ob for ob in outbreaks
            if _lineage_matches(ob.get("lineage"), lineage_id, lineage_name)
        ] or outbreaks  # fall back to all fetched outbreaks for this pathogen if none tag the lineage explicitly

        matching_mol_epi = [
            me for me in molecular_epi
            if _lineage_matches(me.get("lineage"), lineage_id, lineage_name)
        ]

        objects.append(
            EvidenceObject(
                key=f"lineage:{lineage_id or lineage_name}",
                level="lineage",
                lineage=(lineage_row.to_dict() if lineage_row is not None else {"lineage_name": lineage_name}),
                historical_outbreaks=matching_outbreaks,
                molecular_epidemiology=matching_mol_epi,
                epidemiological_context={
                    "transmission": transmission,
                    "demographics": demographics,
                    "surveillance": surveillance,
                },
                interventions=interventions,
                sources=["curated_lineages", "epi_query_engine"],
            )
        )

    return objects
