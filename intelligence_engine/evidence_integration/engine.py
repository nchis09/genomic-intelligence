"""
evidence_integration/engine.py — Core evidence-integration engine.

This module transforms the combined epidemiological, bioinformatics, and curated
reference data into a structured, evidence-based intelligence object.

Phase 1 implements the GenomicSignificanceAnalyzer. Additional analyzers for
molecular epidemiology, comparative genomics, public-health context, knowledge
gaps, evidence-weighted significance, and surveillance prioritization will be
added incrementally.
"""

import csv
import json
import logging
import re
from collections import Counter
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Optional

import pandas as pd

from intelligence_engine.evidence_integration.tree.tree_input import TreeInput, load_tree_input

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Structured intelligence object schema
# ---------------------------------------------------------------------------

@dataclass
class EvidenceRecord:
    """A single, auditable evidence statement from an analyzer."""

    assertion: str
    source: str
    confidence: str
    biological_relevance: str = ""
    epidemiological_relevance: str = ""
    supporting_refs: list[str] = field(default_factory=list)
    # NOTE: this evidence-integration layer intentionally does not carry a
    # public-health-implication field. Interpretive/risk conclusions are the
    # responsibility of the downstream Genomic Intelligence Engine stage.
    record_flagged: bool = False
    finding_type: str = ""


@dataclass
class AnalysisResult:
    """Result of one analyzer."""

    title: str
    findings: list[EvidenceRecord] = field(default_factory=list)
    gaps: list[str] = field(default_factory=list)
    summary: str = ""
    metrics: dict = field(default_factory=dict)


@dataclass
class GenomicIntelligenceObject:
    """Container for all structured intelligence findings."""

    sample_id: str
    sample: dict = field(default_factory=dict)
    variants: list[dict] = field(default_factory=list)
    epi_summary: dict = field(default_factory=dict)
    analyses: dict[str, AnalysisResult] = field(default_factory=dict)
    evidence_records: list[EvidenceRecord] = field(default_factory=list)
    knowledge_gaps: list[str] = field(default_factory=list)
    evidence_weighted_significance: dict = field(default_factory=dict)
    surveillance_priorities: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-serializable dict."""
        return asdict(self)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _parse_pg_array(text: Any) -> list[str]:
    """Parse a Postgres array literal like '{a,b,"c d"}' into a list of strings."""
    if isinstance(text, list):
        return [str(x).strip() for x in text if x is not None]
    if not text or not isinstance(text, str):
        return []
    text = text.strip()
    if text.startswith("{") and text.endswith("}"):
        text = text[1:-1]
    if not text:
        return []
    reader = csv.reader([text], delimiter=",", quotechar='"')
    return [x.strip() for x in next(reader) if x.strip()]


def _safe_str(value: Any) -> str:
    """Return a clean string or empty string for missing/NaN values."""
    if value is None:
        return ""
    if isinstance(value, float) and pd.isna(value):
        return ""
    return str(value).strip()


def _parse_date(value: Any) -> str:
    """Return a normalised YYYY, YYYY-MM or YYYY-MM-DD string."""
    s = _safe_str(value)
    if not s:
        return ""
    m = re.match(r"(\d{4})(?:-\d{2}(?:-\d{2})?)?", s)
    return m.group(0) if m else s


def _to_date(value: Any) -> Optional[pd.Timestamp]:
    """Return a parsed pandas Timestamp or None."""
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return None
    try:
        return pd.to_datetime(value)
    except (ValueError, TypeError):
        return None


def _normalise_lineage_label(value: str) -> str:
    """Normalise lineage labels for comparison across source conventions."""
    value = re.sub(r"[^a-z0-9]+", "-", str(value).strip().lower()).strip("-")
    value = re.sub(r"^(ebov|sudv|bundibugyo|taiforest|reston|bombali)-", "", value)
    return value


def _match_lineage(lineage: str, lineages_df: pd.DataFrame) -> Optional[pd.Series]:
    """Find a lineage row by id, name, or known alias."""
    if lineages_df.empty or not lineage:
        return None

    lineage_norm = lineage.strip().lower()

    for col in ["lineage_id", "lineage_name"]:
        if col in lineages_df.columns:
            values = lineages_df[col].fillna("").astype(str).str.strip().str.lower()
            rows = lineages_df[values == lineage_norm]
            if not rows.empty:
                return rows.iloc[0]

    lineage_label = _normalise_lineage_label(lineage)
    if "known_aliases" in lineages_df.columns:
        for _, row in lineages_df.iterrows():
            aliases = _parse_pg_array(str(row.get("known_aliases", "")))
            candidates = aliases + [str(row.get("lineage_name", "")), str(row.get("lineage_id", ""))]
            if lineage_label in {_normalise_lineage_label(candidate) for candidate in candidates}:
                return row

    return None


# ---------------------------------------------------------------------------
# Analyzers
# ---------------------------------------------------------------------------

class GenomicSignificanceAnalyzer:
    """Evaluate the genomic significance of the detected pathogen/variants."""

    def __init__(
        self,
        associations_df: pd.DataFrame,
        protein_variants_df: pd.DataFrame,
        lineages_df: pd.DataFrame,
        genome_metadata_df: pd.DataFrame,
    ):
        self.associations = associations_df
        self.protein_variants = protein_variants_df
        self.lineages = lineages_df
        self.genome_metadata = genome_metadata_df

        # Normalise columns
        for col in ["protein", "position", "ref_aa", "alt_aa"]:
            if col in self.associations.columns:
                self.associations[col] = self.associations[col].astype(str).str.strip()
            if col in self.protein_variants.columns:
                self.protein_variants[col] = (
                    self.protein_variants[col].astype(str).str.strip()
                )
        for df in [self.associations, self.protein_variants]:
            if "position" in df.columns:
                df["position"] = pd.to_numeric(df["position"], errors="coerce")

    def analyze(
        self,
        stage9: dict,
        variants: list[dict],
        matched_phenotypes: list[dict],
    ) -> AnalysisResult:
        """Run all genomic-significance analyses and return an AnalysisResult."""
        findings = []
        gaps = []

        lineage = stage9.get("lineage")
        species = stage9.get("species")
        sample_id = stage9.get("sample_id")

        # --- Lineage commonness / unusualness --------------------------------
        lineage_row = _match_lineage(lineage, self.lineages) if lineage else None
        if lineage_row is not None:
            lineage_id = lineage_row.get("lineage_id")
            first_detected = _safe_str(lineage_row.get("first_detected"))
            last_detected = _safe_str(lineage_row.get("last_detected"))

            # Count genomes in curated genome_metadata and derive dates if missing
            genome_count = 0
            meta_earliest = meta_latest = None
            if not self.genome_metadata.empty and "lineage_id" in self.genome_metadata.columns:
                meta_mask = self.genome_metadata["lineage_id"].astype(str).str.strip().str.lower() == str(lineage_id).strip().lower()
                meta_rows = self.genome_metadata[meta_mask]
                genome_count = len(meta_rows)
                if genome_count and "collection_date" in meta_rows.columns:
                    dates = pd.to_datetime(meta_rows["collection_date"], errors="coerce")
                    meta_earliest = dates.min() if not dates.empty else None
                    meta_latest = dates.max() if not dates.empty else None

            effective_first = first_detected or (str(meta_earliest.date()) if pd.notna(meta_earliest) else "")
            effective_last = last_detected or (str(meta_latest.date()) if pd.notna(meta_latest) else "")

            if effective_first and effective_last:
                summary = (
                    f"Lineage {lineage} ({lineage_id}) is documented from {effective_first} to {effective_last} "
                    f"across {int(genome_count)} curated genomes."
                )
            else:
                summary = (
                    f"Lineage {lineage} ({lineage_id}) is present in the reference database with "
                    f"{int(genome_count)} curated genomes; first/last detection dates are not curated."
                )

            findings.append(
                EvidenceRecord(
                    assertion=summary,
                    source="PostgreSQL PGIRL database (lineages + genome_metadata)",
                    confidence="high" if genome_count > 0 else "medium",
                    biological_relevance="lineage assignment and commonness",
                    epidemiological_relevance="historical circulation",
                    finding_type="lineage_commonness",
                )
            )
        else:
            findings.append(
                EvidenceRecord(
                    assertion=f"Lineage {lineage} is not found in the curated lineage table; treat as unassigned or uncharacterized.",
                    source="PostgreSQL PGIRL database (lineages)",
                    confidence="medium",
                    biological_relevance="lineage assignment",
                    epidemiological_relevance="unknown historical circulation",
                    finding_type="lineage_unassigned",
                )
            )
            gaps.append(
                f"Curated lineage metadata for {lineage} is missing; lineage commonness and history cannot be assessed."
            )

        # --- Mutation profile typicality --------------------------------------
        typical_variants = []
        atypical_variants = []
        absent_variants = []
        synonymous_variants = []
        variant_contexts: dict[str, dict] = {}

        for v in variants:
            hgvs = v.get("hgvs_p", "")
            gene = str(v.get("gene", "")).strip()
            pos = pd.to_numeric(v.get("position"), errors="coerce")
            ref_aa = str(v.get("ref_aa", "")).strip().upper()
            alt_aa = str(v.get("alt_aa", "")).strip().upper()

            if not gene or pd.isna(pos):
                continue

            if ref_aa == alt_aa:
                synonymous_variants.append(hgvs)
                continue

            if self.protein_variants.empty or "gene" not in self.protein_variants.columns:
                absent_variants.append(hgvs)
                continue

            # Look up curated protein variant frequencies
            mask = (
                (self.protein_variants["gene"].str.upper() == gene.upper())
                & (self.protein_variants["position"] == pos)
                & (self.protein_variants["alt_aa"].str.upper() == alt_aa.upper())
            )
            rows = self.protein_variants[mask]
            if rows.empty:
                absent_variants.append(hgvs)
                continue

            top = rows.iloc[0]
            variant_contexts[hgvs] = {
                "genome_count": int(top.get("genome_count", 0) or 0),
                "total_genomes": int(top.get("species_total_genomes", 0) or 0),
                "first_seen": top.get("first_seen_date"),
                "last_seen": top.get("last_seen_date"),
                "countries_seen": top.get("countries_seen"),
                "lineage_ids": top.get("lineage_ids"),
            }

            # Determine whether the variant is typical for the assigned lineage
            variant_lineages = top.get("lineage_ids")
            if isinstance(variant_lineages, str):
                variant_lineages = _parse_pg_array(variant_lineages)
            if not isinstance(variant_lineages, list):
                variant_lineages = []
            mapped_lineage_id = lineage_row.get("lineage_id") if lineage_row is not None else None
            if mapped_lineage_id and mapped_lineage_id in variant_lineages:
                typical_variants.append(hgvs)
            else:
                atypical_variants.append(hgvs)

        if typical_variants:
            findings.append(
                EvidenceRecord(
                    assertion=f"Variants typical for the assigned lineage are present: {', '.join(typical_variants)}.",
                    source="PostgreSQL PGIRL database (v_variant_summary)",
                    confidence="high",
                    biological_relevance="expected lineage mutation profile",
                    epidemiological_relevance="consistent with known lineage circulation",
                    finding_type="mutation_profile_typical",
                )
            )

        if atypical_variants:
            findings.append(
                EvidenceRecord(
                    assertion=f"Curated variants not typical for the assigned lineage: {', '.join(atypical_variants)}.",
                    source="PostgreSQL PGIRL database (v_variant_summary)",
                    confidence="medium",
                    biological_relevance="atypical mutation profile",
                    epidemiological_relevance="may represent local evolution, laboratory variation, or a new introduction",
                    finding_type="mutation_profile_atypical",
                )
            )

        if absent_variants:
            findings.append(
                EvidenceRecord(
                    assertion=f"Variants absent from curated frequency records: {', '.join(absent_variants)}.",
                    source="PostgreSQL PGIRL database (v_variant_summary)",
                    confidence="medium",
                    biological_relevance="uncharacterized mutation profile",
                    epidemiological_relevance="frequency and lineage distribution are unknown",
                    finding_type="mutation_profile_uncharacterized",
                )
            )

        if synonymous_variants:
            findings.append(
                EvidenceRecord(
                    assertion=f"Synonymous/nonsynonymous-silent amino-acid variants: {', '.join(synonymous_variants)}.",
                    source="bioinformatics output",
                    confidence="high",
                    biological_relevance="likely neutral at protein level; may still affect RNA structure or splicing",
                    epidemiological_relevance="limited unless linked to transmission lineage markers",
                    finding_type="mutation_profile_synonymous",
                )
            )

        # --- Functional domains / hotspots ------------------------------------
        hotspot_variants = [v for v in variants if v.get("known_hotspot")]
        if hotspot_variants:
            hgvs_list = [v.get("hgvs_p", "") for v in hotspot_variants]
            domains = [v.get("domain", "") for v in hotspot_variants]
            findings.append(
                EvidenceRecord(
                    assertion=f"Variant(s) in known functional domains/hotspots: {', '.join(hgvs_list)} (domains: {', '.join(d for d in domains if d)}).",
                    source="bioinformatics output + domain annotations",
                    confidence="high",
                    biological_relevance="functional domain mutation",
                    epidemiological_relevance="potential impact on phenotype",
                    finding_type="functional_hotspot",
                )
            )

        # --- Known phenotype associations -------------------------------------
        if matched_phenotypes:
            for p in matched_phenotypes:
                findings.append(
                    EvidenceRecord(
                        assertion=f"{p.get('genotype_description', p.get('hgvs_p'))} is associated with {p.get('phenotype_category')}: {p.get('phenotype_specific')}.",
                        source="PostgreSQL PGIRL database (genotype_phenotype)",
                        confidence=(p.get("evidence_strength") or "medium").lower(),
                        biological_relevance=p.get("phenotype_category", ""),
                        epidemiological_relevance="known phenotypic effect",
                        supporting_refs=_parse_pg_array(str(p.get("literature_refs", ""))),
                        record_flagged=p.get("record_flagged") == "t" or p.get("record_flagged") is True,
                        finding_type="known_phenotype",
                    )
                )
        else:
            gaps.append(
                "No curated genotype-phenotype associations were found for the detected variants."
            )

        # --- Negative findings: important absences ----------------------------
        negative_findings = self._negative_findings(variants)
        for neg in negative_findings:
            findings.append(neg)

        summary = self._summarize_genomic_significance(
            lineage, lineage_row is not None, typical_variants, atypical_variants, absent_variants, hotspot_variants, matched_phenotypes
        )

        return AnalysisResult(
            title="genomic_significance",
            findings=findings,
            gaps=gaps,
            summary=summary,
        )

    def _negative_findings(self, variants: list[dict]) -> list[EvidenceRecord]:
        """Report absence of known variants of public-health concern."""
        if self.associations.empty:
            return []

        observed = set()
        for v in variants:
            gene = str(v.get("gene", "")).strip().upper()
            pos = v.get("position")
            alt = str(v.get("alt_aa", "")).strip().upper()
            if gene and not pd.isna(pos) and alt:
                observed.add((gene, int(pos), alt))

        findings = []
        categories_of_interest = {
            "vaccine_escape",
            "vaccine_effectiveness",
            "drug_resistance",
            "drug_susceptibility",
            "diagnostic_sensitivity",
            "increased_transmission",
        }

        for cat in categories_of_interest:
            cat_rows = self.associations[
                self.associations["phenotype_category"].str.lower() == cat
            ]
            if cat_rows.empty:
                continue

            present = []
            absent = []
            for _, row in cat_rows.iterrows():
                gene = str(row.get("protein", "")).strip().upper()
                pos = row.get("position")
                alt = str(row.get("alt_aa", "")).strip().upper()
                if not gene or pd.isna(pos) or not alt:
                    continue
                desc = str(row.get("genotype_description", "")).strip()
                if (gene, int(pos), alt) in observed:
                    present.append(desc)
                else:
                    absent.append(desc)

            if absent:
                findings.append(
                    EvidenceRecord(
                        assertion=f"No detected variants with known {cat} association: {', '.join(absent[:5])}{' ...' if len(absent) > 5 else ''}.",
                        source="PostgreSQL PGIRL database (genotype_phenotype)",
                        confidence="high",
                        biological_relevance=f"absence of {cat}",
                        epidemiological_relevance="supports current diagnostic/vaccine/therapeutic effectiveness",
                        finding_type="negative_finding",
                    )
                )

            if present:
                findings.append(
                    EvidenceRecord(
                        assertion=f"Detected variant(s) with known {cat} association: {', '.join(present)}.",
                        source="PostgreSQL PGIRL database (genotype_phenotype)",
                        confidence="high",
                        biological_relevance=cat,
                        epidemiological_relevance="potential intervention impact",
                        finding_type="positive_intervention_concern",
                    )
                )

        return findings

    @staticmethod
    def _summarize_genomic_significance(
        lineage: str,
        lineage_known: bool,
        typical: list[str],
        atypical: list[str],
        absent: list[str],
        hotspots: list[dict],
        matched_phenotypes: list[dict],
    ) -> str:
        parts = []
        if lineage_known:
            parts.append(f"Lineage {lineage} is represented in the reference database.")
        else:
            parts.append(f"Lineage {lineage} is not represented in the reference database.")

        if typical:
            parts.append(f"Typical variants observed: {', '.join(typical)}.")
        if atypical:
            parts.append(f"Curated variants not typical for the lineage: {', '.join(atypical)}.")
        if absent:
            parts.append(f"Variants absent from curated frequency records: {', '.join(absent)}.")
        if hotspots:
            parts.append(
                f"{len(hotspots)} variant(s) fall in known functional domains/hotspots."
            )
        if matched_phenotypes:
            cats = {p.get("phenotype_category") for p in matched_phenotypes if p.get("phenotype_category")}
            if cats:
                parts.append(f"Curated phenotype associations: {', '.join(sorted(cats))}.")
            else:
                parts.append("Curated phenotype associations exist but lack category labels.")
        else:
            parts.append("No curated phenotype associations detected for the sample variants.")

        return " ".join(parts)


# ---------------------------------------------------------------------------
# Molecular epidemiology analyzer
# ---------------------------------------------------------------------------

class MolecularEpidemiologyAnalyzer:
    """Place the detected genome in historical, geographic and host context."""

    def __init__(self, lineages_df: pd.DataFrame, genome_metadata_df: pd.DataFrame):
        self.lineages = lineages_df
        self.genome_metadata = genome_metadata_df

    def analyze(self, stage9: dict) -> AnalysisResult:
        """Run molecular-epidemiology contextualization."""
        findings: list[EvidenceRecord] = []
        gaps: list[str] = []

        lineage = stage9.get("lineage")
        sample_meta = stage9.get("metadata", {})
        sample_country = sample_meta.get("country")
        sample_date = sample_meta.get("collection_date")
        sample_host = sample_meta.get("host") or sample_meta.get("host_species")

        if not lineage:
            gaps.append("No lineage assignment is available; molecular epidemiology cannot be contextualized.")
            return AnalysisResult(
                title="molecular_epidemiology",
                findings=findings,
                gaps=gaps,
                summary="Lineage information is missing; molecular epidemiology analysis cannot be performed.",
            )

        lineage_row = _match_lineage(lineage, self.lineages)
        if lineage_row is None:
            findings.append(
                EvidenceRecord(
                    assertion=f"Lineage {lineage} is not present in the curated lineage table.",
                    source="PostgreSQL PGIRL database (lineages)",
                    confidence="medium",
                    biological_relevance="lineage assignment",
                    epidemiological_relevance="unknown historical circulation",
                    finding_type="lineage_unassigned",
                )
            )
            gaps.append(f"Curated lineage metadata for {lineage} is missing; geographic and temporal context cannot be assessed.")
            return AnalysisResult(
                title="molecular_epidemiology",
                findings=findings,
                gaps=gaps,
                summary=f"Lineage {lineage} is not present in the curated lineage table; molecular epidemiology context is unavailable.",
            )

        lineage_id = lineage_row.get("lineage_id")
        first_detected = _safe_str(lineage_row.get("first_detected"))
        last_detected = _safe_str(lineage_row.get("last_detected"))
        first_country = _safe_str(lineage_row.get("first_country_detected"))
        primary_host = _safe_str(lineage_row.get("primary_host"))
        reservoir = _safe_str(lineage_row.get("reservoir"))
        countries_reported = _parse_pg_array(_safe_str(lineage_row.get("countries_reported")))
        regions_reported = _parse_pg_array(_safe_str(lineage_row.get("regions_reported")))
        endemic_regions = _parse_pg_array(_safe_str(lineage_row.get("endemic_regions")))

        # Derive context from curated genome metadata as a fallback/augmentation
        meta_rows = self._get_metadata_for_lineage(lineage_id)
        if not meta_rows.empty:
            earliest_meta = meta_rows["collection_date"].min()
            latest_meta = meta_rows["collection_date"].max()
            countries_meta = sorted({
                _safe_str(x) for x in meta_rows["collection_country"].dropna().unique()
                if _safe_str(x)
            })
            hosts_meta = sorted({
                _safe_str(x) for x in meta_rows["host"].dropna().unique()
                if _safe_str(x)
            })
            n_genomes = len(meta_rows)
        else:
            earliest_meta = latest_meta = None
            countries_meta = []
            hosts_meta = []
            n_genomes = 0

        effective_countries = list(set(countries_reported) | set(countries_meta)) if (countries_reported or countries_meta) else []
        effective_first = first_detected or (str(earliest_meta.date()) if pd.notna(earliest_meta) else "")
        effective_last = last_detected or (str(latest_meta.date()) if pd.notna(latest_meta) else "")
        effective_hosts = list({h for h in ([primary_host] if primary_host else []) + hosts_meta if h})

        # Historical circulation
        if effective_first and effective_last:
            findings.append(
                EvidenceRecord(
                    assertion=f"Lineage {lineage} ({lineage_id}) has curated genomes spanning {effective_first} to {effective_last} ({n_genomes} genomes).",
                    source="PostgreSQL PGIRL database (lineages + genome_metadata)",
                    confidence="high" if n_genomes > 0 else "medium",
                    biological_relevance="lineage history",
                    epidemiological_relevance="temporal circulation pattern",
                    finding_type="lineage_history",
                )
            )
        else:
            gaps.append(f"First/last detection dates for lineage {lineage} are not curated.")

        # Geographic context
        if effective_countries:
            if sample_country:
                if sample_country in effective_countries:
                    findings.append(
                        EvidenceRecord(
                            assertion=f"{sample_country} has previously reported lineage {lineage}.",
                            source="PostgreSQL PGIRL database (lineages + genome_metadata)",
                            confidence="high",
                            biological_relevance="geographic distribution",
                            epidemiological_relevance="known circulation in country",
                            finding_type="country_previously_reported",
                        )
                    )
                else:
                    findings.append(
                        EvidenceRecord(
                            assertion=f"{sample_country} is a novel geographic detection for lineage {lineage}; previously reported countries include {', '.join(effective_countries[:8])}{' ...' if len(effective_countries) > 8 else ''}.",
                            source="PostgreSQL PGIRL database (lineages + genome_metadata)",
                            confidence="high",
                            biological_relevance="geographic distribution",
                            epidemiological_relevance="potential geographic expansion or new introduction",
                            finding_type="novel_geographic_detection",
                        )
                    )
        else:
            gaps.append(f"No curated country list is available for lineage {lineage}.")

        # Host context
        if effective_hosts:
            findings.append(
                EvidenceRecord(
                    assertion=f"Lineage {lineage} has been associated with hosts: {', '.join(effective_hosts)}.",
                    source="PostgreSQL PGIRL database (lineages + genome_metadata)",
                    confidence="high" if primary_host else "medium",
                    biological_relevance="host range",
                    epidemiological_relevance="transmission ecology",
                    finding_type="host_context",
                )
            )
        else:
            gaps.append(f"No curated host information is available for lineage {lineage}.")

        # Temporal scenario
        scenario, scenario_conf = self._classify_temporal_scenario(
            sample_country, sample_date, effective_countries, effective_first, effective_last
        )
        if scenario:
            findings.append(
                EvidenceRecord(
                    assertion=f"This detection is interpreted as: {scenario}.",
                    source="PostgreSQL PGIRL database (lineages + genome_metadata) + sample metadata",
                    confidence=scenario_conf,
                    biological_relevance="evolutionary and epidemiological interpretation",
                    epidemiological_relevance="surveillance interpretation",
                    finding_type="temporal_scenario",
                )
            )

        summary = self._summarize(
            lineage, lineage_id, sample_country, sample_date, effective_countries, effective_first, effective_last, scenario
        )

        return AnalysisResult(
            title="molecular_epidemiology",
            findings=findings,
            gaps=gaps,
            summary=summary,
        )

    def _get_metadata_for_lineage(self, lineage_id: Any) -> pd.DataFrame:
        """Return genome metadata rows matching the lineage."""
        if self.genome_metadata.empty or "lineage_id" not in self.genome_metadata.columns or not lineage_id:
            return pd.DataFrame()
        mask = self.genome_metadata["lineage_id"].astype(str).str.strip().str.lower() == str(lineage_id).strip().lower()
        df = self.genome_metadata[mask].copy()
        if "collection_date" in df.columns:
            df["collection_date"] = pd.to_datetime(df["collection_date"], errors="coerce")
        return df

    @staticmethod
    def _classify_temporal_scenario(
        sample_country: Optional[str],
        sample_date_str: Optional[str],
        countries_reported: list[str],
        first_detected: str,
        last_detected: str,
    ) -> tuple[str, str]:
        """Classify the detection scenario based on geography and timing."""
        if not sample_date_str or not first_detected or not last_detected:
            return ("", "")

        sample_date = pd.to_datetime(sample_date_str, errors="coerce")
        first = pd.to_datetime(first_detected, errors="coerce")
        last = pd.to_datetime(last_detected, errors="coerce")
        if pd.isna(sample_date) or pd.isna(first) or pd.isna(last):
            return ("", "")

        in_country = bool(sample_country and sample_country in countries_reported)
        after_last = sample_date > last
        year_gap = (sample_date - last).days / 365.25 if after_last else 0

        if in_country:
            if sample_date >= first and sample_date <= last:
                return ("continued circulation", "high")
            if after_last and year_gap <= 2:
                return ("continued circulation", "medium")
            if after_last:
                return ("re-emergence after a reporting gap", "medium")
            return ("continued circulation", "medium")

        if sample_date >= first and sample_date <= last:
            return ("geographic expansion within known circulation period", "medium")
        if after_last:
            return ("new introduction or geographic expansion", "medium")
        return ("uncertain temporal/geographic relationship", "low")

    @staticmethod
    def _summarize(
        lineage: str,
        lineage_id: Any,
        sample_country: Optional[str],
        sample_date: Optional[str],
        countries: list[str],
        first_detected: str,
        last_detected: str,
        scenario: str,
    ) -> str:
        parts = [f"Lineage {lineage} ({lineage_id})"]
        if first_detected and last_detected:
            parts.append(f"has been detected from {first_detected} to {last_detected}.")
        else:
            parts.append("has limited temporal metadata.")

        if countries:
            parts.append(f"Previously reported countries include {', '.join(countries[:5])}{' ...' if len(countries) > 5 else ''}.")
        if sample_country and sample_date:
            parts.append(f"The current sample from {sample_country} ({sample_date}) represents {scenario or 'an undetermined epidemiological scenario'}.")
        return " ".join(parts)


# ---------------------------------------------------------------------------
# Phase 2 decision-oriented analyzers
# ---------------------------------------------------------------------------

class PhylogeographicAnalyzer:
    """Estimate likely origin and dissemination from tree + metadata."""

    def __init__(
        self,
        tree_input: TreeInput,
        lineages_df: pd.DataFrame,
        genome_metadata_df: pd.DataFrame,
    ):
        self.tree = tree_input
        self.lineages = lineages_df
        self.genome_metadata = genome_metadata_df

    def analyze(self, stage9: dict) -> AnalysisResult:
        lineage = _safe_str(stage9.get("lineage"))
        sample_country = _safe_str(stage9.get("metadata", {}).get("country"))
        sample_date = _parse_date(stage9.get("metadata", {}).get("collection_date"))
        lineage_row = _match_lineage(lineage, self.lineages)

        origin, confidence, routes, origin_support, introduction_scenarios = self._infer_origin_and_routes(
            lineage_row, sample_country, sample_date
        )

        findings: list[EvidenceRecord] = []
        gaps: list[str] = []
        metrics: dict[str, Any] = {}

        if origin:
            findings.append(
                EvidenceRecord(
                    assertion=f"Most likely geographic origin of the detected lineage: {origin}.",
                    source="phylogenetic tree + genome_metadata.csv + lineages.csv",
                    confidence=confidence,
                    biological_relevance="phylogeographic origin",
                    epidemiological_relevance="identifies probable source population or location",
                    finding_type="phylogeographic_origin",
                )
            )
        else:
            gaps.append("Insufficient data to infer phylogeographic origin.")

        if routes:
            findings.append(
                EvidenceRecord(
                    assertion=f"Plausible dissemination pathways include: {', '.join(routes)}.",
                    source="phylogenetic tree + genome_metadata.csv",
                    confidence="medium",
                    biological_relevance="dissemination history",
                    epidemiological_relevance="identifies linked regions and potential transmission routes",
                    finding_type="dissemination_pathways",
                )
            )

        if sample_country and origin and sample_country.lower() != origin.lower():
            findings.append(
                EvidenceRecord(
                    assertion=f"Current detection in {sample_country} differs from the inferred origin ({origin}), suggesting importation or expansion.",
                    source="phylogenetic tree + sample metadata",
                    confidence="medium",
                    biological_relevance="geographic mismatch",
                    epidemiological_relevance="possible introduction from origin region",
                    finding_type="geographic_mismatch",
                )
            )

        summary = f"Phylogeographic assessment: most likely origin is {origin or 'unknown'}"
        if routes:
            summary += f"; dissemination routes: {', '.join(routes[:3])}"
        if self.tree.has_tree:
            summary += " (tree-informed)."
        else:
            summary += " (metadata fallback, lower confidence)."

        metrics = {
            "inferred_origin": origin,
            "origin_confidence": confidence,
            "origin_support": origin_support,
            "introduction_scenarios": introduction_scenarios,
            "dissemination_routes": routes,
        }

        return AnalysisResult(
            title="phylogeographic_analysis",
            findings=findings,
            gaps=gaps,
            summary=summary,
            metrics=metrics,
        )

    def _infer_origin_and_routes(
        self,
        lineage_row: Optional[pd.Series],
        sample_country: str,
        sample_date: str,
    ) -> tuple[str, str, list[str], dict[str, Any], dict[str, Any]]:
        """Return (origin_country, confidence, route_list, origin_support, introduction_scenarios)."""
        origin_support: dict[str, Any] = {}
        introduction_scenarios: dict[str, Any] = {}

        # 1. Tree-informed nearest neighbors
        if self.tree.has_tree and self.tree.sample_tip:
            nearest = self.tree.get_nearest_tips(n=10)
            dated = [t for t in nearest if t.country]
            if dated:
                # Count support for each candidate origin country among nearest tips
                country_counts: Counter[str] = Counter(t.country for t in dated)
                total = sum(country_counts.values())
                origin = dated[0].country
                confidence = "high" if len(dated) >= 3 else "medium"
                origin_support = {
                    "method": "nearest_tree_tips",
                    "nearest_tip_countries": dict(country_counts.most_common()),
                    "origin_probability": {c: round(n / total, 3) for c, n in country_counts.most_common()},
                }

                # Build dissemination routes from earliest to most recent nearest tip country
                with_dates = [
                    (t.country, pd.to_datetime(t.date, errors="coerce"))
                    for t in dated
                    if t.country and t.date
                ]
                with_dates = [(c, d) for c, d in with_dates if pd.notna(d)]
                with_dates.sort(key=lambda x: x[1])
                routes = []
                for i in range(1, len(with_dates)):
                    prev, curr = with_dates[i - 1][0], with_dates[i][0]
                    if prev != curr:
                        route = f"{prev} -> {curr}"
                        if route not in routes:
                            routes.append(route)

                # Introduction scenarios
                introduction_scenarios = self._introduction_scenarios(
                    sample_country, sample_date, origin, dated
                )
                return origin, confidence, routes[:5], origin_support, introduction_scenarios

        # 2. Metadata fallback
        lineage_id = _safe_str(lineage_row.get("lineage_id")) if lineage_row is not None else ""
        first_country = _safe_str(lineage_row.get("first_country_detected")) if lineage_row is not None else ""
        countries_reported = _parse_pg_array(
            _safe_str(lineage_row.get("countries_reported"))
        ) if lineage_row is not None else []

        if first_country:
            origin_support = {"method": "first_country_detected_curated", "evidence": first_country}
            introduction_scenarios = self._introduction_scenarios(
                sample_country, sample_date, first_country, []
            )
            return first_country, "medium", [], origin_support, introduction_scenarios

        if countries_reported:
            origin_support = {"method": "countries_reported_list", "evidence": countries_reported}
            introduction_scenarios = self._introduction_scenarios(
                sample_country, sample_date, countries_reported[0], []
            )
            return countries_reported[0], "low", [], origin_support, introduction_scenarios

        if not self.genome_metadata.empty and lineage_id and "lineage_id" in self.genome_metadata.columns:
            mask = self.genome_metadata["lineage_id"].astype(str).str.strip().str.lower() == lineage_id.lower()
            rows = self.genome_metadata[mask]
            if not rows.empty and "collection_country" in rows.columns:
                country_counts = rows["collection_country"].dropna().astype(str).str.strip().value_counts()
                if not country_counts.empty:
                    top_country = country_counts.index[0]
                    total = int(country_counts.sum())
                    origin_support = {
                        "method": "most_common_collection_country",
                        "country_counts": {c: int(n) for c, n in country_counts.items()},
                        "origin_probability": {c: round(n / total, 3) for c, n in country_counts.items()},
                    }
                    introduction_scenarios = self._introduction_scenarios(
                        sample_country, sample_date, top_country, []
                    )
                    return top_country, "low", [], origin_support, introduction_scenarios

        return "", "low", [], origin_support, introduction_scenarios

    def _introduction_scenarios(
        self,
        sample_country: str,
        sample_date: str,
        inferred_origin: str,
        nearest_tips: list[Any],
    ) -> dict[str, Any]:
        """Score plausible introduction scenarios relative to the inferred origin."""
        scenarios: dict[str, Any] = {
            "new_introduction": {"score": 0, "evidence": []},
            "re_emergence": {"score": 0, "evidence": []},
            "local_persistence": {"score": 0, "evidence": []},
        }
        if not sample_country or not inferred_origin:
            return scenarios

        same_country = sample_country.lower() == inferred_origin.lower()

        # New introduction: sample country differs from origin, and nearest tips are not local
        if not same_country:
            scenarios["new_introduction"]["score"] += 3
            scenarios["new_introduction"]["evidence"].append(
                f"Sample country {sample_country} differs from inferred origin {inferred_origin}."
            )

        # Re-emergence: sample date is substantially later than the most recent nearest tip
        if nearest_tips:
            dates = [
                pd.to_datetime(t.date, errors="coerce")
                for t in nearest_tips
                if t.date
            ]
            dates = [d for d in dates if pd.notna(d)]
            if dates and sample_date:
                sample_dt = _to_date(sample_date)
                if sample_dt is not None:
                    most_recent = max(dates)
                    gap_days = (sample_dt - most_recent).days
                    if gap_days > 365:
                        scenarios["re_emergence"]["score"] += 3
                        scenarios["re_emergence"]["evidence"].append(
                            f"Sample date is {gap_days} days after the most recent contextual genome."
                        )
                    else:
                        scenarios["local_persistence"]["score"] += 1
                        scenarios["local_persistence"]["evidence"].append(
                            "Contextual genomes detected within one year of sample date."
                        )

        # Local persistence: same country as origin
        if same_country:
            scenarios["local_persistence"]["score"] += 2
            scenarios["local_persistence"]["evidence"].append(
                "Sample country matches inferred origin."
            )

        return scenarios


class GeneticRelatednessAnalyzer:
    """Report pairwise genetic distances and their uncertainty.

    This analyzer deliberately does **not** infer transmission networks: those
    require contact-tracing, exposure timing and sampling-intensity data that
    genomic sequences alone cannot provide.
    """

    CLOSE_SNP_THRESHOLD = 5  # within ~5 SNPs is consistent with close relatedness

    def __init__(
        self,
        tree_input: TreeInput,
        genome_metadata_df: pd.DataFrame,
        lineages_df: pd.DataFrame,
    ):
        self.tree = tree_input
        self.genome_metadata = genome_metadata_df
        self.lineages = lineages_df

    def analyze(self, stage9: dict) -> AnalysisResult:
        sample_country = _safe_str(stage9.get("metadata", {}).get("country"))
        sample_date = _parse_date(stage9.get("metadata", {}).get("collection_date"))
        lineage = _safe_str(stage9.get("lineage"))

        placement = stage9.get("phylogenetic_placement", {})
        closest_accession = _safe_str(placement.get("closest_genome_accession"))
        closest_distance = placement.get("snps_from_closest")

        findings: list[EvidenceRecord] = []
        gaps: list[str] = []
        metrics: dict[str, Any] = {}

        # Collect nearest contextual genomes from tree or placement metadata
        nearest: list[dict] = []
        if self.tree.has_tree:
            for tip in self.tree.get_nearest_tips(n=5):
                nearest.append({
                    "name": tip.name,
                    "country": tip.country,
                    "date": tip.date,
                    "distance_from_root": tip.distance_from_root,
                })

        if closest_distance is not None:
            try:
                metrics["snps_to_closest_genome"] = int(closest_distance)
                metrics["closest_genome_accession"] = closest_accession
            except (ValueError, TypeError):
                pass

        # Temporal gaps to nearest dated genomes
        temporal_gaps: list[int] = []
        if sample_date and self.tree.has_tree:
            sample_dt = _to_date(sample_date)
            for tip in self.tree.get_nearest_tips(n=10):
                tip_dt = _to_date(tip.date)
                if sample_dt is not None and tip_dt is not None:
                    try:
                        gap = abs((sample_dt - tip_dt).days)
                        if not pd.isna(gap):
                            temporal_gaps.append(int(gap))
                    except (ValueError, TypeError):
                        pass
        if temporal_gaps:
            metrics["min_temporal_gap_days"] = int(min(temporal_gaps))
            median_gap = pd.Series(temporal_gaps).median()
            if not pd.isna(median_gap):
                metrics["median_temporal_gap_days"] = int(median_gap)

        # Report the closest genetic match as relatedness, not transmission
        if metrics.get("snps_to_closest_genome") is not None:
            snps = metrics["snps_to_closest_genome"]
            close = snps <= self.CLOSE_SNP_THRESHOLD
            findings.append(
                EvidenceRecord(
                    assertion=(
                        f"The sample differs by {snps} SNPs from its closest curated genome "
                        f"({closest_accession or 'unknown'}). "
                        f"This is {'consistent with close genetic relatedness' if close else 'not within the close-relatedness window'} "
                        f"(threshold {self.CLOSE_SNP_THRESHOLD} SNPs)."
                    ),
                    source="phylogenetic placement + genome metadata",
                    confidence="medium",
                    biological_relevance="genetic relatedness",
                    epidemiological_relevance="genomic distance to nearest known sequence",
                    finding_type="close_genetic_match" if close else "genetic_relatedness_summary",
                )
            )
        elif nearest:
            names = ", ".join(n.get("name", "") for n in nearest[:3] if n.get("name"))
            findings.append(
                EvidenceRecord(
                    assertion=f"Nearest contextual tree tips: {names}. SNP distances were not supplied.",
                    source="phylogenetic tree",
                    confidence="low",
                    biological_relevance="genetic relatedness",
                    epidemiological_relevance="nearest known sequences in tree topology",
                    finding_type="genetic_relatedness_summary",
                )
            )
        else:
            gaps.append("No phylogenetic placement or tree tips available; genetic relatedness cannot be quantified.")
            findings.append(
                EvidenceRecord(
                    assertion="Genetic relatedness to curated genomes cannot be assessed because neither SNP distances nor a phylogenetic tree were provided.",
                    source="genome metadata",
                    confidence="low",
                    biological_relevance="genetic relatedness",
                    epidemiological_relevance="unknown relationship to other cases",
                    finding_type="no_close_matches",
                )
            )

        if metrics.get("min_temporal_gap_days") is not None:
            gap_days = metrics["min_temporal_gap_days"]
            findings.append(
                EvidenceRecord(
                    assertion=f"The closest dated contextual genome is {gap_days} days from the sample collection date.",
                    source="phylogenetic tree + sample metadata",
                    confidence="medium",
                    biological_relevance="temporal context",
                    epidemiological_relevance="temporal proximity to nearest sampled case",
                    finding_type="temporal_gap",
                )
            )

        # Explicit limitation statement
        findings.append(
            EvidenceRecord(
                assertion="Transmission clusters and networks cannot be reliably inferred from genetic distance alone. Contact tracing, exposure history, and sampling intensity are required.",
                source="methodological limitation",
                confidence="high",
                biological_relevance="genetic relatedness",
                epidemiological_relevance="transmission inference",
                finding_type="transmission_inference_limitation",
            )
        )

        summary_parts = [
            f"Genetic relatedness to nearest curated genome: {metrics.get('snps_to_closest_genome', 'unknown')} SNPs",
        ]
        if metrics.get("min_temporal_gap_days") is not None:
            summary_parts.append(f"minimum temporal gap {metrics['min_temporal_gap_days']} days")
        summary_parts.append("transmission networks cannot be inferred from sequence data alone.")

        return AnalysisResult(
            title="genetic_relatedness_analysis",
            findings=findings,
            gaps=gaps,
            summary=" ".join(summary_parts),
            metrics=metrics,
        )


class MolecularClockAnalyzer:
    """Estimate timing of emergence and introduction from time-scaled tree metadata."""

    def __init__(self, tree_input: TreeInput):
        self.tree = tree_input

    def analyze(self, stage9: dict) -> AnalysisResult:
        time_scaled = (stage9.get("stage6_phylogenetic_tree") or {}).get("time_scaled_tree") or {}
        sample_date = _parse_date(stage9.get("metadata", {}).get("collection_date"))
        lineage = _safe_str(stage9.get("lineage"))

        clock_rate = time_scaled.get("molecular_clock_rate")
        root_age = _safe_str(time_scaled.get("root_age"))
        placement_date = _safe_str(time_scaled.get("sample_placement_date"))

        findings: list[EvidenceRecord] = []
        gaps: list[str] = []

        if clock_rate and (root_age or placement_date):
            try:
                rate = float(clock_rate)
            except (ValueError, TypeError):
                rate = None

            if root_age:
                findings.append(
                    EvidenceRecord(
                        assertion=f"Time-scaled phylogeny estimates the most recent common ancestor (MRCA) of the sampled clade at {root_age} (clock rate {rate} subs/site/year).",
                        source="stage6_phylogenetic_tree.time_scaled_tree",
                        confidence="medium",
                        biological_relevance="evolutionary timing",
                        epidemiological_relevance="estimates when the lineage began circulating",
                        finding_type="molecular_clock_timing",
                    )
                )
            if placement_date:
                findings.append(
                    EvidenceRecord(
                        assertion=f"The sample is placed around {placement_date} on the time-scaled tree.",
                        source="stage6_phylogenetic_tree.time_scaled_tree",
                        confidence="medium",
                        biological_relevance="sample dating",
                        epidemiological_relevance="supports emergence or introduction timing",
                        finding_type="sample_placement_timing",
                    )
                )
            if sample_date and placement_date:
                try:
                    sample_dt = pd.to_datetime(sample_date)
                    place_dt = pd.to_datetime(placement_date)
                    delta = (sample_dt - place_dt).days
                    if delta > 30:
                        findings.append(
                            EvidenceRecord(
                                assertion=f"Sample collection ({sample_date}) is {delta} days after the estimated placement date, consistent with a reporting or transmission chain delay.",
                                source="stage6_phylogenetic_tree.time_scaled_tree + sample metadata",
                                confidence="low",
                                biological_relevance="evolutionary/epidemiological lag",
                                epidemiological_relevance="may indicate undetected intermediate cases",
                                finding_type="reporting_delay",
                            )
                        )
                except (ValueError, TypeError):
                    pass
        else:
            gaps.append(
                "No time-scaled phylogeny is available; molecular-clock estimates of emergence timing cannot be computed."
            )

        summary = "Molecular-clock assessment: "
        if root_age and placement_date:
            summary += f"MRCA estimated at {root_age}; sample placed around {placement_date}."
        elif root_age:
            summary += f"MRCA estimated at {root_age}."
        elif placement_date:
            summary += f"Sample placed around {placement_date}."
        else:
            summary += "insufficient time-scaled tree metadata."

        return AnalysisResult(
            title="molecular_clock_analysis",
            findings=findings,
            gaps=gaps,
            summary=summary,
        )


class EvidenceWeightedThreatAnalyzer:
    """LEGACY / TODO(migrate): combines genomic, epidemiological and
    historical evidence into a quantitative threat assessment.

    This is public-health-conclusion territory and conceptually belongs to a
    downstream Genomic Intelligence Engine stage, not this evidence
    harmonization/cross-evidence-analysis layer. It is kept in place for now
    because intelligence_pipeline.py's risk-tiering and report generation
    still depend on it; extracting it cleanly requires carving `_assess_risk`
    and this analyzer out into that separate downstream stage, which is a
    larger, distinct piece of work.
    """

    def __init__(self, engine: "GenomicIntelligenceEngine"):
        self.engine = engine

    def analyze(
        self,
        stage9: dict,
        all_results: dict[str, AnalysisResult],
        risk: dict,
    ) -> AnalysisResult:
        sample_country = _safe_str(stage9.get("metadata", {}).get("country"))

        # Build an index of evidence so every score point can be traced to a finding.
        finding_types: set[str] = set()
        assertion_text: list[str] = []
        for result in all_results.values():
            for f in result.findings:
                finding_types.add(f.finding_type)
                assertion_text.append(f.assertion.lower())
        all_text = " ".join(assertion_text)

        # Scoring framework: each criterion starts at 0 and accumulates points based on
        # explicit finding types. score_basis records the exact evidence for transparency.
        criteria: dict[str, dict[str, Any]] = {
            "local_transmission": {"rules": [
                ("close_genetic_match", 3, "genetically close genome detected"),
                ("country_previously_reported", 1, "lineage previously reported in sample country"),
                ("continued_circulation", 1, "lineage evidence for continued circulation"),
                ("re-emergence_after_a_reporting_gap", 1, "re-emergence or reporting gap detected"),
                ("r0_estimate", 1, "R0/Rt estimate available"),
            ], "max": 5},
            "geographic_expansion": {"rules": [
                ("novel_geographic_detection", 2, "novel geographic detection"),
                ("geographic_mismatch", 1, "sample country differs from inferred origin"),
                ("spillover_risk", 1, "animal spillover risk"),
            ], "max": 5},
            "international_spread": {"rules": [
                ("dissemination_pathways", 2, "inferred dissemination pathways"),
                ("genetic_relatedness_summary", 1, "distant genetic relatedness supports introduction"),
                ("international_spread_risk", 1, "prior international spread signal"),
                ("geographic_coverage_gap", 1, "geographic coverage gap"),
            ], "max": 5},
            "diagnostic_escape": {"rules": [
                ("known_phenotype", 1, "known phenotype association present"),
            ], "max": 5},
            "vaccine_impact": {"rules": [
                ("known_phenotype", 1, "known phenotype association present"),
                ("functional_hotspot", 1, "mutation in functional hotspot"),
                ("vaccine_escape", 1, "vaccine escape phenotype"),
            ], "max": 5},
            "therapeutic_relevance": {"rules": [
                ("known_phenotype", 1, "known phenotype association present"),
                ("functional_hotspot", 1, "mutation in functional hotspot"),
                ("drug_resistance", 1, "drug resistance phenotype"),
            ], "max": 5},
            "surveillance_priority": {"rules": [
                ("lineage_unassigned", 2, "lineage not in curated database"),
                ("mutation_profile_atypical", 1, "atypical mutation profile"),
                ("mutation_profile_uncharacterized", 1, "uncharacterized mutation profile"),
                ("re-emergence_after_a_reporting_gap", 2, "re-emergence or reporting gap"),
                ("geographic_coverage_gap", 1, "geographic coverage gap"),
                ("undetected_transmission_chain", 1, "possible undetected transmission chain"),
                ("reporting_delay", 1, "reporting delay or lag"),
            ], "max": 5},
        }

        scores: dict[str, int] = {}
        score_basis: dict[str, list[str]] = {}
        for criterion, cfg in criteria.items():
            basis: list[str] = []
            score = 0
            for finding_type, points, reason in cfg["rules"]:
                if finding_type == "continued_circulation":
                    if "continued circulation" in all_text or finding_type in finding_types:
                        score += points
                        basis.append(reason)
                elif finding_type == "re-emergence_after_a_reporting_gap":
                    if "re-emergence" in all_text or finding_type in finding_types:
                        score += points
                        basis.append(reason)
                elif finding_type == "genetic_relatedness_summary":
                    if finding_type in finding_types and "close_genetic_match" not in finding_types:
                        score += points
                        basis.append(reason)
                elif finding_type in finding_types:
                    score += points
                    basis.append(reason)
            scores[criterion] = min(score, cfg["max"])
            score_basis[criterion] = basis

        # Derive an overall label that requires multiple strong signals for the top tiers.
        risk_tier = risk.get("risk_tier", "routine")
        tier_rank = {"routine": 1, "monitor": 2, "investigate": 3, "high_priority": 4, "emergency": 5}
        max_score = max(scores.values())
        tier_value = tier_rank.get(risk_tier, 1)
        # Count how many criteria are elevated (>=3) to avoid a single noisy signal dominating.
        elevated_criteria = sum(1 for v in scores.values() if v >= 3)

        if max_score >= 4 and elevated_criteria >= 2 and tier_value >= 4:
            overall_label = "high_priority"
        elif max_score >= 4 or (elevated_criteria >= 2 and tier_value >= 3):
            overall_label = "investigate"
        elif max_score >= 3 or tier_value >= 2:
            overall_label = "monitor"
        else:
            overall_label = "routine"
        overall = tier_rank[overall_label]

        # Identify top concern criteria and include their evidence basis.
        top_criteria = sorted(scores.items(), key=lambda x: x[1], reverse=True)[:3]
        top_list = [f"{k.replace('_', ' ')} ({v}/5)" for k, v in top_criteria if v > 0]

        # Build transparent findings with traceable scoring basis.
        findings: list[EvidenceRecord] = [
            EvidenceRecord(
                assertion=(
                    f"Evidence-weighted assessment: overall rank {overall_label} "
                    f"(max criterion score {max_score}/5). Top concerns: {', '.join(top_list)}. "
                    "Each score is derived from explicit analyzer findings (see score_basis)."
                ),
                source="transparent scoring of integrated genomic, epidemiological and historical evidence",
                confidence="medium",
                biological_relevance="integrated threat profile",
                epidemiological_relevance="prioritises operational response",
                finding_type="evidence_weighted_assessment",
            )
        ]

        for criterion, cfg in criteria.items():
            if scores[criterion] > 0:
                basis = score_basis.get(criterion, [])
                basis_text = "; ".join(basis[:3]) + (" ..." if len(basis) > 3 else "")
                findings.append(
                    EvidenceRecord(
                        assertion=f"{criterion.replace('_', ' ').title()} score {scores[criterion]}/{cfg['max']}: {basis_text}" if basis_text else f"{criterion.replace('_', ' ').title()} score {scores[criterion]}/{cfg['max']}.",
                        source="transparent scoring framework",
                        confidence="medium",
                        biological_relevance=criterion,
                        epidemiological_relevance="score contribution",
                        finding_type=f"{criterion}_score",
                    )
                )

        gaps = []
        if "molecular_clock_timing" not in finding_types:
            gaps.append("Molecular-clock timing unavailable; temporal inference relies on metadata.")
        if "phylogeographic_origin" not in finding_types:
            gaps.append("Phylogeographic origin unavailable; source investigation is empirical.")

        summary = f"Quantitative assessment rank: {overall_label}. Highest-scoring criteria: {', '.join(top_list)}."

        metrics = {
            "scores": scores,
            "score_basis": {k: v for k, v in score_basis.items()},
            "overall_rank": overall,
            "overall_label": overall_label,
            "top_criteria": [{"criterion": k, "score": v, "basis": score_basis.get(k, [])} for k, v in top_criteria],
            "scoring_framework_version": "1.0",
        }

        return AnalysisResult(
            title="evidence_weighted_threat",
            findings=findings,
            gaps=gaps,
            summary=summary,
            metrics=metrics,
        )


# ---------------------------------------------------------------------------
# Engine orchestrator
# ---------------------------------------------------------------------------

class GenomicIntelligenceEngine:
    """Orchestrate the analyzers and produce the structured intelligence object."""

    def __init__(
        self,
        associations_csv_path: str = "database/exports/genotype_phenotype.csv",
        protein_variants_csv_path: str = "database/exports/protein_variants.csv",
        lineages_csv_path: str = "database/exports/lineages.csv",
        genome_metadata_csv_path: str = "database/exports/genome_metadata.csv",
        epi_output_path: Optional[str] = None,
        pathogen_id: Optional[str] = None,
        species_id: Optional[str] = None,
        db_url: Optional[str] = None,
    ):
        self.associations_csv_path = Path(associations_csv_path)
        self.protein_variants_csv_path = Path(protein_variants_csv_path)
        self.lineages_csv_path = Path(lineages_csv_path)
        self.genome_metadata_csv_path = Path(genome_metadata_csv_path)
        self.epi_output_path = Path(epi_output_path) if epi_output_path else None

        self.pathogen_id = pathogen_id
        self.species_id = species_id
        self.lineage_id: Optional[str] = None
        self.detected_variants: Optional[list[dict]] = None
        self.db_url = db_url

        self._associations: Optional[pd.DataFrame] = None
        self._protein_variants: Optional[pd.DataFrame] = None
        self._lineages: Optional[pd.DataFrame] = None
        self._genome_metadata: Optional[pd.DataFrame] = None
        self._epi_output: Optional[dict] = None

    def set_pathogen_context(
        self,
        pathogen_id: Optional[str] = None,
        species_id: Optional[str] = None,
        db_url: Optional[str] = None,
    ) -> None:
        """Set the pathogen/species context so reference data can be loaded from the DB."""
        self.pathogen_id = pathogen_id or self.pathogen_id
        self.species_id = species_id or self.species_id
        self.db_url = db_url or self.db_url

    def set_lineage_context(self, lineage_id: Optional[str] = None) -> None:
        """Set the sample lineage so contextual metadata queries can be narrowed."""
        self.lineage_id = lineage_id

    def set_detected_variants(self, variants: Optional[list[dict]] = None) -> None:
        """Set the detected variants so variant/phenotype lookups can be narrowed."""
        self.detected_variants = variants or []

    def get_reference_data(self) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """Force-load and return the four curated reference DataFrames."""
        self._load_data()
        return (
            self._associations if self._associations is not None else pd.DataFrame(),
            self._protein_variants if self._protein_variants is not None else pd.DataFrame(),
            self._lineages if self._lineages is not None else pd.DataFrame(),
            self._genome_metadata if self._genome_metadata is not None else pd.DataFrame(),
        )

    def _load_csv(self, path: Path, dtype: Optional[dict] = None) -> pd.DataFrame:
        if not path.exists():
            return pd.DataFrame()
        return pd.read_csv(path, dtype=dtype, low_memory=False)

    def _load_data(self) -> None:
        """Load only the reference rows actually needed for the current sample.

        When pathogen_id/species_id are set, the engine issues targeted DB
        queries: one lineage record, contextual genomes for that lineage,
        protein variants matching detected mutations, and phenotype
        associations for detected variants/motifs plus intervention categories.
        It only falls back to CSV exports when the DB is unreachable.
        """
        from intelligence_engine.evidence_integration.data_sources import db_loader

        if self._associations is None:
            self._associations = self._load_from_db_or_csv(
                lambda: db_loader.load_genotype_phenotype_for_variants(
                    self.pathogen_id, self.species_id, self.detected_variants or [], self.db_url
                )
                if self.detected_variants is not None
                else None,
                fallback=lambda: db_loader.load_genotype_phenotype(
                    self.pathogen_id, self.species_id, self.db_url
                )
                if self.pathogen_id and self.species_id
                else None,
                csv_path=self.associations_csv_path,
                label="genotype_phenotype",
            )
        if self._protein_variants is None:
            self._protein_variants = self._load_from_db_or_csv(
                lambda: db_loader.load_protein_variants_for_variants(
                    self.pathogen_id, self.species_id, self.detected_variants or [], self.db_url
                )
                if self.detected_variants is not None
                else None,
                fallback=lambda: db_loader.load_protein_variants(
                    self.pathogen_id, self.species_id, self.db_url
                )
                if self.pathogen_id and self.species_id
                else None,
                csv_path=self.protein_variants_csv_path,
                label="protein_variants",
            )
        if self._lineages is None:
            self._lineages = self._load_from_db_or_csv(
                lambda: db_loader.load_lineage(
                    self.pathogen_id, self.species_id, self.lineage_id, self.db_url
                )
                if self.pathogen_id and self.species_id
                else None,
                fallback=lambda: db_loader.load_lineages(
                    self.pathogen_id, self.species_id, self.db_url
                )
                if self.pathogen_id and self.species_id
                else None,
                csv_path=self.lineages_csv_path,
                label="lineages",
            )
        if self._genome_metadata is None:
            self._genome_metadata = self._load_from_db_or_csv(
                lambda: db_loader.load_genome_metadata_for_lineage(
                    self.pathogen_id, self.species_id, self.lineage_id, db_url=self.db_url
                )
                if self.pathogen_id and self.species_id
                else None,
                fallback=lambda: db_loader.load_genome_metadata(
                    self.pathogen_id, self.species_id, self.db_url
                )
                if self.pathogen_id and self.species_id
                else None,
                csv_path=self.genome_metadata_csv_path,
                label="genome_metadata",
            )
        if self._epi_output is None and self.epi_output_path and self.epi_output_path.exists():
            with open(self.epi_output_path) as f:
                self._epi_output = json.load(f)
        if self._epi_output is None:
            self._epi_output = {}

    def _load_from_db_or_csv(
        self,
        primary_query,
        fallback,
        csv_path: Path,
        label: str,
    ) -> pd.DataFrame:
        """Run the targeted DB query, then a broader fallback query, then CSV."""
        if self.pathogen_id and self.species_id:
            try:
                df = primary_query()
                if df is not None and not df.empty:
                    log.info(
                        "Loaded %s targeted rows of %s from DB for %s/%s",
                        len(df),
                        label,
                        self.pathogen_id,
                        self.species_id,
                    )
                    return df
            except Exception as exc:
                log.warning("Targeted DB load failed for %s: %s", label, exc)
            # If targeted query returned nothing, try a broader DB load before CSV.
            try:
                df = fallback()
                if df is not None and not df.empty:
                    log.info(
                        "Loaded %s rows of %s from DB (broad query) for %s/%s",
                        len(df),
                        label,
                        self.pathogen_id,
                        self.species_id,
                    )
                    return df
            except Exception as exc:
                log.warning("Broad DB load failed for %s: %s", label, exc)
        return self._load_csv(csv_path)

    def analyze_genomic_significance(
        self,
        stage9: dict,
        variants: list[dict],
        matched_phenotypes: list[dict],
    ) -> AnalysisResult:
        """Run the Phase 1 genomic-significance analyzer."""
        self._load_data()
        analyzer = GenomicSignificanceAnalyzer(
            associations_df=self._associations,
            protein_variants_df=self._protein_variants,
            lineages_df=self._lineages,
            genome_metadata_df=self._genome_metadata,
        )
        return analyzer.analyze(stage9, variants, matched_phenotypes)

    def analyze_molecular_epidemiology(self, stage9: dict) -> AnalysisResult:
        """Run the molecular-epidemiology contextualization analyzer."""
        self._load_data()
        analyzer = MolecularEpidemiologyAnalyzer(
            lineages_df=self._lineages,
            genome_metadata_df=self._genome_metadata,
        )
        return analyzer.analyze(stage9)

    def get_lineage_metadata(self, lineage: str) -> dict:
        """Return curated lineage metadata for a given label, enriched with genome metadata."""
        self._load_data()
        row = _match_lineage(lineage, self._lineages)
        if row is None:
            return {}

        lineage_id = _safe_str(row.get("lineage_id"))
        genome_count = 0
        earliest = latest = None
        derived_countries: set[str] = set()
        derived_regions: set[str] = set()
        derived_hosts: set[str] = set()

        if not self._genome_metadata.empty and "lineage_id" in self._genome_metadata.columns and lineage_id:
            mask = self._genome_metadata["lineage_id"].astype(str).str.strip().str.lower() == lineage_id.lower()
            meta_rows = self._genome_metadata[mask]
            genome_count = len(meta_rows)
            if genome_count and "collection_date" in meta_rows.columns:
                dates = pd.to_datetime(meta_rows["collection_date"], errors="coerce")
                earliest = dates.min() if not dates.empty else None
                latest = dates.max() if not dates.empty else None
            if "collection_country" in meta_rows.columns:
                for c in meta_rows["collection_country"].dropna().astype(str):
                    c = c.strip()
                    if c:
                        derived_countries.add(c)
            if "collection_region" in meta_rows.columns:
                for r in meta_rows["collection_region"].dropna().astype(str):
                    r = r.strip()
                    if r:
                        derived_regions.add(r)
            if "host" in meta_rows.columns:
                for h in meta_rows["host"].dropna().astype(str):
                    h = h.strip()
                    if h:
                        derived_hosts.add(h)

        first_curated = _safe_str(row.get("first_detected"))
        last_curated = _safe_str(row.get("last_detected"))
        first = first_curated or (str(earliest.date()) if pd.notna(earliest) else "")
        last = last_curated or (str(latest.date()) if pd.notna(latest) else "")

        curated_countries = _parse_pg_array(_safe_str(row.get("countries_reported")))
        curated_regions = _parse_pg_array(_safe_str(row.get("regions_reported")))
        curated_endemic = _parse_pg_array(_safe_str(row.get("endemic_regions")))

        # Merge curated and derived metadata, preferring curated values.
        countries_reported = curated_countries or sorted(derived_countries)
        regions_reported = curated_regions or sorted(derived_regions)
        # Endemic regions fall back to derived regions when absent.
        endemic_regions = curated_endemic or regions_reported

        primary_host = _safe_str(row.get("primary_host")) or (sorted(derived_hosts)[0] if derived_hosts else "")
        reservoir = _safe_str(row.get("reservoir"))

        return {
            "lineage_id": lineage_id,
            "lineage_name": _safe_str(row.get("lineage_name")),
            "parent_lineage": _safe_str(row.get("parent_lineage")),
            "first_detected": first,
            "last_detected": last,
            "total_genomes": genome_count,
            "countries_reported": countries_reported,
            "regions_reported": regions_reported,
            "endemic_regions": endemic_regions,
            "primary_host": primary_host,
            "reservoir": reservoir,
        }

    def load_tree_input(
        self,
        stage9: dict,
        tree_path: Optional[str] = None,
    ) -> TreeInput:
        """Load a phylogenetic tree and annotate tips, or build a metadata fallback."""
        self._load_data()

        sample_id = _safe_str(stage9.get("sample_id"))
        sample_accession = _safe_str(stage9.get("genome_accession"))
        sample_metadata = stage9.get("metadata", {})
        lineage = _safe_str(stage9.get("lineage"))

        # Tree metadata from the bioinformatics output if no CLI path given
        stage6_tree = stage9.get("stage6_phylogenetic_tree") or {}
        time_scaled = stage6_tree.get("time_scaled_tree") or {}
        bio_tree_path = stage6_tree.get("tree_file")
        effective_tree_path = tree_path or bio_tree_path

        matched_lineage = _match_lineage(lineage, self._lineages)
        lineage_id = _safe_str(matched_lineage.get("lineage_id")) if matched_lineage is not None else ""

        # Prefer a dedicated tips metadata file next to the tree when available,
        # and fall back to the genome metadata catalogue.
        metadata_df = self._genome_metadata
        tips_metadata_path = None
        if effective_tree_path:
            candidate = Path(effective_tree_path).parent / "tree.tips_metadata.csv"
            if candidate.exists():
                tips_metadata_path = str(candidate)
                try:
                    tips_df = pd.read_csv(candidate, low_memory=False)
                    if not tips_df.empty:
                        # If the genome metadata catalogue has overlapping accessions,
                        # the tips file takes precedence for tree-specific annotations.
                        if not self._genome_metadata.empty and "genome_accession" in self._genome_metadata.columns:
                            catalogue = self._genome_metadata.copy()
                            catalogue["_from_tips"] = 0
                            tips_df["_from_tips"] = 1
                            tips_df = tips_df.reindex(columns=catalogue.columns, fill_value="")
                            combined = pd.concat([catalogue, tips_df], ignore_index=True)
                            metadata_df = combined.sort_values("_from_tips", ascending=False).drop_duplicates(
                                subset=["genome_accession"], keep="first"
                            )
                            metadata_df = metadata_df.drop(columns=["_from_tips"])
                        else:
                            metadata_df = tips_df
                except Exception as exc:
                    log.warning("Could not load tree tips metadata %s: %s", candidate, exc)

        return load_tree_input(
            tree_path=effective_tree_path,
            metadata_df=metadata_df,
            sample_id=sample_id,
            sample_accession=sample_accession,
            sample_metadata=sample_metadata,
            time_scaled=time_scaled,
            lineage_id=lineage_id,
        )

    def analyze_phylogeography(self, stage9: dict, tree_input: TreeInput) -> AnalysisResult:
        """Run the phylogeographic origin and dissemination analyzer."""
        self._load_data()
        analyzer = PhylogeographicAnalyzer(
            tree_input=tree_input,
            lineages_df=self._lineages,
            genome_metadata_df=self._genome_metadata,
        )
        return analyzer.analyze(stage9)

    def analyze_genetic_relatedness(self, stage9: dict, tree_input: TreeInput) -> AnalysisResult:
        """Run the genetic relatedness analyzer.

        This deliberately avoids inferring transmission networks; it only reports
        pairwise genetic distances and temporal gaps to contextual genomes.
        """
        self._load_data()
        analyzer = GeneticRelatednessAnalyzer(
            tree_input=tree_input,
            genome_metadata_df=self._genome_metadata,
            lineages_df=self._lineages,
        )
        return analyzer.analyze(stage9)

    def analyze_molecular_clock(self, stage9: dict, tree_input: TreeInput) -> AnalysisResult:
        """Run the molecular clock timing analyzer."""
        analyzer = MolecularClockAnalyzer(tree_input=tree_input)
        return analyzer.analyze(stage9)

    def analyze_evidence_weighted_threat(
        self,
        stage9: dict,
        all_results: dict[str, AnalysisResult],
        risk: dict,
    ) -> AnalysisResult:
        """Run the integrated evidence-weighted threat assessment."""
        analyzer = EvidenceWeightedThreatAnalyzer(engine=self)
        return analyzer.analyze(stage9, all_results, risk)

    # ------------------------------------------------------------------
    # Phase 3: extended decision-oriented analyses
    # ------------------------------------------------------------------
    def analyze_epidemic_dynamics(self, stage9: dict) -> AnalysisResult:
        """Run the epidemic-dynamics analyzer."""
        self._load_data()
        # Lazy import avoids a circular dependency with extended_analyzers.
        from intelligence_engine.evidence_integration.analyzers.extended_analyzers import (
            EpidemicDynamicsAnalyzer,
        )

        analyzer = EpidemicDynamicsAnalyzer(
            epi_output=self._epi_output,
            lineages_df=self._lineages,
            genome_metadata_df=self._genome_metadata,
        )
        return analyzer.analyze(stage9)

    def analyze_genomic_signal(
        self,
        stage9: dict,
        variants: list[dict],
        matched_phenotypes: list[dict],
    ) -> AnalysisResult:
        """Run the genomic-signal analyzer."""
        self._load_data()
        from intelligence_engine.evidence_integration.analyzers.extended_analyzers import (
            GenomicSignalAnalyzer,
        )

        analyzer = GenomicSignalAnalyzer(
            associations_df=self._associations,
            protein_variants_df=self._protein_variants,
            lineages_df=self._lineages,
            genome_metadata_df=self._genome_metadata,
        )
        return analyzer.analyze(stage9, variants, matched_phenotypes)

    def analyze_surveillance(self, stage9: dict) -> AnalysisResult:
        """Run the surveillance and sampling analyzer."""
        self._load_data()
        from intelligence_engine.evidence_integration.analyzers.extended_analyzers import (
            SurveillanceAnalyzer,
        )

        analyzer = SurveillanceAnalyzer(
            lineages_df=self._lineages,
            genome_metadata_df=self._genome_metadata,
        )
        return analyzer.analyze(stage9)

    def analyze_lineage_behavior(self, stage9: dict) -> AnalysisResult:
        """Run the lineage-behaviour analyzer."""
        self._load_data()
        from intelligence_engine.evidence_integration.analyzers.extended_analyzers import (
            LineageBehaviorAnalyzer,
        )

        analyzer = LineageBehaviorAnalyzer(
            lineages_df=self._lineages,
            genome_metadata_df=self._genome_metadata,
        )
        return analyzer.analyze(stage9)

    def analyze_integrated_inference(
        self,
        stage9: dict,
        tree_input: TreeInput,
        all_results: dict[str, AnalysisResult],
    ) -> AnalysisResult:
        """Run the integrated epi+genomic inference analyzer."""
        from intelligence_engine.evidence_integration.analyzers.extended_analyzers import (
            IntegratedInferenceAnalyzer,
        )

        analyzer = IntegratedInferenceAnalyzer(epi_output=self._epi_output or {})
        return analyzer.analyze(stage9, tree_input, all_results)

    def analyze_evidence_consistency(
        self,
        stage9: dict,
        all_results: dict[str, AnalysisResult],
    ) -> AnalysisResult:
        """Run the evidence consistency analyzer."""
        self._load_data()
        from intelligence_engine.evidence_integration.analyzers.extended_analyzers import (
            EvidenceConsistencyAnalyzer,
        )

        analyzer = EvidenceConsistencyAnalyzer(
            epi_output=self._epi_output or {},
            lineages_df=self._lineages,
            genome_metadata_df=self._genome_metadata,
        )
        return analyzer.analyze(stage9, all_results)

    def analyze_knowledge_gaps(
        self,
        stage9: dict,
        all_results: dict[str, AnalysisResult],
    ) -> AnalysisResult:
        """Run the consolidated knowledge-gap analyzer."""
        from intelligence_engine.evidence_integration.analyzers.extended_analyzers import (
            KnowledgeGapAnalyzer,
        )

        analyzer = KnowledgeGapAnalyzer(epi_output=self._epi_output or {})
        return analyzer.analyze(stage9, all_results)

    def analyze_comparative_outbreak(self, stage9: dict) -> AnalysisResult:
        """Run the comparative outbreak analyzer."""
        self._load_data()
        from intelligence_engine.evidence_integration.analyzers.extended_analyzers import (
            ComparativeOutbreakAnalyzer,
        )

        analyzer = ComparativeOutbreakAnalyzer(
            epi_output=self._epi_output or {},
            lineages_df=self._lineages,
            genome_metadata_df=self._genome_metadata,
        )
        return analyzer.analyze(stage9)
