"""
evidence_integration/cross_evidence.py — Cross-evidence statistical analyses.

Each analyzer below produces quantitative summaries, statistical
associations, and confidence measures from the harmonized evidence objects
and curated reference DataFrames. None of them generate public-health
conclusions or risk interpretations -- that is the responsibility of a
downstream Genomic Intelligence Engine stage. Output here is deliberately
restricted to numbers, statistics, and their supporting evidence.
"""

from __future__ import annotations

import itertools
import logging
from dataclasses import asdict, dataclass, field
from typing import Any, Optional

import pandas as pd
from scipy import stats as scipy_stats

from intelligence_engine.evidence_integration.harmonization import EvidenceObject

log = logging.getLogger(__name__)


@dataclass
class StatisticalFinding:
    """A single quantitative/statistical result, fully traceable to its inputs."""

    metric: str
    value: Any
    method: str = ""
    sample_size: Optional[int] = None
    p_value: Optional[float] = None
    supporting_evidence: list[str] = field(default_factory=list)
    notes: str = ""


@dataclass
class CrossEvidenceResult:
    """Result of one cross-evidence analyzer."""

    title: str
    findings: list[StatisticalFinding] = field(default_factory=list)
    data_gaps: list[str] = field(default_factory=list)
    metrics: dict = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


# ---------------------------------------------------------------------------
# 1. Mutation co-occurrence
# ---------------------------------------------------------------------------

class MutationCooccurrenceAnalyzer:
    """Quantify which detected mutations co-occur within the same genome, and
    report each mutation's independent historical frequency for context."""

    def analyze(self, variants: list[dict]) -> CrossEvidenceResult:
        findings: list[StatisticalFinding] = []
        gaps: list[str] = []
        hgvs_list = [v.get("hgvs_p") for v in variants if v.get("hgvs_p")]

        if len(hgvs_list) >= 2:
            for a, b in itertools.combinations(hgvs_list, 2):
                findings.append(
                    StatisticalFinding(
                        metric="co_occurring_mutation_pair",
                        value={"pair": [a, b]},
                        method="within-genome co-occurrence (single sample)",
                        sample_size=1,
                        supporting_evidence=[a, b],
                        notes="Both mutations detected in the same genome. Historical "
                              "co-occurrence frequency across other genomes is not "
                              "computable: the curated database stores aggregate "
                              "per-mutation statistics, not per-genome mutation sets.",
                    )
                )
        elif len(hgvs_list) == 1:
            gaps.append("Only one mutation detected; no co-occurrence pairs to evaluate.")
        else:
            gaps.append("No mutations detected; co-occurrence analysis not applicable.")

        for v in variants:
            ctx = v.get("_curated_context") or {}
            if ctx.get("frequency"):
                findings.append(
                    StatisticalFinding(
                        metric="mutation_historical_frequency",
                        value=ctx.get("frequency"),
                        method="curated protein_variants aggregate (genome_count / species_total_genomes)",
                        supporting_evidence=[str(v.get("hgvs_p"))],
                        notes=f"first_seen={ctx.get('first_seen_date')}, last_seen={ctx.get('last_seen_date')}, "
                              f"countries_seen={ctx.get('countries_seen')}",
                    )
                )

        return CrossEvidenceResult(
            title="mutation_cooccurrence",
            findings=findings,
            data_gaps=gaps,
            metrics={"n_mutations_detected": len(hgvs_list), "n_pairs": max(0, len(hgvs_list) * (len(hgvs_list) - 1) // 2)},
        )


# ---------------------------------------------------------------------------
# 2. Lineage-phenotype statistical association
# ---------------------------------------------------------------------------

class LineagePhenotypeAssociationAnalyzer:
    """Test whether specific phenotype categories are statistically
    associated with the sample's lineage vs. other lineages in the curated
    genotype_phenotype table, using Fisher's exact test (2x2) or a
    Chi-square test of independence (contingency table)."""

    def __init__(self, associations_df: pd.DataFrame):
        self.associations = associations_df if associations_df is not None else pd.DataFrame()

    def analyze(self, lineage_id: Optional[str]) -> CrossEvidenceResult:
        findings: list[StatisticalFinding] = []
        gaps: list[str] = []
        df = self.associations

        if df.empty or not lineage_id or "lineage_id" not in df.columns or "phenotype_category" not in df.columns:
            gaps.append(
                "Insufficient genotype_phenotype rows (need lineage_id and phenotype_category "
                "across multiple lineages) to run a statistical association test."
            )
            return CrossEvidenceResult(title="lineage_phenotype_association", findings=findings, data_gaps=gaps)

        categories = [c for c in df["phenotype_category"].dropna().unique()]
        n_lineages = df["lineage_id"].dropna().nunique()

        for category in categories:
            table = pd.crosstab(df["lineage_id"] == lineage_id, df["phenotype_category"] == category)
            if table.shape != (2, 2):
                continue
            a, b = table.iloc[1, 1], table.iloc[1, 0]
            c, d = table.iloc[0, 1], table.iloc[0, 0]
            n = a + b + c + d
            if n < 5 or n_lineages < 2:
                continue
            try:
                if n < 30 or min(a, b, c, d) < 5:
                    odds_ratio, p_value = scipy_stats.fisher_exact([[a, b], [c, d]])
                    method = "Fisher's exact test"
                    stat_value = odds_ratio
                else:
                    chi2, p_value, _, _ = scipy_stats.chi2_contingency([[a, b], [c, d]])
                    method = "Chi-square test of independence"
                    stat_value = chi2
            except (ValueError, ZeroDivisionError):
                continue

            findings.append(
                StatisticalFinding(
                    metric=f"lineage_association::{category}",
                    value=stat_value,
                    method=method,
                    sample_size=int(n),
                    p_value=float(p_value),
                    supporting_evidence=[f"lineage_id={lineage_id}", f"phenotype_category={category}"],
                    notes=f"2x2 contingency: lineage-matches x category-matches = [[{a},{b}],[{c},{d}]]",
                )
            )

        if not findings and not gaps:
            gaps.append(f"No phenotype categories had a valid 2x2 contingency table for lineage {lineage_id}.")

        return CrossEvidenceResult(
            title="lineage_phenotype_association",
            findings=findings,
            data_gaps=gaps,
            metrics={"n_lineages_compared": int(n_lineages), "n_categories_tested": len(findings)},
        )


# ---------------------------------------------------------------------------
# 3. Temporal trend analysis
# ---------------------------------------------------------------------------

class TemporalTrendAnalyzer:
    """Fit a Poisson GLM of genome counts (or outbreak case counts) over time
    to quantify whether a lineage's sequencing/reporting rate is increasing,
    decreasing, or stable."""

    def analyze(self, genome_metadata_df: pd.DataFrame, outbreaks: list[dict]) -> CrossEvidenceResult:
        findings: list[StatisticalFinding] = []
        gaps: list[str] = []

        df = genome_metadata_df
        if df is not None and not df.empty and "collection_year" in df.columns:
            yearly = (
                df.dropna(subset=["collection_year"])
                .groupby("collection_year").size()
                .reset_index(name="genome_count")
            )
            findings_added = self._fit_trend(
                yearly, x_col="collection_year", y_col="genome_count",
                metric="genome_collection_trend", unit="genomes/year",
            )
            findings.extend(findings_added)
            if not findings_added:
                gaps.append("Not enough distinct years of genome_metadata to fit a temporal trend (need >= 3).")
        else:
            gaps.append("No genome_metadata with collection_year available for temporal trend analysis.")

        outbreak_rows = [
            {"year": int(str(ob.get("start_date"))[:4]), "cases": ob.get("cases")}
            for ob in outbreaks
            if ob.get("start_date") and str(ob.get("start_date"))[:4].isdigit() and ob.get("cases") is not None
        ]
        if len(outbreak_rows) >= 3:
            odf = pd.DataFrame(outbreak_rows).groupby("year")["cases"].sum().reset_index()
            findings.extend(
                self._fit_trend(odf, x_col="year", y_col="cases", metric="outbreak_case_trend", unit="cases/year")
            )
        else:
            gaps.append("Fewer than 3 dated outbreak records with case counts; cannot fit an outbreak case trend.")

        return CrossEvidenceResult(
            title="temporal_trend",
            findings=findings,
            data_gaps=gaps,
            metrics={"n_trend_findings": len(findings)},
        )

    def _fit_trend(self, df: pd.DataFrame, x_col: str, y_col: str, metric: str, unit: str) -> list[StatisticalFinding]:
        if df.shape[0] < 3:
            return []
        try:
            import statsmodels.api as sm

            x = sm.add_constant(df[x_col].astype(float))
            y = df[y_col].astype(float)
            model = sm.GLM(y, x, family=sm.families.Poisson()).fit()
            slope = float(model.params[x_col])
            p_value = float(model.pvalues[x_col])
            direction = "increasing" if slope > 0 else ("decreasing" if slope < 0 else "flat")
            return [
                StatisticalFinding(
                    metric=metric,
                    value={"slope": slope, "direction": direction, "unit": unit},
                    method="Poisson GLM (statsmodels)",
                    sample_size=int(df.shape[0]),
                    p_value=p_value,
                    supporting_evidence=[f"{x_col}={int(r[x_col])}:{y_col}={int(r[y_col])}" for _, r in df.iterrows()],
                )
            ]
        except Exception as e:  # noqa: BLE001 - degrade gracefully, this is a supporting stat not critical path
            log.warning(f"GLM trend fit failed for {metric}: {e}")
            return []


# ---------------------------------------------------------------------------
# 4. Geographic distribution
# ---------------------------------------------------------------------------

class GeographicDistributionAnalyzer:
    """Summarize the geographic distribution of curated genomes and reported
    outbreaks/surveillance rows for the lineage's pathogen."""

    def analyze(self, genome_metadata_df: pd.DataFrame, outbreaks: list[dict], surveillance: list[dict]) -> CrossEvidenceResult:
        findings: list[StatisticalFinding] = []
        gaps: list[str] = []

        df = genome_metadata_df
        if df is not None and not df.empty and "collection_country" in df.columns:
            counts = df["collection_country"].dropna().value_counts()
            total = int(counts.sum())
            for country, n in counts.items():
                findings.append(
                    StatisticalFinding(
                        metric="genome_country_distribution",
                        value={"country": country, "count": int(n), "pct_of_genomes": round(100 * n / total, 1)},
                        method="curated genome_metadata frequency count",
                        sample_size=total,
                    )
                )
        else:
            gaps.append("No genome_metadata with collection_country available.")

        outbreak_countries = pd.Series([ob.get("country") for ob in outbreaks if ob.get("country")])
        if not outbreak_countries.empty:
            ob_counts = outbreak_countries.value_counts()
            for country, n in ob_counts.items():
                findings.append(
                    StatisticalFinding(
                        metric="outbreak_country_distribution",
                        value={"country": country, "n_reports": int(n)},
                        method="fetched outbreak report frequency count",
                        sample_size=int(ob_counts.sum()),
                    )
                )
        else:
            gaps.append("No outbreak reports with a country field to summarize.")

        surveillance_countries = [s.get("country") for s in surveillance if s.get("country")]
        n_surveilled = len(set(surveillance_countries))

        return CrossEvidenceResult(
            title="geographic_distribution",
            findings=findings,
            data_gaps=gaps,
            metrics={
                "n_countries_with_genomes": int(df["collection_country"].nunique()) if (df is not None and not df.empty and "collection_country" in df.columns) else 0,
                "n_countries_with_outbreak_reports": int(outbreak_countries.nunique()) if not outbreak_countries.empty else 0,
                "n_countries_under_surveillance": n_surveilled,
            },
        )


# ---------------------------------------------------------------------------
# 5. Mutation persistence
# ---------------------------------------------------------------------------

class MutationPersistenceAnalyzer:
    """Quantify how long each detected mutation has been observed in the
    curated record (first_seen -> last_seen span)."""

    def analyze(self, variants: list[dict]) -> CrossEvidenceResult:
        findings: list[StatisticalFinding] = []
        gaps: list[str] = []

        for v in variants:
            ctx = v.get("_curated_context") or {}
            first_seen = ctx.get("first_seen_date")
            last_seen = ctx.get("last_seen_date")
            if not first_seen or not last_seen:
                gaps.append(f"No first/last-seen dates in curated data for {v.get('hgvs_p')}.")
                continue
            try:
                span_days = (pd.to_datetime(last_seen) - pd.to_datetime(first_seen)).days
            except (ValueError, TypeError):
                continue
            findings.append(
                StatisticalFinding(
                    metric="mutation_persistence_span",
                    value={"first_seen": str(first_seen), "last_seen": str(last_seen), "span_days": span_days},
                    method="curated protein_variants first_seen_date/last_seen_date delta",
                    supporting_evidence=[str(v.get("hgvs_p"))],
                )
            )

        return CrossEvidenceResult(
            title="mutation_persistence",
            findings=findings,
            data_gaps=gaps,
            metrics={"n_variants_with_persistence_data": len(findings)},
        )


# ---------------------------------------------------------------------------
# 6. Intervention association
# ---------------------------------------------------------------------------

class InterventionAssociationAnalyzer:
    """Quantify how many phenotype associations fall into
    vaccine/diagnostic/therapeutic-relevant categories, purely as a count
    summary (no effectiveness/impact interpretation)."""

    INTERVENTION_CATEGORIES = {
        "vaccine_escape", "vaccine_effectiveness", "drug_resistance",
        "diagnostic_escape", "immune_escape",
    }

    def analyze(self, matched_phenotypes: list[dict]) -> CrossEvidenceResult:
        findings: list[StatisticalFinding] = []
        gaps: list[str] = []

        if not matched_phenotypes:
            gaps.append("No matched genotype-phenotype associations to evaluate for intervention relevance.")
            return CrossEvidenceResult(title="intervention_association", findings=findings, data_gaps=gaps)

        counts: dict[str, int] = {}
        for p in matched_phenotypes:
            cat = str(p.get("phenotype_category", "")).lower()
            if cat in self.INTERVENTION_CATEGORIES:
                counts[cat] = counts.get(cat, 0) + 1

        for cat, n in counts.items():
            findings.append(
                StatisticalFinding(
                    metric="intervention_relevant_association_count",
                    value={"category": cat, "count": n},
                    method="count of matched genotype-phenotype rows by category",
                    sample_size=len(matched_phenotypes),
                )
            )
        if not counts:
            gaps.append("No matched associations fall into an intervention-relevant phenotype category.")

        return CrossEvidenceResult(
            title="intervention_association",
            findings=findings,
            data_gaps=gaps,
            metrics={"n_intervention_relevant_associations": sum(counts.values()), "n_total_associations": len(matched_phenotypes)},
        )


# ---------------------------------------------------------------------------
# 7. Confidence scoring
# ---------------------------------------------------------------------------

class ConfidenceScoringAnalyzer:
    """Score each EvidenceObject's confidence from 0-1 based purely on the
    quantity and quality of its supporting evidence (source count, curated
    evidence_strength, sample-size/genome_count support, conflict flags).
    No biological or public-health interpretation."""

    STRENGTH_WEIGHT = {"strong": 1.0, "moderate": 0.6, "weak": 0.3}

    def analyze(self, evidence_objects: list[EvidenceObject]) -> CrossEvidenceResult:
        findings: list[StatisticalFinding] = []

        for obj in evidence_objects:
            score = 0.0
            components: dict[str, float] = {}

            n_sources = len(set(obj.sources))
            source_component = min(n_sources / 3.0, 1.0) * 0.3
            components["source_diversity"] = round(source_component, 3)
            score += source_component

            strengths = [
                self.STRENGTH_WEIGHT.get(str(p.get("evidence_strength", "")).lower(), 0.0)
                for p in obj.phenotype_associations
            ]
            strength_component = (max(strengths) if strengths else 0.0) * 0.4
            components["max_evidence_strength"] = round(strength_component, 3)
            score += strength_component

            genome_count = 0
            if obj.variant and obj.variant.get("curated_context"):
                genome_count = obj.variant["curated_context"].get("genome_count") or 0
            sample_component = min(genome_count / 50.0, 1.0) * 0.2
            components["sample_size_support"] = round(sample_component, 3)
            score += sample_component

            flagged = any(p.get("record_flagged") for p in obj.phenotype_associations)
            conflict_penalty = 0.1 if flagged else 0.0
            components["conflict_penalty"] = -conflict_penalty
            score -= conflict_penalty

            consistency_component = 0.1 if (obj.historical_outbreaks or obj.molecular_epidemiology) else 0.0
            components["cross_source_epi_context"] = consistency_component
            score += consistency_component

            score = max(0.0, min(1.0, score))
            findings.append(
                StatisticalFinding(
                    metric="evidence_confidence_score",
                    value={"key": obj.key, "level": obj.level, "score": round(score, 3)},
                    method="weighted sum: source_diversity(0.3) + max_evidence_strength(0.4) + "
                           "sample_size_support(0.2) + cross_source_epi_context(0.1) - conflict_penalty(0.1)",
                    supporting_evidence=list(set(obj.sources)),
                    notes=str(components),
                )
            )

        return CrossEvidenceResult(
            title="confidence_scoring",
            findings=findings,
            metrics={"n_evidence_objects_scored": len(findings)},
        )
