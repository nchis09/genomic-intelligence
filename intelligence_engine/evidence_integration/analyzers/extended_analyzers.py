"""Extended decision-oriented analyzers for the Genomic Intelligence Engine.

Each analyzer consumes the curated reference tables, the epidemiological output,
and the bioinformatics stage9 object. They return structured AnalysisResult
objects that feed the evidence-weighted threat assessment and the R figures.
"""

from __future__ import annotations

from collections import Counter
from datetime import datetime
from typing import Any, Optional

import pandas as pd

from intelligence_engine.evidence_integration.engine import (  # noqa: E402
    AnalysisResult,
    EvidenceRecord,
    _match_lineage,
    _parse_date,
    _parse_pg_array,
    _safe_str,
)


def _to_float(value: Any) -> Optional[float]:
    try:
        return float(value)
    except (ValueError, TypeError):
        return None


def _to_date(value: Any) -> Optional[datetime]:
    s = _parse_date(value)
    if not s:
        return None
    try:
        return pd.to_datetime(s)
    except (ValueError, TypeError):
        return None


def _days_between(a: Any, b: Any) -> Optional[int]:
    da = _to_date(a)
    db = _to_date(b)
    if da is None or db is None:
        return None
    return abs((da - db).days)


# ---------------------------------------------------------------------------
# 1. Epidemic dynamics
# ---------------------------------------------------------------------------
class EpidemicDynamicsAnalyzer:
    """Compute epidemic dynamics metrics from the epi output and curated data."""

    def __init__(
        self,
        epi_output: dict,
        lineages_df: pd.DataFrame,
        genome_metadata_df: pd.DataFrame,
    ):
        self.epi = epi_output or {}
        self.lineages = lineages_df
        self.genome_metadata = genome_metadata_df

    def analyze(self, stage9: dict) -> AnalysisResult:
        findings: list[EvidenceRecord] = []
        gaps: list[str] = []
        metrics: dict[str, Any] = {}

        outbreaks = self.epi.get("outbreaks", [])
        if not outbreaks:
            gaps.append("No outbreak records available in epi output.")

        # Parse outbreaks. OutbreakRecord.cfr is a numeric percentage (0-100);
        # normalise to a fraction (0-1) to match downstream crude_cfr math.
        records = []
        for o in outbreaks:
            cases = _to_float(o.get("cases"))
            deaths = _to_float(o.get("deaths"))
            cfr_pct = _to_float(o.get("cfr"))
            cfr = cfr_pct / 100 if cfr_pct is not None else None
            start_date = _safe_str(o.get("start_date"))
            end_date = _safe_str(o.get("end_date"))
            duration = None
            if start_date and end_date:
                sd = pd.to_datetime(start_date, errors="coerce")
                ed = pd.to_datetime(end_date, errors="coerce")
                if pd.notna(sd) and pd.notna(ed):
                    duration = (ed - sd).days
            records.append(
                {
                    "date": start_date,
                    "country": _safe_str(o.get("country")),
                    "cases": cases,
                    "deaths": deaths,
                    "cfr": cfr,
                    "duration": duration,
                }
            )

        df = pd.DataFrame(records)

        # Basic totals (guard against empty outbreak table)
        if "cases" in df.columns:
            total_cases = df["cases"].sum(skipna=True)
            total_deaths = df["deaths"].sum(skipna=True)
            metrics["total_cases"] = int(total_cases) if pd.notna(total_cases) else 0
            metrics["total_deaths"] = int(total_deaths) if pd.notna(total_deaths) else 0
            if metrics["total_cases"] > 0:
                metrics["crude_cfr"] = round(metrics["total_deaths"] / metrics["total_cases"], 3)
            else:
                metrics["crude_cfr"] = None

            # Outbreak size / duration
            size_rows = df.dropna(subset=["cases"])
            duration_rows = df.dropna(subset=["duration"])
            metrics["outbreaks_with_case_data"] = len(size_rows)
            metrics["outbreaks_with_duration_data"] = len(duration_rows)
            if not size_rows.empty:
                metrics["median_outbreak_size"] = round(size_rows["cases"].median(), 1)
                metrics["max_outbreak_size"] = int(size_rows["cases"].max())
            if not duration_rows.empty:
                metrics["median_outbreak_duration_days"] = round(duration_rows["duration"].median(), 1)

            # Attack rate proxy: cases / outbreak if duration is known
            attack_rows = df.dropna(subset=["cases", "duration"])
            if not attack_rows.empty:
                metrics["median_daily_case_rate"] = round((attack_rows["cases"] / attack_rows["duration"]).median(), 2)
        else:
            metrics["total_cases"] = 0
            metrics["total_deaths"] = 0
            metrics["crude_cfr"] = None
            metrics["outbreaks_with_case_data"] = 0
            metrics["outbreaks_with_duration_data"] = 0

        # Countries affected
        countries = {r["country"] for r in records if r["country"]}
        metrics["countries_with_reported_outbreaks"] = sorted(countries)

        # Transmission parameters (R0, serial interval) if present.
        # 'transmission' is now a single pathogen-level object, not a list.
        transmission = self.epi.get("transmission") or {}
        r0_low = _to_float(transmission.get("r0_low"))
        r0_high = _to_float(transmission.get("r0_high"))
        r0 = None
        if r0_low is not None and r0_high is not None:
            r0 = round((r0_low + r0_high) / 2, 2)
        elif r0_low is not None:
            r0 = r0_low
        elif r0_high is not None:
            r0 = r0_high
        serial_interval = _to_float(transmission.get("serial_interval_days"))

        if r0 is not None:
            metrics["effective_reproduction_number"] = r0
            findings.append(
                EvidenceRecord(
                    assertion=f"Reported reproduction number/R0 estimate is {r0}.",
                    source="transmission",
                    confidence="medium",
                    biological_relevance="transmissibility",
                    epidemiological_relevance="potential for sustained transmission",
                    finding_type="r0_estimate",
                )
            )
        else:
            gaps.append("No R0/Rt estimate found in epi output.")

        if serial_interval is not None:
            metrics["serial_interval_days"] = serial_interval
            findings.append(
                EvidenceRecord(
                    assertion=f"Reported serial interval is {serial_interval} days.",
                    source="transmission",
                    confidence="medium",
                    biological_relevance="generation time",
                    epidemiological_relevance="contact tracing window",
                    finding_type="serial_interval",
                )
            )
        else:
            gaps.append("No serial interval found in epi output.")

        # Doubling time proxy: if we have dates and cumulative cases, estimate growth
        if not df.empty and all(c in df.columns for c in ["date", "cases"]):
            dated = df.dropna(subset=["date", "cases"]).sort_values("date")
            if len(dated) >= 2:
                dated["date_dt"] = pd.to_datetime(dated["date"], errors="coerce")
                dated = dated.dropna(subset=["date_dt"])
                if len(dated) >= 2:
                    dated["cum_cases"] = dated["cases"].cumsum()
                    first = dated.iloc[0]
                    last = dated.iloc[-1]
                    days = (last["date_dt"] - first["date_dt"]).days
                    if days > 0 and first["cum_cases"] > 0 and last["cum_cases"] > 0:
                        growth_rate = (last["cum_cases"] / first["cum_cases"]) ** (1 / days) - 1
                        metrics["cumulative_growth_rate_per_day"] = round(growth_rate, 4)
                        if growth_rate > 0:
                            metrics["doubling_time_days"] = round(0.693 / growth_rate, 1)

        # CFR evidence
        if metrics.get("crude_cfr") is not None:
            if metrics["crude_cfr"] >= 0.5:
                implication = "very high mortality; prepare clinical surge and burial capacity"
            elif metrics["crude_cfr"] >= 0.3:
                implication = "high mortality; prioritise supportive care and contact tracing"
            else:
                implication = "lower observed mortality; still investigate for under-ascertainment"
            findings.append(
                EvidenceRecord(
                    assertion=f"Crude case fatality risk across reported outbreaks is {metrics['crude_cfr']:.1%} ({metrics['total_deaths']} deaths / {metrics['total_cases']} cases).",
                    source="outbreaks",
                    confidence="low" if metrics["outbreaks_with_case_data"] < 5 else "medium",
                    biological_relevance="clinical severity",
                    epidemiological_relevance="burden estimate",
                    finding_type="cfr_estimate",
                )
            )

        summary = (
            f"Epidemic dynamics: {metrics.get('total_cases', 0)} cases and "
            f"{metrics.get('total_deaths', 0)} deaths across {len(records)} records; "
            f"crude CFR {metrics.get('crude_cfr') or 'unknown'}."
        )
        return AnalysisResult(
            title="epidemic_dynamics",
            findings=findings,
            gaps=gaps,
            summary=summary,
            metrics=metrics,
        )


# ---------------------------------------------------------------------------
# 2. Genomic signal
# ---------------------------------------------------------------------------
class GenomicSignalAnalyzer:
    """Quantify genomic novelty, mutation burden, phenotype associations and enrichment."""

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

    def analyze(self, stage9: dict, variants: list[dict], matched_phenotypes: list[dict]) -> AnalysisResult:
        findings: list[EvidenceRecord] = []
        gaps: list[str] = []
        metrics: dict[str, Any] = {}

        total_variants = len(variants)
        metrics["total_amino_acid_variants"] = total_variants

        # Normalise protein variants table
        pv = self.protein_variants.copy()
        if not pv.empty:
            pv["gene_norm"] = pv["gene"].astype(str).str.upper().str.strip()
            pv["position"] = pd.to_numeric(pv["position"], errors="coerce")
            pv["alt_aa_norm"] = pv["alt_aa"].astype(str).str.upper().str.strip()

        # Determine known vs novel variants
        novel_count = 0
        known_count = 0
        gene_counts = Counter()
        hotspot_count = 0
        phenotype_associated = 0

        for v in variants:
            gene = str(v.get("gene", "")).strip().upper()
            pos = _to_float(v.get("position"))
            alt = str(v.get("alt_aa", "")).strip().upper()
            ref = str(v.get("ref_aa", "")).strip().upper()
            hgvs = _safe_str(v.get("hgvs_p"))
            domain = _safe_str(v.get("domain"))
            gene_counts[gene] += 1

            if domain and "hotspot" in domain.lower():
                hotspot_count += 1

            if not gene or pos is None or not alt or ref == alt:
                continue

            if pv.empty:
                novel_count += 1
                continue

            mask = (
                (pv["gene_norm"] == gene)
                & (pv["position"] == pos)
                & (pv["alt_aa_norm"] == alt)
            )
            if pv[mask].empty:
                novel_count += 1
            else:
                known_count += 1

        metrics["known_variants_in_catalogue"] = known_count
        metrics["novel_variants_not_in_catalogue"] = novel_count
        if total_variants:
            metrics["genomic_novelty_score"] = round(novel_count / total_variants, 3)
        else:
            metrics["genomic_novelty_score"] = 0.0

        # Mutation burden vs lineage average from genome metadata
        matched = _match_lineage(_safe_str(stage9.get("lineage")), self.lineages)
        lineage_id = _safe_str(matched.get("lineage_id")) if matched is not None else ""
        metrics["lineage_average_variants"] = None
        metrics["mutation_burden_difference"] = None
        if lineage_id and not self.genome_metadata.empty:
            # Heuristic: use genome length as a proxy if no per-genome variant count exists
            lineage_meta = self.genome_metadata[self.genome_metadata["lineage_id"].astype(str) == lineage_id]
            if not lineage_meta.empty and "genome_length" in lineage_meta.columns:
                avg_length = pd.to_numeric(lineage_meta["genome_length"], errors="coerce").mean()
                if pd.notna(avg_length) and avg_length > 0:
                    # Very rough: assume one variant per ~19,000 nt for a typical genome
                    expected = avg_length / 19000
                    metrics["lineage_average_variants"] = round(expected, 2)
                    metrics["mutation_burden_difference"] = round(total_variants - expected, 2)

        # Gene-level diversity vs protein_variants catalog
        gene_diversity: dict[str, Any] = {}
        if not pv.empty:
            for gene, count in gene_counts.items():
                catalog_positions = pv[pv["gene_norm"] == gene]["position"].nunique()
                gene_diversity[gene] = {
                    "sample_variants": count,
                    "catalog_positions": int(catalog_positions),
                    "density_ratio": round(count / catalog_positions, 3) if catalog_positions else None,
                }
        metrics["gene_level_diversity"] = gene_diversity

        # Phenotype burden by category
        phenotype_burden: dict[str, int] = Counter()
        for mp in matched_phenotypes:
            cat = _safe_str(mp.get("phenotype_category")).lower().replace(" ", "_")
            if cat:
                phenotype_burden[cat] += 1
        metrics["phenotype_burden"] = dict(phenotype_burden)

        # Mutation enrichment in hotspots / functional domains
        if total_variants:
            metrics["hotspot_ratio"] = round(hotspot_count / total_variants, 3)
        else:
            metrics["hotspot_ratio"] = 0.0

        # Build findings
        if metrics["genomic_novelty_score"] and metrics["genomic_novelty_score"] > 0.3:
            findings.append(
                EvidenceRecord(
                    assertion=f"{int(metrics['genomic_novelty_score'] * 100)}% of detected variants are not in the curated protein-variants catalogue.",
                    source="protein_variants.csv",
                    confidence="medium",
                    biological_relevance="genetic novelty",
                    epidemiological_relevance="may represent lineage evolution or unusual variants",
                    finding_type="mutation_profile_uncharacterized",
                )
            )

        if phenotype_burden:
            top_cat = max(phenotype_burden, key=phenotype_burden.get)
            findings.append(
                EvidenceRecord(
                    assertion=f"Detected variants carry {sum(phenotype_burden.values())} known phenotype associations, most commonly {top_cat.replace('_', ' ')}.",
                    source="genotype_phenotype.csv",
                    confidence="medium",
                    biological_relevance="known functional impact",
                    epidemiological_relevance="may affect transmission, severity, or interventions",
                    finding_type="known_phenotype",
                )
            )

        if not findings:
            findings.append(
                EvidenceRecord(
                    assertion="Genomic signal analysis found no strongly divergent or phenotype-associated variants.",
                    source="protein_variants.csv + genotype_phenotype.csv",
                    confidence="medium",
                    biological_relevance="genomic typicality",
                    epidemiological_relevance="limited novel biological signal",
                    finding_type="mutation_profile_typical",
                )
            )

        summary = (
            f"Genomic signal: {total_variants} variants, "
            f"{metrics['genomic_novelty_score']:.1%} novel, "
            f"{sum(phenotype_burden.values())} phenotype associations."
        )
        return AnalysisResult(
            title="genomic_signal",
            findings=findings,
            gaps=gaps,
            summary=summary,
            metrics=metrics,
        )


# ---------------------------------------------------------------------------
# 3. Surveillance and sampling
# ---------------------------------------------------------------------------
class SurveillanceAnalyzer:
    """Assess surveillance gaps, sequencing intensity and geographic coverage."""

    def __init__(
        self,
        lineages_df: pd.DataFrame,
        genome_metadata_df: pd.DataFrame,
    ):
        self.lineages = lineages_df
        self.genome_metadata = genome_metadata_df

    def analyze(self, stage9: dict) -> AnalysisResult:
        findings: list[EvidenceRecord] = []
        gaps: list[str] = []
        metrics: dict[str, Any] = {}

        sample_country = _safe_str(stage9.get("metadata", {}).get("country"))
        sample_date = _parse_date(stage9.get("metadata", {}).get("collection_date"))
        lineage = _safe_str(stage9.get("lineage"))
        matched = _match_lineage(lineage, self.lineages)
        lineage_id = _safe_str(matched.get("lineage_id")) if matched is not None else ""

        if not lineage_id or self.genome_metadata.empty:
            gaps.append("Lineage or genome metadata unavailable for surveillance analysis.")
            return AnalysisResult(
                title="surveillance_metrics",
                findings=findings,
                gaps=gaps,
                summary="Surveillance metrics: insufficient metadata.",
                metrics=metrics,
            )

        lineage_meta = self.genome_metadata[self.genome_metadata["lineage_id"].astype(str) == lineage_id].copy()
        lineage_meta["collection_date_dt"] = pd.to_datetime(lineage_meta["collection_date"], errors="coerce")
        lineage_meta["collection_year"] = pd.to_numeric(
            lineage_meta["collection_date"].astype(str).str[:4], errors="coerce"
        )

        # Surveillance gap: time since last sequenced case in same country/lineage
        country_meta = lineage_meta[lineage_meta["collection_country"].astype(str).str.strip().str.lower() == sample_country.lower()]
        if not country_meta.empty:
            latest_country = country_meta["collection_date_dt"].max()
            metrics["last_sequenced_case_country"] = str(latest_country.date()) if pd.notna(latest_country) else None
            days_gap = _days_between(sample_date, latest_country)
            if days_gap is not None:
                metrics["surveillance_gap_country_days"] = int(days_gap)
                if days_gap > 365 * 2:
                    findings.append(
                        EvidenceRecord(
                            assertion=f"No sequenced cases from {sample_country} in this lineage for {int(days_gap / 365)} years.",
                            source="genome_metadata.csv",
                            confidence="medium",
                            biological_relevance="re-emergence or under-sampling",
                            epidemiological_relevance="surveillance gap",
                            finding_type="re-emergence_after_a_reporting_gap",
                        )
                    )

        # Global surveillance gap
        if not lineage_meta.empty:
            latest_global = lineage_meta["collection_date_dt"].max()
            metrics["last_sequenced_case_global"] = str(latest_global.date()) if pd.notna(latest_global) else None
            days_global = _days_between(sample_date, latest_global)
            if days_global is not None:
                metrics["surveillance_gap_global_days"] = int(days_global)

        # Sequencing intensity: genomes per year for this lineage
        if "collection_year" in lineage_meta.columns:
            year_counts = lineage_meta.dropna(subset=["collection_year"]).groupby("collection_year").size().to_dict()
            metrics["genomes_per_year"] = {int(k): int(v) for k, v in year_counts.items()}
            metrics["median_genomes_per_year"] = round(pd.Series(year_counts).median(), 1) if year_counts else 0

        # Geographic coverage gaps
        if matched is not None:
            reported_countries = set(_parse_pg_array(_safe_str(matched.get("countries_reported"))))
            recent_countries = set()
            if not lineage_meta.empty and "collection_date_dt" in lineage_meta.columns:
                recent = lineage_meta[lineage_meta["collection_date_dt"] >= pd.Timestamp("2023-01-01")]
                recent_countries = set(recent["collection_country"].astype(str).str.strip().unique())
            metrics["lineage_reported_countries"] = sorted(reported_countries)
            metrics["recent_sequenced_countries"] = sorted(recent_countries)
            missing = reported_countries - recent_countries - {sample_country}
            metrics["countries_without_recent_sequences"] = sorted(missing)
            if missing:
                findings.append(
                    EvidenceRecord(
                        assertion=f"{len(missing)} countries with historical lineage reports have no sequences since 2023: {', '.join(sorted(missing)[:5])}.",
                        source="lineages.csv + genome_metadata.csv",
                        confidence="low",
                        biological_relevance="surveillance coverage",
                        epidemiological_relevance="undetected circulation or importation risk",
                        finding_type="geographic_coverage_gap",
                    )
                )

        summary = (
            f"Surveillance: {len(lineage_meta)} genomes for lineage; "
            f"last global sequence {metrics.get('last_sequenced_case_global') or 'unknown'}; "
            f"{len(metrics.get('countries_without_recent_sequences', []))} countries with recent coverage gaps."
        )
        return AnalysisResult(
            title="surveillance_metrics",
            findings=findings,
            gaps=gaps,
            summary=summary,
            metrics=metrics,
        )


# ---------------------------------------------------------------------------
# 4. Lineage behaviour
# ---------------------------------------------------------------------------
class LineageBehaviorAnalyzer:
    """Characterise lineage expansion, persistence and spillover risk."""

    def __init__(
        self,
        lineages_df: pd.DataFrame,
        genome_metadata_df: pd.DataFrame,
    ):
        self.lineages = lineages_df
        self.genome_metadata = genome_metadata_df

    def analyze(self, stage9: dict) -> AnalysisResult:
        findings: list[EvidenceRecord] = []
        gaps: list[str] = []
        metrics: dict[str, Any] = {}

        lineage = _safe_str(stage9.get("lineage"))
        matched = _match_lineage(lineage, self.lineages)
        lineage_id = _safe_str(matched.get("lineage_id")) if matched is not None else ""

        if not lineage_id or self.genome_metadata.empty:
            gaps.append("Lineage or genome metadata unavailable for lineage behaviour analysis.")
            return AnalysisResult(
                title="lineage_behavior",
                findings=findings,
                gaps=gaps,
                summary="Lineage behaviour: insufficient metadata.",
                metrics=metrics,
            )

        lineage_meta = self.genome_metadata[self.genome_metadata["lineage_id"].astype(str) == lineage_id].copy()
        lineage_meta["collection_date_dt"] = pd.to_datetime(lineage_meta["collection_date"], errors="coerce")
        lineage_meta["collection_year"] = pd.to_numeric(
            lineage_meta["collection_date"].astype(str).str[:4], errors="coerce"
        )

        # Temporal span and persistence
        dates = lineage_meta["collection_date_dt"].dropna()
        if not dates.empty:
            first = dates.min()
            last = dates.max()
            metrics["first_detection_date"] = str(first.date())
            metrics["last_detection_date"] = str(last.date())
            metrics["detection_span_years"] = round((last - first).days / 365.25, 2)

            year_counts = lineage_meta.dropna(subset=["collection_year"]).groupby("collection_year").size()
            metrics["years_with_genomes"] = sorted([int(y) for y in year_counts.index])
            metrics["genome_count"] = int(len(lineage_meta))

            # Persistence index
            active_years = set(year_counts.index.astype(int))
            if len(active_years) >= 5:
                persistence = "continuous"
            elif max(active_years) - min(active_years) >= 5:
                persistence = "intermittent"
            else:
                persistence = "limited"
            metrics["persistence_index"] = persistence

            if persistence == "continuous":
                findings.append(
                    EvidenceRecord(
                        assertion=f"Lineage has been detected across {len(active_years)} years ({min(active_years)}-{max(active_years)}), indicating continuous circulation.",
                        source="genome_metadata.csv",
                        confidence="medium",
                        biological_relevance="lineage persistence",
                        epidemiological_relevance="ongoing transmission potential",
                        finding_type="continued_circulation",
                    )
                )
            elif persistence == "intermittent":
                findings.append(
                    EvidenceRecord(
                        assertion=f"Lineage has intermittent detections over {metrics['detection_span_years']} years with reporting gaps, consistent with re-emergence.",
                        source="genome_metadata.csv",
                        confidence="medium",
                        biological_relevance="re-emergence",
                        epidemiological_relevance="sporadic outbreaks or cryptic transmission",
                        finding_type="re-emergence_after_a_reporting_gap",
                    )
                )

        # Geographic expansion rate
        country_years = lineage_meta.dropna(subset=["collection_year"]).groupby("collection_country")["collection_year"].min().to_dict()
        if country_years:
            sorted_countries = sorted(country_years.items(), key=lambda x: x[1])
            metrics["first_detection_by_country"] = {c: int(y) for c, y in sorted_countries}
            if len(sorted_countries) >= 2 and metrics.get("detection_span_years"):
                metrics["geographic_expansion_rate_countries_per_year"] = round(
                    len(sorted_countries) / max(metrics["detection_span_years"], 0.5), 2
                )

        # Spillover risk
        reservoir = _safe_str(matched.get("reservoir")) if matched is not None else ""
        primary_host = _safe_str(matched.get("primary_host")) if matched is not None else ""
        host = _safe_str(stage9.get("metadata", {}).get("host")).lower()

        spillover_score = 0
        if reservoir and reservoir.lower() not in ("human", "unknown"):
            spillover_score += 2
        if primary_host and primary_host.lower() not in ("human", "unknown"):
            spillover_score += 1
        if host and host not in ("human", ""):
            spillover_score += 1

        metrics["spillover_risk_score"] = min(spillover_score, 5)
        metrics["known_reservoir"] = reservoir or None
        metrics["primary_host"] = primary_host or None

        if spillover_score >= 2:
            findings.append(
                EvidenceRecord(
                    assertion=f"Known animal reservoir ({reservoir}) and primary host ({primary_host}) suggest spillover potential.",
                    source="lineages.csv",
                    confidence="low",
                    biological_relevance="host range",
                    epidemiological_relevance="zoonotic introduction risk",
                    finding_type="spillover_risk",
                )
            )

        summary = (
            f"Lineage behaviour: {metrics.get('genome_count', 0)} genomes, "
            f"persistence {metrics.get('persistence_index', 'unknown')}, "
            f"{len(country_years)} countries reported."
        )
        return AnalysisResult(
            title="lineage_behavior",
            findings=findings,
            gaps=gaps,
            summary=summary,
            metrics=metrics,
        )


# ---------------------------------------------------------------------------
# 5. Integrated epi+genomic inference
# ---------------------------------------------------------------------------
class IntegratedInferenceAnalyzer:
    """Combine epi, tree and genomic evidence for actionable inferences."""

    def __init__(self, epi_output: dict):
        self.epi = epi_output or {}

    def analyze(
        self,
        stage9: dict,
        tree_input: Any,
        all_results: dict[str, AnalysisResult],
    ) -> AnalysisResult:
        findings: list[EvidenceRecord] = []
        gaps: list[str] = []
        metrics: dict[str, Any] = {}

        sample_country = _safe_str(stage9.get("metadata", {}).get("country"))
        sample_date = _parse_date(stage9.get("metadata", {}).get("collection_date"))

        # Reporting delay from molecular clock placement
        clock_result = all_results.get("molecular_clock_analysis")
        placement_date = None
        if clock_result:
            for f in clock_result.findings:
                if f.finding_type == "sample_placement_timing":
                    # Extract date from assertion text
                    import re

                    m = re.search(r"(\d{4}-\d{2})", f.assertion)
                    if m:
                        placement_date = m.group(1)
                    break

        if placement_date and sample_date:
            delay_days = _days_between(sample_date, placement_date)
            if delay_days is not None:
                metrics["reporting_delay_days"] = int(delay_days)
                if delay_days > 30:
                    findings.append(
                        EvidenceRecord(
                            assertion=f"Sample collection is {delay_days} days after the molecular-clock placement date, suggesting a reporting or transmission-chain delay.",
                            source="molecular clock + sample metadata",
                            confidence="low",
                            biological_relevance="evolutionary/epidemiological lag",
                            epidemiological_relevance="undetected intermediate cases",
                            finding_type="reporting_delay",
                        )
                    )

        # Undetected transmission chain length from time-scaled tree placement
        if tree_input and tree_input.sample_tip:
            time_scaled = (stage9.get("stage6_phylogenetic_tree") or {}).get("time_scaled_tree") or {}
            placement_date = time_scaled.get("sample_placement_date")
            root_age = time_scaled.get("root_age")
            reference_date = placement_date or root_age
            if reference_date and sample_date:
                delta_days = _days_between(sample_date, reference_date)
                if delta_days is not None and delta_days > 0:
                    # Ebola generation time ~14 days
                    generations = delta_days / 14
                    metrics["estimated_undetected_generations"] = round(generations, 1)
                    if generations > 2:
                        findings.append(
                            EvidenceRecord(
                                assertion=f"Sample collection is {int(delta_days)} days after the time-scaled placement ({reference_date}), implying ~{int(generations)} generations of possible undetected transmission.",
                                source="time-scaled phylogeny + sample metadata",
                                confidence="low",
                                biological_relevance="cryptic transmission",
                                epidemiological_relevance="missed cases in the chain",
                                finding_type="undetected_transmission_chain",
                            )
                        )
                elif delta_days is not None and delta_days <= 0:
                    metrics["estimated_undetected_generations"] = 0

        # Importation vs local transmission probability
        relatedness_result = all_results.get("genetic_relatedness_analysis")
        if relatedness_result:
            close_match = any(
                f.finding_type == "close_genetic_match"
                for f in relatedness_result.findings
            )
            temporal_gap_days = relatedness_result.metrics.get("min_temporal_gap_days")
            # A genetically close match within ~6 months supports local/recent transmission;
            # otherwise an introduction/persistence event is more probable.
            if close_match and temporal_gap_days is not None and temporal_gap_days <= 180:
                metrics["importation_probability"] = 0.2
            else:
                metrics["importation_probability"] = 0.7

        # International spread risk
        spread_score = 0
        geo_result = all_results.get("phylogeographic_analysis")
        if geo_result:
            for f in geo_result.findings:
                if f.finding_type == "dissemination_pathways":
                    spread_score += 2
                if f.finding_type == "geographic_mismatch":
                    spread_score += 1
        threat_result = all_results.get("evidence_weighted_threat")
        if threat_result and threat_result.metrics.get("international_spread"):
            spread_score += threat_result.metrics["international_spread"]
        metrics["international_spread_risk_score"] = min(int(spread_score), 5)

        if metrics.get("international_spread_risk_score", 0) >= 3:
            findings.append(
                EvidenceRecord(
                    assertion="Multiple dissemination signals suggest elevated risk of international spread.",
                    source="phylogeography + threat integration",
                    confidence="medium",
                    biological_relevance="geographic dissemination",
                    epidemiological_relevance="exportation/importation potential",
                    finding_type="international_spread_risk",
                )
            )

        if not findings:
            findings.append(
                EvidenceRecord(
                    assertion="Integrated epi+genomic inference did not identify strong reporting delays or undetected chains.",
                    source="epi output + tree + clock",
                    confidence="low",
                    biological_relevance="na",
                    epidemiological_relevance="limited integrated signal",
                    finding_type="no_major_integrated_signal",
                )
            )

        summary = (
            f"Integrated inference: reporting delay {metrics.get('reporting_delay_days') or 'unknown'} days, "
            f"undetected generations {metrics.get('estimated_undetected_generations') or 'unknown'}, "
            f"international spread risk {metrics.get('international_spread_risk_score') or 0}/5."
        )
        return AnalysisResult(
            title="integrated_inference",
            findings=findings,
            gaps=gaps,
            summary=summary,
            metrics=metrics,
        )


# ---------------------------------------------------------------------------
# 6. Evidence consistency
# ---------------------------------------------------------------------------
class EvidenceConsistencyAnalyzer:
    """Compare genomic, epidemiological and literature evidence for agreement or conflict."""

    def __init__(
        self,
        epi_output: dict,
        lineages_df: pd.DataFrame,
        genome_metadata_df: pd.DataFrame,
    ):
        self.epi = epi_output or {}
        self.lineages = lineages_df
        self.genome_metadata = genome_metadata_df

    def analyze(self, stage9: dict, all_results: dict[str, AnalysisResult]) -> AnalysisResult:
        findings: list[EvidenceRecord] = []
        gaps: list[str] = []
        metrics: dict[str, Any] = {"consistency_checks": []}

        sample_country = _safe_str(stage9.get("metadata", {}).get("country"))
        travel_locs = _parse_pg_array(_safe_str(stage9.get("metadata", {}).get("travel_locations")))
        sample_date = _parse_date(stage9.get("metadata", {}).get("collection_date"))

        # 1. Phylogeographic origin vs travel history
        geo_result = all_results.get("phylogeographic_analysis")
        inferred_origin = ""
        if geo_result and geo_result.metrics:
            inferred_origin = _safe_str(geo_result.metrics.get("inferred_origin"))

        def _normalise_country(name: str) -> str:
            n = name.lower().strip()
            # Common aliases
            if n in {"drc", "democratic republic of the congo", "dr congo", "congo-kinshasa", "congo (drc)"}:
                return "democratic republic of the congo"
            if n in {"roc", "congo", "congo-brazzaville", "republic of the congo"}:
                return "republic of the congo"
            if n in {"usa", "united states", "united states of america", "us"}:
                return "united states"
            if n in {"uk", "united kingdom", "great britain"}:
                return "united kingdom"
            return n

        if inferred_origin and travel_locs:
            inferred_norm = _normalise_country(inferred_origin)
            loc_norms = [_normalise_country(loc) for loc in travel_locs]
            match = inferred_norm in loc_norms or any(inferred_norm in ln or ln in inferred_norm for ln in loc_norms)
            if match:
                msg = f"Inferred phylogeographic origin ({inferred_origin}) is consistent with reported travel history ({', '.join(travel_locs)})."
                findings.append(
                    EvidenceRecord(
                        assertion=msg,
                        source="phylogeographic analysis + sample metadata",
                        confidence="medium",
                        biological_relevance="geographic origin",
                        epidemiological_relevance="travel exposure supports origin inference",
                        finding_type="consistent_evidence",
                    )
                )
                metrics["consistency_checks"].append({"check": "origin_vs_travel", "consistent": True, "detail": msg})
            else:
                msg = f"Inferred phylogeographic origin ({inferred_origin}) does not match reported travel locations ({', '.join(travel_locs)}); other exposure sources may exist."
                findings.append(
                    EvidenceRecord(
                        assertion=msg,
                        source="phylogeographic analysis + sample metadata",
                        confidence="low",
                        biological_relevance="geographic origin",
                        epidemiological_relevance="origin and travel history diverge",
                        finding_type="conflicting_evidence",
                    )
                )
                metrics["consistency_checks"].append({"check": "origin_vs_travel", "consistent": False, "detail": msg})
        elif inferred_origin:
            gaps.append("No travel history available; cannot cross-check phylogeographic origin with epidemiological exposure.")
        else:
            gaps.append("No phylogeographic origin inferred; consistency with travel history cannot be assessed.")

        # 2. Lineage temporal scenario vs molecular clock / genetic relatedness
        relatedness_result = all_results.get("genetic_relatedness_analysis")
        molecular_clock = all_results.get("molecular_clock_analysis")
        lineage_behavior = all_results.get("lineage_behavior")

        has_clock = molecular_clock and any(f.finding_type == "molecular_clock_timing" for f in molecular_clock.findings)
        has_recent_relatedness = (
            relatedness_result
            and relatedness_result.metrics.get("min_temporal_gap_days") is not None
            and relatedness_result.metrics["min_temporal_gap_days"] <= 180
        )
        persistence = ""
        if lineage_behavior and lineage_behavior.metrics:
            persistence = _safe_str(lineage_behavior.metrics.get("persistence_index"))

        if persistence == "continuous" and not has_recent_relatedness and not has_clock:
            msg = "Lineage shows continuous historical circulation, but no recent close genetic match or molecular-clock signal is available; evidence for ongoing local transmission is incomplete."
            findings.append(
                EvidenceRecord(
                    assertion=msg,
                    source="lineage behavior + genetic relatedness + molecular clock",
                    confidence="low",
                    biological_relevance="lineage persistence",
                    epidemiological_relevance="historical continuity without confirming recent linkage",
                    finding_type="insufficient_evidence",
                )
            )
            metrics["consistency_checks"].append({"check": "persistence_vs_recent_relatedness", "consistent": None, "detail": msg})
        elif persistence and has_recent_relatedness:
            msg = "Continuous lineage persistence and a close recent genomic match are consistent with ongoing or recently linked transmission."
            findings.append(
                EvidenceRecord(
                    assertion=msg,
                    source="lineage behavior + genetic relatedness",
                    confidence="medium",
                    biological_relevance="lineage persistence",
                    epidemiological_relevance="genomic and temporal signals align",
                    finding_type="consistent_evidence",
                )
            )
            metrics["consistency_checks"].append({"check": "persistence_vs_recent_relatedness", "consistent": True, "detail": msg})

        # 3. Phenotype evidence vs mutation catalogue
        genomic_signal = all_results.get("genomic_signal")
        if genomic_signal and genomic_signal.metrics:
            known = genomic_signal.metrics.get("known_variants_in_catalogue", 0)
            novel = genomic_signal.metrics.get("novel_variants_not_in_catalogue", 0)
            if known and novel:
                msg = f"Sample contains both catalogue-known ({known}) and catalogue-novel ({novel}) variants; phenotype associations are limited to known mutations."
                findings.append(
                    EvidenceRecord(
                        assertion=msg,
                        source="genomic_signal analyzer",
                        confidence="medium",
                        biological_relevance="variant characterisation",
                        epidemiological_relevance="functional interpretations may be incomplete",
                        finding_type="insufficient_evidence",
                    )
                )
                metrics["consistency_checks"].append({"check": "known_vs_novel_variants", "consistent": None, "detail": msg})
            elif novel and not known:
                msg = "All detected variants are catalogue-novel; curated phenotype evidence is unavailable."
                findings.append(
                    EvidenceRecord(
                        assertion=msg,
                        source="genomic_signal analyzer",
                        confidence="medium",
                        biological_relevance="variant characterisation",
                        epidemiological_relevance="phenotypic impact unknown",
                        finding_type="insufficient_evidence",
                    )
                )
                metrics["consistency_checks"].append({"check": "known_vs_novel_variants", "consistent": None, "detail": msg})

        if not findings:
            findings.append(
                EvidenceRecord(
                    assertion="No clear consistencies or conflicts were identified across the available evidence streams.",
                    source="evidence consistency analyzer",
                    confidence="low",
                    biological_relevance="evidence integration",
                    epidemiological_relevance="limited cross-evidence comparison",
                    finding_type="insufficient_evidence",
                )
            )

        summary = (
            f"Evidence consistency: {sum(1 for c in metrics['consistency_checks'] if c.get('consistent') is True)} consistent, "
            f"{sum(1 for c in metrics['consistency_checks'] if c.get('consistent') is False)} conflicting, "
            f"{sum(1 for c in metrics['consistency_checks'] if c.get('consistent') is None)} indeterminate checks."
        )
        return AnalysisResult(
            title="evidence_consistency",
            findings=findings,
            gaps=gaps,
            summary=summary,
            metrics=metrics,
        )


# ---------------------------------------------------------------------------
# 7. Knowledge gap analyzer
# ---------------------------------------------------------------------------
class KnowledgeGapAnalyzer:
    """Aggregate explicit gaps from all analyzers and add cross-cutting data sufficiency checks."""

    def __init__(self, epi_output: dict):
        self.epi = epi_output or {}

    def analyze(
        self,
        stage9: dict,
        all_results: dict[str, AnalysisResult],
    ) -> AnalysisResult:
        findings: list[EvidenceRecord] = []
        gaps: list[str] = []

        # Collect gaps reported by each analyzer
        for result in all_results.values():
            for g in result.gaps:
                if g and g not in gaps:
                    gaps.append(g)

        # Cross-cutting metadata sufficiency checks
        sample_meta = stage9.get("metadata", {})
        if not _safe_str(sample_meta.get("collection_date")):
            gaps.append("Sample collection date is missing; temporal analyses are limited.")
        if not _safe_str(sample_meta.get("country")):
            gaps.append("Sample country is missing; geographic analyses are limited.")
        if not _safe_str(sample_meta.get("suspected_exposure")):
            gaps.append("Suspected exposure route is missing; transmission hypotheses cannot be refined.")
        if not _safe_str(sample_meta.get("travel_history")):
            gaps.append("Travel history is missing; importation risk is harder to assess.")

        # Epi output sufficiency
        if not self.epi.get("outbreaks"):
            gaps.append("No species outbreak records in epi output; epidemic dynamics are inferred from curated genomes only.")
        if not self.epi.get("transmission"):
            gaps.append("No transmission parameters (R0, serial interval) in epi output; transmission potential cannot be quantified.")

        # Phylogenetic sufficiency
        tree = stage9.get("stage6_phylogenetic_tree", {})
        if not tree.get("tree_file") and not stage9.get("phylogenetic_placement"):
            gaps.append("No phylogenetic tree or placement provided; molecular-clock and relatedness analyses rely on metadata.")
        elif not tree.get("time_scaled_tree"):
            gaps.append("No time-scaled phylogeny available; molecular-clock dating is not possible.")

        # Literature/phenotype sufficiency
        if not stage9.get("closest_reference"):
            gaps.append("No closest reference genome information; assembly context is incomplete.")

        # Summarise
        if gaps:
            findings.append(
                EvidenceRecord(
                    assertion=f"{len(gaps)} knowledge or data gaps were identified across the genomic, epidemiological and literature evidence streams.",
                    source="knowledge gap analyzer",
                    confidence="high",
                    biological_relevance="data sufficiency",
                    epidemiological_relevance="limits confidence in integrated conclusions",
                    finding_type="knowledge_gap_summary",
                )
            )

        summary = f"Knowledge gaps: {len(gaps)} cross-cutting gaps identified."
        return AnalysisResult(
            title="knowledge_gaps",
            findings=findings,
            gaps=gaps,
            summary=summary,
            metrics={"gap_count": len(gaps), "gaps": gaps},
        )


# ---------------------------------------------------------------------------
# 8. Comparative outbreak analysis
# ---------------------------------------------------------------------------
class ComparativeOutbreakAnalyzer:
    """Compare the current detection to historical outbreaks to find the most similar event."""

    def __init__(
        self,
        epi_output: dict,
        lineages_df: pd.DataFrame,
        genome_metadata_df: pd.DataFrame,
    ):
        self.epi = epi_output or {}
        self.lineages = lineages_df
        self.genome_metadata = genome_metadata_df

    def analyze(self, stage9: dict) -> AnalysisResult:
        findings: list[EvidenceRecord] = []
        gaps: list[str] = []
        metrics: dict[str, Any] = {}

        sample_country = _safe_str(stage9.get("metadata", {}).get("country"))
        sample_date = _parse_date(stage9.get("metadata", {}).get("collection_date"))
        sample_lineage = _safe_str(stage9.get("lineage"))
        sample_species = _safe_str(stage9.get("species"))

        outbreaks = self.epi.get("outbreaks", [])
        records = []
        for o in outbreaks:
            cases = _to_float(o.get("cases"))
            deaths = _to_float(o.get("deaths"))
            cfr_pct = _to_float(o.get("cfr"))
            cfr = cfr_pct / 100 if cfr_pct is not None else None
            records.append({
                "date": _safe_str(o.get("start_date")),
                "country": _safe_str(o.get("country")),
                "cases": cases,
                "deaths": deaths,
                "cfr": cfr,
            })

        df = pd.DataFrame(records)

        # Score historical outbreaks by lineage, country, and temporal proximity
        scored = []
        for _, row in df.iterrows():
            score = 0
            reasons = []
            if sample_country and row["country"].lower() == sample_country.lower():
                score += 2
                reasons.append("same country")
            if sample_date and row["date"]:
                try:
                    d = pd.to_datetime(row["date"], errors="coerce")
                    if pd.notna(d):
                        years = abs(pd.to_datetime(sample_date).year - d.year)
                        if years <= 5:
                            score += 2
                            reasons.append("within 5 years")
                        elif years <= 10:
                            score += 1
                            reasons.append("within 10 years")
                except Exception:
                    pass
            if row["cases"] is not None and row["cases"] > 0:
                score += 0.5
            if score > 0:
                scored.append({
                    "date": row["date"],
                    "country": row["country"],
                    "cases": row["cases"],
                    "deaths": row["deaths"],
                    "cfr": row["cfr"],
                    "similarity_score": score,
                    "reasons": reasons,
                })

        scored.sort(key=lambda x: x["similarity_score"], reverse=True)
        metrics["most_similar_outbreaks"] = scored[:5]

        if scored:
            top = scored[0]
            findings.append(
                EvidenceRecord(
                    assertion=(
                        f"Most similar historical outbreak: {top['country']} "
                        f"({top['date'] or 'unknown date'}), "
                        f"similarity score {top['similarity_score']} "
                        f"({', '.join(top['reasons'])})."
                    ),
                    source="outbreaks",
                    confidence="low" if len(scored) < 5 else "medium",
                    biological_relevance="historical comparison",
                    epidemiological_relevance="context for expected scale and dynamics",
                    finding_type="comparative_outbreak_match",
                )
            )
        else:
            gaps.append("No comparable historical outbreaks with country/date overlap were found in the epi output.")
            findings.append(
                EvidenceRecord(
                    assertion="No sufficiently similar historical outbreak could be identified from available epidemiological records.",
                    source="outbreaks",
                    confidence="low",
                    biological_relevance="historical comparison",
                    epidemiological_relevance="insufficient outbreak metadata for comparison",
                    finding_type="no_comparable_outbreak",
                )
            )

        summary = (
            f"Comparative outbreak analysis: {len(scored)} historical outbreak(s) with non-zero similarity; "
            f"top match in {metrics['most_similar_outbreaks'][0]['country'] if metrics['most_similar_outbreaks'] else 'none'}."
        )
        return AnalysisResult(
            title="comparative_outbreak_analysis",
            findings=findings,
            gaps=gaps,
            summary=summary,
            metrics=metrics,
        )
