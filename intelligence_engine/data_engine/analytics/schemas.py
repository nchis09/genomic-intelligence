"""schemas.py — Normalized, entity-based epidemiological master schema.

These Pydantic v2 models define the single validated master object
(``NormalizedEpiObject``) that ``EpiQueryEngine.query`` produces. It replaces
the old question-based (Q1-Q14) output: instead of one array per question,
the LLM (or the deterministic fallback) populates fixed entity sections --
pathogen profile, molecular epidemiology, outbreaks, transmission,
demographics, clinical features, interventions, diagnostics, therapeutics,
vaccines, surveillance, genomic links, knowledge assertions, and references.

Only objective, structured epidemiological evidence lives here. Narrative
reasoning (temporal assessments, novelty judgments, risk interpretation) is
NOT extracted at this stage -- it is generated later by the Genomic
Intelligence Engine, which integrates this structured evidence with genomic
analyses.

``epi_output.json`` is now this object's ``model_dump()``. Downstream
consumers must read entity keys (``outbreaks``, ``transmission``, etc.)
rather than the old ``Q1_species_outbreaks`` / ``Q5_transmission_params``
question keys.
"""
from __future__ import annotations

import hashlib
import re
from typing import Optional

from pydantic import BaseModel, Field, field_validator, model_validator


def _coerce_to_str_field(v):
    """Coerce a value bound for a str/Optional[str] field into a plain string.

    LLMs frequently return a single-element list for fields we asked to be a
    scalar string (e.g. ``source_url: ["https://..."]``) or a bare int/float
    for date-like fields (e.g. ``start_date: 2004``). Rather than dropping
    the whole row for a superficial type mismatch, normalise it here.
    """
    if isinstance(v, list):
        if not v:
            return None
        return ", ".join(str(x) for x in v)
    if isinstance(v, (int, float)):
        return str(v)
    return v


class LenientModel(BaseModel):
    """Base model that coerces list/int/float inputs into strings for any
    field annotated as ``str`` or ``Optional[str]`` before validation.

    This absorbs common LLM formatting quirks (wrapping a scalar in a
    single-element list, returning a year as an int) without silently
    dropping the whole row.
    """

    @model_validator(mode="before")
    @classmethod
    def _coerce_str_fields(cls, data):
        if not isinstance(data, dict):
            return data
        out = dict(data)
        for name, info in cls.model_fields.items():
            annotation = info.annotation
            args = getattr(annotation, "__args__", None)
            is_str_field = annotation is str or (args and str in args)
            if is_str_field and name in out:
                out[name] = _coerce_to_str_field(out[name])
        return out


def _make_outbreak_id(pathogen: Optional[str], country: Optional[str],
                       start_date: Optional[str], cases: Optional[int] = None,
                       deaths: Optional[int] = None) -> str:
    """Deterministic ID so the same outbreak always hashes to the same id
    regardless of which provider reported it.

    Only pathogen, country, and start year are used (not cases/deaths) so
    that multiple providers reporting the same event -- with slightly
    different case/death counts -- collapse to the same outbreak_id. Callers
    are expected to reconcile differing numeric fields during merge.
    """
    year = ""
    if start_date:
        m = re.search(r"\d{4}", str(start_date))
        if m:
            year = m.group(0)
    key = f"{(pathogen or '').lower().strip()}|{(country or '').lower().strip()}|{year}"
    return hashlib.sha1(key.encode("utf-8")).hexdigest()[:16]


def _coerce_int(v):
    if v is None or v == "":
        return None
    if isinstance(v, (int, float)):
        return int(v)
    if isinstance(v, str):
        m = re.search(r"[\d,]+", v)
        if m:
            try:
                return int(m.group(0).replace(",", ""))
            except ValueError:
                return None
    return None


def _coerce_float(v):
    if v is None or v == "":
        return None
    if isinstance(v, (int, float)):
        return float(v)
    if isinstance(v, str):
        m = re.search(r"[\d.]+", v)
        if m:
            try:
                return float(m.group(0))
            except ValueError:
                return None
    return None


class OutbreakRecord(LenientModel):
    """One row per outbreak. This is the primary analytical fact table."""

    outbreak_id: str = Field(default="", description="Deterministic hash ID for deduplication/joins")
    pathogen: Optional[str] = None
    lineage: Optional[str] = None
    country: Optional[str] = None
    admin_region: Optional[str] = None
    latitude: Optional[float] = None
    longitude: Optional[float] = None
    start_date: Optional[str] = None
    end_date: Optional[str] = None
    collection_date: Optional[str] = None
    cases: Optional[int] = None
    deaths: Optional[int] = None
    cfr: Optional[float] = Field(None, description="Case fatality rate, percentage (0-100)")
    incidence: Optional[float] = None
    prevalence: Optional[float] = None
    source_url: Optional[str] = None
    reporting_agency: Optional[str] = None
    confidence_score: Optional[float] = Field(None, ge=0, le=1)
    evidence_level: Optional[str] = Field(
        None, description="'structured_api' | 'llm_extracted' | 'text_regex'"
    )
    retrieval_date: Optional[str] = None

    @field_validator("cases", "deaths", mode="before")
    @classmethod
    def _v_int(cls, v):
        return _coerce_int(v)

    @field_validator("cfr", "incidence", "prevalence", "latitude", "longitude", mode="before")
    @classmethod
    def _v_float(cls, v):
        return _coerce_float(v)

    @field_validator("cfr", mode="after")
    @classmethod
    def _v_cfr_range(cls, v):
        if v is not None and (v < 0 or v > 100):
            return None
        return v

    def model_post_init(self, __context) -> None:
        if not self.outbreak_id:
            self.outbreak_id = _make_outbreak_id(
                self.pathogen, self.country, self.start_date, self.cases, self.deaths
            )
        if self.cfr is None and self.cases and self.deaths is not None:
            try:
                self.cfr = round((self.deaths / self.cases) * 100, 2) if self.cases else None
            except ZeroDivisionError:
                self.cfr = None


class EpiMetadata(LenientModel):
    """Run-level context: what was queried, when, and against what sources."""

    query_timestamp: Optional[str] = None
    species: Optional[str] = None
    lineage: Optional[str] = None
    country: Optional[str] = None
    collection_date: Optional[str] = None
    total_reports_fetched: Optional[int] = None
    api_records_fetched: Optional[int] = None
    provider: Optional[str] = None
    model: Optional[str] = None
    structuring_method: Optional[str] = Field(
        None, description="'llm_entity_extraction' | 'deterministic_extraction'"
    )

    @field_validator("total_reports_fetched", "api_records_fetched", mode="before")
    @classmethod
    def _v_int(cls, v):
        return _coerce_int(v)


class PathogenProfile(LenientModel):
    """Species-level profile, independent of any single outbreak."""

    species: Optional[str] = None
    pathogen_family: Optional[str] = None
    pathogen_genus: Optional[str] = None
    reservoir: Optional[str] = None
    host: Optional[str] = None
    first_documented_year: Optional[int] = None
    first_documented_location: Optional[str] = None
    pathogenic_in_humans: Optional[bool] = None
    confirmed_outbreaks_count: Optional[int] = None
    source_url: Optional[str] = None

    @field_validator("first_documented_year", "confirmed_outbreaks_count", mode="before")
    @classmethod
    def _v_int(cls, v):
        return _coerce_int(v)


class MolecularEpidemiology(LenientModel):
    """Strain/genotype/clade-level epidemiological burden (which circulating
    strain/lineage caused how many cases/deaths, where, and when).

    This is NOT where individual mutations live -- specific amino-acid
    changes belong in ``genomic_links``. This section answers "how did
    strain/genotype/clade X behave epidemiologically", not "what mutations
    does strain X carry".
    """

    lineage: Optional[str] = None
    strain: Optional[str] = None
    genotype: Optional[str] = None
    clade: Optional[str] = None
    country: Optional[str] = None
    year: Optional[str] = None
    cases: Optional[int] = None
    deaths: Optional[int] = None
    cfr: Optional[float] = None
    key_transmission_features: Optional[str] = None
    source_url: Optional[str] = None

    @field_validator("cases", "deaths", mode="before")
    @classmethod
    def _v_int(cls, v):
        return _coerce_int(v)

    @field_validator("cfr", mode="before")
    @classmethod
    def _v_float(cls, v):
        return _coerce_float(v)


class TransmissionParams(LenientModel):
    """Pathogen-level transmission parameters: R0, incubation, serial interval, route."""

    r0_low: Optional[float] = None
    r0_high: Optional[float] = None
    incubation_period_days_low: Optional[float] = None
    incubation_period_days_high: Optional[float] = None
    serial_interval_days: Optional[float] = None
    transmission_route: Optional[str] = None
    host: Optional[str] = None
    reservoir: Optional[str] = None
    vector: Optional[str] = Field(None, description="Arthropod/animal vector, if vector-borne (One Health)")
    source_url: Optional[str] = None

    @field_validator(
        "r0_low", "r0_high", "incubation_period_days_low",
        "incubation_period_days_high", "serial_interval_days", mode="before"
    )
    @classmethod
    def _v_float(cls, v):
        return _coerce_float(v)


class DemographicSummary(LenientModel):
    """Who was affected and how -- age/sex/occupation breakdowns plus the
    routinely-reported variables needed for outbreak investigation:
    population group, risk group, exposure history, infection setting, and
    (for zoonotic pathogens) the affected host species."""

    age_group: Optional[str] = None
    sex: Optional[str] = None
    occupation: Optional[str] = Field(
        None, description="e.g. healthcare worker, farmer, veterinarian"
    )
    population_affected: Optional[str] = Field(
        None, description="e.g. children, adults, elderly, pregnant women, immunocompromised"
    )
    risk_group: Optional[str] = Field(
        None, description="e.g. healthcare workers, household contacts, travelers, refugees"
    )
    exposure_history: Optional[str] = Field(
        None, description="e.g. animal contact, healthcare exposure, funeral attendance, "
                           "travel, vector exposure"
    )
    setting: Optional[str] = Field(
        None, description="Setting of infection: community, healthcare facility, "
                           "household, workplace, school, refugee camp, etc."
    )
    host_species: Optional[str] = Field(
        None, description="Affected host species, for zoonotic pathogens (e.g. cattle, "
                           "primates, bats, humans)"
    )
    case_count: Optional[int] = None
    risk_factor: Optional[str] = None
    notes: Optional[str] = None
    source_url: Optional[str] = None

    @field_validator("case_count", mode="before")
    @classmethod
    def _v_int(cls, v):
        return _coerce_int(v)


class ClinicalFeature(LenientModel):
    """Clinical severity, symptoms, and atypical presentations."""

    feature: Optional[str] = None
    lineage: Optional[str] = None
    severity: Optional[str] = None
    frequency: Optional[str] = None
    notes: Optional[str] = None
    source_url: Optional[str] = None


class InterventionRecord(LenientModel):
    """Public-health measures (contact tracing, isolation, PPE, safe burial, etc.)."""

    intervention_type: Optional[str] = Field(
        None, description="'public_health_measure' | 'clinical_management' | 'ipc'"
    )
    name: Optional[str] = None
    status: Optional[str] = None
    effectiveness: Optional[str] = None
    notes: Optional[str] = None
    source_url: Optional[str] = None


class DiagnosticMethod(LenientModel):
    """Diagnostic and laboratory confirmation methods."""

    method: Optional[str] = None
    method_type: Optional[str] = Field(None, alias="type")
    target: Optional[str] = None
    notes: Optional[str] = None
    source_url: Optional[str] = None

    model_config = {"populate_by_name": True}


class TherapeuticProduct(LenientModel):
    """Therapeutic products and their status/effectiveness."""

    product: Optional[str] = None
    status: Optional[str] = None
    effectiveness: Optional[str] = None
    notes: Optional[str] = None
    source_url: Optional[str] = None


class VaccineProduct(LenientModel):
    """Vaccine products and their status/effectiveness."""

    product: Optional[str] = None
    status: Optional[str] = None
    effectiveness: Optional[str] = None
    notes: Optional[str] = None
    source_url: Optional[str] = None


class SurveillanceSummary(LenientModel):
    """Geographic range and surveillance context, one row per country."""

    country: Optional[str] = None
    region: Optional[str] = None
    first_documented: Optional[int] = None
    reservoir: Optional[str] = None
    seroprevalence_pct: Optional[float] = Field(
        None, description="Seroepidemiology: % of surveyed population with antibodies, where available"
    )
    seroprevalence_population: Optional[str] = Field(
        None, description="Population/group sampled for the seroprevalence study"
    )
    source_url: Optional[str] = None

    @field_validator("first_documented", mode="before")
    @classmethod
    def _v_int(cls, v):
        return _coerce_int(v)

    @field_validator("seroprevalence_pct", mode="before")
    @classmethod
    def _v_float(cls, v):
        return _coerce_float(v)


class GenomicLink(LenientModel):
    """Genomic/phylogenetic context tied to an outbreak."""

    outbreak_id: Optional[str] = None
    lineage: Optional[str] = None
    clade: Optional[str] = None
    mutations: Optional[str] = None
    genomic_accession: Optional[str] = None
    source_url: Optional[str] = None


class KnowledgeAssertion(LenientModel):
    """A single grounded claim extracted by the LLM that doesn't fit a strict
    tabular field, with provenance and a confidence score for downstream
    report generation and evidence-consistency checks."""

    claim: str
    source_url: Optional[str] = None
    reporting_agency: Optional[str] = None
    confidence: Optional[float] = Field(None, ge=0, le=1)
    evidence_level: Optional[str] = Field(
        None, description="'structured_api' | 'llm_extracted' | 'text_regex'"
    )

    @field_validator("confidence", mode="before")
    @classmethod
    def _v_float(cls, v):
        return _coerce_float(v)


class Reference(LenientModel):
    """Consolidated source reference."""

    source_url: str
    title: Optional[str] = None
    reporting_agency: Optional[str] = None
    retrieval_date: Optional[str] = None
    credibility: Optional[int] = None


class NormalizedEpiObject(BaseModel):
    """The single validated master object produced by EpiQueryEngine.query().

    This IS epi_output.json. Every downstream consumer (analytical Parquet
    export, Genomic Intelligence Engine analyzers, report generation) reads
    from this structure.
    """

    metadata: EpiMetadata = Field(default_factory=EpiMetadata)
    pathogen_profile: Optional[PathogenProfile] = None
    molecular_epidemiology: list[MolecularEpidemiology] = Field(default_factory=list)
    outbreaks: list[OutbreakRecord] = Field(default_factory=list)
    transmission: Optional[TransmissionParams] = None
    demographics: list[DemographicSummary] = Field(default_factory=list)
    clinical: list[ClinicalFeature] = Field(default_factory=list)
    interventions: list[InterventionRecord] = Field(default_factory=list)
    diagnostics: list[DiagnosticMethod] = Field(default_factory=list)
    therapeutics: list[TherapeuticProduct] = Field(default_factory=list)
    vaccines: list[VaccineProduct] = Field(default_factory=list)
    surveillance: list[SurveillanceSummary] = Field(default_factory=list)
    genomic_links: list[GenomicLink] = Field(default_factory=list)
    knowledge_assertions: list[KnowledgeAssertion] = Field(default_factory=list)
    references: list[Reference] = Field(default_factory=list)
