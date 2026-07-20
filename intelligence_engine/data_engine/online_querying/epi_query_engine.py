"""
epi_query_engine.py — Epidemiological context query engine.

Reads the bioinformatics output JSON (from the pipeline), fetches real
epidemiological data from trusted sources, and asks the LLM (or a
deterministic fallback) to populate a single normalized, entity-based JSON
object (``NormalizedEpiObject``) -- pathogen profile, molecular epidemiology,
outbreaks, transmission, demographics, clinical features, interventions,
diagnostics, therapeutics, vaccines, surveillance, genomic links, knowledge
assertions, and references. This object IS ``epi_output.json`` and is the
single source of truth for analytical Parquet export, the Genomic
Intelligence Engine analyzers, and report generation.

Usage:
    from intelligence_engine.data_engine.online_querying.epi_query_engine import EpiQueryEngine

    engine = EpiQueryEngine()
    epi_context = engine.query(bioinformatics_output=bio_output_dict)
    # epi_context is a dict matching NormalizedEpiObject's schema.
"""

import argparse
import json
import logging
import re
import sys
from pathlib import Path
from typing import Any, Optional

from pydantic import ValidationError

from intelligence_engine.data_engine.analytics.schemas import (
    ClinicalFeature,
    DemographicSummary,
    DiagnosticMethod,
    GenomicLink,
    InterventionRecord,
    KnowledgeAssertion,
    MolecularEpidemiology,
    NormalizedEpiObject,
    OutbreakRecord,
    PathogenProfile,
    Reference,
    SurveillanceSummary,
    TherapeuticProduct,
    TransmissionParams,
    VaccineProduct,
)
from intelligence_engine.data_engine.llm_querying.llm_client import LLMClient
from intelligence_engine.data_engine.online_querying.ddg_epi_search import DuckDuckGoEpiFetcher, _COUNTRY_LIST

log = logging.getLogger(__name__)


class EpiQueryEngine:
    """Fetches epidemiological data and populates a normalized entity-based schema.

    The engine reads the bioinformatics output JSON, extracts relevant fields,
    fetches real epidemiological reports from trusted sources, and asks the
    LLM (or a deterministic fallback) to populate ``NormalizedEpiObject``'s
    entity sections (pathogen profile, outbreaks, transmission, etc.).

    Every extracted value must be grounded in the fetched trusted-source text
    or the curated source registry -- no training-data leakage. Fields with
    no supporting evidence are left null/empty rather than guessed.
    """

    def __init__(self, llm_client: Optional[LLMClient] = None,
                 data_fetcher: Optional[DuckDuckGoEpiFetcher] = None,
                 use_llm: bool = True):
        """
        Args:
            llm_client: Optional LLMClient. If None and use_llm is True, one is created.
            data_fetcher: Optional EpiDataFetcher.
            use_llm: If True (default), structure epidemiological answers with the LLM
                     using only fetched online data. If False, use deterministic Python extraction.
        """
        self.use_llm = use_llm
        if llm_client:
            self.llm = llm_client
        elif use_llm:
            try:
                self.llm = LLMClient()
            except Exception as e:
                log.warning(
                    f"Local Ollama LLM unavailable ({e}); falling back to deterministic mode."
                )
                self.use_llm = False
                self.llm = None
        else:
            self.llm = None
        self.fetcher = data_fetcher or DuckDuckGoEpiFetcher(
            min_credibility=5,
            max_results_per_query=20,
            max_total_results=30,
        )

    def query(self, bioinformatics_output: dict,
              local_db_results: Optional[dict] = None) -> dict:
        """Fetch real epidemiological data and populate the normalized entity schema.

        The LLM (or the deterministic fallback, in --no-llm mode) is grounded
        strictly in the fetched online text from official/trusted public-health
        sources plus the curated source registry. It must not use training
        data or prior knowledge.

        Args:
            bioinformatics_output: The bioinformatics pipeline output JSON
                                   (stage9_normalised_output or full pipeline output).
            local_db_results: Results from the deterministic SQL queries
                              (layer1-layer5). Optional context.

        Returns:
            A dict matching ``NormalizedEpiObject``'s schema (metadata,
            pathogen_profile, molecular_epidemiology, outbreaks, transmission,
            demographics, clinical, interventions, diagnostics, therapeutics,
            vaccines, surveillance, genomic_links, knowledge_assertions,
            references).
        """
        stage9 = self._extract_stage9(bioinformatics_output)
        db_summary = local_db_results or {}
        context = self._build_context(stage9, db_summary)

        # ── Fetch real data from the internet ──
        print("Searching DuckDuckGo for epidemiological data from trusted sources (WHO, CDC, UN, MoH, etc.)...", flush=True)
        log.info("Fetching epidemiological data via DuckDuckGo search from trusted sources "
                 "(WHO, CDC, UN agencies, national MoH, Africa CDC, ECDC)")
        fetched_data = self.fetcher.fetch_all(
            species_name=context.get("species_name") or "unknown",
            lineage=context.get("lineage"),
            country=context.get("country"),
            collection_date=context.get("collection_date"),
        )

        # ── Keep only official/trusted sources before sending to the LLM ──
        fetched_data = self._filter_trusted_sources(fetched_data)
        trusted_items = fetched_data.get("items", [])
        fetched_context = self.fetcher.format_for_llm(fetched_data, max_chars=12000)

        # Deterministic fallback, built from the fetched/registry items directly.
        fallback = self._build_deterministic_fallback(trusted_items, context)

        if not self.use_llm or self.llm is None:
            print(f"Fetched {len(trusted_items)} reports. Running in free (no-LLM) mode with deterministic extraction...", flush=True)
            entity_dict = fallback
            provider, model, structuring_method = "deterministic (no LLM)", "N/A", "deterministic_extraction"
        else:
            print(f"Fetched {len(trusted_items)} reports from trusted sources. Calling LLM (this can take 30-120 seconds)...", flush=True)
            prompt = self._build_entity_prompt(context, fetched_context)
            log.info(f"Calling LLM for entity-based epidemiological extraction (prompt: {len(prompt)} chars)")
            try:
                llm_response = self.llm.query_epidemiology(prompt)
                if not isinstance(llm_response, dict):
                    llm_response = {}
                print("LLM response received. Merging with deterministic fallback...", flush=True)
                entity_dict = self._merge_entity_dicts(llm_response, fallback)
                provider, model, structuring_method = self.llm.provider, self.llm._get_model(), "llm_entity_extraction"
            except Exception as e:
                log.warning(f"LLM call failed ({e}); using deterministic fallback.")
                print(f"LLM call failed ({e}); using deterministic fallback.", flush=True)
                entity_dict = fallback
                provider, model, structuring_method = "deterministic (no LLM)", "N/A", "deterministic_extraction"

        entity_dict["metadata"] = {
            "query_timestamp": fetched_data.get("_summary", {}).get("fetch_timestamp"),
            "species": context.get("species_name"),
            "lineage": context.get("lineage"),
            "country": context.get("country"),
            "collection_date": context.get("collection_date"),
            "total_reports_fetched": len(trusted_items),
            "api_records_fetched": fetched_data.get("_summary", {}).get("api_records_fetched"),
            "provider": provider,
            "model": model,
            "structuring_method": structuring_method,
        }
        entity_dict["references"] = [r.model_dump() for r in self._build_references(fetched_data)]

        return self._validate_entity_dict(entity_dict).model_dump()

    def query_from_file(self, input_path: str,
                        local_db_results: Optional[dict] = None) -> dict:
        """Convenience method: read the epi engine input from a JSON file.

        Accepts two shapes:
          1. A workflow input file with:
               { "bioinformatics_output_path": "..." }
             or:
               { "bioinformatics_output": {...} }
          2. A raw bioinformatics output JSON.

        Args:
            input_path: Path to the workflow input or bioinformatics output JSON.
            local_db_results: Optional DB query results dict.

        Returns:
            A dict matching ``NormalizedEpiObject``'s schema.
        """
        path = Path(input_path)
        if not path.exists():
            raise FileNotFoundError(f"Input not found: {path}")

        with open(path) as f:
            data = json.load(f)

        if "bioinformatics_output" in data or "bioinformatics_output_path" in data:
            bio_output = data.get("bioinformatics_output")
            if not bio_output:
                bio_path = path.parent / data["bioinformatics_output_path"]
                if not bio_path.exists():
                    raise FileNotFoundError(f"bioinformatics_output_path not found: {bio_path}")
                with open(bio_path) as f:
                    bio_output = json.load(f)
            return self.query(bio_output, local_db_results)

        return self.query(data, local_db_results)

    def _extract_stage9(self, bio_output: dict) -> dict:
        """Extract stage9_normalised_output from full pipeline output,
        or use the dict directly if it's already stage9."""
        if "stage9_normalised_output" in bio_output:
            return bio_output["stage9_normalised_output"]
        if "intelligence_engine_queries" in bio_output:
            return bio_output.get("stage9_normalised_output", bio_output)
        return bio_output

    def _build_context(self, stage9: dict, db_summary: dict) -> dict:
        """Build the context dict used for template substitution and condition evaluation."""
        metadata = stage9.get("metadata", {})

        gq = stage9.get("genome_quality", {})
        if isinstance(gq, dict):
            quality_flag = gq.get("flag")
        else:
            quality_flag = gq

        context = {
            "sample_id": stage9.get("sample_id"),
            "pathogen_id": stage9.get("pathogen", stage9.get("pathogen_id")),
            "species_id": stage9.get("species_id"),
            "species_name": stage9.get("species", stage9.get("species_name")),
            "lineage": stage9.get("lineage"),
            "clade": stage9.get("clade"),
            "country": metadata.get("country"),
            "collection_date": metadata.get("collection_date"),
            "genome_quality_flag": quality_flag,
            "mutations_detected": stage9.get("mutations", []),
            "negative_findings": stage9.get("negative_findings", []),
            "pathogen_family": stage9.get("pathogen_family"),
            "pathogen_genus": stage9.get("pathogen_genus"),
            "lineage_last_seen": db_summary.get("lineage_last_seen_date"),
            "lineage_countries_in_db": db_summary.get("lineage_countries_in_db", []),
            "lineage_not_in_db_for": db_summary.get("lineage_not_in_db_for", []),
            "variant_frequencies": db_summary.get("variant_frequencies", {}),
            "phenotype_associations": db_summary.get("phenotype_associations", {}),
        }

        for k, v in context.items():
            if v == "" or v == "unknown" or v == "null":
                context[k] = None

        return context

    def _build_entity_prompt(self, ctx: dict, fetched_context: str = "") -> str:
        """Assemble the full LLM prompt asking for the entity-based schema.

        Args:
            ctx: Bioinformatics + DB context dict
            fetched_context: Formatted string of real data fetched from
                             trusted sources (WHO, ProMED, CDC, etc.)
        """
        bio_lines = [
            "## BIOINFORMATICS SUMMARY",
            f"Sample ID: {ctx.get('sample_id', 'unknown')}",
            f"Pathogen: {ctx.get('pathogen_id', 'unknown')}",
            f"Species: {ctx.get('species_name', 'unknown')} ({ctx.get('species_id', 'unknown')})",
        ]
        if ctx.get("lineage"):
            bio_lines.append(f"Lineage: {ctx['lineage']}")
        if ctx.get("clade"):
            bio_lines.append(f"Clade: {ctx['clade']}")
        if ctx.get("country"):
            bio_lines.append(f"Collection country: {ctx['country']}")
        if ctx.get("collection_date"):
            bio_lines.append(f"Collection date: {ctx['collection_date']}")
        if ctx.get("genome_quality_flag"):
            bio_lines.append(f"Genome quality: {ctx['genome_quality_flag']}")

        mutations = ctx.get("mutations_detected", [])
        if mutations:
            if isinstance(mutations, list):
                mut_str = ", ".join(
                    m.get("hgvs_p", str(m)) if isinstance(m, dict) else str(m)
                    for m in mutations
                )
            else:
                mut_str = str(mutations)
            bio_lines.append(f"Mutations detected: {mut_str}")

        neg = ctx.get("negative_findings", [])
        if neg:
            if isinstance(neg, list):
                neg_str = ", ".join(
                    n.get("mutation", str(n)) if isinstance(n, dict) else str(n)
                    for n in neg
                )
            else:
                neg_str = str(neg)
            bio_lines.append(f"Negative findings (NOT detected): {neg_str}")

        schema_lines = [
            "## OUTPUT SCHEMA",
            "Return a single JSON object with EXACTLY these top-level keys. Each key's",
            "value must be an ARRAY of objects (or null for singleton keys), using ONLY",
            "the exact field names listed. Leave a field null, or return an empty array,",
            "when the fetched text has no supporting evidence for it -- do not guess.",
            "",
            "- pathogen_profile (single object, not array): species, pathogen_family, "
            "pathogen_genus, reservoir, host, first_documented_year, first_documented_location, "
            "pathogenic_in_humans, confirmed_outbreaks_count, source_url",
            "- molecular_epidemiology (array): lineage, strain, genotype, clade, country, year, "
            "cases, deaths, cfr, key_transmission_features, source_url "
            "(strain/genotype/clade-level case-death burden -- NOT individual mutations, "
            "those go in genomic_links)",
            "- outbreaks (array): pathogen, lineage, country, admin_region, start_date, end_date, "
            "cases, deaths, cfr, source_url, reporting_agency",
            "- transmission (single object, not array): r0_low, r0_high, "
            "incubation_period_days_low, incubation_period_days_high, serial_interval_days, "
            "transmission_route, host, reservoir, vector, source_url",
            "- demographics (array): age_group, sex, occupation, population_affected, "
            "risk_group, exposure_history, setting, host_species, case_count, risk_factor, "
            "notes, source_url",
            "- clinical (array): feature, lineage, severity, frequency, notes, source_url",
            "- interventions (array): intervention_type, name, status, effectiveness, notes, source_url",
            "- diagnostics (array): method, type, target, notes, source_url",
            "- therapeutics (array): product, status, effectiveness, notes, source_url",
            "- vaccines (array): product, status, effectiveness, notes, source_url",
            "- surveillance (array): country, region, first_documented, reservoir, "
            "seroprevalence_pct, seroprevalence_population, source_url",
            "- genomic_links (array): lineage, clade, mutations, genomic_accession, source_url",
            "- knowledge_assertions (array): claim, source_url, reporting_agency, confidence "
            "(0-1), evidence_level",
        ]

        # Assemble full prompt: bio summary + fetched data + schema instructions
        sections = [
            "\n".join(bio_lines),
        ]

        if fetched_context:
            sections.append(fetched_context)
        else:
            sections.append("## FETCHED EPIDEMIOLOGICAL DATA\n(No data could be fetched from trusted sources. "
                            "Return empty arrays/nulls for all sections.)")

        sections.append("\n".join(schema_lines))
        sections.append(
            "IMPORTANT:\n"
            "1. Extract and populate the schema above ONLY from the fetched epidemiological text.\n"
            "2. Do NOT use your training data, memory, prior knowledge, or general world facts.\n"
            "3. Verify that each source you cite is an official or trusted public-health page "
            "(WHO, CDC, Africa CDC, ECDC, GOARN, national public-health agencies, UNICEF, or recognised international organisations). "
            "Ignore any claim from an unofficial, unknown, or low-credibility source.\n"
            "4. Populate ONLY the sections and fields listed above. Do not add extra top-level keys "
            "or invent new field names.\n"
            "5. Do NOT include narrative reasoning, temporal assessments of the sample's collection "
            "date, novelty judgments about the pathogen, or risk interpretation -- only objective, "
            "sourced epidemiological facts.\n"
            "6. Be thorough: extract every relevant fact from the fetched reports into the "
            "matching section. Only leave a section empty/null when the fetched text genuinely "
            "contains no information for it.\n"
            "7. Do NOT include prose, markdown, or explanation outside the JSON.\n"
            "8. Extract numbers from the text when present (cases, deaths, CFR, dates, R0, "
            "incubation period, serial interval).\n"
            "9. Cite the source URL in each row/object when one is available in the fetched data.\n"
            "10. The fetched text may contain pipe-separated tables (e.g. 'Country | Cases | Deaths'). "
            "Read each row and extract the numbers.\n"
            "11. If a number is approximate (e.g. 'about 150', '>1000'), record the number and note "
            "'approximate' in a notes field if one exists.\n"
            "12. For any claim that is well-supported by the text but does not fit a strict tabular "
            "field above, add it to knowledge_assertions with its source_url and a confidence score.\n"
            "13. Capture the full epidemiological picture, not just outbreak counts: classical "
            "epidemiology (who/what/when/where), molecular epidemiology (strain/genotype/clade "
            "burden), demographics and exposure history, seroepidemiology (seroprevalence_pct in "
            "surveillance, where reported), One Health signals (animal hosts, reservoirs, vectors), "
            "and environmental/ecological drivers (climate, seasonality, land use) whenever the "
            "fetched text supports them."
        )

        return "\n\n".join(sections)

    def _merge_entity_dicts(self, llm_dict: dict, fallback_dict: dict) -> dict:
        """Merge the LLM's entity dict with the deterministic fallback.

        For list-valued sections, prefer the LLM's list when it contains at
        least one row with a non-null, non-source_url value; otherwise use
        the fallback. The 'outbreaks' section is merged (not replaced) since
        the fallback often has more complete case/death numbers from the
        curated registry than the LLM's free-text extraction.
        """
        merged = {}
        list_sections = (
            "molecular_epidemiology", "demographics", "clinical", "interventions",
            "diagnostics", "therapeutics", "vaccines", "surveillance",
            "genomic_links", "knowledge_assertions",
        )

        def _has_signal(rows):
            return isinstance(rows, list) and any(
                isinstance(r, dict) and any(
                    k != "source_url" and v not in (None, "")
                    for k, v in r.items()
                )
                for r in rows
            )

        for section in list_sections:
            llm_rows = llm_dict.get(section)
            merged[section] = llm_rows if _has_signal(llm_rows) else fallback_dict.get(section, [])

        # Outbreaks: merge LLM + fallback candidates and deduplicate.
        llm_outbreaks = llm_dict.get("outbreaks") if isinstance(llm_dict.get("outbreaks"), list) else []
        fallback_outbreaks = fallback_dict.get("outbreaks", [])
        merged["outbreaks"] = self._deduplicate_outbreak_rows(
            [r for r in llm_outbreaks if isinstance(r, dict)] + fallback_outbreaks
        )

        # Singleton sections: prefer LLM's if it has any non-null value, else fallback.
        for section in ("pathogen_profile", "transmission"):
            llm_val = llm_dict.get(section)
            if isinstance(llm_val, dict) and any(v not in (None, "") for v in llm_val.values()):
                merged[section] = llm_val
            else:
                merged[section] = fallback_dict.get(section)

        return merged

    def _validate_entity_dict(self, entity_dict: dict) -> NormalizedEpiObject:
        """Validate each entity section independently, dropping invalid rows
        rather than failing the whole object on one bad row."""
        def _safe_list(model_cls, rows):
            out = []
            for row in rows or []:
                if not isinstance(row, dict):
                    continue
                try:
                    out.append(model_cls(**row).model_dump())
                except ValidationError as e:
                    log.warning(f"Dropping invalid {model_cls.__name__} row: {e}")
            return out

        def _safe_singleton(model_cls, value):
            if not isinstance(value, dict):
                return None
            try:
                return model_cls(**value).model_dump()
            except ValidationError as e:
                log.warning(f"Dropping invalid {model_cls.__name__} object: {e}")
                return None

        safe = {
            "metadata": entity_dict.get("metadata", {}),
            "pathogen_profile": _safe_singleton(PathogenProfile, entity_dict.get("pathogen_profile")),
            "molecular_epidemiology": _safe_list(MolecularEpidemiology, entity_dict.get("molecular_epidemiology")),
            "outbreaks": _safe_list(OutbreakRecord, entity_dict.get("outbreaks")),
            "transmission": _safe_singleton(TransmissionParams, entity_dict.get("transmission")),
            "demographics": _safe_list(DemographicSummary, entity_dict.get("demographics")),
            "clinical": _safe_list(ClinicalFeature, entity_dict.get("clinical")),
            "interventions": _safe_list(InterventionRecord, entity_dict.get("interventions")),
            "diagnostics": _safe_list(DiagnosticMethod, entity_dict.get("diagnostics")),
            "therapeutics": _safe_list(TherapeuticProduct, entity_dict.get("therapeutics")),
            "vaccines": _safe_list(VaccineProduct, entity_dict.get("vaccines")),
            "surveillance": _safe_list(SurveillanceSummary, entity_dict.get("surveillance")),
            "genomic_links": _safe_list(GenomicLink, entity_dict.get("genomic_links")),
            "knowledge_assertions": _safe_list(KnowledgeAssertion, entity_dict.get("knowledge_assertions")),
            "references": _safe_list(Reference, entity_dict.get("references")),
        }
        return NormalizedEpiObject(**safe)

    def _all_text(self, ob: dict) -> str:
        """Return all fetched text for an outbreak as a single string."""
        parts = [
            ob.get("title", ""),
            ob.get("description", ""),
            ob.get("content_excerpt", ""),
            ob.get("summary", ""),
            ob.get("epidemiology", ""),
            ob.get("overview", ""),
            ob.get("content", ""),
        ]
        return " ".join(str(p) for p in parts if p)

    @staticmethod
    def _is_duplicate_outbreak(existing_rows: list, start_year: Any, country: Any,
                               cases: Any, deaths: Any) -> bool:
        """Return True if a matching outbreak is already in *existing_rows*.

        Two records are considered duplicates when they share the same start year
        and country and their reported case/death numbers are within a small
        tolerance.  This merges near-identical reports from different providers
        while preserving genuinely separate outbreaks in the same country/year.
        """
        if not start_year or not country:
            return False
        try:
            year = int(start_year)
            c = int(cases) if cases is not None else None
            d = int(deaths) if deaths is not None else None
        except (ValueError, TypeError):
            return False
        for r in existing_rows:
            if r.get("country") != country:
                continue
            try:
                ry = int(str(r.get("start_date"))[:4]) if r.get("start_date") else None
                if ry is None:
                    continue
                if ry != year:
                    continue
                rc = int(r.get("cases")) if r.get("cases") is not None else None
                rd = int(r.get("deaths")) if r.get("deaths") is not None else None
            except (ValueError, TypeError):
                continue
            if rc is None or c is None:
                # If numbers are missing, treat same year/country as duplicate only when exact
                if rc == c and rd == d:
                    return True
                continue
            # Tolerance: within 2% of the larger outbreak.  This keeps genuinely
            # separate small outbreaks distinct while merging near-identical reports
            # from different providers.
            c_diff = abs(rc - c)
            d_diff = abs(rd - d) if rd is not None and d is not None else 0
            c_tol = max(0, int(rc * 0.02)) if rc is not None else 0
            d_tol = max(0, int(rd * 0.02)) if rd is not None else 0
            if c_diff <= c_tol and d_diff <= d_tol:
                return True
        return False

    def _deduplicate_outbreak_rows(self, rows: list) -> list:
        """Return *rows* with near-duplicate same-year/country outbreaks merged."""
        out = []
        for r in rows:
            start_year = r.get("start_date")
            if isinstance(start_year, str) and start_year:
                m = re.search(r"(\d{4})", start_year)
                start_year = m.group(1) if m else start_year[:4]
            if self._is_duplicate_outbreak(out, start_year, r.get("country"),
                                           r.get("cases"), r.get("deaths")):
                continue
            out.append(r)
        return out

    def _extract_outbreak_records_from_text(self, text: str, url: str) -> list:
        """Parse CDC/WHO-style outbreak pages into individual outbreak records.

        Splits by year headings (### 2022) and country headings (#### Country)
        and extracts cases, deaths, and CFR from each section.
        """
        if not text:
            return []
        text = re.sub(r"[*_]+", " ", text)
        records = []
        # Split by year headings (### 2022, ### 2021, etc.)
        sections = re.split(r"\n(?=###\s+\d{4}\b)", text)
        for section in sections:
            year_match = re.search(r"###\s+(\d{4})\b", section)
            year = year_match.group(1) if year_match else None
            # Split by country headings (#### Country)
            subsections = re.split(r"\n(?=####\s+)", section)
            for sub in subsections:
                heading_match = re.search(r"####\s+(.+)", sub)
                if not heading_match:
                    continue
                heading = heading_match.group(1).strip()
                heading = re.sub(r"\s*\([^)]*\)$", "", heading)
                country = self._extract_country_from_text(heading)
                if not country:
                    country = self._extract_country_from_text(sub[:300])
                cases = self._extract_number_from_text(sub, "cases")
                deaths = self._extract_number_from_text(sub, "deaths")
                if not country and not cases and not deaths:
                    continue
                cfr = self._calculate_cfr(sub)
                records.append({
                    "date": year,
                    "country": country,
                    "cases": cases,
                    "deaths": deaths,
                    "cfr": cfr,
                    "source_url": url,
                })
        return records

    @staticmethod
    def _item_data(ob: dict, *keys) -> Any:
        """Return the first matching key from the structured ``data`` payload.

        The source registry stores normalised fields in ``item["data"]``.
        This helper lets the deterministic extractor use those values first and
        fall back to ``None`` only when absent.
        """
        data = ob.get("data") if isinstance(ob.get("data"), dict) else {}
        for k in keys:
            if k in data and data[k] not in (None, ""):
                return data[k]
        return None

    def _build_deterministic_fallback(self, outbreaks: list, ctx: dict) -> dict:
        """Build the entity-based fallback dict from Python-extracted fetched data.

        This deterministic fallback reuses the fetched/registry outbreak
        reports and maps them directly onto ``NormalizedEpiObject``'s entity
        sections. Values that are not present in the fetched text are left
        null/empty rather than guessed.
        """
        pathogen = ctx.get("species_name") or "unknown"
        lineage = ctx.get("lineage")

        # Keep only fetched items that match the detected species before any
        # downstream extraction. This prevents a corrupted upstream registry
        # from polluting outbreaks, surveillance, demographics, etc.
        filtered_outbreaks = [
            ob for ob in outbreaks
            if self._matches_target_species(ob, ctx)
        ]
        if len(filtered_outbreaks) < len(outbreaks):
            log.info(
                "Filtered %d fetched item(s) that do not match target species '%s'",
                len(outbreaks) - len(filtered_outbreaks),
                pathogen,
            )
        outbreaks = filtered_outbreaks

        # ── outbreaks: one row per individual outbreak report ──
        outbreak_rows = []
        for ob in outbreaks:
            if self._item_data(ob, "item_type") == "strain_profile":
                continue
            text = self._all_text(ob)
            date = ob.get("date", "")
            if isinstance(date, str) and "T" in date:
                date = date.split("T")[0]
            cases = self._item_data(ob, "cases", "total_cases") or self._extract_number_from_text(text, "cases")
            deaths = self._item_data(ob, "deaths", "total_deaths") or self._extract_number_from_text(text, "deaths")
            cfr = self._item_data(ob, "cfr") or self._calculate_cfr(text)
            country = self._item_data(ob, "country") or self._extract_country_from_text(text)
            start_year = self._item_data(ob, "year_start") or (date[:4] if date else None)
            if self._is_duplicate_outbreak(outbreak_rows, start_year, country, cases, deaths):
                continue
            outbreak_rows.append({
                "pathogen": pathogen,
                "lineage": lineage,
                "country": country,
                "start_date": date or None,
                "cases": cases,
                "deaths": deaths,
                "cfr": cfr,
                "source_url": ob.get("url", ""),
                "evidence_level": "structured_api" if self._item_data(ob, "cases", "total_cases") else "text_regex",
            })

        # ── molecular_epidemiology: lineage-specific profile ──
        molecular_epidemiology = []
        if lineage:
            lineage_lower = lineage.lower()
            for ob in outbreaks:
                text = self._all_text(ob)
                strain = self._item_data(ob, "strain", "abbreviation", "species")
                matches = lineage_lower in text.lower()
                if not matches and strain:
                    matches = lineage_lower in str(strain).lower()
                if not matches:
                    continue
                date = ob.get("date", "")
                if isinstance(date, str) and "T" in date:
                    date = date.split("T")[0]
                cases = self._item_data(ob, "cases", "total_cases") or self._extract_number_from_text(text, "cases")
                deaths = self._item_data(ob, "deaths", "total_deaths") or self._extract_number_from_text(text, "deaths")
                cfr = self._item_data(ob, "cfr") or self._calculate_cfr(text)
                country = self._item_data(ob, "country") or self._extract_country_from_text(text)
                genotype = self._item_data(ob, "genotype")
                clade = self._item_data(ob, "clade")
                molecular_epidemiology.append({
                    "lineage": lineage,
                    "strain": strain if strain else None,
                    "genotype": genotype,
                    "clade": clade,
                    "country": country,
                    "year": date[:4] if date else None,
                    "cases": cases,
                    "deaths": deaths,
                    "cfr": cfr,
                    "source_url": ob.get("url", ""),
                })

        # ── pathogen_profile: aggregated from strain_profile registry items ──
        pathogen_profile = None
        for ob in outbreaks:
            if self._item_data(ob, "item_type") != "strain_profile":
                continue
            pathogen_profile = {
                "species": self._item_data(ob, "species") or pathogen,
                "pathogen_family": ctx.get("pathogen_family"),
                "pathogen_genus": ctx.get("pathogen_genus"),
                "pathogenic_in_humans": self._item_data(ob, "pathogenic_in_humans"),
                "first_documented_year": self._item_data(ob, "discovered_year"),
                "first_documented_location": self._item_data(ob, "discovered_location"),
                "confirmed_outbreaks_count": self._item_data(ob, "confirmed_outbreaks_count"),
                "source_url": ob.get("url", ""),
            }
            break

        # ── entity sections derived from text-mining the fetched reports ──
        demographics = self._extract_demographics(outbreaks)
        transmission_rows = self._extract_transmission_params(outbreaks)
        transmission = self._fold_transmission_params(transmission_rows)
        interventions = self._extract_interventions(outbreaks)
        clinical = self._extract_clinical_features(outbreaks, lineage or "")
        surveillance = self._extract_surveillance(outbreaks)
        diagnostics = self._extract_diagnostics(outbreaks)
        vaccines, therapeutics = self._extract_vaccine_therapeutics(outbreaks)

        return {
            "pathogen_profile": pathogen_profile,
            "molecular_epidemiology": molecular_epidemiology,
            "outbreaks": outbreak_rows,
            "transmission": transmission,
            "demographics": demographics,
            "clinical": clinical,
            "interventions": interventions,
            "diagnostics": diagnostics,
            "therapeutics": therapeutics,
            "vaccines": vaccines,
            "surveillance": surveillance,
            "genomic_links": [],
            "knowledge_assertions": [],
        }

    def _extract_number_from_text(self, text: str, kind: str) -> Optional[int]:
        """Extract case/death counts from fetched text, preferring totals."""
        if not text:
            return None
        text = text.lower()
        # Remove markdown bold/italic so '**164**' or '__164__' can be parsed
        text = re.sub(r"[*_]+", " ", text)

        def _parse_int(s: str) -> Optional[int]:
            try:
                return int(s.replace(",", "").replace(" ", ""))
            except (ValueError, AttributeError):
                return None

        if kind == "cases":
            patterns = [
                r"cumulative\s+total\s+of\s+([\d, ]{1,11})\s+(?:clinical\s+|confirmed\s+|probable\s+)?cases",
                r"total\s+(?:number\s+)?of\s+([\d, ]{1,11})\s+(?:clinical\s+|confirmed\s+|probable\s+)?cases",
                r"total\s+of\s+([\d, ]{1,11})\s+(?:clinical\s+|confirmed\s+|probable\s+|new\s+|probable\s+and\s+confirmed\s+)?cases",
                r"([\d, ]{1,11})\s+(?:probable\s+and\s+confirmed\s+|confirmed\s+|probable\s+|suspect\s+|clinical\s+|new\s+)?cases",
                r"([\d, ]{1,11})\s+cases\s*\(",
                r"cases\s*[:\-]\s*([\d, ]{1,11})",
                r"case\s+count\s*[:\-]\s*([\d, ]{1,11})",
                r"\b([\d, ]{1,11})\s+confirmed\s+cases?\b",
            ]
        else:
            patterns = [
                r"deaths\s*[:\-]\s*([\d, ]{1,11})",
                r"(?:including|and)\s+([\d, ]{1,11})\s*deaths",
                r"([\d, ]{1,11})\s*deaths",
                r"([\d, ]{1,11})\s*have\s*been\sfatal",
                r"([\d, ]{1,11})\s*fatal",
            ]

        for pattern in patterns:
            m = re.search(pattern, text, re.IGNORECASE)
            if m:
                val = _parse_int(m.group(1))
                if val is not None:
                    return val
        return None

    def _calculate_cfr(self, text: str) -> Optional[str]:
        """Return a CFR string if both cases and deaths are found, else None."""
        if not text:
            return None
        cases = self._extract_number_from_text(text, "cases")
        deaths = self._extract_number_from_text(text, "deaths")
        if cases and deaths and cases > 0:
            try:
                return f"{round(deaths / cases * 100, 1)}%"
            except Exception:
                return None
        return None

    def _extract_country_from_text(self, text: str) -> Optional[str]:
        """Return the first recognised country mentioned in the text."""
        if not text:
            return None
        text_lower = text.lower()
        earliest = None
        earliest_pos = None
        for country in _COUNTRY_LIST:
            if len(country) <= 3:
                # Short aliases need word boundaries so 'us' doesn't match 'virus'
                pattern = rf"\b{re.escape(country.lower())}\b"
            else:
                # Longer names should not be part of a longer word
                pattern = rf"(?<![a-zA-Z]){re.escape(country.lower())}(?![a-zA-Z])"
            for m in re.finditer(pattern, text_lower):
                pos = m.start()
                if earliest_pos is None or pos < earliest_pos:
                    earliest_pos = pos
                    if country == "DRC":
                        earliest = "Democratic Republic of the Congo"
                    elif country == "Congo, The Democratic Republic of the":
                        earliest = "Democratic Republic of the Congo"
                    else:
                        earliest = country
                break
        return earliest

    def _extract_demographics(self, outbreaks: list) -> list:
        """Extract demographic groups mentioned in the fetched text."""
        rows = []
        seen = set()
        for ob in outbreaks:
            text = self._all_text(ob).lower()
            url = ob.get("url", "")
            checks = [
                (("health care worker", "healthcare worker", "health worker", "health workers"), "healthcare_workers", "occupational exposure", None),
                (("child", "children", "under 5", "under five"), "children", None, "children"),
                (("pregnant", "pregnancy"), "pregnant_women", None, None),
                (("elderly", "older adult", "old age"), "elderly", None, "elderly"),
            ]
            for terms, key, risk, age_group in checks:
                if any(term in text for term in terms) and key not in seen:
                    seen.add(key)
                    rows.append({
                        "age_group": age_group,
                        "sex": "female" if key == "pregnant_women" else None,
                        "occupation": "healthcare worker" if key == "healthcare_workers" else None,
                        "population_affected": age_group,
                        "risk_group": "healthcare workers" if key == "healthcare_workers" else None,
                        "risk_factor": risk,
                        "notes": f"{key.replace('_', ' ')} mentioned in fetched source",
                        "source_url": url,
                    })

            # Exposure history
            exposure_checks = [
                (("animal contact", "contact with animal", "bushmeat", "handling of infected animals"), "animal contact"),
                (("funeral", "burial"), "funeral/burial attendance"),
                (("nosocomial", "healthcare exposure", "hospital-acquired"), "healthcare exposure"),
                (("travel history", "recent travel", "traveler", "travellers"), "travel"),
                (("mosquito bite", "tick bite", "vector exposure", "insect bite"), "vector exposure"),
            ]
            for terms, exposure in exposure_checks:
                key = f"exposure_{exposure}"
                if any(term in text for term in terms) and key not in seen:
                    seen.add(key)
                    rows.append({
                        "age_group": None, "sex": None, "occupation": None,
                        "population_affected": None, "risk_group": None,
                        "exposure_history": exposure, "setting": None, "host_species": None,
                        "risk_factor": None,
                        "notes": f"{exposure} reported in fetched source",
                        "source_url": url,
                    })

            # Setting of infection
            setting_checks = [
                (("healthcare facility", "hospital setting", "clinic setting"), "healthcare facility"),
                (("household transmission", "household contact", "within the household"), "household"),
                (("refugee camp", "displacement camp", "idp camp"), "refugee/displacement camp"),
                (("community transmission", "community setting"), "community"),
            ]
            for terms, setting in setting_checks:
                key = f"setting_{setting}"
                if any(term in text for term in terms) and key not in seen:
                    seen.add(key)
                    rows.append({
                        "age_group": None, "sex": None, "occupation": None,
                        "population_affected": None, "risk_group": None,
                        "exposure_history": None, "setting": setting, "host_species": None,
                        "risk_factor": None,
                        "notes": f"{setting} setting reported in fetched source",
                        "source_url": url,
                    })

            # Affected host species (One Health / zoonotic spillover)
            host_species_checks = [
                (("fruit bat", "fruit bats"), "fruit bats"),
                (("nonhuman primate", "non-human primate", "chimpanzee", "gorilla"), "non-human primates"),
                (("cattle", "livestock"), "cattle/livestock"),
                (("pig", "swine"), "pigs"),
                (("poultry", "chicken", "bird flu"), "poultry"),
                (("camel",), "camels"),
            ]
            for terms, species in host_species_checks:
                key = f"host_species_{species}"
                if any(term in text for term in terms) and key not in seen:
                    seen.add(key)
                    rows.append({
                        "age_group": None, "sex": None, "occupation": None,
                        "population_affected": None, "risk_group": None,
                        "exposure_history": None, "setting": None, "host_species": species,
                        "risk_factor": None,
                        "notes": f"{species} reported as affected host species in fetched source",
                        "source_url": url,
                    })

            # Male/female require whole-word boundaries so "female" does not match "male"
            if re.search(r"\bmales?\b|\bmen\b", text) and "males" not in seen:
                seen.add("males")
                rows.append({
                    "age_group": None,
                    "sex": "male",
                    "occupation": None,
                    "risk_factor": None,
                    "notes": "mentioned in fetched source",
                    "source_url": url,
                })
            if re.search(r"\bfemales?\b|\bwomen\b", text) and "females" not in seen:
                seen.add("females")
                rows.append({
                    "age_group": None,
                    "sex": "female",
                    "occupation": None,
                    "risk_factor": None,
                    "notes": "mentioned in fetched source",
                    "source_url": url,
                })
            # Look for explicit male/female counts
            m = re.search(r"(\d+)\s*(?:male|men)\s*(?:and|&)?\s*(\d+)\s*(?:female|women)", text)
            if m:
                key = f"gender_split_{m.group(1)}_{m.group(2)}"
                if key not in seen:
                    seen.add(key)
                    rows.append({
                        "age_group": None,
                        "sex": "male/female",
                        "occupation": None,
                        "risk_factor": None,
                        "notes": f"{m.group(1)} male, {m.group(2)} female",
                        "source_url": url,
                    })
        return rows

    def _extract_transmission_params(self, outbreaks: list) -> list:
        """Extract incubation period, serial interval, R0 and transmission routes."""
        rows = []
        seen = set()
        for ob in outbreaks:
            text = self._all_text(ob)
            url = ob.get("url", "")
            lower = text.lower()

            # Prefer structured registry fields
            r0_range = self._item_data(ob, "r0_range")
            if r0_range and r0_range not in seen:
                seen.add(r0_range)
                rows.append({"parameter": "R0", "value": r0_range, "unit": None, "notes": None, "source_url": url})

            transmission = self._item_data(ob, "transmission")
            if transmission:
                key = f"route_{transmission}"
                if key not in seen:
                    seen.add(key)
                    rows.append({"parameter": "main transmission route", "value": transmission, "unit": None, "notes": None, "source_url": url})

            # Incubation period
            m = re.search(r"incubation\s+period(?:\s+is|\s+of)?\s+(?:about|approximately|up\s+to)?\s*(\d+\s*(?:to|-|–)\s*\d+)\s*days?", lower)
            if m:
                key = f"incubation_{m.group(1)}"
                if key not in seen:
                    seen.add(key)
                    rows.append({"parameter": "incubation period", "value": m.group(1).replace("–", "-").strip(), "unit": "days", "notes": None, "source_url": url})
            else:
                m = re.search(r"incubation\s+period(?:\s+is|\s+of)?\s+up\s+to\s+(\d+)\s*days?", lower)
                if m:
                    key = f"incubation_up_to_{m.group(1)}"
                    if key not in seen:
                        seen.add(key)
                        rows.append({"parameter": "incubation period", "value": f"up to {m.group(1)}", "unit": "days", "notes": None, "source_url": url})

            # Serial interval
            m = re.search(r"serial\s+interval(?:\s+is|\s+of)?\s+(?:about|approximately)?\s*(\d+\.?\d*)\s*days?", lower)
            if m:
                key = f"serial_{m.group(1)}"
                if key not in seen:
                    seen.add(key)
                    rows.append({"parameter": "serial interval", "value": m.group(1), "unit": "days", "notes": None, "source_url": url})

            # R0 / reproductive number (text fallback)
            m = re.search(r"(?:r0|reproductive\s+number|reproduction\s+number|basic\s+reproduction\s+number|r_0)\s*[:\s-]\s*(\d+\.?\d*)", lower)
            if m:
                key = f"r0_{m.group(1)}"
                if key not in seen:
                    seen.add(key)
                    rows.append({"parameter": "R0", "value": m.group(1), "unit": None, "notes": None, "source_url": url})

            # Transmission routes (text fallback)
            routes = []
            if "person-to-person" in lower or "direct contact" in lower or "close contact" in lower:
                routes.append("person-to-person via direct/close contact")
            if "bodily fluid" in lower or "body fluid" in lower or "body fluids" in lower:
                routes.append("contact with bodily fluids")
            if "sexual" in lower:
                routes.append("sexual transmission")
            if "airborne" in lower:
                routes.append("airborne (not primary)")
            if "needle" in lower or "sharps" in lower or "inject" in lower:
                routes.append("needlestick / contaminated equipment")
            for route in routes:
                key = f"route_{route}"
                if key not in seen:
                    seen.add(key)
                    rows.append({"parameter": "main transmission route", "value": route, "unit": None, "notes": None, "source_url": url})

            # Vector-borne transmission (One Health)
            m = re.search(r"\b(mosquito(?:es)?|tick(?:s)?|sandfl(?:y|ies)|midges?|fle(?:a|as))\b\s*(?:bite|bites|vector)?", lower)
            if m and ("vector" in lower or "bite" in lower or "borne" in lower):
                vector = m.group(1)
                key = f"vector_{vector}"
                if key not in seen:
                    seen.add(key)
                    rows.append({"parameter": "vector", "value": vector, "unit": None, "notes": None, "source_url": url})
        return rows

    def _fold_transmission_params(self, rows: list) -> Optional[dict]:
        """Fold raw {parameter, value, unit, notes, source_url} rows (as
        produced by ``_extract_transmission_params``) into a single
        TransmissionParams-shaped dict, since transmission is a pathogen-level
        singleton section, not a per-row list."""
        if not rows:
            return None
        tp = {
            "r0_low": None, "r0_high": None,
            "incubation_period_days_low": None, "incubation_period_days_high": None,
            "serial_interval_days": None, "transmission_route": None,
            "host": None, "reservoir": None, "vector": None, "source_url": None,
        }
        for row in rows:
            param = str(row.get("parameter", "")).lower()
            value = row.get("value")
            if "r0" in param or "reproductive" in param:
                m = re.findall(r"[\d.]+", str(value or ""))
                if len(m) >= 2:
                    tp["r0_low"], tp["r0_high"] = float(m[0]), float(m[1])
                elif len(m) == 1:
                    tp["r0_low"] = tp["r0_high"] = float(m[0])
            elif "incubation" in param:
                m = re.findall(r"[\d.]+", str(value or ""))
                if len(m) >= 2:
                    tp["incubation_period_days_low"], tp["incubation_period_days_high"] = float(m[0]), float(m[1])
                elif len(m) == 1:
                    tp["incubation_period_days_low"] = tp["incubation_period_days_high"] = float(m[0])
            elif "serial interval" in param:
                m = re.findall(r"[\d.]+", str(value or ""))
                if m:
                    tp["serial_interval_days"] = float(m[0])
            elif "route" in param and not tp["transmission_route"]:
                tp["transmission_route"] = str(value) if value else None
            elif param == "vector" and not tp["vector"]:
                tp["vector"] = str(value) if value else None
            if row.get("source_url") and not tp["source_url"]:
                tp["source_url"] = row.get("source_url")
        return tp

    def _extract_interventions(self, outbreaks: list) -> list:
        """Extract public health interventions and medical countermeasures."""
        rows = []
        seen = set()
        keyword_map = [
            ("vaccine", "vaccination", "vaccine"),
            ("vaccination", "vaccination", "vaccine"),
            ("ring vaccination", "ring vaccination", "vaccine"),
            ("contact tracing", "contact tracing", "public_health_measure"),
            ("isolation", "isolation", "clinical_management"),
            ("quarantine", "quarantine", "public_health_measure"),
            ("personal protective equipment", "PPE", "PPE"),
            ("ppe", "PPE", "PPE"),
            ("safe burial", "safe burial", "public_health_measure"),
            ("monoclonal", "monoclonal antibody", "treatment"),
            ("therapeutic", "therapeutics", "treatment"),
            ("antiviral", "antiviral", "treatment"),
            ("infection prevention", "infection prevention and control", "ipc"),
            ("hand hygiene", "hand hygiene", "public_health_measure"),
            ("surveillance", "surveillance", "public_health_measure"),
            ("risk communication", "risk communication", "public_health_measure"),
        ]
        for ob in outbreaks:
            text = self._all_text(ob).lower()
            url = ob.get("url", "")
            for keyword, label, itype in keyword_map:
                if keyword in text and label not in seen:
                    seen.add(label)
                    rows.append({
                        "intervention_type": itype,
                        "name": label,
                        "effectiveness": None,
                        "status": "used",
                        "source_url": url,
                    })
        return rows

    def _extract_clinical_features(self, outbreaks: list, lineage: str) -> list:
        """Extract clinical features and severity information."""
        rows = []
        seen = set()
        for ob in outbreaks:
            text = self._all_text(ob).lower()
            url = ob.get("url", "")
            for feature, severity in [
                ("asymptomatic", "mild"),
                ("severe", "severe"),
                ("hospitalized", "severe"),
                ("fever", None),
                ("fatigue", None),
                ("muscle pain", None),
                ("headache", None),
                ("vomiting", None),
                ("diarrhoea", None),
                ("diarrhea", None),
                ("haemorrhage", "severe"),
                ("hemorrhage", "severe"),
                ("bleeding", "severe"),
            ]:
                if feature in text and feature not in seen:
                    seen.add(feature)
                    rows.append({
                        "feature": feature,
                        "lineage": lineage or None,
                        "severity": severity,
                        "frequency": None,
                        "notes": "clinical feature mentioned in fetched source",
                        "source_url": url,
                    })
        return rows

    def _region_for_country(self, country: str) -> Optional[str]:
        """Map a country to a broad WHO region."""
        africa = {"Uganda", "Gabon", "Democratic Republic of Congo", "Congo", "Sierra Leone", "Liberia", "Guinea", "Nigeria", "Senegal", "Mali", "Sudan", "South Sudan", "Kenya", "Tanzania", "Cameroon", "Ghana", "Ivory Coast", "Côte d'Ivoire", "Zaire", "Angola", "Zambia", "Congo Republic"}
        europe = {"Spain", "Italy", "Germany", "France", "United Kingdom", "UK", "Russia"}
        asia = {"India", "China", "Japan", "South Korea", "Vietnam", "Thailand", "Philippines", "Indonesia", "Jordan", "Turkey", "Cambodia", "Palestine", "Gaza"}
        americas = {"United States", "USA", "Brazil", "Mexico", "Colombia", "Peru", "Haiti"}
        if country in africa:
            return "Africa"
        if country in europe:
            return "Europe"
        if country in asia:
            return "Asia"
        if country in americas:
            return "Americas"
        if country == "Australia":
            return "Oceania"
        return None

    def _extract_surveillance(self, outbreaks: list) -> list:
        """Extract geographic range, region, reservoir and earliest documented year per country."""
        rows = []
        countries = set()
        # per-country earliest year and URL
        earliest_by_country = {}
        reservoir = None
        seroprevalence_by_country = {}

        for ob in outbreaks:
            text = self._all_text(ob)
            url = ob.get("url", "")
            country = self._item_data(ob, "country") or self._extract_country_from_text(text)
            if country:
                countries.add(country)
            # seroepidemiology: seroprevalence percentage + sampled population
            m = re.search(
                r"seroprevalence\s+(?:of\s+|was\s+)?([\d.]+)\s*%\s*(?:among|in)?\s*([a-z ,]{0,60})",
                text.lower(),
            )
            if m and country and country not in seroprevalence_by_country:
                seroprevalence_by_country[country] = {
                    "pct": float(m.group(1)),
                    "population": m.group(2).strip() or None,
                    "url": url,
                }
            # earliest year per country
            year = self._item_data(ob, "year_start", "discovered_year")
            if year is None:
                date = ob.get("date", "")
                m = re.search(r"(\d{4})", str(date))
                if m:
                    year = int(m.group(1))
            if year and country:
                try:
                    y = int(year)
                    if country not in earliest_by_country or y < earliest_by_country[country]["year"]:
                        earliest_by_country[country] = {"year": y, "url": url}
                except (ValueError, TypeError):
                    pass
            # reservoir clues - only use explicit statements
            if not reservoir:
                lower = text.lower()
                if "fruit bat" in lower or "pteropus" in lower:
                    reservoir = "fruit bats (Pteropodidae)"
                elif re.search(r"\breservoir\b.{0,40}\b(bats?|bat\b)", lower):
                    reservoir = "bats"
                elif re.search(r"\bnatural\s+host\b.{0,40}\b(bats?|bat\b)", lower):
                    reservoir = "bats"

        for country in sorted(countries):
            info = earliest_by_country.get(country, {})
            sero = seroprevalence_by_country.get(country, {})
            rows.append({
                "country": country,
                "region": self._region_for_country(country),
                "first_documented": info.get("year"),
                "reservoir": reservoir,
                "seroprevalence_pct": sero.get("pct"),
                "seroprevalence_population": sero.get("population"),
                "source_url": info.get("url") or sero.get("url"),
            })
        return rows

    def _extract_diagnostics(self, outbreaks: list) -> list:
        """Extract diagnostic and laboratory methods from fetched text."""
        rows = []
        seen = set()
        methods = [
            (("rt-pcr", "rt pcr"), "RT-PCR", "molecular", "viral RNA"),
            (("pcr",), "PCR", "molecular", "viral RNA"),
            (("antigen", "rapid diagnostic"), "antigen detection", "rapid", "viral antigen"),
            (("serology", "serological"), "serology", "immunological", "antibodies"),
            (("elisa",), "ELISA", "immunological", "antibodies"),
            (("sequencing", "genome", "genomic"), "genomic sequencing", "molecular", "whole genome"),
            (("viral culture", "virus isolation"), "viral culture", "isolation", "live virus"),
            (("biosafety", "bsl-4", "bsl4"), "BSL-4/biosafety", "laboratory", "containment"),
            (("biopsy",), "biopsy", "pathology", "tissue sample"),
            (("blood sample", "blood test"), "blood testing", "laboratory", "blood sample"),
        ]
        for ob in outbreaks:
            text = self._all_text(ob).lower()
            url = ob.get("url", "")
            for terms, method, mtype, target in methods:
                if any(term in text for term in terms) and method not in seen:
                    # avoid double-counting PCR when RT-PCR is already present in the same text
                    if method == "PCR" and "RT-PCR" in seen:
                        continue
                    seen.add(method)
                    rows.append({
                        "method": method,
                        "type": mtype,
                        "target": target,
                        "notes": "mentioned in fetched source",
                        "source_url": url,
                    })
        return rows

    def _extract_vaccine_therapeutics(self, outbreaks: list) -> tuple[list, list]:
        """Extract vaccine and therapeutic product names/status from registry data and text.

        Returns:
            (vaccines, therapeutics) -- two separate lists matching the
            VaccineProduct and TherapeuticProduct schemas.
        """
        vaccines = []
        therapeutics = []
        seen = set()
        products = [
            (("ervebo", "rvsv-zebov", "rvsv-zebov-gp"), "Ervebo (rVSV-ZEBOV)", "vaccine", "used", "highly effective"),
            (("rvsv", "vsv-zebov"), "rVSV-based vaccine", "vaccine", "used", None),
            (("mab114", "mab 114"), "mAb114", "therapeutic", "used", "effective"),
            (("regn-eb3", "regn eb3"), "REGN-EB3", "therapeutic", "used", "effective"),
            (("zmapp", "z-mapp"), "ZMapp", "therapeutic", "used", "effective"),
            (("remdesivir",), "Remdesivir", "therapeutic", "used", "limited efficacy"),
            (("favipiravir",), "Favipiravir", "therapeutic", "investigational", None),
            (("monoclonal antibody", "monoclonal antibodies"), "monoclonal antibodies", "therapeutic", "used", "effective"),
            (("vaccine", "vaccination"), "vaccination (unspecified)", "vaccine", "used", None),
            (("therapeutic", "treatment"), "therapeutics (unspecified)", "therapeutic", "used", None),
        ]
        for ob in outbreaks:
            text = self._all_text(ob).lower()
            url = ob.get("url", "")

            # Structured coverage flags from source registry
            vaccine_flag = self._item_data(ob, "vaccine_coverage")
            if vaccine_flag is True and "vaccine coverage" not in seen:
                seen.add("vaccine coverage")
                vaccines.append({
                    "product": "vaccine coverage available",
                    "status": "used",
                    "effectiveness": None,
                    "notes": "supported by species/strain profile",
                    "source_url": url,
                })
            therapeutic_flag = self._item_data(ob, "therapeutic_coverage")
            if therapeutic_flag is True and "therapeutic coverage" not in seen:
                seen.add("therapeutic coverage")
                therapeutics.append({
                    "product": "therapeutic coverage available",
                    "status": "used",
                    "effectiveness": None,
                    "notes": "supported by species/strain profile",
                    "source_url": url,
                })

            for terms, product, ptype, status, effectiveness in products:
                if any(term in text for term in terms) and product not in seen:
                    seen.add(product)
                    row = {
                        "product": product,
                        "status": status,
                        "effectiveness": effectiveness,
                        "notes": "mentioned in fetched source",
                        "source_url": url,
                    }
                    (vaccines if ptype == "vaccine" else therapeutics).append(row)
        return vaccines, therapeutics

    def _build_references(self, fetched_data: dict) -> list:
        """Return the consolidated reference list for the output."""
        refs = []
        for r in self.fetcher.get_reference_urls(fetched_data):
            try:
                refs.append(Reference(
                    source_url=r.get("url"),
                    title=r.get("title"),
                    reporting_agency=r.get("source"),
                ))
            except ValidationError as e:
                log.warning(f"Skipping invalid reference: {e}")
        return refs

    # Official / trusted public-health source keys. Others are dropped before LLM processing.
    MIN_SOURCE_CREDIBILITY = 5

    # Ebolavirus species keywords used to keep fetched items aligned with the
    # sample's detected species. Generic terms like "ebola" or "ebolavirus" are
    # intentionally neutral; only explicit species names/abbreviations drive
    # filtering, so we do not accidentally discard cross-species comparisons.
    _EBOLAVIRUS_SPECIES = {
        "bundibugyo": {"bundibugyo", "bdbv"},
        "zaire": {"zaire", "ebov"},
        "sudan": {"sudan", "sudv"},
        "taï forest": {"taï forest", "tai forest", "tafv"},
        "reston": {"reston", "restv"},
        "bombali": {"bombali", "bomv"},
    }

    def _detect_species(self, text: str) -> set[str]:
        """Return the ebolavirus species explicitly named in *text*."""
        if not text:
            return set()
        lower = str(text).lower()
        matched = set()
        for species, keywords in self._EBOLAVIRUS_SPECIES.items():
            if any(kw in lower for kw in keywords):
                matched.add(species)
        return matched

    def _matches_target_species(self, ob: dict, ctx: dict) -> bool:
        """Return False if the fetched item clearly belongs to a different ebolavirus species.

        True means "keep": either the item matches the target species, is generic,
        or we cannot confidently assign it to another species.
        """
        target = (ctx.get("species_name") or ctx.get("species_id") or "").strip().lower()
        all_keywords = {
            kw for keywords in self._EBOLAVIRUS_SPECIES.values() for kw in keywords
        }
        if not target or not any(kw in target for kw in all_keywords):
            return True  # not ebolavirus or unknown target, keep everything

        # Resolve target to a canonical species key
        target_species = None
        for species, keywords in self._EBOLAVIRUS_SPECIES.items():
            if any(kw in target for kw in keywords):
                target_species = species
                break
        if not target_species:
            return True

        # Prefer structured data fields, then title/description
        data = ob.get("data") or {}
        text_parts = [
            data.get("species"),
            data.get("pathogen"),
            data.get("strain"),
            ob.get("title"),
            ob.get("description"),
        ]
        text = " ".join(str(p) for p in text_parts if p)
        matched = self._detect_species(text)

        # EBOV abbreviation is ambiguous; do not let it alone drive a decision.
        if target_species == "zaire" and matched == set():
            # If target is Zaire and no other species is named, keep the item.
            return True

        if target_species in matched:
            return True
        if matched:
            # Item names one or more other species but not the target
            return False
        return True

    def _filter_trusted_sources(self, fetched_data: dict) -> dict:
        threshold = getattr(self.fetcher, "min_credibility", self.MIN_SOURCE_CREDIBILITY)
        filtered = {"_summary": fetched_data.get("_summary", {})}
        trusted_items = [
            item for item in fetched_data.get("items", [])
            if item.get("credibility", 0) >= threshold
        ]
        filtered["items"] = trusted_items
        log.info(f"Kept {len(trusted_items)} results above credibility threshold {threshold}")
        return filtered


def _find_db_query_results_for_sample(sample_id: str, db_results_arg: Optional[str]) -> Optional[Path]:
    """Return the most specific db_query_results.json available for a sample."""
    if not db_results_arg:
        return None
    db_path = Path(db_results_arg)
    if db_path.is_file():
        return db_path
    # If a directory was passed, look for per-sample result file
    candidate = db_path / sample_id / "db_query_results.json"
    if candidate.exists():
        return candidate
    return None


def _run_one(
    bio_path: Path,
    db_query_path: Optional[Path],
    output_path: Path,
    no_llm: bool,
    export_normalized: bool,
) -> None:
    """Run the epi query engine for a single bio_output.json file."""
    with open(bio_path) as f:
        bio_output = json.load(f)

    sample_id = (
        bio_output.get("sample", {}).get("sample_id")
        or (bio_output.get("stage9_normalised_output") or {}).get("sample_id")
        or bio_path.parent.name
    )

    local_db_results = None
    db_full: dict[str, Any] = {}
    if db_query_path:
        if db_query_path.exists():
            with open(db_query_path) as f:
                db_full = json.load(f)
            local_db_results = db_full.get("local_db_results") or db_full
        else:
            print(f"DB query results not found: {db_query_path}", file=sys.stderr)

    engine = EpiQueryEngine(use_llm=not no_llm)
    result = engine.query(
        bioinformatics_output=bio_output,
        local_db_results=local_db_results,
    )

    # Preserve the deterministic DB query layers so the intelligence pipeline
    # can match phenotypes, retrieve literature references, and enrich outputs
    # even when the LLM-driven epi step is the input.
    if local_db_results:
        result["local_db_results"] = local_db_results
    for layer_key in (
        "layer1_variant_lookup",
        "layer2_lineage_context",
        "layer3_geographic_temporal",
        "layer4_gene_function_context",
        "layer5_species_surveillance",
    ):
        if layer_key in db_full:
            result[layer_key] = db_full[layer_key]

    output_path.parent.mkdir(parents=True, exist_ok=True)
    print(f"[{sample_id}] Writing {output_path}...", flush=True)
    with open(output_path, "w") as f:
        json.dump(result, f, indent=2, default=str, ensure_ascii=False)
    print(f"[{sample_id}] Wrote {output_path}")

    if export_normalized:
        from intelligence_engine.data_engine.analytics.storage import write_parquet_tables

        print(f"[{sample_id}] Exporting normalized analytical tables (Parquet)...", flush=True)
        dataset = NormalizedEpiObject(**result)
        written = write_parquet_tables(dataset, str(output_path.parent))
        for table_name, path in written.items():
            value = getattr(dataset, table_name)
            count = len(value) if isinstance(value, list) else (1 if value is not None else 0)
            print(f"  {table_name}: {count} row(s) -> {path}")


def main():
    """CLI entry point for the epidemiological query engine."""
    parser = argparse.ArgumentParser(
        description="Fetch epidemiological context for one or more bioinformatics outputs."
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--bio-output",
        default=None,
        help="Path to a single bio_output.json.",
    )
    group.add_argument(
        "--bioinformatics-dir",
        default=None,
        help="Directory containing per-sample folders, each with bio_output.json. "
             "Processes every sample found and writes per-sample epi_output.json.",
    )
    parser.add_argument(
        "--db-query-results",
        default=None,
        help="Path to a single db_query_results.json, or a directory of per-sample results "
             "(used with --bioinformatics-dir).",
    )
    parser.add_argument(
        "--output",
        default=None,
        help="Path to write a single epi_output.json (required with --bio-output). "
             "Ignored when --bioinformatics-dir is used.",
    )
    parser.add_argument(
        "--output-dir",
        default="output/data_query",
        help="Directory to write per-sample epi_output.json files when --bioinformatics-dir is used "
             "(default: output/data_query).",
    )
    parser.add_argument(
        "--no-llm",
        action="store_true",
        help="Run in free mode with deterministic extraction only (no LLM required).",
    )
    parser.add_argument(
        "--export-normalized",
        action="store_true",
        help="Additionally export a normalized analytical dataset (Parquet tables) "
             "alongside the standard epi_output.json. Does not modify epi_output.json.",
    )
    args = parser.parse_args()

    if args.bio_output:
        if not args.output:
            parser.error("--output is required when using --bio-output")
        _run_one(
            bio_path=Path(args.bio_output),
            db_query_path=_find_db_query_results_for_sample(
                Path(args.bio_output).parent.name, args.db_query_results
            ),
            output_path=Path(args.output),
            no_llm=args.no_llm,
            export_normalized=args.export_normalized,
        )
        return

    # Directory mode: discover all per-sample bio_output.json files
    bioinformatics_dir = Path(args.bioinformatics_dir)
    if not bioinformatics_dir.is_dir():
        print(f"Bioinformatics directory not found: {bioinformatics_dir}", file=sys.stderr)
        sys.exit(1)

    output_dir = Path(args.output_dir)
    processed = 0
    for sample_dir in sorted(bioinformatics_dir.iterdir()):
        if not sample_dir.is_dir():
            continue
        bio_path = sample_dir / "bio_output.json"
        if not bio_path.exists():
            log.warning("Skipping %s — no bio_output.json", sample_dir.name)
            continue

        sample_output_dir = output_dir / sample_dir.name
        output_path = sample_output_dir / "epi_output.json"
        db_query_path = _find_db_query_results_for_sample(sample_dir.name, args.db_query_results)

        try:
            _run_one(
                bio_path=bio_path,
                db_query_path=db_query_path,
                output_path=output_path,
                no_llm=args.no_llm,
                export_normalized=args.export_normalized,
            )
            processed += 1
        except Exception as exc:
            print(f"[{sample_dir.name}] Epi query failed: {exc}", file=sys.stderr)

    print(f"Processed {processed} sample(s); outputs in {output_dir}")


if __name__ == "__main__":
    main()
