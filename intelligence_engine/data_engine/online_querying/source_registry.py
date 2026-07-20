"""Curated source-driven API registry for epidemiological data.

This module replaces the generic DuckDuckGo API-discovery with a curated registry
of trusted epidemiological data providers plus a lightweight discovery path over
public API directories (APIs.guru, OpenAPI registries).  Each provider entry
describes the source, its endpoints, the fields it supplies, and provenance
metadata.  Records are normalised into the same ``item`` schema used by the
epi fetcher, with a structured ``data`` payload so the query engine can extract
numbers deterministically without fragile text parsing.
"""
from __future__ import annotations

import json
import logging
import re
import urllib.request
from datetime import datetime, timezone
from typing import Any, Optional

try:
    import pycountry

    _country_names = sorted(
        {c.name for c in pycountry.countries}
        | {getattr(c, "official_name", "") for c in pycountry.countries}
        | {"DRC", "Congo", "Democratic Republic of the Congo", "Democratic Republic of Congo", "USA", "United States", "United Kingdom", "South Korea", "North Korea", "Russia", "Czechia", "Türkiye", "Turkey", "Cote d'Ivoire", "Côte d'Ivoire"},
        key=len,
        reverse=True,
    )
except Exception:
    _country_names = []

log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Curated registry of trusted epidemiological data providers
# ---------------------------------------------------------------------------

SOURCE_REGISTRY: list[dict[str, Any]] = [
    {
        "name": "EbolaIntel",
        "domains": ["ebolaintel.com"],
        "pathogen_terms": ["ebola", "ebolavirus", "ebov", "sudv", "bdbv", "tafv"],
        "endpoints": {
            "outbreaks": "https://ebolaintel.com/api/outbreaks.json",
            "strains": "https://ebolaintel.com/api/strains.json",
        },
        "supported_fields": [
            "outbreaks", "cases", "deaths", "cfr", "country", "year", "location",
            "strain", "active", "notable", "primary_source", "species_profile",
            "vaccine_coverage", "therapeutic_coverage", "confirmed_outbreaks_count",
        ],
        "credibility": 9,
        "reporting_agency": "EbolaIntel (aggregates CDC, WHO, peer-reviewed sources)",
        "update_frequency": "daily",
        "provenance": "curated aggregation",
        "api_version": "v1",
    },
    {
        "name": "OutbreakTinder",
        "domains": ["outbreaktinder.pages.dev", "cdn.jsdelivr.net", "github.com"],
        "pathogen_terms": ["all"],
        "endpoints": {
            "outbreaks": "https://cdn.jsdelivr.net/gh/ByteWorthyLLC/outbreaktinder@main/data/outbreaks.json",
        },
        "supported_fields": [
            "outbreaks", "cases", "deaths", "cfr_range", "transmission",
            "r0_range", "geographic_origin", "summary", "tags", "citations",
            "pathogen", "pathogen_type", "deaths_estimated",
        ],
        "credibility": 7,
        "reporting_agency": "OutbreakTinder (public historical outbreak dataset)",
        "update_frequency": "periodic",
        "provenance": "aggregated historical",
        "api_version": "1",
    },
    {
        # Placeholder / discovery-only entries for authoritative sources.  These
        # will be skipped until concrete endpoints are configured, but they keep
        # the registry aligned with the requested source framework.
        "name": "WHO",
        "domains": ["who.int"],
        "pathogen_terms": ["all"],
        "endpoints": {},
        "supported_fields": [],
        "credibility": 10,
        "reporting_agency": "World Health Organization",
        "update_frequency": "event-driven",
        "provenance": "official",
        "api_version": None,
    },
    {
        "name": "CDC",
        "domains": ["cdc.gov", "data.cdc.gov"],
        "pathogen_terms": ["all"],
        "endpoints": {},
        "supported_fields": [],
        "credibility": 10,
        "reporting_agency": "Centers for Disease Control and Prevention",
        "update_frequency": "weekly",
        "provenance": "official",
        "api_version": None,
    },
    {
        "name": "Africa CDC",
        "domains": ["africacdc.org"],
        "pathogen_terms": ["all"],
        "endpoints": {},
        "supported_fields": [],
        "credibility": 9,
        "reporting_agency": "Africa Centres for Disease Control and Prevention",
        "update_frequency": "event-driven",
        "provenance": "official",
        "api_version": None,
    },
    {
        "name": "ECDC",
        "domains": ["ecdc.europa.eu"],
        "pathogen_terms": ["all"],
        "endpoints": {},
        "supported_fields": [],
        "credibility": 9,
        "reporting_agency": "European Centre for Disease Prevention and Control",
        "update_frequency": "weekly",
        "provenance": "official",
        "api_version": None,
    },
    {
        "name": "ProMED",
        "domains": ["promedmail.org"],
        "pathogen_terms": ["all"],
        "endpoints": {},
        "supported_fields": [],
        "credibility": 6,
        "reporting_agency": "ProMED-mail (ISID)",
        "update_frequency": "daily",
        "provenance": "surveillance network",
        "api_version": None,
    },
    {
        "name": "HealthMap",
        "domains": ["healthmap.org"],
        "pathogen_terms": ["all"],
        "endpoints": {},
        "supported_fields": [],
        "credibility": 6,
        "reporting_agency": "HealthMap (Boston Children's Hospital)",
        "update_frequency": "real-time",
        "provenance": "automated surveillance",
        "api_version": None,
    },
    {
        "name": "ReliefWeb",
        "domains": ["reliefweb.int"],
        "pathogen_terms": ["all"],
        "endpoints": {},
        "supported_fields": [],
        "credibility": 7,
        "reporting_agency": "UN OCHA ReliefWeb",
        "update_frequency": "daily",
        "provenance": "humanitarian reports",
        "api_version": None,
    },
    {
        "name": "Our World in Data",
        "domains": ["ourworldindata.org"],
        "pathogen_terms": ["all"],
        "endpoints": {},
        "supported_fields": [],
        "credibility": 8,
        "reporting_agency": "Our World in Data (Oxford)",
        "update_frequency": "daily",
        "provenance": "curated data",
        "api_version": None,
    },
    {
        "name": "NCBI",
        "domains": ["ncbi.nlm.nih.gov"],
        "pathogen_terms": ["all"],
        "endpoints": {},
        "supported_fields": [],
        "credibility": 8,
        "reporting_agency": "National Center for Biotechnology Information",
        "update_frequency": "daily",
        "provenance": "genomic surveillance",
        "api_version": None,
    },
]


# ---------------------------------------------------------------------------
# Utility helpers
# ---------------------------------------------------------------------------

def _fetch_json(url: str, timeout: int = 25) -> Optional[Any]:
    """Fetch a JSON URL and return parsed data."""
    try:
        req = urllib.request.Request(
            url,
            headers={
                "User-Agent": "PGIRL-EpiFetcher/1.0",
                "Accept": "application/json, */*",
            },
        )
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            raw = resp.read().decode("utf-8", errors="replace")
            return json.loads(raw)
    except Exception as e:
        log.debug(f"JSON fetch failed for {url}: {e}")
        return None


def _terms_match(terms: list[str], provider: dict[str, Any]) -> bool:
    """Return True if any search term matches a provider's pathogen terms."""
    if "all" in provider.get("pathogen_terms", []):
        return True
    provider_terms = [t.lower() for t in provider.get("pathogen_terms", [])]
    for term in terms:
        t = term.lower()
        for pt in provider_terms:
            if pt in t or t in pt:
                return True
    return False


def _extract_country(text: str) -> Optional[str]:
    """Return the first known country name found in *text*."""
    if not text or not _country_names:
        return None
    lower = text.lower()
    for name in _country_names:
        if not name:
            continue
        pattern = rf"(?<![a-zA-Z]){re.escape(name.lower())}(?![a-zA-Z])"
        if re.search(pattern, lower):
            # Normalise common aliases
            if name in ("Congo", "DRC"):
                return "Democratic Republic of the Congo"
            if name in ("Cote d'Ivoire",):
                return "Côte d'Ivoire"
            return name
    return None


def _format_content(data: dict[str, Any]) -> str:
    """Serialize a structured data payload into a readable excerpt for the LLM."""
    parts = []
    for key, value in sorted(data.items()):
        if value is None or value == "":
            continue
        if isinstance(value, (list, dict)):
            value = json.dumps(value, ensure_ascii=False)
        parts.append(f"{key.replace('_', ' ').title()}: {value}")
    return " | ".join(parts)


def _make_item(
    title: str,
    url: str,
    source_name: str,
    provider: dict[str, Any],
    data: dict[str, Any],
    date: Optional[str] = None,
) -> dict[str, Any]:
    """Normalise a fetched record into the fetcher item schema."""
    domain = re.sub(r"^www\.", "", urllib.parse.urlparse(url).netloc.lower())
    content = _format_content(data)
    return {
        "source": f"API ({source_name})",
        "source_key": source_name.lower().replace(" ", "_"),
        "domain": domain,
        "credibility": provider.get("credibility", 5),
        "title": title,
        "url": url,
        "date": date or "",
        "description": content[:500],
        "content_excerpt": content[:4000],
        "data": data,
        "source_metadata": {
            "reporting_agency": provider.get("reporting_agency", source_name),
            "update_frequency": provider.get("update_frequency", "unknown"),
            "provenance": provider.get("provenance", "unknown"),
            "api_version": provider.get("api_version"),
            "provider_name": provider.get("name", source_name),
            "reporting_date": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
        },
    }


# ---------------------------------------------------------------------------
# Provider-specific fetchers
# ---------------------------------------------------------------------------

def _fetch_ebolaintel(provider: dict[str, Any], terms: list[str]) -> list[dict[str, Any]]:
    """Fetch EbolaIntel outbreak and strain records."""
    items: list[dict[str, Any]] = []
    endpoints = provider.get("endpoints", {})

    outbreaks = _fetch_json(endpoints.get("outbreaks", "")) or {}
    for record in outbreaks.get("outbreaks", []):
        country = ""
        countries = record.get("countries") or [record.get("country")]
        if isinstance(countries, list) and countries:
            country = countries[0]
        elif isinstance(countries, str):
            country = countries

        # Normalise historical DRC labels
        if country and ("zaire" in country.lower() or country in ("DRC", "Congo, The Democratic Republic of the")):
            country = "Democratic Republic of the Congo"
        elif country == "Congo":
            # Preserve the distinction; EbolaIntel uses 'Republic of the Congo' explicitly elsewhere
            country = "Republic of the Congo"

        year_start = record.get("yearStart")
        year_end = record.get("yearEnd")
        date = str(year_start) if year_start else ""
        if year_end and year_end != year_start:
            date = f"{year_start}-{year_end}"

        data = {
            "item_type": "outbreak",
            "cases": record.get("cases"),
            "deaths": record.get("deaths"),
            "cfr": record.get("cfrPct"),
            "country": country,
            "year_start": year_start,
            "year_end": year_end,
            "location": record.get("location"),
            "strain": record.get("strain"),
            "active": record.get("active"),
            "notable": record.get("notable"),
            "primary_source": record.get("primarySource"),
        }

        url = record.get("primarySource") or endpoints.get("outbreaks", "")
        title = f"Outbreak {country} {year_start or ''} {record.get('strain') or ''}".strip()
        items.append(_make_item(title, url, provider["name"], provider, data, date=date))

    strains = _fetch_json(endpoints.get("strains", "")) or {}
    for record in strains.get("strains", []):
        name = record.get("scientificName") or record.get("commonName") or record.get("abbreviation")
        cfr_range = record.get("cfrRangePct") or []
        totals = record.get("totals") or {}
        data = {
            "item_type": "strain_profile",
            "species": record.get("scientificName"),
            "abbreviation": record.get("abbreviation"),
            "common_name": record.get("commonName"),
            "discovered_year": record.get("discoveredYear"),
            "discovered_location": record.get("discoveredLocation"),
            "pathogenic_in_humans": record.get("pathogenicInHumans"),
            "cfr_range": "-".join(str(v) for v in cfr_range) if cfr_range else None,
            "confirmed_outbreaks_count": record.get("confirmedHumanOutbreaks"),
            "vaccine_coverage": record.get("approvedVaccineCoverage"),
            "therapeutic_coverage": record.get("approvedTherapeuticCoverage"),
            "total_cases": totals.get("cases"),
            "total_deaths": totals.get("deaths"),
            "total_outbreaks": totals.get("count"),
        }
        title = f"Species profile: {name}"
        url = endpoints.get("strains", "")
        date = str(record.get("discoveredYear") or "")
        items.append(_make_item(title, url, provider["name"], provider, data, date=date))

    return items


def _fetch_outbreaktinder(provider: dict[str, Any], terms: list[str]) -> list[dict[str, Any]]:
    """Fetch OutbreakTinder records and keep those matching the pathogen terms."""
    items: list[dict[str, Any]] = []
    endpoints = provider.get("endpoints", {})
    data = _fetch_json(endpoints.get("outbreaks", ""))
    if not isinstance(data, list):
        return items

    terms_lower = [t.lower() for t in terms]
    for record in data:
        pathogen = str(record.get("pathogen", "")).lower()
        name = str(record.get("name", "")).lower()
        tags = [str(t).lower() for t in record.get("tags", [])]
        if not any(
            t in pathogen or t in name or any(t in tag for tag in tags)
            for t in terms_lower
        ):
            continue

        cfr_range = record.get("cfr_range") or {}
        r0_range = record.get("r0_range") or {}
        deaths_est = record.get("deaths_estimated") or {}
        years = record.get("years", "")
        start_year = None
        end_year = None
        if years and isinstance(years, str):
            m = re.findall(r"\d{4}", years)
            if m:
                start_year = m[0]
                end_year = m[-1] if len(m) > 1 else start_year

        # Try to parse explicit cases/deaths from the estimate note, e.g.
        # "318 reported cases, 280 deaths."
        cases = None
        deaths = None
        note = str(deaths_est.get("note", ""))
        if note:
            m_cases = re.search(r"([\d, ]+)\s+reported\s+cases", note, re.IGNORECASE)
            if m_cases:
                try:
                    cases = int(m_cases.group(1).replace(",", "").replace(" ", ""))
                except (ValueError, AttributeError):
                    pass
            m_deaths = re.search(r"([\d, ]+)\s+deaths", note, re.IGNORECASE)
            if m_deaths:
                try:
                    deaths = int(m_deaths.group(1).replace(",", "").replace(" ", ""))
                except (ValueError, AttributeError):
                    pass

        # Derive a clean country from the geographic origin / name
        origin = str(record.get("geographic_origin", ""))
        country = _extract_country(origin) or _extract_country(record.get("name", ""))

        data_payload = {
            "item_type": "outbreak",
            "name": record.get("name"),
            "years": years,
            "year_start": start_year,
            "year_end": end_year,
            "pathogen": record.get("pathogen"),
            "pathogen_type": record.get("pathogen_type"),
            "country": country,
            "transmission": record.get("transmission"),
            "r0_range": f"{r0_range.get('min')}-{r0_range.get('max')}" if r0_range else None,
            "cfr_range": f"{cfr_range.get('min')}-{cfr_range.get('max')}" if cfr_range else None,
            "cases": cases,
            "deaths": deaths,
            "deaths_estimated_low": deaths_est.get("low"),
            "deaths_estimated_high": deaths_est.get("high"),
            "deaths_estimate_note": deaths_est.get("note"),
            "summary": record.get("summary"),
            "geographic_origin": record.get("geographic_origin"),
            "tags": record.get("tags"),
            "first_reported_month_day": record.get("first_reported_month_day"),
        }

        # Build citations into a provenance string
        citations = record.get("citations") or []
        primary_source = ""
        if citations:
            primary_source = citations[0].get("source_url", "")
        data_payload["primary_sources"] = [c.get("source_url") for c in citations if c.get("source_url")]

        title = f"Historical outbreak: {record.get('name', '')}"
        url = primary_source or endpoints.get("outbreaks", "")
        items.append(_make_item(title, url, provider["name"], provider, data_payload, date=start_year))

    return items


# ---------------------------------------------------------------------------
# Public API discovery (APIs.guru and OpenAPI registries)
# ---------------------------------------------------------------------------

def _discover_apis_guru(terms: list[str]) -> list[dict[str, Any]]:
    """Query APIs.guru registry for APIs matching the pathogen terms."""
    candidates: list[dict[str, Any]] = []
    try:
        listing = _fetch_json("https://api.apis.guru/v2/list.json", timeout=20)
        if not isinstance(listing, dict):
            return candidates
        for api_key, api_info in listing.items():
            info = (api_info.get("versions") or {}).get(api_info.get("preferred"), {}).get("info", {})
            text = " ".join([
                str(info.get("title", "")),
                str(info.get("description", "")),
                str(api_key),
            ]).lower()
            if any(term.lower() in text for term in terms):
                candidates.append({
                    "name": info.get("title", api_key),
                    "domain": api_info.get("preferredScheme", "") + "://" + api_key,
                    "endpoints": {api_info.get("preferred", ""): api_info.get("swaggerUrl") or api_info.get("openapiUrl")},
                    "credibility": 5,
                    "reporting_agency": "APIs.guru discovery",
                    "update_frequency": "unknown",
                    "provenance": "discovered",
                    "api_version": api_info.get("preferred"),
                })
    except Exception as e:
        log.debug(f"APIs.guru discovery failed: {e}")
    return candidates


# ---------------------------------------------------------------------------
# Source registry class
# ---------------------------------------------------------------------------

class SourceRegistry:
    """Curated registry + discovery of epidemiological data sources."""

    def __init__(self, registry: Optional[list[dict[str, Any]]] = None):
        self.registry = registry or SOURCE_REGISTRY

    def fetch_items(self, species_name: str, search_terms: list[str]) -> list[dict[str, Any]]:
        """Fetch all curated and discovered API items for the given species."""
        items: list[dict[str, Any]] = []
        for provider in self.registry:
            if not _terms_match(search_terms, provider):
                continue
            name = provider.get("name", "").lower()
            if name == "ebolaintel":
                items.extend(_fetch_ebolaintel(provider, search_terms))
            elif name == "outbreaktinder":
                items.extend(_fetch_outbreaktinder(provider, search_terms))
            # Placeholder providers with no endpoints are skipped until configured.

        # Discovery path: also look for new APIs in public registries.  These are
        # returned as low-credibility source cards so the pipeline can decide
        # whether to try them.
        for candidate in _discover_apis_guru(search_terms):
            if candidate.get("endpoints"):
                items.extend(self._try_candidate(candidate, search_terms))
        return items

    def _try_candidate(self, candidate: dict[str, Any], terms: list[str]) -> list[dict[str, Any]]:
        """Attempt to fetch a discovered API candidate if it looks like a JSON endpoint."""
        items: list[dict[str, Any]] = []
        for endpoint in candidate.get("endpoints", {}).values():
            if not endpoint or not endpoint.lower().endswith(".json"):
                continue
            data = _fetch_json(endpoint, timeout=15)
            if not data:
                continue
            records = data if isinstance(data, list) else data.get("outbreaks") or data.get("results") or data.get("data") or []
            if not isinstance(records, list):
                continue
            for record in records[:50]:
                if not isinstance(record, dict):
                    continue
                title = record.get("name") or record.get("title") or f"Discovered API record for {terms[0]}"
                items.append(_make_item(
                    title,
                    endpoint,
                    candidate["name"],
                    candidate,
                    record,
                    date=str(record.get("year") or record.get("date") or ""),
                ))
        return items
