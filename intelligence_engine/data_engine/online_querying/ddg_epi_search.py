"""ddg_epi_search.py — DuckDuckGo-backed epidemiological data fetcher.

Searches DuckDuckGo broadly for epidemiological reports, ranks the returned
pages by dynamic source-credibility heuristics (WHO, CDC, UN, government
domains, Ministries of Health, etc.), fetches the pages, and passes clean
text to the LLM for structured extraction.

No hard-coded per-site HTML parsers; new pathogens only need an entry in
SPECIES_SEARCH_TERMS and new official domains are picked up automatically by
the credibility scorer.

Recommended dependency (auto-detected, falls back to a lightweight HTML
scraper if unavailable):
    pip install ddgs
"""

import html
import json
import logging
import re
import urllib.parse
import urllib.request
from datetime import datetime, timedelta, timezone
from itertools import groupby
from typing import Any, Optional

from .source_registry import SourceRegistry

try:
    import pycountry
    names = {c.name for c in pycountry.countries}
    official = {getattr(c, "official_name", "") for c in pycountry.countries}
    aliases = {"DRC", "Congo", "Democratic Republic of the Congo", "Democratic Republic of Congo", "USA", "United States", "United Kingdom", "South Korea", "North Korea", "Russia", "Venezuela", "Cote d'Ivoire", "Côte d'Ivoire", "Czechia", "Türkiye", "Turkey"}
    _COUNTRY_LIST = sorted(
        ((names | official | aliases) - {""}),
        key=len, reverse=True,
    )
except Exception:
    _COUNTRY_LIST = [
        "Afghanistan", "Albania", "Algeria", "Angola", "Argentina", "Armenia", "Australia", "Austria",
        "Bangladesh", "Belarus", "Belgium", "Benin", "Bolivia", "Bosnia and Herzegovina", "Botswana",
        "Brazil", "Bulgaria", "Burkina Faso", "Burundi", "Cambodia", "Cameroon", "Canada", "Central African Republic",
        "Chad", "Chile", "China", "Colombia", "Congo", "Costa Rica", "Croatia", "Cuba", "Cyprus", "Czech Republic",
        "Democratic Republic of Congo", "Denmark", "Djibouti", "Dominican Republic", "DRC", "Ecuador", "Egypt",
        "El Salvador", "Equatorial Guinea", "Eritrea", "Estonia", "Eswatini", "Ethiopia", "Fiji", "Finland",
        "France", "Gabon", "Gambia", "Georgia", "Germany", "Ghana", "Greece", "Guatemala", "Guinea", "Guinea-Bissau",
        "Guyana", "Haiti", "Honduras", "Hungary", "Iceland", "India", "Indonesia", "Iran", "Iraq", "Ireland", "Israel",
        "Italy", "Ivory Coast", "Jamaica", "Japan", "Jordan", "Kazakhstan", "Kenya", "Kuwait", "Kyrgyzstan", "Laos",
        "Latvia", "Lebanon", "Lesotho", "Liberia", "Libya", "Lithuania", "Luxembourg", "Madagascar", "Malawi", "Malaysia",
        "Maldives", "Mali", "Malta", "Mauritania", "Mauritius", "Mexico", "Moldova", "Mongolia", "Montenegro", "Morocco",
        "Mozambique", "Myanmar", "Namibia", "Nepal", "Netherlands", "New Zealand", "Nicaragua", "Niger", "Nigeria",
        "North Korea", "North Macedonia", "Norway", "Oman", "Pakistan", "Panama", "Papua New Guinea", "Paraguay", "Peru",
        "Philippines", "Poland", "Portugal", "Qatar", "Romania", "Russia", "Rwanda", "Saudi Arabia", "Senegal",
        "Serbia", "Sierra Leone", "Singapore", "Slovakia", "Slovenia", "Somalia", "South Africa", "South Korea",
        "South Sudan", "Spain", "Sri Lanka", "Sudan", "Suriname", "Sweden", "Switzerland", "Syria", "Taiwan",
        "Tajikistan", "Tanzania", "Thailand", "Timor-Leste", "Togo", "Trinidad and Tobago", "Tunisia", "Turkey",
        "Turkmenistan", "Uganda", "Ukraine", "United Arab Emirates", "United Kingdom", "United States", "Uruguay",
        "Uzbekistan", "Venezuela", "Vietnam", "Yemen", "Zambia", "Zimbabwe",
    ]

log = logging.getLogger(__name__)

HAVE_DDG = False
try:
    from ddgs import DDGS
    HAVE_DDG = True
except ImportError:
    pass

DDG_LITE_URL = "https://html.duckduckgo.com/html/"

SPECIES_SEARCH_TERMS = {
    "Zaire ebolavirus": ["ebola", "EVD", "Zaire ebolavirus", "EBOV"],
    "Sudan ebolavirus": ["sudan ebolavirus", "SUDV", "ebola sudan"],
    "Bundibugyo ebolavirus": ["bundibugyo", "BDBV"],
    "Marburg marburgvirus": ["marburg", "MARV"],
    "Dengue virus": ["dengue", "DENV"],
    "Influenza A virus": ["influenza", "flu", "H5N1", "H1N1", "H3N2"],
    "Rift Valley fever virus": ["rift valley fever", "RVF"],
    "Monkeypox virus": ["mpox", "monkeypox", "MPXV"],
}

_TITLE_REQUIRED = {
    "Zaire ebolavirus": ["ebola", "EVD", "Zaire ebolavirus", "EBOV"],
    "Sudan ebolavirus": ["ebola", "SUDV", "Sudan ebolavirus", "Sudan"],
    "Bundibugyo ebolavirus": ["ebola", "BDBV", "Bundibugyo"],
    "Marburg marburgvirus": ["marburg", "MARV"],
    "Dengue virus": ["dengue", "DENV"],
    "Influenza A virus": ["influenza", "flu", "H5N1", "H1N1", "H3N2"],
    "Rift Valley fever virus": ["rift valley fever", "RVF"],
    "Monkeypox virus": ["mpox", "monkeypox", "MPXV"],
}

# (domain substring, credibility bonus, display name, source key)
OFFICIAL_DOMAIN_PATTERNS = [
    ("who.int", 10, "WHO", "who"),
    ("cdc.gov", 9, "CDC", "cdc"),
    ("un.org", 9, "UN", "un"),
    ("unicef.org", 9, "UNICEF", "unicef"),
    ("africacdc.org", 8, "Africa CDC", "africa_cdc"),
    ("ecdc.europa.eu", 8, "ECDC", "ecdc"),
    ("reliefweb.int", 7, "ReliefWeb", "reliefweb"),
    ("gov.uk", 8, "GOV.UK", "gov_uk"),
    ("nih.gov", 7, "NIH/PMC", "nih_pmc"),
]

OFFICIAL_TLD_SUFFIXES = [
    ".gov", ".gov.uk", ".gov.au", ".gov.ng", ".gov.sl", ".gov.lr", ".gov.gh",
    ".gouv.cd", ".gouv.ci", ".gouv.sn", ".gouv.ml", ".gouv.ga", ".gouv.gn",
    ".go.ug", ".go.ke", ".go.tz", ".go.zm", ".go.zw", ".go.rw", ".go.za",
    ".gob.mx", ".gob.ar", ".gob.es", ".govt.nz", ".mil",
]

UNTRUSTED_DOMAINS = [
    "wikipedia.org", "youtube.com", "facebook.com", "twitter.com", "x.com",
    "reddit.com", "tiktok.com", "instagram.com", "pinterest.com",
]



class DuckDuckGoSearchClient:
    def __init__(self, timeout: int = 20):
        self.timeout = timeout

    def search(self, query: str, max_results: int = 10) -> list:
        if HAVE_DDG:
            return self._search_with_library(query, max_results)
        return self._search_html(query, max_results)

    def _search_with_library(self, query: str, max_results: int) -> list:
        with DDGS() as ddgs:
            results = ddgs.text(query, max_results=max_results)
            return [
                {"title": r.get("title", ""), "href": r.get("href", ""), "body": r.get("body", "")}
                for r in results
            ]

    def _search_html(self, query: str, max_results: int) -> list:
        encoded = urllib.parse.quote_plus(query)
        url = f"{DDG_LITE_URL}?q={encoded}&kl=us-en"
        headers = {
            "User-Agent": (
                "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 "
                "(KHTML, like Gecko) Chrome/115.0.0.0 Safari/537.36"
            ),
            "Accept": "text/html,application/xhtml+xml",
        }
        req = urllib.request.Request(url, headers=headers, method="GET")
        try:
            with urllib.request.urlopen(req, timeout=self.timeout) as resp:
                raw = resp.read().decode("utf-8", errors="replace")
        except Exception as e:
            raise RuntimeError(f"DuckDuckGo HTML search failed: {e}")

        results = []
        blocks = re.findall(r'<div class="result__body">(.*?)</div>\s*</div>', raw, re.DOTALL)
        for block in blocks[:max_results]:
            title_match = re.search(
                r'<a class="result__a" href="([^"]+)"[^>]*>(.*?)</a>',
                block,
                re.DOTALL,
            )
            if not title_match:
                continue
            href = self._decode_ddg_redirect(title_match.group(1))
            title = re.sub(r"<[^>]+", "", title_match.group(2)).strip()
            snippet_match = re.search(
                r'<div class="result__snippet">(.*?)</div>',
                block,
                re.DOTALL,
            )
            snippet = ""
            if snippet_match:
                snippet = re.sub(r"<[^>]+", "", snippet_match.group(1)).strip()
            if href:
                results.append({
                    "title": html.unescape(title),
                    "href": href,
                    "body": html.unescape(snippet),
                })
        return results

    @staticmethod
    def _decode_ddg_redirect(href: str) -> str:
        if href.startswith("http") and "duckduckgo.com" not in href:
            return href
        parsed = urllib.parse.urlparse(href)
        qs = urllib.parse.parse_qs(parsed.query)
        if "uddg" in qs:
            return urllib.parse.unquote(qs["uddg"][0])
        return href


class DuckDuckGoEpiFetcher:
    def __init__(
        self,
        timeout: int = 20,
        min_credibility: int = 5,
        max_results_per_query: int = 20,
        max_total_results: int = 30,
    ):
        self.timeout = timeout
        self.min_credibility = min_credibility
        self.max_results_per_query = max_results_per_query
        self.max_total_results = max_total_results
        self.ddg = DuckDuckGoSearchClient(timeout=timeout)

    def fetch_all(
        self,
        species_name: str,
        lineage: Optional[str] = None,
        country: Optional[str] = None,
        collection_date: Optional[str] = None,
    ) -> dict:
        self.species_name = species_name
        search_terms = SPECIES_SEARCH_TERMS.get(species_name, [species_name])

        # Fetch from the curated source-driven API registry first. This is more
        # reliable and maintainable than DuckDuckGo scraping/API discovery.
        registry = SourceRegistry()
        api_items = registry.fetch_items(species_name, search_terms)

        queries = self._build_queries(species_name, lineage, country)

        items: list[dict[str, Any]] = list(api_items)
        seen_urls: set[str] = {it["url"] for it in items}

        # If the API already returned solid outbreak records, skip the expensive
        # DDG page scraping for this species and only use it for the odd page
        # that is already returned as a primary source in the API records.
        if not api_items:
            for query in queries:
                try:
                    ddg_results = self.ddg.search(query, max_results=self.max_results_per_query)
                except Exception as e:
                    log.warning(f"DuckDuckGo search failed for query '{query}': {e}")
                    continue

                for sr in ddg_results:
                    url = sr.get("href", "")
                    if not url or url in seen_urls:
                        continue
                    title = sr.get("title", "")
                    snippet = sr.get("body", "")
                    combined = f"{title} {snippet}"
                    if not self._matches_terms(combined, search_terms):
                        continue
                    if not self._is_relevant_species(title, combined, species_name):
                        continue

                    credibility, source_name, source_key = self._classify_result(url, title, snippet)
                    if credibility < self.min_credibility:
                        continue

                    seen_urls.add(url)
                    page_text = self._fetch_page_text(url)
                    content = page_text[:4000] if page_text else snippet
                    items.append({
                        "source": source_name,
                        "source_key": source_key,
                        "domain": urllib.parse.urlparse(url).netloc,
                        "credibility": credibility,
                        "title": self._clean_text(title, 300),
                        "url": url,
                        "date": self._extract_date_from_text(page_text or snippet),
                        "description": self._clean_text(snippet, 500),
                        "content_excerpt": self._clean_text(content, 2000),
                    })

        items.sort(key=lambda x: (-x.get("credibility", 0), x.get("source", "")))
        items = items[: self.max_total_results + len(api_items)]

        results = {
            "items": items,
            "_summary": {
                "total_reports_fetched": len(items),
                "api_records_fetched": len(api_items),
                "queries": queries,
                "fetch_timestamp": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
                "search_terms": search_terms,
                "species": species_name,
                "lineage": lineage,
                "country": country,
            },
        }
        log.info(f"Fetched {len(items)} epidemiological reports (API records: {len(api_items)})")
        return results

    def _build_queries(
        self, species_name: str, lineage: Optional[str], country: Optional[str]
    ) -> list[str]:
        terms = SPECIES_SEARCH_TERMS.get(species_name, [species_name])
        term_str = " ".join(f'"{t}"' if " " in t else t for t in terms)
        base = f"{term_str} {country or ''}".strip()
        queries = [
            # Global epidemiology for the pathogen
            f"{term_str} outbreak cases deaths".strip(),
            f"{term_str} epidemiology statistics".strip(),
            f"{term_str} situation report".strip(),
            f"{term_str} outbreak history cases deaths".strip(),
            # Official health-agency historical outbreak pages
            f"{term_str} \"History of Ebola Outbreaks\" site:cdc.gov".strip(),
            f"{term_str} \"Cases and outbreaks by year\" site:cdc.gov".strip(),
            f"{term_str} \"Disease Outbreak News\" site:who.int".strip(),
            f"{term_str} outbreak site:who.int".strip(),
            # Country-specific context when known
            f"{base} outbreak cases deaths".strip(),
            # Molecular epidemiology: strain/genotype/clade burden
            f"{term_str} strain genotype clade cases".strip(),
            # Seroepidemiology
            f"{term_str} seroprevalence antibody survey".strip(),
            # One Health / zoonotic spillover
            f"{term_str} reservoir host animal spillover".strip(),
            # Environmental / ecological epidemiology
            f"{term_str} environmental risk factors seasonality".strip(),
        ]
        return queries

    def _classify_result(self, url: str, title: str, snippet: str) -> tuple[int, str, str]:
        parsed = urllib.parse.urlparse(url)
        domain = parsed.netloc.lower()
        if domain.startswith("www."):
            domain = domain[4:]
        path = parsed.path.lower()
        text = f"{domain} {path} {title} {snippet}".lower()
        score = 0
        source_name: Optional[str] = None
        source_key = "unknown"

        for pat, bonus, name, key in OFFICIAL_DOMAIN_PATTERNS:
            if pat in domain:
                score += bonus
                source_name = name
                source_key = key
                return score, source_name, source_key

        if any(domain.endswith(suffix) for suffix in OFFICIAL_TLD_SUFFIXES):
            score += 5
            source_name = "Government"
            source_key = "government"
            health_signals = ["health", "moh", "sante", "cdc", "medic", "epidem", "doh"]
            if any(s in domain or s in path for s in health_signals):
                score += 3
                source_name = "Government Health Agency"
                source_key = "government_health"

        if source_name is None:
            if "ministry of health" in text or "ministère de la santé" in text:
                score += 2
                source_name = "Ministry of Health"
                source_key = "moh"
            elif "health" in domain or "sante" in domain:
                score += 1
                source_name = "Health-related source"
                source_key = "health_source"

        if "who" in text:
            score += 1
        if "cdc" in text:
            score += 1
        if "unicef" in text or "un " in text or "united nations" in text:
            score += 1

        for bad in UNTRUSTED_DOMAINS:
            if bad in domain:
                score -= 10
                source_name = bad.split(".")[0].capitalize()
                source_key = bad.split(".")[0]

        if source_name is None:
            source_name = domain.split(".")[0].upper() or "Web source"
        return score, source_name, source_key

    def format_for_llm(self, fetched_data: dict, max_chars: int = 8000) -> str:
        lines = ["## FETCHED EPIDEMIOLOGICAL DATA (FROM TRUSTED SOURCES)"]
        lines.append("The data below was fetched from official public health sources.")
        lines.append("Use ONLY this data to answer the epidemiological questions.")
        lines.append("If the fetched data does not contain an answer, return null.")
        lines.append("")

        summary = fetched_data.get("_summary", {})
        lines.append(f"Fetched at: {summary.get('fetch_timestamp', 'unknown')}")
        lines.append(f"Search terms: {', '.join(summary.get('search_terms', []))}")
        lines.append("")

        items = sorted(
            fetched_data.get("items", []),
            key=lambda x: (-x.get("credibility", 0), x.get("source", "")),
        )
        for source_name, group in groupby(items, key=lambda x: x.get("source", "Unknown")):
            group_items = list(group)
            lines.append(f"### {source_name} ({len(group_items)} reports)")
            for i, item in enumerate(group_items[:10], 1):
                lines.append(f"  [{i}] {item.get('title', 'untitled')}")
                if item.get("date"):
                    lines.append(f"      Date: {item['date']}")
                if item.get("description"):
                    lines.append(f"      Details: {item['description'][:300]}")
                if item.get("content_excerpt"):
                    lines.append(f"      Content: {item['content_excerpt'][:2000]}")
                if item.get("url"):
                    lines.append(f"      Source URL: {item['url']}")
                lines.append("")

            current = len("\n".join(lines))
            if current > max_chars * 0.8:
                lines.append("(Remaining sources truncated to fit context window)")
                break

        result = "\n".join(lines)
        if len(result) > max_chars:
            result = result[:max_chars] + "\n... (truncated)"
        return result

    def get_reference_urls(self, fetched_data: dict) -> list:
        refs = []
        seen: set[str] = set()
        for item in fetched_data.get("items", []):
            url = item.get("url")
            if not url or url in seen:
                continue
            seen.add(url)
            refs.append({
                "source": item.get("source", "Web source"),
                "title": item.get("title", "")[:80],
                "url": url,
            })
        return refs

    def structure_fetched_data(
        self,
        fetched_data: dict,
        species_name: str = "",
        lineage: str = "",
        country: str = "",
    ) -> dict:
        outbreaks = []
        countries_affected = set()
        source_counts: dict[str, int] = {}

        for item in fetched_data.get("items", []):
            source_key = item.get("source_key", "unknown")
            source_counts[source_key] = source_counts.get(source_key, 0) + 1
            outbreak = {
                "source": item.get("source", source_key),
                "source_key": source_key,
                "credibility": item.get("credibility", 0),
                "title": self._clean_text(item.get("title", ""), 300),
                "date": item.get("date", ""),
                "url": item.get("url", ""),
            }
            if item.get("description"):
                outbreak["description"] = self._clean_text(item["description"], 800)
            if item.get("content_excerpt"):
                outbreak["content"] = self._clean_text(item["content_excerpt"], 1200)

            text = f"{outbreak.get('title','')} {outbreak.get('description','')} {outbreak.get('content','')}"
            found_countries = self._extract_countries(text)
            if found_countries:
                outbreak["countries"] = found_countries
                countries_affected.update(found_countries)
            outbreaks.append(outbreak)

        outbreaks = self._deduplicate_outbreaks(outbreaks)
        outbreaks.sort(key=lambda x: self._parse_date(x.get("date", "")), reverse=True)

        cutoff = datetime.now(timezone.utc) - timedelta(days=730)
        recent = [
            ob for ob in outbreaks
            if self._parse_date(ob.get("date", "")) > cutoff
        ]

        return {
            "pathogen": species_name,
            "lineage": lineage,
            "query_country": country,
            "total_reports": len(outbreaks),
            "recent_outbreaks_2yr": len(recent),
            "countries_affected": sorted(countries_affected),
            "source_counts": source_counts,
            "latest_report_date": outbreaks[0]["date"] if outbreaks else None,
            "outbreaks": outbreaks,
            "recent_outbreaks": recent,
            "fetch_timestamp": fetched_data.get("_summary", {}).get("fetch_timestamp", ""),
        }

    def _fetch_page_text(self, url: str) -> Optional[str]:
        raw = self._make_request(url, timeout=self.timeout)
        if not raw:
            return None
        text = re.sub(r"<script[^>]*>.*?</script>", "", raw, flags=re.DOTALL)
        text = re.sub(r"<style[^>]*>.*?</style>", "", text, flags=re.DOTALL)
        # Preserve table row/cell structure so the LLM can read case/death numbers
        text = re.sub(r"</tr\s*>", "\n", text, flags=re.IGNORECASE)
        text = re.sub(r"<t[dh]\b[^>]*>", " | ", text, flags=re.IGNORECASE)
        text = re.sub(r"</t[dh]\s*>", " ", text, flags=re.IGNORECASE)
        text = re.sub(r"<[^>]+>", " ", text)
        text = html.unescape(text)
        text = re.sub(r"\s+\|\s+", " | ", text)
        text = re.sub(r"\n\s*\n+", "\n", text)
        text = re.sub(r"[ \t]+", " ", text).strip()
        return text[:8000]

    @staticmethod
    def _extract_date_from_text(text: Optional[str]) -> str:
        if not text:
            return ""
        text = str(text)
        m = re.search(r"\b(\d{4}-\d{2}-\d{2})\b", text)
        if m:
            return m.group(1)
        m = re.search(
            r"\b((?:January|February|March|April|May|June|July|August|September|October|November|December)\s+\d{1,2},?\s+\d{4})\b",
            text,
            re.IGNORECASE,
        )
        if m:
            return m.group(1)
        m = re.search(
            r"\b(\d{1,2}\s+(?:January|February|March|April|May|June|July|August|September|October|November|December)\s+\d{4})\b",
            text,
            re.IGNORECASE,
        )
        if m:
            return m.group(1)
        return ""

    def _make_request(
        self, url: str, headers: Optional[dict] = None, timeout: Optional[int] = None
    ) -> Optional[str]:
        default_headers = {
            "User-Agent": "PGIRL-EpiFetcher/1.0 (Genomic Intelligence System)",
            "Accept": "application/json, text/html, application/xml, */*",
        }
        if headers:
            default_headers.update(headers)
        try:
            req = urllib.request.Request(url, headers=default_headers, method="GET")
            with urllib.request.urlopen(req, timeout=timeout or self.timeout) as resp:
                return resp.read().decode("utf-8", errors="replace")
        except Exception as e:
            log.warning(f"Direct request failed for {url}: {e}")
            return self._fetch_via_jina(url, timeout=timeout or self.timeout)

    def _fetch_via_jina(self, url: str, timeout: int = 20) -> Optional[str]:
        """Use r.jina.ai to extract markdown text when a site blocks direct access."""
        try:
            encoded = urllib.parse.quote(url, safe="")
            proxy_url = f"https://r.jina.ai/http://{encoded}"
            req = urllib.request.Request(proxy_url, headers={
                "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36",
                "Accept": "text/plain",
            })
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read().decode("utf-8", errors="replace")
        except Exception as e:
            log.warning(f"Proxy request failed for {url}: {e}")
            return None

    @staticmethod
    def _matches_terms(text: str, terms: list) -> bool:
        text_lower = text.lower()
        return any(term.lower() in text_lower for term in terms)

    def _is_relevant_species(self, title: str, text: str, species_name: str) -> bool:
        title_lower = (title or "").lower()
        species = species_name or ""
        required = _TITLE_REQUIRED.get(species, [species.lower()])
        return any(req.lower() in title_lower for req in required)

    @staticmethod
    def _clean_text(text: Any, max_chars: int = 3000) -> str:
        if not text:
            return ""
        text = str(text)
        text = re.sub(r"<[^>]+", " ", text)
        text = html.unescape(text)
        text = re.sub(r"\s+", " ", text).strip()
        return text[:max_chars]

    def _extract_countries(self, text: str) -> list:
        found = []
        text_lower = text.lower()
        for country in _COUNTRY_LIST:
            if country.lower() in text_lower:
                norm = "Democratic Republic of Congo" if country in ("DRC", "Congo") else country
                if norm not in found:
                    found.append(norm)
        return found

    def _parse_date(
        self, date_str: str, default: datetime = datetime.min.replace(tzinfo=timezone.utc)
    ) -> datetime:
        if not date_str:
            return default
        date_str = date_str.strip()
        formats = [
            "%Y-%m-%dT%H:%M:%SZ",
            "%Y-%m-%dT%H:%M:%S.%fZ",
            "%Y-%m-%dT%H:%M:%S",
            "%Y-%m-%d",
            "%a, %d %b %Y %H:%M:%S %Z",
            "%a, %d %b %Y",
            "%d %b %Y",
        ]
        for fmt in formats:
            try:
                return datetime.strptime(date_str[:25], fmt).replace(tzinfo=timezone.utc)
            except (ValueError, TypeError):
                continue
        return default

    def _deduplicate_outbreaks(self, outbreaks: list) -> list:
        best = {}
        for ob in outbreaks:
            key = self._outbreak_key(ob)
            if key not in best or self._outbreak_score(ob) > self._outbreak_score(best[key]):
                best[key] = ob
        return list(best.values())

    def _outbreak_key(self, outbreak: dict) -> str:
        date = (outbreak.get("date") or "")[:10]
        title = outbreak.get("title", "")
        title = re.sub(r"\d{4}", "", title)
        title = re.sub(r"[^a-zA-Z\s]", " ", title)
        title = re.sub(r"\s+", " ", title).lower().strip()
        title = re.sub(r"\b(update|en|ebola|virus|disease|outbreak|item|html|\d+)\b", "", title).strip()
        country = (outbreak.get("countries") or [""])[0] if outbreak.get("countries") else ""
        return f"{date}|{title}|{country}"

    def _outbreak_score(self, outbreak: dict) -> int:
        score = outbreak.get("credibility", 5) * 10
        text = f"{outbreak.get('summary','')} {outbreak.get('epidemiology','')} {outbreak.get('overview','')} {outbreak.get('description','')}"
        score += len(text) / 50
        for field in ("cases", "deaths", "cfr", "summary", "overview", "epidemiology", "description"):
            if outbreak.get(field) is not None and outbreak.get(field) != "":
                score += 5
        return score


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Smoke-test the DDG epi fetcher")
    parser.add_argument("--species", default="Zaire ebolavirus")
    parser.add_argument("--lineage", default=None)
    parser.add_argument("--country", default="Uganda")
    parser.add_argument("--output", default="ddg_epi_fetch.json")
    args = parser.parse_args()

    fetcher = DuckDuckGoEpiFetcher()
    data = fetcher.fetch_all(
        species_name=args.species,
        lineage=args.lineage,
        country=args.country,
    )
    with open(args.output, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, default=str, ensure_ascii=False)
    print(f"Wrote {args.output}")
