#!/usr/bin/env python3
"""
Search Europe PMC for per-species, per-domain literature and save results as JSON/TSV.

One Europe PMC query is built per domain using the species synonyms from the YAML.
The species group, domain group, and a 10-year publication date filter are AND-ed.
"""

import argparse
import csv
import html
import json
import sys
import time
import datetime
import re
from pathlib import Path

import requests
import yaml

BASE_URL = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"


def quote_term(term: str) -> str:
    """Quote multi-word or special phrases for Europe PMC."""
    if " " in term or "|" in term or "," in term or "+" in term:
        term = term.replace('"', '\\"')
        return f'"{term}"'
    return term


def build_query_value(terms: list) -> str:
    """Join terms with the Europe PMC OR operator."""
    return " OR ".join(quote_term(t) for t in terms)


def build_query(species_terms: list, domain_terms: list, from_date: str, to_date: str) -> str:
    """
    Build a Europe PMC `query` string.
    The species and domain term groups are AND-ed with a 10-year date range.
    """
    species_group = build_query_value(species_terms)
    domain_group = build_query_value(domain_terms)
    return f"({species_group}) AND ({domain_group}) AND FIRST_PDATE:[{from_date} TO {to_date}]"


def work_text(work: dict) -> str:
    """Combine title and abstract for local term matching."""
    title = work.get("display_name", "")
    abstract = work.get("abstract", "")
    if not isinstance(abstract, str):
        abstract = ""
    return f"{title} {abstract}".lower()


def normalize_work(hit: dict) -> dict:
    """Map a Europe PMC hit to the common schema used for filtering and saving."""
    raw_abstract = hit.get("abstractText")
    if isinstance(raw_abstract, str):
        abstract = re.sub(r"<[^>]+>", "", raw_abstract)
        abstract = html.unescape(abstract)
    else:
        abstract = ""
    author = hit.get("authorString") or ""
    first_author = author.split(",")[0] if author else ""
    return {
        "id": hit.get("pmid") or hit.get("id", ""),
        "display_name": hit.get("title", ""),
        "doi": hit.get("doi", ""),
        "publication_year": hit.get("pubYear", ""),
        "publication_date": hit.get("firstPublicationDate", ""),
        "abstract": abstract,
        "cited_by_count": int(hit.get("citedByCount") or 0),
        "first_author": first_author,
        "source": hit.get("journalTitle", ""),
        "is_oa": hit.get("isOpenAccess", ""),
    }


def term_hits(text: str, terms: list) -> int:
    """Count how many terms appear in text as whole phrases (case-insensitive)."""
    if not terms:
        return 0
    return sum(
        1
        for term in terms
        if re.search(r"\b" + re.escape(term) + r"\b", text, re.IGNORECASE)
    )


def filter_and_score(
    works: list,
    exact_terms: list,
    broad_terms: list,
    shared_terms: list,
    excluded_terms: list,
    required_terms: list,
    related_terms: list,
    domain_excluded_terms: list = None,
) -> list:
    """Drop off-topic works and score the rest by species/domain relevance."""
    domain_excluded_terms = domain_excluded_terms or []
    filtered = []
    for work in works:
        text = work_text(work)

        # Drop if it mentions another species or Marburg.
        if term_hits(text, excluded_terms) > 0:
            continue

        # Drop if it matches a domain-level exclusion (e.g. animal-model-only
        # studies in the human "clinical" domain).
        if term_hits(text, domain_excluded_terms) > 0:
            continue

        exact_hits = term_hits(text, exact_terms)
        broad_hits = term_hits(text, broad_terms)
        shared_hits = term_hits(text, shared_terms)
        # Require at least one species-specific (exact or broad) term.
        # Shared terms like "orthoebolavirus" are too generic on their own.
        if exact_hits + broad_hits == 0:
            continue

        required_hits = term_hits(text, required_terms)
        related_hits = term_hits(text, related_terms)
        if required_hits == 0 and related_hits == 0:
            # No domain signal at all.
            continue

        score = (
            10 * exact_hits
            + 5 * broad_hits
            + 2 * shared_hits
            + 4 * required_hits
            + 1 * related_hits
            + 0.01 * work.get("cited_by_count", 0)
        )
        work["_relevance_score"] = score
        filtered.append(work)

    return filtered


MAX_RETRY_WAIT = 60.0  # seconds; never honor a Retry-After longer than this


def fetch_page(params: dict, timeout: int = 60, max_retries: int = 5) -> dict:
    """Fetch one page from the Europe PMC API, with 429 backoff."""
    for attempt in range(max_retries + 1):
        time.sleep(1.0)  # Throttle to stay within Europe PMC per-second limits
        response = requests.get(BASE_URL, params=params, timeout=timeout)
        if response.status_code == 429:
            if attempt == max_retries:
                response.raise_for_status()
            retry_after = response.headers.get("Retry-After")
            if retry_after:
                try:
                    # Increase backoff on each 429 even when the server suggests a value
                    wait = min(float(retry_after) + (attempt * 15), MAX_RETRY_WAIT)
                except ValueError:
                    wait = min((2 ** attempt) * 10, MAX_RETRY_WAIT)
            else:
                wait = min((2 ** attempt) * 10, MAX_RETRY_WAIT)
            print(
                f"Europe PMC 429: rate limited. Waiting {wait}s (retry {attempt + 1}/{max_retries}) ...",
                file=sys.stderr,
            )
            time.sleep(wait)
            continue
        response.raise_for_status()
        return response.json()
    raise requests.exceptions.HTTPError("Max retries exceeded for Europe PMC request.")


def search_europepmc(
    species_terms: list,
    domain_terms: list,
    max_results: int,
    per_page: int = 100,
    sleep: float = 1.0,
    max_retries: int = 2,
) -> list:
    """Paginate Europe PMC `search` requests up to `max_results`."""
    results = []
    cursor = "*"

    min_date = (datetime.date.today() - datetime.timedelta(days=10 * 365)).isoformat()
    max_date = datetime.date.today().isoformat()

    while len(results) < max_results:
        query = build_query(species_terms, domain_terms, min_date, max_date)
        params = {
            "query": query,
            "format": "json",
            "resultType": "core",
            "pageSize": min(per_page, max_results),
            "cursorMark": cursor,
        }

        data = fetch_page(params, max_retries=max_retries)
        page_results = data.get("resultList", {}).get("result", [])

        if not page_results:
            break

        for hit in page_results:
            results.append(normalize_work(hit))

        next_cursor = data.get("nextCursorMark")
        if not next_cursor or next_cursor == cursor or len(page_results) < per_page:
            break

        cursor = next_cursor
        time.sleep(sleep)

    return results[:max_results]


def save_domain(outdir: Path, domain: str, works: list) -> None:
    """Write `results.json` and `results.tsv` for a single domain."""
    domain_dir = Path(outdir) / domain
    domain_dir.mkdir(parents=True, exist_ok=True)

    json_path = domain_dir / "results.json"
    with open(json_path, "w", encoding="utf-8") as fh:
        json.dump({"count": len(works), "results": works}, fh, ensure_ascii=False, indent=2)

    tsv_path = domain_dir / "results.tsv"
    with open(tsv_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow([
            "id",
            "title",
            "doi",
            "year",
            "publication_date",
            "first_author",
            "source",
            "is_oa",
            "cited_by_count",
        ])
        for work in works:
            writer.writerow([
                work.get("id", ""),
                work.get("display_name", ""),
                work.get("doi", ""),
                work.get("publication_year", ""),
                work.get("publication_date", ""),
                work.get("first_author", ""),
                work.get("source", ""),
                work.get("is_oa", ""),
                work.get("cited_by_count", 0),
            ])


def main():
    parser = argparse.ArgumentParser(
        description="Search Europe PMC per species and domain and save results."
    )
    parser.add_argument("--species", required=True, help="Species key, e.g. ebov")
    parser.add_argument("--terms-yaml", required=True, help="Path to literature_search_terms.yml")
    parser.add_argument("--max-results", type=int, default=1000, help="Max results per domain")
    parser.add_argument("--mailto", default=None, help="Unused legacy OpenAlex email argument")
    parser.add_argument("--api-key", default=None, help="Unused legacy OpenAlex API key argument")
    parser.add_argument("--outdir", required=True, help="Output directory")
    args = parser.parse_args()

    with open(args.terms_yaml, "r", encoding="utf-8") as fh:
        config = yaml.safe_load(fh)

    species_cfg = config.get("species_synonyms", {}).get(args.species, {})
    exact_terms = species_cfg.get("exact", [])
    broad_terms = species_cfg.get("broad", [])
    shared_terms = config.get("shared_terms", [])
    excluded_terms = config.get("excluded_terms", {}).get(args.species, [])

    species_search_terms = list(dict.fromkeys(exact_terms + broad_terms + shared_terms))
    if not species_search_terms:
        print(
            f"Error: no species terms for '{args.species}' in the terms YAML.",
            file=sys.stderr,
        )
        sys.exit(1)

    domains = config.get("domains", {})
    if not domains:
        print("Error: no domains defined in terms YAML.", file=sys.stderr)
        sys.exit(1)

    all_domain_excluded_terms = config.get("domain_excluded_terms", {})

    outdir = Path(args.outdir)
    for domain_key, domain_info in domains.items():
        required_terms = domain_info.get("required", [])
        related_terms = domain_info.get("related", [])
        domain_search_terms = list(dict.fromkeys(required_terms + related_terms))
        if not domain_search_terms:
            continue
        domain_excluded_terms = all_domain_excluded_terms.get(domain_key, [])

        json_path = outdir / domain_key / "results.json"
        if json_path.is_file() and json_path.stat().st_size > 0:
            try:
                with open(json_path, "r", encoding="utf-8") as fh:
                    existing = json.load(fh)
                print(
                    f"[{args.species}/{domain_key}] Found existing {len(existing.get('results', []))} results, skipping.",
                    file=sys.stderr,
                )
                continue
            except json.JSONDecodeError:
                # Corrupt/partial file; re-download.
                pass

        print(f"[{args.species}/{domain_key}] Searching Europe PMC ...", file=sys.stderr)
        try:
            raw = search_europepmc(
                species_terms=species_search_terms,
                domain_terms=domain_search_terms,
                max_results=args.max_results,
            )
            works = filter_and_score(
                raw,
                exact_terms=exact_terms,
                broad_terms=broad_terms,
                shared_terms=shared_terms,
                excluded_terms=excluded_terms,
                required_terms=required_terms,
                related_terms=related_terms,
                domain_excluded_terms=domain_excluded_terms,
            )
            works.sort(key=lambda w: w.get("_relevance_score", 0), reverse=True)
            for w in works:
                w.pop("_relevance_score", None)
                w.pop("abstract", None)
            works = works[: args.max_results]
        except requests.exceptions.RequestException as exc:
            print(
                f"[{args.species}/{domain_key}] Search failed ({exc}). Writing empty results.",
                file=sys.stderr,
            )
            works = []
        save_domain(outdir, domain_key, works)
        print(f"[{args.species}/{domain_key}] Saved {len(works)} results.", file=sys.stderr)
        time.sleep(3.0)  # Space out Europe PMC domain requests


if __name__ == "__main__":
    main()
