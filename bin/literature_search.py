#!/usr/bin/env python3
"""
Search OpenAlex for per-species, per-domain literature and save results as JSON/TSV.

One OpenAlex query is built per domain using the species synonyms from the YAML
and a separate `title_and_abstract.search` filter is used for the domain terms.
The two filters are AND-ed by repeating the `title_and_abstract.search` filter.
"""

import argparse
import csv
import json
import sys
import time
import datetime
from pathlib import Path

import requests
import yaml

BASE_URL = "https://api.openalex.org/works"
SELECT_FIELDS = (
    "id,display_name,doi,publication_year,publication_date,"
    "abstract_inverted_index,authorships,primary_location,locations,"
    "type,open_access,cited_by_count"
)


def quote_term(term: str) -> str:
    """Quote multi-word or special terms for OpenAlex filter values."""
    if " " in term or "|" in term or "," in term or "+" in term:
        term = term.replace('"', '\\"')
        return f'"{term}"'
    return term


def build_pipe_group(terms: list) -> str:
    """Join terms with the OpenAlex OR pipe."""
    return "|".join(quote_term(t) for t in terms)


def build_filter_value(species_terms: list, domain_terms: list, from_date: str = None) -> str:
    """
    Build the OpenAlex `filter` value.
    Two `title_and_abstract.search` filters are AND-ed by repeating the filter
    key and separating with a comma. Optionally AND with a `from_publication_date`.
    """
    species_group = build_pipe_group(species_terms)
    domain_group = build_pipe_group(domain_terms)
    filter_value = (
        f"title_and_abstract.search:{species_group},"
        f"title_and_abstract.search:{domain_group}"
    )
    if from_date:
        filter_value += f",from_publication_date:{from_date}"
    return filter_value


def abstract_from_inverted_index(inv_index: dict) -> str:
    """Reconstruct plain text from OpenAlex abstract_inverted_index."""
    if not inv_index:
        return ""
    positions = {}
    for word, indices in inv_index.items():
        for pos in indices:
            positions[pos] = word
    if not positions:
        return ""
    max_pos = max(positions.keys())
    words = [positions.get(i, "") for i in range(max_pos + 1)]
    return " ".join(words)


MAX_RETRY_WAIT = 30.0  # seconds; never honor a Retry-After longer than this


def fetch_page(params: dict, timeout: int = 60, max_retries: int = 5) -> dict:
    """Fetch one page from the OpenAlex API, with 429 backoff."""
    for attempt in range(max_retries + 1):
        response = requests.get(BASE_URL, params=params, timeout=timeout)
        if response.status_code == 429:
            if attempt == max_retries:
                response.raise_for_status()
            retry_after = response.headers.get("Retry-After")
            if retry_after:
                try:
                    wait = min(float(retry_after), MAX_RETRY_WAIT)
                except ValueError:
                    wait = min(2 ** attempt, MAX_RETRY_WAIT)
            else:
                wait = min(2 ** attempt, MAX_RETRY_WAIT)
            print(
                f"OpenAlex 429: rate limited. Waiting {wait}s (retry {attempt + 1}/{max_retries}) ...",
                file=sys.stderr,
            )
            time.sleep(wait)
            continue
        response.raise_for_status()
        return response.json()
    raise requests.exceptions.HTTPError("Max retries exceeded for OpenAlex request.")


def search_openalex(
    species_terms: list,
    domain_terms: list,
    max_results: int,
    mailto: str = None,
    api_key: str = None,
    per_page: int = 100,
    sleep: float = 1.0,
    max_retries: int = 2,
) -> list:
    """Paginate OpenAlex `filter` searches up to `max_results`."""
    results = []
    cursor = "*"

    ten_years_ago = (datetime.date.today() - datetime.timedelta(days=10 * 365)).isoformat()

    while len(results) < max_results:
        params = {
            "filter": build_filter_value(species_terms, domain_terms, from_date=ten_years_ago),
            "per_page": per_page,
            "cursor": cursor,
            "select": SELECT_FIELDS,
        }
        if mailto:
            params["mailto"] = mailto
        if api_key:
            params["api_key"] = api_key

        data = fetch_page(params, max_retries=max_retries)
        page_results = data.get("results", [])
        meta = data.get("meta", {})

        if not page_results:
            break

        results.extend(page_results)

        next_cursor = meta.get("next_cursor")
        if not next_cursor or len(page_results) < per_page:
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
            "openalex_id",
            "title",
            "doi",
            "year",
            "publication_date",
            "first_author",
            "source",
            "is_oa",
            "cited_by_count",
            "abstract",
        ])
        for work in works:
            openalex_id = work.get("id", "")
            title = work.get("display_name", "")
            doi = work.get("doi", "")
            year = work.get("publication_year", "")
            pub_date = work.get("publication_date", "")

            first_author = ""
            authorships = work.get("authorships", [])
            if authorships:
                author = authorships[0].get("author", {})
                first_author = author.get("display_name", "")

            source = ""
            primary = work.get("primary_location") or {}
            if primary and primary.get("source"):
                source = primary["source"].get("display_name", "")

            is_oa = ""
            open_access = work.get("open_access", {})
            if open_access:
                is_oa = open_access.get("is_oa", "")

            cited = work.get("cited_by_count", "")
            abstract = abstract_from_inverted_index(work.get("abstract_inverted_index", {}))

            writer.writerow([
                openalex_id,
                title,
                doi,
                year,
                pub_date,
                first_author,
                source,
                is_oa,
                cited,
                abstract,
            ])


def main():
    parser = argparse.ArgumentParser(
        description="Search OpenAlex per species and domain and save results."
    )
    parser.add_argument("--species", required=True, help="Species key, e.g. ebov")
    parser.add_argument("--terms-yaml", required=True, help="Path to literature_search_terms.yml")
    parser.add_argument("--max-results", type=int, default=1000, help="Max results per domain")
    parser.add_argument("--mailto", default=None, help="Email for OpenAlex polite pool")
    parser.add_argument("--api-key", default=None, help="OpenAlex API key")
    parser.add_argument("--outdir", required=True, help="Output directory")
    args = parser.parse_args()

    with open(args.terms_yaml, "r", encoding="utf-8") as fh:
        config = yaml.safe_load(fh)

    species_synonyms = config.get("species_synonyms", {}).get(args.species, [])
    shared_terms = config.get("shared_terms", [])

    if not species_synonyms:
        print(
            f"Warning: no synonyms for species '{args.species}', using shared terms only.",
            file=sys.stderr,
        )
        species_terms = list(shared_terms)
    else:
        species_terms = list(dict.fromkeys(species_synonyms + shared_terms))

    domains = config.get("domains", {})
    if not domains:
        print("Error: no domains defined in terms YAML.", file=sys.stderr)
        sys.exit(1)

    outdir = Path(args.outdir)
    for domain_key, domain_info in domains.items():
        domain_terms = domain_info.get("terms", [])
        if not domain_terms:
            continue

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

        print(f"[{args.species}/{domain_key}] Searching OpenAlex ...", file=sys.stderr)
        try:
            works = search_openalex(
                species_terms=species_terms,
                domain_terms=domain_terms,
                max_results=args.max_results,
                mailto=args.mailto,
                api_key=args.api_key,
            )
        except requests.exceptions.RequestException as exc:
            print(
                f"[{args.species}/{domain_key}] Search failed ({exc}). Writing empty results.",
                file=sys.stderr,
            )
            works = []
        save_domain(outdir, domain_key, works)
        print(f"[{args.species}/{domain_key}] Saved {len(works)} results.", file=sys.stderr)


if __name__ == "__main__":
    main()
