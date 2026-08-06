#!/usr/bin/env python3
"""Fetch PubMed metadata for a single domain's Europe PMC results."""
import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any, Dict, List

from Bio import Entrez


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Fetch PubMed metadata for one LITERATURE_SEARCH domain."
    )
    parser.add_argument(
        "--results-json", required=True, help="Path to LITERATURE_SEARCH results.json"
    )
    parser.add_argument("--species", required=True, help="Species key")
    parser.add_argument("--domain", required=True, help="Domain key")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument(
        "--email", default=None, help="NCBI Entrez contact email"
    )
    parser.add_argument(
        "--batch-size", type=int, default=100, help="PMIDs per EFetch request"
    )
    parser.add_argument(
        "--sleep", type=float, default=0.5, help="Seconds to sleep between requests"
    )
    parser.add_argument(
        "--retries", type=int, default=3, help="Retries per batch"
    )
    return parser.parse_args()


def extract_text(obj: Any) -> str:
    """Recursively extract text from a Biopython Entrez string/list/dict."""
    if obj is None:
        return ""
    if isinstance(obj, str):
        return obj
    if isinstance(obj, list):
        return " ".join(extract_text(item) for item in obj).strip()
    if isinstance(obj, dict):
        if "#text" in obj:
            return str(obj["#text"])
        return " ".join(extract_text(v) for v in obj.values()).strip()
    return str(obj)


def load_results(path: Path) -> List[Dict[str, Any]]:
    with open(path, "r", encoding="utf-8") as fh:
        data = json.load(fh)
    return data.get("results", [])


def parse_article(article: Dict[str, Any]) -> Dict[str, Any]:
    medline = article.get("MedlineCitation", {})
    pmid = str(medline.get("PMID", ""))
    article_data = medline.get("Article", {})

    title = extract_text(article_data.get("ArticleTitle", ""))

    abstract_node = article_data.get("Abstract", {})
    abstract = extract_text(abstract_node.get("AbstractText", ""))

    journal = article_data.get("Journal", {})
    journal_title = extract_text(journal.get("Title", ""))
    journal_issue = journal.get("JournalIssue", {})
    pub_date = journal_issue.get("PubDate", {})
    year = pub_date.get("Year", "")
    if not year:
        medline_date = pub_date.get("MedlineDate", "")
        if medline_date:
            year = str(medline_date)[:4]
    if not year:
        article_date = article_data.get("ArticleDate", {})
        year = article_date.get("Year", "")

    authors = []
    for author in article_data.get("AuthorList", []):
        if not isinstance(author, dict):
            continue
        name = ""
        if "CollectiveName" in author:
            name = extract_text(author["CollectiveName"])
        else:
            parts = []
            last = author.get("LastName", "")
            if last:
                parts.append(last)
            fore = author.get("ForeName", "") or author.get("Initials", "")
            if fore:
                parts.append(fore)
            name = " ".join(parts)
        if name:
            authors.append(name)

    article_ids: Dict[str, str] = {}
    for aid in article_data.get("ArticleIdList", []):
        if not isinstance(aid, dict):
            continue
        id_type = aid.get("IdType", "")
        value = str(aid.get("value", ""))
        if id_type and value:
            article_ids[id_type] = value

    doi = article_ids.get("doi", "")
    pmcid = article_ids.get("pmc", "")

    keywords = []
    for k in medline.get("KeywordList", []):
        keywords.append(extract_text(k))
    keywords = [k for k in keywords if k]

    return {
        "pmid": pmid,
        "title": title,
        "abstract": abstract,
        "authors": authors,
        "year": year,
        "doi": doi,
        "journal": journal_title,
        "keywords": keywords,
        "pmcid": pmcid,
    }


def fetch_pubmed(
    pmids: List[str],
    email: str,
    batch_size: int,
    sleep: float,
    retries: int,
) -> Dict[str, Dict[str, Any]]:
    Entrez.email = email
    results: Dict[str, Dict[str, Any]] = {}

    for i in range(0, len(pmids), batch_size):
        batch = pmids[i : i + batch_size]
        for attempt in range(retries):
            try:
                handle = Entrez.efetch(db="pubmed", id=batch, retmode="xml")
                record = Entrez.read(handle)
                for article in record.get("PubmedArticle", []):
                    parsed = parse_article(article)
                    if parsed["pmid"]:
                        results[parsed["pmid"]] = parsed
                handle.close()
                break
            except Exception as exc:
                print(
                    f"[fetch] Batch {i} failed (attempt {attempt + 1}/{retries}): {exc}",
                    file=sys.stderr,
                )
                if attempt == retries - 1:
                    print(f"[fetch] Giving up on batch {i}", file=sys.stderr)
                time.sleep(sleep)
        time.sleep(sleep)

    return results


def main() -> None:
    args = parse_args()
    if not args.email:
        print(
            "Error: NCBI Entrez requires a contact email. Set --email or params.pubmed_email.",
            file=sys.stderr,
        )
        sys.exit(1)

    works = load_results(Path(args.results_json))
    if not works:
        print("[fetch] No publications in results.json", file=sys.stderr)
        Path(args.outdir).mkdir(parents=True, exist_ok=True)
        return

    pmids = []
    original: Dict[str, Dict[str, Any]] = {}
    for w in works:
        pmid = str(w.get("id", "")).strip()
        if not pmid:
            continue
        pmids.append(pmid)
        original[pmid] = w

    numeric_pmids = [p for p in pmids if p.isdigit()]
    pubmed: Dict[str, Dict[str, Any]] = {}
    if numeric_pmids:
        pubmed = fetch_pubmed(
            numeric_pmids,
            args.email,
            args.batch_size,
            args.sleep,
            args.retries,
        )
    else:
        print("[fetch] No numeric PMIDs to query PubMed", file=sys.stderr)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    for pmid in pmids:
        original_work = original.get(pmid, {})

        if pmid in pubmed:
            record = pubmed[pmid]
        else:
            record = {
                "pmid": pmid,
                "title": original_work.get("display_name", ""),
                "abstract": "",
                "authors": [original_work.get("first_author", "")] if original_work.get("first_author") else [],
                "year": original_work.get("publication_year", ""),
                "doi": original_work.get("doi", ""),
                "journal": original_work.get("source", ""),
                "keywords": [],
                "pmcid": "",
            }

        record["species"] = args.species
        record["domain"] = args.domain
        record["publication_date"] = original_work.get("publication_date", "")
        record["cited_by_count"] = original_work.get("cited_by_count", 0)
        record["is_oa"] = original_work.get("is_oa", "")

        out_path = outdir / f"{pmid}.json"
        with open(out_path, "w", encoding="utf-8") as fh:
            json.dump(record, fh, ensure_ascii=False, indent=2)
        print(f"[fetch] Wrote {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
