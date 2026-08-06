#!/usr/bin/env python3
"""Download full-text PDFs for screened literature records.

Tries three resolvers in order for each paper:
1. Europe PMC Open Access (using PMCID if available, otherwise searching by PMID)
2. Unpaywall (using DOI and a contact email)
3. https://doi.org/ DOI redirect (Accept: application/pdf)

Only files whose content-type is application/pdf (or whose body starts with the
%PDF magic bytes) are written. Papers that cannot be retrieved are recorded in
pdf_download_summary.json with a failed status.
"""
import argparse
import json
import os
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional

import requests


HEADERS = {
    "User-Agent": "Mozilla/5.0 (compatible; pgirl-genomic-intelligence/1.0; +https://github.com/pgirl/genomic-intelligence)",
    "Accept": "application/pdf, application/octet-stream, */*",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Download full-text PDFs for screened literature records.")
    parser.add_argument("--input-dir", required=True, help="Directory with screened paper JSON files.")
    parser.add_argument("--outdir", required=True, help="Directory for PDF output and summary.")
    parser.add_argument("--species", required=True, help="Species key.")
    parser.add_argument("--domain", required=True, help="Domain key.")
    parser.add_argument("--email", default=None, help="Contact email for Unpaywall API.")
    parser.add_argument("--timeout", type=int, default=30, help="HTTP request timeout in seconds.")
    parser.add_argument("--sleep", type=float, default=1.0, help="Seconds to sleep between requests.")
    parser.add_argument("--retries", type=int, default=3, help="Retries per resolver.")
    return parser.parse_args()


def _is_pdf_response(response: requests.Response) -> bool:
    content_type = (response.headers.get("Content-Type") or "").lower()
    if "application/pdf" in content_type:
        return True
    # Some publishers serve application/octet-stream; check magic bytes.
    if response.content[:4] == b"%PDF":
        return True
    return False


def _safe_get(
    url: str, timeout: int, retries: int, sleep: float
) -> Optional[requests.Response]:
    for attempt in range(retries):
        try:
            r = requests.get(url, headers=HEADERS, timeout=timeout, allow_redirects=True)
            r.raise_for_status()
            return r
        except Exception:
            if attempt == retries - 1:
                return None
            time.sleep(sleep)
    return None


def _write_pdf(outdir: Path, pmid: str, content: bytes) -> Path:
    pdf_path = outdir / f"{pmid}.pdf"
    with open(pdf_path, "wb") as fh:
        fh.write(content)
    return pdf_path


def _try_europepmc(
    record: Dict[str, Any], timeout: int, retries: int, sleep: float
) -> tuple[Optional[bytes], Optional[str]]:
    pmid = record.get("pmid", "")
    pmcid = record.get("pmcid", "")

    if not pmcid and pmid:
        # Try to resolve PMCID from Europe PMC search.
        search_url = (
            f"https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=ext_id:{pmid}%20src:med"
            f"&format=json&resultType=core&pageSize=1"
        )
        r = _safe_get(search_url, timeout, retries, sleep)
        if r is not None:
            try:
                data = r.json()
                results = data.get("resultList", {}).get("result", [])
                if results:
                    pmcid = results[0].get("pmcid", "")
            except Exception:
                pass

    if not pmcid:
        return None, "No PMCID available for Europe PMC OA lookup"

    # OA PDF endpoint: works for OA articles in PubMed Central.
    pdf_url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/pdf/oa"
    r = _safe_get(pdf_url, timeout, retries, sleep)
    if r is None:
        return None, "Europe PMC OA PDF request failed"
    if not _is_pdf_response(r):
        return None, "Europe PMC response was not a PDF"
    return r.content, None


def _try_unpaywall(
    record: Dict[str, Any], email: str, timeout: int, retries: int, sleep: float
) -> tuple[Optional[bytes], Optional[str]]:
    doi = record.get("doi", "")
    if not doi:
        return None, "No DOI for Unpaywall lookup"
    if not email:
        return None, "No email provided for Unpaywall API"

    url = f"https://api.unpaywall.org/v2/{doi}?email={email}"
    r = _safe_get(url, timeout, retries, sleep)
    if r is None:
        return None, "Unpaywall API request failed"

    try:
        data = r.json()
    except Exception:
        return None, "Unpaywall returned invalid JSON"

    if not data.get("is_oa"):
        return None, "Unpaywall reports no open access location"

    location = data.get("best_oa_location") or {}
    pdf_url = location.get("url_for_pdf") or ""
    if not pdf_url:
        return None, "Unpaywall has no url_for_pdf"

    r = _safe_get(pdf_url, timeout, retries, sleep)
    if r is None:
        return None, "Unpaywall PDF download failed"
    if not _is_pdf_response(r):
        return None, "Unpaywall response was not a PDF"
    return r.content, None


def _try_doi_resolver(
    record: Dict[str, Any], timeout: int, retries: int, sleep: float
) -> tuple[Optional[bytes], Optional[str]]:
    doi = record.get("doi", "")
    if not doi:
        return None, "No DOI for resolver"

    url = f"https://doi.org/{doi}"
    headers = dict(HEADERS)
    headers["Accept"] = "application/pdf"
    for attempt in range(retries):
        try:
            r = requests.get(url, headers=headers, timeout=timeout, allow_redirects=True)
            r.raise_for_status()
            if _is_pdf_response(r):
                return r.content, None
            return None, "DOI resolver did not return a PDF"
        except Exception as exc:
            if attempt == retries - 1:
                return None, f"DOI resolver failed: {exc}"
            time.sleep(sleep)
    return None, "DOI resolver failed"


def _download_pdf(
    record: Dict[str, Any],
    email: Optional[str],
    timeout: int,
    sleep: float,
    retries: int,
) -> Dict[str, Any]:
    pmid = record.get("pmid", "") or ""
    doi = record.get("doi", "") or ""
    pmcid = record.get("pmcid", "") or ""

    result: Dict[str, Any] = {
        "pmid": pmid,
        "doi": doi,
        "pmcid": pmcid,
        "status": "failed",
        "source": None,
        "error": None,
    }

    sources = [
        ("europepmc", lambda rec: _try_europepmc(rec, timeout, retries, sleep)),
        ("unpaywall", lambda rec: _try_unpaywall(rec, email, timeout, retries, sleep)),
        ("doi", lambda rec: _try_doi_resolver(rec, timeout, retries, sleep)),
    ]

    for source, fn in sources:
        content, err = fn(record)
        if content is not None:
            result["source"] = source
            result["status"] = "success"
            result["pdf_bytes"] = len(content)
            return result, content
        result["error"] = err
        time.sleep(sleep)

    return result, None


def load_records(input_dir: Path) -> List[Dict[str, Any]]:
    records: List[Dict[str, Any]] = []
    for path in sorted(input_dir.glob("*.json")):
        if path.name == "screening_summary.json" or path.name.endswith("_summary.json"):
            continue
        try:
            with open(path, "r", encoding="utf-8") as fh:
                rec = json.load(fh)
            if rec.get("pmid"):
                records.append(rec)
        except Exception as exc:
            print(f"[warn] Could not read {path}: {exc}", file=sys.stderr)
    return records


def main() -> None:
    args = parse_args()
    input_dir = Path(args.input_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    records = load_records(input_dir)
    if not records:
        print("[fetch] No paper JSONs found to download.", file=sys.stderr)

    summary: List[Dict[str, Any]] = []
    for rec in records:
        result, content = _download_pdf(
            rec,
            email=args.email,
            timeout=args.timeout,
            sleep=args.sleep,
            retries=args.retries,
        )
        if content is not None:
            _write_pdf(outdir, result["pmid"], content)
            print(f"[fetch] Downloaded PDF for PMID {result['pmid']} from {result['source']}", file=sys.stderr)
        else:
            print(f"[fetch] No PDF for PMID {result['pmid']}: {result['error']}", file=sys.stderr)
        summary.append(result)
        time.sleep(args.sleep)

    summary_record = {
        "species": args.species,
        "domain": args.domain,
        "input_count": len(records),
        "success_count": sum(1 for s in summary if s["status"] == "success"),
        "failed_count": sum(1 for s in summary if s["status"] != "success"),
        "results": summary,
    }
    summary_path = outdir / "pdf_download_summary.json"
    with open(summary_path, "w", encoding="utf-8") as fh:
        json.dump(summary_record, fh, ensure_ascii=False, indent=2)
    print(f"[fetch] Wrote {summary_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
