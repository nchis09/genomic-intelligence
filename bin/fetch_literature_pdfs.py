#!/usr/bin/env python3
"""Download full-text PDFs for screened literature records.

Tries three resolvers in order for each paper:
1. Europe PMC Open Access via pyeuropepmc (FullTextClient with built-in
   render→backend→ZIP→Unpaywall cascade; PMCID resolved via SearchClient
   when not already present in the metadata)
2. Unpaywall direct lookup via pyeuropepmc (UnpaywallClient, for papers
   with a DOI but no resolvable PMCID)
3. https://doi.org/ DOI redirect with landing-page PDF link scraping
   (custom fallback for publishers not covered by the above)

Only files whose content-type is application/pdf (or whose body starts with the
%PDF magic bytes) are written. Papers that cannot be retrieved are recorded in
pdf_download_summary.json with a failed status.
"""
import argparse
import json
import os
import re
import sys
import tempfile
import time
from pathlib import Path
from typing import Any, Dict, List, Optional
from urllib.parse import urljoin

import requests

from pyeuropepmc.clients.search import SearchClient
from pyeuropepmc.clients.fulltext import FullTextClient
from pyeuropepmc.clients.unpaywall_client import UnpaywallClient


HEADERS = {
    "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/124.0.0.0 Safari/537.36",
    "Accept": "application/pdf, application/octet-stream, */*",
    "Accept-Language": "en-US,en;q=0.9",
}

# Publishers whose landing pages predictably expose a PDF at <landing_url>/pdf
# even when it wasn't discoverable via the generic HTML link scrape (e.g. due to
# anti-bot blocking or JS-rendered link markup).
PDF_SUFFIX_DOMAINS = ("mdpi.com",)


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
    url: str,
    timeout: int,
    retries: int,
    sleep: float,
    headers: Optional[Dict[str, str]] = None,
    referer: Optional[str] = None,
) -> Optional[requests.Response]:
    req_headers = dict(headers if headers is not None else HEADERS)
    if referer:
        req_headers["Referer"] = referer
    for attempt in range(retries):
        try:
            r = requests.get(url, headers=req_headers, timeout=timeout, allow_redirects=True)
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


def _resolve_pmcid(
    pmid: str, doi: str, search_client: SearchClient
) -> str:
    """Resolve a PMCID from PMID or DOI using the pyeuropepmc SearchClient."""
    for query, label in [
        (f"ext_id:{pmid} src:med", f"PMID {pmid}") if pmid else (None, None),
        (f"DOI:{doi}", f"DOI {doi}") if doi else (None, None),
    ]:
        if query is None:
            continue
        try:
            data = search_client.search(query, resultType="lite", pageSize=1, format="json")
            results = data.get("resultList", {}).get("result", []) if isinstance(data, dict) else []
            if results:
                pmcid = results[0].get("pmcid", "")
                if pmcid:
                    print(f"[fetch] {label} resolved to PMCID {pmcid} via Europe PMC search", file=sys.stderr)
                    return pmcid
        except Exception as exc:
            print(f"[fetch] Europe PMC search for {label} failed: {exc}", file=sys.stderr)
    return ""


def _try_europepmc(
    record: Dict[str, Any],
    timeout: int,
    retries: int,
    sleep: float,
    search_client: SearchClient,
    fulltext_client: FullTextClient,
) -> tuple[Optional[bytes], Optional[str]]:
    pmid = record.get("pmid", "")
    doi = record.get("doi", "")
    pmcid = record.get("pmcid", "")

    if not pmcid:
        pmcid = _resolve_pmcid(pmid, doi, search_client)

    if not pmcid:
        return None, "No PMCID available for Europe PMC OA lookup"

    # FullTextClient.download_pdf_by_pmcid tries (in order):
    #   1. Europe PMC render endpoint
    #   2. Backend render service (ptpmcrender.fcgi)
    #   3. ZIP archive from OA bulk collection
    #   4. Unpaywall via DOI lookup
    # It downloads to a file and validates PDF magic bytes.
    with tempfile.NamedTemporaryFile(suffix=".pdf", delete=False) as tmp:
        tmp_path = Path(tmp.name)
    try:
        result_path = fulltext_client.download_pdf_by_pmcid(pmcid, output_path=tmp_path)
        if result_path is not None and result_path.exists() and result_path.stat().st_size > 0:
            content = result_path.read_bytes()
            if content[:4] == b"%PDF":
                return content, None
            return None, "Europe PMC downloaded file is not a valid PDF"
        return None, "Europe PMC FullTextClient returned no PDF"
    except Exception as exc:
        return None, f"Europe PMC download error: {exc}"
    finally:
        tmp_path.unlink(missing_ok=True)


def _try_unpaywall(
    record: Dict[str, Any],
    email: str,
    timeout: int,
    retries: int,
    sleep: float,
    unpaywall_client: UnpaywallClient,
) -> tuple[Optional[bytes], Optional[str]]:
    doi = record.get("doi", "")
    if not doi:
        return None, "No DOI for Unpaywall lookup"

    try:
        pdf_url = unpaywall_client.get_pdf_url(doi)
    except Exception as exc:
        return None, f"Unpaywall API error: {exc}"

    if not pdf_url:
        return None, "Unpaywall has no PDF URL for this DOI"

    r = _safe_get(pdf_url, timeout, retries, sleep)
    if r is None:
        return None, "Unpaywall PDF download failed"
    if not _is_pdf_response(r):
        return None, "Unpaywall response was not a PDF"
    return r.content, None


def _extract_pdf_urls(base_url: str, html: str) -> List[str]:
    """Extract likely PDF URLs from a publisher landing page."""
    candidates: List[str] = []

    # citation_pdf_url meta tag
    for match in re.finditer(
        r'<meta[^>]+(?:name|property)=["\']citation_pdf_url["\'][^>]+content=["\']([^"\']+)["\']',
        html,
        re.IGNORECASE,
    ):
        candidates.append(match.group(1))
    for match in re.finditer(
        r'<meta[^>]+content=["\']([^"\']+\.pdf[^"\']*)["\'][^>]+(?:name|property)=["\']citation_pdf_url["\']',
        html,
        re.IGNORECASE,
    ):
        candidates.append(match.group(1))

    # Any link ending in .pdf
    for match in re.finditer(r'href=["\']([^"\']*\.pdf[^"\']*)["\']', html, re.IGNORECASE):
        candidates.append(match.group(1))

    # Any link whose text/aria-label says PDF
    for match in re.finditer(r'href=["\']([^"\']+)["\'][^>]*>([^<]*(?:PDF|pdf)[^<]*)', html, re.IGNORECASE):
        candidates.append(match.group(1))

    # <link rel="alternate" type="application/pdf" href="...">
    for match in re.finditer(
        r'<link[^>]+rel=["\']alternate["\'][^>]+type=["\']application/pdf["\'][^>]+href=["\']([^"\']+)["\']',
        html,
        re.IGNORECASE,
    ):
        candidates.append(match.group(1))
    for match in re.finditer(
        r'<link[^>]+href=["\']([^"\']+)["\'][^>]+rel=["\']alternate["\'][^>]+type=["\']application/pdf["\']',
        html,
        re.IGNORECASE,
    ):
        candidates.append(match.group(1))

    # og:pdf-url style meta tag used by some publisher CMSes
    for match in re.finditer(
        r'<meta[^>]+(?:name|property)=["\'][^"\']*pdf[^"\']*["\'][^>]+content=["\']([^"\']+\.pdf[^"\']*)["\']',
        html,
        re.IGNORECASE,
    ):
        candidates.append(match.group(1))

    # Deduplicate and resolve relative URLs
    seen: set = set()
    resolved: List[str] = []
    for raw in candidates:
        raw = raw.strip()
        if not raw:
            continue
        full = urljoin(base_url, raw)
        if full not in seen:
            seen.add(full)
            resolved.append(full)
    return resolved


def _try_doi_resolver(
    record: Dict[str, Any], timeout: int, retries: int, sleep: float
) -> tuple[Optional[bytes], Optional[str]]:
    doi = record.get("doi", "")
    if not doi:
        return None, "No DOI for resolver"

    url = f"https://doi.org/{doi}"
    html_headers = dict(HEADERS)
    html_headers["Accept"] = "text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8"
    # Some publishers (MDPI, CDC) block direct/bot-like requests with no Referer.
    # Pretending the request came from the DOI resolver itself mimics a normal
    # click-through and avoids a chunk of false 403s.
    html_headers["Referer"] = "https://doi.org/"

    landing_url = None
    landing_html = None
    last_error = None

    for attempt in range(retries):
        try:
            r = requests.get(url, headers=html_headers, timeout=timeout, allow_redirects=True)
            r.raise_for_status()
            if _is_pdf_response(r):
                return r.content, None

            landing_url = r.url
            landing_html = r.text
            break
        except requests.exceptions.HTTPError as exc:
            last_error = f"DOI resolver failed: {exc}"
            # Capture the final (blocked) landing URL even though raise_for_status
            # raised, so the MDPI-style /pdf suffix fallback below can still
            # target the real publisher URL rather than the bare doi.org link.
            if exc.response is not None:
                landing_url = exc.response.url
            # 403s are frequently anti-bot rules on the publisher's landing page
            # rather than a real access restriction; don't retry the same
            # blocked request repeatedly.
            break
        except Exception as exc:
            last_error = f"DOI resolver failed: {exc}"
            if attempt == retries - 1:
                break
            time.sleep(sleep)

    if landing_url is not None and landing_html is not None:
        candidates = _extract_pdf_urls(landing_url, landing_html)
        for pdf_url in candidates:
            pdf_r = _safe_get(pdf_url, timeout, 2, sleep, headers=html_headers, referer=landing_url)
            if pdf_r is not None and _is_pdf_response(pdf_r):
                return pdf_r.content, None
        if candidates:
            last_error = "DOI landing page PDF links did not return a PDF"
        else:
            last_error = "DOI landing page did not contain a PDF link"

    # Predictable-URL fallback for publishers whose landing pages block scraping
    # or use JS-rendered links our regex can't see (e.g. MDPI's <article>/pdf).
    probe_url = landing_url or url
    domain = probe_url.split("/")[2] if "://" in probe_url else ""
    if any(domain.endswith(d) for d in PDF_SUFFIX_DOMAINS):
        suffix_url = probe_url.rstrip("/") + "/pdf"
        pdf_r = _safe_get(suffix_url, timeout, retries, sleep, headers=html_headers, referer=probe_url)
        if pdf_r is not None and _is_pdf_response(pdf_r):
            return pdf_r.content, None

    return None, last_error or "DOI resolver failed"


def _download_pdf(
    record: Dict[str, Any],
    email: Optional[str],
    timeout: int,
    sleep: float,
    retries: int,
    search_client: SearchClient,
    fulltext_client: FullTextClient,
    unpaywall_client: Optional[UnpaywallClient],
) -> tuple[Dict[str, Any], Optional[bytes]]:
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

    attempts: List[Dict[str, str]] = []

    def _ep(rec: Dict[str, Any]) -> tuple[Optional[bytes], Optional[str]]:
        return _try_europepmc(rec, timeout, retries, sleep, search_client, fulltext_client)

    def _up(rec: Dict[str, Any]) -> tuple[Optional[bytes], Optional[str]]:
        if unpaywall_client is None:
            return None, "No email provided for Unpaywall API"
        return _try_unpaywall(rec, email, timeout, retries, sleep, unpaywall_client)

    def _dr(rec: Dict[str, Any]) -> tuple[Optional[bytes], Optional[str]]:
        return _try_doi_resolver(rec, timeout, retries, sleep)

    sources = [
        ("europepmc", _ep),
        ("unpaywall", _up),
        ("doi", _dr),
    ]

    for source, fn in sources:
        result["source"] = source
        content, err = fn(record)
        attempts.append({"source": source, "status": "success" if content is not None else "failed", "error": err or ""})
        if content is not None:
            result["status"] = "success"
            result["pdf_bytes"] = len(content)
            result["attempts"] = attempts
            return result, content
        result["error"] = err
        time.sleep(sleep)

    result["attempts"] = attempts
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

    # Instantiate pyeuropepmc clients once for reuse across all records.
    search_client = SearchClient()
    fulltext_client = FullTextClient()
    unpaywall_client = None
    if args.email:
        unpaywall_client = UnpaywallClient(email=args.email)

    summary: List[Dict[str, Any]] = []
    for rec in records:
        result, content = _download_pdf(
            rec,
            email=args.email,
            timeout=args.timeout,
            sleep=args.sleep,
            retries=args.retries,
            search_client=search_client,
            fulltext_client=fulltext_client,
            unpaywall_client=unpaywall_client,
        )
        if content is not None:
            _write_pdf(outdir, result["pmid"], content)
            print(f"[fetch] Downloaded PDF for PMID {result['pmid']} from {result['source']}", file=sys.stderr)
        else:
            print(f"[fetch] No PDF for PMID {result['pmid']}: {result['error']}", file=sys.stderr)
        summary.append(result)
        time.sleep(args.sleep)

    # Clean up client resources.
    try:
        search_client.close()
    except Exception:
        pass
    try:
        fulltext_client.close()
    except Exception:
        pass
    if unpaywall_client is not None:
        try:
            unpaywall_client.close()
        except Exception:
            pass

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
