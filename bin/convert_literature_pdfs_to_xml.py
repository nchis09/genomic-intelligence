#!/usr/bin/env python3
"""Convert downloaded PDFs to structured TEI XML using GROBID.

Uses the official grobid-client-python library when possible, falling back to
direct HTTP POST calls to the GROBID /api/processFulltextDocument endpoint.
"""
import argparse
import json
import shutil
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import requests


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Convert PDFs to TEI XML with GROBID.")
    parser.add_argument("--input-dir", required=True, help="Directory with PDF files.")
    parser.add_argument("--outdir", required=True, help="Directory for TEI XML output and summary.")
    parser.add_argument("--grobid-url", required=True, help="GROBID server URL, e.g. http://localhost:8070")
    parser.add_argument("--timeout", type=int, default=60, help="HTTP timeout in seconds.")
    parser.add_argument("--sleep", type=float, default=1.0, help="Seconds between GROBID requests.")
    parser.add_argument("--retries", type=int, default=3, help="Retries per PDF.")
    parser.add_argument("--batch-size", type=int, default=10, help="Ignored; kept for future batching.")
    return parser.parse_args()


def _read_pdf(path: Path) -> bytes:
    with open(path, "rb") as fh:
        return fh.read()


def _looks_tei(content: bytes) -> bool:
    text = content[:200].decode("utf-8", errors="ignore").strip()
    return text.startswith("<?xml") or text.startswith("<TEI") or text.startswith("<tei:")


def _pmid_from_filename(pdf_path: Path) -> str:
    # {pmid}.pdf
    return pdf_path.stem


def _find_client_output(outdir: Path, pmid: str) -> Optional[Path]:
    for path in outdir.glob(f"{pmid}*.tei.xml"):
        return path
    return None


def _try_grobid_client(
    pdfs: List[Path], outdir: Path, grobid_input: Path, grobid_url: str
) -> Tuple[Dict[str, Any], Optional[Dict[str, Any]]]:
    try:
        from grobid_client.grobid_client import GrobidClient
    except Exception as exc:
        return {}, {"source": "grobid_client", "error": f"import failed: {exc}"}

    try:
        client = GrobidClient(grobid_server=grobid_url)
        client.process(
            "processFulltextDocument",
            str(grobid_input),
            output=str(outdir),
            n=1,
            force=True,
        )
    except Exception as exc:
        return {}, {"source": "grobid_client", "error": str(exc)}

    results: Dict[str, Any] = {}
    for pdf in pdfs:
        pmid = _pmid_from_filename(pdf)
        src = _find_client_output(outdir, pmid)
        target = outdir / f"{pmid}.tei.xml"
        if src and src != target:
            shutil.move(src, target)
        if target.is_file():
            results[pmid] = {"status": "success", "source": "grobid_client"}
        else:
            results[pmid] = {"status": "failed", "source": "grobid_client", "error": "no TEI output"}
    return results, None


def _process_with_requests(
    pdf: Path, grobid_url: str, timeout: int, retries: int, sleep: float
) -> Tuple[Optional[bytes], Optional[str]]:
    url = f"{grobid_url.rstrip('/')}/api/processFulltextDocument"
    for attempt in range(retries):
        try:
            with open(pdf, "rb") as fh:
                files = {"input": (pdf.name, fh, "application/pdf")}
                data = {
                    "consolidateCitations": "0",
                    "includeRawCitations": "0",
                    "teiCoordinates": "0",
                    "segmentSentences": "0",
                }
                r = requests.post(url, files=files, data=data, timeout=timeout)
            if r.status_code == 200 and _looks_tei(r.content):
                return r.content, None
            if r.status_code == 200:
                return None, "GROBID response was not TEI XML"
            return None, f"GROBID returned {r.status_code}"
        except Exception as exc:
            if attempt == retries - 1:
                return None, f"requests failed: {exc}"
            time.sleep(sleep)
    return None, "exhausted retries"


def _setup_client_input(input_dir: Path, outdir: Path) -> Tuple[Path, Path]:
    """GROBID client needs an input directory, so symlink/copy PDFs there."""
    workdir = outdir / "grobid_work"
    workdir.mkdir(parents=True, exist_ok=True)
    grobid_input = workdir / "grobid_input"
    grobid_input.mkdir(parents=True, exist_ok=True)
    for pdf in sorted(input_dir.glob("*.pdf")):
        link = grobid_input / pdf.name
        if not link.exists():
            try:
                link.symlink_to(pdf.resolve())
            except OSError:
                shutil.copy2(pdf, link)
    return workdir, grobid_input


def main() -> None:
    args = parse_args()
    input_dir = Path(args.input_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    pdfs = sorted(input_dir.glob("*.pdf"))
    if not pdfs:
        print("[grobid] No PDFs found to process.", file=sys.stderr)

    # Try grobid-client-python first.
    workdir, grobid_input = _setup_client_input(input_dir, outdir)
    client_results, client_error = _try_grobid_client(pdfs, outdir, grobid_input, args.grobid_url)

    summary: List[Dict[str, Any]] = []
    for pdf in pdfs:
        pmid = _pmid_from_filename(pdf)
        target = outdir / f"{pmid}.tei.xml"

        if pmid in client_results and client_results[pmid]["status"] == "success":
            summary.append({
                "pmid": pmid,
                "status": "success",
                "source": "grobid_client",
                "error": None,
            })
            continue

        # Fallback to direct HTTP requests.
        content, err = _process_with_requests(
            pdf, args.grobid_url, args.timeout, args.retries, args.sleep
        )
        if content is not None:
            with open(target, "wb") as fh:
                fh.write(content)
            summary.append({
                "pmid": pmid,
                "status": "success",
                "source": "requests",
                "error": None,
            })
            print(f"[grobid] Converted PMID {pmid} via requests", file=sys.stderr)
        else:
            client_err = client_results.get(pmid, {}).get("error") or (client_error or {}).get("error")
            summary.append({
                "pmid": pmid,
                "status": "failed",
                "source": "grobid_client",
                "error": err or client_err or "unknown error",
            })
            print(f"[grobid] Failed PMID {pmid}: {err or client_err}", file=sys.stderr)

        time.sleep(args.sleep)

    summary_record = {
        "species": None,
        "domain": None,
        "grobid_url": args.grobid_url,
        "input_count": len(pdfs),
        "success_count": sum(1 for s in summary if s["status"] == "success"),
        "failed_count": sum(1 for s in summary if s["status"] != "success"),
        "results": summary,
    }
    summary_path = outdir / "grobid_summary.json"
    with open(summary_path, "w", encoding="utf-8") as fh:
        json.dump(summary_record, fh, ensure_ascii=False, indent=2)
    print(f"[grobid] Wrote {summary_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
