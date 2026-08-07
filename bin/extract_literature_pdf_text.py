#!/usr/bin/env python3
"""Extract plain text from downloaded PDFs using pymupdf or pdfplumber.

No external server (e.g. GROBID) is required. Output is one .txt file per PMID
plus a pdf_text_summary.json file.
"""
import argparse
import json
import os
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Extract text from PDFs with pymupdf or pdfplumber.")
    parser.add_argument("--input-dir", required=True, help="Directory with PDF files.")
    parser.add_argument("--outdir", required=True, help="Directory for text output and summary.")
    parser.add_argument("--species", required=True, help="Species key for the summary.")
    parser.add_argument("--domain", required=True, help="Domain key for the summary.")
    parser.add_argument("--format", choices=["txt", "md"], default="txt", help="Output format.")
    parser.add_argument("--sleep", type=float, default=0.0, help="Seconds between files.")
    return parser.parse_args()


def _pmid_from_filename(pdf_path: Path) -> str:
    return pdf_path.stem


def _has_pymupdf() -> bool:
    try:
        import fitz  # noqa: F401
        return True
    except Exception:
        return False


def _has_pdfplumber() -> bool:
    try:
        import pdfplumber  # noqa: F401
        return True
    except Exception:
        return False


def _extract_with_pymupdf(pdf_path: Path) -> Tuple[Optional[str], Optional[str]]:
    try:
        import fitz
    except Exception as exc:
        return None, f"pymupdf not available: {exc}"

    try:
        text_parts = []
        with fitz.open(str(pdf_path)) as doc:
            for page_num, page in enumerate(doc, start=1):
                page_text = page.get_text("text")
                if page_text and page_text.strip():
                    text_parts.append(f"--- Page {page_num} ---\n{page_text.strip()}")
        full_text = "\n\n".join(text_parts)
        if not full_text.strip():
            return None, "pymupdf produced empty text"
        return full_text, None
    except Exception as exc:
        return None, f"pymupdf extraction failed: {exc}"


def _extract_with_pdfplumber(pdf_path: Path) -> Tuple[Optional[str], Optional[str]]:
    try:
        import pdfplumber
    except Exception as exc:
        return None, f"pdfplumber not available: {exc}"

    try:
        text_parts = []
        with pdfplumber.open(str(pdf_path)) as pdf:
            for page_num, page in enumerate(pdf.pages, start=1):
                page_text = page.extract_text()
                if page_text and page_text.strip():
                    text_parts.append(f"--- Page {page_num} ---\n{page_text.strip()}")
        full_text = "\n\n".join(text_parts)
        if not full_text.strip():
            return None, "pdfplumber produced empty text"
        return full_text, None
    except Exception as exc:
        return None, f"pdfplumber extraction failed: {exc}"


def _extract_text(pdf_path: Path) -> Tuple[Optional[str], Optional[str], str]:
    """Try pymupdf first, then pdfplumber. Return (text, error, source)."""
    if _has_pymupdf():
        text, err = _extract_with_pymupdf(pdf_path)
        if text is not None:
            return text, None, "pymupdf"
        if _has_pdfplumber():
            text2, err2 = _extract_with_pdfplumber(pdf_path)
            if text2 is not None:
                return text2, None, "pdfplumber"
            return None, f"pymupdf: {err}; pdfplumber: {err2}", "pdfplumber"
        return None, f"pymupdf: {err}; no pdfplumber", "pymupdf"

    if _has_pdfplumber():
        text, err = _extract_with_pdfplumber(pdf_path)
        if text is not None:
            return text, None, "pdfplumber"
        return None, f"pdfplumber: {err}", "pdfplumber"

    return None, " neither pymupdf nor pdfplumber is installed", "none"


def _pdf_page_count(pdf_path: Path) -> int:
    try:
        import fitz
        with fitz.open(str(pdf_path)) as doc:
            return len(doc)
    except Exception:
        pass
    try:
        import pdfplumber
        with pdfplumber.open(str(pdf_path)) as pdf:
            return len(pdf.pages)
    except Exception:
        return 0


def main() -> None:
    args = parse_args()
    input_dir = Path(args.input_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    pdfs = sorted(input_dir.glob("*.pdf"))
    if not pdfs:
        print("[extract_text] No PDFs found to process.", file=sys.stderr)

    summary: List[Dict[str, Any]] = []
    for pdf in pdfs:
        pmid = _pmid_from_filename(pdf)
        target = outdir / f"{pmid}.{args.format}"

        text, err, source = _extract_text(pdf)
        pages = _pdf_page_count(pdf) if text is not None else 0
        chars = len(text) if text is not None else 0

        if text is not None:
            with open(target, "w", encoding="utf-8") as fh:
                fh.write(text)
            print(f"[extract_text] Extracted PMID {pmid} ({chars} chars, {pages} pages) via {source}", file=sys.stderr)
            summary.append({
                "pmid": pmid,
                "status": "success",
                "source": source,
                "error": None,
                "pages": pages,
                "chars": chars,
            })
        else:
            print(f"[extract_text] Failed PMID {pmid}: {err}", file=sys.stderr)
            summary.append({
                "pmid": pmid,
                "status": "failed",
                "source": source,
                "error": err,
                "pages": 0,
                "chars": 0,
            })

        if args.sleep:
            time.sleep(args.sleep)

    summary_record = {
        "species": args.species,
        "domain": args.domain,
        "input_count": len(pdfs),
        "success_count": sum(1 for s in summary if s["status"] == "success"),
        "failed_count": sum(1 for s in summary if s["status"] != "success"),
        "results": summary,
    }
    summary_path = outdir / "pdf_text_summary.json"
    with open(summary_path, "w", encoding="utf-8") as fh:
        json.dump(summary_record, fh, ensure_ascii=False, indent=2)
    print(f"[extract_text] Wrote {summary_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
