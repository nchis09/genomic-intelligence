"""
PubMed ID Verification Tool
============================

Validates PMIDs cited in curated CSV files against the PubMed E-utilities API.
This tool was created after discovering that LLM-generated PMIDs were hallucinated
(wrong digits, non-existent papers, wrong author attributions).

Usage:
    python pubmed_verifier.py <csv_file> [--column literature_refs] [--output report.json]

The tool:
1. Extracts all PMID:XXXXX patterns from the specified column
2. Queries PubMed E-utilities esummary for each PMID
3. Reports verified vs unverified PMIDs with citation details
4. Exits with code 1 if any PMIDs are unverified (for CI/CD integration)

This is a CRITICAL quality control tool. No curated CSV should be published
without running this verification first.
"""

import re
import csv
import json
import time
import argparse
import sys
from typing import List, Dict, Tuple, Optional
import urllib.request
import urllib.error

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
ESUMMARY_URL = f"{EUTILS_BASE}/esummary.fcgi"
ESEARCH_URL = f"{EUTILS_BASE}/esearch.fcgi"

PMID_PATTERN = re.compile(r"PMID:(\d+)")


def extract_pmids_from_text(text: str) -> List[str]:
    """Extract all PMID:XXXXX patterns from a text string."""
    if not text or text == "UNVERIFIED":
        return []
    return PMID_PATTERN.findall(text)


def extract_pmids_from_csv(csv_path: str, column_name: str = "literature_refs") -> Dict[str, List[str]]:
    """Extract PMIDs from a CSV file's specified column.
    
    Returns:
        Dict mapping row_id -> list of PMID strings
    """
    pmids_by_row = {}
    with open(csv_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        id_col = None
        for col in reader.fieldnames:
            if col.lower() in ("association_id", "outbreak_id", "lineage_id", "id"):
                id_col = col
                break
        if not id_col:
            id_col = reader.fieldnames[0]
        
        for row in reader:
            row_id = row.get(id_col, "unknown")
            refs_text = row.get(column_name, "")
            pmids = extract_pmids_from_text(refs_text)
            if pmids:
                pmids_by_row[row_id] = pmids
    return pmids_by_row


def verify_pmid(pmid: str, retries: int = 3) -> Optional[Dict]:
    """Verify a single PMID against PubMed E-utilities.
    
    Returns:
        Dict with citation details if verified, None if not found.
    """
    params = (
        f"?db=pubmed&id={pmid}&retmode=json"
    )
    url = ESUMMARY_URL + params
    
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url)
            req.add_header("User-Agent", "GenomicIntelligence/1.0 PMID Verifier")
            with urllib.request.urlopen(req, timeout=15) as resp:
                data = json.loads(resp.read().decode("utf-8"))
            
            result = data.get("result", {})
            if pmid in result and "title" in result[pmid]:
                entry = result[pmid]
                authors = entry.get("authors", [])
                first_author = authors[0]["name"] if authors else "Unknown"
                return {
                    "pmid": pmid,
                    "title": entry.get("title", ""),
                    "first_author": first_author,
                    "journal": entry.get("fulljournalname", entry.get("source", "")),
                    "pubdate": entry.get("pubdate", ""),
                    "verified": True,
                }
            elif "error" in result.get(pmid, {}):
                return None
            else:
                return None
        except urllib.error.URLError as e:
            if attempt < retries - 1:
                time.sleep(2 ** attempt)
                continue
            return {"pmid": pmid, "verified": False, "error": str(e)}
        except Exception as e:
            if attempt < retries - 1:
                time.sleep(2 ** attempt)
                continue
            return {"pmid": pmid, "verified": False, "error": str(e)}
    
    return {"pmid": pmid, "verified": False, "error": "Max retries exceeded"}


def verify_all_pmids(csv_path: str, column_name: str = "literature_refs") -> Dict:
    """Verify all PMIDs in a CSV file.
    
    Returns:
        Dict with verification report
    """
    pmids_by_row = extract_pmids_from_csv(csv_path, column_name)
    
    all_pmids = set()
    for pmids in pmids_by_row.values():
        all_pmids.update(pmids)
    
    print(f"Found {len(all_pmids)} unique PMIDs across {len(pmids_by_row)} rows")
    print(f"Verifying against PubMed E-utilities...")
    
    results = {}
    verified_count = 0
    unverified_count = 0
    unverified_pmids = []
    
    for i, pmid in enumerate(sorted(all_pmids)):
        if pmid in results:
            continue
        
        info = verify_pmid(pmid)
        if info and info.get("verified"):
            results[pmid] = info
            verified_count += 1
            print(f"  [{i+1}/{len(all_pmids)}] PMID:{pmid} ✓ — {info['first_author']} et al., {info['journal']} ({info['pubdate']})")
        else:
            results[pmid] = {"pmid": pmid, "verified": False}
            unverified_count += 1
            unverified_pmids.append(pmid)
            print(f"  [{i+1}/{len(all_pmids)}] PMID:{pmid} ✗ — NOT FOUND IN PUBMED")
        
        time.sleep(0.34)
    
    report = {
        "csv_file": csv_path,
        "column_checked": column_name,
        "total_unique_pmids": len(all_pmids),
        "verified": verified_count,
        "unverified": unverified_count,
        "unverified_pmids": unverified_pmids,
        "results": results,
        "rows_with_pmids": {k: v for k, v in pmids_by_row.items()},
    }
    
    return report


def main():
    parser = argparse.ArgumentParser(
        description="Verify PMIDs in a curated CSV file against PubMed"
    )
    parser.add_argument("csv_file", help="Path to the CSV file to verify")
    parser.add_argument(
        "--column", default="literature_refs",
        help="Column name containing PMID references (default: literature_refs)"
    )
    parser.add_argument(
        "--output", default=None,
        help="Output JSON report path (default: print to stdout)"
    )
    args = parser.parse_args()
    
    report = verify_all_pmids(args.csv_file, args.column)
    
    print(f"\n{'='*60}")
    print(f"PMID VERIFICATION REPORT")
    print(f"{'='*60}")
    print(f"File: {args.csv_file}")
    print(f"Column: {args.column}")
    print(f"Total unique PMIDs: {report['total_unique_pmids']}")
    print(f"Verified: {report['verified']}")
    print(f"Unverified: {report['unverified']}")
    
    if report["unverified_pmids"]:
        print(f"\nUNVERIFIED PMIDs:")
        for pmid in report["unverified_pmids"]:
            rows = [r for r, ps in report["rows_with_pmids"].items() if pmid in ps]
            print(f"  PMID:{pmid} — used in rows: {', '.join(rows)}")
    
    if args.output:
        with open(args.output, "w", encoding="utf-8") as f:
            json.dump(report, f, indent=2)
        print(f"\nReport saved to: {args.output}")
    
    if report["unverified"] > 0:
        sys.exit(1)
    else:
        print(f"\nAll PMIDs verified successfully!")
        sys.exit(0)


if __name__ == "__main__":
    main()
