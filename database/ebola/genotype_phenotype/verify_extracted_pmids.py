#!/usr/bin/env python3
"""
Verify PMIDs in the genotype_phenotype table.

Reads unverified genotype-phenotype associations, checks each PMID against
PubMed E-utilities, and updates the table:
  - verification_status = 'verified' when all PMIDs resolve.
  - verification_status = 'invalid' when at least one PMID does not resolve.
  - flag_reason is updated to show the PubMed title(s) found for review.

Usage:
    python verify_extracted_pmids.py [--batch-size 50]
"""

from __future__ import annotations

import argparse
import logging
import re
import sys
import time
from datetime import date
from pathlib import Path
from typing import Optional

import psycopg2

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from config import DB_URL
from database.ebola.genotype_phenotype.pubmed_verifier import verify_pmid

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s", datefmt="%H:%M:%S")
log = logging.getLogger(__name__)

PMID_PATTERN = re.compile(r"PMID:(\d+)")


def extract_pmids_from_refs(refs: list[str] | str) -> list[str]:
    """Extract PMID strings from a Postgres text array or a plain string."""
    if not refs:
        return []
    if isinstance(refs, list):
        text = " ".join(str(r) for r in refs)
    else:
        text = str(refs)
    return PMID_PATTERN.findall(text)


def verify_extracted_pmids(batch_size: int = 50, dry_run: bool = False) -> dict:
    conn = psycopg2.connect(DB_URL)
    cur = conn.cursor()
    try:
        cur.execute(
            """
            SELECT association_id, literature_refs, notes
            FROM genotype_phenotype
            WHERE verification_status = 'unverified'
            ORDER BY association_id
            LIMIT %s
            """,
            (batch_size,),
        )
        rows = cur.fetchall()
    except Exception as e:
        log.error("Database query failed: %s", e)
        conn.close()
        return {}

    verified_map = {}
    results = {"verified": 0, "invalid": 0, "no_pmids": 0, "pmids_checked": {}, "details": []}

    for association_id, refs, notes in rows:
        pmids = extract_pmids_from_refs(refs)
        if not pmids:
            results["no_pmids"] += 1
            if not dry_run:
                cur.execute(
                    "UPDATE genotype_phenotype SET verification_status = 'invalid', flag_reason = 'No PMIDs found', last_updated = %s WHERE association_id = %s",
                    (date.today(), association_id),
                )
            continue

        all_ok = True
        titles = []
        for pmid in pmids:
            if pmid not in verified_map:
                verified_map[pmid] = verify_pmid(pmid)
                time.sleep(0.34)
            info = verified_map[pmid]
            results["pmids_checked"][pmid] = info
            if info and info.get("verified"):
                titles.append(f"{pmid}: {info['title'][:60]}...")
            else:
                all_ok = False

        status = "verified" if all_ok else "invalid"
        reason = "All PMIDs verified. " + "; ".join(titles) if all_ok else "At least one PMID not found: " + ", ".join(pmids)

        results[status] += 1
        results["details"].append({"association_id": association_id, "status": status, "pmids": pmids, "reason": reason})

        if not dry_run:
            cur.execute(
                """
                UPDATE genotype_phenotype
                SET verification_status = %s,
                    flag_reason = %s,
                    last_updated = %s
                WHERE association_id = %s
                """,
                (status, reason, date.today(), association_id),
            )

    if not dry_run:
        conn.commit()
    else:
        log.info("DRY RUN: no DB updates applied")

    cur.close()
    conn.close()
    return results


def main():
    parser = argparse.ArgumentParser(description="Verify PMIDs in genotype_phenotype table")
    parser.add_argument("--batch-size", type=int, default=50)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    results = verify_extracted_pmids(batch_size=args.batch_size, dry_run=args.dry_run)
    print(f"PMID verification report:")
    print(f"  verified: {results.get('verified', 0)}")
    print(f"  invalid: {results.get('invalid', 0)}")
    print(f"  no_pmids: {results.get('no_pmids', 0)}")
    print(f"  unique PMIDs checked: {len(results.get('pmids_checked', {}))}")


if __name__ == "__main__":
    main()
