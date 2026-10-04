#!/usr/bin/env python3
"""Fetch FAO EMPRES-i avian influenza outbreak data.

FAO exposes EMPRES-i disease events through a public BigQuery endpoint:
    https://api.data.apps.fao.org/api/v2/bigquery

This script queries confirmed avian influenza events for a date range and
writes the original wide-format rows (with DB-friendly column names) into
epi_data/<species>/empresi_avian.csv.
"""
import argparse
import csv
import io
import os
import re
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, Optional

import requests

BASE_URL = "https://api.data.apps.fao.org/api/v2/bigquery"
SQL_URL = (
    "https://data.apps.fao.org/catalog/dataset/"
    "3ff164cd-8b44-46d3-8f88-92b0361c7878/resource/"
    "137a69a0-ad5f-48c3-b927-3a61d2c9a2ce/download/"
    "animal-major-disease-parameterized-query.sql"
)
DISEASE_NAME = "Influenza - Avian"
DIAGNOSIS_STATUS = "Confirmed"
ANIMAL_TYPE = "all"
COUNTRY = "all"

EMPRESI_COLUMN_MAP: Dict[str, str] = {
    "global_id": "global_id",
    "lat": "lat",
    "lon": "lon",
    "locality": "locality",
    "region": "region",
    "location": "location",
    "observation_date": "observation_date",
    "report_date": "report_date",
    "display_date": "display_date",
    "species_overview_list": "species_affected",
    "humans_affected": "humans_affected",
    "humans_deaths": "humans_deaths",
    "diagnosis_source": "diagnosis_source",
    "diagnosis_status": "diagnosis_status",
    "animal_type_list": "animal_type",
    "disease": "disease",
    "country": "country",
}

DEFAULT_MIN_YEAR = datetime.now().year - 5
DEFAULT_MAX_YEAR = datetime.now().year


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Fetch FAO EMPRES-i avian influenza events for a detected avian H5 species."
    )
    parser.add_argument("--species", required=True, help="Species short code (e.g. h5nx, h5n1).")
    parser.add_argument("--min-year", type=int, default=DEFAULT_MIN_YEAR, help="Earliest year to fetch.")
    parser.add_argument("--max-year", type=int, default=DEFAULT_MAX_YEAR, help="Latest year to fetch.")
    parser.add_argument("--outdir", required=True, help="Output directory.")
    parser.add_argument("--timeout", type=int, default=300, help="HTTP timeout in seconds.")
    return parser.parse_args()


def _int(value: Optional[str]) -> Optional[int]:
    if value is None:
        return None
    value = str(value).strip().replace(",", "")
    if value in ("", "NA", "N/A", "None", "null", "NULL", "no_value"):
        return None
    try:
        return int(float(value))
    except ValueError:
        return None


def _float(value: Optional[str]) -> Optional[float]:
    if value is None:
        return None
    value = str(value).strip().replace(",", "")
    if value in ("", "NA", "N/A", "None", "null", "NULL", "no_value"):
        return None
    try:
        return float(value)
    except ValueError:
        return None


def _normalize_date(raw: Optional[str]) -> str:
    if not raw:
        return ""
    m = re.match(r"(\d{4})-(\d{2})-(\d{2})", str(raw))
    if m:
        return m.group(0)
    return str(raw).strip()


def _clean_value(db_name: str, value: Optional[str]) -> str:
    if value is None:
        return ""
    value = str(value).strip()
    if value in ("", "NA", "N/A", "None", "null", "NULL", "no_value"):
        return ""
    if db_name in ("observation_date", "report_date", "display_date"):
        return _normalize_date(value)
    if db_name in ("humans_affected", "humans_deaths"):
        v = _int(value)
        return str(v) if v is not None else ""
    if db_name in ("lat", "lon"):
        v = _float(value)
        return str(v) if v is not None else ""
    return value


def main() -> int:
    args = parse_args()
    species = args.species
    min_year = args.min_year
    max_year = args.max_year
    if min_year > max_year:
        print(f"[fetch_empresi] WARNING: min_year {min_year} > max_year {max_year}, swapping", file=sys.stderr)
        min_year, max_year = max_year, min_year
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    species_dir = outdir / "epi_data" / species
    species_dir.mkdir(parents=True, exist_ok=True)

    start_date = f"{min_year}-01-01"
    end_date = f"{max_year}-12-31"

    params = {
        "sql_url": SQL_URL,
        "start_date": start_date,
        "end_date": end_date,
        "country": COUNTRY,
        "disease": DISEASE_NAME,
        "diagnosis_status": DIAGNOSIS_STATUS,
        "animal_type": ANIMAL_TYPE,
        "download": "true",
    }

    print(f"[fetch_empresi] {species}: {BASE_URL} (start={start_date} end={end_date})", file=sys.stderr)

    try:
        r = requests.get(BASE_URL, params=params, timeout=args.timeout, stream=True)
        r.raise_for_status()
    except requests.RequestException as exc:
        print(f"[fetch_empresi] API request failed for {species}: {exc}", file=sys.stderr)
        _write_empty(species_dir, outdir, species, min_year, max_year, "", str(exc))
        return 0

    out_csv = species_dir / "empresi_avian.csv"
    reader = csv.DictReader(io.StringIO(r.text))
    if not reader.fieldnames:
        print(f"[fetch_empresi] Empty CSV for {species}", file=sys.stderr)
        _write_empty(species_dir, outdir, species, min_year, max_year, r.url, "empty csv")
        return 0

    out_fieldnames = [EMPRESI_COLUMN_MAP.get(fn, fn) for fn in reader.fieldnames]

    with open(out_csv, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=out_fieldnames)
        writer.writeheader()
        for row in reader:
            mapped = {}
            for api_name, value in row.items():
                db_name = EMPRESI_COLUMN_MAP.get(api_name, api_name)
                mapped[db_name] = _clean_value(db_name, value)
            writer.writerow({k: mapped.get(k, "") for k in out_fieldnames})

    rows = reader.line_num - 1
    _write_summary(outdir, species, min_year, max_year, r.url, rows, "success")

    print(f"[fetch_empresi] {species}: wrote {out_csv} ({rows} rows)", file=sys.stderr)
    return 0


def _write_empty(species_dir: Path, outdir: Path, species: str, min_year: int, max_year: int, url: str, message: str) -> None:
    out_csv = species_dir / "empresi_avian.csv"
    with open(out_csv, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(EMPRESI_COLUMN_MAP.values()))
        writer.writeheader()
    _write_summary(outdir, species, min_year, max_year, url, 0, f"empty: {message}")


def _write_summary(outdir: Path, species: str, min_year: int, max_year: int, url: str, rows: int, status: str) -> None:
    summary_path = outdir / "empresi_search_summary.tsv"
    with open(summary_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=[
            "species", "min_year", "max_year", "url", "rows", "source", "status",
        ], delimiter="\t")
        writer.writeheader()
        writer.writerow({
            "species": species,
            "min_year": min_year,
            "max_year": max_year,
            "url": url,
            "rows": rows,
            "source": "FAO EMPRES-i",
            "status": status,
        })


if __name__ == "__main__":
    sys.exit(main())
