#!/usr/bin/env python3
"""Fetch WHO FluNet data and preserve every original column.

WHO publishes FluNet data through a public xMart endpoint:
    https://xmart-api-public.who.int/FLUMART/VIW_FNT?%24format=csv

This script filters by ISO year range and by the subtype column(s) relevant
for the detected species, then writes the original wide-format rows (with
DB-friendly column names) into epi_data/<species>/flunet.csv.
"""
import argparse
import csv
import io
import json
import os
import re
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import requests

FLUMART_URL = "https://xmart-api-public.who.int/FLUMART/VIW_FNT"

# Map the species short code to the FluNet subtype column(s) that represent it.
# Only rows where at least one of these columns is positive are downloaded.
SPECIES_COLUMNS: Dict[str, List[str]] = {
    "h1n1pdm": ["AH1N12009"],
    "h3n2": ["AH3"],
    "vic": ["BVIC_2DEL", "BVIC_3DEL", "BVIC_NODEL", "BVIC_DELUNK"],
    "yam": ["BYAM"],
    "h5nx": ["AH5"],
    "h5n1": ["AH5"],
}

# Map WHO API column names to database column names.  The CSV produced by this
# script uses the database names directly, so build_knowledge_db.py can insert
# without further mapping.
FLUNET_COLUMN_MAP: Dict[str, str] = {
    "WHOREGION": "who_region",
    "FLUSEASON": "fluseason",
    "HEMISPHERE": "hemisphere",
    "ITZ": "influenza_transmission_zone",
    "COUNTRY_CODE": "location_code",
    "COUNTRY_AREA_TERRITORY": "country",
    "ISO_WEEKSTARTDATE": "record_date",
    "ISO_YEAR": "year",
    "ISO_WEEK": "week",
    "MMWR_WEEKSTARTDATE": "mmwr_weekstartdate",
    "MMWR_YEAR": "mmwr_year",
    "MMWR_WEEK": "mmwr_week",
    "ORIGIN_SOURCE": "origin_source",
    "SPEC_PROCESSED_NB": "specimens_processed",
    "SPEC_RECEIVED_NB": "specimens_received",
    "AH1N12009": "ah1n1pdm09",
    "AH1": "ah1",
    "AH3": "ah3",
    "AH5": "ah5",
    "AH7N9": "ah7n9",
    "ANOTSUBTYPED": "a_not_subtyped",
    "ANOTSUBTYPABLE": "a_not_subtypable",
    "AOTHER_SUBTYPE": "a_other_subtype",
    "AOTHER_SUBTYPE_DETAILS": "a_other_subtype_details",
    "INF_A": "inf_a",
    "BVIC_2DEL": "b_vic_2del",
    "BVIC_3DEL": "b_vic_3del",
    "BVIC_NODEL": "b_vic_nodel",
    "BVIC_DELUNK": "b_vic_delunk",
    "BYAM": "b_yam",
    "BNOTDETERMINED": "b_not_determined",
    "INF_B": "inf_b",
    "INF_ALL": "inf_all",
    "INF_NEGATIVE": "inf_negative",
    "ILI_ACTIVITY": "ili_activity",
    "ADENO": "adeno",
    "BOCA": "boca",
    "HUMAN_CORONA": "human_corona",
    "METAPNEUMO": "metapneumo",
    "PARAINFLUENZA": "parainfluenza",
    "RHINO": "rhino",
    "RSV_PROCESSED": "rsv_processed",
    "RSV": "rsv",
    "OTHERRESPVIRUS": "other_resp_virus",
    "OTHER_RESPVIRUS_DETAILS": "other_resp_virus_details",
    "LAB_RESULT_COMMENT": "lab_result_comment",
    "WCR_COMMENT": "wcr_comment",
    "ISO2": "iso2",
    "ISOYW": "isoyw",
    "MMWRYW": "mmwryw",
    "PSOURCE_SUBTYPE_INF": "psource_subtype_inf",
    "PSOURCE_PPOS_INF": "psource_ppos_inf",
    "PSOURCE_RSV": "psource_rsv",
}

INTEGER_COLUMNS = {
    "year", "week", "mmwr_year", "mmwr_week", "specimens_processed",
    "specimens_received", "ah1n1pdm09", "ah1", "ah3", "ah5", "ah7n9",
    "a_not_subtyped", "a_not_subtypable", "inf_a", "b_vic_2del", "b_vic_3del",
    "b_vic_nodel", "b_vic_delunk", "b_yam", "b_not_determined", "inf_b",
    "inf_all", "inf_negative", "adeno", "boca", "human_corona", "metapneumo",
    "parainfluenza", "rhino", "rsv_processed", "rsv", "other_resp_virus",
}

DEFAULT_MIN_YEAR = datetime.now().year - 5
DEFAULT_MAX_YEAR = datetime.now().year


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Fetch WHO FluNet data for a detected influenza species.")
    parser.add_argument("--species", required=True, help="Species short code (e.g. h3n2, h1n1pdm, vic, h5nx).")
    parser.add_argument("--min-year", type=int, default=DEFAULT_MIN_YEAR, help="Earliest ISO year to fetch.")
    parser.add_argument("--max-year", type=int, default=DEFAULT_MAX_YEAR, help="Latest ISO year to fetch.")
    parser.add_argument("--outdir", required=True, help="Output directory.")
    parser.add_argument("--timeout", type=int, default=180, help="HTTP timeout in seconds.")
    return parser.parse_args()


def _int(value: Optional[str]) -> Optional[int]:
    if value is None:
        return None
    value = str(value).strip().replace(",", "")
    if value in ("", "NA", "N/A", "None", "null", "NULL"):
        return None
    try:
        return int(float(value))
    except ValueError:
        return None


def _build_filter(min_year: int, max_year: int, columns: List[str]) -> str:
    """Build the OData $filter string for the FluMart API."""
    year_filter = f"ISO_YEAR ge {min_year} and ISO_YEAR le {max_year}"
    col_filters = " or ".join(f"{c} gt 0" for c in columns)
    if len(columns) > 1:
        col_filters = f"({col_filters})"
    return f"{year_filter} and {col_filters}"


def _get_columns_for_species(species: str) -> Optional[List[str]]:
    if species in SPECIES_COLUMNS:
        return SPECIES_COLUMNS[species]
    for key in SPECIES_COLUMNS:
        if key in species:
            return SPECIES_COLUMNS[key]
    return None


def _normalize_date(raw: Optional[str]) -> str:
    if not raw:
        return ""
    m = re.match(r"(\d{4})-(\d{2})-(\d{2})", str(raw))
    if m:
        return m.group(0)
    return str(raw).strip()


def main() -> int:
    args = parse_args()
    species = args.species
    min_year = args.min_year
    max_year = args.max_year
    if min_year > max_year:
        print(f"[fetch_flunet] WARNING: min_year {min_year} > max_year {max_year}, swapping", file=sys.stderr)
        min_year, max_year = max_year, min_year
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    species_dir = outdir / "epi_data" / species
    species_dir.mkdir(parents=True, exist_ok=True)

    columns = _get_columns_for_species(species)
    if not columns:
        print(f"[fetch_flunet] No FluNet column mapping for species '{species}'", file=sys.stderr)
        _write_empty(species_dir, species, min_year, max_year, outdir, f"unknown species {species}")
        return 0

    filter_str = _build_filter(min_year, max_year, columns)
    url = f"{FLUMART_URL}?%24format=csv&%24filter={requests.utils.quote(filter_str)}"

    print(f"[fetch_flunet] {species}: {url}", file=sys.stderr)

    try:
        r = requests.get(url, timeout=args.timeout)
        r.raise_for_status()
    except requests.RequestException as exc:
        print(f"[fetch_flunet] API request failed for {species}: {exc}", file=sys.stderr)
        _write_empty(species_dir, species, min_year, max_year, outdir, str(exc))
        return 0

    reader = csv.DictReader(io.StringIO(r.text))
    if not reader.fieldnames:
        print(f"[fetch_flunet] Empty CSV for {species}", file=sys.stderr)
        _write_empty(species_dir, species, min_year, max_year, outdir, "empty csv")
        return 0

    out_csv = species_dir / "flunet.csv"
    out_fieldnames = [FLUNET_COLUMN_MAP.get(fn, fn) for fn in reader.fieldnames]

    with open(out_csv, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=out_fieldnames)
        writer.writeheader()
        for row in reader:
            mapped = {}
            for api_name, value in row.items():
                db_name = FLUNET_COLUMN_MAP.get(api_name, api_name)
                mapped[db_name] = _clean_value(db_name, value)
            writer.writerow({k: mapped.get(k, "") for k in out_fieldnames})

    _write_summary(outdir, species, min_year, max_year, url, reader.line_num - 1, "success")

    print(f"[fetch_flunet] {species}: wrote {out_csv} ({reader.line_num - 1} rows)", file=sys.stderr)
    return 0


def _clean_value(db_name: str, value: Optional[str]) -> str:
    if value is None:
        return ""
    value = str(value).strip()
    if value in ("", "NA", "N/A", "None", "null", "NULL"):
        return ""
    if db_name in ("record_date", "mmwr_weekstartdate"):
        return _normalize_date(value)
    if db_name in INTEGER_COLUMNS:
        v = _int(value)
        return str(v) if v is not None else ""
    return value


def _write_empty(species_dir: Path, species: str, min_year: int, max_year: int, outdir: Path, message: str) -> None:
    (species_dir).mkdir(parents=True, exist_ok=True)
    out_csv = species_dir / "flunet.csv"
    with open(out_csv, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(FLUNET_COLUMN_MAP.values()))
        writer.writeheader()
    _write_summary(outdir, species, min_year, max_year, "", 0, f"empty: {message}")


def _write_summary(outdir: Path, species: str, min_year: int, max_year: int, url: str, rows: int, status: str) -> None:
    summary_path = outdir / "flunet_search_summary.tsv"
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
            "source": "WHO FluNet",
            "status": status,
        })


if __name__ == "__main__":
    sys.exit(main())
