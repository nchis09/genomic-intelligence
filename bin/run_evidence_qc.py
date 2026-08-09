#!/usr/bin/env python3
"""Validate per-PMID evidence extraction JSONs and split into clean/failed.

Uses Pandera for schema validation of the flattened extraction DataFrame and
computes a per-paper QC score (fraction of expected fields that are present
and type-valid). Papers whose QC score is below the configured threshold or
that fail schema validation are written to the failed/ directory with a
qc_failure annotation; passing papers are copied to clean/.
"""

import argparse
import json
import sys
from pathlib import Path

import pandas as pd
import pandera as pa
import yaml
from pandera import Check, Column, DataFrameSchema


TOP_STATUSES = {"success", "no_evidence", "unknown"}
CONFIDENCES = {"high", "medium", "low"}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Evidence extraction QC")
    p.add_argument("--input-dir", required=True, help="Directory containing *.json extraction files")
    p.add_argument("--outdir", required=True, help="Output directory for clean/failed/qc_report")
    p.add_argument("--species", required=True, help="Species being processed")
    p.add_argument("--domain", required=True, help="Domain being processed")
    p.add_argument("--templates-yml", required=True, help="Path to evidence_templates.yml")
    p.add_argument("--min-completeness", type=float, default=0.0,
                   help="Minimum QC score (fraction of expected fields) for a paper to be clean")
    return p.parse_args()


def load_templates(templates_yml: str, domain: str) -> tuple[dict, list]:
    data = yaml.safe_load(open(templates_yml, encoding="utf-8"))
    fields = data.get(domain, [])
    if not fields:
        raise ValueError(f"No template fields for domain '{domain}' in {templates_yml}")
    return {f["field"]: f for f in fields}, [f["field"] for f in fields]


def build_pandera_schema(expected_fields: list) -> DataFrameSchema:
    field_set = set(expected_fields)
    return DataFrameSchema(
        {
            "pmid": Column(pa.String, nullable=False),
            "field": Column(pa.String, Check.isin(field_set), nullable=False),
            "present": Column(pa.Bool, nullable=False),
            "value": Column(object, nullable=True),
            "quote": Column(object, nullable=True),
            "confidence": Column(pa.String, Check.isin(CONFIDENCES), nullable=True),
        },
        strict="filter",
        coerce=True,
    )


def python_type_matches(value, ftype: str, present: bool) -> bool:
    if not present:
        return value is None
    if value is None:
        return False
    if ftype == "string":
        return isinstance(value, str)
    if ftype == "number":
        return isinstance(value, (int, float)) and not isinstance(value, bool)
    if ftype == "boolean":
        return isinstance(value, bool)
    if ftype == "list":
        return isinstance(value, list)
    return True


def validate_paper(path: Path, schema: DataFrameSchema, field_map: dict, expected_fields: list):
    failures: list = []
    raw_text = path.read_text(encoding="utf-8", errors="replace")

    try:
        data = json.loads(raw_text)
    except json.JSONDecodeError as e:
        failures.append(f"JSON decode error: {e}")
        return None, failures, 0.0

    pmid = data.get("pmid")
    species = data.get("species")
    domain = data.get("domain")
    status = data.get("status", "unknown")

    if not isinstance(pmid, str) or not pmid:
        failures.append("pmid is missing or not a non-empty string")
    if not isinstance(species, str) or not species:
        failures.append("species is missing or not a non-empty string")
    if not isinstance(domain, str) or not domain:
        failures.append("domain is missing or not a non-empty string")
    if status not in TOP_STATUSES:
        failures.append(f"status '{status}' not in {TOP_STATUSES}")
    if not isinstance(data.get("extraction"), list):
        failures.append("extraction must be a list")
    if not isinstance(data.get("not_found"), list):
        failures.append("not_found must be a list")

    extraction = data.get("extraction", [])
    not_found = set(data.get("not_found", []))
    fields_in_extraction: set = set()
    present_valid = 0
    rows = []

    for item in extraction:
        if not isinstance(item, dict):
            failures.append("extraction contains a non-object item")
            continue

        field = item.get("field")
        fields_in_extraction.add(field)
        fdef = field_map.get(field, {})
        ftype = fdef.get("type", "string")
        present = item.get("present")
        value = item.get("value")

        if isinstance(present, bool) and present:
            if value is None:
                failures.append(f"{field}: present=True but value is null")
            elif not python_type_matches(value, ftype, present):
                failures.append(
                    f"{field}: value type {type(value).__name__} does not match template type '{ftype}'"
                )
            else:
                present_valid += 1
        elif present is False:
            if value is not None:
                failures.append(f"{field}: present=False but value is not null")
        else:
            failures.append(f"{field}: 'present' must be a boolean")

        rows.append({
            "pmid": pmid if isinstance(pmid, str) else str(path.stem),
            "field": field if isinstance(field, str) else None,
            "present": bool(present) if isinstance(present, bool) else False,
            "value": value,
            "quote": item.get("quote"),
            "confidence": item.get("confidence"),
        })

    missing = set(expected_fields) - fields_in_extraction
    extra = fields_in_extraction - set(expected_fields)
    if missing:
        failures.append(f"Missing extraction entries for fields: {sorted(missing)}")
    if extra:
        failures.append(f"Unexpected fields in extraction: {sorted(extra)}")

    for field in expected_fields:
        item = next((i for i in extraction if isinstance(i, dict) and i.get("field") == field), None)
        present = item.get("present") if isinstance(item, dict) else None
        in_not_found = field in not_found
        if in_not_found and present is True:
            failures.append(f"{field} listed in not_found but has present=True")
        elif not in_not_found and (present is False or item is None):
            failures.append(f"{field} not in not_found but present=False or missing")

    if rows:
        try:
            df = pd.DataFrame(rows)
            schema.validate(df, lazy=True)
        except pa.errors.SchemaErrors as e:
            for err in e.schema_errors:
                failures.append(f"Pandera schema error: {err}")
        except Exception as e:
            failures.append(f"Pandera validation error: {e}")

    qc_score = present_valid / len(expected_fields) if expected_fields else 0.0
    return data, failures, qc_score


def write_clean_or_failed(outdir: Path, data, score: float, failures: list, min_score: float) -> bool:
    data = dict(data) if data else {"pmid": None}
    data["qc_score"] = round(score, 4)

    pmid = data.get("pmid")
    if not pmid:
        pmid = "unknown"

    if score >= min_score and not failures:
        clean_dir = outdir / "clean"
        clean_dir.mkdir(exist_ok=True)
        (clean_dir / f"{pmid}.json").write_text(json.dumps(data, indent=2), encoding="utf-8")
        return True
    else:
        failed_dir = outdir / "failed"
        failed_dir.mkdir(exist_ok=True)
        data["qc_failure"] = failures
        (failed_dir / f"{pmid}.json").write_text(json.dumps(data, indent=2), encoding="utf-8")
        return False


def build_report(clean: list, failed: list, expected_fields: list, species: str, domain: str,
                 min_completeness: float) -> dict:
    all_rows = []
    for data in clean:
        for item in data.get("extraction", []):
            if not isinstance(item, dict):
                continue
            row = dict(item)
            row["pmid"] = data.get("pmid")
            all_rows.append(row)

    if all_rows:
        df = pd.DataFrame(all_rows)
        conf_counts = df["confidence"].value_counts().to_dict() if "confidence" in df.columns else {}
        present_rows = df[df["present"] == True] if "present" in df.columns else pd.DataFrame()
        quote_total = len(present_rows)
        quote_present = int(present_rows["quote"].notna().sum()) if "quote" in present_rows.columns else 0
        field_coverage = {f: int((df["field"] == f).any()) if "field" in df.columns else 0 for f in expected_fields}
        summary = {
            "confidence_distribution": conf_counts,
            "quote_presence_rate": round(quote_present / quote_total, 4) if quote_total else None,
            "field_coverage": field_coverage,
            "duplicate_pmids": int(df["pmid"].duplicated().sum()) if "pmid" in df.columns else 0,
        }
    else:
        summary = {}

    per_paper = []
    for data in clean:
        per_paper.append({
            "pmid": data.get("pmid"),
            "qc_score": data.get("qc_score"),
            "status": "clean",
            "qc_failure": [],
        })
    for data in failed:
        per_paper.append({
            "pmid": data.get("pmid"),
            "qc_score": data.get("qc_score"),
            "status": "failed",
            "qc_failure": data.get("qc_failure", []),
        })

    return {
        "species": species,
        "domain": domain,
        "total_papers": len(clean) + len(failed),
        "clean_count": len(clean),
        "failed_count": len(failed),
        "min_completeness": min_completeness,
        "expected_fields": expected_fields,
        "summary": summary,
        "per_paper": per_paper,
    }


def make_html(report: dict, clean: list, failed: list) -> str:
    rows = []
    for p in report["per_paper"]:
        rows.append(
            "<tr>"
            f"<td>{p.get('pmid')}</td>"
            f"<td>{p.get('status')}</td>"
            f"<td>{p.get('qc_score')}</td>"
            f"<td>{' | '.join(p.get('qc_failure', []))}</td>"
            "</tr>"
        )
    table = (
        "<table border='1' cellspacing='0' cellpadding='4'>"
        "<tr><th>PMID</th><th>Status</th><th>QC Score</th><th>Failures</th></tr>"
        + "".join(rows)
        + "</table>"
    )
    summary = json.dumps(report.get("summary", {}), indent=2)
    return f"""<!DOCTYPE html>
<html>
<head>
  <title>Evidence QC - {report['species']} / {report['domain']}</title>
  <style>
    body {{ font-family: sans-serif; margin: 2em; }}
    h1, h2 {{ color: #333; }}
    pre {{ background: #f5f5f5; padding: 1em; overflow-x: auto; }}
    table {{ border-collapse: collapse; margin-top: 1em; }}
    th {{ background: #eee; }}
  </style>
</head>
<body>
  <h1>Evidence QC Report — {report['species']} / {report['domain']}</h1>
  <p>
    Total: {report['total_papers']} |
    Clean: {report['clean_count']} |
    Failed: {report['failed_count']} |
    Min completeness: {report['min_completeness']}
  </p>
  <h2>Summary</h2>
  <pre>{summary}</pre>
  <h2>Per-paper</h2>
  {table}
</body>
</html>"""


def main() -> int:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    field_map, expected_fields = load_templates(args.templates_yml, args.domain)
    schema = build_pandera_schema(expected_fields)

    clean: list = []
    failed: list = []

    for path in sorted(Path(args.input_dir).glob("*.json")):
        data, failures, score = validate_paper(path, schema, field_map, expected_fields)
        if data is None:
            data = {"pmid": path.stem}
            score = 0.0
        data["qc_score"] = round(score, 4)

        if write_clean_or_failed(outdir, data, score, failures, args.min_completeness):
            clean.append(data)
        else:
            failed.append(data)

    report = build_report(clean, failed, expected_fields, args.species, args.domain, args.min_completeness)
    (outdir / "qc_report.json").write_text(json.dumps(report, indent=2, default=str), encoding="utf-8")
    (outdir / "qc_report.html").write_text(make_html(report, clean, failed), encoding="utf-8")

    for data in clean:
        print(outdir / "clean" / f"{data['pmid']}.json")

    return 0


if __name__ == "__main__":
    sys.exit(main())
