#!/usr/bin/env python3
"""Extract concrete, quotable evidence claims from paper text files using Ollama.

Reads plain-text outputs of LITERATURE_TEXT and the matching PubMed metadata
JSONs, calls a local Ollama model, and emits a single TSV of extracted claims
(evidence_extracted.tsv) plus an extraction log (extraction_log.json).
"""
import argparse
import csv
import datetime
import io
import json
import os
import re
import sys
import urllib.request
from pathlib import Path
from typing import Any, Dict, List, Optional

import yaml

BASE_COLUMNS = [
    "pmid",
    "pmcid",
    "species",
    "domain",
    "topic",
    "claim_type",
    "finding",
    "quote",
    "evidence_level",
    "population",
    "geography",
    "source_file",
]


def expected_columns(details_fields: List[Dict[str, Any]]) -> List[str]:
    """TSV header for the current domain: base columns + domain fields + details JSON."""
    detail_names = [f["field"] for f in details_fields if isinstance(f, dict) and f.get("field")]
    return BASE_COLUMNS + detail_names + ["details"]

VALID_TOPICS = {
    "vaccine", "therapeutic", "diagnostic", "surveillance", "transmission",
    "clinical", "reservoir", "policy", "unknown",
}
VALID_CLAIM_TYPES = {
    "efficacy", "safety", "coverage", "availability", "resistance",
    "pharmacokinetics", "unknown",
}
VALID_EVIDENCE_LEVELS = {
    "RCT", "observational", "case_report", "in_silico", "review",
    "expert_opinion", "abstract_only", "unknown",
}

THINKING_FAMILIES = {"qwen3", "deepseek", "deepseek-r1", "gpt-oss"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Extract evidence claims from paper text with Ollama.")
    parser.add_argument("--input-dir", required=True, help="Directory with .txt paper files.")
    parser.add_argument("--metadata-dir", required=True, help="Directory with matching .json metadata files.")
    parser.add_argument("--outdir", required=True, help="Directory for evidence_extracted.tsv and extraction_log.json.")
    parser.add_argument("--species", required=True, help="Species key.")
    parser.add_argument("--domain", required=True, help="Domain key.")
    parser.add_argument("--templates-yml", required=True, help="Path to evidence_templates.yml (used for domain context).")
    parser.add_argument("--ollama-host", default=os.environ.get("PG_OLLAMA_HOST", "http://localhost:11434"), help="Ollama base URL.")
    parser.add_argument("--ollama-model", default=os.environ.get("PG_OLLAMA_MODEL", ""), help="Ollama model name.")
    parser.add_argument("--n-ctx", type=int, default=4096, help="Ollama context size (token budget).")
    parser.add_argument("--temperature", type=float, default=0.1, help="Ollama temperature.")
    parser.add_argument("--timeout", type=int, default=240, help="HTTP timeout in seconds per call.")
    return parser.parse_args()


def load_templates(path: Path) -> Dict[str, Any]:
    with open(path, "r", encoding="utf-8") as fh:
        return yaml.safe_load(fh) or {}


def _is_thinking_model(m: Dict[str, Any]) -> bool:
    details = m.get("details") or {}
    families = set(details.get("families") or [details.get("family", "")])
    return bool(families & THINKING_FAMILIES) or any(f in (m.get("name") or "").lower() for f in THINKING_FAMILIES)


def resolve_ollama_model(host: str, model_arg: str) -> str:
    if model_arg and model_arg.strip():
        return model_arg.strip()
    try:
        req = urllib.request.Request(f"{host}/api/tags", method="GET")
        with urllib.request.urlopen(req, timeout=10) as resp:
            data = json.loads(resp.read().decode("utf-8"))
        models = data.get("models", [])
        if not models:
            raise RuntimeError("No Ollama models found. Run 'ollama pull <model>'.")

        # 1) explicit llama3 if present
        for m in models:
            name = (m.get("name") or m.get("model") or "").lower()
            if "llama3" in name:
                return m["name"]

        # 2) prefer a non-thinking model, ideally qwen2.5 (prefer 14b over 7b if available)
        non_thinking = [m for m in models if not _is_thinking_model(m)]
        qwen = [m for m in non_thinking if (m.get("name") or "").startswith("qwen2.5")]
        preferred = [m for m in qwen if ":14b" in (m.get("name") or "")]
        if (preferred):
            return preferred[0]["name"]
        if qwen:
            return qwen[0]["name"]
        if non_thinking:
            return non_thinking[0]["name"]

        # 3) fall back to the first model and warn that it may think
        return models[0]["name"]
    except Exception as exc:
        raise RuntimeError(f"Could not resolve Ollama model at {host}: {exc}")


def _call_ollama_generate(host: str, model: str, prompt: str, n_ctx: int, temperature: float, timeout: int) -> str:
    body = {
        "model": model,
        "prompt": prompt,
        "stream": False,
        "options": {
            "temperature": temperature,
            "num_ctx": n_ctx,
            "num_predict": 8192,
        },
    }
    data = json.dumps(body).encode("utf-8")
    req = urllib.request.Request(
        f"{host}/api/generate",
        data=data,
        headers={"Content-Type": "application/json"},
        method="POST",
    )
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            raw = json.loads(resp.read().decode("utf-8"))
    except urllib.error.HTTPError as exc:
        raise RuntimeError(f"Ollama HTTP {exc.code}: {exc.read().decode('utf-8', errors='ignore')}")
    except Exception as exc:
        raise RuntimeError(f"Ollama request failed: {exc}")

    if not raw.get("response"):
        raise RuntimeError("Ollama returned an empty response")
    return raw["response"].strip()


BASE_FIELDS = set(BASE_COLUMNS) | {"details"}


def _details_fields(templates: Dict[str, Any], domain: str) -> List[Dict[str, Any]]:
    """Return per-domain detail fields not already captured at the top level."""
    fields: List[Dict[str, Any]] = []
    counts = templates.get("counts", {})
    n_field = counts.get("n") if isinstance(counts, dict) else None
    if n_field:
        fields.append({**n_field, "field": "n"})
    domain_fields = templates.get(domain, []) or []
    for f in domain_fields:
        if isinstance(f, dict) and f.get("field") not in BASE_FIELDS:
            fields.append(f)
    return fields


def _details_prompt_block(fields: List[Dict[str, Any]]) -> str:
    if not fields:
        return "For this domain, the 'details' object is optional and may be left empty if no domain-specific data is present.\n"
    lines = ["For this domain, the 'details' object must include these keys when stated in the text (use null or '' if not available):\n"]
    for f in fields:
        key = f.get("field", "")
        question = f.get("question", "")
        ftype = f.get("type", "string")
        lines.append(f'  "{key}" ({ftype}): {question}\n')
    return "".join(lines)


def _build_prompt(domain: str, details_fields: List[Dict[str, Any]], metadata: Dict[str, Any], text: str, n_ctx: int) -> str:
    details_block = _details_prompt_block(details_fields)

    instruction = (
        "You are an exact biomedical evidence extractor. Read the paper text and return\n"
        "a JSON list of 1 to 5 relevant, explicitly stated statements. Each item must be an object with exactly\n"
        "these top-level fields and no others:\n\n"
        "  pmid, pmcid, species, domain, topic, claim_type,\n"
        "  finding, quote, evidence_level, population, geography, details\n\n"
        "Rules:\n"
        "- topic must be one of: vaccine, therapeutic, diagnostic, surveillance, transmission, clinical, reservoir, policy, unknown.\n"
        "- claim_type must be one of: efficacy, safety, coverage, availability, resistance, pharmacokinetics, unknown.\n"
        "- evidence_level must be one of: RCT, observational, case_report, in_silico, review, expert_opinion, abstract_only, unknown.\n"
        "- product_name: a named product, drug, test, vaccine, or intervention if mentioned; otherwise ''.\n"
        "- finding: a one-sentence summary of the statement. This field must never be empty.\n"
        "- quote: the exact sentence or clause from the text that supports the finding. This field is required; do not leave it empty if a source sentence exists.\n"
        "- details field types: number = exact numeric value; string = exact text; boolean = true/false; list = semicolon-separated items.\n"
        "- population: e.g. humans, nonhuman_primate, cell_culture, or ''.\n"
        "- geography: country/region mentioned, or ''.\n"
        "- details: a nested object for the domain of interest. See the domain-specific list below.\n"
        "- Use only explicit statements from the text. Do not invent, infer, or synthesise facts.\n"
        "- The 'domain' is the search keyword that found the paper; extract claims from any relevant topic.\n"
        "- Return the exact PMID, PMCID, species and domain values given below.\n"
        "- Return valid JSON only. No markdown, no explanation, no preamble.\n\n"
        + details_block +
        "\n"
        "Example (a diagnostic test):\n"
        "[\n"
        "  {\n"
        "    \"pmid\": \"87654321\",\n"
        "    \"pmcid\": \"\",\n"
        "    \"species\": \"ebov\",\n"
        "    \"domain\": \"diagnostic\",\n"
        "    \"topic\": \"diagnostic\",\n"
        "    \"claim_type\": \"unknown\",\n"
        "    \"finding\": \"OraQuick demonstrated 84.0% sensitivity in archived EVD patient venous whole-blood samples.\",\n"
        "    \"quote\": \"OraQuick Ebola demonstrated clinical sensitivity of 84.0% in archived EVD patient venous whole-blood samples.\",\n"
        "    \"evidence_level\": \"observational\",\n"
        "    \"population\": \"humans\",\n"
        "    \"geography\": \"\",\n"
        "    \"details\": {\n"
        "      \"test_name\": \"OraQuick Ebola Rapid Antigen Test\",\n"
        "      \"test_type\": \"rapid antigen RDT\",\n"
        "      \"sample_type\": \"venous whole-blood\",\n"
        "      \"sensitivity\": \"84.0%\",\n"
        "      \"specificity\": \"98.0%\",\n"
        "      \"n\": 482\n"
        "    }\n"
        "  }\n"
        "]\n\n"
    )

    header = (
        f"pmid: {metadata.get('pmid', '')}\n"
        f"pmcid: {metadata.get('pmcid', '')}\n"
        f"species: {metadata.get('species', 'unknown')}\n"
        f"domain: {domain}\n"
        f"title: {metadata.get('title', '')}\n"
        f"abstract: {metadata.get('abstract', '')}\n\n"
        "Text to extract from:\n"
    )
    footer = "\n\nJSON list of claims:\n"

    # Approx 1 token ~ 4 characters. Leave a 1k-token margin for the model answer.
    budget = (n_ctx - 1024) * 4
    current = len(instruction) + len(header) + len(footer)
    if current + len(text) > budget:
        keep = max(budget - current, 0)
        text = text[:keep]

    return instruction + header + text + footer


def _parse_llm_json(raw: str) -> List[Dict[str, Any]]:
    """Parse a JSON list or {'claims': [...]} from the Ollama response."""
    raw = raw.strip()
    if raw.startswith("```"):
        raw = re.sub(r"^```(?:json)?\s*", "", raw)
        raw = re.sub(r"\s*```$", "", raw)
    raw = raw.strip()
    if not raw:
        return []
    try:
        parsed = json.loads(raw)
    except Exception as exc:
        # Try to extract the first JSON array that looks like a list of objects
        m = re.search(r"(\[\s*\{.*?\}\s*\])", raw, re.DOTALL)
        if m:
            try:
                parsed = json.loads(m.group(1))
            except Exception:
                raise RuntimeError(f"JSON parse failed: {exc}")
        else:
            raise RuntimeError(f"JSON parse failed: {exc}")

    if isinstance(parsed, list):
        return parsed
    if isinstance(parsed, dict):
        return parsed.get("claims") or parsed.get("extraction") or parsed.get("results") or []
    return []


def _normalise_value(v: Any) -> str:
    if v is None:
        return ""
    if isinstance(v, list):
        return "; ".join(str(x) for x in v)
    return str(v).strip()


def _clean_row(row: Dict[str, Any], pmid: str, pmcid: str, species: str, domain: str,
               source_file: str, details_fields: List[Dict[str, Any]]) -> Optional[Dict[str, str]]:
    finding = _normalise_value(row.get("finding"))
    if not finding or finding.lower() in ("na", "n/a", "none"):
        return None

    raw_details = row.get("details", {})
    if not isinstance(raw_details, dict):
        raw_details = {}

    field_names = [f["field"] for f in details_fields if isinstance(f, dict) and f.get("field")]
    filtered_details = {k: v for k, v in raw_details.items() if k in field_names}

    clean: Dict[str, str] = {
        "pmid": _normalise_value(pmid),
        "pmcid": _normalise_value(pmcid),
        "species": _normalise_value(species),
        "domain": _normalise_value(domain),
        "topic": _normalise_value(row.get("topic", "unknown")).lower(),
        "claim_type": _normalise_value(row.get("claim_type", "unknown")).lower(),
        "product_name": _normalise_value(row.get("product_name")),
        "finding": finding,
        "quote": _normalise_value(row.get("quote")),
        "evidence_level": _normalise_value(row.get("evidence_level", "unknown")),
        "population": _normalise_value(row.get("population")),
        "geography": _normalise_value(row.get("geography")),
        "source_file": source_file,
    }
    # Domain-specific fields go to their own columns in the TSV.
    for fname in field_names:
        clean[fname] = _normalise_value(filtered_details.get(fname))

    clean["details"] = json.dumps(filtered_details, ensure_ascii=False, sort_keys=True)

    if clean["topic"] not in VALID_TOPICS:
        clean["topic"] = "unknown"
    if clean["claim_type"] not in VALID_CLAIM_TYPES:
        clean["claim_type"] = "unknown"
    if clean["evidence_level"] not in VALID_EVIDENCE_LEVELS:
        clean["evidence_level"] = "unknown"

    # Anti-hallucination: quote should appear in the text if present and not too long
    quote = clean["quote"]
    if quote and (quote not in text_for_quote_check or len(quote) > 2000):
        clean["quote"] = ""

    return clean


def _load_metadata(metadata_dir: Path, stem: str) -> Dict[str, Any]:
    meta_path = metadata_dir / f"{stem}.json"
    if not meta_path.is_file():
        return {}
    try:
        return json.loads(meta_path.read_text(encoding="utf-8"))
    except Exception:
        return {}


def main() -> None:
    args = parse_args()
    input_dir = Path(args.input_dir)
    metadata_dir = Path(args.metadata_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    templates = load_templates(Path(args.templates_yml))
    details_fields = _details_fields(templates, args.domain)

    host = args.ollama_host.rstrip("/")
    model = resolve_ollama_model(host, args.ollama_model)

    out_tsv = outdir / "evidence_extracted.tsv"
    out_log = outdir / "extraction_log.json"

    txt_files = sorted(input_dir.glob("*.txt"))
    if not txt_files:
        with open(out_tsv, "w", encoding="utf-8", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=EXPECTED_COLUMNS, delimiter="\t")
            writer.writeheader()
        log = {
            "species": args.species,
            "domain": args.domain,
            "model": model,
            "input_count": 0,
            "claim_count": 0,
            "success_count": 0,
            "failed_count": 0,
            "results": [],
        }
        out_log.write_text(json.dumps(log, ensure_ascii=False, indent=2), encoding="utf-8")
        print(f"[extract_evidence] No .txt files; wrote empty {out_tsv}", file=sys.stderr)
        return

    all_rows: List[Dict[str, str]] = []
    results = []
    global text_for_quote_check

    for txt in txt_files:
        stem = txt.stem
        metadata = _load_metadata(metadata_dir, stem)
        pmid = metadata.get("pmid") or stem
        pmcid = metadata.get("pmcid") or ""
        metadata["pmid"] = pmid
        metadata["pmcid"] = pmcid
        metadata["species"] = args.species

        text = txt.read_text(encoding="utf-8")
        text_for_quote_check = text
        if not text.strip():
            results.append({"pmid": pmid, "status": "empty_text", "claims": 0, "error": None})
            continue

        try:
            prompt = _build_prompt(args.domain, details_fields, metadata, text, args.n_ctx)
            raw = _call_ollama_generate(host, model, prompt, args.n_ctx, args.temperature, args.timeout)
            parsed = _parse_llm_json(raw)
            kept = []
            for row in parsed:
                if not isinstance(row, dict):
                    continue
                clean = _clean_row(row, pmid, pmcid, args.species, args.domain, txt.name, details_fields)
                if clean:
                    kept.append(clean)
            all_rows.extend(kept)
            results.append({
                "pmid": pmid,
                "status": "success",
                "claims": len(kept),
                "error": None,
            })
            print(f"[extract_evidence] PMID {pmid}: {len(kept)} claims", file=sys.stderr)
        except Exception as exc:
            print(f"[extract_evidence] PMID {pmid} failed: {exc}", file=sys.stderr)
            results.append({"pmid": pmid, "status": "failed", "claims": 0, "error": str(exc)})

    fieldnames = expected_columns(details_fields)
    with open(out_tsv, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in all_rows:
            writer.writerow({k: row.get(k, "") for k in fieldnames})

    log = {
        "species": args.species,
        "domain": args.domain,
        "model": model,
        "input_count": len(txt_files),
        "claim_count": len(all_rows),
        "success_count": sum(1 for r in results if r["status"] == "success"),
        "failed_count": sum(1 for r in results if r["status"] == "failed"),
        "results": results,
    }
    with open(out_log, "w", encoding="utf-8") as fh:
        json.dump(log, fh, ensure_ascii=False, indent=2)

    print(f"[extract_evidence] Wrote {out_tsv} ({len(all_rows)} claims)", file=sys.stderr)


if __name__ == "__main__":
    main()
