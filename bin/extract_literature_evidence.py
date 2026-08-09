#!/usr/bin/env python3
"""Extract structured, quotable, confidence-scored evidence from paper text files.

Reads plain-text outputs of `LITERATURE_TEXT` and uses a local llama-cpp model
to attempt a per-domain evidence template defined in `database/evidence_templates.yml`.
"""
import argparse
import datetime
import json
import os
import re
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional

import yaml

try:
    from llama_cpp import Llama
    _HAS_LLAMA = True
except Exception:
    _HAS_LLAMA = False


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Extract structured evidence from paper text.")
    parser.add_argument("--input-dir", required=True, help="Directory with .txt paper files.")
    parser.add_argument("--outdir", required=True, help="Directory for JSON output and summary.")
    parser.add_argument("--species", required=True, help="Species key.")
    parser.add_argument("--domain", required=True, help="Domain key.")
    parser.add_argument("--templates-yml", required=True, help="Path to evidence_templates.yml.")
    parser.add_argument("--gguf", default=None, help="Path to local .gguf model file.")
    parser.add_argument("--n-ctx", type=int, default=4096, help="LLM context size.")
    parser.add_argument("--temperature", type=float, default=0.1, help="LLM temperature.")
    parser.add_argument("--seed", type=int, default=42, help="LLM seed.")
    return parser.parse_args()


def load_templates(path: Path) -> Dict[str, Any]:
    with open(path, "r", encoding="utf-8") as fh:
        data = yaml.safe_load(fh)
    return data


def pmid_from_filename(p: Path) -> str:
    return p.stem


def _model_name(gguf: Optional[str]) -> Optional[str]:
    if not gguf:
        return None
    return Path(gguf).name


def _truncate_text(text: str, model: Llama, max_tokens: int) -> str:
    """Keep the start of the text that fits within max_tokens."""
    try:
        tokens = model.tokenize(text.encode("utf-8"), add_bos=False)
        if len(tokens) <= max_tokens:
            return text
        kept = tokens[:max_tokens]
        return model.detokenize(kept).decode("utf-8", errors="ignore")
    except Exception:
        # crude fallback: ~4 chars per token
        return text[: max_tokens * 4]


def _build_prompt(domain: str, fields: List[Dict[str, Any]], text: str) -> str:
    field_descriptions = []
    for f in fields:
        optional = "optional" if f.get("optional") else "required if present in text"
        field_descriptions.append(f"- {f['field']} ({f.get('type','string')}, {optional}): {f['question']}")

    instructions = (
        "You are an exact evidence extractor. Read the paper text below and extract "
        "only the facts that match the requested variables. "
        "Return a JSON object with a single key 'extraction' whose value is a list of objects. "
        "Each object in the list must have: field, present (true or false), value (string, null if absent), "
        "quote (exact verbatim sentence or phrase from the text, null if absent), and confidence (high/medium/low). "
        "For any numeric, rate, or percentage value, also include these optional keys when present in the text: "
        "numerator (number), denominator (number), ci_95_lower (number), ci_95_upper (number), "
        "p_value (string), unit (string), group (string), n (number). "
        "If the text does not contain the requested fact, set present to false and all other keys to null. "
        "Do not infer or synthesise information. Use only explicit statements from the text.\n\n"
        "Variables to extract:\n" + "\n".join(field_descriptions) + "\n\n"
        f"Paper text ({domain} domain):\n{text}\n\n"
        "Respond with valid JSON only."
    )
    return instructions


def _call_llm(prompt: str, model: Llama, n_ctx: int, temperature: float, seed: int) -> Optional[str]:
    try:
        prompt_tokens = model.tokenize(prompt.encode("utf-8"), add_bos=False)
        max_tokens = max(64, n_ctx - len(prompt_tokens) - 64)
        messages = [
            {"role": "system", "content": "You are a precise structured evidence extractor that outputs only valid JSON."},
            {"role": "user", "content": prompt},
        ]
        # Try structured JSON mode; fall back to plain chat completion if unsupported.
        try:
            resp = model.create_chat_completion(
                messages,
                response_format={"type": "json_object"},
                max_tokens=max_tokens,
                temperature=temperature,
                seed=seed,
            )
        except TypeError:
            resp = model.create_chat_completion(
                messages,
                max_tokens=max_tokens,
                temperature=temperature,
                seed=seed,
            )
        return resp["choices"][0]["message"]["content"].strip()
    except Exception as exc:
        print(f"[extract_evidence] LLM call failed: {exc}", file=sys.stderr)
        return None


def _parse_json(raw: str) -> Optional[Dict[str, Any]]:
    raw = raw.strip()
    # Remove possible markdown fences
    if raw.startswith("```"):
        raw = re.sub(r"^```(?:json)?\s*", "", raw)
        raw = re.sub(r"\s*```$", "", raw)
    try:
        return json.loads(raw)
    except Exception as exc:
        print(f"[extract_evidence] JSON parse failed: {exc}", file=sys.stderr)
        return None


def _coerce_output(parsed: Optional[Dict[str, Any]], expected_fields: List[str], text: str) -> List[Dict[str, Any]]:
    if not parsed:
        return []
    extr = parsed.get("extraction") or parsed.get("fields") or parsed.get("results") or []
    by_field = {}
    for item in extr:
        if not isinstance(item, dict):
            continue
        fname = item.get("field")
        if not fname or fname not in expected_fields:
            continue
        rec = {
            "field": fname,
            "present": bool(item.get("present", False)),
            "value": item.get("value") if item.get("value") not in (None, "", "null") else None,
            "quote": item.get("quote") if item.get("quote") not in (None, "", "null") else None,
            "confidence": item.get("confidence", "low"),
            "numerator": _to_num(item.get("numerator")) if "numerator" in item else None,
            "denominator": _to_num(item.get("denominator")) if "denominator" in item else None,
            "ci_95_lower": _to_num(item.get("ci_95_lower")) if "ci_95_lower" in item else None,
            "ci_95_upper": _to_num(item.get("ci_95_upper")) if "ci_95_upper" in item else None,
            "p_value": item.get("p_value") if "p_value" in item else None,
            "unit": item.get("unit") if "unit" in item else None,
            "group": item.get("group") if "group" in item else None,
            "n": _to_num(item.get("n")) if "n" in item else None,
        }
        # Anti-hallucination: quote must be in the original text
        if rec["present"] and rec["quote"]:
            if rec["quote"] not in text:
                rec["confidence"] = "low"
                rec["quote"] = None
        by_field[fname] = rec

    # Fill missing expected fields as not found
    results = []
    for f in expected_fields:
        if f in by_field:
            results.append(by_field[f])
        else:
            results.append({
                "field": f,
                "present": False,
                "value": None,
                "quote": None,
                "confidence": None,
                "numerator": None,
                "denominator": None,
                "ci_95_lower": None,
                "ci_95_upper": None,
                "p_value": None,
                "unit": None,
                "group": None,
                "n": None,
            })
    return results


def _to_num(v: Any) -> Optional[float]:
    if v is None:
        return None
    try:
        if isinstance(v, (int, float)):
            return float(v)
        s = str(v).replace(",", "")
        return float(s)
    except Exception:
        return None


def _write_evidence_json(
    outdir: Path,
    pmid: str,
    species: str,
    domain: str,
    model_name: Optional[str],
    status: str,
    extraction: List[Dict[str, Any]],
) -> None:
    not_found = [e["field"] for e in extraction if not e["present"]]
    record = {
        "pmid": pmid,
        "species": species,
        "domain": domain,
        "model": model_name,
        "extraction_timestamp": datetime.datetime.utcnow().isoformat() + "Z",
        "status": status,
        "extraction": extraction,
        "not_found": not_found,
    }
    out_path = outdir / f"{pmid}.json"
    with open(out_path, "w", encoding="utf-8") as fh:
        json.dump(record, fh, ensure_ascii=False, indent=2)


def main() -> None:
    args = parse_args()
    input_dir = Path(args.input_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    templates = load_templates(Path(args.templates_yml))
    domain_fields = templates.get(args.domain, [])
    expected_fields = [f["field"] for f in domain_fields]

    txt_files = sorted(input_dir.glob("*.txt"))
    if not txt_files:
        print("[extract_evidence] No .txt files found to process.", file=sys.stderr)

    model = None
    model_name = _model_name(args.gguf)
    if _HAS_LLAMA and args.gguf and Path(args.gguf).is_file():
        try:
            print(f"[extract_evidence] Loading model {model_name}...", file=sys.stderr)
            model = Llama(
                model_path=args.gguf,
                n_ctx=args.n_ctx,
                n_batch=512,
                verbose=False,
            )
        except Exception as exc:
            print(f"[extract_evidence] Failed to load model: {exc}", file=sys.stderr)
            model = None

    results = []
    for txt in txt_files:
        pmid = pmid_from_filename(txt)
        text = txt.read_text(encoding="utf-8")

        if not model:
            print(f"[extract_evidence] Skipping PMID {pmid}: no model available.", file=sys.stderr)
            _write_evidence_json(outdir, pmid, args.species, args.domain, model_name, "skipped_model_missing", [])
            results.append({"pmid": pmid, "status": "skipped_model_missing", "present_count": 0, "not_found_count": len(expected_fields), "error": "no model available"})
            continue

        prompt = _build_prompt(args.domain, domain_fields, text)
        # Ensure prompt fits in context; keep system+instructions and truncate text
        prompt_tokens = model.tokenize(prompt.encode("utf-8"), add_bos=False)
        if len(prompt_tokens) > args.n_ctx - 256:
            # Over-estimation; rebuild with truncated text
            fixed_intro = _build_prompt(args.domain, domain_fields, "")
            intro_tokens = model.tokenize(fixed_intro.encode("utf-8"), add_bos=False)
            allowed_text = args.n_ctx - len(intro_tokens) - 256
            text = _truncate_text(text, model, allowed_text)
            prompt = _build_prompt(args.domain, domain_fields, text)

        raw = _call_llm(prompt, model, args.n_ctx, args.temperature, args.seed)
        if raw is None:
            _write_evidence_json(outdir, pmid, args.species, args.domain, model_name, "failed", [])
            results.append({"pmid": pmid, "status": "failed", "present_count": 0, "not_found_count": 0, "error": "LLM call failed"})
            continue

        parsed = _parse_json(raw)
        extraction = _coerce_output(parsed, expected_fields, text)
        present_count = sum(1 for e in extraction if e["present"])

        status = "success" if present_count > 0 else "no_evidence"
        _write_evidence_json(outdir, pmid, args.species, args.domain, model_name, status, extraction)
        results.append({
            "pmid": pmid,
            "status": status,
            "present_count": present_count,
            "not_found_count": len(expected_fields) - present_count,
            "error": None,
        })
        print(f"[extract_evidence] PMID {pmid}: {present_count}/{len(expected_fields)} fields present", file=sys.stderr)

    summary = {
        "species": args.species,
        "domain": args.domain,
        "model": model_name,
        "input_count": len(txt_files),
        "success_count": sum(1 for r in results if r["status"] == "success"),
        "no_evidence_count": sum(1 for r in results if r["status"] == "no_evidence"),
        "failed_count": sum(1 for r in results if r["status"] in ("failed", "skipped_model_missing")),
        "results": results,
    }
    with open(outdir / "evidence_summary.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, ensure_ascii=False, indent=2)
    print(f"[extract_evidence] Wrote {outdir/'evidence_summary.json'}", file=sys.stderr)


if __name__ == "__main__":
    main()
