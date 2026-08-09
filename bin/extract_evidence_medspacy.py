#!/usr/bin/env python3
"""MedSpaCy-based evidence extraction from full-text .txt files.

Replaces AgentSLR + LLM approach with deterministic, rule-based NLP extraction
using MedSpaCy (spaCy + clinical NLP) with:
  - TargetRule matching for domain-specific concepts
  - ConText algorithm for negation/uncertainty filtering
  - Section detection for confidence scoring
  - Regex value extraction for numeric metrics

Produces per-PMID JSON files with the same schema as the previous AgentSLR output.
"""

import argparse
import json
import os
import re
import sys
from pathlib import Path

import yaml


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract structured evidence from full-text papers using MedSpaCy."
    )
    parser.add_argument("--input-dir", required=True, help="Directory containing .txt files")
    parser.add_argument("--outdir", required=True, help="Output directory for per-paper JSON files")
    parser.add_argument("--species", required=True, help="Pipeline species key (e.g. bdbv, ebov)")
    parser.add_argument("--domain", required=True, help="Extraction domain (e.g. diagnostic, seroprevalence)")
    parser.add_argument("--rules-yml", required=True, help="Path to evidence_rules.yml")
    parser.add_argument("--metadata-dir", default=None, help="Directory containing PubMed metadata JSON files")
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Metadata loading (same as run_agentslr_extraction.py)
# ---------------------------------------------------------------------------

METADATA_FIELDS = [
    "title", "authors", "year", "doi", "journal",
    "publication_date", "keywords", "pmcid",
    "cited_by_count", "is_oa",
]


def load_metadata_map(metadata_dir: Path) -> dict:
    """Load all metadata JSON files from a directory, keyed by PMID (filename stem)."""
    metadata_map = {}
    if metadata_dir is None or not metadata_dir.exists():
        return metadata_map
    for json_file in metadata_dir.glob("*.json"):
        try:
            with open(json_file, "r", encoding="utf-8") as f:
                data = json.load(f)
            metadata_map[json_file.stem] = data
        except (json.JSONDecodeError, OSError) as e:
            print(f"[extract_evidence] Warning: could not read metadata file {json_file}: {e}")
    return metadata_map


# ---------------------------------------------------------------------------
# Rules loading
# ---------------------------------------------------------------------------

def load_rules(rules_path: Path) -> dict:
    """Load evidence_rules.yml and return a dict of domain -> list of field rules."""
    with open(rules_path, "r", encoding="utf-8") as f:
        rules = yaml.safe_load(f)
    return rules


# ---------------------------------------------------------------------------
# MedSpaCy pipeline setup
# ---------------------------------------------------------------------------

def setup_nlp_pipeline():
    """Create a MedSpaCy NLP pipeline with target matcher, context, and sectionizer."""
    import medspacy

    nlp = medspacy.load(
        enable=["medspacy_target_matcher", "medspacy_context", "medspacy_sectionizer"]
    )

    return nlp


def build_target_rules(field_rules: list):
    """Convert YAML rule definitions into MedSpaCy TargetRule objects."""
    from medspacy.ner import TargetRule

    target_rules = []
    for field_def in field_rules:
        field_name = field_def["field"]
        for tr in field_def.get("target_rules", []):
            literal = tr["literal"]
            label = tr["label"]
            pattern = tr.get("pattern")
            if pattern:
                rule = TargetRule(literal, label, pattern=pattern)
            else:
                rule = TargetRule(literal, label)
            target_rules.append(rule)
    return target_rules


# ---------------------------------------------------------------------------
# Extraction logic
# ---------------------------------------------------------------------------

# ConText modifier categories to treat as "not present"
NEGATIVE_MODIFIERS = {"NEGATED_EXISTENCE", "UNCERTAIN", "HISTORICAL"}

# Section titles considered high-confidence evidence sources
HIGH_CONF_SECTIONS = {"abstract", "results", "result", "conclusion", "conclusions", "findings", "principal findings", "key findings"}
LOW_CONF_SECTIONS = {"references", "acknowledgments", "acknowledgements", "author contributions", "author information", "funding", "data availability", "conflicts of interest", "supplementary material"}

# Known symptom terms for normalization/deduplication
SYMPTOM_SYNONYMS = {
    "bleeding": ["hemorrhage", "haemorrhage"],
    "fever": [],
    "fatigue": ["tiredness", "exhaustion"],
    "headache": [],
    "vomiting": ["emesis"],
    "diarrhea": ["diarrhoea"],
    "abdominal pain": ["stomach pain"],
    "muscle pain": ["myalgia", "muscle ache"],
    "sore throat": [],
    "rash": [],
    "conjunctivitis": [],
    "hiccups": [],
    "internal bleeding": [],
    "ocular pain": [],
}


def is_negated_or_uncertain(span) -> bool:
    """Check if a span is modified by negation, uncertainty, or historical context."""
    modifiers = span._.modifiers
    if not modifiers:
        return False
    for mod in modifiers:
        mod_category = str(mod.category).upper()
        if mod_category in NEGATIVE_MODIFIERS:
            return True
    return False


def detect_section_by_regex(text: str) -> str:
    """Infer section title from common markdown-style headers in the text."""
    lines = text.splitlines()
    # Build a mapping of line start positions to section titles
    section_map = []
    pos = 0
    for line in lines:
        stripped = line.strip()
        if stripped.startswith("---") and stripped.endswith("---"):
            stripped = stripped.strip("-").strip()
        if re.match(r"^[A-Z][A-Za-z0-9\s/:–-]{2,60}$", stripped) or re.match(r"^#+\s+\w+", stripped):
            section_map.append((pos, stripped.lower().strip("# ")))
        pos += len(line) + 1
    return section_map


def get_section_title(span, section_map: list | None = None) -> str:
    """Get the section title for a given span, using MedSpaCy or regex fallback."""
    try:
        section = span._.section
        if section is not None:
            return str(section.title).lower() if hasattr(section, 'title') else ""
    except (AttributeError, KeyError):
        pass
    # Fallback: use regex-detected section headers based on span character offset
    if section_map and hasattr(span, "start_char"):
        char_offset = span.start_char
        current = ""
        for pos, title in section_map:
            if pos <= char_offset:
                current = title
            else:
                break
        return current
    return ""


def is_high_confidence_section(section_title: str) -> bool:
    return section_title and any(s in section_title for s in HIGH_CONF_SECTIONS)


def is_low_confidence_section(section_title: str) -> bool:
    return section_title and any(s in section_title for s in LOW_CONF_SECTIONS)


def context_matches(sentence_text: str, context_regex: str | None) -> bool:
    """Check whether the sentence contains required context keywords."""
    if not context_regex:
        return True
    return bool(re.search(context_regex, sentence_text, re.IGNORECASE))


def extract_value_with_regex(text: str, value_regex: str) -> str | None:
    """Apply a regex to the surrounding text and return the first non-empty group."""
    if not value_regex:
        return None
    for match in re.finditer(value_regex, text, re.IGNORECASE):
        if match.groups():
            for g in match.groups():
                if g:
                    return g.strip()
        else:
            return match.group(0).strip()
    return None


def get_sentence_text(span, doc) -> str:
    """Get the sentence containing the span."""
    try:
        sent = span.sent
        return sent.text
    except (AttributeError, ValueError):
        start = max(0, span.start - 50)
        end = min(len(doc), span.end + 50)
        return doc[start:end].text


def get_context_window(span, doc, window: int = 100) -> str:
    """Get a text window around the span for regex value extraction."""
    start = max(0, span.start - window)
    end = min(len(doc), span.end + window)
    return doc[start:end].text


def determine_confidence(section_title: str, has_value: bool) -> str:
    """Determine confidence level based on section context and value presence."""
    if is_low_confidence_section(section_title):
        return "low"
    if is_high_confidence_section(section_title):
        if has_value:
            return "high"
        return "medium"
    if has_value:
        return "medium"
    return "low"


def extract_field_value(field_def, span, doc):
    """Extract the value for a field from a matched span.

    Returns (value, quote) tuple. value is None when regex extraction fails and
    the field requires a value beyond the matched span text.
    """
    value_type = field_def.get("value_type", "string")
    value_regex = field_def.get("value_regex")
    sentence_text = get_sentence_text(span, doc)

    # Try regex value extraction first
    if value_regex:
        context_text = get_context_window(span, doc, window=140)
        regex_value = extract_value_with_regex(context_text, value_regex)
        if regex_value:
            return regex_value, sentence_text.strip()
        # If a value_regex is defined but returns nothing, do not fall back to the span text
        return None, sentence_text.strip()

    # For boolean fields: presence of the concept = true
    if value_type == "boolean":
        return "true", sentence_text.strip()

    # For number fields: try to find a number near the match
    if value_type == "number":
        context_text = get_context_window(span, doc, window=60)
        num_match = re.search(r"(\d[\d,]*)", context_text)
        if num_match:
            return num_match.group(1).replace(",", ""), sentence_text.strip()
        return None, sentence_text.strip()

    # For list and string types: use the matched span text
    if value_type == "list":
        return span.text, sentence_text.strip()

    # Default: use the matched literal text
    return span.text, sentence_text.strip()


def validate_percentage(value: str, field_name: str) -> bool:
    """Return True if value looks like a percentage between 0 and 100."""
    if value is None:
        return False
    try:
        num = float(value)
        return 0 <= num <= 100
    except (ValueError, TypeError):
        return field_name not in {"case_fatality_rate", "hospitalization_rate", "efficacy_value", "seroprevalence_value", "attack_rate", "sensitivity", "specificity"}


def normalize_symptoms(values: list) -> list:
    """Normalize symptom strings: lowercase, strip, merge synonyms/deduplicate."""
    canonical_map = {}
    for canonical, synonyms in SYMPTOM_SYNONYMS.items():
        canonical_map[canonical.lower()] = canonical
        for syn in synonyms:
            canonical_map[syn.lower()] = canonical
    seen = set()
    out = []
    for v in values:
        key = v.strip().lower().rstrip("s")
        # Map plural/synonym to canonical
        canonical = canonical_map.get(key, v.strip().lower())
        if canonical not in seen:
            seen.add(canonical)
            out.append(canonical)
    return out


def deduplicate_string_values(values: list) -> list:
    """Case-insensitive deduplication preserving first-seen casing."""
    seen = set()
    out = []
    for v in values:
        key = v.strip().lower()
        if key and key not in seen:
            seen.add(key)
            out.append(v.strip())
    return out


def is_keyword_only(value, field_name: str) -> bool:
    """Reject values that are just trigger words without real content."""
    if value is None:
        return True
    s = str(value).strip().lower()
    keyword_triggers = {
        "geographic_location": {"outbreak in", "conducted in", "located in", "reported in", "identified in", "implemented in", "applicable to", "study site"},
        "efficacy_value": {"efficacy", "effectiveness", "survival rate"},
        "effectiveness": {"effectiveness", "impact", "reduction in", "decreased transmission"},
        "coverage": {"coverage", "adherence", "compliance"},
        "severity_score": {"disease severity", "severity score", "severity classification"},
        "hospitalization_rate": {"hospitalization rate", "admission rate", "hospitalized"},
        "risk_factors": {"risk factor", "associated with", "contact with"},
        "age_group": {"age group"},
    }
    return s in keyword_triggers.get(field_name, set())


def validate_field_value(field_name: str, value, allowed_values: list | None) -> bool:
    """Run field-specific validation rules."""
    if value is None:
        return False
    if is_keyword_only(value, field_name):
        return False
    if field_name == "severity_score" and allowed_values:
        if str(value).strip().lower() not in {v.lower() for v in allowed_values}:
            return False
    if field_name in {"case_fatality_rate", "hospitalization_rate", "efficacy_value", "seroprevalence_value", "attack_rate", "sensitivity", "specificity"}:
        return validate_percentage(value, field_name)
    return True


def post_process_extraction(extraction: list, field_rules: list) -> list:
    """Apply validation, normalization, deduplication, and cross-field consistency."""
    field_map = {e["field"]: e for e in extraction}
    field_def_map = {fr["field"]: fr for fr in field_rules}

    for entry in extraction:
        field = entry["field"]
        value = entry["value"]
        field_def = field_def_map.get(field, {})
        allowed = field_def.get("allowed_values")
        value_type = field_def.get("value_type", "string")

        # Normalize list values
        if value_type == "list" and isinstance(value, list):
            if field == "symptoms":
                value = normalize_symptoms(value)
            else:
                value = deduplicate_string_values(value)
            entry["value"] = value if value else None
            entry["present"] = bool(value)

        # Validate non-list values
        elif value_type != "list":
            if not validate_field_value(field, value, allowed):
                entry["value"] = None
                entry["present"] = False
                entry["confidence"] = "low"

    # Cross-field consistency: hospitalization_rate should not equal case_fatality_rate
    cfr = field_map.get("case_fatality_rate", {}).get("value")
    hosp = field_map.get("hospitalization_rate", {}).get("value")
    if cfr and hosp and str(cfr) == str(hosp):
        # Reject hospitalization_rate unless it was explicitly supported by hospital context
        # (value_regex already enforces that), but if it equals CFR it is likely a duplicate.
        field_map["hospitalization_rate"]["value"] = None
        field_map["hospitalization_rate"]["present"] = False
        field_map["hospitalization_rate"]["confidence"] = "low"

    return extraction


def process_paper(txt_path: Path, nlp, field_rules: list, label_to_field: dict) -> dict:
    """Process a single paper and extract evidence for all fields.

    Args:
        txt_path: Path to the .txt file
        nlp: MedSpaCy NLP pipeline
        field_rules: List of field definitions from evidence_rules.yml
        label_to_field: Dict mapping label -> list of field_names

    Returns:
        Dict with extraction results
    """
    article_id = txt_path.stem
    text = txt_path.read_text(encoding="utf-8", errors="replace")

    if not text.strip():
        return {
            "article_id": article_id,
            "extraction": [],
            "not_found": [fr["field"] for fr in field_rules],
            "status": "no_evidence",
        }

    # Detect section headers before NLP to handle PDF-extracted text
    section_map = detect_section_by_regex(text)

    # Run NLP pipeline
    doc = nlp(text)

    # Collect matches per field
    field_results = {}
    for fr in field_rules:
        field_results[fr["field"]] = {
            "field": fr["field"],
            "present": False,
            "value": None,
            "quote": None,
            "confidence": "low",
            "value_type": fr.get("value_type", "string"),
            "value_regex": fr.get("value_regex"),
            "context_regex": fr.get("context_regex"),
            "context_filters": fr.get("context_filters", []),
            "section_boost": fr.get("section_boost", []),
            "allowed_values": fr.get("allowed_values"),
            "matches": [],
        }

    def _process_span(span):
        label = span.label_
        if label not in label_to_field:
            return
        # Check ConText modifiers
        if is_negated_or_uncertain(span):
            return
        section_title = get_section_title(span, section_map)
        # Skip matches in low-confidence sections
        if is_low_confidence_section(section_title):
            return
        sentence_text = get_sentence_text(span, doc)
        for field_name in label_to_field[label]:
            field_def = next(fr for fr in field_rules if fr["field"] == field_name)
            # Validate sentence context
            if not context_matches(sentence_text, field_def.get("context_regex")):
                continue
            value, quote = extract_field_value(field_def, span, doc)
            if value is None:
                continue
            confidence = determine_confidence(section_title, value is not None)
            field_results[field_name]["matches"].append((value, quote, confidence))

    # Process entities from doc.ents
    for ent in doc.ents:
        _process_span(ent)

    # Also check span groups (MedSpaCy may put matches in span groups)
    try:
        spans = doc.spans.get("medspacy_spans", [])
        for span in spans:
            _process_span(span)
    except (KeyError, AttributeError):
        pass

    # Aggregate matches per field
    extraction = []
    not_found = []

    for field_name, info in field_results.items():
        matches = info["matches"]
        if not matches:
            not_found.append(field_name)
            extraction.append({
                "field": field_name,
                "present": False,
                "value": None,
                "quote": None,
                "confidence": "low",
            })
            continue

        value_type = info["value_type"]

        if value_type == "list":
            values = []
            quotes = []
            for v, q, c in matches:
                if v and v not in values:
                    values.append(v)
                    quotes.append(q)
            value = values if values else None
            quote = " | ".join(quotes) if quotes else None
        elif value_type == "boolean":
            value = "true"
            quote = matches[0][1]
        elif value_type == "number":
            value = matches[0][0] if matches else None
            quote = matches[0][1] if matches else None
        else:
            # For string fields, prefer high/medium confidence and deduplicate
            best = None
            seen = set()
            for v, q, c in matches:
                if v and v.lower() not in seen:
                    seen.add(v.lower())
                    if best is None or c == "high" or (best[2] == "low" and c == "medium"):
                        best = (v, q, c)
            value = best[0] if best else None
            quote = best[1] if best else None

        confidences = [c for _, _, c in matches]
        if "high" in confidences:
            confidence = "high"
        elif "medium" in confidences:
            confidence = "medium"
        else:
            confidence = "low"

        extraction.append({
            "field": field_name,
            "present": True,
            "value": value,
            "quote": quote,
            "confidence": confidence,
        })

    # Apply post-processing
    extraction = post_process_extraction(extraction, field_rules)

    present_count = sum(1 for e in extraction if e["present"])
    status = "success" if present_count > 0 else "no_evidence"

    return {
        "article_id": article_id,
        "extraction": extraction,
        "not_found": not_found,
        "status": status,
    }


# ---------------------------------------------------------------------------
# Output formatting (same schema as AgentSLR output)
# ---------------------------------------------------------------------------

def build_output_json(result: dict, species: str, domain: str, metadata_map: dict) -> dict:
    """Build the final per-PMID JSON output, merging metadata."""
    article_id = result["article_id"]
    output = {
        "pmid": article_id,
        "species": species,
        "domain": domain,
        "extraction": result["extraction"],
        "not_found": result["not_found"],
        "status": result["status"],
    }

    meta = metadata_map.get(str(article_id))
    if meta:
        for field in METADATA_FIELDS:
            if field in meta:
                output[field] = meta[field]
        if "pmid" in meta:
            output["pmid"] = meta["pmid"]
    else:
        if metadata_map:
            print(f"[extract_evidence] Warning: no metadata found for PMID {article_id}")
        for field in METADATA_FIELDS:
            output[field] = None

    return output


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()

    input_dir = Path(args.input_dir).resolve()
    outdir = Path(args.outdir).resolve()
    rules_path = Path(args.rules_yml).resolve()
    metadata_dir = Path(args.metadata_dir).resolve() if args.metadata_dir else None

    # Load rules
    all_rules = load_rules(rules_path)
    domain_rules = all_rules.get(args.domain)
    if domain_rules is None:
        print(f"[extract_evidence] ERROR: domain '{args.domain}' not found in {rules_path}")
        sys.exit(1)

    print(f"[extract_evidence] Loaded {len(domain_rules)} field rules for domain '{args.domain}'")

    # Setup MedSpaCy pipeline
    print("[extract_evidence] Initializing MedSpaCy pipeline...")
    nlp = setup_nlp_pipeline()

    # Build and add target rules
    target_rules = build_target_rules(domain_rules)
    target_matcher = nlp.get_pipe("medspacy_target_matcher")
    target_matcher.add(target_rules)
    print(f"[extract_evidence] Added {len(target_rules)} target rules to pipeline")

    # Build label-to-field mapping
    label_to_field = {}
    for fr in domain_rules:
        for tr in fr.get("target_rules", []):
            label_to_field.setdefault(tr["label"], []).append(fr["field"])

    # Load metadata
    metadata_map = load_metadata_map(metadata_dir) if metadata_dir else {}
    if metadata_map:
        print(f"[extract_evidence] Loaded {len(metadata_map)} metadata records from {metadata_dir}")

    # Find .txt files
    txt_files = sorted(input_dir.glob("*.txt"))
    if not txt_files:
        print("[extract_evidence] No .txt files found. Writing empty output.")
        outdir.mkdir(parents=True, exist_ok=True)
        return

    print(f"[extract_evidence] Found {len(txt_files)} .txt files for species={args.species}, domain={args.domain}")

    # Process each paper
    outdir.mkdir(parents=True, exist_ok=True)
    results = []

    for txt_file in txt_files:
        print(f"[extract_evidence] Processing: {txt_file.name}")
        result = process_paper(txt_file, nlp, domain_rules, label_to_field)

        output = build_output_json(result, args.species, args.domain, metadata_map)
        out_path = outdir / f"{result['article_id']}.json"
        with open(out_path, "w", encoding="utf-8") as f:
            json.dump(output, f, indent=2, ensure_ascii=False)

        present_count = sum(1 for e in result["extraction"] if e["present"])
        total = len(result["extraction"])
        print(f"  -> {present_count}/{total} fields extracted (status: {result['status']})")

        results.append(output)

    # Summary
    total_papers = len(results)
    success_count = sum(1 for r in results if r["status"] == "success")
    total_fields = sum(len(r["extraction"]) for r in results)
    present_fields = sum(1 for r in results for e in r["extraction"] if e["present"])

    print(f"\n[extract_evidence] Summary:")
    print(f"  Papers processed: {total_papers}")
    print(f"  Papers with evidence: {success_count}")
    print(f"  Fields extracted: {present_fields}/{total_fields} ({100*present_fields/max(total_fields,1):.1f}%)")
    print(f"[extract_evidence] Done. Output: {outdir}")


if __name__ == "__main__":
    main()
