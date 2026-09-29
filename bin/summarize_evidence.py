#!/usr/bin/env python3
"""summarize_evidence.py

Generate a plain-language narrative per literature domain (and one
species-level overview) from the extracted evidence TSVs, using the same
local Ollama model the extractor uses. Output is a single
domain_summaries.tsv the knowledge warehouse loads, so the dashboard can
display text without calling the LLM itself.

Usage:
    python3 summarize_evidence.py \
        --input-dir . --species ebov --outdir . \
        [--ollama-host http://localhost:11434] [--ollama-model NAME]
"""

import argparse
import csv
import json
import re
import sys
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

# Reasoning-model families burn num_predict on hidden <think> tokens —
# same guard as bin/extract_literature_evidence.py.
THINKING_FAMILIES = ("qwen3", "deepseek", "gpt-oss")

UNIVERSAL_COLS = {
    "pmid", "pmcid", "species", "domain", "topic", "claim_type", "finding",
    "quote", "evidence_level", "population", "geography", "source_file", "n",
    "details",
}

DOMAIN_LABELS = {
    "diagnostic": "diagnostic tests",
    "vaccine_therapeutic": "vaccines and therapeutics",
    "intervention": "interventions",
    "policy": "policies and recommendations",
    "clinical": "clinical presentation",
    "epidemiological_context": "epidemiological context",
    "transmission": "transmission",
    "seroprevalence": "seroprevalence",
    "reservoir": "animal reservoirs",
    "genomic_surveillance": "genomic surveillance",
    "variant_phenotype": "variant phenotypes",
}

FREQUENCY_FIELDS = {
    "policy": ["issuing_body", "policy_type", "intervention", "travel_measure",
               "policy_outcome", "target_population", "geographic_scope"],
    "intervention": ["intervention_name", "intervention_type",
                     "implementation_setting", "target_population"],
    "reservoir": ["host_species", "detection_method", "geographic_scope"],
    "clinical": ["clinical_manifestations", "symptoms", "disease_severity",
                 "treatment_setting"],
    "transmission": ["transmission_route", "transmission_setting",
                     "risk_factors"],
    "genomic_surveillance": ["sequencing_platform", "lineage", "clade",
                             "genotype", "geographic_distribution"],
    "variant_phenotype": ["gene", "mutation", "phenotype"],
    "diagnostic": ["test_name", "test_type", "sample_type"],
    "vaccine_therapeutic": ["product_name", "product_type", "study_phase"],
    "epidemiological_context": ["outbreak_name", "geographic_scope"],
    "seroprevalence": ["host_species", "population"],
}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input-dir", required=True,
                   help="Directory containing staged evidence_*.tsv files")
    p.add_argument("--species", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--ollama-host", default="http://localhost:11434")
    p.add_argument("--ollama-model", default="")
    p.add_argument("--n-ctx", type=int, default=8192)
    p.add_argument("--temperature", type=float, default=0.3)
    p.add_argument("--timeout", type=int, default=60)
    return p.parse_args()


def ollama_generate(host, model, prompt, n_ctx, temperature, timeout):
    url = host.rstrip("/") + "/api/generate"
    body = {
        "model": model,
        "prompt": prompt,
        "stream": False,
        "options": {"temperature": temperature, "num_predict": 400,
                    "num_ctx": n_ctx},
    }
    req = urllib.request.Request(
        url, data=json.dumps(body).encode("utf-8"),
        headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        raw = json.loads(resp.read().decode("utf-8"))
    text = (raw.get("response") or "").strip()
    text = re.sub(r"<think>.*?</think>", "", text, flags=re.S).strip()
    return text


def resolve_model(host, model):
    if model:
        return model
    try:
        with urllib.request.urlopen(host.rstrip("/") + "/api/tags",
                                    timeout=5) as resp:
            tags = json.loads(resp.read().decode("utf-8"))
        names = [(m.get("name") or m.get("model") or "")
                 for m in tags.get("models", [])]
        names = [n for n in names if n]
        if not names:
            return ""
        non_thinking = [n for n in names
                        if not any(f in n.lower() for f in THINKING_FAMILIES)]
        for n in non_thinking:
            if n.startswith("qwen2.5"):
                return n
        return (non_thinking or names)[0]
    except Exception:
        return ""


def read_rows(path):
    with open(path, encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def nonempty(v):
    return bool(v) and str(v).strip() not in ("", "null", "NA")


def top_values(rows, field, k=5):
    counts = {}
    for r in rows:
        v = (r.get(field) or "").strip()
        if not nonempty(v):
            continue
        for part in v.replace(";", "\n").split("\n"):
            part = part.strip().strip("[]\"'")
            if part and part.lower() not in ("null", "na", "none"):
                counts[part] = counts.get(part, 0) + 1
    return sorted(counts.items(), key=lambda kv: kv[1], reverse=True)[:k]


def build_facts(domain, rows):
    pmids = {r.get("pmid", "") for r in rows if nonempty(r.get("pmid"))}
    facts = {
        "n_claims": len(rows),
        "n_papers": len(pmids),
        "frequencies": {},
        "findings": [],
    }
    for f in FREQUENCY_FIELDS.get(domain, []):
        tops = top_values(rows, f)
        if tops:
            facts["frequencies"][f] = dict(tops)
    seen = set()
    for r in rows:
        f_ = (r.get("finding") or "").strip()
        if nonempty(f_) and f_ not in seen:
            seen.add(f_)
            facts["findings"].append(f_)
        if len(facts["findings"]) >= 40:
            break
    return facts


def template_summary(domain, facts, species):
    label = DOMAIN_LABELS.get(domain, domain.replace("_", " "))
    if not facts["n_claims"]:
        return (f"No extracted claims are available for {label} in "
                f"{species.upper()}.")
    parts = [f"The {label} literature for {species.upper()} yields "
             f"{facts['n_claims']} claims across {facts['n_papers']} papers."]
    for f, counts in list(facts["frequencies"].items())[:3]:
        tops = "; ".join(sorted(counts, key=counts.get, reverse=True)[:3])
        parts.append(f"Most-mentioned {f.replace('_', ' ')}: {tops}.")
    return " ".join(parts)


def llm_summary(domain, facts, species, host, model, n_ctx, temperature,
                timeout):
    if not model:
        return None
    label = DOMAIN_LABELS.get(domain, domain.replace("_", " "))
    if not facts["n_claims"]:
        return template_summary(domain, facts, species)
    prompt = (
        "You are an epidemiologist summarising extracted literature evidence "
        f"on {label} for {species.upper()}. Below is JSON with per-field "
        "frequencies and the raw findings (up to 40). Write 3-5 plain "
        "sentences for a decision-maker: what the evidence actually says, "
        "the most common themes and any clear gaps. Use only what is given; "
        "no invented numbers; no markdown.\n\nFACTS:\n"
        + json.dumps(facts, ensure_ascii=False))
    try:
        return ollama_generate(host, model, prompt, n_ctx, temperature,
                               timeout) or None
    except Exception as exc:
        print(f"[summarize_evidence] Ollama failed for {domain}: {exc}",
              file=sys.stderr)
        return None


def overview_facts(per_domain):
    out = {}
    for domain, (rows, facts) in per_domain.items():
        out[domain] = {
            "n_claims": facts["n_claims"],
            "n_papers": facts["n_papers"],
            "top": {f: list(c)[:3] for f, c in facts["frequencies"].items()},
        }
    return out


# -- Evidence highlights (brief) ------------------------------------------------

PCT_METRICS = {
    "vaccine_therapeutic": ("efficacy", ("product_name", "product_type",
                                        "intervention_name", "target")),
    "diagnostic":          ("sensitivity", ("test_name", "product_name",
                                            "test_type",
                                            "target_gene_or_antigen")),
    "intervention":        ("effectiveness", ("intervention_name",
                                              "product_name",
                                              "intervention_type")),
}


def parse_pct(v):
    """Parse a proportion/percent field to 0-100."""
    if not nonempty(v):
        return None
    s = str(v).strip().rstrip("%")
    m = re.search(r"-?\d+(?:\.\d+)?", s)
    if not m:
        return None
    x = float(m.group(0))
    if x <= 1.0:
        x *= 100
    return x


def group_label(r, fields):
    for f in fields:
        v = (r.get(f) or "").strip()
        if nonempty(v):
            return v
    return ""


def top_products(per_domain, k=6):
    """Median effectiveness metric per concrete product/test/intervention."""
    out = []
    for domain, (metric, groups) in PCT_METRICS.items():
        rows = per_domain.get(domain, ([], {}))[0]
        agg = {}
        for r in rows:
            g = group_label(r, groups)
            v = parse_pct(r.get(metric))
            if g and v is not None:
                agg.setdefault(g, []).append(v)
        for g, vs in agg.items():
            vs.sort()
            med = vs[len(vs) // 2] if len(vs) % 2 else (vs[len(vs) // 2 - 1]
                                                      + vs[len(vs) // 2]) / 2
            out.append({"product": g, "domain": domain,
                        "metric": metric, "value": round(med, 1),
                        "n": len(vs)})
    out.sort(key=lambda d: d["value"], reverse=True)
    return out[:k]


def highlights_prompt(facts, species):
    return (
        "You are an epidemiologist writing field-brief highlights on "
        f"{species.upper()} for a colleague about to deploy. From the JSON "
        "below — per-domain literature summaries, claim coverage, gaps and "
        "top countermeasure products — write 4-6 plain sentences: what is "
        "known about this pathogen from published evidence, which concrete "
        "diagnostics/vaccines/therapeutics exist and their key metrics, and "
        "the most important knowledge gap. Name products and numbers; "
        "summarise only what is given — no outside knowledge, no markdown.\n\n"
        "FACTS:\n" + json.dumps(facts, ensure_ascii=False))


def write_highlights(outdir, species, facts, text, source, model, generated_at):
    payload = {
        "species": species, "text": text, "source": source,
        "model": model if source == "ollama" else "",
        "generated_at": generated_at,
        "domain_claims": facts["domain_claims"],
        "top_products": facts["top_products"],
        "papers_with_claims": facts["papers_with_claims"],
        "domains_with_evidence": facts["domains_with_evidence"],
        "domains_with_gap": facts["domains_with_gap"],
    }
    with open(Path(outdir) / "evidence_highlights.json", "w",
              encoding="utf-8") as f:
        json.dump(payload, f, indent=2, ensure_ascii=False)


def overview_prompt(facts, species):
    return (
        "You are an epidemiologist writing a short 'known vs unknown' "
        f"briefing on {species.upper()} from per-domain literature evidence "
        "counts and top values. Write 3-5 sentences: what is well evidenced, "
        "what is thin or missing, and the most important gap for a "
        "decision-maker. Plain language, no markdown, no invented numbers.\n\n"
        "FACTS:\n" + json.dumps(facts, ensure_ascii=False))


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out_tsv = outdir / "domain_summaries.tsv"

    host = args.ollama_host.rstrip("/")
    model = resolve_model(host, args.ollama_model)
    if not model:
        print("[summarize_evidence] no Ollama model; template summaries only",
              file=sys.stderr)

    tsvs = sorted(Path(args.input_dir).glob("*.tsv"))
    per_domain = {}
    for p in tsvs:
        rows = read_rows(p)
        domains = {r.get("domain", "") for r in rows}
        domain = domains.pop() if len(domains) == 1 else p.stem
        per_domain[domain] = (rows, build_facts(domain, rows))

    results = []
    now = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    for domain in sorted(per_domain):
        rows, facts = per_domain[domain]
        text = llm_summary(domain, facts, args.species, host, model,
                           args.n_ctx, args.temperature, args.timeout)
        source = "ollama" if text else "template"
        results.append({
            "species": args.species, "domain": domain,
            "n_claims": facts["n_claims"], "n_papers": facts["n_papers"],
            "summary": text or template_summary(domain, facts, args.species),
            "source": source, "model": model if source == "ollama" else "",
            "generated_at": now,
        })

    ov_facts = overview_facts(per_domain)
    if model:
        try:
            text = ollama_generate(host, model, overview_prompt(ov_facts, args.species),
                                   args.n_ctx, args.temperature, args.timeout)
        except Exception as exc:
            print(f"[summarize_evidence] overview Ollama call failed: {exc}",
                  file=sys.stderr)
            text = None
    else:
        text = None
    n_claims = sum(f["n_claims"] for _, f in per_domain.values())
    n_papers = sum(f["n_papers"] for _, f in per_domain.values())
    present = sorted(d for d, (_, f) in per_domain.items() if f["n_claims"])
    absent = sorted(set(DOMAIN_LABELS) - set(per_domain))
    fallback = (
        f"{args.species.upper()} literature produced {n_claims} claims across "
        f"{n_papers} papers. Domains with extracted evidence: "
        + (", ".join(present) if present else "none")
        + (". Domains with no extraction: " + ", ".join(absent) if absent else ".")
    )
    results.append({
        "species": args.species, "domain": "species_overview",
        "n_claims": n_claims, "n_papers": n_papers,
        "summary": text or fallback,
        "source": "ollama" if text else "template",
        "model": model if text else "", "generated_at": now,
    })

    # -- Intelligence-brief highlights ------------------------------------------
    products = top_products(per_domain)
    hl_facts = {
        "domain_summaries": {r["domain"]: r["summary"] for r in results
                             if r["domain"] != "species_overview"},
        "domain_claims": {d: f["n_claims"] for d, (_, f)
                          in per_domain.items()},
        "top_products": products,
        "papers_with_claims": n_papers,
        "domains_with_evidence": present,
        "domains_with_gap": [d for d in sorted(DOMAIN_LABELS)
                             if d not in present],
    }
    hl_text = None
    if model:
        try:
            hl_text = ollama_generate(host, model,
                                      highlights_prompt(hl_facts, args.species),
                                      args.n_ctx, args.temperature, args.timeout)
        except Exception as exc:
            print(f"[summarize_evidence] highlights Ollama failed: {exc}",
                  file=sys.stderr)
    if not hl_text:
        prod_txt = ("; ".join(f"{p['product']} ({p['metric']} "
                              f"{p['value']}%, n={p['n']})"
                              for p in products[:3])
                    if products else "none quantified")
        hl_text = (fallback
                   + f" Countermeasure products with metrics: {prod_txt}.")
    write_highlights(outdir, args.species, hl_facts, hl_text,
                     "ollama" if hl_text and model else "template",
                     model, now)

    cols = ["species", "domain", "n_claims", "n_papers", "summary", "source",
            "model", "generated_at"]
    with open(out_tsv, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(results)
    print(f"[summarize_evidence] wrote {out_tsv} "
          f"({len(results)} summaries, model={model or 'none'})",
          file=sys.stderr)


if __name__ == "__main__":
    main()
