#!/usr/bin/env python3
"""summarize_mutation_profile.py

Generate two LLM-assisted outputs for one species' mutation profile at
pipeline time, so the dashboard renders text without calling Ollama:

1. mutation_summary.json  — a plain-language 'mutation landscape' note:
   query-vs-background separation (PLS-DA + burden test), proteins with
   elevated mutation load in the query genomes, phenotype-linked mutations
   and what they plausibly affect, and query-enriched positions.
2. protein_summaries.tsv  — one short digest per protein condensed from the
   UniProt function/domain text attached to mutated positions (the raw text
   is rich but repetitive and hard to follow in the UI).

Fallback is a deterministic template summary — both files are always written.

Usage:
    python3 summarize_mutation_profile.py \
        --profile-dir mutation_profile --species ebov --outdir . \
        [--ollama-host http://localhost:11434] [--ollama-model NAME]
"""

import argparse
import csv
import json
import re
import sys
import urllib.request
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

THINKING_FAMILIES = ("qwen3", "deepseek", "gpt-oss")
MAX_PROTEIN_DIGESTS = 12


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--profile-dir", required=True,
                   help="Directory containing the 01_*.tsv mutation profile tables")
    p.add_argument("--species", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--ollama-host", default="http://localhost:11434")
    p.add_argument("--ollama-model", default="")
    p.add_argument("--temperature", type=float, default=0.3)
    p.add_argument("--timeout", type=int, default=60)
    return p.parse_args()


def ollama_generate(host, model, prompt, temperature, timeout):
    url = host.rstrip("/") + "/api/generate"
    body = {
        "model": model, "prompt": prompt, "stream": False,
        "options": {"temperature": temperature, "num_predict": 800},
    }
    req = urllib.request.Request(url, data=json.dumps(body).encode("utf-8"),
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


def read_tsv(path):
    if not path.exists():
        return []
    with open(path, encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def fnum(v):
    try:
        return float(str(v).strip())
    except (TypeError, ValueError):
        return None


def truthy(v):
    return str(v).strip().lower() in ("true", "1", "yes")


def build_facts(profile_dir):
    facts = {}

    # -- Overall burden: query vs background ---------------------------------
    rows = read_tsv(profile_dir / "01_mutation_burden_summary.tsv")
    for r in rows:
        side = "query" if truthy(r.get("is_query")) else "background"
        facts.setdefault("burden", {})[side] = {
            "n": fnum(r.get("n")),
            "mean_mutation_positions": fnum(r.get("mean_burden")),
            "test": r.get("test", ""),
            "p_value": fnum(r.get("p_value")),
            "effect": fnum(r.get("estimate")),
        }

    # -- Per-protein burden: query vs background ------------------------------
    prot_rows = read_tsv(profile_dir / "01_protein_burden_summary.tsv")
    by_prot = defaultdict(dict)
    for r in prot_rows:
        p = r.get("protein_name", "").strip()
        if not p:
            continue
        side = "query" if truthy(r.get("is_query")) else "background"
        by_prot[p][side] = {
            "n": fnum(r.get("n")),
            "mean_mutated_positions": fnum(r.get("mean_mutation_positions")),
            "p_value": fnum(r.get("p_value")),
            "effect": fnum(r.get("estimate")),
        }
    diffs = []
    for p, d in by_prot.items():
        q, b = d.get("query"), d.get("background")
        if not q or not b or q["mean_mutated_positions"] is None:
            continue
        dms = q["mean_mutated_positions"] - (b["mean_mutated_positions"] or 0)
        diffs.append({"protein": p,
                      "query_mean_positions": q["mean_mutated_positions"],
                      "background_mean_positions": b["mean_mutated_positions"],
                      "difference": round(dms, 3),
                      "p_value": d.get("query", {}).get("p_value")
                      or by_prot[p].get("p_value")})
    diffs.sort(key=lambda d: abs(d["difference"]), reverse=True)
    facts["protein_burden"] = diffs[:8]

    # -- PLS-DA: query/background separation + drivers -------------------------
    scores = read_tsv(profile_dir / "01_plsda_scores.tsv")
    if scores:
        q = [r for r in scores if truthy(r.get("is_query"))]
        b = [r for r in scores if not truthy(r.get("is_query"))]
        def centroid(rs):
            xs = [fnum(r.get("PC1")) for r in rs if fnum(r.get("PC1")) is not None]
            ys = [fnum(r.get("PC2")) for r in rs if fnum(r.get("PC2")) is not None]
            return (sum(xs) / len(xs), sum(ys) / len(ys)) if xs else None
        qc, bc = centroid(q), centroid(b)
        dist = None
        if qc and bc:
            dist = round(((qc[0] - bc[0]) ** 2 + (qc[1] - bc[1]) ** 2) ** 0.5, 3)
        facts["plsda"] = {"n_query": len(q), "n_background": len(b),
                          "query_centroid": qc, "background_centroid": bc,
                          "centroid_distance": dist}
    vip = read_tsv(profile_dir / "01_plsda_vip.tsv")
    vip = [(r.get("protein_name", ""), fnum(r.get("vip")))
           for r in vip if fnum(r.get("vip")) is not None]
    vip.sort(key=lambda t: t[1], reverse=True)
    facts["plsda_top_drivers"] = [{"protein": p, "vip": v}
                                  for p, v in vip[:6] if p]

    # -- Phenotype-linked mutations -------------------------------------------
    pheno = read_tsv(profile_dir / "01_mutation_phenotypes.tsv")
    seen = set()
    hits = []
    for r in pheno:
        key = (r.get("mutation_label", ""), r.get("protein_name", ""),
               r.get("phenotype", ""))
        if key in seen or not r.get("mutation_label"):
            continue
        seen.add(key)
        hits.append({"mutation": r.get("mutation_label"),
                     "protein": r.get("protein_name"),
                     "phenotype": r.get("phenotype"),
                     "effect": r.get("effect"),
                     "source": r.get("source")})
        if len(hits) >= 12:
            break
    facts["phenotype_mutations"] = hits

    # -- Query-enriched positions ---------------------------------------------
    pos = read_tsv(profile_dir / "01_position_summary.tsv")
    enriched = []
    for r in pos:
        nq, nb = fnum(r.get("n_query_mutated")), fnum(r.get("n_bg_mutated"))
        if not truthy(r.get("has_query_mutation")) or nq is None:
            continue
        enriched.append({
            "protein": r.get("protein_name"),
            "position": r.get("position"),
            "reference_aa": r.get("reference_aa"),
            "query_mutated": int(nq),
            "background_mutated": int(nb or 0),
            "entropy": fnum(r.get("shannon_entropy")),
        })
    enriched.sort(key=lambda d: d["query_mutated"], reverse=True)
    facts["query_enriched_positions"] = enriched[:10]
    return facts


def template_summary(facts, species):
    parts = [f"Mutation landscape for {species.upper()}:"]
    bur = facts.get("burden", {})
    q, b = bur.get("query"), bur.get("background")
    if q and b and q.get("mean_mutation_positions") is not None:
        p = q.get("p_value")
        parts.append(
            f"query genomes carry on average {q['mean_mutation_positions']:.1f} "
            f"mutated positions vs {b.get('mean_mutation_positions', 0):.1f} in "
            f"background isolates" +
            (f" (p={p:.3g})" if p is not None else "") + ".")
    pls = facts.get("plsda", {})
    if pls.get("centroid_distance") is not None:
        parts.append(
            f"query and background separate on the PLS-DA axes "
            f"(centroid distance {pls['centroid_distance']}).")
    tops = facts.get("protein_burden") or []
    if tops:
        parts.append("Proteins most divergent in the query: "
                     + ", ".join(d["protein"] for d in tops[:4]) + ".")
    pm = facts.get("phenotype_mutations") or []
    if pm:
        parts.append("Phenotype-linked mutations include "
                     + ", ".join(h["mutation"] for h in pm[:5]) + ".")
    return " ".join(parts)


LANDSCAPE_PROMPT = (
    "You are a genomic epidemiologist describing the mutation landscape of "
    "query {species} genomes against background isolates for a public-health "
    "audience. From the JSON facts below, write 4-6 sentences covering: "
    "(1) whether the query genomes form a distinct cluster from the "
    "background (PLS-DA centroid distance and per-protein burden tests); "
    "(2) which proteins show an elevated mutation load in the query and "
    "what that may imply; (3) key mutations linked to a known biological "
    "effect or phenotype — name them and say what the mutation may do, "
    "including its likely evolutionary origin (recently arisen in the query "
    "vs shared with background); (4) one caveat about the analysis. Use only "
    "the numbers given; plain language; no markdown.\n\nFACTS:\n{facts}")


def protein_digest_prompt(protein, texts):
    joined = "\n---\n".join(texts[:3])
    return (
        f"Summarise the biological function of the {protein} protein for a "
        "public-health genomics audience in 1-2 plain sentences. The UniProt "
        "text below is long and repetitive — distil only what the protein "
        "does and why it matters for the virus. No markdown.\n\nTEXT:\n"
        + joined)


def main():
    args = parse_args()
    profile_dir = Path(args.profile_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    host = args.ollama_host.rstrip("/")
    model = resolve_model(host, args.ollama_model)
    if not model:
        print("[mutation_summary] no Ollama model; template only",
              file=sys.stderr)

    facts = build_facts(profile_dir)
    now = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

    # -- Landscape narrative --------------------------------------------------
    text, source, error = None, "template", None
    if model and any(facts.values()):
        try:
            text = ollama_generate(
                host, model,
                LANDSCAPE_PROMPT.format(species=args.species.upper(),
                                        facts=json.dumps(facts)),
                args.temperature, args.timeout)
        except Exception as exc:
            error = str(exc)
            print(f"[mutation_summary] Ollama failed: {exc}", file=sys.stderr)
    payload = {
        "species": args.species,
        "text": text or template_summary(facts, args.species),
        "source": "ollama" if text else "template",
        "model": model if text else "", "error": error,
        "generated_at": now,
    }
    with open(outdir / "mutation_summary.json", "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)

    # -- Per-protein UniProt function digests ----------------------------------
    ctx = read_tsv(profile_dir / "01_mutation_context.tsv")
    func_by_prot = defaultdict(list)
    for r in ctx:
        if (r.get("context_type") or "").strip().lower() != "function":
            continue
        p = (r.get("protein_name") or "").strip()
        d = (r.get("context_description") or "").strip()
        if p and d and d not in func_by_prot[p]:
            func_by_prot[p].append(d)

    digest_rows = []
    for prot, texts in sorted(func_by_prot.items())[:MAX_PROTEIN_DIGESTS]:
        digest = None
        if model:
            try:
                digest = ollama_generate(
                    host, model, protein_digest_prompt(prot, texts),
                    args.temperature, args.timeout)
            except Exception as exc:
                print(f"[mutation_summary] digest failed for {prot}: {exc}",
                      file=sys.stderr)
        src = "ollama"
        if not digest:
            first = texts[0].split(".")[0].strip()
            digest = (first + ".") if first else ""
            src = "template"
        digest_rows.append({"protein_name": prot, "summary": digest,
                            "source": src,
                            "model": model if src == "ollama" else ""})

    with open(outdir / "protein_summaries.tsv", "w", encoding="utf-8",
              newline="") as f:
        w = csv.DictWriter(f, fieldnames=["protein_name", "summary", "source",
                                          "model"], delimiter="\t")
        w.writeheader()
        w.writerows(digest_rows)

    print(f"[mutation_summary] wrote mutation_summary.json + "
          f"{len(digest_rows)} protein digests (model={model or 'none'})",
          file=sys.stderr)


if __name__ == "__main__":
    main()
