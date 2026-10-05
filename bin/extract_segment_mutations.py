#!/usr/bin/env python3
"""
extract_segment_mutations.py

Produce a per-segment mutation table (``mutations_{segment}.tsv``) for every
genome segment a Nextstrain build actually produced, covering query AND
background tips. Two producers, chosen per segment:

  A. nextclade  — when a staged reference.fasta + genemap.gff exist for the
     segment (seasonal flu: config/{lineage}/{segment}/ inside the build
     work dir; internal segments get theirs from the Nextclade datasets
     NEXTSTRAIN_FLU fetches). Runs `nextclade run` on the build's own
     per-segment aligned.fasta.
  B. auspice    — when no nextclade-usable reference exists (avian flu is
     augur-based), derive ref-relative AA mutations per tip by walking the
     segment's Auspice JSON and accumulating branch_attrs.mutations.aa
     (and .nuc for the count) along each root→tip path.

Output contract (read by extract_query_proteins.py):
  seqName, aaSubstitutions ("GENE:X123Y,GENE:..."), totalSubstitutions

Per-segment failures are logged and skipped — a broken segment never kills
the group.

Inputs:
  --results_dir  Nextstrain build results/work directory
  --pathogen     meta.pathogen (influenza | avian_influenza | ...)
  --lineage      meta.species (h1n1pdm | h3n2 | vic | h5nx | h5n1 | ...)
  --auspice      Zero or more Auspice JSONs (per-segment trees)
  --outdir       Output directory for mutations_{segment}.tsv
  --report       Log/report file path
"""

import argparse
import csv
import glob
import json
import os
import re
import subprocess
import sys

SEG_ORDER = ["pb2", "pb1", "pa", "ha", "np", "na", "mp", "ns"]
MUT_RE = re.compile(r"^([A-Za-z0-9_\-]+):([A-Z*])(\d+)([A-Z*])$")
AA_MUT_RE = re.compile(r"^([A-Z*])(\d+)([A-Z*])$")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--results_dir", required=True)
    p.add_argument("--pathogen", default="")
    p.add_argument("--lineage", default="")
    p.add_argument("--auspice", nargs="*", default=[])
    p.add_argument("--outdir", default=".")
    p.add_argument("--report", default="segment_mutations_report.txt")
    return p.parse_args()


class Reporter:
    def __init__(self, path):
        self.fh = open(path, "w")

    def log(self, msg):
        print(msg, file=sys.stderr)
        self.fh.write(msg + "\n")

    def close(self):
        self.fh.close()


def segment_from_auspice_name(path):
    """Infer the segment token from an Auspice JSON filename.

    seasonal: {lineage}_pgirl_{segment}.json        -> last '_' token
    avian:    avian-flu_{subtype}_{segment}_{time}  -> second-to-last token
    Tip-frequencies sidecars and unrecognized names return None.
    """
    stem = os.path.basename(path)
    if "tip-frequencies" in stem or not stem.endswith(".json"):
        return None
    stem = stem[:-5]
    parts = stem.split("_")
    if stem.startswith("avian-flu") and len(parts) >= 2:
        cand = parts[-2]
    else:
        cand = parts[-1]
    return cand if cand in SEG_ORDER else None


def segment_from_path(path):
    """Extract the segment component from a build path containing
    /{segment}/ between build-level dirs (builds/*/pb2/..., results/*/ha/...).
    """
    parts = path.split(os.sep)
    for i, part in enumerate(parts):
        if part in SEG_ORDER and i > 0:
            return part
    return None


def find_segment_reference(results_dir, segment, lineage, rep):
    """Locate a staged reference.fasta + genemap.gff for a segment.

    Seasonal-flu work dirs keep per-subtype copies at
    config/{subtype}/{segment}/ for every subtype the build supports, so a
    bare glob picks the alphabetically first subtype (h1n1 before h3n2) and
    produces all-failed alignments. Only accept a candidate under
    /{lineage}/{segment}/; a subtype mismatch means no usable nextclade
    reference exists and the auspice producer (subtype-agnostic) takes over.
    Full datasets ds_{segment}/ are subtype-specific already — always safe.
    Returns dict with 'dataset_dir' or 'reference'+'genemap', or None.
    """
    # Full nextclade dataset fetched by the build (internal segments).
    for ds in sorted(glob.glob(os.path.join(results_dir, "**", f"ds_{segment}"), recursive=True)):
        if os.path.exists(os.path.join(ds, "pathogen.json")):
            return {"dataset_dir": ds}
    # Staged reference + genemap (vendored ha/na; converted internals).
    candidates = []
    for ref in sorted(glob.glob(os.path.join(results_dir, "**", segment, "reference.fasta"), recursive=True)):
        gff = os.path.join(os.path.dirname(ref), "genemap.gff")
        if os.path.exists(gff):
            candidates.append(ref)
    if not candidates:
        return None
    if lineage:
        wanted = f"{os.sep}{lineage}{os.sep}{segment}{os.sep}"
        for ref in candidates:
            if wanted in ref + os.sep:
                return {"reference": ref, "genemap": os.path.join(os.path.dirname(ref), "genemap.gff")}
        rep.log(f"  [{segment}] staged reference(s) exist only under other subtype configs "
                f"({', '.join(candidates)}); lineage '{lineage}' not matched — "
                f"skipping nextclade producer")
        return None
    return {"reference": candidates[0], "genemap": os.path.join(os.path.dirname(candidates[0]), "genemap.gff")}


def nextclade_bin():
    """nextclade 3.x ships both `nextclade` and `nextclade3` entrypoints
    depending on the package; use whichever is on PATH."""
    for cand in ("nextclade3", "nextclade"):
        try:
            subprocess.run([cand, "--version"], capture_output=True, timeout=30)
            return cand
        except (OSError, subprocess.TimeoutExpired):
            continue
    return None


def run_nextclade(alignment, ref_info, segment, out_tsv, rep):
    """Producer A: `nextclade run` on the build's per-segment sequences."""
    binary = nextclade_bin()
    if not binary:
        rep.log(f"  [{segment}] nextclade not on PATH; skipping producer A")
        return False
    cmd = [binary, "run", alignment]
    if ref_info.get("dataset_dir"):
        cmd += ["--input-dataset", ref_info["dataset_dir"]]
    else:
        cmd += ["-r", ref_info["reference"], "-m", ref_info["genemap"]]
    cmd += ["--output-tsv", out_tsv, "--include-reference", "--silent"]
    rep.log(f"  [{segment}] nextclade: {' '.join(cmd)}")
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
    except (OSError, subprocess.TimeoutExpired) as e:
        rep.log(f"  [{segment}] nextclade FAILED to run: {e}")
        return False
    if proc.returncode != 0 or not os.path.exists(out_tsv):
        rep.log(f"  [{segment}] nextclade FAILED: {proc.stderr.strip()[-500:]}")
        return False
    rep.log(f"  [{segment}] nextclade -> {out_tsv}")
    return True


def auspice_tip_mutations(auspice_path, rep, segment):
    """Producer B: ref-relative AA (and nuc-count) mutations per tip.

    Walks the Auspice v2 tree; each node's branch_attrs.mutations.aa holds
    the mutations on its incoming branch as "GENE:REFposALT". Accumulating
    position->(ref,alt) along the path collapses reversions: if the final
    alt equals the first-recorded ref at that position, the entry is dropped.
    Returns {tip_name: ({"GENE":["X123Y"]}, nuc_count)}.
    """
    with open(auspice_path) as fh:
        auspice = json.load(fh)
    root = auspice.get("tree")
    if not root:
        rep.log(f"  [{segment}] no 'tree' in {os.path.basename(auspice_path)}")
        return {}

    results = {}

    def branch_muts(node):
        """Auspice stores mutations as {label: [mutstr]} where 'nuc' holds
        nucleotide mutations and every other key is a gene/annotation name
        whose entries are bare 'X123Y' strings (no gene prefix). Some exports
        may instead carry an 'aa' key with 'GENE:X123Y' entries — handle both.
        Returns ([(gene, ref, pos, alt)], nuc_list)."""
        muts = node.get("branch_attrs", {}).get("mutations", {})
        aa = []
        for label, entries in muts.items():
            if label == "nuc" or not isinstance(entries, list):
                continue
            for entry in entries:
                if ":" in entry:
                    m = MUT_RE.match(entry)
                    if m:
                        aa.append((m.group(1), m.group(2), int(m.group(3)), m.group(4)))
                else:
                    m = AA_MUT_RE.match(entry)
                    if m:
                        aa.append((label, m.group(1), int(m.group(2)), m.group(3)))
        return aa, muts.get("nuc") or []

    def walk(node, aa_map, nuc_count):
        aa_here = {g: dict(v) for g, v in aa_map.items()}
        aa_list, nuc_list = branch_muts(node)
        nuc_here = nuc_count + len(nuc_list)
        for gene, ref, pos, alt in aa_list:
            slot = aa_here.setdefault(gene, {})
            if pos in slot:
                slot[pos] = (slot[pos][0], alt)  # keep original ref, latest alt
            else:
                slot[pos] = (ref, alt)
        name = node.get("name", "")
        children = node.get("children") or []
        if not children and name:
            out = {}
            nuc = nuc_here
            for gene, positions in aa_here.items():
                muts = [f"{r}{p}{a}" for p, (r, a) in sorted(positions.items()) if a != r]
                if muts:
                    out[gene] = muts
            results[name] = (out, nuc)
        for child in children:
            walk(child, aa_here, nuc_here)

    walk(root, {}, 0)
    rep.log(f"  [{segment}] auspice walk: {len(results)} tips from {os.path.basename(auspice_path)}")
    return results


def write_tsv(rows, out_tsv, segment):
    """Emit the contract TSV. `rows`: [(seqName, {gene:[muts]}, nuc_count)]"""
    with open(out_tsv, "w") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["seqName", "aaSubstitutions", "totalSubstitutions"])
        for name, gene_muts, nuc in sorted(rows):
            flat = ",".join(f"{g}:{m}" for g in sorted(gene_muts) for m in gene_muts[g])
            w.writerow([name, flat, nuc])


def main():
    args = parse_args()
    rep = Reporter(args.report)
    rep.log(f"extract_segment_mutations: results_dir={args.results_dir} "
            f"pathogen={args.pathogen} lineage={args.lineage}")
    os.makedirs(args.outdir, exist_ok=True)

    # --- Segment -> artifacts discovery -------------------------------------
    seg_auspice = {}
    for js in args.auspice:
        seg = segment_from_auspice_name(js)
        if seg:
            seg_auspice[seg] = js
        else:
            rep.log(f"  ignoring non-segment Auspice file: {js}")

    # Also discover Auspice JSONs inside results_dir (avian keeps them in a
    # sibling auspice/ dir, so the explicit --auspice args normally cover it;
    # seasonal work dirs contain them under {species}/auspice/).
    for js in glob.glob(os.path.join(args.results_dir, "**", "auspice", "*.json"), recursive=True):
        seg = segment_from_auspice_name(js)
        if seg and seg not in seg_auspice:
            seg_auspice[seg] = js

    # nextclade input per segment: prefer the build's unaligned merged
    # sequences.fasta (nextclade re-aligns anyway); aligned.fasta fallback.
    seg_alignment = {}
    for pat in ("**/sequences.fasta", "**/aligned.fasta"):
        for aln in sorted(glob.glob(os.path.join(args.results_dir, pat), recursive=True)):
            seg = segment_from_path(aln)
            if seg and seg not in seg_alignment:
                seg_alignment[seg] = aln

    segments = [s for s in SEG_ORDER if s in set(seg_auspice) | set(seg_alignment)]
    rep.log(f"  segments discovered: {segments or 'NONE'} "
            f"(auspice={sorted(seg_auspice)}, alignments={sorted(seg_alignment)})")
    if not segments:
        rep.log("  no segment artifacts found; nothing to emit")
        rep.close()
        return

    # --- Per-segment production ---------------------------------------------
    emitted = []
    for seg in segments:
        out_tsv = os.path.join(args.outdir, f"mutations_{seg}.tsv")
        try:
            ref_info = find_segment_reference(args.results_dir, seg, args.lineage, rep)
            aln = seg_alignment.get(seg)
            done = False
            if ref_info and aln:
                done = run_nextclade(aln, ref_info, seg, out_tsv, rep)
            elif ref_info and not aln:
                rep.log(f"  [{seg}] reference found but no aligned.fasta; trying auspice")
            if not done:
                js = seg_auspice.get(seg)
                if not js:
                    rep.log(f"  [{seg}] no auspice JSON either; segment skipped")
                    continue
                tip_muts = auspice_tip_mutations(js, rep, seg)
                if not tip_muts:
                    rep.log(f"  [{seg}] auspice produced no tip mutations; skipped")
                    continue
                write_tsv([(n, g, c) for n, (g, c) in tip_muts.items()], out_tsv, seg)
                rep.log(f"  [{seg}] auspice-derived -> {out_tsv} ({len(tip_muts)} tips)")
            emitted.append(seg)
        except Exception as e:  # per-segment isolation — never kill the group
            rep.log(f"  [{seg}] ERROR: {e}; segment skipped")

    rep.log(f"  emitted mutations TSVs for: {emitted or 'NONE'}")
    rep.close()


if __name__ == "__main__":
    main()
