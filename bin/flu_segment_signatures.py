#!/usr/bin/env python3
"""SEGMENT_SIGNATURES — influenza downstream signature tables.

Derives five per-group analysis tables from the per-segment mutation TSVs
(SEGMENT_MUTATIONS), the classification species assignments, per-record
Nextclade QC tables, GenoFLU results (avian) and the curated
database/flu_markers.yml catalogue:

  constellation.tsv      isolate x segment clade/genotype labels + the
                         PB2|PB1|PA|HA|NP|NA|MP|NS signature string and a
                         novelty flag vs the group modal signature.
  reassortment_flags.tsv isolates flagged by species_assignments or by
                         per-segment label discordance vs the modal
                         constellation.
  di_candidates.tsv      defective-interfering screen from Nextclade QC
                         fields (deletions / frameshifts / stop codons /
                         failed CDS / coverage) — a consensus-level screen;
                         true DI confirmation needs read-level coverage.
  markers.tsv            query-isolate mutations joined against the curated
                         resistance/virulence catalogue, with coordinate
                         translation for HA sub-CDS labels and a caveat flag
                         for cross-numbering-system matches.
  group_diversity.tsv    per group_id x segment shared/private mutation
                         counts and mean pairwise difference between
                         isolates sharing a group_id (intra-host / same-site
                         diversity proxy from consensus genomes).

The script deliberately uses only stdlib + pandas + yaml so it can run in
the existing pgirl_nextstrain conda env.
"""

import argparse
import csv
import json
import re
import sys
from pathlib import Path

try:
    import yaml
except ImportError:  # pragma: no cover
    yaml = None

SEGMENT_ORDER = ["pb2", "pb1", "pa", "ha", "np", "na", "mp", "ns"]

# Mutation labels that sit inside a larger polyprotein CDS. For influenza the
# only split CDS is HA (SigPep/HA1/HA2 -> HA0); the marker catalogue's
# segment_cds coordinates are per-protein, so unsplit genes need no
# translation at all.
SUBCDS_TO_PARENT = {"sigpep": "HA", "ha1": "HA", "ha2": "HA"}

SEQNAME_CAND_RE = re.compile(r"_Seg\d+$", re.IGNORECASE)


def log(msg):
    print(f"[flu_segment_signatures] {msg}", file=sys.stderr)


# ---------------------------------------------------------------------------
# Input loading
# ---------------------------------------------------------------------------

def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def seg_from_filename(path):
    m = re.search(r"mutations_([a-zA-Z0-9]+)\.tsv$", Path(path).name)
    return m.group(1).lower() if m else None


def load_assignments(path, pathogen=None, species=None):
    """species_assignments.tsv -> {record: {isolate, segment, best_dataset_file,
    reassortment_suspected, species}}. The file covers every routed group, so
    rows are scoped to this group: exact species match first (keeps same-
    pathogen subtype groups apart, e.g. h1n1pdm vs h3n2), falling back to
    pathogen-only when the resolved species renamed upstream (h5nx -> h5n1
    after GenoFLU) would otherwise empty the set."""
    rows = {}
    n_species_match = 0
    for r in read_tsv(path):
        record = (r.get("sample") or "").strip()
        if not record:
            continue
        if pathogen and (r.get("pathogen") or "").strip() != pathogen:
            continue
        if species and (r.get("species") or "").strip() != species:
            continue
        n_species_match += 1
        rows[record] = {
            "isolate": (r.get("isolate") or record).strip(),
            "segment": (r.get("segment") or "").strip().lower(),
            "record_species": (r.get("record_species") or "").strip(),
            "best_dataset_file": (r.get("best_dataset_file") or "").strip(),
            "reassortment_suspected": str(
                r.get("reassortment_suspected", "")
            ).lower() in ("true", "1", "yes"),
            "species": (r.get("species") or "").strip(),
        }
    if species and not n_species_match and pathogen:
        # Resolved species renamed after assignments were written — retry
        # scoped by pathogen only.
        return load_assignments(path, pathogen=pathogen, species=None)
    return rows


def load_nextclade_tables(paths):
    """Classification nextclade TSVs -> {tsv_name: {seqName: row}}."""
    tables = {}
    for p in paths:
        p = Path(p)
        if not p.exists() or p.name.startswith("NO_FILE"):
            continue
        try:
            rows = read_tsv(p)
        except Exception:
            continue
        if not rows or "seqName" not in rows[0]:
            continue
        tables[p.name] = {(r.get("seqName") or "").strip(): r for r in rows}
    return tables


def _to_float(v):
    try:
        return float(str(v).strip())
    except (TypeError, ValueError):
        return None


def _to_int(v):
    f = _to_float(v)
    return int(f) if f is not None else None


def nc_row(tables, record, seg=None, best_file=None):
    """Row for a record, preferring (a) the record's winning dataset — every
    record appears in every dataset TSV, so a sibling dataset's row would
    carry the wrong clade/QC — then (b) a dataset TSV whose filename names
    the segment, then (c) any TSV holding it."""
    def hit(name):
        t = tables.get(name)
        return t.get(record) if t else None

    if best_file:
        r = hit(best_file)
        if r:
            return r, best_file
    if seg:
        seg_hits = [n for n in tables if re.search(rf"[_./-]{re.escape(seg)}[_.-]", n, re.I)]
        for name in seg_hits:
            r = hit(name)
            if r:
                return r, name
    for name, t in tables.items():
        r = t.get(record)
        if r:
            return r, name
    return None, None


def parse_aa_muts(s):
    """'HA1:D122N,HA2:I32T' or 'HA:K3N' -> [(gene, ref, pos, alt)]."""
    out = []
    for tok in (s or "").split(","):
        tok = tok.strip()
        m = re.match(r"^([^:]+):([A-Z*X-]?)(\d+)([A-Z*X-]*)$", tok)
        if not m:
            continue
        gene, ref, pos, alt = m.group(1), m.group(2), int(m.group(3)), m.group(4)
        if alt.upper() in ("X", "*", "-"):
            # 'X' tails in Auspice-derived calls are alignment-edge
            # unknown/stop artifacts, not clean AA variants.
            continue
        out.append((gene, ref, pos, alt))
    return out


def load_mutations(paths):
    """mutations_{seg}.tsv -> {seg: [rows]} keeping seqName, clade, aa list."""
    segs = {}
    for p in paths:
        p = Path(p)
        if not p.exists() or p.name.startswith("NO_FILE"):
            continue
        seg = seg_from_filename(p)
        if not seg:
            continue
        try:
            rows = read_tsv(p)
        except Exception:
            continue
        if not rows or "seqName" not in rows[0]:
            continue
        parsed = []
        for r in rows:
            name = (r.get("seqName") or "").strip()
            if not name:
                continue
            parsed.append({
                "seqName": name,
                "clade": (r.get("clade") or "").strip() or None,
                "aa": parse_aa_muts(r.get("aaSubstitutions")),
            })
        segs[seg] = parsed
    return segs


def load_metadata_groups(path):
    """metadata.tsv -> {strain -> group_id} (tolerates missing column/file)."""
    groups = {}
    p = Path(path) if path else None
    if not p or not p.exists() or p.name.startswith("NO_FILE"):
        return groups
    try:
        for r in read_tsv(p):
            name = (r.get("strain") or r.get("sample") or r.get("sample_id") or "").strip()
            gid = (r.get("group_id") or "").strip()
            if name:
                groups[name] = gid or None
    except Exception:
        pass
    return groups


def load_genoflu(path):
    """genoflu_results.tsv -> {isolate: {seg: genotype}} — per-segment labels
    parsed from 'Genotype Sample Title List' entries 'geno:strain:SEGMENT'."""
    out = {}
    p = Path(path) if path else None
    if not p or not p.exists() or p.name.startswith("NO_FILE"):
        return out
    for r in read_tsv(p):
        strain = (r.get("Strain") or "").strip()
        titles = (r.get("Genotype Sample Title List") or "").split(",")
        seg_map = {}
        for t in titles:
            parts = t.strip().rsplit(":", 2)
            if len(parts) == 3:
                seg_map[parts[2].lower()] = parts[0]
        if strain and seg_map:
            out[strain] = seg_map
        if strain:
            out.setdefault(strain, seg_map)
    return out


def load_build_clades(build_dir):
    """<seg>/clades.json in the build results dir -> {seg: {name: clade}}."""
    out = {}
    bd = Path(build_dir) if build_dir else None
    if not bd or not bd.exists() or bd.name.startswith("NO_FILE"):
        return out
    for f in bd.rglob("clades.json"):
        try:
            data = json.loads(f.read_text())
        except Exception:
            continue
        seg = f.parent.name.lower()
        m = {}
        for name, attrs in (data.get("nodes") or {}).items():
            cl = (attrs or {}).get("clade_membership") if isinstance(attrs, dict) else None
            if cl:
                m[SEQNAME_CAND_RE.sub("", name)] = cl
        if m:
            out[seg] = m
    return out


def load_ha_offsets(build_dir):
    """Segment genemap.gff -> per sub-CDS aa offsets needed to translate
    SigPep/HA1/HA2 mutation labels onto the HA0 protein and mature-HA
    (h_numbering) coordinates.

    Returns {gene_lower: {'sigpep_start': int, 'sigpep_len': int,
             'cds_offset': int, 'h_offset': int}} where a mutation `GENE:x`
    maps to HA0 residue x + cds_offset and h-number residue x + h_offset."""
    bd = Path(build_dir) if build_dir else None
    if not bd or not bd.exists():
        return None
    for gff in bd.rglob("genemap.gff"):
        parts = [x.lower() for x in gff.parts]
        if "ha" not in parts:
            continue
        feats = []
        for line in gff.read_text().splitlines():
            if line.startswith("#") or not line.strip():
                continue
            f = line.split("\t")
            if len(f) < 9 or f[2] not in ("gene", "CDS"):
                continue
            m = re.search(r'gene_name="?([^";]+)', f[8])
            if m:
                feats.append((m.group(1), int(f[3]), int(f[4])))
        if not feats:
            continue
        sig = next((x for x in feats if x[0].lower() in ("sigpep", "signalpeptide", "sp")), None)
        anchor_start = sig[1] if sig else feats[0][1]
        sig_len = ((sig[2] - sig[1] + 1) // 3) if sig else 0
        out = {}
        for name, start, _end in feats:
            out[name.lower()] = {
                "sigpep_len": sig_len,
                "cds_offset": (start - anchor_start) // 3,
                "h_offset": (start - anchor_start) // 3 - sig_len,
            }
        return out
    return None


# ---------------------------------------------------------------------------
# Analyses
# ---------------------------------------------------------------------------

def isolate_of(seqname, isolates):
    if seqname in isolates:
        return seqname
    base = SEQNAME_CAND_RE.sub("", seqname)
    return base if base in isolates else None


def build_constellation(assignments, mutations, genoflu, nc_tables, build_clades):
    """Per isolate x segment label. Priority: genoflu genotype (segment-exact)
    > clade in the record's own winning dataset TSV > a segment-named dataset
    TSV > build clades.json > mutation TSV clade column."""
    isolates = sorted({a["isolate"] for a in assignments.values()})
    segs = sorted(
        {a["segment"] for a in assignments.values() if a["segment"] and a["segment"] != "genome"}
        | set(mutations.keys()),
        key=lambda s: SEGMENT_ORDER.index(s) if s in SEGMENT_ORDER else 99,
    )

    rec_for = {}
    rec_species = {}
    for rec, a in assignments.items():
        if a["segment"] and a["segment"] != "genome":
            rec_for[(a["isolate"], a["segment"])] = rec
            rec_species[(a["isolate"], a["segment"])] = a["record_species"]

    rows = []
    sig_parts = {iso: {} for iso in isolates}
    for iso in isolates:
        for seg in segs:
            label = None
            source = None
            geno = (genoflu.get(iso) or {}).get(seg)
            if geno:
                label, source = geno, "genoflu"
            rec = rec_for.get((iso, seg))
            clade, dataset = None, None
            if rec:
                a = assignments[rec]
                r, ds = nc_row(nc_tables, rec, seg=seg, best_file=a["best_dataset_file"])
                if r:
                    clade = (r.get("clade") or "").strip() or None
                    dataset = ds
            if label is None and clade:
                label, source = clade, "nextclade_clade"
            if label is None:
                bc = (build_clades.get(seg) or {}).get(iso)
                if bc:
                    label, source = bc, "build_clade"
            rs = rec_species.get((iso, seg)) or None
            if label is None and rs:
                label, source = rs, "record_species"
            if label is None:
                for mrow in mutations.get(seg, []):
                    if isolate_of(mrow["seqName"], {iso}) and mrow["clade"]:
                        label, source = mrow["clade"], "mutations_tsv"
                        break
            label = label or "unknown"
            sig_parts[iso][seg] = label
            rows.append({
                "isolate": iso, "segment": seg, "label": label,
                "clade": clade or "", "genoflu_genotype": geno or "",
                "record_species": rs or "",
                "dataset": dataset or "", "label_source": source or "none",
            })

    def signature(iso):
        return "|".join(sig_parts[iso].get(s, "NA") for s in segs)

    sig_counts = {}
    for iso in isolates:
        s = signature(iso)
        sig_counts[s] = sig_counts.get(s, 0) + 1
    modal_sig = max(sig_counts.items(), key=lambda kv: kv[1])[0] if sig_counts else ""

    out_rows = []
    for r in rows:
        iso = r["isolate"]
        sig = signature(iso)
        out_rows.append({**r,
                         "constellation_signature": sig,
                         "modal_signature": modal_sig,
                         "is_novel": str(sig != modal_sig).lower()})
    return out_rows, segs


def build_reassortment(assignments, const_rows):
    per_iso = {}
    for rec, a in assignments.items():
        iso = a["isolate"]
        d = per_iso.setdefault(iso, {"segments": set(), "reasons": []})
        if a["reassortment_suspected"]:
            d["reasons"].append(f"species_assignments flagged record {rec} ({a['segment']})")
        # Per-segment subtype typing (record_species) disagreeing with the
        # isolate's consensus species is direct reassortment evidence —
        # e.g. an h5nx isolate whose ns typed against an h3n2 dataset.
        rs = a.get("record_species")
        if a["segment"] not in ("", "genome") and rs and rs != a["species"]:
            d["reasons"].append(
                f"{a['segment']} typed as {rs} (isolate species {a['species']})")
        if a["segment"]:
            d["segments"].add(a["segment"])

    discord = {}
    for r in const_rows:
        if r["label"] in ("unknown", ""):
            continue
        discord.setdefault(r["segment"], {}).setdefault(r["label"], []).append(r["isolate"])
    modal_label = {}
    for seg, lab_isos in discord.items():
        modal_label[seg] = max(lab_isos.items(), key=lambda kv: len(kv[1]))[0]

    for r in const_rows:
        seg, lab = r["segment"], r["label"]
        if lab in ("unknown", ""):
            continue
        m = modal_label.get(seg)
        if m and lab != m:
            iso_d = per_iso[r["isolate"]]
            reason = f"{seg} label '{lab}' differs from group modal '{m}'"
            if reason not in iso_d["reasons"]:
                iso_d["reasons"].append(reason)

    rows = []
    for iso, d in sorted(per_iso.items()):
        flagged = bool(d["reasons"])
        rows.append({
            "isolate": iso,
            "reassortment_suspected": str(flagged).lower(),
            "reasons": "; ".join(d["reasons"]),
        })
    return rows


def build_di_candidates(assignments, nc_tables, min_deletions, max_coverage):
    rows = []
    for rec, a in sorted(assignments.items()):
        seg = a["segment"] or "untyped"
        # QC must come from the record's own winning dataset — every record
        # appears in every dataset TSV, so a segment-name match would report
        # QC from an alignment the record was not typed against.
        r, ds = nc_row(nc_tables, rec, seg=None,
                       best_file=a["best_dataset_file"])
        if not r:
            continue
        deletions = _to_int(r.get("totalDeletions")) or 0
        frameshifts = (_to_int(r.get("qc.frameShifts.totalFrameShifts"))
                       or _to_int(r.get("totalFrameShifts")) or 0)
        stop_codons = (_to_int(r.get("qc.stopCodons.totalStopCodons"))
                       or _to_int(r.get("totalStopCodons")) or 0)
        failed = (r.get("failedCdses") or "").strip()
        failed_list = [] if failed in ("", "[]", "[]") else [x for x in re.split(r"[,;]", failed) if x]
        coverage = _to_float(r.get("coverage"))
        missing = _to_int(r.get("totalMissing")) or 0

        reasons = []
        if deletions >= min_deletions:
            reasons.append(f"large_deletion_count={deletions}")
        if frameshifts > 0:
            reasons.append(f"frameshifts={frameshifts}")
        if stop_codons > 0:
            reasons.append(f"premature_stop_codons={stop_codons}")
        if failed_list:
            reasons.append(f"failed_cdses={','.join(failed_list)}")
        if coverage is not None and coverage < max_coverage and deletions > 0:
            reasons.append(f"low_coverage={coverage:.2f}_with_deletions")
        rows.append({
            "isolate": a["isolate"], "record": rec, "segment": seg,
            "dataset": ds or "", "coverage": coverage if coverage is not None else "",
            "total_deletions": deletions, "total_frameshifts": frameshifts,
            "total_stop_codons": stop_codons, "failed_cdses": ",".join(failed_list),
            "total_missing": missing,
            "di_candidate": str(bool(reasons)).lower(),
            "reasons": "; ".join(reasons),
        })
    return rows


def match_markers(catalog, mutations, assignments, ha_offsets, ha_sigpep_fallback):
    """Join query-isolate AA mutations against the marker catalogue.

    Coordinate handling:
      - segment_cds markers match GENE:pos directly (gene positions are
        already per-protein for unsplit genes).
      - HA sub-CDS labels (SigPep/HA1/HA2) are translated to HA0/mature-HA
        coordinates via the genemap offsets; 'HA:' labels (avian) translate
        to mature-HA numbering by subtracting the signal peptide.
      - n1_numbering NA markers: exact match, then +/-1 window match flagged
        coordinate_caveat (subtype numbering can differ by indels).
      - motif markers are skipped here (they are not point substitutions);
        cleavage-motif reporting comes from the build's cleavage-site.json.
    """
    isolates = {a["isolate"] for a in assignments.values()}
    species = next(iter({a["species"] for a in assignments.values() if a["species"]}), "")
    point_markers = [m for m in catalog if m.get("type", "point") == "point"]

    by_gene = {}
    for m in point_markers:
        by_gene.setdefault(str(m.get("gene", "")).upper(), []).append(m)

    hits = []
    for seg, rows in sorted(mutations.items()):
        for mrow in rows:
            iso = isolate_of(mrow["seqName"], isolates)
            if not iso:
                continue
            for gene, ref, pos, alt in mrow["aa"]:
                parent = SUBCDS_TO_PARENT.get(gene.lower(), gene.upper())
                if parent == "HA" and gene.upper() not in ("HA",):
                    cds_pos = pos + (ha_offsets or {}).get(gene.lower(), {}).get("cds_offset", 0)
                    h_pos = pos + (ha_offsets or {}).get(gene.lower(), {}).get("h_offset", 0)
                    conv_ok = bool(ha_offsets)
                elif gene.upper() == "HA":
                    sig_len = (ha_offsets or {}).get("ha1", {}).get("sigpep_len", ha_sigpep_fallback)
                    cds_pos, h_pos = pos, pos - sig_len
                    conv_ok = True
                else:
                    cds_pos = h_pos = pos
                    conv_ok = True

                for mk in by_gene.get(parent, []):
                    cs = (mk.get("coordinate_system") or "segment_cds").lower()
                    caveat = ""
                    if cs == "h_numbering":
                        if parent != "HA":
                            continue
                        if not conv_ok and gene.upper() != "HA":
                            continue
                        if not conv_ok:
                            caveat = "no genemap: signal-peptide offset assumed 0"
                        target = h_pos
                        expected = int(mk["position"])
                        if target != expected:
                            continue
                    elif cs == "n1_numbering":
                        expected = int(mk["position"])
                        equiv = (mk.get("equivalents") or {})
                        eq_pos = equiv.get(species)
                        if eq_pos is not None and pos == int(eq_pos):
                            pass
                        elif pos == expected:
                            caveat = "N1-numbering position matched directly"
                        elif abs(pos - expected) == 1:
                            caveat = "within +/-1 of N1-numbered position — verify subtype numbering"
                        else:
                            continue
                    else:  # segment_cds
                        expected = int(mk["position"])
                        if cds_pos != expected:
                            continue
                    if alt.upper() not in [str(a).upper() for a in (mk.get("alt") or [])]:
                        continue
                    scope = mk.get("subtype_scope") or []
                    if scope and species and species not in scope:
                        caveat = (caveat + "; " if caveat else "") + \
                            f"marker characterized in {'/'.join(scope)}, this run is {species}"
                    hits.append({
                        "sample": iso, "segment": seg, "gene": parent,
                        "mutation": f"{gene}:{ref}{pos}{alt}",
                        "marker_position": expected,
                        "category": mk.get("category", ""),
                        "phenotype": mk.get("phenotype", ""),
                        "drug": mk.get("drug", ""),
                        "effect": mk.get("effect", ""),
                        "confidence": mk.get("confidence", ""),
                        "coordinate_system": cs,
                        "coordinate_caveat": caveat,
                        "source": mk.get("source", ""),
                        "refs": ";".join(mk.get("refs") or []),
                    })
    return hits


def build_group_diversity(groups, mutations, isolates):
    iso_group = {i: groups.get(i) for i in isolates}
    by_group = {}
    for iso, g in iso_group.items():
        if g:
            by_group.setdefault(g, []).append(iso)

    mut_sets = {}
    for seg, rows in mutations.items():
        for mrow in rows:
            iso = isolate_of(mrow["seqName"], isolates)
            if iso:
                mut_sets.setdefault((iso, seg), set()).update(
                    f"{g}:{r}{p}{a}" for g, r, p, a in mrow["aa"])

    rows = []
    for gid, members in sorted(by_group.items()):
        members = sorted(set(members))
        segs = sorted({s for (i, s) in mut_sets if i in members},
                      key=lambda s: SEGMENT_ORDER.index(s) if s in SEGMENT_ORDER else 99)
        if len(members) < 2:
            for seg in segs:
                rows.append({
                    "group_id": gid, "segment": seg, "members": ",".join(members),
                    "n_isolates": 1, "n_shared": "", "n_private": "",
                    "mean_pairwise": "", "shared_mutations": "",
                    "private_mutations": "",
                    "note": "single isolate in group — no within-group diversity",
                })
            continue
        for seg in segs:
            sets = [mut_sets.get((i, seg), set()) for i in members]
            shared = set.intersection(*sets) if sets else set()
            union = set.union(*sets) if sets else set()
            private = {i: sorted(s - shared) for i, s in zip(members, sets)}
            pairs = [len(a ^ b) for i, a in enumerate(sets) for b in sets[i + 1:]]
            rows.append({
                "group_id": gid, "segment": seg, "members": ",".join(members),
                "n_isolates": len(members), "n_shared": len(shared),
                "n_private": len(union - shared),
                "mean_pairwise": round(sum(pairs) / len(pairs), 2) if pairs else 0,
                "shared_mutations": ",".join(sorted(shared)),
                "private_mutations": "; ".join(f"{i}:[{','.join(v)}]" for i, v in private.items() if v),
                "note": "",
            })
    return rows


# ---------------------------------------------------------------------------

def write_tsv(path, rows, fields):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mutations", nargs="*", default=[],
                    help="mutations_{seg}.tsv files from SEGMENT_MUTATIONS")
    ap.add_argument("--assignments", required=True, help="species_assignments.tsv")
    ap.add_argument("--nextclade", nargs="*", default=[],
                    help="classification nextclade TSVs (per-record QC + clade)")
    ap.add_argument("--metadata", default=None, help="input metadata.tsv (group_id)")
    ap.add_argument("--genoflu", default=None, help="genoflu_results.tsv (avian)")
    ap.add_argument("--build-dir", default=None,
                    help="nextstrain build results dir (clades.json / genemap.gff)")
    ap.add_argument("--markers", required=True, help="database/flu_markers.yml")
    ap.add_argument("--species", default="")
    ap.add_argument("--pathogen", default="",
                    help="meta.pathogen — scopes species_assignments rows to this group")
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--report", default="signatures_report.txt")
    ap.add_argument("--di-min-deletion", type=int, default=15,
                    help="total deletions threshold for the DI screen")
    ap.add_argument("--di-max-coverage", type=float, default=0.7,
                    help="coverage below which deletions also count toward DI")
    ap.add_argument("--ha-signal-peptide", type=int, default=16,
                    help="fallback HA signal peptide length when no genemap is available")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    assignments = load_assignments(args.assignments, pathogen=args.pathogen or None,
                                   species=args.species or None)
    isolates = {a["isolate"] for a in assignments.values()}
    mutations = load_mutations(
        [p for p in args.mutations if not Path(p).name.startswith("NO_FILE")]
    )
    nc_tables = load_nextclade_tables(args.nextclade)
    groups = load_metadata_groups(args.metadata)
    genoflu = load_genoflu(args.genoflu)
    build_clades = load_build_clades(args.build_dir)
    ha_offsets = load_ha_offsets(args.build_dir)

    catalog = []
    if yaml is not None:
        cat = yaml.safe_load(Path(args.markers).read_text())
        catalog = (cat or {}).get("markers") or []
    else:
        log("WARNING: pyyaml not available — marker screening skipped")

    log(f"{len(isolates)} isolates, {len(mutations)} segments, "
        f"{len(nc_tables)} nextclade tables, {len(catalog)} catalog markers")

    const_rows, segs = build_constellation(assignments, mutations, genoflu, nc_tables, build_clades)
    write_tsv(outdir / "constellation.tsv", const_rows,
              ["isolate", "segment", "label", "clade", "genoflu_genotype",
               "record_species", "dataset", "label_source",
               "constellation_signature", "modal_signature", "is_novel"])

    re_rows = build_reassortment(assignments, const_rows)
    write_tsv(outdir / "reassortment_flags.tsv", re_rows,
              ["isolate", "reassortment_suspected", "reasons"])

    di_rows = build_di_candidates(assignments, nc_tables,
                                  args.di_min_deletion, args.di_max_coverage)
    write_tsv(outdir / "di_candidates.tsv", di_rows,
              ["isolate", "record", "segment", "dataset", "coverage",
               "total_deletions", "total_frameshifts", "total_stop_codons",
               "failed_cdses", "total_missing", "di_candidate", "reasons"])

    marker_rows = match_markers(catalog, mutations, assignments, ha_offsets,
                                args.ha_signal_peptide)
    write_tsv(outdir / "markers.tsv", marker_rows,
              ["sample", "segment", "gene", "mutation", "marker_position",
               "category", "phenotype", "drug", "effect", "confidence",
               "coordinate_system", "coordinate_caveat", "source", "refs"])

    gd_rows = build_group_diversity(groups, mutations, isolates)
    write_tsv(outdir / "group_diversity.tsv", gd_rows,
              ["group_id", "segment", "members", "n_isolates", "n_shared",
               "n_private", "mean_pairwise", "shared_mutations",
               "private_mutations", "note"])

    n_di = sum(1 for r in di_rows if r["di_candidate"] == "true")
    n_re = sum(1 for r in re_rows if r["reassortment_suspected"] == "true")
    n_nov = sum(1 for r in const_rows if r["is_novel"] == "true" and r["segment"] == segs[0])
    report = f"""flu_segment_signatures report
=============================
species:            {args.species or '(unset)'}
isolates:           {len(isolates)}
segments:           {len(segs)} ({', '.join(segs)})
constellation rows: {len(const_rows)}  (novel constellations: {n_nov})
reassortment flags: {n_re} / {len(re_rows)}
DI candidates:      {n_di} / {len(di_rows)}
marker hits:        {len(marker_rows)}
group diversity:    {len(gd_rows)} group x segment rows ({len({r['group_id'] for r in gd_rows})} groups)
"""
    (outdir / args.report).write_text(report)
    log("done — outputs in " + str(outdir))


if __name__ == "__main__":
    main()
