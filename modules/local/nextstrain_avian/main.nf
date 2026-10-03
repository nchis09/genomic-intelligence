/*
 * Local module: NEXTSTRAIN_AVIAN
 *
 * Run the Nextstrain avian-flu segment-focused Snakemake workflow which:
 *   1. Merges locally-ingested NCBI background sequences (config `inputs`)
 *   2. Injects the group's query sequences as `additional_inputs`, one
 *      FASTA per segment to build. Segment identity comes from
 *      species_assignments.tsv (the Nextclade dataset each record scored
 *      against), never from record names. Records are renamed to their
 *      ISOLATE id — one strain name per isolate per segment tree.
 *   3. 'genome'-typed records (untypable at screening — usually the NA,
 *      which has no avian Nextclade dataset) become NA-build candidates
 *      for isolates lacking a typed NA. Duplicates get _gN suffixes;
 *      augur align + min_length drops impostor segments.
 *   4. Runs the augur pipeline per segment (all typed segments + na
 *      candidates) through to Auspice export — the NA tree is what
 *      resolves the `x` in `h5nx` phylogenetically.
 */

process NEXTSTRAIN_AVIAN {
    tag "$meta.id"
    label 'process_high'

    conda "${projectDir}/envs/pgirl_nextstrain.yml"

    input:
    tuple val(meta), path(fasta), path(metadata), path(aligned), path(bg_metadata), path(bg_fastas), path(assignments)

    output:
    tuple val(meta), path("${meta.species}/auspice/*.json")   , emit: auspice
    tuple val(meta), path("${meta.species}/results/")         , emit: results_dir

    when:
    task.ext.when == null || task.ext.when

    script:
    def subtype      = task.ext.subtype ?: meta.species
    def time         = task.ext.time    ?: 'all-time'
    def nextstrain_dir = "${projectDir}/data/avian_flu"
    """
    # Snakemake rules invoke `python3 scripts/*.py` — pin python3 to the
    # task env's interpreter so it sees augur/pandas (system python lacks
    # them). Same trick as the ingest modules.
    PY_BIN=\$(dirname \$(which python))
    if [ ! -e "\${PY_BIN}/python3" ]; then
        ln -sf "\${PY_BIN}/python" "\${PY_BIN}/python3"
    fi
    export PATH="\${PY_BIN}:\${PATH}"

    NEXTSTRAIN_DIR="${nextstrain_dir}"

    # --- Per-segment query extraction keyed on species_assignments.
    # ha records come from the winning dataset's aligned output (aligned to
    # the HA reference); every other segment is taken from the raw group
    # fasta — augur aligns each to its own segment reference.
    python3 - "${assignments}" "${meta.query_samples}" "${fasta}" "${aligned}" <<'PYEOF'
import csv, os, sys
from collections import defaultdict

assignments_tsv, group_samples, group_fasta, aligned_fasta = sys.argv[1:5]
SEG_ORDER = ["pb2", "pb1", "pa", "ha", "np", "na", "mp", "ns"]

group_ids = set(group_samples.split(','))
rows = []
with open(assignments_tsv) as fh:
    for r in csv.DictReader(fh, delimiter="\\t"):
        if r["sample"] in group_ids:
            rows.append(r)
isolates = sorted({r["isolate"] for r in rows})

best = {}   # (isolate, seg) -> (score, sample)
for r in rows:
    seg = r["segment"]
    key = (r["isolate"], seg)
    try:
        score = float(r["qc_score"])
    except ValueError:
        score = float("inf")
    if key not in best or score < best[key][0]:
        best[key] = (score, r["sample"])

seg_out = defaultdict(list)   # seg -> [(tip_name, record_id)]
for (iso, seg), (score, sample) in best.items():
    if seg in SEG_ORDER:
        seg_out[seg].append((iso, sample))

# 'genome'-typed records (screening couldn't tether them — typically the
# NA, which has no avian Nextclade dataset) feed the NA build for isolates
# without a typed NA. Extras get _gN suffixes; augur's align + min_length
# drop impostor segments.
by_iso = defaultdict(list)
for r in rows:
    by_iso[r["isolate"]].append(r)
for iso, rs in by_iso.items():
    if not any(seg == "na" for (i, seg) in best if i == iso):
        extras = [r["sample"] for r in rs if r["segment"] == "genome"]
        for i, sample in enumerate(extras):
            tip = iso if i == 0 else f"{iso}_g{i+1}"
            seg_out["na"].append((tip, sample))

def read_fasta(path):
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line)
    if name is not None:
        yield name, "".join(seq)

built = []
for seg in SEG_ORDER:
    if seg not in seg_out:
        continue
    src = aligned_fasta if seg == "ha" else group_fasta
    if not os.path.exists(src) or os.path.getsize(src) == 0:
        sys.stderr.write(f"NEXTSTRAIN_AVIAN: no source fasta for segment {seg}\\n")
        continue
    want = dict((rec, tip) for tip, rec in seg_out[seg])
    n = 0
    with open(f"query_{seg}.fasta", "w") as out:
        for name, seq in read_fasta(src):
            if name in want:
                out.write(f">{want[name]}\\n{seq}\\n")
                n += 1
    if n:
        built.append(seg)
    sys.stderr.write(f"NEXTSTRAIN_AVIAN: query_{seg}.fasta <- {n} records\\n")

with open("built_segments.txt", "w") as fh:
    fh.write(" ".join(built) + "\\n")
with open("isolates.txt", "w") as fh:
    fh.write("\\n".join(isolates) + "\\n")
# All emitted tip names (isolates plus any _gN genome-candidate aliases) —
# each needs a metadata row and a force-include entry to survive the build.
# Dedupe: an isolate appears once per segment but augur metadata requires
# unique ids.
with open("tips.txt", "w") as fh:
    seen = set()
    for seg in built:
        for tip, _rec in seg_out[seg]:
            if tip not in seen:
                seen.add(tip)
                fh.write(tip + "\\n")
PYEOF

    BUILT_SEGMENTS=\$(cat built_segments.txt)

    # Force-include the query tips in the build's subsampling step:
    # --include is a union applied after all filtering, so it protects
    # query rows from exclude_where / min-date / min-length dropping them
    # when the query metadata is sparse.
    cat "${nextstrain_dir}/config/${subtype}/include_strains_${subtype}_${time}.txt" \
        tips.txt | sort -u > include_strains.txt

    # --- Query metadata: one row per emitted tip name. augur merge keys
    # inputs on `strain` (the avian merge rules don't support a per-input
    # id_field), so `strain` must equal the FASTA header — the isolate id
    # (or its _gN genome-candidate alias, which inherits the isolate's row).
    python3 - "${metadata}" tips.txt query_metadata.tsv <<'PYEOF'
import csv, datetime, sys

today = datetime.date.today().isoformat()
with open(sys.argv[2]) as fh:
    tips = [l.strip() for l in fh if l.strip()]
with open(sys.argv[1]) as fh:
    rows = list(csv.DictReader(fh, delimiter="\\t"))
by_iso = {r.get("strain") or r.get("accession") or "": r for r in rows}
out_rows = []
for tip in tips:
    iso = tip.split("_g")[0] if "_g" in tip else tip
    src = by_iso.get(iso)
    if src is None:
        continue
    r = dict(src)
    r["query_strain"] = iso
    r["accession"] = tip
    r["strain"] = tip
    r.setdefault("is_query", "true")
    if not (r.get("date") or "").strip():
        r["date"] = today
    out_rows.append(r)
header = list(rows[0].keys()) if rows else ["strain", "is_query", "date"]
for extra in ("query_strain", "accession", "strain", "is_query", "date"):
    if extra not in header:
        header.append(extra)
with open(sys.argv[3], "w") as fh:
    w = csv.DictWriter(fh, fieldnames=header, delimiter="\\t", extrasaction="ignore")
    w.writeheader()
    w.writerows(out_rows)
PYEOF

    # Background sequences live beside the staged metadata in this work dir
    # (bg_fastas are staged next to it as sequences_<segment>.fasta).
    BG_DIR=\$(dirname \$(realpath ${bg_metadata}))

    SEG_YAML=""
    QUERY_SEQ_YAML=""
    for SEG in \${BUILT_SEGMENTS}; do
        SEG_YAML="\${SEG_YAML}      - \${SEG}
"
        if [ -s "query_\${SEG}.fasta" ]; then
            QUERY_SEQ_YAML="\${QUERY_SEQ_YAML}      \${SEG}: \$(realpath query_\${SEG}.fasta)
"
        fi
    done

    # The subtype_query keeps every background strain of the subtype's H5
    # window AND every query row (queries carry no 'subtype' column, so
    # without the is_query clause they would be filtered out before the
    # tree step). Broad h5nx covers all H5 serotypes; a GenoFLU-resolved
    # subtype narrows the window to itself. Augur coerces the is_query
    # flag to a boolean, so compare with the bare True literal (verified:
    # 'true'/'True' string forms don't match).
    if [ "${subtype}" = "h5nx" ]; then
        SUBTYPE_LIST="'h5n1', 'h5n2', 'h5n3', 'h5n4', 'h5n5', 'h5n6', 'h5n7', 'h5n8', 'h5n9'"
    else
        SUBTYPE_LIST="'${subtype}'"
    fi

    cat > pgirl_override.yaml <<EOF
builds:
  - subtype:
      - ${subtype}
    segment:
\${SEG_YAML}
    time:
      - ${time}

inputs:
  - name: ncbi_local
    metadata: \$(realpath ${bg_metadata})
    sequences: \${BG_DIR}/sequences_{segment}.fasta

additional_inputs:
  - name: pgirl_query
    metadata: \$(realpath query_metadata.tsv)
    sequences:
\${QUERY_SEQ_YAML}

subtype_query:
  "${subtype}/*/*": "subtype in [\${SUBTYPE_LIST}] or is_query == True"

include_strains: \$(realpath include_strains.txt)
EOF

    mkdir -p ${meta.species}

    snakemake \\
        --snakefile \${NEXTSTRAIN_DIR}/segment-focused/Snakefile \\
        --cores ${task.cpus} \\
        --directory ${meta.species} \\
        --configfile \$(realpath pgirl_override.yaml) \\
        --rerun-incomplete \\
        --nolock
    """
}
