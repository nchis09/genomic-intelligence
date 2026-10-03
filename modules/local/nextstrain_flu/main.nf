/*
 * Local module: NEXTSTRAIN_FLU
 *
 * Run the vendored nextstrain/seasonal-flu Snakemake workflow for one
 * influenza lineage group (h1n1pdm, h3n2 or vic — vic also covers yam and
 * generic B queries, which have no dedicated ingest source).
 *
 *   1. GenSpectrum/LAPIS background + query sequences are merged via the
 *      workflow's `inputs`/`additional_inputs` mechanism into
 *      data/{lineage}/metadata.tsv + data/{lineage}/{segment}.fasta
 *   2. Per-segment query FASTAs are extracted from the lineage's Nextclade
 *      aligned outputs. Segment identity comes from species_assignments.tsv
 *      (the dataset each record scored against), never from record names.
 *      Records are renamed to their ISOLATE id so the same strain name is
 *      shared across all segment trees — required by the vendored
 *      sanitize_trees step which intersects leaf names across segments.
 *   3. The "{lineage}_pgirl" build runs every segment shared by all query
 *      isolates (sanitize_trees would drop an isolate absent from any
 *      segment tree), with a 6-year temporal subsample, force-including
 *      every query isolate.
 *   4. Segments missing vendored config assets (internals: pb2…ns) get
 *      reference.fasta/genemap.gff fetched from their Nextclade dataset.
 *
 * Auspice output: {species}/auspice/{lineage}_pgirl_{segment}.json
 *
 * Vendored repo: nextstrain/seasonal-flu @ a1e704c77c894cb868b9e162da6ec5220de5efc7
 */

process NEXTSTRAIN_FLU {
    tag "$meta.id"
    label 'process_high'

    conda "${projectDir}/envs/pgirl_nextstrain.yml"

    input:
    tuple val(meta), path(fasta), path(metadata), path(aligned_fastas), val(aligned_segs), path(bg_metadata), path(bg_fastas), path(assignments)

    output:
    // Optional: a group whose records extract no usable segment sequence
    // skips the build — the auspice glob is then empty.
    tuple val(meta), path("${meta.species}/auspice/*.json", optional: true), emit: auspice
    tuple val(meta), path("${meta.species}")               , emit: results_dir

    when:
    task.ext.when == null || task.ext.when

    script:
    def lineage = task.ext.lineage ?: meta.species   // h1n1pdm | h3n2 | vic
    def nextstrain_dir = "${projectDir}/data/seasonal_flu"
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
    WORK_DIR="${meta.species}"
    mkdir -p "\${WORK_DIR}"

    # The workflow's helper rules run `python3 scripts/*.py` and read
    # `config/...` assets relative to --directory. The clades/subclades
    # download rules also write into `config/{lineage}/{segment}/` — copy
    # (not symlink) so downloads land in the work dir, not the vendored repo.
    cp -r "\${NEXTSTRAIN_DIR}/config"  "\${WORK_DIR}/config"
    cp -r "\${NEXTSTRAIN_DIR}/scripts" "\${WORK_DIR}/scripts"

    # --- Per-segment query extraction.
    # The parallel `aligned_segs`/`aligned_fastas` inputs give a segment ->
    # aligned-fasta map; assignments give record -> (isolate, segment,
    # qc_score). One record per (isolate, segment) is kept (best qc_score);
    # 'genome'-typed records (untethered junk hits) feed nothing.
    printf '%s\\n' ${aligned_segs.join(' ')} > "\${WORK_DIR}/aligned_segs.txt"
    printf '%s\\n' ${aligned_fastas.join(' ')} > "\${WORK_DIR}/aligned_files.txt"
    paste "\${WORK_DIR}/aligned_segs.txt" "\${WORK_DIR}/aligned_files.txt" > "\${WORK_DIR}/aligned_map.tsv"

    python3 - "${assignments}" "${meta.query_samples}" "\${WORK_DIR}/aligned_map.tsv" "\${WORK_DIR}" <<'PYEOF'
import csv, os, sys
from collections import defaultdict

assignments_tsv, group_samples, map_tsv, out_dir = sys.argv[1:5]
SEG_ORDER = ["pb2", "pb1", "pa", "ha", "np", "na", "mp", "ns"]

group_ids = set(group_samples.split(','))
rows = []
with open(assignments_tsv) as fh:
    for r in csv.DictReader(fh, delimiter="\\t"):
        if r["sample"] in group_ids:
            rows.append(r)
isolates = sorted({r["isolate"] for r in rows})
typed = defaultdict(set)            # isolate -> {segments}
best  = {}                          # (isolate, seg) -> (score, sample)
for r in rows:
    seg = r["segment"]
    if seg in ("genome", ""):
        continue
    key = (r["isolate"], seg)
    try:
        score = float(r["qc_score"])
    except ValueError:
        score = float("inf")
    typed[r["isolate"]].add(seg)
    if key not in best or score < best[key][0]:
        if key in best:
            sys.stderr.write(
                f"NEXTSTRAIN_FLU: duplicate {seg} records for {r['isolate']}; "
                f"keeping {r['sample']} (qc={score})\\n")
        best[key] = (score, r["sample"])

seg_file = {}
with open(map_tsv) as fh:
    for line in fh:
        parts = line.rstrip("\\n").split("\\t")
        if len(parts) == 2 and parts[0] and parts[1]:
            seg_file[parts[0]] = parts[1]

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

# Segments every isolate covers — sanitize_trees intersects leaf names
# across all built segments, so a segment missing for ANY isolate would
# silently drop that isolate from every tree. ha must always be present:
# the workflow's clade machinery annotates from the ha tree.
shared = set.intersection(*[typed[i] for i in isolates]) if isolates else set()
if "ha" not in shared:
    shared.add("ha")
built = [s for s in SEG_ORDER if s in shared]

for seg in built:
    src = seg_file.get(seg)
    if not src or not os.path.exists(src):
        sys.stderr.write(f"NEXTSTRAIN_FLU: no aligned fasta for segment {seg}\\n")
        continue
    want = {}   # record_id -> isolate
    for iso in isolates:
        rec = best.get((iso, seg))
        if rec:
            want[rec[1]] = iso
    n = 0
    out_path = os.path.join(out_dir, f"query_{seg}.fasta")
    with open(out_path, "w") as out:
        for name, seq in read_fasta(src):
            if name in want:
                out.write(f">{want[name]}\\n{seq}\\n")
                n += 1
    sys.stderr.write(f"NEXTSTRAIN_FLU: query_{seg}.fasta <- {n} records\\n")

with open(os.path.join(out_dir, "built_segments.txt"), "w") as fh:
    fh.write(" ".join(built) + "\\n")
with open(os.path.join(out_dir, "isolates.txt"), "w") as fh:
    fh.write("\\n".join(isolates) + "\\n")
PYEOF

    # Internal segments (pb2…ns) have no vendored reference/genemap — fetch
    # them from the segment's own Nextclade dataset so augur align/translate
    # get per-segment assets.
    BUILT_SEGMENTS=\$(cat "\${WORK_DIR}/built_segments.txt")
    for SEG in \${BUILT_SEGMENTS}; do
        REF_DIR="\${WORK_DIR}/config/${lineage}/\${SEG}"
        if [ ! -s "\${REF_DIR}/reference.fasta" ]; then
            mkdir -p "\${REF_DIR}"
            nextclade dataset get --name "nextstrain/flu/${lineage}/\${SEG}" --output-dir "\${WORK_DIR}/ds_\${SEG}" >/dev/null 2>&1
            if [ -s "\${WORK_DIR}/ds_\${SEG}/reference.fasta" ]; then
                cp "\${WORK_DIR}/ds_\${SEG}/reference.fasta" "\${REF_DIR}/reference.fasta"
                if [ -s "\${WORK_DIR}/ds_\${SEG}/genome_annotation.gff3" ]; then
                    # augur's GFF reader only accepts 'gene'/'source'/'region'
                    # rows (augur/io/sequences.py load_features), while the
                    # Nextclade dataset ships 'CDS' rows keyed by Name=. Emit
                    # flat `gene` features with `gene_name` like the vendored
                    # genemaps. Spliced ORFs (multi-CDS: M2, NEP, PA-X) can't
                    # be represented flat — they are dropped and must stay out
                    # of GENES in core.smk.
                    python3 - "\${WORK_DIR}/ds_\${SEG}/genome_annotation.gff3" "\${REF_DIR}/genemap.gff" <<'PYEOF'
import sys
src, dst = sys.argv[1:3]
seqid, region, cds = None, None, []
count = {}
with open(src) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.rstrip("\\n").split("\\t")
        if len(f) < 9:
            continue
        if f[2] == "region":
            region = (f[0], f[3], f[4])
            seqid = f[0]
        elif f[2] == "CDS":
            seqid = seqid or f[0]
            name = next((kv[5:] for kv in f[8].split(";") if kv.startswith("Name=")), None)
            if name:
                count[name] = count.get(name, 0) + 1
                cds.append((f[3], f[4], f[6], name))
if seqid is None or region is None:
    sys.exit("could not parse nextclade gff3 " + src)
with open(dst, "w") as out:
    out.write("\\t".join([seqid, "feature", "region", region[1], region[2], ".", "+", ".", 'gene_name="nuc"']) + "\\n")
    for start, end, strand, name in cds:
        if count[name] == 1:
            out.write("\\t".join([seqid, "feature", "gene", start, end, ".", strand, ".", f'gene_name="{name}"']) + "\\n")
PYEOF
                fi
            else
                echo "NEXTSTRAIN_FLU [${meta.id}]: could not fetch reference for segment \${SEG}" >&2
            fi
        fi
    done

    # --- Query metadata: one row per isolate. augur merge keys the query
    # input on `accession`, which must equal the FASTA header — the isolate
    # id. The metadata rows are already fanned out per record with
    # `strain` = isolate, so dedupe on that key.
    python3 - "${metadata}" "\${WORK_DIR}/isolates.txt" "\${WORK_DIR}/query_metadata.tsv" <<'PYEOF'
import csv, datetime, sys

today = datetime.date.today().isoformat()
with open(sys.argv[2]) as fh:
    isolates = {l.strip() for l in fh if l.strip()}
with open(sys.argv[1]) as fh:
    rows = list(csv.DictReader(fh, delimiter="\\t"))
seen, out_rows = set(), []
for r in rows:
    iso = r.get("strain") or r.get("accession") or ""
    if iso not in isolates or iso in seen:
        continue
    seen.add(iso)
    # augur merge emits `strain` as the shared id column — a query-side
    # column of the same name would collide, so the isolate goes to
    # `query_strain` and `accession` (the merge key) only.
    r["query_strain"] = iso
    r["accession"] = iso
    r.pop("strain", None)
    r.setdefault("is_query", "true")
    if not (r.get("date") or "").strip():
        r["date"] = today
    out_rows.append(r)
header = [c for c in (list(rows[0].keys()) if rows else []) if c != "strain"]
if not header:
    header = ["accession", "query_strain", "is_query", "date"]
for extra in ("query_strain", "accession", "is_query", "date"):
    if extra not in header:
        header.append(extra)
with open(sys.argv[3], "w") as fh:
    w = csv.DictWriter(fh, fieldnames=header, delimiter="\\t", extrasaction="ignore")
    w.writeheader()
    w.writerows(out_rows)
PYEOF

    # Force-include: every query isolate plus the lineage's reference strains.
    cat "\${WORK_DIR}/config/${lineage}/reference_strains.txt" \
        "\${WORK_DIR}/isolates.txt" | sort -u > "\${WORK_DIR}/include_strains.txt"

    # A group can legitimately contain zero records with extracted segment
    # sequence. Typing still succeeded; skip the build rather than emit a
    # background-only tree.
    N_QUERY=\$(cat "\${WORK_DIR}"/query_*.fasta 2>/dev/null | grep -c '^>' || true)
    if [ "\${N_QUERY:-0}" -eq 0 ]; then
        echo "NEXTSTRAIN_FLU [${meta.id}]: group has no extracted segment records — skipping seasonal-flu build" >&2
        mkdir -p "\${WORK_DIR}/auspice"
        echo "No query segment records extracted for lineage ${lineage}; build skipped." \
            > "\${WORK_DIR}/auspice/NO_QUERY_SEGMENTS.txt"
        exit 0
    fi

    BG_DIR=\$(dirname \$(realpath "${bg_metadata}"))

    # 6-year temporal resolution (user decision). Build a literal filter
    # string — the workflow .format()s filters with build_params, and a
    # brace-free string passes through untouched.
    MIN_DATE=\$(python3 -c "import datetime; print((datetime.date.today()-datetime.timedelta(days=6*365)).isoformat())")

    QUERY_SEQ_YAML=""
    for SEG in \${BUILT_SEGMENTS}; do
        if [ -s "\${WORK_DIR}/query_\${SEG}.fasta" ]; then
            QUERY_SEQ_YAML="\${QUERY_SEQ_YAML}        \${SEG}: \$(realpath \${WORK_DIR}/query_\${SEG}.fasta)
"
        fi
    done

    SEG_YAML=""
    for SEG in \${BUILT_SEGMENTS}; do
        SEG_YAML="\${SEG_YAML}  - \${SEG}
"
    done

    cat > pgirl_builds.yaml <<EOF
segments:
\${SEG_YAML}
lat-longs: config/lat_longs.tsv

inputs:
  - name: genspectrum
    lineage: ${lineage}
    metadata: \$(realpath ${bg_metadata})
    id_field: strain
    sequences: \${BG_DIR}/{segment}.fasta

additional_inputs:
  - name: pgirl_query
    lineage: ${lineage}
    metadata: \$(realpath "\${WORK_DIR}/query_metadata.tsv")
    id_field: accession
    sequences:
\${QUERY_SEQ_YAML}

builds:
  ${lineage}_pgirl:
    lineage: ${lineage}
    reference: "config/${lineage}/{segment}/reference.fasta"
    annotation: "config/${lineage}/{segment}/genemap.gff"
    # No subclades key: subclade TSVs only exist for ha/na, and the
    # download_subclades rule would curl a None URL for internal segments.
    # Subclade labels fall back to the ha clades.json via import_clades.
    clades: "config/${lineage}/ha/clades.tsv"
    auspice_config: "config/${lineage}/auspice_config.json"
    include: include_strains.txt
    subsamples:
      global:
        filters: "--min-date \${MIN_DATE} --include include_strains.txt --exclude config/${lineage}/outliers.txt --group-by region year month --subsample-max-sequences 3000"
    min_date: "\${MIN_DATE}"
EOF

    snakemake \\
        --snakefile "\${NEXTSTRAIN_DIR}/Snakefile" \\
        --cores ${task.cpus} \\
        --directory "\${WORK_DIR}" \\
        --configfile \$(realpath pgirl_builds.yaml) \\
        --rerun-incomplete \\
        --nolock
    """
}
