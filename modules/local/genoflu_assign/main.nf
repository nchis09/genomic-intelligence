/*
 * Local module: GENOFLU_ASSIGN
 *
 * Run the vendored GenoFLU-multi genotyper on an avian influenza species
 * group BEFORE ingest/build, so a resolved subtype narrows the background
 * fetch and the build subtype (h5nx -> h5n1 when cleanly assigned).
 *
 * GenoFLU is an HPAI H5N1 constellation genotyper: a clean `Genotype`
 * call implies the isolate is H5N1; "Not assigned" stays h5nx (the NA
 * segment tree remains the evidence for other N subtypes).
 *
 * genoflu-multi groups FASTA records by header id into per-strain record
 * sets, so records are renamed to their ISOLATE id first (group metadata
 * carries accession=record id / strain=isolate id per row).
 */

process GENOFLU_ASSIGN {
    tag "$meta.id"
    label 'process_medium'

    conda "${projectDir}/envs/pgirl_nextstrain.yml"

    input:
    tuple val(meta), path(fasta), path(metadata)

    output:
    tuple val(meta), path("genoflu_results.tsv"), emit: results

    when:
    task.ext.when == null || task.ext.when

    script:
    def genoflu_dir = "${projectDir}/data/avian_flu/ingest/vendored-GenoFLU-multi"
    """
    mkdir -p genoflu_input

    # Rename each record's header to its isolate id so genoflu-multi treats
    # the isolate's segments as one strain. Mapping comes from the group
    # metadata (accession -> strain), not the record-name suffix.
    python3 - "${fasta}" "${metadata}" genoflu_input/queries.fasta <<'PYEOF'
import csv, sys

fasta_path, metadata_path, out_path = sys.argv[1:4]
acc_to_iso = {}
with open(metadata_path) as fh:
    for row in csv.DictReader(fh, delimiter="\\t"):
        acc = (row.get("accession") or "").strip()
        iso = (row.get("strain") or acc).strip()
        if acc:
            acc_to_iso[acc] = iso

name = None
with open(fasta_path) as fh, open(out_path, "w") as out:
    for line in fh:
        if line.startswith(">"):
            rid = line[1:].split()[0]
            out.write(f">{acc_to_iso.get(rid, rid)}\\n")
        else:
            out.write(line)
PYEOF

    # -i: genotype incomplete genomes too (partial submissions still get
    # any-call evidence). GenoFLU failures never kill the pipeline — the
    # group just stays h5nx.
    python "${genoflu_dir}/bin/genoflu-multi.py" -f genoflu_input -i || true

    if [ -s genoflu_input/results/results.tsv ]; then
        cp genoflu_input/results/results.tsv genoflu_results.tsv
    else
        : > genoflu_results.tsv
    fi
    """
}
