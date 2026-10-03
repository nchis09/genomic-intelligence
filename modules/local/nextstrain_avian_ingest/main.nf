/*
 * Local module: NEXTSTRAIN_AVIAN_INGEST
 *
 * Run the Nextstrain avian-flu ingest Snakemake workflow (NCBI GenBank
 * source) to download and curate background sequences/metadata. Unlike the
 * Ebola ingest (per-species), this fetches ALL configured Influenza A
 * serotypes in one shot -- the H5 filters below widen upstream's default
 * H5N1/North-America/since-2023 slice to global coverage back to 1996 so
 * the phylogenetic background covers the full diversity of H5Nx.
 *
 * Output is per-segment: background/sequences_{pb2,pb1,pa,ha,np,na,mp,ns}.fasta
 * plus a shared background/metadata.tsv. This process runs on every
 * pipeline invocation so the background data stays current (same pattern
 * as NEXTSTRAIN_EBOLA_INGEST).
 */

process NEXTSTRAIN_AVIAN_INGEST {
    tag "$meta.id"
    label 'process_medium'

    conda "${projectDir}/envs/pgirl_nextstrain.yml"

    input:
    val(meta)

    output:
    tuple val(meta), path('background/metadata.tsv'), path('background/sequences_*.fasta'), emit: background

    when:
    !params.skip_nextstrain_ingest

    script:
    def nextstrain_dir = "${projectDir}/data/avian_flu"
    // GenoFLU may have narrowed the group to a resolved subtype (h5n1) —
    // fetch background for just that serotype. Unresolved h5nx keeps the
    // full H5N1..H5N9 window.
    def serotypes = meta.species ==~ /h5n[1-9]/
        ? '"' + meta.species.toUpperCase() + '"'
        : '"H5N1" "H5N2" "H5N3" "H5N4" "H5N5" "H5N6" "H5N7" "H5N8" "H5N9"'
    """
    # Ensure the ingest Python helpers are available in the activated conda env.
    PY_BIN=\$(dirname \$(which python))
    if [ ! -e "\${PY_BIN}/python3" ]; then
        ln -sf "\${PY_BIN}/python" "\${PY_BIN}/python3"
    fi
    export PATH="\${PY_BIN}:\${PATH}"

    INGEST_DIR="${nextstrain_dir}/ingest"
    RESULTS_DIR="\${INGEST_DIR}/ncbi/results"

    # Widen the upstream defaults: no region restriction and all-time
    # history (h5nx builds start at 1996). Serotype window narrows to the
    # resolved subtype when GenoFLU assigned one.
    cat > pgirl_ingest_override.yaml <<EOF
ncbi_virus_min_collection_date: "1996"
ncbi_virus_filters:
  - 'Serotype_s:(${serotypes})'
EOF

    # Remove stale results so Snakemake fetches fresh data (mirrors the
    # Ebola ingest's rm -rf of its per-species data dir).
    rm -rf "\${RESULTS_DIR}"

    # Run from the ingest dir: the Snakefile's relative `configfile:` and
    # include paths resolve against the working directory, and --configfile
    # is variadic (one flag, multiple files; later files win).
    cd "\${INGEST_DIR}"

    snakemake \
        --snakefile Snakefile \
        --cores ${task.cpus} \
        --configfile build-configs/ncbi/defaults/config.yaml \$OLDPWD/pgirl_ingest_override.yaml \
        --rerun-incomplete \
        --nolock \
        ingest_ncbi

    cd "\$OLDPWD"
    mkdir -p background
    cp "\${RESULTS_DIR}/metadata.tsv" background/metadata.tsv
    cp "\${RESULTS_DIR}"/sequences_*.fasta background/
    """
}
