/*
 * Local module: NEXTSTRAIN_EBOLA_INGEST
 *
 * Run the Nextstrain Ebola ingest Snakemake workflow to download and curate
 * background sequences/metadata for a given species. This process runs on
 * every pipeline invocation so the background data stays current.
 */

process NEXTSTRAIN_EBOLA_INGEST {
    tag "$meta.id"
    label 'process_medium'

    conda "${projectDir}/envs/pgirl_nextstrain.yml"

    input:
    val(meta)

    output:
    tuple val(meta), path('background/sequences.fasta'), path('background/metadata.tsv'), emit: background

    when:
    !params.skip_nextstrain_ingest

    script:
    def species        = (meta.species ?: meta.id).replaceAll(/^.*_/, '')
    def nextstrain_dir = "${projectDir}/data/nextstrain_ebola"
    """
    # Ensure the ingest Python helpers are available in the activated conda env.
    # Some conda envs expose only 'python', but the upstream ingest scripts call
    # 'python3' in their shebang, so make sure python3 resolves to the env Python.
    PY_BIN=\$(dirname \$(which python))
    if [ ! -e "\${PY_BIN}/python3" ]; then
        ln -sf "\${PY_BIN}/python" "\${PY_BIN}/python3"
    fi
    export PATH="\${PY_BIN}:\${PATH}"
    python3 -c "import pandas" 2>/dev/null || python3 -m pip install --quiet pandas
    python3 -c "import bio" 2>/dev/null || python3 -m pip install --quiet bio

    INGEST_DIR="${nextstrain_dir}/ingest"
    SPP_DATA="\${INGEST_DIR}/data/${species}"

    # Remove stale local ingest files so the Snakemake workflow fetches and
    # curates the latest background sequences/metadata.
    rm -rf "\${SPP_DATA}"

    snakemake \
        --snakefile \${INGEST_DIR}/Snakefile \
        --cores ${task.cpus} \
        --config species=[${species}] \
        --rerun-incomplete \
        --nolock \
        data/${species}/sequences.fasta data/${species}/metadata.tsv

    mkdir -p background
    cp "\${SPP_DATA}/sequences.fasta" background/sequences.fasta
    cp "\${SPP_DATA}/metadata.tsv" background/metadata.tsv
    """
}
