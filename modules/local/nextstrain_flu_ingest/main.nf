/*
 * Local module: NEXTSTRAIN_FLU_INGEST
 *
 * Run the Nextstrain seasonal-flu ingest-open Snakemake workflow
 * (GenSpectrum/LAPIS public source) for ONE influenza lineage
 * (h1n1pdm, h3n2 or vic — the only lineages ingest-open supports).
 * Fetches all eight genome segments so the staged background covers the
 * full genome even though the build layer uses only HA + NA.
 *
 * Yamagata and generic Flu B have no independent public ingest source —
 * those groups are routed onto the vic lineage (see
 * params.nextclade_species_overrides) and share this vic background.
 *
 * Staged results stay under data/seasonal_flu/ingest-open/results/{lineage}/
 * so --skip_nextstrain_ingest can reuse them (same pattern as
 * NEXTSTRAIN_AVIAN_INGEST).
 */

process NEXTSTRAIN_FLU_INGEST {
    tag "$meta.id"
    label 'process_medium'

    conda "${projectDir}/envs/pgirl_nextstrain.yml"

    input:
    val(meta)

    output:
    tuple val(meta), path('background/metadata.tsv'), path('background/*.fasta'), emit: background

    when:
    !params.skip_nextstrain_ingest

    script:
    def nextstrain_dir = "${projectDir}/data/seasonal_flu"
    def lineage = meta.species   // h1n1pdm | h3n2 | vic (see overrides)
    """
    # Ensure the ingest Python helpers are available in the activated conda env.
    PY_BIN=\$(dirname \$(which python))
    if [ ! -e "\${PY_BIN}/python3" ]; then
        ln -sf "\${PY_BIN}/python" "\${PY_BIN}/python3"
    fi
    export PATH="\${PY_BIN}:\${PATH}"

    INGEST_DIR="${nextstrain_dir}/ingest-open"
    RESULTS_DIR="\${INGEST_DIR}/results/${lineage}"

    # Restrict the fetch to this group's lineage. `segments` stays at the
    # upstream default (all 8) — full-genome scope per project design.
    cat > pgirl_ingest_override.yaml <<EOF
lineages:
  - ${lineage}
EOF

    # Remove stale results so Snakemake fetches fresh data (mirrors the
    # Ebola/avian ingest pattern).
    rm -rf "\${RESULTS_DIR}"

    # Run from the ingest dir: the Snakefile reads an optional config.yaml
    # from the working directory — we use --configfile instead.
    cd "\${INGEST_DIR}"

    snakemake \
        --snakefile Snakefile \
        --cores ${task.cpus} \
        --configfile defaults/config.yaml "\$OLDPWD/pgirl_ingest_override.yaml" \
        --rerun-incomplete \
        --nolock

    cd "\$OLDPWD"
    mkdir -p background
    cp "\${RESULTS_DIR}/metadata.tsv" background/metadata.tsv
    cp "\${RESULTS_DIR}"/*.fasta background/
    """
}
