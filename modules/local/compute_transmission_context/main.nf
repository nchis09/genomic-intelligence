/*
 * Local module: COMPUTE_TRANSMISSION_CONTEXT
 *
 * Pre-compute Transmission & Spread tables (epidemiological burden, spatial
 * hotspots, anomaly alerts and transmission potential) for one species from
 * the exported DuckDB knowledge warehouse.
 */

process COMPUTE_TRANSMISSION_CONTEXT {
    tag "$meta.id"
    label 'process_low'

    conda "${projectDir}/envs/pgirl_knowledge.yml"

    input:
    tuple val(meta), path(duckdb_file)

    output:
    tuple val(meta), path("transmission_context/${meta.species}"), emit: tsv

    script:
    def species = meta.species ?: 'ebov'
    """
    export PATH="\${CONDA_PREFIX:+\$CONDA_PREFIX/bin:}\$PATH"

    python3 ${projectDir}/bin/compute_transmission_context.py \
        --db ${duckdb_file} \
        --species ${species} \
        --outdir transmission_context \
        --si-mean ${params.si_mean} \
        --si-sd ${params.si_sd} \
        --cluster-threshold ${params.cluster_threshold}

    python3 ${projectDir}/bin/model_transmission_potential.py \
        --input transmission_context/${species}/transmission_potential.tsv \
        --outdir transmission_context/${species}
    """
}
