/*
 * Local module: PATHOGEN_MUTATION_PROFILE
 *
 * For each detected species, compute the overall mutation profile
 * (query vs background) from a portable DuckDB export of the knowledge
 * warehouse. Produces TSVs and MultiQC custom-content tables only;
 * no PNGs — the Shiny dashboard will create the visualisations.
 */

process PATHOGEN_MUTATION_PROFILE {
    tag "$meta.id"
    label 'process_low'

    conda "${projectDir}/envs/pgirl_pathogen_genomics.yml"

    input:
    tuple val(meta), path(duckdb_file)

    output:
    path "mutation_profile/*.tsv",         emit: tsv
    path "mutation_profile/mqc/*_mqc.tsv", emit: mqc_tsv

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args   = task.ext.args   ?: ''
    """
    export PATH="\$CONDA_PREFIX/bin:\$PATH"

    Rscript ${projectDir}/bin/pathogen_mutation_profile.R \
        --duckdb ${duckdb_file} \
        --species ${meta.species} \
        --run-id ${prefix} \
        --outdir . \
        ${args}
    """
}
