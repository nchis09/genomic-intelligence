/*
 * Local module: PATHOGEN_IDENTIFICATION
 *
 * For each detected species, run an R script (species identification
 * analysis) against a portable DuckDB export of the knowledge warehouse
 * (see EXPORT_KNOWLEDGE_DB) to produce per-species statistical TSV tables.
 * Reading a local file instead of holding a live connection to the shared
 * Postgres server means this stage has a single, plain file-based input --
 * no db_host/db_port channel plumbing needed -- and the shared server can
 * already be stopped by the time this runs.
 */

process PATHOGEN_IDENTIFICATION {
    tag "$meta.id"
    label 'process_low'

    conda "${projectDir}/envs/pgirl_pathogen_id.yml"

    input:
    tuple val(meta), path(duckdb_file)

    output:
    path "species_identification/*.tsv", emit: tsv
    path "species_identification/mqc/*_mqc.tsv", emit: mqc_tsv

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    export PATH="\$CONDA_PREFIX/bin:\$PATH"

    Rscript ${projectDir}/bin/pathogen_identification.R \\
        --species ${meta.species} \\
        --run-id ${meta.id} \\
        --duckdb-file ${duckdb_file} \\
        --outdir .
    """
}
