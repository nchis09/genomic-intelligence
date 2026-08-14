/*
 * Subworkflow: PATHOGEN_IDENTIFICATION_WF
 *
 * Runs the species identification statistical analysis against a portable
 * DuckDB export of the knowledge warehouse, for each detected species.
 * Wrapped in its own subworkflow (rather than calling the
 * PATHOGEN_IDENTIFICATION process directly from the main workflow) so
 * Nextflow's own DAG exporter gives it a native subgraph boundary.
 *
 * Input: ch_pathogen_id - channel: [ val(meta), path(duckdb_file) ]
 *                          A plain file-based input (the portable DuckDB
 *                          export of the knowledge warehouse -- see
 *                          EXPORT_KNOWLEDGE_DB) rather than a live
 *                          db_host/db_port connection.
 * Output: tsv - channel: path(*.tsv)
 */

include { PATHOGEN_IDENTIFICATION } from '../../../modules/local/pathogen_identification/main'

workflow PATHOGEN_IDENTIFICATION_WF {
    take:
    ch_pathogen_id // channel: [ meta, duckdb_file ]

    main:
    PATHOGEN_IDENTIFICATION(ch_pathogen_id)

    emit:
    tsv = PATHOGEN_IDENTIFICATION.out.tsv
}
