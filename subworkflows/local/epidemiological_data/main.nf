/*
 * Subworkflow: EPIDEMIOLOGICAL_DATA
 *
 * Placeholder subworkflow for future epidemiological analyses
 * (e.g. case mapping, transmission dynamics, growth rate estimation).
 * Visible in the DAG via a placeholder process until populated.
 */

include { EPI_DATA_PLACEHOLDER } from '../../../modules/local/epi_data_placeholder/main'

workflow EPIDEMIOLOGICAL_DATA {
    // no inputs yet

    main:
    EPI_DATA_PLACEHOLDER()

    // no outputs yet
}
