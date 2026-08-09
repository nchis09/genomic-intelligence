/*
 * Local module: EPI_DATA_PLACEHOLDER
 *
 * Lightweight placeholder process so the empty EPIDEMIOLOGICAL_DATA
 * subworkflow is visible in the nf-metro DAG until real epidemiological
 * processes are added.
 */

process EPI_DATA_PLACEHOLDER {
    tag 'epi_data_placeholder'
    label 'process_low'

    output:
    val(true), emit: ready

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    echo "Epidemiological data placeholder (to be populated)"
    """
}
