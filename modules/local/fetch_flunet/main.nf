/*
 * Local module: FETCH_FLUNET
 *
 * Download strain-specific influenza virological surveillance from the WHO
 * FluMart public API and reshape it into a long-format CSV for the
 * knowledge warehouse.
 */
process FETCH_FLUNET {
    tag "$meta.id"
    label 'process_low'

    conda '/Users/christianndekezi/anaconda3/envs/pgirl_literature'
    container null

    input:
    tuple val(meta), val(species)

    output:
    tuple val(meta), path("epi_data"),               emit: epi_raw
    tuple val(meta), path("flunet_search_summary.tsv"), emit: search_summary

    when:
    task.ext.when == null || task.ext.when

    script:
    def min_arg = (params.epi_min_year != null) ? "--min-year ${params.epi_min_year}" : ""
    def max_arg = (params.epi_max_year != null) ? "--max-year ${params.epi_max_year}" : ""
    """
    /Users/christianndekezi/anaconda3/envs/pgirl_literature/bin/python3 ${projectDir}/bin/fetch_flunet.py \
        --species "${species}" \
        ${min_arg} \
        ${max_arg} \
        --outdir .
    """
}
