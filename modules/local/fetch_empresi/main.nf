/*
 * Local module: FETCH_EMPRESI
 *
 * Download FAO EMPRES-i avian influenza outbreak events and add them to the
 * existing epi_data directory for an avian H5 species.
 */
process FETCH_EMPRESI {
    tag "$meta.id"
    label 'process_low'

    conda '/Users/christianndekezi/anaconda3/envs/pgirl_literature'
    container null

    input:
    tuple val(meta), val(species), path(epi_raw)

    output:
    tuple val(meta), path("epi_data"), emit: epi_raw
    tuple val(meta), path("empresi_search_summary.tsv"), emit: search_summary

    when:
    task.ext.when == null || task.ext.when

    script:
    def min_arg = (params.epi_min_year != null) ? "--min-year ${params.epi_min_year}" : ""
    def max_arg = (params.epi_max_year != null) ? "--max-year ${params.epi_max_year}" : ""
    """
    mkdir -p epi_data
    if [ -d "${epi_raw}" ]; then
        cp -r ${epi_raw}/* epi_data/ 2>/dev/null || true
    fi

    /Users/christianndekezi/anaconda3/envs/pgirl_literature/bin/python3 ${projectDir}/bin/fetch_empresi.py \
        --species "${species}" \
        ${min_arg} \
        ${max_arg} \
        --outdir .
    """
}
