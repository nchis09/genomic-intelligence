process EPI_QUERY {
    tag "${sample_id}"

    input:
    tuple val(sample_id), val(bio_dir), val(db_query_dir)
    val outdir
    val use_llm

    output:
    tuple val(sample_id), val("${projectDir}/${outdir}/data_query/${sample_id}"), emit: data_query_tuple

    script:
    def llm_flag = use_llm ? "" : "--no-llm"
    """
    export PYTHONPATH="${projectDir}"
    ${params.python} -m intelligence_engine.data_engine.online_querying.epi_query_engine \
        --bio-output "${bio_dir}/bio_output.json" \
        --db-query-results "${db_query_dir}/db_query_results.json" \
        --output "${projectDir}/${outdir}/data_query/${sample_id}/epi_output.json" \
        ${llm_flag}
    """
}
