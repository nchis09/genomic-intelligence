process DB_QUERY {
    tag "${sample_id}"

    input:
    tuple val(sample_id), val(bio_parent_dir)
    val outdir
    val db_url

    output:
    tuple val(sample_id), val("${projectDir}/${outdir}/data_query/${sample_id}"), emit: data_query_tuple

    script:
    """
    export PYTHONPATH="${projectDir}"
    python3 -m intelligence_engine.data_engine.sql_querying.bioinformatics_query \
        --bioinformatics-dir "${bio_parent_dir}" \
        --output-dir "${projectDir}/${outdir}/data_query" \
        --sample "${sample_id}" \
        --db-url "${db_url}"
    """
}
