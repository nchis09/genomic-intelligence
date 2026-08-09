/*
 * Local module: START_KNOWLEDGE_DB
 *
 * Starts (or reuses) the single, pipeline-lifetime shared PostgreSQL
 * instance used by every species' BUILD_KNOWLEDGE_DB call and the FIGURES
 * stage's QUERY_KNOWLEDGE_DB call. Runs once, before any species-specific
 * work begins. The server is started as an independent daemon bound to a
 * fixed data directory outside this task's own work dir, so it keeps
 * running after this (quick) task completes -- STOP_KNOWLEDGE_DB tears it
 * down once, at the very end of the run.
 */

process START_KNOWLEDGE_DB {
    label 'process_single'

    conda "${projectDir}/envs/pgirl_knowledge.yml"

    input:
    val(db_host)
    val(db_port)
    val(_gate)

    output:
    val(true), emit: ready

    script:
    def data_dir = "${params.outdir}/knowledge_warehouse/_shared_pg_data"
    def log_file = "${params.outdir}/knowledge_warehouse/_shared_postgres.log"
    """
    export PATH="\${CONDA_PREFIX:+\$CONDA_PREFIX/bin:}\$PATH"
    mkdir -p "${params.outdir}/knowledge_warehouse"
    python3 ${projectDir}/bin/start_shared_db.py \\
        --data-dir ${data_dir} \\
        --log-file ${log_file} \\
        --host ${db_host} \\
        --port ${db_port}
    """
}
