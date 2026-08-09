/*
 * Local module: STOP_KNOWLEDGE_DB
 *
 * Dumps the shared knowledge warehouse to a single, portable .sql file and
 * stops the server. Gated (via a collected value channel input) on every
 * species' BUILD_KNOWLEDGE_DB and QUERY_KNOWLEDGE_DB call having finished,
 * so nothing is still writing/reading when the server goes down.
 */

process STOP_KNOWLEDGE_DB {
    label 'process_single'

    conda "${projectDir}/envs/pgirl_knowledge.yml"

    input:
    val(db_host)
    val(db_port)
    val(_gate)

    output:
    path("*_genomic_intelligence.sql"), emit: dump

    script:
    def data_dir = "${params.outdir}/knowledge_warehouse/_shared_pg_data"
    """
    export PATH="\${CONDA_PREFIX:+\$CONDA_PREFIX/bin:}\$PATH"
    python3 ${projectDir}/bin/stop_shared_db.py \\
        --data-dir ${data_dir} \\
        --host ${db_host} \\
        --port ${db_port} \\
        --outdir . \\
        --prefix run
    """
}
