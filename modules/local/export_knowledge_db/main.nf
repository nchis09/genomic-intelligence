/*
 * Local module: EXPORT_KNOWLEDGE_DB
 *
 * Copies every table in the shared PostgreSQL knowledge warehouse into a
 * single, portable DuckDB file, so PATHOGEN_IDENTIFICATION can query a
 * self-contained local file instead of connecting to the live shared
 * server. Gated (via a collected value channel input) on every species'
 * BUILD_KNOWLEDGE_DB call having finished.
 */

process EXPORT_KNOWLEDGE_DB {
    label 'process_low'
    cache false

    conda "${projectDir}/envs/pgirl_knowledge.yml"

    input:
    val(db_host)
    val(db_port)
    val(_gate)

    output:
    path("knowledge_warehouse.duckdb"), emit: duckdb_dump

    script:
    """
    export PATH="\${CONDA_PREFIX:+\$CONDA_PREFIX/bin:}\$PATH"
    python3 ${projectDir}/bin/export_knowledge_db.py \\
        --host ${db_host} \\
        --port ${db_port} \\
        --output knowledge_warehouse.duckdb
    """
}
