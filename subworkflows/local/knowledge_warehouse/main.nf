/*
 * Subworkflow: KNOWLEDGE_WAREHOUSE
 *
 * Collect key pipeline outputs and build a per-run PostgreSQL knowledge warehouse.
 * This subworkflow starts the shared PostgreSQL instance internally, so the
 * main workflow does not need to expose the DB lifecycle node to nf-metro.
 *
 * Once every species' ingestion is done, the whole warehouse is exported to
 * a single, portable DuckDB file (EXPORT_KNOWLEDGE_DB) and the shared
 * Postgres server is stopped immediately afterwards (STOP_KNOWLEDGE_DB) --
 * downstream stages like PATHOGEN_IDENTIFICATION read that DuckDB file
 * instead of holding a live connection open for the rest of the run.
 *
 * Input:  ch_warehousable - channel: [ val(meta), path(species_assignments),
 *                            path(metadata_tsv), path(epi_raw_dir),
 *                            path(epi_search_summary), path(bioinformatics_results),
 *                            path(uniprotr_results), path(extractr_results),
 *                            path(rbioapi_results), path(auspice_json), path(iqtree),
 *                            path(query_data_files), path(hmm_files),
 *                            val(evidence_qc_ready), val(db_host), val(db_port) ]
 *                            Pre-joined by the caller (see ebola_workflow/main.nf),
 *                            including the shared Postgres host/port to be used,
 *                            so this subworkflow's DAG entry has a single inbound
 *                            edge -- keeping db_host/db_port as separate `take:`
 *                            params would give this subworkflow multiple inbound
 *                            edges and can crash nf-metro's layout with a
 *                            CurveInvariantError.
 * Output: knowledge_db    - channel: [ val(meta), path(knowledge_warehouse) ]
 *         duckdb_dump     - channel: path(knowledge_warehouse.duckdb)
 */

include { BUILD_KNOWLEDGE_DB  } from '../../../modules/local/build_knowledge_db/main'
include { START_KNOWLEDGE_DB  } from '../../../modules/local/start_knowledge_db/main'
include { EXPORT_KNOWLEDGE_DB } from '../../../modules/local/export_knowledge_db/main'
include { STOP_KNOWLEDGE_DB   } from '../../../modules/local/stop_knowledge_db/main'

workflow KNOWLEDGE_WAREHOUSE {
    take:
    ch_warehousable // channel: [ meta, assignments, metadata, epi_dir, epi_summary, bioinfo_dir, uniprotr_dir, extractr_dir, rbioapi_dir, auspice_json, iqtree, query_data, hmm, evidence_qc_ready, db_host, db_port ]

    main:
    ch_schema = Channel.fromPath("${projectDir}/database/knowledge_schema.sql").first()
    ch_views  = Channel.fromPath("${projectDir}/database/knowledge_views.sql").first()
    ch_warehousable = ch_warehousable.combine(ch_schema).combine(ch_views)

    // Start the shared DB once every species' warehousable input is ready
    ch_db_gate = ch_warehousable.collect().map { true }
    ch_db_host = Channel.value(params.kw_db_host)
    ch_db_port = Channel.value(params.kw_db_port)
    START_KNOWLEDGE_DB(ch_db_host, ch_db_port, ch_db_gate)
    ch_db_ready = START_KNOWLEDGE_DB.out.ready

    ch_build_input = ch_warehousable.map { meta, assignments, metadata, epi_dir, epi_summary, bioinfo_dir, uniprotr_dir, extractr_dir, rbioapi_dir, auspice_json, iqtree, query_data, hmm_files, evidence_qc_ready, _host, _port, schema, views ->
        [ meta, assignments, metadata, epi_dir, epi_summary, bioinfo_dir, uniprotr_dir, extractr_dir, rbioapi_dir, auspice_json, iqtree, query_data, hmm_files, evidence_qc_ready, schema, views ]
    }
    ch_build_host = ch_warehousable.map { meta, assignments, metadata, epi_dir, epi_summary, bioinfo_dir, uniprotr_dir, extractr_dir, rbioapi_dir, auspice_json, iqtree, query_data, hmm_files, evidence_qc_ready, host, port, schema, views -> host }
    ch_build_port = ch_warehousable.map { meta, assignments, metadata, epi_dir, epi_summary, bioinfo_dir, uniprotr_dir, extractr_dir, rbioapi_dir, auspice_json, iqtree, query_data, hmm_files, evidence_qc_ready, _host, port, schema, views -> port }

    // Gate the actual ingestion on the DB being ready
    ch_build_input_gated = ch_build_input.combine(ch_db_ready).map { it -> it[0..-2] }
    ch_build_host_gated  = ch_build_host.combine(ch_db_ready).map { it -> it[0] }
    ch_build_port_gated  = ch_build_port.combine(ch_db_ready).map { it -> it[0] }

    BUILD_KNOWLEDGE_DB(ch_build_input_gated, ch_build_host_gated, ch_build_port_gated)

    // Once every species has been ingested, export the whole warehouse to a
    // portable DuckDB file, then stop the shared server immediately -- no
    // downstream stage needs a live connection to it any more.
    ch_export_gate = BUILD_KNOWLEDGE_DB.out.knowledge_db.collect().map { true }
    EXPORT_KNOWLEDGE_DB(ch_db_host, ch_db_port, ch_export_gate)
    STOP_KNOWLEDGE_DB(ch_db_host, ch_db_port, EXPORT_KNOWLEDGE_DB.out.duckdb_dump.collect().map { true })

    emit:
    knowledge_db = BUILD_KNOWLEDGE_DB.out.knowledge_db     // channel: [ meta, path ]
    mqc_summary  = BUILD_KNOWLEDGE_DB.out.mqc_summary      // channel: [ meta, path ]
    duckdb_dump  = EXPORT_KNOWLEDGE_DB.out.duckdb_dump      // channel: path(knowledge_warehouse.duckdb)
}
