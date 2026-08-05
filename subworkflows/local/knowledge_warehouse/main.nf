/*
 * Subworkflow: KNOWLEDGE_WAREHOUSE
 *
 * Collect key pipeline outputs and build a per-run PostgreSQL knowledge warehouse.
 *
 * Input:  ch_warehousable - channel: [ val(meta), path(species_assignments),
 *                            path(metadata_tsv), path(epi_raw_dir),
 *                            path(epi_search_summary), path(bioinformatics_results),
 *                            path(uniprotr_results), path(extractr_results),
 *                            path(rbioapi_results), path(auspice_json), path(iqtree),
 *                            path(query_data_files), path(hmm_files), val(db_host), val(db_port) ]
 *                            Pre-joined by the caller (see ebola_workflow/main.nf),
 *                            including the shared Postgres host/port (pointing
 *                            at the instance START_KNOWLEDGE_DB started once at
 *                            pipeline launch), so this subworkflow's DAG entry
 *                            has a single inbound edge -- keeping db_host/
 *                            db_port as separate `take:` params would give
 *                            this subworkflow multiple inbound edges and can
 *                            crash nf-metro's layout with a CurveInvariantError.
 * Output: knowledge_db    - channel: [ val(meta), path(knowledge_warehouse) ]
 */

include { BUILD_KNOWLEDGE_DB } from '../../../modules/local/build_knowledge_db/main'

workflow KNOWLEDGE_WAREHOUSE {
    take:
    ch_warehousable // channel: [ meta, assignments, metadata, epi_dir, epi_summary, bioinfo_dir, uniprotr_dir, extractr_dir, rbioapi_dir, auspice_json, iqtree, query_data, hmm, db_host, db_port ]

    main:
    ch_build_input = ch_warehousable.map { meta, assignments, metadata, epi_dir, epi_summary, bioinfo_dir, uniprotr_dir, extractr_dir, rbioapi_dir, auspice_json, iqtree, query_data, hmm_files, _host, _port ->
        [ meta, assignments, metadata, epi_dir, epi_summary, bioinfo_dir, uniprotr_dir, extractr_dir, rbioapi_dir, auspice_json, iqtree, query_data, hmm_files ]
    }
    ch_db_host = ch_warehousable.map { meta, assignments, metadata, epi_dir, epi_summary, bioinfo_dir, uniprotr_dir, extractr_dir, rbioapi_dir, auspice_json, iqtree, query_data, hmm_files, host, _port -> host }
    ch_db_port = ch_warehousable.map { meta, assignments, metadata, epi_dir, epi_summary, bioinfo_dir, uniprotr_dir, extractr_dir, rbioapi_dir, auspice_json, iqtree, query_data, hmm_files, _host, port -> port }

    BUILD_KNOWLEDGE_DB(ch_build_input, ch_db_host, ch_db_port)

    emit:
    knowledge_db = BUILD_KNOWLEDGE_DB.out.knowledge_db  // channel: [ meta, path ]
    mqc_summary  = BUILD_KNOWLEDGE_DB.out.mqc_summary   // channel: [ meta, path ]
}
