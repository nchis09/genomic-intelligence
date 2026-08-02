/*
 * Subworkflow: KNOWLEDGE_WAREHOUSE
 *
 * Collect key pipeline outputs and build a per-run PostgreSQL knowledge warehouse.
 *
 * Input:  ch_warehousable - channel: [ val(meta), path(species_assignments),
 *                            path(metadata_tsv), path(epi_raw_dir),
 *                            path(epi_search_summary), path(bioinformatics_results),
 *                            path(uniprotr_results), path(extractr_results),
 *                            path(rbioapi_results) ]
 *                            Pre-joined by the caller (see ebola_workflow/main.nf)
 *                            so this subworkflow's DAG entry has a single
 *                            inbound edge.
 * Output: knowledge_db    - channel: [ val(meta), path(knowledge_warehouse) ]
 */

include { BUILD_KNOWLEDGE_DB } from '../../../modules/local/build_knowledge_db/main'

workflow KNOWLEDGE_WAREHOUSE {
    take:
    ch_warehousable // channel: [ meta, assignments, metadata, epi_dir, epi_summary, bioinfo_dir, uniprotr_dir, extractr_dir, rbioapi_dir ]

    main:
    BUILD_KNOWLEDGE_DB(ch_warehousable)

    emit:
    knowledge_db = BUILD_KNOWLEDGE_DB.out.knowledge_db  // channel: [ meta, path ]
    mqc_summary  = BUILD_KNOWLEDGE_DB.out.mqc_summary   // channel: [ meta, path ]
}
