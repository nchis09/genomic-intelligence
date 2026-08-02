/*
 * Subworkflow: KNOWLEDGE_WAREHOUSE
 *
 * Collect key pipeline outputs and build a per-run PostgreSQL knowledge warehouse.
 *
 * Input:  ch_species_data        - channel: [ val(meta), path(fasta), path(metadata) ]
 *         ch_species_assignments - path: species_assignments.tsv (broadcast)
 *         ch_epi_raw             - channel: [ val(meta), path(epi_data_dir) ]
 *         ch_epi_search_summary  - channel: [ val(meta), path(rhdx_search_results.tsv) ]
 * Output: knowledge_db           - channel: [ val(meta), path(knowledge_warehouse) ]
 */

include { BUILD_KNOWLEDGE_DB } from '../../../modules/local/build_knowledge_db/main'

workflow KNOWLEDGE_WAREHOUSE {
    take:
    ch_species_data        // channel: [ val(meta), path(fasta), path(metadata) ]
    ch_species_assignments // path: species_assignments.tsv
    ch_epi_raw             // channel: [ val(meta), path(epi_data_dir) ]
    ch_epi_search_summary  // channel: [ val(meta), path(tsv) ]

    main:
    // Re-key inputs by pathogen so epidemiological data (one entry per pathogen)
    // can be broadcast to every species group for that pathogen.
    def no_file_meta    = file('NO_FILE_metadata')
    def no_file_epi     = file('NO_FILE_epi')
    def no_file_summary = file('NO_FILE_epi_summary')

    ch_species_by_pathogen = ch_species_data
        .map { meta, fasta, metadata ->
            def meta_tsv = (metadata && !metadata.name.startsWith('NO_FILE')) ? metadata : no_file_meta
            [ meta.pathogen, meta, meta_tsv ]
        }

    ch_epi_by_pathogen = ch_epi_raw
        .map { meta, epi_dir ->
            def epi = (epi_dir && !epi_dir.name.startsWith('NO_FILE')) ? epi_dir : no_file_epi
            [ meta.pathogen, epi ]
        }

    ch_epi_summary_by_pathogen = ch_epi_search_summary
        .map { meta, tsv ->
            def summary = (tsv && !tsv.name.startsWith('NO_FILE')) ? tsv : no_file_summary
            [ meta.pathogen, summary ]
        }

    // Join metadata with epidemiological data by pathogen, then broadcast the
    // single species_assignments file to every species group.
    ch_warehousable = ch_species_by_pathogen
        .join(ch_epi_by_pathogen, by: 0, remainder: true)
        .join(ch_epi_summary_by_pathogen, by: 0, remainder: true)
        .map { pathogen, meta, metadata, epi_dir, epi_summary ->
            [ meta, metadata ?: no_file_meta, epi_dir ?: no_file_epi, epi_summary ?: no_file_summary ]
        }
        .combine(ch_species_assignments)
        .map { meta, metadata, epi_dir, epi_summary, assignments ->
            [ meta, assignments, metadata, epi_dir, epi_summary ]
        }

    BUILD_KNOWLEDGE_DB(ch_warehousable)

    emit:
    knowledge_db = BUILD_KNOWLEDGE_DB.out.knowledge_db  // channel: [ meta, path ]
}
