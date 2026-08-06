/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { CLASSIFICATION         } from '../subworkflows/local/classification/main'
include { PATHOGEN_ROUTER        } from '../subworkflows/local/pathogen_router/main'
include { KNOWLEDGE_WAREHOUSE    } from '../subworkflows/local/knowledge_warehouse/main'
include { START_KNOWLEDGE_DB     } from '../modules/local/start_knowledge_db/main'
include { STOP_KNOWLEDGE_DB      } from '../modules/local/stop_knowledge_db/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow GENOMIC_INTELLIGENCE {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:


    //
    // Shared PostgreSQL connection details. They are set once the DB is
    // started (after per-species analysis outputs are ready) or default to
    // configured values if the knowledge warehouse is skipped.
    //
    def db_host
    def db_port

    //
    // SUBWORKFLOW: Nextclade classification + species assignment
    //
    // Nextclade datasets to screen against (configurable via params)
    ch_datasets = channel.fromList(params.nextclade_datasets)

    // Pass full samplesheet — CLASSIFICATION extracts fasta internally,
    // runs Nextclade, assigns species, and emits per-species groups
    CLASSIFICATION(ch_samplesheet, ch_datasets)

    //
    // SUBWORKFLOW: Pathogen router (dispatches each species group to its workflow)
    //
    // CLASSIFICATION auto-detected pathogen/species → route to the matching
    // pathogen-specific workflow (currently only Ebola is registered). Each
    // pathogen workflow runs the per-species analysis modules (bioinformatics,
    // phenotype annotation, epidemiological data) it needs. The shared DB is
    // then built in the main workflow (see below). Groups whose
    // pathogen has no registered workflow are skipped with a warning; see
    // PATHOGEN_ROUTER.out.unsupported for the summary file.
    //
    // All Nextclade JSON outputs (one per sample x dataset run) and the
    // species_assignments.tsv (sample -> winning dataset) are broadcast to
    // every species group below (collect()/single-path outputs behave as
    // value channels) for use by each pathogen workflow's own bioinformatics
    // analysis stage.
    ch_nextclade_json_all = CLASSIFICATION.out.json
        .map { _meta, json -> json }
        .collect()

    PATHOGEN_ROUTER(
        CLASSIFICATION.out.species_groups,
        ch_nextclade_json_all,
        CLASSIFICATION.out.assignments
    )

    //
    // Shared knowledge-warehouse: start the DB once the pathogen analyses
    // have produced per-species data, ingest the data, and stop the DB once
    // ingestion is complete.
    //
    def ch_knowledge_db = channel.empty()
    def ch_knowledge_db_summary = channel.empty()

    if (!params.skip_knowledge_warehouse) {
        START_KNOWLEDGE_DB(
            channel.value(params.kw_db_host),
            channel.value(params.kw_db_port),
            PATHOGEN_ROUTER.out.kw_input.collect()
        )
        db_host = START_KNOWLEDGE_DB.out.ready.map { params.kw_db_host }
        db_port = START_KNOWLEDGE_DB.out.ready.map { params.kw_db_port }

        KNOWLEDGE_WAREHOUSE(
            PATHOGEN_ROUTER.out.kw_input.combine(db_host).combine(db_port)
        )
        ch_knowledge_db = KNOWLEDGE_WAREHOUSE.out.knowledge_db
        ch_knowledge_db_summary = KNOWLEDGE_WAREHOUSE.out.mqc_summary

        def ch_db_done = KNOWLEDGE_WAREHOUSE.out.knowledge_db.collect()
        STOP_KNOWLEDGE_DB(db_host, db_port, ch_db_done)
    } else {
        db_host = channel.value(params.kw_db_host)
        db_port = channel.value(params.kw_db_port)
    }

    //
    // No downstream reporting or figures: the pipeline stops at the knowledge warehouse.
    //

    emit:
    unsupported    = PATHOGEN_ROUTER.out.unsupported
    phenotype_mutations = PATHOGEN_ROUTER.out.mutations
    phenotype_summary   = PATHOGEN_ROUTER.out.query_summary
    uniprotr_results    = PATHOGEN_ROUTER.out.uniprotr_results
    extractr_results    = PATHOGEN_ROUTER.out.extractr_results
    rbioapi_results     = PATHOGEN_ROUTER.out.rbioapi_results
    epi_raw             = PATHOGEN_ROUTER.out.epi_raw
    epi_search_summary  = PATHOGEN_ROUTER.out.epi_search_summary
    lit_results         = PATHOGEN_ROUTER.out.lit_results
    knowledge_db        = ch_knowledge_db
    knowledge_db_summary = ch_knowledge_db_summary
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
