/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { CLASSIFICATION         } from '../subworkflows/local/classification/main'
include { PATHOGEN_ROUTER        } from '../subworkflows/local/pathogen_router/main'
include { KNOWLEDGE_WAREHOUSE    } from '../subworkflows/local/knowledge_warehouse/main'
include { PATHOGEN_IDENTIFICATION_WF } from '../subworkflows/local/pathogen_identification/main'
include { MUTATION_PROFILE_WF    } from '../subworkflows/local/mutation_profile/main'

include { REPORTING              } from '../subworkflows/local/reporting/main'

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
    ch_datasets = nextflow.Channel.fromList(params.nextclade_datasets)

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
    // Shared knowledge-warehouse: the subworkflow starts the shared PostgreSQL
    // instance internally, ingests per-species data, exports the whole
    // warehouse to a portable DuckDB file, then stops the server -- all
    // before this main workflow finishes, so PATHOGEN_IDENTIFICATION reads
    // that file instead of needing a live DB connection for the rest of the
    // run. A workflow.onComplete safety net (see main.nf) only acts if the
    // server is somehow still running (e.g. the pipeline failed midway).
    //
    def ch_knowledge_db = nextflow.Channel.empty()
    def ch_knowledge_db_summary = nextflow.Channel.empty()
    def ch_duckdb_dump = nextflow.Channel.empty()
    def ch_pathogen_id_tsv = nextflow.Channel.empty()
    def ch_pathogen_id_mqc = nextflow.Channel.empty()
    def ch_mutation_profile_tsv = nextflow.Channel.empty()
    def ch_mutation_profile_mqc = nextflow.Channel.empty()
    def ch_multiqc_report = nextflow.Channel.empty()

    db_host = nextflow.Channel.value(params.kw_db_host)
    db_port = nextflow.Channel.value(params.kw_db_port)

    if (!params.skip_knowledge_warehouse) {
        KNOWLEDGE_WAREHOUSE(
            PATHOGEN_ROUTER.out.kw_input.combine(db_host).combine(db_port)
        )
        ch_knowledge_db = KNOWLEDGE_WAREHOUSE.out.knowledge_db
        ch_knowledge_db_summary = KNOWLEDGE_WAREHOUSE.out.mqc_summary
        ch_duckdb_dump = KNOWLEDGE_WAREHOUSE.out.duckdb_dump

        if (!params.skip_pathogen_identification) {
            ch_knowledge_db
                .combine(ch_duckdb_dump)
                .map { meta, kw_dir, duckdb_file -> tuple(meta, duckdb_file) }
                .set { ch_pathogen_id }
            PATHOGEN_IDENTIFICATION_WF(ch_pathogen_id)
            ch_pathogen_id_tsv = PATHOGEN_IDENTIFICATION_WF.out.tsv
            ch_pathogen_id_mqc = PATHOGEN_IDENTIFICATION_WF.out.mqc_tsv
        }

        if (!params.skip_mutation_profile) {
            ch_knowledge_db
                .combine(ch_duckdb_dump)
                .map { meta, kw_dir, duckdb_file -> tuple(meta, duckdb_file) }
                .set { ch_mutation_profile }

            // Derive translations directory from Nextstrain/Nextclade outputs
            ch_translations = ch_knowledge_db
                .map { meta, kw_dir ->
                    def trans = file("${params.outdir}/nextstrain_ebola/${meta.species}/results/${meta.species}/translations")
                    trans.exists() ? trans : file("NO_TRANSLATIONS")
                }

            MUTATION_PROFILE_WF(ch_mutation_profile, ch_translations)
            ch_mutation_profile_tsv = MUTATION_PROFILE_WF.out.tsv
            ch_mutation_profile_mqc = MUTATION_PROFILE_WF.out.mqc_tsv
        }

    }

    //
    // Downstream reporting: MultiQC, fed with the per-species pathogen
    // identification tables and the knowledge-warehouse load summary.
    //
    if (!params.skip_multiqc) {
        def ch_multiqc_files = ch_pathogen_id_mqc.mix(ch_knowledge_db_summary.map { _meta, mqc_path -> mqc_path })

        REPORTING(
            ch_multiqc_files,
            multiqc_config,
            multiqc_logo
        )
        ch_multiqc_report = REPORTING.out.multiqc_report
    }

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
    lit_evidence        = PATHOGEN_ROUTER.out.lit_evidence
    knowledge_db        = ch_knowledge_db
    knowledge_db_summary = ch_knowledge_db_summary
    identification_tsv   = ch_pathogen_id_tsv
    multiqc_report       = ch_multiqc_report
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
