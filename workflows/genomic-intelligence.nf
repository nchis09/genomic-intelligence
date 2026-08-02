/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { CLASSIFICATION         } from '../subworkflows/local/classification/main'
include { PATHOGEN_ROUTER        } from '../subworkflows/local/pathogen_router/main'
include { REPORTING              } from '../subworkflows/local/reporting/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_genomic-intelligence_pipeline'

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

    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()

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
    // pathogen workflow runs its own phenotype annotation (UniprotR +
    // UniProtExtractR + rbioapi) and tree annotation (PHYLO_ANNOTATE)
    // internally, since those depend on that pathogen's own auspice/results
    // (see subworkflows/local/ebola_workflow/main.nf). Groups whose pathogen
    // has no registered workflow are skipped with a warning; see
    // PATHOGEN_ROUTER.out.unsupported for the summary file.
    //
    // All Nextclade JSON outputs (one per sample x dataset run) and the
    // species_assignments.tsv (sample -> winning dataset) are needed to
    // resolve each tip's own Nextclade call for direct-vs-reference AA
    // mutations (see PHYLO_ANNOTATE). Both are broadcast to every species
    // group below (collect()/single-path outputs behave as value channels).
    ch_nextclade_json_all = CLASSIFICATION.out.json
        .map { _meta, json -> json }
        .collect()

    PATHOGEN_ROUTER(
        CLASSIFICATION.out.species_groups,
        ch_nextclade_json_all,
        CLASSIFICATION.out.assignments
    )

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name:  'genomic-intelligence_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // SUBWORKFLOW: Reporting (MultiQC)
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        PATHOGEN_ROUTER.out.epi_search_summary
            .filter { _meta, file -> file.name != 'NO_FILE' }
            .map { _meta, file -> file }
    )
    ch_multiqc_files = ch_multiqc_files.mix(
        PATHOGEN_ROUTER.out.epi_raw
            .filter { _meta, file -> file.name != 'NO_FILE' }
            .map { _meta, file -> file }
    )
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))

    REPORTING(
        ch_multiqc_files,
        multiqc_config,
        multiqc_logo
    )

    emit:
    multiqc_report = REPORTING.out.multiqc_report.map { _meta, report -> [report] }.toList()
    figures        = PATHOGEN_ROUTER.out.figures
    unsupported    = PATHOGEN_ROUTER.out.unsupported
    versions       = ch_versions
    phenotype_mutations = PATHOGEN_ROUTER.out.mutations
    phenotype_summary   = PATHOGEN_ROUTER.out.query_summary
    uniprotr_results    = PATHOGEN_ROUTER.out.uniprotr_results
    extractr_results    = PATHOGEN_ROUTER.out.extractr_results
    rbioapi_results     = PATHOGEN_ROUTER.out.rbioapi_results
    epi_raw             = PATHOGEN_ROUTER.out.epi_raw
    epi_search_summary  = PATHOGEN_ROUTER.out.epi_search_summary
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
