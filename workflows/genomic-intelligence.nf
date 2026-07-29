/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { CLASSIFICATION         } from '../subworkflows/local/classification/main'
include { PATHOGEN_ROUTER        } from '../subworkflows/local/pathogen_router/main'
include { PHENOTYPE_ANNOTATION   } from '../subworkflows/local/phenotype_annotation/main'
include { PHYLO_ANNOTATE         } from '../modules/local/phylo_annotate/main'
include { PHYLO_VISUALIZE        } from '../modules/local/phylo_visualize/main'
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
    // pathogen-specific workflow (currently only Ebola is registered).
    // Groups whose pathogen has no registered workflow are skipped with a
    // warning; see PATHOGEN_ROUTER.out.unsupported for the summary file.
    PATHOGEN_ROUTER(CLASSIFICATION.out.species_groups)

    //
    // SUBWORKFLOW: Phenotype annotation (UniprotR + UniProtExtractR + rbioapi)
    //
    // Uses Auspice JSON + results from PATHOGEN_ROUTER to extract query sample
    // mutations and retrieve species-specific UniProt functional annotations,
    // mutation-level phenotype effects, STRING interactions, and Reactome pathways.
    //
    if (!params.skip_phenotype_annotation) {
        ch_phenotype_input = PATHOGEN_ROUTER.out.auspice
            .join(PATHOGEN_ROUTER.out.results, by: [0])
            .map { meta, auspice_files, results_dir ->
                def auspice = auspice_files instanceof List ? auspice_files[0] : auspice_files
                [ meta, auspice, results_dir ]
            }

        PHENOTYPE_ANNOTATION(ch_phenotype_input)
    }

    //
    // MODULES: Tree visualization (PHYLO_ANNOTATE + PHYLO_VISUALIZE)
    //
    // Deliberately run AFTER phenotype annotation (not inside the
    // pathogen-specific workflow) — the tree figure incorporates
    // phenotype-derived annotations (e.g. mutations with documented
    // phenotypic evidence from rbioapi_results), so joining on
    // PHENOTYPE_ANNOTATION.out.rbioapi_results enforces that ordering.
    // When phenotype annotation is skipped, fall back to a NO_FILE
    // placeholder per group so the tree still gets built without waiting
    // on a subworkflow that was never run.
    //
    ch_rbioapi_for_tree = params.skip_phenotype_annotation
        ? PATHOGEN_ROUTER.out.auspice.map { meta, _auspice -> [ meta, file('NO_FILE') ] }
        : PHENOTYPE_ANNOTATION.out.rbioapi_results

    // All Nextclade JSON outputs (one per sample x dataset run) and the
    // species_assignments.tsv (sample -> winning dataset) are needed to
    // resolve each tip's own Nextclade call for direct-vs-reference AA
    // mutations (see PHYLO_ANNOTATE). Both are broadcast to every species
    // group below (collect()/single-path outputs behave as value channels).
    ch_nextclade_json_all = CLASSIFICATION.out.json
        .map { _meta, json -> json }
        .collect()

    ch_annotate_input = PATHOGEN_ROUTER.out.auspice
        .join(PATHOGEN_ROUTER.out.results, by: [0])
        .join(ch_rbioapi_for_tree, by: [0])
        .map { meta, auspice_files, results_dir, rbioapi_dir ->
            def auspice = auspice_files instanceof List ? auspice_files[0] : auspice_files
            [ meta, auspice, results_dir, rbioapi_dir ]
        }

    PHYLO_ANNOTATE(
        ch_annotate_input,
        ch_nextclade_json_all,
        CLASSIFICATION.out.assignments
    )

    ch_plot_base = PHYLO_ANNOTATE.out.tip_metadata
        .join(PATHOGEN_ROUTER.out.results, by: [0])
        .map { meta, tip_meta, results_dir ->
            [ meta, results_dir, tip_meta ]
        }

    ch_plot_input = ch_plot_base
        .join(PHYLO_ANNOTATE.out.mutation_matrix, by: [0], remainder: true)
        .join(PHYLO_ANNOTATE.out.mutation_legend, by: [0], remainder: true)
        .join(PHYLO_ANNOTATE.out.protein_burden, by: [0], remainder: true)
        .map { meta, results, tip_meta, mut_matrix=null, mut_legend=null, protein_burden=null ->
            [ meta, results, tip_meta, mut_matrix ?: file('NO_FILE'), mut_legend ?: file('NO_FILE'), protein_burden ?: file('NO_FILE') ]
        }

    PHYLO_VISUALIZE(ch_plot_input)

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
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'genomic-intelligence'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )

    emit:
    multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList()
    figures        = PHYLO_VISUALIZE.out.tree_heatmap.mix(PHYLO_VISUALIZE.out.geo_map)
    auspice        = PATHOGEN_ROUTER.out.auspice
    unsupported    = PATHOGEN_ROUTER.out.unsupported
    versions       = ch_versions
    phenotype_mutations = params.skip_phenotype_annotation ? channel.empty() : PHENOTYPE_ANNOTATION.out.mutations
    phenotype_summary   = params.skip_phenotype_annotation ? channel.empty() : PHENOTYPE_ANNOTATION.out.query_summary
    rbioapi_results     = params.skip_phenotype_annotation ? channel.empty() : PHENOTYPE_ANNOTATION.out.rbioapi_results
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
