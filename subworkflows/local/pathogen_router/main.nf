/*
 * Subworkflow: PATHOGEN_ROUTER
 *
 * Routes each species group (output of CLASSIFICATION) to the correct
 * pathogen-specific workflow based on `meta.pathogen`.
 *
 * WORKFLOW_REGISTRY below is the single place to register which pathogen
 * families are handled by which workflow. Multiple pathogen families can
 * share the same workflow entry (e.g. future flu subtypes could all map to
 * a shared 'FLU' workflow). Modules used by a workflow (e.g. NEXTSTRAIN_EBOLA,
 * QUERY_KNOWLEDGE_DB, PHYLO_VISUALIZE) stay generic/reusable so a future
 * workflow can call them again instead of duplicating logic.
 *
 * To add a new pathogen workflow:
 *   1. Add its pathogen family key(s) to WORKFLOW_REGISTRY below.
 *   2. Add a new `.branch{}` case routing that family to your new subworkflow.
 *   3. Call your subworkflow and mix its outputs into the emits below.
 *
 * Any pathogen family NOT present in WORKFLOW_REGISTRY is skipped: an
 * immediate warning is logged, and the group is recorded in
 * `unsupported_pathogens.tsv` for an end-of-run summary (see main.nf
 * workflow.onComplete).
 *
 * Input:  ch_species_data       - channel of [ meta(pathogen, species), fasta, metadata ]
 *         ch_nextclade_json_all  - path: all Nextclade JSONs (broadcast)
 *         ch_species_assignments - path: species_assignments.tsv (broadcast)
 * Output: kw_input, mutations, query_summary, rbioapi_results,
 *         uniprotr_results, extractr_results, epi_raw, epi_search_summary,
 *         lit_results, unsupported
 */

include { ROUTE_PATHOGEN } from '../../../modules/local/route_pathogen/main'
include { EBOLA_WORKFLOW } from '../ebola_workflow/main'

workflow PATHOGEN_ROUTER {
    take:
    ch_species_data        // channel: [ val(meta), path(fasta), path(metadata) ]
    ch_nextclade_json_all  // path: all Nextclade JSONs (broadcast/value channel)
    ch_species_assignments // path: species_assignments.tsv (broadcast/value channel)

    main:
    //
    // Registry: workflow name -> list of pathogen families it handles.
    // This is the single place to register which pathogen families are
    // handled by which workflow. Multiple families can share one workflow
    // (e.g. future flu subtypes could all map to a shared 'FLU' entry).
    //
    def WORKFLOW_REGISTRY = [
        EBOLA: ['orthoebolavirus'],
        // Future pathogens, e.g.:
        // FLU: ['flu'],
        // RSV: ['rsv'],
    ]

    //
    // MODULE: Router marker — gives the pathogen-dispatch decision below a
    // real node in the Nextflow DAG. `.branch{}` alone is a pure channel
    // operator and leaves no trace in the execution graph, so without this
    // the fork point is invisible to DAG-visualization tools (e.g. nf-metro).
    //
    ROUTE_PATHOGEN(ch_species_data)

    //
    // Branch species groups by pathogen family into known workflows vs unsupported
    //
    ch_branched = ROUTE_PATHOGEN.out.routed.branch { meta, fasta, metadata ->
        ebola       : WORKFLOW_REGISTRY.EBOLA.contains(meta.pathogen)
        unsupported : true
    }

    //
    // Log an immediate warning for each unsupported pathogen group, and
    // record it for the end-of-run summary written in main.nf's onComplete.
    //
    ch_branched.unsupported.subscribe { meta, fasta, metadata ->
        log.warn "Workflow for pathogen '${meta.pathogen}' does not exist yet — skipping sample group '${meta.id}' (species: ${meta.species})."
    }

    ch_unsupported = ch_branched.unsupported
        .map { meta, fasta, metadata -> "${meta.id}\t${meta.pathogen}\t${meta.species}" }
        .collectFile(
            name: 'unsupported_pathogens.tsv',
            newLine: true,
            seed: 'sample_id\tpathogen\tspecies',
            storeDir: "${params.outdir}/pipeline_info"
        )

    //
    // Run the Ebola workflow for orthoebolavirus species groups. Ebola now
    // owns only the per-species analysis modules (bioinformatics, phenotype
    // annotation, epidemiological data). The shared DB and figures are handled
    // in the main workflow so they can be used across pathogen-specific
    // analyses.
    //
    EBOLA_WORKFLOW(ch_branched.ebola, ch_nextclade_json_all, ch_species_assignments)

    //
    // As more pathogen workflows are registered above, mix their outputs
    // into these emits (e.g. `.mix(FLU_WORKFLOW.out.auspice)`) — only one
    // pathogen workflow fires per species group, so `.mix()` (not `.join()`)
    // is the correct way to recombine them.
    //
    emit:
    kw_input         = EBOLA_WORKFLOW.out.kw_input         // channel: [ meta, kw inputs ]
    mutations        = EBOLA_WORKFLOW.out.mutations        // channel: [ meta, tsv ]
    query_summary    = EBOLA_WORKFLOW.out.query_summary    // channel: [ meta, json ]
    uniprotr_results = EBOLA_WORKFLOW.out.uniprotr_results // channel: [ meta, dir ]
    extractr_results = EBOLA_WORKFLOW.out.extractr_results // channel: [ meta, dir ]
    rbioapi_results  = EBOLA_WORKFLOW.out.rbioapi_results  // channel: [ meta, dir ]
    epi_raw          = EBOLA_WORKFLOW.out.epi_raw          // channel: [ meta, epi_data.csv ]
    epi_search_summary = EBOLA_WORKFLOW.out.epi_search_summary // channel: [ meta, rhdx_search_results.tsv ]
    lit_results      = EBOLA_WORKFLOW.out.lit_results      // channel: [ meta, [ literature result files ] ]
    lit_evidence     = EBOLA_WORKFLOW.out.lit_evidence     // channel: [ meta, [ clean/*.json files ] ]
    unsupported      = ch_unsupported                        // path: unsupported_pathogens.tsv
}
