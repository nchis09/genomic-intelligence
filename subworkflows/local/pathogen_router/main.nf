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
 * PHYLO_ANNOTATE, PHYLO_VISUALIZE) stay generic/reusable so a future workflow
 * can call them again instead of duplicating logic.
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
 * Input:  ch_species_data - channel of [ meta(pathogen, species), fasta, metadata ]
 * Output: figures, tip_metadata, auspice, results, unsupported
 */

include { EBOLA_WORKFLOW } from '../ebola_workflow/main'

workflow PATHOGEN_ROUTER {
    take:
    ch_species_data  // channel: [ val(meta), path(fasta), path(metadata) ]

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
    // Branch species groups by pathogen family into known workflows vs unsupported
    //
    ch_branched = ch_species_data.branch { meta, fasta, metadata ->
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
    // Run the Ebola workflow for orthoebolavirus species groups
    //
    EBOLA_WORKFLOW(ch_branched.ebola)

    emit:
    figures      = EBOLA_WORKFLOW.out.figures       // channel: [ meta, png ]
    tip_metadata = EBOLA_WORKFLOW.out.tip_metadata   // channel: [ meta, tsv ]
    auspice      = EBOLA_WORKFLOW.out.auspice        // channel: [ meta, json ]
    results      = EBOLA_WORKFLOW.out.results        // channel: [ meta, dir ]
    unsupported  = ch_unsupported                    // path: unsupported_pathogens.tsv
}
