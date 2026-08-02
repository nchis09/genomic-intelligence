/*
 * Subworkflow: EBOLA_WORKFLOW
 *
 * Pathogen-specific workflow for Ebola (orthoebolavirus) species groups:
 *   1. BIOINFORMATICS_ANALYSIS — Nextstrain Ebola build + internal stage markers
 *   2. PHENOTYPE_ANNOTATION    — UniprotR + UniProtExtractR + rbioapi annotation
 *   3. EPIDEMIOLOGICAL_DATA    — placeholder for future epidemiology steps
 *   4. EVIDENCE_SYNTHESIS      — PHYLO_ANNOTATE + PHYLO_VISUALIZE
 *
 * Input:  ch_species_data       - channel of [ meta, fasta, metadata ]
 *         ch_nextclade_json_all  - path: all Nextclade JSONs (broadcast)
 *         ch_species_assignments - path: species_assignments.tsv (broadcast)
 * Output: figures, mutations, query_summary, rbioapi_results,
 *         uniprotr_results, extractr_results
 */

include { BIOINFORMATICS_ANALYSIS } from '../bioinformatics_analysis/main'
include { PHENOTYPE_ANNOTATION    } from '../phenotype_annotation/main'
include { EPIDEMIOLOGICAL_DATA    } from '../epidemiological_data/main'
include { EVIDENCE_SYNTHESIS      } from '../evidence_synthesis/main'
include { KNOWLEDGE_WAREHOUSE     } from '../knowledge_warehouse/main'

workflow EBOLA_WORKFLOW {
    take:
    ch_species_data        // channel: [ val(meta), path(fasta), path(metadata) ]
    ch_nextclade_json_all  // path: all Nextclade JSONs (broadcast/value channel)
    ch_species_assignments // path: species_assignments.tsv (broadcast/value channel)

    main:
    //
    // SUBWORKFLOW: Bioinformatics analysis (Nextstrain Ebola build)
    //
    BIOINFORMATICS_ANALYSIS(ch_species_data)

    ch_auspice_results = BIOINFORMATICS_ANALYSIS.out.auspice
        .join(BIOINFORMATICS_ANALYSIS.out.results, by: [0])
        .map { meta, auspice_files, results_dir ->
            def auspice = auspice_files instanceof List ? auspice_files[0] : auspice_files
            [ meta, auspice, results_dir ]
        }

    //
    // SUBWORKFLOW: Phenotype annotation (UniprotR + UniProtExtractR + rbioapi)
    //
    if (!params.skip_phenotype_annotation) {
        PHENOTYPE_ANNOTATION(ch_auspice_results)
    }

    // When phenotype annotation is skipped, fall back to NO_FILE placeholders
    // so downstream steps still get channels to join on.
    ch_rbioapi_for_synthesis = params.skip_phenotype_annotation
        ? BIOINFORMATICS_ANALYSIS.out.auspice.map { meta, _auspice -> [ meta, file('NO_FILE') ] }
        : PHENOTYPE_ANNOTATION.out.rbioapi_results

    ch_phenotype_mutations = params.skip_phenotype_annotation
        ? channel.empty()
        : PHENOTYPE_ANNOTATION.out.mutations

    ch_phenotype_summary = params.skip_phenotype_annotation
        ? channel.empty()
        : PHENOTYPE_ANNOTATION.out.query_summary

    ch_uniprotr_results = params.skip_phenotype_annotation
        ? channel.empty()
        : PHENOTYPE_ANNOTATION.out.uniprotr_results

    ch_extractr_results = params.skip_phenotype_annotation
        ? channel.empty()
        : PHENOTYPE_ANNOTATION.out.extractr_results

    ch_rbioapi_results = params.skip_phenotype_annotation
        ? channel.empty()
        : PHENOTYPE_ANNOTATION.out.rbioapi_results

    //
    // SUBWORKFLOW: Epidemiological data (HDX via rhdx)
    //
    EPIDEMIOLOGICAL_DATA(ch_species_data)
    ch_epi_raw        = EPIDEMIOLOGICAL_DATA.out.epi_raw
    ch_epi_search_summary = EPIDEMIOLOGICAL_DATA.out.search_summary

    //
    // SUBWORKFLOW: Knowledge warehouse (PostgreSQL)
    //
    if (!params.skip_knowledge_warehouse) {
        KNOWLEDGE_WAREHOUSE(
            ch_species_data,
            ch_species_assignments,
            ch_epi_raw,
            ch_epi_search_summary
        )
        ch_knowledge_db = KNOWLEDGE_WAREHOUSE.out.knowledge_db
    } else {
        ch_knowledge_db = channel.empty()
    }

    //
    // SUBWORKFLOW: Evidence synthesis (tree annotation + visualization)
    //
    EVIDENCE_SYNTHESIS(
        ch_auspice_results,
        ch_rbioapi_for_synthesis,
        ch_nextclade_json_all,
        ch_species_assignments
    )

    emit:
    figures          = EVIDENCE_SYNTHESIS.out.figures          // channel: [ meta, png ]
    mutations        = ch_phenotype_mutations                  // channel: [ meta, tsv ]
    query_summary    = ch_phenotype_summary                    // channel: [ meta, json ]
    uniprotr_results = ch_uniprotr_results                     // channel: [ meta, dir ]
    extractr_results = ch_extractr_results                     // channel: [ meta, dir ]
    rbioapi_results  = ch_rbioapi_results                       // channel: [ meta, dir ]
    epi_raw          = ch_epi_raw                              // channel: [ meta, epi_data.csv ]
    epi_search_summary = ch_epi_search_summary                 // channel: [ meta, rhdx_search_results.tsv ]
    knowledge_db     = ch_knowledge_db                         // channel: [ meta, knowledge_warehouse ]
}
