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

    // Per-meta, NO_FILE-safe fallbacks for the knowledge warehouse (needs a
    // real entry per species group to join on, unlike the emit-only channels
    // above which can stay empty when phenotype annotation is skipped).
    ch_uniprotr_for_kw = params.skip_phenotype_annotation
        ? BIOINFORMATICS_ANALYSIS.out.auspice.map { meta, _auspice -> [ meta, file('NO_FILE_uniprotr_kw') ] }
        : PHENOTYPE_ANNOTATION.out.uniprotr_results

    ch_extractr_for_kw = params.skip_phenotype_annotation
        ? BIOINFORMATICS_ANALYSIS.out.auspice.map { meta, _auspice -> [ meta, file('NO_FILE_extractr_kw') ] }
        : PHENOTYPE_ANNOTATION.out.extractr_results

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
        // Derive a local copy of the broadcast species_assignments channel so
        // the DAG traces it back to a node within this workflow's own region
        // rather than skipping straight back to SPECIES_ASSIGN in Classification
        // (purely cosmetic: keeps the metro-map connecting line intact).
        ch_species_assignments_kw = ch_species_assignments.map { it }

        // Bundle every knowledge-warehouse input into a single per-meta tuple
        // channel before the call, so KNOWLEDGE_WAREHOUSE's DAG entry node has
        // exactly one inbound edge (matching how single-argument subworkflow
        // calls like EPIDEMIOLOGICAL_DATA render in the metro map).
        def no_file_meta     = file('NO_FILE_metadata')
        def no_file_epi      = file('NO_FILE_epi')
        def no_file_summary  = file('NO_FILE_epi_summary')
        def no_file_bioinfo  = file('NO_FILE_bioinfo')
        def no_file_uniprotr = file('NO_FILE_uniprotr')
        def no_file_extractr = file('NO_FILE_extractr')
        def no_file_rbioapi  = file('NO_FILE_rbioapi')

        // Epidemiological data is fetched once per pathogen, so it needs to be
        // re-keyed by pathogen and broadcast to every species group sharing it.
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

        // Bioinformatics/phenotype outputs are already per-species (per-meta).
        ch_bioinfo_safe  = BIOINFORMATICS_ANALYSIS.out.results
            .map { meta, dir -> [ meta, (dir && !dir.name.startsWith('NO_FILE')) ? dir : no_file_bioinfo ] }
        ch_uniprotr_safe = ch_uniprotr_for_kw
            .map { meta, dir -> [ meta, (dir && !dir.name.startsWith('NO_FILE')) ? dir : no_file_uniprotr ] }
        ch_extractr_safe = ch_extractr_for_kw
            .map { meta, dir -> [ meta, (dir && !dir.name.startsWith('NO_FILE')) ? dir : no_file_extractr ] }
        ch_rbioapi_safe  = ch_rbioapi_for_synthesis
            .map { meta, dir -> [ meta, (dir && !dir.name.startsWith('NO_FILE')) ? dir : no_file_rbioapi ] }

        ch_kw_input = ch_species_by_pathogen
            .join(ch_epi_by_pathogen, by: 0, remainder: true)
            .join(ch_epi_summary_by_pathogen, by: 0, remainder: true)
            .map { pathogen, meta, metadata, epi_dir, epi_summary ->
                [ meta, metadata ?: no_file_meta, epi_dir ?: no_file_epi, epi_summary ?: no_file_summary ]
            }
            .join(ch_bioinfo_safe, by: 0)
            .join(ch_uniprotr_safe, by: 0)
            .join(ch_extractr_safe, by: 0)
            .join(ch_rbioapi_safe, by: 0)
            .combine(ch_species_assignments_kw)
            .map { meta, metadata, epi_dir, epi_summary, bioinfo_dir, uniprotr_dir, extractr_dir, rbioapi_dir, assignments ->
                [ meta, assignments, metadata, epi_dir, epi_summary, bioinfo_dir, uniprotr_dir, extractr_dir, rbioapi_dir ]
            }

        KNOWLEDGE_WAREHOUSE(ch_kw_input)
        ch_knowledge_db = KNOWLEDGE_WAREHOUSE.out.knowledge_db
        ch_knowledge_db_summary = KNOWLEDGE_WAREHOUSE.out.mqc_summary
    } else {
        ch_knowledge_db = channel.empty()
        ch_knowledge_db_summary = ch_species_data.map { meta, _fasta, _metadata -> [ meta, file('NO_FILE') ] }
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
    knowledge_db_summary = ch_knowledge_db_summary             // channel: [ meta, mqc_summary.tsv ]
}
