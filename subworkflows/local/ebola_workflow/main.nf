/*
 * Subworkflow: EBOLA_WORKFLOW
 *
 * Pathogen-specific workflow for Ebola (orthoebolavirus) species groups:
 *   1. BIOINFORMATICS_AND_EPIDEMIOLOGICAL — Nextstrain build + HDX epidemiology
 *   2. PHENOTYPE_ANNOTATION               — UniprotR + UniProtExtractR + rbioapi
 *
 * Input:  ch_species_data       - channel of [ meta, fasta, metadata ]
 *         ch_nextclade_json_all  - path: all Nextclade JSONs (broadcast)
 *         ch_species_assignments - path: species_assignments.tsv (broadcast/value channel)
 * Output: kw_input, mutations, query_summary, rbioapi_results,
 *         uniprotr_results, extractr_results, epi_raw, epi_search_summary
 */

include { BIOINFORMATICS_AND_EPIDEMIOLOGICAL } from '../bioinformatics_and_epi/main'
include { PHENOTYPE_ANNOTATION               } from '../phenotype_annotation/main'

workflow EBOLA_WORKFLOW {
    take:
    ch_species_data        // channel: [ val(meta), path(fasta), path(metadata) ]
    ch_nextclade_json_all  // path: all Nextclade JSONs (broadcast/value channel)
    ch_species_assignments // path: species_assignments.tsv (broadcast/value channel)

    main:
    //
    // SUBWORKFLOW: Bioinformatics + Epidemiological data
    //
    BIOINFORMATICS_AND_EPIDEMIOLOGICAL(ch_species_data)
    ch_auspice_results = BIOINFORMATICS_AND_EPIDEMIOLOGICAL.out.auspice_results
    ch_epi_raw        = BIOINFORMATICS_AND_EPIDEMIOLOGICAL.out.epi_raw
    ch_epi_search_summary = BIOINFORMATICS_AND_EPIDEMIOLOGICAL.out.epi_search_summary

    //
    // SUBWORKFLOW: Phenotype annotation (UniprotR + UniProtExtractR + rbioapi)
    //
    if (!params.skip_phenotype_annotation) {
        PHENOTYPE_ANNOTATION(ch_auspice_results)
    }

    // When phenotype annotation is skipped, fall back to NO_FILE placeholders
    // so downstream steps still get channels to join on.
    ch_rbioapi_for_synthesis = params.skip_phenotype_annotation
        ? BIOINFORMATICS_AND_EPIDEMIOLOGICAL.out.auspice.map { meta, _auspice -> [ meta, file('NO_FILE') ] }
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
        ? BIOINFORMATICS_AND_EPIDEMIOLOGICAL.out.auspice.map { meta, _auspice -> [ meta, file('NO_FILE_uniprotr_kw') ] }
        : PHENOTYPE_ANNOTATION.out.uniprotr_results

    ch_extractr_for_kw = params.skip_phenotype_annotation
        ? BIOINFORMATICS_AND_EPIDEMIOLOGICAL.out.auspice.map { meta, _auspice -> [ meta, file('NO_FILE_extractr_kw') ] }
        : PHENOTYPE_ANNOTATION.out.extractr_results

    // Bundle every EXTRACT_QUERY_PROTEINS output (discovery.tsv, accessions.txt,
    // uniprot_download.tsv, query_mutations.tsv, query_proteins.fasta,
    // discovery_summary.json) into one per-meta list for the knowledge
    // warehouse's provenance registration -- no images.
    ch_query_data_for_kw = params.skip_phenotype_annotation
        ? BIOINFORMATICS_AND_EPIDEMIOLOGICAL.out.auspice.map { meta, _auspice -> [ meta, [ file('NO_FILE_query_data') ] ] }
        : PHENOTYPE_ANNOTATION.out.discovery
            .join(PHENOTYPE_ANNOTATION.out.accessions, by: 0)
            .join(PHENOTYPE_ANNOTATION.out.uniprot_tsv, by: 0)
            .join(PHENOTYPE_ANNOTATION.out.mutations, by: 0)
            .join(PHENOTYPE_ANNOTATION.out.query_proteins, by: 0)
            .join(PHENOTYPE_ANNOTATION.out.query_summary, by: 0)
            .map { meta, discovery, accessions, uniprot_tsv, mutations, proteins, summary ->
                [ meta, [ discovery, accessions, uniprot_tsv, mutations, proteins, summary ] ]
            }

    //
    // Bundle all knowledge-warehouse inputs into a single per-meta tuple.
    // The actual KNOWLEDGE_WAREHOUSE subworkflow is called in the main
    // workflow so the shared DB is not coupled to one pathogen's subworkflow.
    //
    if (!params.skip_knowledge_warehouse) {
        // Derive a local copy of the broadcast species_assignments channel so
        // the DAG traces it back to a node within this workflow's own region
        // rather than skipping straight back to SPECIES_ASSIGN in Classification
        // (purely cosmetic: keeps the metro-map connecting line intact).
        ch_species_assignments_kw = ch_species_assignments.map { it }

        def no_file_meta     = file('NO_FILE_metadata')
        def no_file_epi      = file('NO_FILE_epi')
        def no_file_summary  = file('NO_FILE_epi_summary')
        def no_file_bioinfo  = file('NO_FILE_bioinfo')
        def no_file_uniprotr = file('NO_FILE_uniprotr')
        def no_file_extractr = file('NO_FILE_extractr')
        def no_file_rbioapi  = file('NO_FILE_rbioapi')
        def no_file_auspice  = file('NO_FILE_auspice')
        def no_file_tree     = file('NO_FILE_tree')
        def no_file_hmm      = file('NO_FILE_hmm')

        // Epidemiological data is now fetched per species, so it is keyed by meta.
        ch_species_for_kw = ch_species_data
            .map { meta, fasta, metadata ->
                def meta_tsv = (metadata && !metadata.name.startsWith('NO_FILE')) ? metadata : no_file_meta
                [ meta, meta_tsv ]
            }

        ch_epi_for_kw = ch_epi_raw
            .map { meta, epi_dir ->
                def epi = (epi_dir && !epi_dir.name.startsWith('NO_FILE')) ? epi_dir : no_file_epi
                [ meta, epi ]
            }

        ch_epi_summary_for_kw = ch_epi_search_summary
            .map { meta, tsv ->
                def summary = (tsv && !tsv.name.startsWith('NO_FILE')) ? tsv : no_file_summary
                [ meta, summary ]
            }

        // Bioinformatics/phenotype outputs are already per-species (per-meta).
        ch_bioinfo_safe  = BIOINFORMATICS_AND_EPIDEMIOLOGICAL.out.results
            .map { meta, dir -> [ meta, (dir && !dir.name.startsWith('NO_FILE')) ? dir : no_file_bioinfo ] }
        ch_uniprotr_safe = ch_uniprotr_for_kw
            .map { meta, dir -> [ meta, (dir && !dir.name.startsWith('NO_FILE')) ? dir : no_file_uniprotr ] }
        ch_extractr_safe = ch_extractr_for_kw
            .map { meta, dir -> [ meta, (dir && !dir.name.startsWith('NO_FILE')) ? dir : no_file_extractr ] }
        ch_rbioapi_safe  = ch_rbioapi_for_synthesis
            .map { meta, dir -> [ meta, (dir && !dir.name.startsWith('NO_FILE')) ? dir : no_file_rbioapi ] }

        ch_auspice_for_kw = ch_auspice_results
            .join(BIOINFORMATICS_AND_EPIDEMIOLOGICAL.out.tree, by: 0, remainder: true)
            .map { meta, auspice, _results, tree ->
                [ meta, (auspice && !auspice.name.startsWith('NO_FILE')) ? auspice : no_file_auspice, (tree && !tree.name.startsWith('NO_FILE')) ? tree : no_file_tree ]
            }

        ch_hmm_for_kw = params.skip_phenotype_annotation
            ? BIOINFORMATICS_AND_EPIDEMIOLOGICAL.out.auspice.map { meta, _auspice -> [ meta, no_file_hmm ] }
            : PHENOTYPE_ANNOTATION.out.hmm_results

        ch_kw_input = ch_species_for_kw
            .join(ch_epi_for_kw, by: 0, remainder: true)
            .join(ch_epi_summary_for_kw, by: 0, remainder: true)
            .map { meta, metadata, epi_dir, epi_summary ->
                [ meta, metadata ?: no_file_meta, epi_dir ?: no_file_epi, epi_summary ?: no_file_summary ]
            }
            .join(ch_bioinfo_safe, by: 0)
            .join(ch_uniprotr_safe, by: 0)
            .join(ch_extractr_safe, by: 0)
            .join(ch_rbioapi_safe, by: 0)
            .join(ch_auspice_for_kw, by: 0)
            .join(ch_query_data_for_kw, by: 0)
            .join(ch_hmm_for_kw, by: 0)
            .combine(ch_species_assignments_kw)
            .map { meta, metadata, epi_dir, epi_summary, bioinfo_dir, uniprotr_dir, extractr_dir, rbioapi_dir, auspice, tree, query_data, hmm, assignments ->
                [ meta, assignments, metadata, epi_dir, epi_summary, bioinfo_dir, uniprotr_dir, extractr_dir, rbioapi_dir, auspice, tree, query_data, hmm ]
            }
    } else {
        ch_kw_input = channel.empty()
    }

    emit:
    kw_input         = ch_kw_input                             // channel: [ meta, assignments, metadata, epi_dir, epi_summary, bioinfo_dir, uniprotr_dir, extractr_dir, rbioapi_dir, auspice, tree, query_data, hmm ]
    mutations        = ch_phenotype_mutations                  // channel: [ meta, tsv ]
    query_summary    = ch_phenotype_summary                    // channel: [ meta, json ]
    uniprotr_results = ch_uniprotr_results                     // channel: [ meta, dir ]
    extractr_results = ch_extractr_results                     // channel: [ meta, dir ]
    rbioapi_results  = ch_rbioapi_results                       // channel: [ meta, dir ]
    epi_raw          = ch_epi_raw                              // channel: [ meta, epi_data.csv ]
    epi_search_summary = ch_epi_search_summary                 // channel: [ meta, rhdx_search_results.tsv ]
}
