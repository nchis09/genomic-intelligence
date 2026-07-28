/*
 * Subworkflow: PHENOTYPE_ANNOTATION
 *
 * Extract query-sample mutations and protein info from Nextstrain outputs,
 * then annotate using UniprotR, UniProtExtractR, and rbioapi (cloned repos).
 *
 * Input:  ch_auspice_results - channel of [ meta, auspice_json, results_dir ]
 * Output: mutations, uniprotr_results, extractr_results, rbioapi_results, query_summary
 */

include { EXTRACT_QUERY_PROTEINS } from '../../../modules/local/extract_query_proteins/main'
include { UNIPROT_ANNOTATE       } from '../../../modules/local/uniprot_annotate/main'

workflow PHENOTYPE_ANNOTATION {
    take:
    ch_auspice_results  // channel: [ val(meta), path(auspice_json), path(results_dir) ]

    main:
    //
    // MODULE: Extract query sample mutations + download UniProtKB TSV
    //
    EXTRACT_QUERY_PROTEINS(ch_auspice_results)

    //
    // MODULE: Run UniprotR + UniProtExtractR + rbioapi on extracted data
    //
    ch_annotate_input = EXTRACT_QUERY_PROTEINS.out.accessions
        .join(EXTRACT_QUERY_PROTEINS.out.uniprot_tsv, by: [0])
        .join(EXTRACT_QUERY_PROTEINS.out.mutations, by: [0])

    UNIPROT_ANNOTATE(ch_annotate_input)

    emit:
    discovery        = EXTRACT_QUERY_PROTEINS.out.discovery        // [ meta, tsv ]
    mutations        = EXTRACT_QUERY_PROTEINS.out.mutations        // [ meta, tsv ]
    query_summary    = EXTRACT_QUERY_PROTEINS.out.summary          // [ meta, json ]
    query_proteins   = EXTRACT_QUERY_PROTEINS.out.proteins         // [ meta, fasta ]
    uniprotr_results = UNIPROT_ANNOTATE.out.uniprotr_results      // [ meta, dir ]
    extractr_results = UNIPROT_ANNOTATE.out.extractr_results      // [ meta, dir ]
    rbioapi_results  = UNIPROT_ANNOTATE.out.rbioapi_results        // [ meta, dir ]
}
