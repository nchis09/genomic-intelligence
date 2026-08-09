/*
 * Subworkflow: PHENOTYPE_ANNOTATION
 *
 * Extract query-sample mutations and protein info from Nextstrain outputs,
 * then annotate using UniprotR, UniProtExtractR, and rbioapi (cloned repos).
 *
 * Input:  ch_auspice_results - channel of [ meta, auspice_json, results_dir ]
 * Output: mutations, uniprotr_results, extractr_results, rbioapi_results, query_summary
 */

include { EXTRACT_QUERY_PROTEINS      } from '../../../modules/local/extract_query_proteins/main'
include { UNIPROT_ANNOTATE            } from '../../../modules/local/uniprot_annotate/main'
include { UNIPROT_EXTRACTR_ANNOTATE   } from '../../../modules/local/uniprotextractr_annotate/main'
include { RBIOAPI_ANNOTATE            } from '../../../modules/local/rbioapi_annotate/main'
include { PREPARE_HMMDB               } from '../../../modules/local/prepare_hmmdb/main'
include { HMM_ANNOTATE                } from '../../../modules/local/hmm_annotate/main'

workflow PHENOTYPE_ANNOTATION {
    take:
    ch_auspice_results  // channel: [ val(meta), path(auspice_json), path(results_dir) ]

    main:
    //
    // MODULE: Extract query sample mutations + download UniProtKB TSV
    //
    EXTRACT_QUERY_PROTEINS(ch_auspice_results)

    //
    // MODULES: HMM database annotation (optional)
    //
    ch_hmm_annotations = channel.empty()
    ch_hmm_results     = ch_auspice_results.map { meta, _auspice, _results -> [ meta, file('NO_FILE_hmm') ] }
    if (!params.skip_hmm_annotation) {
        PREPARE_HMMDB(params.hmm_db_url)
        HMM_ANNOTATE(EXTRACT_QUERY_PROTEINS.out.proteins.combine(PREPARE_HMMDB.out.hmm_db_dir))
        ch_hmm_annotations = HMM_ANNOTATE.out.hmm_annotations
        ch_hmm_results = HMM_ANNOTATE.out.hmm_annotations
            .join(HMM_ANNOTATE.out.hmm_sequence_table, by: 0)
            .join(HMM_ANNOTATE.out.hmm_domain_table, by: 0)
            .join(HMM_ANNOTATE.out.hmm_pfam_table, by: 0)
            .join(HMM_ANNOTATE.out.hmm_report, by: 0)
            .map { meta, ann, seq, dom, pfam, rep ->
                [ meta, [ ann, seq, dom, pfam, rep ] ]
            }
    }

    //
    // MODULES: Run UniprotR, UniProtExtractR, and rbioapi on extracted data
    //
    ch_annotate_input = EXTRACT_QUERY_PROTEINS.out.accessions
        .join(EXTRACT_QUERY_PROTEINS.out.uniprot_tsv, by: [0])
        .join(EXTRACT_QUERY_PROTEINS.out.mutations, by: [0])

    UNIPROT_ANNOTATE(ch_annotate_input)
    UNIPROT_EXTRACTR_ANNOTATE(ch_annotate_input)
    RBIOAPI_ANNOTATE(ch_annotate_input)

    emit:
    discovery        = EXTRACT_QUERY_PROTEINS.out.discovery        // [ meta, tsv ]
    accessions       = EXTRACT_QUERY_PROTEINS.out.accessions       // [ meta, txt ]
    uniprot_tsv      = EXTRACT_QUERY_PROTEINS.out.uniprot_tsv      // [ meta, tsv ]
    mutations        = EXTRACT_QUERY_PROTEINS.out.mutations        // [ meta, tsv ]
    query_summary    = EXTRACT_QUERY_PROTEINS.out.summary          // [ meta, json ]
    query_proteins   = EXTRACT_QUERY_PROTEINS.out.proteins         // [ meta, fasta ]
    uniprotr_results = UNIPROT_ANNOTATE.out.uniprotr_results      // [ meta, dir ]
    extractr_results = UNIPROT_EXTRACTR_ANNOTATE.out.extractr_results      // [ meta, dir ]
    rbioapi_results  = RBIOAPI_ANNOTATE.out.rbioapi_results        // [ meta, dir ]
    hmm_annotations  = ch_hmm_annotations                          // [ meta, tsv ]
    hmm_results      = ch_hmm_results                              // [ meta, [files] ]
}
