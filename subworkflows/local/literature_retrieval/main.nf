/*
 * Subworkflow: LITERATURE_RETRIEVAL
 *
 * Runs the Europe PMC literature search for each detected species and for every
 * domain defined in the literature search terms YAML, then fetches PubMed
 * metadata (title, abstract, authors, journal, year, DOI, keywords, PMCID)
 * for each individual paper.
 */

include { LITERATURE_SEARCH } from '../../../modules/local/literature_search/main'
include { PUBMED_METADATA } from '../../../modules/local/pubmed_metadata/main'
include { LITERATURE_DEDUPLICATE } from '../../../modules/local/literature_deduplicate/main'
include { LITERATURE_SCREEN } from '../../../modules/local/literature_screen/main'
include { LITERATURE_PDF } from '../../../modules/local/literature_pdf/main'
include { LITERATURE_TEXT } from '../../../modules/local/literature_text/main'
include { LITERATURE_EVIDENCE } from '../../../modules/local/literature_evidence/main'

workflow LITERATURE_RETRIEVAL {
    take:
    ch_species_data  // channel: [ val(meta), path(fasta), path(metadata) ]

    main:
    if (!params.skip_literature_search) {
        LITERATURE_SEARCH(
            ch_species_data,
            file(params.literature_search_terms),
            file("${projectDir}/bin/literature_search.py")
        )
        ch_lit_results = LITERATURE_SEARCH.out.results

        if (!params.skip_pubmed_metadata) {
            ch_pubmed_input = ch_lit_results.flatMap { meta, files ->
                files.findAll { it.name == 'results.json' }.collect { json ->
                    def tsv = files.find { it.name == 'results.tsv' && it.parent == json.parent }
                    [meta + [domain: json.parent.name], json, tsv]
                }
            }
            PUBMED_METADATA(ch_pubmed_input)
            ch_pubmed_results = PUBMED_METADATA.out.metadata

            if (!params.skip_literature_deduplication) {
                LITERATURE_DEDUPLICATE(ch_pubmed_results)
                ch_pubmed_results = LITERATURE_DEDUPLICATE.out.deduplicated
            }

            if (!params.skip_literature_screening) {
                LITERATURE_SCREEN(ch_pubmed_results, file(params.literature_search_terms))
                ch_pubmed_results = LITERATURE_SCREEN.out.screened
            }

            if (!params.skip_literature_pdf) {
                LITERATURE_PDF(ch_pubmed_results)
                if (!params.skip_literature_text) {
                    LITERATURE_TEXT(LITERATURE_PDF.out.pdfs)
                    if (!params.skip_literature_evidence) {
                        LITERATURE_EVIDENCE(LITERATURE_TEXT.out.text)
                    }
                }
            }
        } else {
            ch_pubmed_results = channel.empty()
        }
    } else {
        ch_lit_results = channel.empty()
        ch_pubmed_results = channel.empty()
    }

    emit:
    lit_results = ch_lit_results  // channel: [ meta, [ */results.* files ] ]
    pubmed_metadata = ch_pubmed_results  // channel: [ meta, [ *.json files ] ]
}
