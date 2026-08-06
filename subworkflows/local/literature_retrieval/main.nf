/*
 * Subworkflow: LITERATURE_RETRIEVAL
 *
 * Runs the OpenAlex literature search for each detected species and for every
 * domain defined in the literature search terms YAML.
 */

include { LITERATURE_SEARCH } from '../../../modules/local/literature_search/main'

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
    } else {
        ch_lit_results = channel.empty()
    }

    emit:
    lit_results = ch_lit_results  // channel: [ meta, [ */results.* files ] ]
}
