/*
 * Subworkflow: EPIDEMIOLOGICAL_DATA
 *
 * Search HDX for a disease term using the rhdx R package and download
 * the first usable resource per species group.
 */

include { FETCH_EPIDEMIOLOGICAL_DATA } from '../../../modules/local/fetch_epidemiological_data/main'

workflow EPIDEMIOLOGICAL_DATA {
    take:
    ch_species_data  // channel: [ val(meta), path(fasta), path(metadata) ]

    main:
    // Map detected pathogen family to a canonical HDX search term.
    // params.epi_search_term can override the mapping.
    def disease_map = [
        orthoebolavirus: 'ebola'
    ]

    ch_search_input = ch_species_data
        .map { meta, fasta, metadata ->
            def search_term = params.epi_search_term ?: disease_map.get(meta.pathogen, meta.pathogen)
            [ meta, search_term, meta.species ]
        }

    if (!params.skip_epi_data) {
        FETCH_EPIDEMIOLOGICAL_DATA(ch_search_input)
        ch_epi_raw = FETCH_EPIDEMIOLOGICAL_DATA.out.epi_raw
        ch_search_summary = FETCH_EPIDEMIOLOGICAL_DATA.out.search_summary
    } else {
        ch_epi_raw = ch_search_input.map { meta, _term, _species -> [ meta, file('NO_FILE') ] }
        ch_search_summary = ch_search_input.map { meta, _term, _species -> [ meta, file('NO_FILE') ] }
    }

    emit:
    epi_raw        = ch_epi_raw        // channel: [ meta, epi_data.csv ]
    search_summary = ch_search_summary // channel: [ meta, rhdx_search_results.tsv ]
}
