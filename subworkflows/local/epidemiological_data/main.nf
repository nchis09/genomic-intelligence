/*
 * Subworkflow: EPIDEMIOLOGICAL_DATA
 *
 * Route each species group to the right epidemiological data source:
 *   - orthoebolavirus -> HDX (Ebola datasets)
 *   - influenza / avian_flu (human) -> WHO FluMart (FluNet)
 *   - avian h5 (animal)             -> FAO EMPRES-i
 */

include { FETCH_EPIDEMIOLOGICAL_DATA } from '../../../modules/local/fetch_epidemiological_data/main'
include { FETCH_FLUNET               } from '../../../modules/local/fetch_flunet/main'
include { FETCH_EMPRESI              } from '../../../modules/local/fetch_empresi/main'

workflow EPIDEMIOLOGICAL_DATA {
    take:
    ch_species_data  // channel: [ val(meta), path(fasta), path(metadata) ]

    main:
    if (params.skip_epi_data) {
        ch_epi_raw = ch_species_data
            .map { meta, fasta, metadata -> [ meta, file('NO_FILE_epi') ] }
        ch_search_summary = ch_species_data
            .map { meta, fasta, metadata -> [ meta, file('NO_FILE_epi_summary') ] }
    } else {
        // Split the incoming channel by pathogen family.
        ch_species_data
            .branch { meta, fasta, metadata ->
                ebola : meta.pathogen == 'orthoebolavirus'
                flu   : meta.pathogen == 'influenza' || meta.pathogen == 'avian_influenza'
                other : true
            }
            .set { ch_branched }

        // Ebola/orthoebolavirus -> HDX
        def disease_map = [ orthoebolavirus: 'ebola' ]
        ch_ebola_input = ch_branched.ebola
            .map { meta, fasta, metadata ->
                def search_term = params.epi_search_term ?: disease_map.get(meta.pathogen, meta.pathogen)
                [ meta, search_term, meta.species ]
            }
        FETCH_EPIDEMIOLOGICAL_DATA(ch_ebola_input)

        // Influenza / avian_flu -> WHO FluMart for all flu species
        ch_flu_input = ch_branched.flu
            .map { meta, fasta, metadata ->
                [ meta, meta.species ]
            }
        FETCH_FLUNET(ch_flu_input)

        // Avian H5 -> also fetch FAO EMPRES-i animal outbreaks
        ch_empresi_input = FETCH_FLUNET.out.epi_raw
            .filter { meta, epi_dir ->
                meta.pathogen == 'avian_influenza'
            }
            .map { meta, epi_dir ->
                [ meta, meta.species, epi_dir ]
            }
        FETCH_EMPRESI(ch_empresi_input)

        def no_file_epi     = file('NO_FILE_epi')
        def no_file_summary = file('NO_FILE_epi_summary')

        ch_other_epi = ch_branched.other
            .map { meta, fasta, metadata -> [ meta, no_file_epi ] }
        ch_other_summary = ch_branched.other
            .map { meta, fasta, metadata -> [ meta, no_file_summary ] }

        // For avian H5 use the EMPRES-augmented epi_data, otherwise use FluNet
        ch_flu_epi = FETCH_FLUNET.out.epi_raw
            .filter { meta, epi_dir ->
                meta.pathogen != 'avian_influenza'
            }
            .mix(FETCH_EMPRESI.out.epi_raw)

        ch_epi_raw = FETCH_EPIDEMIOLOGICAL_DATA.out.epi_raw
            .mix(ch_flu_epi)
            .mix(ch_other_epi)
        ch_search_summary = FETCH_EPIDEMIOLOGICAL_DATA.out.search_summary
            .mix(FETCH_FLUNET.out.search_summary)
            .mix(ch_other_summary)
    }

    emit:
    epi_raw        = ch_epi_raw        // channel: [ meta, epi_data_dir ]
    search_summary = ch_search_summary // channel: [ meta, summary_tsv ]
}
