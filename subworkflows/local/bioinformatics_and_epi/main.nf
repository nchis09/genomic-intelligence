/*
 * Subworkflow: BIOINFORMATICS_AND_EPIDEMIOLOGICAL
 *
 * Wraps the per-species Nextstrain Ebola build and the HDX epidemiology
 * download so the parent EBOLA_WORKFLOW only fans to one child for these
 * two stages. This keeps the nf-metro fan to at most two.
 *
 * Input:  ch_species_data - channel of [ meta, fasta, metadata ]
 * Output: auspice, results, auspice_results, epi_raw, epi_search_summary
 */

include { BIOINFORMATICS_ANALYSIS } from '../bioinformatics_analysis/main'
include { EPIDEMIOLOGICAL_DATA    } from '../epidemiological_data/main'

workflow BIOINFORMATICS_AND_EPIDEMIOLOGICAL {
    take:
    ch_species_data  // channel: [ val(meta), path(fasta), path(metadata) ]

    main:
    BIOINFORMATICS_ANALYSIS(ch_species_data)
    EPIDEMIOLOGICAL_DATA(ch_species_data)

    ch_auspice_results = BIOINFORMATICS_ANALYSIS.out.auspice
        .join(BIOINFORMATICS_ANALYSIS.out.results, by: [0])
        .map { meta, auspice_files, results_dir ->
            def auspice = auspice_files instanceof List ? auspice_files[0] : auspice_files
            [ meta, auspice, results_dir ]
        }

    emit:
    auspice        = BIOINFORMATICS_ANALYSIS.out.auspice      // channel: [ meta, json ]
    results        = BIOINFORMATICS_ANALYSIS.out.results      // channel: [ meta, dir ]
    alignment      = BIOINFORMATICS_ANALYSIS.out.alignment    // channel: [ meta, fasta ]
    tree           = BIOINFORMATICS_ANALYSIS.out.tree         // channel: [ meta, newick ]
    auspice_results = ch_auspice_results                      // channel: [ meta, json, dir ]
    epi_raw        = EPIDEMIOLOGICAL_DATA.out.epi_raw         // channel: [ meta, csv ]
    epi_search_summary = EPIDEMIOLOGICAL_DATA.out.search_summary  // channel: [ meta, tsv ]
}
