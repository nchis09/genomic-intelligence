/*
 * Subworkflow: BIOINFORMATICS_ANALYSIS
 *
 * Pathogen-specific bioinformatics stage. Currently wraps the Nextstrain
 * Ebola Snakemake workflow. Additional tools can be added here in future
 * pathogen workflows.
 *
 * Input:  ch_species_data - channel of [ meta, fasta, metadata ]
 * Output: auspice, results
 */

include { NEXTSTRAIN_EBOLA             } from '../../../modules/local/nextstrain_ebola/main'

workflow BIOINFORMATICS_ANALYSIS {
    take:
    ch_species_data  // channel: [ val(meta), path(fasta), path(metadata) ]

    main:
    //
    // Real bioinformatics workhorse: Nextstrain Ebola Snakemake workflow
    //
    NEXTSTRAIN_EBOLA(ch_species_data)

    emit:
    auspice = NEXTSTRAIN_EBOLA.out.auspice      // channel: [ meta, json ]
    results = NEXTSTRAIN_EBOLA.out.results_dir  // channel: [ meta, dir ]
}
