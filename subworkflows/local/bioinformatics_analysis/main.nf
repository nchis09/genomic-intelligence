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
include { MAFFT_ALIGN                  } from '../../../modules/local/mafft_align/main'
include { IQTREE2                      } from '../../../modules/local/iqtree2/main'

workflow BIOINFORMATICS_ANALYSIS {
    take:
    ch_species_data  // channel: [ val(meta), path(fasta), path(metadata) ]

    main:
    //
    // Real bioinformatics workhorse: Nextstrain Ebola Snakemake workflow
    //
    NEXTSTRAIN_EBOLA(ch_species_data)

    //
    // Build a model-aware ML tree from the Nextstrain subsampled sequences
    //
    MAFFT_ALIGN(NEXTSTRAIN_EBOLA.out.results_dir)
    IQTREE2(MAFFT_ALIGN.out.alignment)

    emit:
    auspice   = NEXTSTRAIN_EBOLA.out.auspice      // channel: [ meta, json ]
    results   = NEXTSTRAIN_EBOLA.out.results_dir  // channel: [ meta, dir ]
    alignment = MAFFT_ALIGN.out.alignment         // channel: [ meta, fasta ]
    tree      = IQTREE2.out.tree                  // channel: [ meta, newick ]
}
