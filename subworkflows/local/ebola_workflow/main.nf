/*
 * Subworkflow: EBOLA_WORKFLOW
 *
 * Pathogen-specific workflow for Ebola (orthoebolavirus) species groups:
 *   1. NEXTSTRAIN_EBOLA — run the Nextstrain Ebola Snakemake workflow per species
 *
 * Tree visualization (PHYLO_ANNOTATE + PHYLO_VISUALIZE) is deliberately NOT
 * run here anymore — it now runs in the main workflow AFTER
 * PHENOTYPE_ANNOTATION completes, since the tree figure incorporates
 * phenotype-derived annotations (see workflows/genomic-intelligence.nf).
 * Those two modules are pathogen-agnostic, so this subworkflow only needs
 * to build the tree and hand off auspice/results.
 *
 * Input:  ch_species_data - channel of [ meta, fasta, metadata ]
 *         where meta includes: id, pathogen, species
 *         (species is already "bdbv"/"ebov"/"sudv", matching the Nextstrain
 *         Ebola ingest data directory names — no renaming needed here)
 * Output: auspice, results
 */

include { NEXTSTRAIN_EBOLA } from '../../../modules/local/nextstrain_ebola/main'

workflow EBOLA_WORKFLOW {
    take:
    ch_species_data  // channel: [ val(meta), path(fasta), path(metadata) ]

    main:
    //
    // MODULE: Build Nextstrain Ebola phylogeny per species (bdbv/ebov/sudv)
    //
    NEXTSTRAIN_EBOLA(ch_species_data)

    emit:
    auspice = NEXTSTRAIN_EBOLA.out.auspice        // channel: [ meta, json ]
    results = NEXTSTRAIN_EBOLA.out.results_dir    // channel: [ meta, dir ]
}
