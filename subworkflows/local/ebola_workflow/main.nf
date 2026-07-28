/*
 * Subworkflow: EBOLA_WORKFLOW
 *
 * Pathogen-specific workflow for Ebola (orthoebolavirus) species groups:
 *   1. NEXTSTRAIN_EBOLA — run the Nextstrain Ebola Snakemake workflow per species
 *   2. PHYLO_ANNOTATE   — extract tip metadata & mutation matrix
 *   3. PHYLO_VISUALIZE  — static annotated tree + per-genome heatmap + map
 *
 * Input:  ch_species_data - channel of [ meta, fasta, metadata ]
 *         where meta includes: id, pathogen, species
 *         (species is already "bdbv"/"ebov"/"sudv", matching the Nextstrain
 *         Ebola ingest data directory names — no renaming needed here)
 * Output: figures, tip_metadata, auspice, results
 */

include { NEXTSTRAIN_EBOLA } from '../../../modules/local/nextstrain_ebola/main'
include { PHYLO_ANNOTATE   } from '../../../modules/local/phylo_annotate/main'
include { PHYLO_VISUALIZE  } from '../../../modules/local/phylo_visualize/main'

workflow EBOLA_WORKFLOW {
    take:
    ch_species_data  // channel: [ val(meta), path(fasta), path(metadata) ]

    main:
    //
    // MODULE: Build Nextstrain Ebola phylogeny per species (bdbv/ebov/sudv)
    //
    NEXTSTRAIN_EBOLA(ch_species_data)

    //
    // MODULE: Extract annotations from Nextstrain JSONs
    //
    // Pass the auspice JSON and results directory to PHYLO_ANNOTATE.
    // It will find tree.nwk, muts.json, branch_lengths.json inside results_dir.
    ch_annotate_input = NEXTSTRAIN_EBOLA.out.auspice
        .join(NEXTSTRAIN_EBOLA.out.results_dir, by: [0])
        .map { meta, auspice_files, results_dir ->
            def auspice = auspice_files instanceof List ? auspice_files[0] : auspice_files
            [ meta, auspice, results_dir ]
        }

    PHYLO_ANNOTATE(ch_annotate_input)

    //
    // MODULE: Static annotated tree + per-genome heatmap + geographic map
    //
    // Combine tip_metadata with results_dir and optional mutation_matrix
    ch_plot_base = PHYLO_ANNOTATE.out.tip_metadata
        .join(NEXTSTRAIN_EBOLA.out.results_dir, by: [0])
        .map { meta, tip_meta, results_dir ->
            [ meta, results_dir, tip_meta ]
        }

    // Add mutation_matrix if available, otherwise use NO_FILE placeholder
    ch_mut_matrix = PHYLO_ANNOTATE.out.mutation_matrix

    ch_plot_input = ch_plot_base
        .join(ch_mut_matrix, by: [0], remainder: true)
        .map { meta, results, tip_meta, mut_matrix=null ->
            [ meta, results, tip_meta, mut_matrix ?: file('NO_FILE') ]
        }

    PHYLO_VISUALIZE(ch_plot_input)

    emit:
    figures      = PHYLO_VISUALIZE.out.tree_heatmap
                       .mix(PHYLO_VISUALIZE.out.geo_map) // channel: [ meta, png ]
    tip_metadata = PHYLO_ANNOTATE.out.tip_metadata     // channel: [ meta, tsv ]
    auspice      = NEXTSTRAIN_EBOLA.out.auspice        // channel: [ meta, json ]
    results      = NEXTSTRAIN_EBOLA.out.results_dir    // channel: [ meta, dir ]
}
