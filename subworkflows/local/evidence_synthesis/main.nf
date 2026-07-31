/*
 * Subworkflow: EVIDENCE_SYNTHESIS
 *
 * Synthesise pathogen-specific outputs into interpretable evidence:
 *   1. PHYLO_ANNOTATE — extract tip/node metadata and mutation matrices
 *   2. PHYLO_VISUALIZE — render the annotated tree + geographic map
 *
 * Input:  ch_auspice_results       - channel of [ meta, auspice_json, results_dir ]
 *         ch_rbioapi_results       - channel of [ meta, rbioapi_dir ]
 *         ch_nextclade_json_all    - path: all Nextclade JSONs (broadcast)
 *         ch_species_assignments   - path: species_assignments.tsv (broadcast)
 * Output: figures
 */

include { PHYLO_ANNOTATE   } from '../../../modules/local/phylo_annotate/main'
include { PHYLO_VISUALIZE  } from '../../../modules/local/phylo_visualize/main'

workflow EVIDENCE_SYNTHESIS {
    take:
    ch_auspice_results      // channel: [ val(meta), path(auspice_json), path(results_dir) ]
    ch_rbioapi_results      // channel: [ val(meta), path(rbioapi_dir) ]
    ch_nextclade_json_all   // path: all Nextclade JSONs (broadcast/value channel)
    ch_species_assignments  // path: species_assignments.tsv (broadcast/value channel)

    main:
    //
    // MODULE: Tree annotation (tip/node metadata + mutation matrices)
    //
    ch_annotate_input = ch_auspice_results
        .join(ch_rbioapi_results, by: [0])
        .map { meta, auspice_files, results_dir, rbioapi_dir ->
            def auspice = auspice_files instanceof List ? auspice_files[0] : auspice_files
            [ meta, auspice, results_dir, rbioapi_dir ]
        }

    PHYLO_ANNOTATE(
        ch_annotate_input,
        ch_nextclade_json_all,
        ch_species_assignments
    )

    //
    // MODULE: Tree and map visualization
    //
    ch_plot_base = PHYLO_ANNOTATE.out.tip_metadata
        .join(ch_auspice_results, by: [0])
        .map { meta, tip_meta, auspice_files, results_dir ->
            [ meta, results_dir, tip_meta ]
        }

    ch_plot_input = ch_plot_base
        .join(PHYLO_ANNOTATE.out.mutation_matrix, by: [0], remainder: true)
        .join(PHYLO_ANNOTATE.out.mutation_legend, by: [0], remainder: true)
        .join(PHYLO_ANNOTATE.out.protein_burden, by: [0], remainder: true)
        .map { meta, results, tip_meta, mut_matrix=null, mut_legend=null, protein_burden=null ->
            [ meta, results, tip_meta, mut_matrix ?: file('NO_FILE'), mut_legend ?: file('NO_FILE'), protein_burden ?: file('NO_FILE') ]
        }

    PHYLO_VISUALIZE(ch_plot_input)

    emit:
    figures = PHYLO_VISUALIZE.out.tree_heatmap.mix(PHYLO_VISUALIZE.out.geo_map)  // channel: [ meta, png ]
}
