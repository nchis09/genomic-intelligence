/*
 * Local module: PHYLO_VISUALIZE
 *
 * Rich static visualisation of Nextstrain phylogenetic outputs per species:
 *   1. Annotated tree + per-genome heatmap using the published plotTree R
 *      function (katholt/plotTree, wrapped in bin/plot_tree_with_heatmap.R)
 *   2. Static geographic distribution map using ggplot2 + maps
 *      (bin/plot_geo_map.R)
 */

process PHYLO_VISUALIZE {
    tag "$meta.id"
    label 'process_low'

    conda "${projectDir}/envs/pgirl_ggtree.yml"

    input:
    tuple val(meta), path(results_dir), path(tip_metadata), path(mutation_matrix)

    output:
    tuple val(meta), path("*_tree_heatmap.png"), emit: tree_heatmap, optional: true
    tuple val(meta), path("*_geo_map.png"),     emit: geo_map,     optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def colour_by = task.ext.colour_nodes_by ?: 'country'
    def width     = task.ext.width ?: 2400
    def height    = task.ext.height ?: 3000
    def mut_arg   = mutation_matrix.name != 'NO_FILE' ? "--mutation-matrix ${mutation_matrix}" : ''
    """
    # Locate the refined Nextstrain tree inside the species results directory
    TREE=\$(find -L ${results_dir} -name "tree.nwk" -path "*/all-outbreaks/*" | head -1)
    if [ -z "\${TREE}" ]; then
        TREE=\$(find -L ${results_dir} -name "tree.nwk" | head -1)
    fi

    if [ -z "\${TREE}" ]; then
        echo "ERROR: no tree.nwk found in ${results_dir}" >&2
        exit 1
    fi

    # 1) Annotated phylogeny + per-genome mutation heatmap
    plot_tree_with_heatmap.R \\
        --tree \${TREE} \\
        --metadata ${tip_metadata} \\
        ${mut_arg} \\
        --output ${prefix}_tree_heatmap.png \\
        --prefix ${prefix} \\
        --colour-nodes-by ${colour_by} \\
        --width ${width} \\
        --height ${height} \\
        ${args}

    # 2) Geographic distribution map
    plot_geo_map.R \\
        --metadata ${tip_metadata} \\
        --output ${prefix}_geo_map.png \\
        --prefix ${prefix} \\
        --colour-by ${colour_by} \\
        --width 2000 \\
        --height 1400
    """
}
