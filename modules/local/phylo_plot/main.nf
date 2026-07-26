/*
 * Local module: PHYLO_PLOT
 *
 * Generate an annotated phylogenetic tree figure using ggtree/ggtreeExtra.
 * Branches are colored by country (majority rule), query tips are highlighted,
 * and a binary heatmap of top amino-acid mutations is appended.
 */

process PHYLO_PLOT {
    tag "$meta.id"
    label 'process_low'

    conda "${projectDir}/envs/pgirl_ggtree.yml"

    input:
    tuple val(meta), path(results_dir), path(tip_metadata), path(mutation_matrix)

    output:
    tuple val(meta), path("*.png"), emit: figure

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def title  = meta.title ?: "${meta.id} annotated phylogeny"
    def layout = meta.layout ?: 'rectangular'
    def mut_arg = mutation_matrix.name != 'NO_FILE' ? "--mutation-matrix=${mutation_matrix}" : ''
    """
    # Find tree.nwk inside results directory (follow symlinks with -L)
    TREE=\$(find -L ${results_dir} -name "tree.nwk" -path "*/all-outbreaks/*" | head -1)
    if [ -z "\${TREE}" ]; then
        TREE=\$(find -L ${results_dir} -name "tree.nwk" | head -1)
    fi

    plot_ggtreeExtra.R \\
        --tree=\${TREE} \\
        --metadata=${tip_metadata} \\
        ${mut_arg} \\
        --output=${prefix}_ggtreeExtra.png \\
        --title="${title}" \\
        --layout=${layout} \\
        ${args}
    """
}
