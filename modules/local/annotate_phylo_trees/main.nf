/*
 * Local module: ANNOTATE_PHYLO_TREES
 *
 * Produce fully-annotated phylogenetic tree figures (PNG + SVG) for a species
 * using ggtree/ggtreeExtra. Reads trees and tip metadata from the knowledge
 * warehouse DuckDB export. Generates:
 *   1. Augur evolutionary tree (if nextstrain/auspice tree exists)
 *   2. IQ-TREE bootstrap tree (if iqtree2 tree exists)
 *
 * Both light and dark mode variants are produced.
 */

process ANNOTATE_PHYLO_TREES {
    tag "$meta.id"
    label 'process_low'

    conda "${projectDir}/envs/pgirl_ggtree.yml"

    input:
    tuple val(meta), path(duckdb_file)

    output:
    tuple val(meta), path("*.png"), emit: png, optional: true
    tuple val(meta), path("*.svg"), emit: svg, optional: true

    script:
    def species = meta.species ?: meta.id
    """
    export PATH="\$CONDA_PREFIX/bin:\$PATH"

    # Light mode
    Rscript ${projectDir}/bin/annotate_phylo_trees.R \\
        --duckdb ${duckdb_file} \\
        --species ${species} \\
        --outdir .

    # Dark mode variant
    Rscript ${projectDir}/bin/annotate_phylo_trees.R \\
        --duckdb ${duckdb_file} \\
        --species ${species} \\
        --outdir . \\
        --dark-mode
    """
}
