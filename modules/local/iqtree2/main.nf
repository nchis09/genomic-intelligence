/*
 * Local module: IQTREE2
 *
 * Build a maximum-likelihood phylogeny with IQ-TREE2.
 */

process IQTREE2 {
    tag "$meta.id"
    label 'process_high'

    conda "${projectDir}/envs/pgirl_iqtree.yml"

    input:
    tuple val(meta), path(alignment)

    output:
    tuple val(meta), path("${meta.id}.treefile"), emit: tree
    tuple val(meta), path("${meta.id}.iqtree"),   emit: report, optional: true
    tuple val(meta), path("${meta.id}.log"),      emit: log,    optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix  = "${meta.id}"
    def args    = task.ext.args ?: ''
    """
    iqtree -s ${alignment} -pre ${prefix} -nt ${task.cpus} -m GTR+I+G ${args}
    """
}
