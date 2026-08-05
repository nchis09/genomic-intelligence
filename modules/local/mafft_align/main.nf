/*
 * Local module: MAFFT_ALIGN
 *
 * Align the Nextstrain subsampled FASTA with MAFFT to prepare an MSA for IQ-TREE2.
 */

process MAFFT_ALIGN {
    tag "$meta.id"
    label 'process_medium'

    conda "${projectDir}/envs/pgirl_iqtree.yml"

    input:
    tuple val(meta), path(results_dir)

    output:
    tuple val(meta), path("${meta.id}_aligned.fasta"), emit: alignment

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix  = "${meta.id}"
    def args    = task.ext.args ?: '--retree 2'
    """
    SUBSAMPLED=\$(find -L ${results_dir} -name "subsampled.fasta" -print -quit)
    [[ -n "\$SUBSAMPLED" ]] || { echo "subsampled.fasta not found in ${results_dir}"; exit 1; }
    mafft $args --thread ${task.cpus} "\$SUBSAMPLED" > ${prefix}_aligned.fasta
    """
}
