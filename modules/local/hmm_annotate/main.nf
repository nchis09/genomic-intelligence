/*
 * Local module: HMM_ANNOTATE
 *
 * Run hmmscan with the pressed VOGDB HMM library against the
 * reconstructed query-protein FASTA produced by EXTRACT_QUERY_PROTEINS.
 */

process HMM_ANNOTATE {
    tag "$meta.id"
    label 'process_low'

    conda "${projectDir}/envs/pgirl_hmmer.yml"
    container null

    input:
    tuple val(meta), path(query_fasta), path(vogdb_dir)

    output:
    tuple val(meta), path("*_hmm_annotations.tsv"), emit: hmm_annotations

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def evalue  = task.ext.evalue ?: params.hmm_evalue ?: 1e-5
    """
    hmmscan --noali --tblout ${prefix}_hmm_tblout.txt \
        -E ${evalue} --domE ${evalue} \
        ${vogdb_dir}/vog.hmm ${query_fasta} > /dev/null

    python3 ${projectDir}/bin/parse_hmmscan.py \
        --tblout ${prefix}_hmm_tblout.txt \
        --prefix ${prefix} \
        --evalue ${evalue}
    """
}
