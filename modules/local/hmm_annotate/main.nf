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
    tuple val(meta), path(query_fasta), path(hmm_db_dir)

    output:
    tuple val(meta), path("*_hmm_annotations.tsv"), emit: hmm_annotations
    path("*_hmm_sequence_table.txt"),     emit: hmm_sequence_table
    path("*_hmm_domain_table.txt"),       emit: hmm_domain_table
    path("*_hmm_pfam_table.txt"),         emit: hmm_pfam_table
    path("*_hmm_report.txt"),             emit: hmm_report

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def evalue  = task.ext.evalue ?: params.hmm_evalue ?: 1e-5
    """
    hmmscan --noali \
        --tblout ${prefix}_hmm_sequence_table.txt \
        --domtblout ${prefix}_hmm_domain_table.txt \
        --pfamtblout ${prefix}_hmm_pfam_table.txt \
        -E ${evalue} --domE ${evalue} \
        ${hmm_db_dir}/hmm_db.hmm ${query_fasta} > ${prefix}_hmm_report.txt

    python3 ${projectDir}/bin/parse_hmmscan.py \
        --tblout ${prefix}_hmm_sequence_table.txt \
        --prefix ${prefix} \
        --evalue ${evalue}
    """
}
