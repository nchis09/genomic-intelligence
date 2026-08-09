/*
 * Local module: UNIPROT_EXTRACTR_ANNOTATE
 *
 * Run UniProtExtractR structured field extraction on the UniProtKB TSV
 * downloaded by EXTRACT_QUERY_PROTEINS.
 *
 * Uses cloned repo at tools/UniProtExtractR.
 */

process UNIPROT_EXTRACTR_ANNOTATE {
    tag "$meta.id"
    label 'process_low'

    conda "${projectDir}/envs/pgirl_uniprotr.yml"
    container null

    input:
    tuple val(meta), path(accessions_txt), path(uniprot_tsv), path(mutations_tsv)

    output:
    tuple val(meta), path("uniprotextractr_results/") , emit: extractr_results

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix       = task.ext.prefix ?: "${meta.id}"
    def extractr_dir = "${projectDir}/tools/UniProtExtractR"
    """
    mkdir -p uniprotextractr_results

    echo "=== Running UniProtExtractR annotation ==="
    Rscript ${projectDir}/bin/annotate_uniprotextractr.R \\
        --input ${uniprot_tsv} \\
        --extractr_dir ${extractr_dir} \\
        --outdir uniprotextractr_results \\
        --prefix ${prefix}
    """
}
