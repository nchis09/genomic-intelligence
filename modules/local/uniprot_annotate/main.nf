/*
 * Local module: UNIPROT_ANNOTATE
 *
 * Run UniprotR (broad functional annotation), UniProtExtractR
 * (structured field extraction), and rbioapi (mutation-level phenotype,
 * STRING interactions, Reactome pathways) on query-sample-specific
 * accessions and UniProtKB TSV downloaded by EXTRACT_QUERY_PROTEINS.
 *
 * Uses cloned repos at tools/UniprotR, tools/UniProtExtractR, tools/rbioapi.
 */

process UNIPROT_ANNOTATE {
    tag "$meta.id"
    label 'process_low'

    conda "${projectDir}/envs/pgirl_uniprotr.yml"
    container null

    input:
    tuple val(meta), path(accessions_txt), path(uniprot_tsv), path(mutations_tsv)

    output:
    tuple val(meta), path("uniprotr_results/")       , emit: uniprotr_results
    tuple val(meta), path("uniprotextractr_results/") , emit: extractr_results
    tuple val(meta), path("rbioapi_results/")         , emit: rbioapi_results

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix       = task.ext.prefix ?: "${meta.id}"
    def species      = meta.species ?: meta.id.replaceAll(/^.*_/, '')
    def uniprotr_dir = "${projectDir}/tools/UniprotR"
    def extractr_dir = "${projectDir}/tools/UniProtExtractR"
    def rbioapi_dir  = "${projectDir}/tools/rbioapi"
    """
    mkdir -p uniprotr_results uniprotextractr_results rbioapi_results

    echo "=== Running UniprotR annotation ==="
    Rscript ${projectDir}/bin/annotate_uniprotr.R \\
        --accessions ${accessions_txt} \\
        --uniprotr_dir ${uniprotr_dir} \\
        --outdir uniprotr_results \\
        --prefix ${prefix}

    echo "=== Running UniProtExtractR annotation ==="
    Rscript ${projectDir}/bin/annotate_uniprotextractr.R \\
        --input ${uniprot_tsv} \\
        --extractr_dir ${extractr_dir} \\
        --outdir uniprotextractr_results \\
        --prefix ${prefix}

    echo "=== Running rbioapi annotation (UniProt mutagenesis + STRING + Reactome) ==="
    Rscript ${projectDir}/bin/annotate_rbioapi.R \\
        --mutations ${mutations_tsv} \\
        --accessions ${accessions_txt} \\
        --rbioapi_dir ${rbioapi_dir} \\
        --species ${species} \\
        --outdir rbioapi_results \\
        --prefix ${prefix}
    """
}
