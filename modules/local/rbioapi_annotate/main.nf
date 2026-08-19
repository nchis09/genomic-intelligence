/*
 * Local module: RBIOAPI_ANNOTATE
 *
 * Run rbioapi for mutation-level phenotype evidence, STRING interactions,
 * and Reactome pathways on query-sample-specific accessions and mutations.
 *
 * Uses cloned repo at tools/rbioapi.
 */

process RBIOAPI_ANNOTATE {
    tag "$meta.id"
    label 'process_low'

    conda "${projectDir}/envs/pgirl_uniprotr.yml"
    container null

    input:
    tuple val(meta), path(accessions_txt), path(uniprot_tsv), path(mutations_tsv)

    output:
    tuple val(meta), path("rbioapi_results/") , emit: rbioapi_results

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def species     = meta.species ?: meta.id.replaceAll(/^.*_/, '')
    def rbioapi_dir = "${projectDir}/tools/rbioapi"
    """
    # Auto-download rbioapi on first run if the local clone is empty/missing.
    if [ ! -f "${projectDir}/tools/rbioapi/DESCRIPTION" ]; then
        rm -rf "${projectDir}/tools/rbioapi"
        git clone https://github.com/moosa-r/rbioapi.git "${projectDir}/tools/rbioapi"
    fi

    mkdir -p rbioapi_results

    echo "=== Running rbioapi annotation (UniProt mutagenesis + STRING + Reactome) ==="
    Rscript ${projectDir}/bin/annotate_rbioapi.R \\
        --mutations ${mutations_tsv} \\
        --accessions ${accessions_txt} \\
        --uniprot_tsv ${uniprot_tsv} \\
        --rbioapi_dir ${rbioapi_dir} \\
        --species ${species} \\
        --outdir rbioapi_results \\
        --prefix ${prefix}
    """
}
