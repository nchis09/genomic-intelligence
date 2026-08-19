/*
 * Local module: UNIPROT_ANNOTATE
 *
 * Run UniprotR broad functional annotation on query-sample-specific
 * accessions downloaded by EXTRACT_QUERY_PROTEINS.
 *
 * Uses cloned repo at tools/UniprotR.
 */

process UNIPROT_ANNOTATE {
    tag "$meta.id"
    label 'process_low'

    conda "${projectDir}/envs/pgirl_uniprotr.yml"
    container null

    input:
    tuple val(meta), path(accessions_txt), path(uniprot_tsv), path(mutations_tsv)

    output:
    tuple val(meta), path("uniprotr_results/") , emit: uniprotr_results

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix       = task.ext.prefix ?: "${meta.id}"
    def uniprotr_dir = "${projectDir}/tools/UniprotR"
    """
    # Auto-download UniprotR on first run if the local clone is empty/missing.
    if [ ! -f "${projectDir}/tools/UniprotR/DESCRIPTION" ]; then
        rm -rf "${projectDir}/tools/UniprotR"
        git clone https://github.com/Proteomicslab57357/UniprotR.git "${projectDir}/tools/UniprotR"
    fi

    mkdir -p uniprotr_results

    echo "=== Running UniprotR annotation ==="
    Rscript ${projectDir}/bin/annotate_uniprotr.R \\
        --accessions ${accessions_txt} \\
        --uniprotr_dir ${uniprotr_dir} \\
        --outdir uniprotr_results \\
        --prefix ${prefix}
    """
}
