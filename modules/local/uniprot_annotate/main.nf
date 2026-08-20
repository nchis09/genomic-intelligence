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
    # Ensure UniprotR is cloned once in a concurrency-safe way.
    UNIPROTR_DIR="\$(${projectDir}/bin/setup_tool.sh UniprotR https://github.com/Proteomicslab57357/UniprotR.git)"

    mkdir -p uniprotr_results

    echo "=== Running UniprotR annotation ==="
    \$CONDA_PREFIX/bin/Rscript ${projectDir}/bin/annotate_uniprotr.R \\
        --accessions ${accessions_txt} \\
        --uniprotr_dir ${uniprotr_dir} \\
        --outdir uniprotr_results \\
        --prefix ${prefix}
    """
}
