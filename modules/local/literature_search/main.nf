/*
 * Local module: LITERATURE_SEARCH
 *
 * Search OpenAlex for each detected species and for every domain defined in
 * database/literature_search_terms.yml. One JSON and one TSV are written per
 * domain under <species>/<domain>/.
 */

process LITERATURE_SEARCH {
    tag "$meta.id - $meta.species"
    label 'process_low'

    // Use the env YAML; Nextflow builds and activates the conda env.
    // Call 'python' (not 'python3') because 'python3' on this system resolves
    // to the Homebrew interpreter before the activated conda env.
    conda "${projectDir}/envs/pgirl_literature.yml"

    input:
    tuple val(meta), path(fasta), path(metadata)

    output:
    tuple val(meta), path("*/results.*"), emit: results

    when:
    task.ext.when == null || task.ext.when

    script:
    def max_results = params.literature_search_max_results ?: 1000
    def mailto = params.literature_search_mailto ? "--mailto ${params.literature_search_mailto}" : ""
    def api_key = params.literature_search_api_key ? "--api-key ${params.literature_search_api_key}" : ""
    """
    python ${projectDir}/bin/literature_search.py \
        --species "${meta.species}" \
        --terms-yaml "${params.literature_search_terms}" \
        --max-results ${max_results} \
        ${mailto} \
        ${api_key} \
        --outdir .
    """
}
