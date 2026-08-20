/*
 * Local module: LITERATURE_SEARCH
 *
 * Search Europe PMC for each detected species and for every domain defined in
 * database/literature_search_terms.yml. One JSON and one TSV are written per
 * domain under <species>/<domain>/.
 */

process LITERATURE_SEARCH {
    tag "$meta.id - $meta.species"
    label 'process_low'
    maxForks 1

    // Use the env YAML; Nextflow builds and activates the conda env.
    // Call 'python' (not 'python3') because 'python3' on this system resolves
    // to the Homebrew interpreter before the activated conda env.
    conda "${projectDir}/envs/pgirl_literature.yml"

    input:
    tuple val(meta), path(fasta), path(metadata)
    path terms_yaml
    path script

    output:
    tuple val(meta), path("*/results.*"), emit: results

    when:
    task.ext.when == null || task.ext.when

    script:
    def max_results = params.literature_search_max_results ?: 1000
    """
    [ -n "\${CONDA_PREFIX}" ] && export PATH="\${CONDA_PREFIX}/bin:\${PATH}"

    # Fallback: if the conda env did not expose python, locate the cached
    # literature env directly and prepend its bin directory to PATH.
    if ! command -v python >/dev/null 2>&1; then
        LIT_ENV=\$(find ${projectDir}/work/conda -maxdepth 1 -type d -name 'pgirl_literature-*' | head -n1)
        if [ -n "\${LIT_ENV}" ]; then
            export PATH="\${LIT_ENV}/bin:\${PATH}"
        fi
    fi

    python ${script} \
        --species "${meta.species}" \
        --terms-yaml "${terms_yaml}" \
        --max-results ${max_results} \
        --outdir .
    """
}
