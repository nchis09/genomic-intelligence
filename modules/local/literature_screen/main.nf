/*
 * Local module: LITERATURE_SCREEN
 *
 * Runs a fully automated asreview simulation on the deduplicated paper JSONs.
 * Uses domain-specific keywords, a 10-year publication window, and species
 * synonyms/exclusions to provide prior labels for the asreview model.
 */

process LITERATURE_SCREEN {
    tag "$meta.id - $meta.species - $meta.domain"
    label 'process_low'

    conda "${projectDir}/envs/pgirl_asreview.yml"

    input:
    tuple val(meta), path("*.json")
    path terms_file

    output:
    tuple val(meta), path("*.json"), emit: screened

    when:
    !params.skip_literature_screening && (task.ext.when == null || task.ext.when)

    script:
    def n_prior_inc = params.asreview_n_prior_included ?: 5
    def n_prior_exc = params.asreview_n_prior_excluded ?: 5
    def n_stop = params.asreview_n_stop ?: 10
    def top_n = params.asreview_top_n ?: 50
    def min_year_arg = params.asreview_min_year ? "--min-year ${params.asreview_min_year}" : ""
    """
    [ -n "\${CONDA_PREFIX}" ] && export PATH="\${CONDA_PREFIX}/bin:\${PATH}"

    python -m pip install -q asreview

    Rscript ${projectDir}/bin/screen_literature.R \
        --input-dir . \
        --outdir . \
        --species "${meta.species}" \
        --domain "${meta.domain}" \
        --terms-file ${terms_file} \
        --n-prior-included ${n_prior_inc} \
        --n-prior-excluded ${n_prior_exc} \
        --n-stop ${n_stop} \
        --top-n ${top_n} \
        ${min_year_arg}
    """
}
