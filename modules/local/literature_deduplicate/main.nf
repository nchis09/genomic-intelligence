/*
 * Local module: LITERATURE_DEDUPLICATE
 *
 * Removes papers without abstracts and deduplicates the per-paper JSON files
 * produced by PUBMED_METADATA. Writes one cleaned JSON per retained paper and
 * a small summary file for each species/domain.
 */

process LITERATURE_DEDUPLICATE {
    tag "$meta.id - $meta.species - $meta.domain"
    label 'process_low'

    conda "${projectDir}/envs/pgirl_revtools.yml"

    input:
    tuple val(meta), path("*.json")

    output:
    tuple val(meta), path("*.json"), emit: deduplicated

    when:
    !params.skip_literature_deduplication && (task.ext.when == null || task.ext.when)

    script:
    """
    [ -n "\${CONDA_PREFIX}" ] && export PATH="\${CONDA_PREFIX}/bin:\${PATH}"
    Rscript ${projectDir}/bin/deduplicate_literature.R \
        --input-dir . \
        --outdir . \
        --species "${meta.species}" \
        --domain "${meta.domain}"
    """
}
