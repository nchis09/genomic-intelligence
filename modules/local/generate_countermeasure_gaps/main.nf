/*
 * Local module: GENERATE_COUNTERMEASURE_GAPS
 *
 * Aggregates per-domain evidence_extracted.tsv files for a species into
 * countermeasure_readiness.tsv and knowledge_gaps.tsv.
 */

process GENERATE_COUNTERMEASURE_GAPS {
    tag "$meta.id"
    label 'process_low'

    conda "${projectDir}/envs/pgirl_dashboard.yml"

    input:
    tuple val(meta), path("evidence/*")

    output:
    tuple val(meta), path("countermeasure_readiness.tsv"), emit: countermeasures
    tuple val(meta), path("knowledge_gaps.tsv"), emit: gaps

    when:
    !params.skip_countermeasure_gaps && !params.skip_literature_evidence && (task.ext.when == null || task.ext.when)

    script:
    """
    [ -n "\${CONDA_PREFIX}" ] && export PATH="\${CONDA_PREFIX}/bin:\$PATH"
    set -e

    if ls evidence/*.tsv 1>/dev/null 2>&1; then
        Rscript ${projectDir}/bin/aggregate_countermeasures_gaps.R \
            --species "${meta.species}" \
            --evidence-dir evidence \
            --outdir .
    else
        echo "No evidence TSVs found; writing empty countermeasure/gap tables"
        touch countermeasure_readiness.tsv
        touch knowledge_gaps.tsv
    fi
    """
}
