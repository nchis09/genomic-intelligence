/*
 * Local module: EVIDENCE_QC
 *
 * Validates per-PMID extraction JSONs produced by LITERATURE_EVIDENCE using
 * Pandera schema validation and a per-paper completeness score. Passing
 * JSONs are copied to clean/ and failing JSONs (with qc_failure annotated)
 * are copied to failed/; an aggregate qc_report.json and qc_report.html are
 * also generated per species/domain.
 */

process EVIDENCE_QC {
    tag "$meta.id - $meta.species - $meta.domain"
    label 'process_low'

    conda "${projectDir}/envs/pgirl_qc.yml"

    input:
    tuple val(meta), path("*.json")
    path templates_yml

    output:
    tuple val(meta), path("clean/*.json"), optional: true, emit: clean
    tuple val(meta), path("qc_report.json"), emit: report_json
    path "qc_report.html", emit: report_html

    when:
    !params.skip_evidence_qc && !params.skip_literature_evidence && !params.skip_literature_text && !params.skip_literature_pdf

    script:
    """
    [ -n "\${CONDA_PREFIX}" ] && export PATH="\${CONDA_PREFIX}/bin:\${PATH}"

    set -e
    python ${projectDir}/bin/run_evidence_qc.py \
        --input-dir . \
        --outdir . \
        --species "${meta.species}" \
        --domain "${meta.domain}" \
        --templates-yml "${templates_yml}" \
        --min-completeness ${params.evidence_qc_min_completeness}
    """
}
