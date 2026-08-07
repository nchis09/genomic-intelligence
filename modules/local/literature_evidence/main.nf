/*
 * Local module: LITERATURE_EVIDENCE
 *
 * Extracts structured, quotable, confidence-scored evidence from full-text
 * .txt files using a local llama-cpp model and a per-domain template.
 */

process LITERATURE_EVIDENCE {
    tag "$meta.id - $meta.species - $meta.domain"
    label 'process_low'

    conda "${projectDir}/envs/pgirl_agentslr.yml"

    input:
    tuple val(meta), path("*.txt")

    output:
    tuple val(meta), path("*.json"), optional: true, emit: evidence
    path "evidence_summary.json", emit: summary

    when:
    !params.skip_literature_evidence && !params.skip_literature_text && (task.ext.when == null || task.ext.when)

    script:
    def gguf_arg = params.llm_gguf ? "--gguf ${params.llm_gguf}" : ""
    def n_ctx_arg = params.llm_n_ctx ? "--n-ctx ${params.llm_n_ctx}" : ""
    def temp_arg = params.llm_temperature ? "--temperature ${params.llm_temperature}" : ""
    def seed_arg = params.llm_seed ? "--seed ${params.llm_seed}" : ""
    """
    [ -n "\${CONDA_PREFIX}" ] && export PATH="\${CONDA_PREFIX}/bin:\${PATH}"

    python ${projectDir}/bin/extract_literature_evidence.py \
        --input-dir . \
        --outdir . \
        --species "${meta.species}" \
        --domain "${meta.domain}" \
        --templates-yml ${projectDir}/database/evidence_templates.yml \
        ${gguf_arg} \
        ${n_ctx_arg} \
        ${temp_arg} \
        ${seed_arg}
    """
}
