/*
 * Local module: EVIDENCE_SUMMARY
 *
 * Generates a plain-language narrative per literature domain (plus a
 * species-level overview) from the extracted evidence TSVs, using the local
 * Ollama model. Output is domain_summaries.tsv — stored in the knowledge
 * warehouse so the dashboard displays it without calling the LLM itself.
 */

process EVIDENCE_SUMMARY {
    tag "$meta.id - $meta.species"
    label 'process_low'

    conda "${projectDir}/envs/pgirl_evidence.yml"

    input:
    tuple val(meta), path(evidence_files, stageAs: "evidence_?.tsv")

    output:
    tuple val(meta), path("domain_summaries.tsv"), emit: summaries
    tuple val(meta), path("evidence_highlights.json"), emit: highlights

    when:
    !params.skip_literature_evidence && !params.skip_literature_text && (task.ext.when == null || task.ext.when)

    script:
    def ollama_host  = task.ext.ollama_host  ?: params.ollama_host  ?: 'http://localhost:11434'
    def ollama_model = task.ext.ollama_model ?: params.ollama_model ?: ''
    def ollama_n_ctx = task.ext.ollama_n_ctx ?: params.ollama_n_ctx ?: '4096'
    def ollama_temp  = task.ext.ollama_temperature ?: params.ollama_temperature ?: '0.1'
    """
    [ -n "\${CONDA_PREFIX}" ] && export PATH="\${CONDA_PREFIX}/bin:\$PATH"

    python3 ${projectDir}/bin/summarize_evidence.py \
        --input-dir . \
        --outdir . \
        --species "${meta.species}" \
        --ollama-host "${ollama_host}" \
        --ollama-model "${ollama_model}" \
        --n-ctx ${ollama_n_ctx} \
        --temperature ${ollama_temp}
    """
}
