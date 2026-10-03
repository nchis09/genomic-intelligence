/*
 * Local module: LITERATURE_EVIDENCE
 *
 * Extracts concrete, quotable evidence claims from full-text .txt paper files
 * using a local Ollama model. Outputs a per-domain TSV (evidence_extracted.tsv)
 * and an extraction log (extraction_log.json) instead of the previous
 * per-paper JSONs.
 */

process LITERATURE_EVIDENCE {
    tag "$meta.id - $meta.species - $meta.domain"
    label 'process_low'

    conda "${projectDir}/envs/pgirl_evidence.yml"

    input:
    tuple val(meta), path("*.txt"), path("*.json")
    path templates_yml

    output:
    tuple val(meta), path("evidence_extracted.tsv"), emit: tsv
    tuple val(meta), path("extraction_log.json"), emit: log

    when:
    !params.skip_literature_evidence && !params.skip_literature_text && (task.ext.when == null || task.ext.when)

    script:
    def ollama_host  = task.ext.ollama_host  ?: params.ollama_host  ?: 'http://localhost:11434'
    def ollama_model = task.ext.ollama_model ?: params.ollama_model ?: ''
    def ollama_n_ctx = task.ext.ollama_n_ctx ?: params.ollama_n_ctx ?: '4096'
    def ollama_temp  = task.ext.ollama_temperature ?: params.ollama_temperature ?: '0.1'
    """
    [ -n "\${CONDA_PREFIX}" ] && export PATH="\${CONDA_PREFIX}/bin:\$PATH"

    # --- Stage metadata JSONs into metadata/ subdirectory ---
    mkdir -p metadata
    mv *.json metadata/ 2>/dev/null || true

    # --- Run Ollama evidence extraction ---
    set -e
    python3 ${projectDir}/bin/extract_literature_evidence.py \
        --input-dir . \
        --metadata-dir metadata \
        --outdir . \
        --species "${meta.species}" \
        --domain "${meta.domain}" \
        --templates-yml "${templates_yml}" \
        --ollama-host "${ollama_host}" \
        --ollama-model "${ollama_model}" \
        --n-ctx ${ollama_n_ctx} \
        --temperature ${ollama_temp}
    """
}
