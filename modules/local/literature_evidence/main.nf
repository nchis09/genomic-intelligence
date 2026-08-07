/*
 * Local module: LITERATURE_EVIDENCE
 *
 * Extracts structured, quotable, confidence-scored evidence from full-text
 * .txt files using AgentSLR's data_extraction_custom stage with per-domain
 * tool-call schemas derived from evidence_templates.yml.
 * Merges PubMed metadata (title, year, authors, journal, DOI) into output JSON.
 * Starts a local Ollama server for LLM inference if not already running.
 */

process LITERATURE_EVIDENCE {
    tag "$meta.id - $meta.species - $meta.domain"
    label 'process_low'

    conda "${projectDir}/envs/pgirl_agentslr.yml"

    input:
    tuple val(meta), path("*.txt"), path("*.json")
    path templates_yml

    output:
    tuple val(meta), path("*.json"), optional: true, emit: evidence

    when:
    !params.skip_literature_evidence && !params.skip_literature_text && (task.ext.when == null || task.ext.when)

    script:
    def config_json_arg = params.agentslr_config_json ? "--config-json ${params.agentslr_config_json}" : ""
    def metadata_arg = params.skip_pubmed_metadata ? "" : "--metadata-dir metadata"
    """
    [ -n "\${CONDA_PREFIX}" ] && export PATH="\${CONDA_PREFIX}/bin:\${PATH}"

    # --- Stage metadata JSONs into metadata/ subdirectory ---
    mkdir -p metadata
    mv *.json metadata/ 2>/dev/null || true

    # --- Start Ollama server for local LLM inference ---
    # Use system-level ollama (not conda) since conda's ollama lacks llama-server
    # Search common system paths first, bypassing conda's PATH
    OLLAMA_BIN=""
    for candidate in /opt/homebrew/bin/ollama /usr/local/bin/ollama /usr/bin/ollama; do
        if [ -x "\${candidate}" ]; then
            OLLAMA_BIN="\${candidate}"
            break
        fi
    done
    if [ -z "\${OLLAMA_BIN}" ]; then
        echo "[LITERATURE_EVIDENCE] ERROR: system ollama binary not found."
        echo "Install with: brew install ollama  (or download from ollama.com)"
        exit 1
    fi
    echo "[LITERATURE_EVIDENCE] Using ollama: \${OLLAMA_BIN}"
    export OLLAMA_MODELS="${params.ollama_models_dir}"
    mkdir -p "\${OLLAMA_MODELS}"

    # Check if Ollama is already running; if not, start it in background
    if ! curl -s http://localhost:11434/api/tags > /dev/null 2>&1; then
        echo "[LITERATURE_EVIDENCE] Starting Ollama server..."
        "\${OLLAMA_BIN}" serve > ollama_server.log 2>&1 &
        OLLAMA_PID=\$!
        # Wait for server to be ready (up to 30 seconds)
        for i in \$(seq 1 30); do
            if curl -s http://localhost:11434/api/tags > /dev/null 2>&1; then
                echo "[LITERATURE_EVIDENCE] Ollama server is ready."
                break
            fi
            sleep 1
        done
    else
        echo "[LITERATURE_EVIDENCE] Ollama server already running."
    fi

    # Pull model if not already cached
    if ! "\${OLLAMA_BIN}" list 2>/dev/null | grep -q "${params.agentslr_model_name}"; then
        echo "[LITERATURE_EVIDENCE] Pulling model ${params.agentslr_model_name}..."
        "\${OLLAMA_BIN}" pull "${params.agentslr_model_name}"
    fi

    # --- Run evidence extraction ---
    python ${projectDir}/bin/run_agentslr_extraction.py \
        --input-dir . \
        --outdir . \
        --species "${meta.species}" \
        --domain "${meta.domain}" \
        --agentslr-dir ${projectDir}/tools/AgentSLR \
        --templates-yml "${templates_yml}" \
        --model-name "${params.agentslr_model_name}" \
        --base-url "${params.agentslr_base_url}" \
        --api-key "${params.agentslr_api_key}" \
        --concurrency ${params.agentslr_extraction_concurrency} \
        --max-completion-tokens ${params.agentslr_max_completion_tokens} \
        ${config_json_arg} \
        ${metadata_arg}

    # --- Cleanup: stop Ollama server if we started it ---
    if [ -n "\${OLLAMA_PID:-}" ]; then
        echo "[LITERATURE_EVIDENCE] Stopping Ollama server (PID \${OLLAMA_PID})..."
        kill "\${OLLAMA_PID}" 2>/dev/null || true
        wait "\${OLLAMA_PID}" 2>/dev/null || true
    fi
    """
}
