#!/usr/bin/env bash
# Run the PGIRL Nextflow pipeline with a local Ollama LLM.
# This helper starts Ollama, ensures the model is available, runs the pipeline,
# and stops Ollama afterwards.
#
# Usage:
#   scripts/run_with_ollama.sh [extra_nextflow_args...]
#
# Example:
#   scripts/run_with_ollama.sh --db_url postgresql://localhost:5432/pgirl

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

MODEL="qwen2.5:7b"

# Download Nextflow if it is not present (e.g., after cloning the repo).
if [ ! -x "${PROJECT_DIR}/nextflow" ]; then
    echo "Nextflow not found; downloading..."
    cd "$PROJECT_DIR"
    curl -fsSL https://get.nextflow.io | bash
fi

# Parse Nextflow's -log option so it can be placed before the `run` keyword.
NF_ARGS=()
NF_LOG=""
while [ $# -gt 0 ]; do
    case "$1" in
        -log)
            shift
            NF_LOG="${1:-}"
            shift
            ;;
        -log=*)
            NF_LOG="${1#*=}"
            shift
            ;;
        *)
            NF_ARGS+=("$1")
            shift
            ;;
    esac
done

# Check Ollama is installed.
if ! command -v ollama &>/dev/null; then
    echo "ERROR: ollama is not installed. Install it first:"
    echo "  curl -fsSL https://ollama.com/install.sh | sh"
    exit 1
fi

# Start Ollama server in the background only if it is not already running.
STARTED_OLLAMA=false
if curl -fsSL http://localhost:11434/api/tags &>/dev/null; then
    echo "Ollama is already running; using existing server."
else
    echo "Starting Ollama server..."
    ollama serve &
    OLLAMA_PID=$!
    STARTED_OLLAMA=true
fi

# Stop only the server this script started.
cleanup() {
    if [ "$STARTED_OLLAMA" = true ]; then
        echo "Stopping Ollama server (PID: $OLLAMA_PID)..."
        kill "$OLLAMA_PID" 2>/dev/null || true
        wait "$OLLAMA_PID" 2>/dev/null || true
    fi
}
trap cleanup EXIT

# Wait for the API to become ready.
for i in {1..30}; do
    if curl -fsSL http://localhost:11434/api/tags &>/dev/null; then
        break
    fi
    sleep 1
done

# Pull the model if it is not already present.
if ! ollama list | grep -q "^${MODEL}"; then
    echo "Pulling model ${MODEL} (this may take a few minutes)..."
    ollama pull "$MODEL"
fi

echo "Running Nextflow pipeline with LLM enabled..."
if [ -n "$NF_LOG" ]; then
    "${PROJECT_DIR}/nextflow" -log "$NF_LOG" run main.nf --use_llm true "${NF_ARGS[@]}"
else
    "${PROJECT_DIR}/nextflow" run main.nf --use_llm true "${NF_ARGS[@]}"
fi
