#!/usr/bin/env bash
# Single entry point for running the PGIRL main.nf Nextflow workflow.
#
# Usage:
#   ./run_pipeline.sh [nextflow_args...]
#
# If --use_llm true is requested, this script automatically starts a local
# Ollama server (if not already running), pulls the configured model, runs
# Nextflow, and stops the Ollama server it started. Otherwise it runs
# Nextflow directly.
#
# Example:
#   ./run_pipeline.sh \
#       --input_fasta input/input_FASTA.fasta \
#       --input_metadata input/input_metadata.csv \
#       --db_url postgresql://localhost:5432/pgirl \
#       --use_llm true

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Download Nextflow if it is not present (e.g., after cloning the repo).
if [ ! -x "${SCRIPT_DIR}/nextflow" ]; then
    echo "Nextflow not found; downloading..."
    curl -fsSL https://get.nextflow.io | bash
fi

# Ensure the active Conda environment (or a local base install) is first in
# PATH so the correct Python and tools are used. This avoids accidentally picking
# up a Homebrew/system Python that lacks the project dependencies.
if [ -n "${CONDA_PREFIX:-}" ] && [ -d "${CONDA_PREFIX}/bin" ]; then
    export PATH="${CONDA_PREFIX}/bin:${PATH}"
else
    for conda_bin in "$HOME/anaconda3/bin" "$HOME/miniconda3/bin" "/opt/miniconda3/bin"; do
        if [ -d "$conda_bin" ]; then
            export PATH="$conda_bin:$PATH"
            break
        fi
    done
fi

# Defaults
USE_LLM=false
OUTDIR="output"

args=("$@")
i=0
while [ $i -lt ${#args[@]} ]; do
    case "${args[$i]}" in
        --use_llm)
            next="${args[$((i+1))]:-}"
            if [ "$next" = "true" ]; then
                USE_LLM=true
            fi
            ;;
        --use_llm=true)
            USE_LLM=true
            ;;
        --outdir)
            OUTDIR="${args[$((i+1))]:-output}"
            ;;
        --outdir=*)
            OUTDIR="${args[$i]#*=}"
            ;;
    esac
    i=$((i+1))
done

mkdir -p "$OUTDIR"
LOG_FILE="$OUTDIR/.nextflow.log"

if [ "$USE_LLM" = true ]; then
    export PGIRL_USE_LLM=true
    exec scripts/run_with_ollama.sh -log "$LOG_FILE" "$@"
else
    export PGIRL_USE_LLM=false
    exec ./nextflow -log "$LOG_FILE" run main.nf "$@"
fi
