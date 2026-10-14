#!/usr/bin/env bash
# run_with_ollama.sh
#
# Wrapper around `nextflow run` that starts a local Ollama server if one is
# not already listening, then stops only the server it started after the
# pipeline finishes (success or failure). This keeps the user's normal
# Nextflow command unchanged while giving the pipeline a local LLM to call.
#
# Usage:
#   ./bin/run_with_ollama.sh main.nf -profile conda ...
#
# The script passes all arguments through to nextflow unchanged.

set -euo pipefail

OLLAMA_HOST="${PG_OLLAMA_HOST:-http://127.0.0.1:11434}"
PID_FILE="${PWD}/.ollama.pid"

ollama_reachable() {
  curl -s -o /dev/null -w '%{http_code}' "${OLLAMA_HOST}" 2>/dev/null || true
}

start_ollama() {
  if command -v ollama >/dev/null 2>&1; then
    echo "Starting local Ollama server at ${OLLAMA_HOST} ..."
    nohup ollama serve > "${PWD}/.ollama.log" 2>&1 &
    echo $! > "${PID_FILE}"
    # Wait for the server to respond (up to 30s)
    for _ in $(seq 1 30); do
      if [[ "$(ollama_reachable)" == "200" ]]; then
        echo "Ollama server is ready."
        return 0
      fi
      sleep 1
    done
    echo "Ollama server did not become reachable; LLM notes may fall back to templates." >&2
  else
    echo "Ollama binary not found in PATH; LLM notes will fall back to templates." >&2
  fi
}

stop_ollama() {
  if [[ -f "${PID_FILE}" ]]; then
    pid=$(cat "${PID_FILE}")
    echo "Stopping local Ollama server (pid ${pid}) ..."
    kill "${pid}" 2>/dev/null || true
    rm -f "${PID_FILE}"
  fi
}

# Clean up the server we started on exit.
trap stop_ollama EXIT

if [[ "$(ollama_reachable)" != "200" ]]; then
  start_ollama
else
  echo "Ollama already reachable at ${OLLAMA_HOST}; not starting a new server."
fi

# Run Nextflow with the supplied arguments.
nextflow "$@"
