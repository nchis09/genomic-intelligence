#!/usr/bin/env bash
# Convenience wrapper for setup.nf.
#
# Checks the runtime environment and sets up the PostgreSQL reference database
# in one step. Uses the active Conda/base Python and PostgreSQL binaries.
#
# Usage:
#   ./setup.sh
#   PGIRL_DB_URL=postgresql://localhost:5432/other_db ./setup.sh

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

exec ./nextflow run setup.nf "$@"
