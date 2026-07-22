#!/usr/bin/env bash
# =============================================================================
# PGIRL Setup Script
# =============================================================================
# Creates the 3 conda environments needed to run the pipeline and verifies
# that all required tools are available.
#
# Usage:
#   ./setup.sh                  # Install all 3 envs
#   ./setup.sh --env nextstrain # Install only the nextstrain env
#   ./setup.sh --check          # Just verify existing installation
#
# Environments:
#   pgirl_nextstrain  — Nextclade classification + nextstrain/ebola phylogenetics
#   pgirl_db          — Database queries, epi data gathering (future)
#   pgirl_analysis    — Evidence integration, intelligence synthesis (future)
#
# Prerequisites:
#   - conda or mamba installed
#   - Java 11+ (for Nextflow)
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENVS_DIR="${SCRIPT_DIR}/envs"

# --- Find conda/mamba ---
if command -v mamba &>/dev/null; then
    CONDA_CMD="mamba"
elif command -v conda &>/dev/null; then
    CONDA_CMD="conda"
else
    # Try common install locations
    for candidate in "$HOME/anaconda3/bin/conda" "$HOME/miniconda3/bin/conda" "/opt/miniconda3/bin/conda"; do
        if [ -x "$candidate" ]; then
            CONDA_CMD="$candidate"
            break
        fi
    done
    if [ -z "${CONDA_CMD:-}" ]; then
        echo "ERROR: conda or mamba not found. Please install conda first:"
        echo "  https://docs.conda.io/en/latest/miniconda.html"
        exit 1
    fi
fi

echo "============================================"
echo "PGIRL Setup"
echo "============================================"
echo "Using: ${CONDA_CMD}"
echo ""

# --- Parse arguments ---
TARGET_ENV="all"
CHECK_ONLY=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        --env)
            TARGET_ENV="$2"
            shift 2
            ;;
        --check)
            CHECK_ONLY=true
            shift
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# --- Environment definitions (compatible with bash 3.x / macOS) ---
ALL_ENVS="nextstrain db analysis"

env_file() {
    case "$1" in
        nextstrain) echo "${ENVS_DIR}/pgirl_nextstrain.yml" ;;
        db)         echo "${ENVS_DIR}/pgirl_db.yml" ;;
        analysis)   echo "${ENVS_DIR}/pgirl_analysis.yml" ;;
        *)          echo "" ;;
    esac
}

env_name() {
    case "$1" in
        nextstrain) echo "pgirl_nextstrain" ;;
        db)         echo "pgirl_db" ;;
        analysis)   echo "pgirl_analysis" ;;
        *)          echo "" ;;
    esac
}

# --- Create or update environments ---
create_env() {
    local key="$1"
    local yml; yml=$(env_file "$key")
    local name; name=$(env_name "$key")

    if [ -z "$yml" ] || ! [ -f "$yml" ]; then
        echo "  ERROR: ${yml} not found"
        return 1
    fi

    # Check if env already exists
    if ${CONDA_CMD} env list 2>/dev/null | grep -q "/${name}$\|/${name} "; then
        echo "  [${name}] Already exists — updating..."
        ${CONDA_CMD} env update -f "$yml" --prune -q
    else
        echo "  [${name}] Creating..."
        ${CONDA_CMD} env create -f "$yml" -q
    fi
    echo "  [${name}] Done ✓"
}

if [ "$CHECK_ONLY" = false ]; then
    echo "--- Creating conda environments ---"
    echo ""

    if [ "$TARGET_ENV" = "all" ]; then
        for key in $ALL_ENVS; do
            create_env "$key"
            echo ""
        done
    else
        if [ -z "$(env_file "$TARGET_ENV")" ]; then
            echo "ERROR: Unknown env '${TARGET_ENV}'. Options: nextstrain, db, analysis"
            exit 1
        fi
        create_env "$TARGET_ENV"
        echo ""
    fi
fi

# --- Clone and patch nextstrain/ebola ---
if [ "$CHECK_ONLY" = false ]; then
    echo "--- Setting up nextstrain/ebola pipeline ---"
    echo ""

    EBOLA_DIR="${SCRIPT_DIR}/intelligence_engine/bioinformatics/nextstrain_ebola"

    if [ -d "${EBOLA_DIR}/.git" ]; then
        echo "  [nextstrain_ebola] Already cloned — pulling latest..."
        git -C "${EBOLA_DIR}" pull --quiet 2>/dev/null || true
    else
        echo "  [nextstrain_ebola] Cloning from GitHub..."
        git clone --quiet https://github.com/nextstrain/ebola.git "${EBOLA_DIR}"
    fi

    # Apply patch: augur renamed _resolve_filepath → resolve_filepath
    PATCH_FILE="${EBOLA_DIR}/shared/vendored/snakemake/config.smk"
    if [ -f "${PATCH_FILE}" ] && grep -q "_resolve_filepath" "${PATCH_FILE}" 2>/dev/null; then
        echo "  [nextstrain_ebola] Applying augur compatibility patch..."
        sed -i'' -e 's/_resolve_filepath/resolve_filepath/g' "${PATCH_FILE}"
    fi

    echo "  [nextstrain_ebola] Done ✓"
    echo ""

    # Install nextstrain ebola pathogen (downloads to ~/.nextstrain/pathogens/ebola/)
    if [ -n "${NEXTSTRAIN_PREFIX:-}" ] && [ -x "${NEXTSTRAIN_PREFIX}/bin/nextstrain" ]; then
        NEXTSTRAIN_BIN="${NEXTSTRAIN_PREFIX}/bin/nextstrain"
    elif command -v nextstrain &>/dev/null; then
        NEXTSTRAIN_BIN="nextstrain"
    fi

    if [ -n "${NEXTSTRAIN_BIN:-}" ]; then
        if [ ! -d "${HOME}/.nextstrain/pathogens/ebola" ]; then
            echo "  Installing nextstrain ebola pathogen..."
            "${NEXTSTRAIN_BIN}" setup ebola 2>/dev/null || echo "  (manual setup may be needed: nextstrain setup ebola)"
        else
            echo "  [nextstrain ebola pathogen] ✓ already installed"
        fi

        # Apply same patch to installed pathogen
        INSTALLED_PATCH=$(find "${HOME}/.nextstrain/pathogens/ebola" -path "*/vendored/snakemake/config.smk" 2>/dev/null | head -1)
        if [ -n "${INSTALLED_PATCH}" ] && grep -q "_resolve_filepath" "${INSTALLED_PATCH}" 2>/dev/null; then
            echo "  [nextstrain ebola pathogen] Applying augur compatibility patch..."
            sed -i'' -e 's/_resolve_filepath/resolve_filepath/g' "${INSTALLED_PATCH}"
        fi
    fi

    echo ""
fi

# --- Verify installation ---
echo "--- Verifying installation ---"
echo ""

# Check each env exists
for key in $ALL_ENVS; do
    name=$(env_name "$key")
    if ${CONDA_CMD} env list 2>/dev/null | grep -q "/${name}$\|/${name} "; then
        echo "  [${name}] ✓ installed"
    else
        echo "  [${name}] ✗ NOT installed"
    fi
done

echo ""

# Check key tools in nextstrain env
NEXTSTRAIN_PREFIX=$(${CONDA_CMD} env list 2>/dev/null | grep "pgirl_nextstrain" | awk '{print $NF}')

if [ -n "${NEXTSTRAIN_PREFIX}" ]; then
    echo "--- Key tools (pgirl_nextstrain) ---"
    echo ""
    for tool in nextflow nextclade augur snakemake mafft iqtree; do
        bin="${NEXTSTRAIN_PREFIX}/bin/${tool}"
        if [ -x "$bin" ]; then
            version=$("$bin" --version 2>&1 | head -1 || echo "?")
            echo "  [${tool}] ✓ ${version}"
        else
            echo "  [${tool}] ✗ not found"
        fi
    done
    echo ""
fi

echo "============================================"
echo "Setup complete"
echo "============================================"
echo ""
echo "To run the bioinformatics pipeline:"
echo "  conda activate pgirl_nextstrain"
echo "  nextflow run main.nf"
echo ""
echo "To run with custom inputs:"
echo "  nextflow run main.nf --input_fasta input/input_FASTA.fasta --input_metadata input/metadata.tsv"
echo ""
