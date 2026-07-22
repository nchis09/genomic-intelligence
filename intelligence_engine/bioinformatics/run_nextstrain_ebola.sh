#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# PGIRL — Run Nextstrain Ebola phylogenetic workflow
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJ_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"
NEXTSTRAIN_DIR="${SCRIPT_DIR}/nextstrain_ebola"

# ── Defaults ─────────────────────────────────────────────────────────────────
CONDA_ENV="${PGIRL_NEXTSTRAIN_ENV:-pgirl_nextstrain}"
INPUT_FASTA="${1:-}"
INPUT_METADATA="${2:-}"
OUTDIR="${3:-${PROJ_ROOT}/output/nextstrain}"
CONFIGFILE="${4:-}"  # optional: custom snakemake config

# ── Usage ────────────────────────────────────────────────────────────────────
if [ -z "$INPUT_FASTA" ] || [ -z "$INPUT_METADATA" ]; then
    echo "Usage: $0 <sequences.fasta> <metadata.tsv> [outdir] [config.yaml]"
    echo ""
    echo "  sequences.fasta : consensus sequences in FASTA format"
    echo "  metadata.tsv    : sample metadata (strain, date, country, etc.)"
    echo "  outdir          : output directory (default: output/nextstrain)"
    echo "  config.yaml     : optional Snakemake config override"
    echo ""
    echo "Environment variables:"
    echo "  PGIRL_NEXTSTRAIN_ENV : conda env name (default: pgirl_nextstrain)"
    exit 1
fi

# ── Activate conda env ───────────────────────────────────────────────────────
echo ">>> Activating conda environment: ${CONDA_ENV}"
eval "$(conda shell.bash hook)"
conda activate "${CONDA_ENV}"

# ── Prepare input ────────────────────────────────────────────────────────────
# Nextstrain/ebola expects input in phylogenetic/data/ or via config
PHYLO_DIR="${NEXTSTRAIN_DIR}/phylogenetic"
DATA_DIR="${PHYLO_DIR}/data"
mkdir -p "${DATA_DIR}"

# Symlink or copy input files
echo ">>> Linking input data..."
ln -sf "$(realpath "${INPUT_FASTA}")" "${DATA_DIR}/sequences.fasta"
ln -sf "$(realpath "${INPUT_METADATA}")" "${DATA_DIR}/metadata.tsv"

# ── Build command ────────────────────────────────────────────────────────────
CMD="snakemake --snakefile ${PHYLO_DIR}/Snakefile \
    --cores all \
    --directory ${OUTDIR} \
    --config \
        sequences=${DATA_DIR}/sequences.fasta \
        metadata=${DATA_DIR}/metadata.tsv"

if [ -n "$CONFIGFILE" ]; then
    CMD="${CMD} --configfile ${CONFIGFILE}"
fi

# ── Run ──────────────────────────────────────────────────────────────────────
echo ">>> Running Nextstrain Ebola phylogenetic workflow..."
echo "    Command: ${CMD}"
echo ""
eval "${CMD}"

echo ""
echo ">>> Nextstrain complete. Results in: ${OUTDIR}"
echo "    View with: auspice view --datasetDir ${OUTDIR}/auspice"
