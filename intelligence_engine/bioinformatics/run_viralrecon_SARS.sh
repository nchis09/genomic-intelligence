#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# PGIRL — Run nf-core/viralrecon for QC, variant calling, and annotation
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJ_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"

# ── Defaults (override via CLI args or env vars) ─────────────────────────────
CONDA_ENV="${PGIRL_VIRALRECON_ENV:-pgirl_viralrecon}"
INPUT_SAMPLESHEET="${1:-}"
OUTDIR="${2:-${PROJ_ROOT}/output/viralrecon}"
GENOME="${3:-}"                        # e.g. NC_002549.1 for Ebola Zaire
FASTA_REF="${4:-}"                     # path to reference FASTA (if no --genome)
PLATFORM="${PGIRL_PLATFORM:-illumina}" # illumina or nanopore
PROTOCOL="${PGIRL_PROTOCOL:-metagenomic}" # amplicon or metagenomic
PROFILE="${PGIRL_NF_PROFILE:-docker}"  # docker, singularity, or conda

# ── Usage ────────────────────────────────────────────────────────────────────
if [ -z "$INPUT_SAMPLESHEET" ]; then
    echo "Usage: $0 <samplesheet.csv> [outdir] [genome_id] [fasta_ref]"
    echo ""
    echo "  samplesheet.csv  : nf-core samplesheet (sample,fastq_1,fastq_2)"
    echo "  outdir           : output directory (default: output/viralrecon)"
    echo "  genome_id        : NCBI accession for --genome (e.g. NC_002549.1)"
    echo "  fasta_ref        : path to reference FASTA (alternative to genome_id)"
    echo ""
    echo "Environment variables:"
    echo "  PGIRL_VIRALRECON_ENV  : conda env name (default: pgirl_viralrecon)"
    echo "  PGIRL_PLATFORM        : illumina|nanopore (default: illumina)"
    echo "  PGIRL_PROTOCOL        : amplicon|metagenomic (default: metagenomic)"
    echo "  PGIRL_NF_PROFILE      : docker|singularity|conda (default: docker)"
    exit 1
fi

# ── Activate conda env ───────────────────────────────────────────────────────
echo ">>> Activating conda environment: ${CONDA_ENV}"
eval "$(conda shell.bash hook)"
conda activate "${CONDA_ENV}"

# ── Build command ────────────────────────────────────────────────────────────
CMD="nextflow run ${SCRIPT_DIR}/viralrecon \
    --input ${INPUT_SAMPLESHEET} \
    --outdir ${OUTDIR} \
    --platform ${PLATFORM} \
    --protocol ${PROTOCOL} \
    -profile ${PROFILE}"

# Add genome reference
if [ -n "$GENOME" ]; then
    CMD="${CMD} --genome ${GENOME}"
elif [ -n "$FASTA_REF" ]; then
    CMD="${CMD} --fasta ${FASTA_REF}"
else
    echo "WARNING: No --genome or --fasta specified. Viralrecon may fail."
fi

# ── Run ──────────────────────────────────────────────────────────────────────
echo ">>> Running viralrecon..."
echo "    Command: ${CMD}"
echo ""
eval "${CMD}"

echo ""
echo ">>> viralrecon complete. Results in: ${OUTDIR}"
