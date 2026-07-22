#!/usr/bin/env bash
# =============================================================================
# PGIRL Nextclade Classification & Analysis Pipeline
# =============================================================================
# Screens unknown viral sequences against all Nextclade datasets,
# assigns pathogen/species, runs full Nextclade analysis per sample,
# then auto-routes to pathogen-specific phylogenetic pipelines.
#
# Usage:
#   ./run_nextclade.sh <input.fasta> <metadata.tsv> <output_dir> [--datasets d1,d2,...] [--skip-phylo]
#   ./run_nextclade.sh <input.fasta> <metadata.tsv> <output_dir> --step screen|assign|analyze|phylo
#
# If --datasets is omitted, screens against all representative datasets.
# =============================================================================

set -euo pipefail

NEXTSTRAIN_BIN="${HOME}/.nextstrain/runtimes/conda/env/bin"
NEXTCLADE="${NEXTSTRAIN_BIN}/nextclade"
SNAKEMAKE="${NEXTSTRAIN_BIN}/snakemake"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EBOLA_REPO="${SCRIPT_DIR}/nextstrain_ebola"

# --- Representative screening panel (one per pathogen/species) ---
DEFAULT_DATASETS=(
    "nextstrain/orthoebolavirus/ebov"
    "nextstrain/orthoebolavirus/bdbv"
    "nextstrain/orthoebolavirus/sudv"
    "nextstrain/sars-cov-2/wuhan-hu-1/orfs"
    "nextstrain/flu/h1n1pdm/ha/MW626062"
    "nextstrain/flu/h3n2/ha/EPI1857216"
    "nextstrain/flu/vic/ha/KX058884"
    "nextstrain/mpox/all-clades"
    "nextstrain/dengue/all"
    "nextstrain/measles/genome/WHO-2012"
    "nextstrain/rsv/a/EPI_ISL_412866"
    "nextstrain/rsv/b/EPI_ISL_1653999"
    "nextstrain/yellow-fever/prM-E"
    "nextstrain/hmpv/all-clades/NC_039199"
    "nextstrain/wnv/all-lineages"
)

# --- Parse arguments ---
if [[ $# -lt 3 ]]; then
    echo "Usage: $0 <input.fasta> <metadata.tsv> <output_dir> [--datasets dataset1,dataset2,...] [--skip-phylo]"
    exit 1
fi

INPUT_FASTA="$1"
INPUT_METADATA="$2"
OUTPUT_DIR="$3"
shift 3

DATASETS=("${DEFAULT_DATASETS[@]}")
SKIP_PHYLO=false
RUN_STEP="all"
RUN_SAMPLE=""
RUN_SPECIES=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --datasets)
            IFS=',' read -ra DATASETS <<< "$2"
            shift 2
            ;;
        --skip-phylo)
            SKIP_PHYLO=true
            shift
            ;;
        --step)
            RUN_STEP="$2"
            shift 2
            ;;
        --sample)
            RUN_SAMPLE="$2"
            shift 2
            ;;
        --species)
            RUN_SPECIES="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# --- Setup output directories ---
NCLADE_DIR="${OUTPUT_DIR}/nextclade_classification"
PHYLO_DIR="${OUTPUT_DIR}/nextstrain_ebola"
SCREEN_DIR="${NCLADE_DIR}/.screening"
mkdir -p "${SCREEN_DIR}" "${NCLADE_DIR}" "${PHYLO_DIR}"

echo "============================================"
echo "PGIRL Nextclade Classification Pipeline"
echo "============================================"
echo "Input:    ${INPUT_FASTA}"
echo "Metadata: ${INPUT_METADATA}"
echo "Output:   ${OUTPUT_DIR}"
echo "Datasets: ${#DATASETS[@]}"
echo ""

# =============================================================================
# STEP 1: Screen against all datasets
# =============================================================================
if [[ "${RUN_STEP}" == "all" || "${RUN_STEP}" == "screen" ]]; then
echo "--- Step 1: Screening against ${#DATASETS[@]} datasets ---"
echo ""

for dataset in "${DATASETS[@]}"; do
    short_name=$(echo "$dataset" | tr '/' '_')
    echo -n "  Screening: ${dataset} ... "

    "${NEXTCLADE}" run \
        --dataset-name "${dataset}" \
        --output-tsv "${SCREEN_DIR}/${short_name}.tsv" \
        "${INPUT_FASTA}" 2>/dev/null && echo "done" || echo "failed (skipping)"
done

echo ""
fi  # end step screen

# =============================================================================
# STEP 2: Assign pathogen/species per sequence (sample-centric output)
# =============================================================================
if [[ "${RUN_STEP}" == "all" || "${RUN_STEP}" == "assign" ]]; then
echo "--- Step 2: Assigning pathogen/species per sample ---"
echo ""

python3 - "${SCREEN_DIR}" "${NCLADE_DIR}/assignments.tsv" <<'PYTHON'
import sys, csv
from pathlib import Path

screen_dir = Path(sys.argv[1])
output_file = sys.argv[2]

# Dataset name to short species label mapping
SPECIES_LABELS = {
    'nextstrain_orthoebolavirus_ebov': 'ebov',
    'nextstrain_orthoebolavirus_bdbv': 'bdbv',
    'nextstrain_orthoebolavirus_sudv': 'sudv',
    'nextstrain_sars-cov-2_wuhan-hu-1_orfs': 'sars-cov-2',
    'nextstrain_flu_h1n1pdm_ha_MW626062': 'flu-h1n1pdm',
    'nextstrain_flu_h3n2_ha_EPI1857216': 'flu-h3n2',
    'nextstrain_flu_vic_ha_KX058884': 'flu-vic',
    'nextstrain_mpox_all-clades': 'mpox',
    'nextstrain_dengue_all': 'dengue',
    'nextstrain_measles_genome_WHO-2012': 'measles',
    'nextstrain_rsv_a_EPI_ISL_412866': 'rsv-a',
    'nextstrain_rsv_b_EPI_ISL_1653999': 'rsv-b',
    'nextstrain_yellow-fever_prM-E': 'yellow-fever',
    'nextstrain_hmpv_all-clades_NC_039199': 'hmpv',
    'nextstrain_wnv_all-lineages': 'wnv',
}

# Pathogen family mapping (for routing to phylogenetic pipelines)
PATHOGEN_FAMILY = {
    'ebov': 'orthoebolavirus',
    'bdbv': 'orthoebolavirus',
    'sudv': 'orthoebolavirus',
}

scores = {}  # {seqName: [(dataset_key, score, status, coverage, clade, outbreak)]}

for tsv_file in sorted(screen_dir.glob("*.tsv")):
    dataset_key = tsv_file.stem
    try:
        with open(tsv_file) as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                seq = row.get('seqName', '')
                if not seq:
                    continue
                try:
                    score = float(row.get('qc.overallScore', '999999999') or '999999999')
                except ValueError:
                    score = 999999999.0
                status = row.get('qc.overallStatus', 'unknown')
                try:
                    coverage = float(row.get('coverage', '0') or '0')
                except ValueError:
                    coverage = 0.0
                clade = row.get('clade', '') or ''
                outbreak = row.get('outbreak', '') or ''

                if seq not in scores:
                    scores[seq] = []
                scores[seq].append((dataset_key, score, status, coverage, clade, outbreak))
    except Exception as e:
        pass  # silently skip unparseable datasets

assignments = []
for seq, results in sorted(scores.items()):
    best = min(results, key=lambda x: x[1])
    dataset_key, score, status, coverage, clade, outbreak = best
    species_label = SPECIES_LABELS.get(dataset_key, dataset_key)
    family = PATHOGEN_FAMILY.get(species_label, species_label)
    # Sample-centric folder name: SAMPLENAME_species
    folder_name = f"{seq}_{species_label}"

    assignments.append({
        'seqName': seq,
        'assigned_dataset': dataset_key,
        'species_label': species_label,
        'pathogen_family': family,
        'folder_name': folder_name,
        'qc_score': score,
        'qc_status': status,
        'coverage': coverage,
        'clade': clade,
        'outbreak': outbreak,
    })
    print(f"  {seq} → {species_label} (QC={status}, score={score:.2f}, coverage={coverage:.4f}, clade={clade}, outbreak={outbreak})")

with open(output_file, 'w', newline='') as f:
    writer = csv.DictWriter(f,
        fieldnames=['seqName', 'assigned_dataset', 'species_label', 'pathogen_family',
                    'folder_name', 'qc_score', 'qc_status', 'coverage', 'clade', 'outbreak'],
        delimiter='\t')
    writer.writeheader()
    writer.writerows(assignments)

print(f"\n  Assignments written to: {output_file}")
PYTHON

echo ""
fi  # end step assign

# =============================================================================
# STEP 3: Full Nextclade analysis per sample (sample-centric folders)
# =============================================================================
if [[ "${RUN_STEP}" == "all" || "${RUN_STEP}" == "analyze" ]]; then
echo "--- Step 3: Full Nextclade analysis per sample ---"
echo ""

while IFS=$'\t' read -r seqName assigned_dataset species_label pathogen_family folder_name qc_score qc_status coverage clade outbreak; do
    # Skip header
    [[ "$seqName" == "seqName" ]] && continue

    # If --sample is set, only process that sample
    if [[ -n "${RUN_SAMPLE}" && "${seqName}" != "${RUN_SAMPLE}" ]]; then
        continue
    fi

    sample_dir="${NCLADE_DIR}/${folder_name}"
    mkdir -p "${sample_dir}"

    # Convert dataset key back to dataset name
    dataset_name=$(echo "$assigned_dataset" | sed 's/_/\//1; s/_/\//1; s/_/\//1')

    # Extract this sample's sequence
    temp_fasta="${sample_dir}/input_sequence.fasta"
    python3 -c "
seq_name = '${seqName}'
writing = False
with open('${INPUT_FASTA}') as f, open('${temp_fasta}', 'w') as out:
    for line in f:
        if line.startswith('>'):
            name = line[1:].strip().split()[0]
            writing = (name == seq_name)
        if writing:
            out.write(line)
"

    echo "  ${seqName} → ${species_label}"
    echo "    Dataset: ${dataset_name}"

    "${NEXTCLADE}" run \
        --dataset-name "${dataset_name}" \
        --output-all "${sample_dir}/" \
        "${temp_fasta}" 2>/dev/null

    echo "    Output:  ${sample_dir}/"
    echo ""

done < "${NCLADE_DIR}/assignments.tsv"
fi  # end step analyze

# =============================================================================
# STEP 4a: Generate nextstrain/ebola config + input files
# =============================================================================
if [[ "${RUN_STEP}" == "all" || "${RUN_STEP}" == "phylo" || "${RUN_STEP}" == "phylo_config" ]]; then
if [[ "${SKIP_PHYLO}" == "true" ]]; then
    echo "--- Step 4a: Phylogenetic config SKIPPED (--skip-phylo) ---"
    echo ""
else
    echo "--- Step 4a: Generating nextstrain/ebola config ---"
    echo ""

    # Collect ebolavirus samples and group by species
    EBOLA_SPECIES=$(awk -F'\t' 'NR>1 && $4=="orthoebolavirus" {print $3}' "${NCLADE_DIR}/assignments.tsv" | sort -u)

    if [[ -n "${EBOLA_SPECIES}" ]]; then
        echo "  Ebolavirus samples detected. Generating nextstrain/ebola config..."
        echo ""

        # Build the config.yaml for nextstrain/ebola
        python3 - "${NCLADE_DIR}/assignments.tsv" "${PHYLO_DIR}/config.yaml" "${NCLADE_DIR}" "${INPUT_METADATA}" <<'PYCONFIG'
import sys, csv, os
from collections import defaultdict

assignments_file = sys.argv[1]
config_file = sys.argv[2]
results_dir = sys.argv[3]
metadata_file = sys.argv[4]

# S3 public data URLs per species
PUBLIC_DATA = {
    'ebov': {
        'name': 'ppx_open_ebov',
        'metadata': 's3://nextstrain-data/files/workflows/ebola/ebov/metadata_open.tsv.zst',
        'sequences': 's3://nextstrain-data/files/workflows/ebola/ebov/sequences_open.fasta.zst',
    },
    'bdbv': {
        'name': 'ppx_open_bdbv',
        'metadata': 's3://nextstrain-data/files/workflows/ebola/bdbv/metadata_open.tsv.zst',
        'sequences': 's3://nextstrain-data/files/workflows/ebola/bdbv/sequences_open.fasta.zst',
    },
    'sudv': {
        'name': 'ppx_open_sudv',
        'metadata': 's3://nextstrain-data/files/workflows/ebola/sudv/metadata_open.tsv.zst',
        'sequences': 's3://nextstrain-data/files/workflows/ebola/sudv/sequences_open.fasta.zst',
    },
}

# Read assignments
species_samples = defaultdict(list)
with open(assignments_file) as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        if row['pathogen_family'] == 'orthoebolavirus':
            species_samples[row['species_label']].append(row)

# Build config
config = {'inputs': [], 'builds': []}

for species, samples in sorted(species_samples.items()):
    # Add public data for this species
    if species in PUBLIC_DATA:
        pub = PUBLIC_DATA[species]
        config['inputs'].append({
            'name': pub['name'],
            'species': species,
            'metadata': pub['metadata'],
            'sequences': pub['sequences'],
        })

    # Add local PGIRL samples for this species
    # Combine all samples of same species into one FASTA + metadata
    config['inputs'].append({
        'name': f'pgirl_{species}',
        'species': species,
        'metadata': f'input/{species}_metadata.tsv',
        'sequences': f'input/{species}_sequences.fasta',
    })

    # Add build
    config['builds'].append(f'{species}/all-outbreaks')

    # Write combined FASTA and metadata for this species
    input_dir = os.path.join(os.path.dirname(config_file), 'input')
    os.makedirs(input_dir, exist_ok=True)

    fasta_out = os.path.join(input_dir, f'{species}_sequences.fasta')
    meta_out = os.path.join(input_dir, f'{species}_metadata.tsv')

    with open(fasta_out, 'w') as fout:
        for sample in samples:
            src = os.path.join(results_dir, sample['folder_name'], 'input_sequence.fasta')
            if os.path.exists(src):
                with open(src) as fin:
                    fout.write(fin.read())

    # Read user-provided metadata and extract rows for this species' samples
    sample_names = {s['seqName'] for s in samples}
    meta_rows = []
    with open(metadata_file) as mf:
        reader = csv.DictReader(mf, delimiter='\t')
        meta_header = reader.fieldnames
        for row in reader:
            if row.get('accession', row.get('strain', '')) in sample_names:
                meta_rows.append(row)

    with open(meta_out, 'w') as fout:
        fout.write('\t'.join(meta_header) + '\n')
        for row in meta_rows:
            fout.write('\t'.join(row.get(h, '') for h in meta_header) + '\n')
        # If any samples are missing from metadata, add minimal rows
        found = {row.get('accession', row.get('strain', '')) for row in meta_rows}
        for name in sorted(sample_names - found):
            fout.write(f'{name}\t{name}\t\t\t\t\t\t\n')

# Write YAML config
with open(config_file, 'w') as f:
    # Write manually for clean formatting
    f.write('# Auto-generated by PGIRL classification pipeline\n')
    f.write('# Species assignments from Nextclade screening\n\n')
    f.write('inputs:\n')
    for inp in config['inputs']:
        f.write(f'  - name: {inp["name"]}\n')
        f.write(f'    species: {inp["species"]}\n')
        f.write(f'    metadata: "{inp["metadata"]}"\n')
        f.write(f'    sequences: "{inp["sequences"]}"\n')
    f.write('\nbuilds:\n')
    for build in config['builds']:
        f.write(f'  - {build}\n')

print(f"    Config:   {config_file}")
print(f"    Builds:   {', '.join(config['builds'])}")
for species in sorted(species_samples):
    samples = species_samples[species]
    print(f"    {species}: {', '.join(s['seqName'] for s in samples)}")
PYCONFIG

    else
        echo "  No ebolavirus samples found — skipping nextstrain/ebola."
        echo "  (Other pathogen phylogenetic pipelines not yet configured)"
        echo ""
    fi
fi
fi  # end step phylo_config

# =============================================================================
# STEP 4b: Run nextstrain/ebola phylogenetic build (per species if --species)
# =============================================================================
if [[ "${RUN_STEP}" == "all" || "${RUN_STEP}" == "phylo" || "${RUN_STEP}" == "phylo_build" ]]; then
if [[ "${SKIP_PHYLO}" == "true" ]]; then
    echo "--- Step 4b: Phylogenetic build SKIPPED (--skip-phylo) ---"
    echo ""
elif [[ ! -f "${PHYLO_DIR}/config.yaml" ]]; then
    echo "--- Step 4b: No config.yaml found — no ebolavirus samples to build ---"
    echo ""
else
    # Determine which builds to run
    if [[ -n "${RUN_SPECIES}" ]]; then
        BUILD_TARGET="${RUN_SPECIES}/all-outbreaks"
        echo "--- Step 4b: Building nextstrain/ebola — ${RUN_SPECIES} ---"
    else
        BUILD_TARGET=""
        echo "--- Step 4b: Building nextstrain/ebola — all species ---"
    fi
    echo ""

    # Auto-discover the nextstrain ebola pathogen directory (hash changes on updates)
    EBOLA_PATHOGEN_BASE=$(find "${HOME}/.nextstrain/pathogens/ebola" -maxdepth 1 -type d -name "main=*" 2>/dev/null | head -1)
    EBOLA_PHYLO_DIR="${EBOLA_PATHOGEN_BASE}/phylogenetic"

    if [[ -d "${EBOLA_PHYLO_DIR}" ]]; then
        # Copy input files into Snakemake's working directory
        cp -r "${PHYLO_DIR}/input" "${EBOLA_PHYLO_DIR}/"

        # Build a per-species or full config
        if [[ -n "${RUN_SPECIES}" ]]; then
            # Create a species-specific config from the full config
            python3 - "${PHYLO_DIR}/config.yaml" "${EBOLA_PHYLO_DIR}/pgirl_config.yaml" "${RUN_SPECIES}" <<'PYFILTER'
import sys, yaml

full_config = sys.argv[1]
out_config = sys.argv[2]
species = sys.argv[3]

with open(full_config) as f:
    lines = f.readlines()

# Simple filter: keep only inputs and builds for this species
with open(out_config, 'w') as f:
    f.write(f'# Auto-generated for species: {species}\n\n')
    f.write('inputs:\n')
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.strip().startswith('- name:'):
            # Read the next 3 lines (name, species, metadata, sequences)
            block = [lines[i+j] for j in range(4) if i+j < len(lines)]
            if any(f'species: {species}' in b for b in block):
                for b in block:
                    f.write(b)
            i += 4
        else:
            i += 1
    f.write(f'\nbuilds:\n')
    f.write(f'  - {species}/all-outbreaks\n')
PYFILTER
        else
            cp "${PHYLO_DIR}/config.yaml" "${EBOLA_PHYLO_DIR}/pgirl_config.yaml"
        fi

        echo "  Running nextstrain/ebola phylogenetic pipeline..."
        echo ""

        eval "$("${HOME}/.nextstrain/cli-standalone/nextstrain" init-shell zsh 2>/dev/null || true)"
        nextstrain build --cpus 4 "${EBOLA_PHYLO_DIR}" --configfile pgirl_config.yaml --forceall 2>&1 | tail -10

        echo ""

        # Move results to our output directory
        mkdir -p "${PHYLO_DIR}/auspice" "${PHYLO_DIR}/results"
        if [[ -n "${RUN_SPECIES}" ]]; then
            # Move only this species' results
            if [[ -d "${EBOLA_PHYLO_DIR}/results/${RUN_SPECIES}" ]]; then
                mv "${EBOLA_PHYLO_DIR}/results/${RUN_SPECIES}" "${PHYLO_DIR}/results/" 2>/dev/null || true
            fi
            mv "${EBOLA_PHYLO_DIR}/auspice/ebola_${RUN_SPECIES}_all-outbreaks.json" "${PHYLO_DIR}/auspice/" 2>/dev/null || true
        else
            for species_dir in ${EBOLA_PHYLO_DIR}/results/*/; do
                species=$(basename "${species_dir}")
                if grep -q "${species}" "${PHYLO_DIR}/config.yaml" 2>/dev/null; then
                    mv "${species_dir}" "${PHYLO_DIR}/results/" 2>/dev/null || true
                fi
            done
            mv "${EBOLA_PHYLO_DIR}"/auspice/ebola_*_all-outbreaks.json "${PHYLO_DIR}/auspice/" 2>/dev/null || true
        fi

        # Clean up copied files from pipeline directory
        rm -rf "${EBOLA_PHYLO_DIR}/input" "${EBOLA_PHYLO_DIR}/pgirl_config.yaml"

        echo "  Phylogenetic results saved to: ${PHYLO_DIR}/"
    else
        echo "  WARNING: nextstrain/ebola not found at ${EBOLA_PHYLO_DIR}"
    fi
    echo ""
fi
fi  # end step phylo_build

# =============================================================================
# Summary
# =============================================================================
echo "============================================"
echo "COMPLETE"
echo "============================================"
echo ""
echo "Results:"
echo "  Assignments:  ${NCLADE_DIR}/assignments.tsv"
echo "  Per-sample:   ${NCLADE_DIR}/<SAMPLE>_<species>/"
echo ""
echo "Per-sample outputs include:"
echo "  - nextclade.tsv    (mutations, clade, QC)"
echo "  - nextclade.aligned.fasta (aligned sequences)"
echo "  - nextclade.cds_translation.*.fasta (protein translations)"
echo "  - nextclade.auspice.json (phylogenetic placement)"
echo "  - nextclade.nwk (Newick tree with placement)"
echo ""
if [[ -f "${PHYLO_DIR}/config.yaml" ]]; then
    echo "Phylogenetic config generated at:"
    echo "  ${PHYLO_DIR}/config.yaml"
    echo ""
fi
