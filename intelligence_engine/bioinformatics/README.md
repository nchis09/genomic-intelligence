# PGIRL Bioinformatics — Established Workflow Integration

This directory integrates two established, community-maintained tools for
viral genomic analysis, replacing the custom module-by-module approach.

## Architecture

```
INPUT: Unknown consensus FASTA sequences
                    │
    ┌───────────────▼───────────────┐
    │   Nextclade (Classification)  │
    │                               │
    │ Step 1: Screen all datasets   │
    │ Step 2: Assign pathogen/species│
    │ Step 3: Full analysis:        │
    │   • Alignment                 │
    │   • Translation               │
    │   • Mutation calling           │
    │   • Clade/outbreak assignment │
    │   • QC                        │
    │   • Phylogenetic placement    │
    └───────────────┬───────────────┘
                    │
    ┌───────────────▼───────────────┐
    │   nextstrain/ebola            │
    │   (Snakemake)                 │
    │                               │
    │   • Subsampling               │
    │   • Tree (IQ-TREE)            │
    │   • Timetree                  │
    │   • Phylogeography            │
    │   • Ancestral reconstruction  │
    │   • Auspice export            │
    └───────────────┬───────────────┘
                    │
    ┌───────────────▼───────────────┐
    │  PGIRL Intelligence Layer     │
    │  (downstream)                 │
    └───────────────────────────────┘
```

## Pipeline Logic

1. **Nextclade** — Pathogen identification, classification, QC, and mutation calling
   - Screens input sequences against all supported pathogen datasets
   - Assigns species/pathogen based on best QC score + coverage
   - Runs full analysis (alignment, translation, mutations, clade) using the correct reference
   - Supported pathogens: Ebola (EBOV, BDBV, SUDV), SARS-CoV-2, Influenza (H1N1, H3N2, Vic),
     Mpox, Dengue, Measles, RSV, Yellow Fever, Marburg, West Nile, HMPV, and more

2. **nextstrain/ebola** — Phylogenetics and transmission dynamics
   - Builds de novo phylogenetic tree (IQ-TREE)
   - Infers timetree with molecular clock (TreeTime)
   - Reconstructs ancestral sequences and mutations per gene
   - Infers phylogeography (transmission routes)
   - Exports interactive Auspice visualization

## Directory Layout

```
bioinformatics/
├── nextclade/               # Cloned: nextstrain/nextclade (reference)
├── nextstrain_ebola/        # Cloned: nextstrain/ebola
├── viralrecon/              # Cloned: nf-core/viralrecon (reserved for FASTQ workflows)
├── envs/
│   ├── viralrecon.yml       # Conda env for Nextflow + nf-core (future use)
│   └── nextstrain.yml       # Conda env for Augur + Snakemake + tools
├── run_nextclade.sh         # Wrapper: classification + full analysis
└── README.md                # This file
```

## Setup

### Prerequisites

Nextstrain CLI must be installed (provides nextclade + augur + snakemake):

```bash
# Install Nextstrain CLI (already done)
nextstrain --version

# Nextclade binary location:
# /Users/christianndekezi/.nextstrain/runtimes/conda/env/bin/nextclade
```

### Verify installation

```bash
nextclade --version
nextstrain --version
```

## Running the Pipeline

### Nextclade — Classification + Full Analysis

Accepts **any viral consensus FASTA** (unknown pathogen):

```bash
./run_nextclade.sh input/sequences.fasta output/nextclade
```

This will:
1. Screen against all Nextclade datasets
2. Identify pathogen and species per sequence
3. Run full analysis with correct reference
4. Output: TSV (mutations, clade, QC), aligned FASTA, translations, tree placement

### Nextstrain — Phylogenetics & Genomic Epidemiology

After classification, run phylogenetics for the identified pathogen:

```bash
nextstrain build phylogenetic --cores 4
```

View results interactively:
```bash
nextstrain view phylogenetic/auspice/
```

## Key Differences from Old Pipeline

| Aspect | Old (bioinformatics_old/) | New |
|--------|--------------------------|-----|
| Classification | Manual species assignment | Auto-screening all Nextclade datasets |
| Scope | Ebola only | All Nextclade-supported pathogens (100+) |
| Mutation calling | Custom scripts | Nextclade (validated, published) |
| Phylogenetics | Nextclade placement tree | Full IQ-TREE + timetree + phylogeography |
| Visualization | Static PNGs | Interactive Auspice |
| Maintenance | Us | Community-maintained (nextstrain) |

## Notes

- The **PGIRL intelligence layer** (functional annotation, novelty detection,
  outbreak comparison, diagnostic primer checks, surveillance trends) remains
  downstream and consumes outputs from these pipelines.
- **viralrecon** is kept for future use with FASTQ input (raw reads → consensus).
- For new pathogens, Nextstrain has separate pathogen repos that can be cloned
  similarly (e.g., `nextstrain/dengue`, `nextstrain/avian-flu`).
- The old pipeline is preserved at `bioinformatics_old/` for reference.
