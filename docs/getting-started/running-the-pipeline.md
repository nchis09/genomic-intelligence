# Running the pipeline

Everything the dashboard shows comes out of one Nextflow run. This page walks
through the command used for the worked example in this guide (an Ebola /
influenza mixed FASTA with metadata), then the options that matter most for
the visualizations.

## The example run

```bash
./bin/run_with_ollama.sh main.nf -profile conda \
  --fasta input/all/all_fasta.fasta \
  --metadata input/all/all_metadata.tsv \
  --outdir ./results \
  --epi_min_year 2015 --epi_max_year 2026 \
  --skip_literature true --skip_iqtree true \
  -resume
```

Step by step:

1. **`./bin/run_with_ollama.sh main.nf ...`** — the wrapper starts a local
   Ollama server (if one isn't already running), polls until it answers, runs
   Nextflow with the remaining arguments, then stops only the server it
   started. Without it, LLM narrative notes fall back to template text.
2. **`-profile conda`** — per-process conda environments are created
   automatically.
3. **`--fasta` / `--metadata`** — the query genomes and their sample metadata.
   The pipeline classifies every sequence (Nextclade), groups them by
   detected species, and runs the per-species branches (Ebola via
   nextstrain-ebola, seasonal flu via pgirl builds, avian flu via
   segment-focused builds).
4. **`--outdir ./results`** — all dashboard data lives under here. Keep it
   stable; the dashboard points at this directory.
5. **`--epi_min_year 2015 --epi_max_year 2026`** — the epidemiological window.
   !!! warning "Don't anchor this to your query's year"
       The historical-comparison panels only work when background genomes
       from **earlier years** are in scope. If your query is from 2026 but you
       fetch only 2025–2026 epi data, the "how did it behave historically"
       sections have nothing to compare against. A decade or more is a good
       default.
6. **`--skip_literature true --skip_iqtree true`** — optional speed-ups used in
   this example.
7. **`-resume`** — cached tasks are reused; only steps whose inputs changed
   re-run.

## What the run produces for the dashboard

| Output | Feeds |
|---|---|
| `results/pathogen_identification/<species>/` | Biological Threat tables & charts |
| `results/knowledge_warehouse/knowledge_warehouse.duckdb` | Trees, tips, samples, epi records — most analytical panels |
| `results/transmission_context/<species>/*.tsv` | Transmission & Spread, geographic panels |
| `results/pathogen_mutation_profile/<species>/` | Mutation landscape, PLS-DA |
| `results/intelligence_brief/<species>/intelligence_brief.json` | Overview narrative + Intelligence Brief |
| `results/*/auspice/*.json` | Phylogenetic trees, tip metadata (clades, dates, geography) |
| `results/epi_raw/` + epi tables | FluNet / HDX / EMPRES-i historical context |

## Typical durations

A mixed Ebola + influenza run of a few genomes takes roughly 45 minutes to
2 hours on a laptop, dominated by Nextstrain builds and sequence annotation.
With `-resume`, re-runs only redo what changed (usually
`BUILD_KNOWLEDGE_DB` → `COMPUTE_TRANSMISSION_CONTEXT` onward).
