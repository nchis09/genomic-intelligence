# Installation

Two environments are involved: the **pipeline environment** (Nextflow + conda,
runs the analyses) and the **dashboard environment** (R + Shiny, displays the
results). You need the pipeline env to produce results, and the dashboard env
to view them.

## 1. Pipeline environment

```bash
conda create -n nf-genomic-intelligence nextflow conda
conda activate nf-genomic-intelligence
```

The pipeline itself pulls its per-process dependencies via the `-profile
conda` flag — Nextflow creates the process environments on first run, so the
first execution is slower than subsequent ones.

## 2. Dashboard environment

The dashboard needs its own R environment (`pgirl_dashboard`):

```bash
conda env create -f envs/pgirl_dashboard.yml
conda activate pgirl_dashboard
```

!!! note "First launch auto-installs missing packages"
    If a package is missing, `app.R` installs it from CRAN / Bioconductor on
    first launch (shiny, bs4Dash, DT, plotly, leaflet, duckdb, ggtree,
    treeio, ...). That first start can take a few minutes; later launches are
    fast.

## 3. Optional: local LLM summaries

Some panels show AI-written narrative notes (tree notes, spread assessment,
intelligence brief). These are generated **locally** by
[Ollama](https://ollama.com) — no data leaves your machine.

```bash
# macOS
brew install ollama
ollama pull qwen2.5:14b     # default model configured in nextflow.config
```

The pipeline default is `qwen2.5:14b` at `http://127.0.0.1:11434` (see
`params.ollama_host` / `params.ollama_model` in `nextflow.config`). Run the
pipeline through `bin/run_with_ollama.sh` and Ollama is started/stopped for
you. If Ollama is unavailable, the pipeline still writes every narrative —
just as a deterministic **template** note instead of an LLM note (the
dashboard marks which is which).

## What gets installed where

| Component | Where it runs | Environment |
|---|---|---|
| Nextflow pipeline | your machine or HPC | `nf-genomic-intelligence` |
| Shiny dashboard | your machine, in a browser | `pgirl_dashboard` |
| Ollama (optional) | your machine (`localhost:11434`) | system install |
| This docs site | GitHub Pages (edit locally with `mkdocs serve`) | any Python ≥3.9 |
