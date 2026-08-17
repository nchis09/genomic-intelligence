# Genomic Intelligence Framework — Dashboard (v3, Intelligence Overview)

A standalone R Shiny app, built on **bs4Dash** + **fresh**, that reads the flat
TSV outputs already produced by `bin/pathogen_identification.R` and presents
them as an assessment-first "Intelligence Overview," alongside a roadmap of
the other intelligence objectives that don't have a wired data source yet.

It is **not** wired into the Nextflow pipeline — run the pipeline first, then
launch this app manually to browse the results.

## Layout

- **Intelligence Overview** (Home) — species selector (when more than one
  species is found), an assessment header for the selected species —
  **MONITOR / INVESTIGATE / ESCALATE**, Confidence, Evidence availability,
  and 3–5 Key Intelligence Signals explaining the assessment (see
  `modules/assessment.R`) — two CTAs ("Explore Evidence", "Generate Genomic
  Intelligence Brief"), and a de-emphasized roadmap strip for the other
  objectives. Only **Biological Threat** has real content today.
- **Sidebar** — Intelligence Overview, then **Biological Threat** (formerly
  "Pathogen Identification"; expands into a submenu of species discovered
  under `--outdir`, e.g. `BDBV`, `SUDV`), then the other 6 intelligence
  objectives (Transmission & Spread, Geographic & Temporal Context, Health
  Impact, Populations & Settings at Risk, Countermeasure Readiness, Risk &
  Action), then Evidence & Knowledge Gaps and Intelligence Brief. Everything
  except Biological Threat opens a clean "Planned" roadmap page (no data
  source wired in yet). Intelligence Brief additionally has a UI-only
  "Preview & Edit Brief" modal scaffold for its future export flow. *Compare
  is a distinct future capability and isn't in the nav yet.*
- **Footer** (every page) — a one-line "About the Genomic Intelligence
  Framework" caption (reusing MultiQC's `intro_text`), plus the
  supporting-institution logos from `assets/institution_logo/`.

Code layout:

```
dashboard/
  app.R                             # bs4Dash page, sidebar/body wiring, theme
  modules/
    home.R                          # Intelligence Overview UI + GIF intro-text helper
    assessment.R                    # MONITOR/INVESTIGATE/ESCALATE + confidence heuristic
    pathogen_identification.R       # Biological Threat: interactive comparison workspace (selectize + plotly charts + DT table)
    placeholder.R                   # roadmap objectives + "Planned" placeholder cards + Brief preview modal
```

## Requirements

Create the conda env (once):

```bash
conda env create -f ../envs/pgirl_dashboard.yml
conda activate pgirl_dashboard
```

## Launch

From the repo root:

```bash
Rscript -e "shiny::runApp('dashboard')"
```

This starts a local R web server and opens your browser automatically
(typically `http://127.0.0.1:<port>`). Leave the terminal running — closing it
ends the session. Re-run the command any time you want to view it again.

Alternatively, open `dashboard/app.R` in RStudio and click **Run App**.

## Usage

1. In the sidebar, set **Pipeline --outdir** to the `--outdir` you used for the
   Nextflow run (defaults to `../results`, relative to `dashboard/`).
2. On **Intelligence Overview**, pick a species (if more than one is found)
   to see its assessment, confidence, and Key Intelligence Signals. Click a
   signal, or the "Explore Evidence" button, to jump straight into that
   species' Biological Threat tables.
3. Under **Biological Threat** in the sidebar, select samples from the
   multi-select picker (up to 5) to see a summary strip and interactive
   plotly comparison charts (coverage, identity, genetic distance, QC score,
   clade distance, MSA identity). The full identification results table
   below has checkboxes synced with the picker.
4. The other sidebar items (Transmission & Spread, Geographic & Temporal
   Context, etc.) are roadmap pages for now — they describe what each
   objective will show once its data source is wired in. **Intelligence
   Brief** additionally has a "Preview & Edit Brief" button that opens a
   placeholder editable preview, previewing the intended export flow.

## Data source

Reads directly from:

```
<outdir>/pathogen_identification/<species>/species_identification/identification_summary.tsv
<outdir>/pathogen_identification/<species>/species_identification/unresolved_samples.tsv
```

No database connection — these are the same plain TSVs written by
`PATHOGEN_IDENTIFICATION` alongside the `*_mqc.tsv` MultiQC copies (see
`bin/pathogen_identification.R`). MultiQC output is untouched; this app is a
separate, parallel way to view the same data.

## Future upgrade path

When genuinely dynamic/filterable tables are needed (beyond what's
precomputed), the plan is to have the analysis stage write result tables
*into the DuckDB export* itself, and have this app query that file directly
(read-only) for fast, live filter/sort/drill-down — see
`.windsurf/plans/shiny-dashboard-v1-529106.md` for details.
