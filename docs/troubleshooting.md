# Troubleshooting

Panels that look wrong almost always trace to one of a handful of causes.
Check these first.

## "No species found" on every page

The app can't see `<outdir>/pathogen_identification/`. Fix the
**Pipeline --outdir** field in the sidebar, or run the pipeline. Details:
[Launching the dashboard](getting-started/launching-the-dashboard.md).

## Narrative notes look generic ("template" source)

The LLM wasn't reachable when `GENERATE_TREE_NOTES` /
`GENERATE_SPREAD_ASSESSMENT` / `GENERATE_INTELLIGENCE_BRIEF` ran — the
pipeline writes a deterministic template instead of skipping the section.

- **Fix**: make sure Ollama is installed with a pulled model, then run
  through `./bin/run_with_ollama.sh` (it starts/stops the server for you).
- Verify locally: `ollama list` and `curl http://127.0.0.1:11434/api/tags`.

## Transmission & Spread: "How did it behave historically?" is empty

The section needs background genomes with a resolvable **strain** (clade →
lineage → outbreak), a **date**, and **coordinates**. Symptom-to-cause map:

- **No strains at all** → background tips have no clade in the warehouse.
  Influenza JSONs store clade under `clade_membership` (seasonal) or
  `label_clade` (avian), not `clade` — fixed in `load_background_from_auspice`;
  rebuild the warehouse.
- **Timeline collapses to one month** → tip dates were missing. Auspice
  stores dates as `num_date` (decimal year), not `date` — the loader now
  converts both; rebuild.
- **Map draws nothing / spread stats blank** → background samples had no
  `sample_geo_location` link or coordinates. Now populated from
  `database/country_centroids.csv` at load; rebuild.

After any of these fixes, `-resume` re-runs `BUILD_KNOWLEDGE_DB` and
`COMPUTE_TRANSMISSION_CONTEXT` and the section fills in.

## Mutation landscape is flat (all entropy ≈ 0)

The landscape needs **translated protein FASTAs**
(`results/<build>/translations/*.fasta`) to compute per-position amino-acid
entropy. The seasonal-flu build produces them; some avian builds currently
don't, so the plot falls back to mutation counts with entropy = 0. The
underlying mutation tables are still correct — only the landscape plot is
starved.

## Geographic map is empty or the slider spans one month

Same root cause as the historical section — background genomes need dates +
coordinates (see above). The slider is built from `min(date) → max(date)`;
a collapsed slider means only the query had a date.

## A panel says "No transmission data for this species"

Check whether
`results/transmission_context/<species>/strain_transmission_profile.tsv`
exists and has rows beyond its header. A header-only TSV means the upstream
warehouse query found no qualifying background genomes — usually missing
strain/date/coordinate data.

## R or Python dependency errors at launch

The dashboard auto-installs missing CRAN/Bioconductor packages on first run
but needs an internet connection to do it. For air-gapped use, pre-create the
`pgirl_dashboard` conda env.
