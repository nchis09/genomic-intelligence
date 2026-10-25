# Limitations

The dashboard is built to show *what the data supports* — and to say so when
the support is thin. Keep these caveats in mind when quoting results.

## Evidence labels

Panels are labelled **observed history**, **estimated** (model fitted to
observed data, e.g. R(t), growth-rate fits), or **tree-derived**. Nothing on
the dashboard forecasts future cases — that capability was removed because
the pipeline does not collect enough covariates to project responsibly.

## Structural limitations

- **Background sampling bias** — the tree's "history" is whatever genomes the
  builds downloaded. Unsampled transmission is invisible to every panel.
- **Assessment heuristic is triage, not risk assessment** — MONITOR /
  INVESTIGATE / ESCALATE uses only identification signals today (see
  [Intelligence Overview](../dashboard/intelligence-overview.md)).
- **Country-level coordinates** — background genomes map at country
  centroids unless finer locations exist; "origin" means *earliest sampled*,
  not true geographic origin.
- **Segmented-genome nuance (influenza)** — the same isolate can carry
  different clade labels per segment (reassortment); HA labels take priority
  in the warehouse.
- **LLM narratives** — Ollama-written text is grounded on the computed facts
  but can phrase them loosely; the numbers to trust are in the tables.
- **Epi coverage varies** — FluNet/HDX/EMPRES-i reporting is uneven across
  countries and years; absence of a signal is not evidence of absence.
