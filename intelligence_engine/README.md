# Intelligence Engine (Tier 3)

Code that turns a new genome + the `reference_library` (PGIRL) into an Intelligence Brief.
The engine performs **lookups and rule application only** — no AI, no biological inference.

## Pipeline stages (folders)

1. **`bioinformatics/`**
   Take a new consensus genome (+ metadata). Sub-divided by task:
   - `qc/` — input validation, quality metrics, and metadata checks.
   - `taxonomic_classification/` — species / lineage / clade assignment (Nextclade, BLAST, Kraken, etc.).
   - `templates/` — staged bioinformatics output templates.
   - `validation/` — metadata schema definitions.
   May wrap existing tools (Nextclade, pangolin-style, BLAST, Kraken) — we do **not** reinvent
   genomic analysis.

2. **`data_engine/`**
   Query the `reference_library` and external sources:
   - `sql_querying/` — curated database queries (PostgreSQL).
     - `bioinformatics_query.py` — reads per-sample `bio_output.json` and queries the local DB for variants, lineages, phenotype associations, and surveillance context.
   - `online_querying/` — DuckDuckGo text search for epidemiological reports.
     Returned pages are ranked dynamically by source credibility (WHO, CDC, UN,
     national Ministries of Health / government health agencies, Africa CDC,
     ECDC, GOV.UK, ReliefWeb, etc.). Requires `pip install ddgs` for best results
     (falls back to a lightweight HTML scraper).
   - `llm_querying/` — structured LLM-based summarisation of fetched evidence.
   - Has this species/lineage been seen in this location before? (`geography`, `outbreaks`)
   - Is the lineage known or novel? (`molecular_epidemiology`)
   - Closest historical outbreak? (`outbreaks`, `reference_genomes`)
   - Are called mutations in the library? What associations + evidence? (`mutations` → `literature`)

3. **`evidence_integration/`**
   Harmonize molecular + epidemiological evidence into unified `EvidenceObject`s and run
   cross-evidence statistical analyses (mutation co-occurrence, lineage-phenotype association
   via Fisher's exact/Chi-square, temporal trend via Poisson GLM, geographic distribution,
   mutation persistence, intervention association, confidence scoring). Produces quantitative
   evidence summaries only -- no risk tiers or public-health conclusions.
   - `engine.py` — core analyzer suite (genomic significance, molecular epidemiology, phylogeography, etc.).
   - `harmonization.py` — builds `EvidenceObject`s linking variants/lineages/phenotypes/outbreaks/epi context.
   - `cross_evidence.py` — statistical cross-evidence analyzers.
   - `visualization.py` — evidence-network (networkx) and geo/temporal charts (matplotlib).
   - `analyzers/` — extended decision-oriented analyzers.
   - `pipeline/` — single-sample and multi-sample batch pipelines; writes `evidence_package.json`.
   - `tree/` — phylogenetic tree input and synthetic tree generation.
   - `figures/` — R scripts that generate decision-oriented figures for the report.
   - `examples/` — example manifests for batch mode.

4. **`genomic_intelligence/`**
   Acting as an expert genomic intelligence analyst, contextualize and synthesize the evidence
   produced by `evidence_integration/` into a coherent, evidence-grounded narrative via an LLM
   (`data_engine/llm_querying/llm_client.py`), reasoning ONLY over the provided evidence --
   never the model's own trained knowledge -- and citing the specific finding/metric/source
   behind every statement. Does not run new statistics or make public-health recommendations.
   - `context_builder.py` — curates intelligence_object.json + evidence_package.json into the
     grounding context passed to the LLM (excludes risk-tier/recommendation content).
   - `synthesize.py` — CLI entry point; writes `genomic_intelligence_assessment.md` and the
     exact `context_used.json` given to the model, for full traceability.

## Contract
- The engine may output "**Not enough information**" — a valid, expected result.
- Every claim in the brief traces to a curated `reference_library` entry (or a flagged
  external source pending curation).
