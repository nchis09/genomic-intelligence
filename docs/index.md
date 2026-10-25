# Genomic Intelligence Framework — Dashboard Guide

The **Genomic Intelligence Framework (PGIF)** is a Nextflow pipeline that turns
consensus pathogen genomes into an epidemic-intelligence picture: what the
pathogen is, how it differs from known references, how it is spreading, where
it has been, and what we can do about it.

This site documents the **results-visualization dashboard** — the R Shiny app
that lets you browse everything the pipeline produced. It is written as a
guided tour: for each dashboard page we show the visualization, explain what
it means, note its limitations, and describe the public-health question it
helps answer.

!!! info "How this guide is organized"
    Each dashboard page gets one chapter, written step-by-step with
    screenshots from a real pipeline run. Where a number or plot comes from,
    the chapter says so — the dashboard only shows what the pipeline computed;
    it never invents results.

## Where to start

- **New to PGIF?** → [Getting started](getting-started/installation.md)
  (install → run the pipeline → launch the dashboard in ~10 minutes).
- **Already have results?** → Jump to
  [Navigation & species selection](dashboard/navigation.md), then follow the
  Dashboard guide in order, or open the page matching the panel you are
  looking at.
- **Looking at a number you don't trust?** →
  [Interpreting results](interpretation/genomic-results.md) and
  [Troubleshooting](troubleshooting.md).

## The dashboard at a glance

| Dashboard page | Question it answers |
|---|---|
| [Intelligence Overview](dashboard/intelligence-overview.md) | What did the run find, and how confident should we be? |
| [Biological Threat](dashboard/biological-threat.md) | What is this genome — species, clade, how close to known references? |
| [Pathogen Genomics](dashboard/pathogen-genomics.md) | Where does it sit on the tree, and which mutations distinguish it? |
| [Transmission & Spread](dashboard/transmission-spread.md) | How did it spread historically, and how intense was transmission? |
| [Geographic & Temporal Context](dashboard/geographic-temporal.md) | Where and when was it (and its relatives) circulating? |
| [Countermeasure Readiness](dashboard/countermeasure-readiness.md) | Do vaccines/therapeutics exist and how well do they work? |
| [Evidence & Knowledge Gaps](dashboard/evidence-knowledge-gaps.md) | What does the literature say — and what do we still not know? |
| [Intelligence Brief](dashboard/intelligence-brief.md) | The shareable written summary for decision-makers. |

## Pipeline documentation

This site also carries the pipeline's own reference docs:

- [Usage](usage.md) — how to run the pipeline, inputs, and parameters
- [Output](output.md) — every file the pipeline produces
- [Contributing](CONTRIBUTING.md) — development guidelines

## A note on evidence

PGIF labels every visualization by the kind of evidence behind it —
**observed history** (measured data), **estimated** (e.g. R(t) fitted to
reported cases), or **model-derived**. Figures that forecast the future were
deliberately removed: the pipeline does not collect enough covariates to
project cases responsibly, so the dashboard reports what the data supports —
nothing more.
