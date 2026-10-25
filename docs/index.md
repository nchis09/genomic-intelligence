# Pathogen Genomic Intelligence Framework — Dashboard Guide

Pathogen genomic data can provide valuable insights into emerging health
threats, outbreak dynamics, and the evolution of infectious diseases. When
interpreted alongside epidemiological and public health information, these
data can help public health professionals understand outbreak patterns,
identify changes in circulating pathogens, and inform decisions about
surveillance and response.

The **Pathogen Genomic Intelligence Framework (PGIF)** translates pathogen
genomic findings into accessible, public health-oriented intelligence. This
guide supports public health professionals, epidemiologists, and outbreak
response teams in **exploring the dashboard, understanding its
visualizations, interpreting findings in their epidemiological context, and
recognizing their implications for public health action**.

Each section provides a practical guide to the results presented in the
dashboard, explaining what the findings may indicate, what conclusions can
reasonably be drawn, the limitations of the available evidence, and what
additional information may be needed to support public health decisions.
Throughout, findings are labelled by the strength of evidence behind them —
*observed*, *estimated*, or *model-derived* — and the dashboard deliberately
avoids forecasting where the underlying data is too sparse to support it.

!!! info "How this guide is organized"
    Each dashboard page gets one chapter, written step-by-step with
    screenshots from a real pipeline run. Where a number or plot comes from,
    the chapter says so — the dashboard only shows what the analysis
    computed; it never invents results. Technical terms are explained where
    they appear, and the [genomic results glossary](interpretation/genomic-results.md)
    covers the rest.

## Where to start

- **New to PGIF?** → [Getting started](getting-started/installation.md)
  (install → run the analysis → open the dashboard).
- **Already have results?** → Jump to
  [Navigation & species selection](dashboard/navigation.md), then follow the
  Dashboard guide in order, or open the chapter matching the panel you are
  looking at.
- **Looking at a finding you don't trust?** →
  [Interpreting results](interpretation/genomic-results.md) and
  [Troubleshooting](troubleshooting.md).

## The dashboard at a glance

| Dashboard page | Question it answers |
|---|---|
| [Intelligence Overview](dashboard/intelligence-overview.md) | What did the analysis find, and how confident should we be? |
| [Biological Threat](dashboard/biological-threat.md) | Have we seen this pathogen before — and how close is it to known references? |
| [Pathogen Genomics](dashboard/pathogen-genomics.md) | Where does it sit among known viruses, and which mutations distinguish it? |
| [Transmission & Spread](dashboard/transmission-spread.md) | How did this lineage spread historically, and how intense was transmission? |
| [Geographic & Temporal Context](dashboard/geographic-temporal.md) | Where and when was it — and its relatives — circulating? |
| [Countermeasure Readiness](dashboard/countermeasure-readiness.md) | Do vaccines, therapeutics, or diagnostics exist — and how strong is the evidence? |
| [Evidence & Knowledge Gaps](dashboard/evidence-knowledge-gaps.md) | What does the literature say — and what do we still not know? |
| [Intelligence Brief](dashboard/intelligence-brief.md) | The shareable written summary for decision-makers. |

## Pipeline documentation

Technical reference for the underlying analysis pipeline:

- [Usage](usage.md) — how to run the pipeline, inputs, and parameters
- [Output](output.md) — every file the pipeline produces
- [Contributing](CONTRIBUTING.md) — development guidelines
