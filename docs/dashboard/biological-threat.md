# Biological Threat

*Have we seen this virus before — and how confident is the identification?*

This page compares your submitted genomes against known references. Pick up
to five samples and it builds a side-by-side profile: how completely each
genome was sequenced, how closely it matches its nearest known relative, and
how far it sits from its assigned family group.

## Step 1 — Pick samples to compare

The bar at the top lists every sample the pipeline confirmed for this
pathogen.

<!-- screenshot: images/biological-threat/01-sample-picker.png — selectize bar + badge "3 samples selected" -->

- Type or click to add samples — **up to 5** keeps the charts readable.
- The badge counts your selection; **Reset selection** clears it.
- Nothing appears until you pick at least one sample — everything on this
  page is a comparison view.

## Step 2 — Read the summary strip

Once samples are selected, three cards summarize the group:

<!-- screenshot: images/biological-threat/02-summary-strip.png — pathogen / samples compared / closest outbreak cards -->

- **Pathogen** — the species all selected samples were identified as.
- **Samples compared** — your selection out of the total confirmed.
- **Closest outbreak(s)** — the past outbreak(s) these genomes resemble most
  closely, or *N/A*.

!!! info "Where it comes from"
    The pipeline's identification results
    (`identification_summary.tsv`). "Closest" here means *closest in genetic
    similarity* — the genomes look alike — not proof of an epidemiological
    link to that outbreak.

## Step 3 — Compare the genomic metrics

With **two or more** samples selected, a grouped bar chart compares six
identification metrics — one coloured bar per sample, grouped by metric:

<!-- screenshot: images/biological-threat/03-metric-bars.png — grouped bars, six metric sections on x-axis -->

| Metric | What it means | Direction |
|---|---|---|
| Coverage (%) | How much of the genome was successfully sequenced | Higher is better |
| % Identity to closest reference | How similar the genome is to its nearest known match | Higher is better |
| Genetic distance to closest reference | How many mutations separate it from that match | **Lower is better** |
| QC overall score | The classifier's overall quality flag for the sequence | **Lower is better** |
| Distance to assigned clade | How far it sits from the centre of its assigned family group | **Lower is better** |
| MSA Identity (%) | How similar the query is to the rest of the alignment, on average | Higher is better |

!!! warning "A taller bar is not always better"
    The six metrics live on different scales, so each bar is drawn as a
    percentage of the highest value in its section. Read each group against
    its direction label: a query standing out on a *lower-is-better* metric
    (genetic distance, clade distance, QC score) is the anomaly worth
    attention — it suggests a novel or unusual genome, or a quality problem.

## Step 4 — The identification results table

Below the charts, the full table lists every confirmed sample: coverage, QC
score, the closest reference genome (with where and when it was collected),
% identity, genetic and clade distances, overall alignment similarity, and
closest matching outbreak.

<!-- screenshot: images/biological-threat/04-results-table.png — DT table with checkbox column -->

- The **checkbox column is synced with the picker** — ticking rows here is
  the same as selecting samples at the top.
- Sort and filter any column; "Not Provided" in the location field means the
  reference record simply didn't carry one.

!!! note "Only confirmed identifications appear"
    The table is restricted to samples that passed the identification
    quality check. Rejected samples aren't missing — they're in the last
    card on this page.

## Step 5 — Genomic variation profile (radar)

The next card has its own picker and draws a **radar chart** of how each
selected genome differs from the reference: counts of substitutions,
deletions, insertions, frameshifts, missing bases, ambiguous characters, and
amino-acid changes.

<!-- screenshot: images/biological-threat/05-radar.png — radar chart, several axes -->

- Each arm is one kind of variation; each coloured shape is one sample.
- Use it to compare *how* genomes differ, not just *how much*: two samples
  can have identical identification stats but very different profiles — for
  example, one dominated by genuine substitutions versus another dominated
  by missing data (a sequencing-quality difference, not a biological one).
- The table below covers **all** query samples, not only confirmed ones.

## Step 6 — Unresolved samples

The last card (collapsed by default) lists samples the classifier **could not
assign with confidence** — ambiguous or rejected calls, with their coverage
and QC status.

<!-- screenshot: images/biological-threat/06-unresolved.png — unresolved samples table -->

Treat this as a data-quality queue: poor sequencing or mixed samples explain
most rejections, and they're candidates for re-sequencing rather than
signals of something new.

## What to watch for

Three reads matter most on this page:

1. **The assignment and outbreak match** — does the genome map onto a known
   threat we already understand?
2. **Coverage and QC** — can you trust the call? Low coverage degrades every
   distance metric; read novelty cautiously on poorly sequenced samples.
3. **The distance columns** — a well-sequenced genome sitting far from every
   reference is a possible novel lineage and a candidate for confirmatory
   work.

!!! warning "Limitations"
    - "Closest outbreak" is genetic similarity, not epidemiological linkage —
      an old isolate can match a new query purely by relatedness.
    - Identification is against a *reference set*: a genuinely novel pathogen
      far from all references may appear as unresolved rather than as a new
      finding.
    - This page says what the genome *is*, not what it will *do* — spread and
      impact are covered by the Transmission and Geographic pages.
