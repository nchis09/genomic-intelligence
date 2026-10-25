# Intelligence Overview (home)

The Intelligence Overview is the assessment-first landing page for the
selected pathogen. In one screen it answers: *what did the run find, how
confident should we be, and where do I look next?*

This walkthrough uses an Ebola run (SUDV) as the example.

## Step 1 — Read the assessment badges

At the top of the card, a row of badges compresses the whole run into a
verdict:

<!-- screenshot: images/home/01-assessment-badges.png — badge row: MONITOR / Confidence: High / Transmission risk / Narrative: ollama -->

- **Assessment** — `MONITOR`, `INVESTIGATE`, or `ESCALATE`. This is a
  deliberately conservative, rule-based heuristic, **not** a predictive
  model. It escalates only on *novelty signals* (elevated genetic distance to
  the closest reference, low sequence identity, distant placement from the
  assigned clade) and only when data confidence is adequate.
- **Confidence** — `High` / `Medium` / `Low` / `Insufficient evidence`. This
  is the *data-quality* axis: driven by genome coverage and whether the
  identification table produced confirmed rows. Low confidence never pushes
  the assessment toward ESCALATE — it just means "don't over-read this."
- **Transmission risk** — labels aggregated from the per-query transmission
  assessments (e.g. *sustained_human_transmission*).
- **Narrative source** — `ollama` when the paragraph below was written by the
  local LLM, `template` when it used the deterministic fallback (see
  [Troubleshooting](../troubleshooting.md) if you expected Ollama).

!!! note "Two separate axes"
    Assessment = *severity/novelty*; Confidence = *data quality*. A
    well-sequenced but ordinary genome reads `MONITOR · High confidence`; a
    possibly-novel genome on thin data reads `INVESTIGATE · Medium
    confidence`. The heuristic almost never reaches ESCALATE on
    identification data alone — that is intentional.

## Step 2 — Scan the situation strip

Four headline numbers summarize the query batch:

<!-- screenshot: images/home/02-situation-strip.png — Query genomes / Countries / Collection dates / Median coverage -->

- **Query genomes** — how many query samples resolved to this species.
- **Countries** — where they were collected (names up to three).
- **Collection dates** — the span between earliest and latest query sample.
- **Median coverage** — median genome coverage across queries; the main
  input to the Confidence badge.

## Step 3 — Read the narrative

The blue-bordered paragraph is the generated situation narrative: a plain
-language synthesis of identification, genomic, transmission and evidence
facts.

<!-- screenshot: images/home/03-narrative.png — narrative block -->

!!! info "Where it comes from"
    Written at pipeline time by `GENERATE_INTELLIGENCE_BRIEF` —
    `intelligence_brief.json`. With Ollama running it's an LLM-written
    paragraph grounded on the computed facts; otherwise a deterministic
    template. Either way, **every number in it traces to a pipeline table** —
    it never contains a figure you can't find elsewhere in the dashboard.

## Step 4 — Walk the four evidence cards

Below the narrative, four clickable cards summarize each evidence layer —
each one jumps to the full page for that objective.

<!-- screenshot: images/home/04-evidence-cards.png — the four section cards stacked -->

**What is it — identification.**
Closest reference and % identity, median/max genetic distance, closest
outbreak match. Click → [Biological Threat](biological-threat.md).

**What's different about it — genomics.**
Query vs background mutation burden (with p-value), phenotype-linked
mutations, top discriminating proteins — plus two inline mini-plots:

- *Divergence boxplot* — query genomes (red dots) against the background
  distribution of root-to-tip divergence. A query far above the box is
  unusually divergent.
- *Top proteins by VIP* — the proteins contributing most to separating query
  from background in the PLS-DA mutation profile.

Click → [Pathogen Genomics](pathogen-genomics.md).

**How it's spreading — transmission evidence.**
Number of genetic clusters the queries fall in (with largest cluster size),
likely origin, estimated introductions, **R at sampling** — plus an inline
map of where the query's relatives were sampled.

<!-- screenshot: images/home/05-transmission-card.png — transmission stats + mini map -->

Click → [Transmission & Spread](transmission-spread.md).

**What we know & what we can do — evidence & countermeasures.**
Papers screened that yielded claims, how many of the 11 evidence domains have
coverage, the top countermeasure product — plus an "Evidence highlights"
AI/auto summary and two small charts (claims per domain; top products).

Click → [Evidence & Knowledge Gaps](evidence-knowledge-gaps.md).

## Step 5 — Check the caveats

The amber **Caveats** strip (when present) lists the run's own data warnings —
e.g. sparse epi coverage, small background sets. Read it before quoting any
number downstream.

## Step 6 — The objectives strip

Beneath the cards, a de-emphasized row of tiles mirrors the sidebar's
objective list — a reminder of what else is available without competing with
the assessment above.

## Public-health relevance

This page is built for the 60-second read: a duty officer should be able to
say *"the run found SUDV, confidence is high, transmission evidence is
moderate, here is the narrative"* and know exactly which card to click for
the supporting evidence. Everything on it is **observed or estimated from
data** — the verdict heuristic is transparent rules, not a black box.

!!! warning "Limitations"
    - The assessment only uses Pathogen Identification signals today; it is a
      first-pass triage, not a risk assessment.
    - "Insufficient evidence" means *missing*, not *safe*.
    - Template vs Ollama narratives differ in phrasing, not in facts.
