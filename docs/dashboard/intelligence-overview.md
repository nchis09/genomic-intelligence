# Intelligence Overview (home)

The Intelligence Overview is the landing page for the selected pathogen — a
one-screen summary designed for a quick read: *what did the analysis find,
how much confidence does it deserve, and where do I look next?*

This walkthrough uses an Ebola run (SUDV) as the example.

## Step 1 — Read the assessment badges

At the top of the card, a row of badges compresses the whole analysis into a
verdict:

<!-- screenshot: images/home/01-assessment-badges.png — badge row: MONITOR / Confidence: High / Transmission risk / Narrative: ollama -->

- **Assessment** — `MONITOR`, `INVESTIGATE`, or `ESCALATE`: how unusual the
  genome looks compared with everything we already know about. It is a
  cautious, transparent rule set — **not a model making predictions** — and
  it only escalates on *novelty signals* (the genome sitting unusually far
  from its closest known relative or its assigned family group), and only
  when the underlying data is good enough to trust.
- **Confidence** — `High` / `Medium` / `Low` / `Insufficient evidence`:
  how solid the underlying data is, driven mainly by how much of each genome
  was successfully sequenced (coverage). Low confidence never pushes the
  assessment toward ESCALATE — it simply means "don't over-read this."
- **Transmission risk** — short labels summarizing what the transmission
  analysis found for each query genome (e.g. *sustained human transmission*).
- **Narrative source** — whether the paragraph below was drafted by the local
  AI assistant (`ollama`) or assembled from a fixed template (`template`).
  See [Troubleshooting](../troubleshooting.md) if you expected AI text.

!!! note "Two separate axes"
    Assessment answers *"how worrying is this genome?"*; Confidence answers
    *"how good is the data behind that call?"* A clean, ordinary genome on
    solid data reads `MONITOR · High confidence`; a possibly-novel genome on
    thin data reads `INVESTIGATE · Medium confidence`. The rules almost never
    reach ESCALATE from identification evidence alone — that is deliberate.

## Step 2 — Scan the situation strip

Four headline numbers summarize what was submitted:

<!-- screenshot: images/home/02-situation-strip.png — Query genomes / Countries / Collection dates / Median coverage -->

- **Query genomes** — how many submitted samples were identified as this
  pathogen.
- **Countries** — where those samples were collected.
- **Collection dates** — earliest to latest sampling date.
- **Median coverage** — how completely the genomes were sequenced, on
  average; the main input to the Confidence badge.

## Step 3 — Read the narrative

The blue-bordered paragraph is the generated situation summary: a
plain-language synthesis of what the analysis found across identification,
genomics, transmission, and the published evidence.

<!-- screenshot: images/home/03-narrative.png — narrative block -->

!!! info "Where it comes from"
    Written during the analysis run and stored in
    `intelligence_brief.json`. When the local AI assistant is available it
    writes the paragraph; otherwise a fixed template fills in the same facts.
    Either way, **every number in it traces back to a result table** — it
    never contains a figure you can't find elsewhere in the dashboard.

## Step 4 — Walk the four evidence cards

Below the narrative, four clickable cards summarize each layer of evidence —
clicking one jumps straight to the full page for that question.

<!-- screenshot: images/home/04-evidence-cards.png — the four section cards stacked -->

**What is it — identification.**
The closest known reference genome and how similar the query is to it, plus
the closest matching past outbreak. Click → [Biological Threat](biological-threat.md).

**What's different about it — genomics.**
Whether the query carries more mutations than its background relatives (with
a statistical significance check), which mutations are linked to known
traits, and which proteins most distinguish it — plus two mini-plots:

- *Divergence boxplot* — how far each query sits from the root of the family
  tree (red dots) compared with the background distribution. A query far
  above the box is unusually divergent.
- *Top proteins* — the proteins whose mutations most separate the query from
  background genomes.

Click → [Pathogen Genomics](pathogen-genomics.md).

**How it's spreading — transmission evidence.**
Whether the queries belong to known genetic clusters, where the lineage most
likely originated, how many separate introductions the data suggests, and the
estimated **R at sampling** — the average number of people each case infected
at the time the samples were taken — plus a mini-map of where the query's
relatives were collected.

<!-- screenshot: images/home/05-transmission-card.png — transmission stats + mini map -->

Click → [Transmission & Spread](transmission-spread.md).

**What we know & what we can do — evidence & countermeasures.**
How much published research yielded usable findings, how many evidence areas
are covered, and the leading countermeasure product — plus an AI/auto
highlights summary and two small charts.

Click → [Evidence & Knowledge Gaps](evidence-knowledge-gaps.md).

## Step 5 — Check the caveats

The amber **Caveats** strip (when present) lists the analysis' own warnings —
e.g. thin epidemiological coverage or a small background set. Read it before
quoting any number further.

## Step 6 — The objectives strip

Beneath the cards, a muted row of tiles mirrors the sidebar's objective list —
a reminder of what else is available without competing with the assessment
above.

## Public-health relevance

This page is built for the 60-second read: a duty officer should be able to
say *"the analysis found SUDV, confidence is high, transmission evidence is
moderate, here is the summary"* — and know exactly which card to click for
the supporting evidence. Everything shown is **observed or estimated from
data**; the verdict is a transparent rule set, not a black box.

!!! warning "Limitations"
    - The assessment currently uses identification signals only; it is
      first-pass triage, not a full risk assessment.
    - "Insufficient evidence" means *missing data*, not *low threat*.
    - AI-drafted and template narratives differ in phrasing, not in facts.
