# Navigation & species selection

How the dashboard is organized and how to move between pathogens. Every other
chapter assumes you can do the three things on this page.

## The landing screen: pathogen cards

When the app loads, the home page is a **pathogen-family selector** — one
card per family the run detected, so you only ever see dashboards for
pathogens that were actually found in your genomes.

![Pathogen-family landing page: Ebola and Influenza cards](../images/navigation/02-family-cards.png)

Each card shows:

- the **family name** and a one-line summary of its biology (for example,
  a segmented genome means the dashboard will offer per-segment views);
- the **species detected** in your run as badges — here *BDBV, EBOV, SUDV*
  under Ebola and *H3N2, H5NX* under Influenza;
- an **"Open … dashboard"** button that takes you into that family's
  Intelligence Overview.

!!! note "Why families, not just a species list"
    Different pathogens need different views — a segmented influenza genome
    gets per-segment analyses that a single-genome filovirus doesn't. You
    pick the family first; the dashboard then shows only what applies to
    that biology.

The button row also includes a reset — the **"change pathogen"** flow
returns you to these cards at any time.

## The global pathogen switcher

Once a family is selected, a dropdown is pinned to the **centre of the top
navbar**, next to a family badge:

<!-- screenshot: images/navigation/03-navbar-switcher.png — centred dropdown showing SUDV + EBOLA badge -->

- It lists every detected species grouped by family, so you can jump between,
  say, `SUDV` and `H3N2` without going back to the cards.
- The pill-shaped badge beside it names the family of the current species.
- It stays visible on **every** tab — the current pathogen is never ambiguous.

## The sidebar

The sidebar is organized around the eight intelligence objectives. Sections
expand per species, and only objectives with real data get full pages — the
rest show a short "Planned" roadmap note.

<!-- screenshot: images/navigation/04-sidebar.png — expanded sidebar, one species under each objective -->

| Sidebar group | Contains |
|---|---|
| **Intelligence Overview** | The assessment-first home page for the selected species |
| **INTELLIGENCE OBJECTIVES** | Biological Threat → per-species pages; Pathogen Genomics → *Phylotree* and *Mutation* submenus |
| **INFLUENZA ANALYSIS** | *(influenza runs only)* Segment-signatures pages |
| Objectives continued | Transmission & Spread, Geographic & Temporal Context, Countermeasure Readiness — one page per species |
| **EVIDENCE & REPORTING** | Evidence & Knowledge Gaps, Intelligence Brief |

!!! tip
    Expanding one objective auto-collapses the others — only one tree stays
    open at a time, so the sidebar never scrolls far.

## The --outdir field

Top of the sidebar: the **Pipeline --outdir** text box. Change it and the
whole app reloads against a different results directory — useful for
comparing two runs. See [Launching the
dashboard](../getting-started/launching-the-dashboard.md#point-it-at-your-results).

## How pages relate

Most visualizations deep-link into each other:

- Clicking a section card on the **Intelligence Overview** jumps straight to
  the corresponding objective page (e.g. *"What's different about it"* →
  Pathogen Genomics).
- The Pathogen Genomics landing page offers per-species **Phylotree** and
  **Mutation** shortcuts.
- The **Intelligence Brief** aggregates facts from every section into one
  exportable document.

The flow mirrors an analyst's read: *what is it → what's different → how is
it spreading → where/when → what can we do → what's the brief*.
