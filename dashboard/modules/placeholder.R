# Roadmap module: the 8 intelligence objectives/pages that don't have a wired
# data source yet (everything except Biological Threat, which now carries the
# former Pathogen Identification content). Deliberately neutral/quiet styling
# -- no "coming soon" badge clutter -- since these are an intentional, visible
# roadmap rather than broken features.
#
# NOTE: Compare (nearest isolates / lineages / previous outbreaks) is a
# distinct future capability and is NOT included here -- it has no nav item
# or roadmap tile in this pass.

ROADMAP_OBJECTIVES <- list(
  list(id = "transmission_spread", label = "Transmission & Spread",
       icon = "share-nodes",
       desc = "Estimated transmissibility and spread indicators from phylogenetic and epidemiological linkage."),
  list(id = "geographic_temporal", label = "Geographic & Temporal Context",
       icon = "earth-africa",
       desc = "Where and when related cases or sequences have been detected."),
  list(id = "health_impact", label = "Health Impact",
       icon = "heart-pulse",
       desc = "Case counts, severity and burden metrics associated with this pathogen or lineage."),
  list(id = "populations_at_risk", label = "Populations & Settings at Risk",
       icon = "users",
       desc = "Demographic and setting-level risk factors relevant to further spread."),
  list(id = "countermeasure_readiness", label = "Countermeasure Readiness",
       icon = "syringe",
       desc = "Diagnostic, therapeutic and vaccine relevance/readiness for this pathogen or lineage."),
  list(id = "risk_action", label = "Risk & Action",
       icon = "triangle-exclamation",
       desc = "Composite risk rating and recommended public-health actions."),
  list(id = "evidence_knowledge_gaps", label = "Evidence & Knowledge Gaps",
       icon = "magnifying-glass",
       desc = "Supporting literature, evidence strength, and open knowledge gaps behind the current assessment."),
  list(id = "intelligence_brief", label = "Intelligence Brief",
       icon = "file-lines",
       desc = "A concise, standardized brief communicating the assessment, confidence, key findings and priority actions to decision-makers.")
)

roadmap_tab_name <- function(id) paste0("roadmap_", id)

# Full-page roadmap card for an objective's own tab.
roadmap_ui <- function(objective) {
  bs4Dash::bs4Card(
    title = objective$label,
    width = 12,
    status = "secondary",
    div(
      style = "text-align: center; padding: 48px 0; color: #6c757d;",
      icon(objective$icon, class = "fa-3x"),
      h4(objective$label, style = "margin-top: 16px;"),
      p(style = "max-width: 560px; margin: 8px auto 16px;", objective$desc),
      span(class = "badge badge-light", style = "border: 1px solid #ced4da; padding: 4px 10px;", "Planned")
    )
  )
}

# Small, quiet tile used in the Intelligence Overview's de-emphasized
# roadmap strip. Clicking it jumps the sidebar to the objective's own tab.
roadmap_tile_ui <- function(objective) {
  div(
    style = paste(
      "cursor: pointer; border: 1px solid #e9ecef; border-radius: 6px;",
      "padding: 10px 12px; text-align: center; color: #6c757d;",
      "background: #fff; height: 100%;"
    ),
    onclick = sprintf(
      "Shiny.setInputValue('overview_nav_click', '%s', {priority: 'event'});",
      roadmap_tab_name(objective$id)
    ),
    icon(objective$icon),
    div(style = "font-size: 0.78rem; margin-top: 4px;", objective$label),
    span(class = "badge badge-light", style = "font-size: 0.65rem; margin-top: 2px;", "Planned")
  )
}

# Intelligence Brief's own roadmap tab: same "Planned" framing as the other
# objectives, but with an added UI-only scaffold for the intended "preview &
# edit before export" flow -- no real brief content or PDF rendering yet
# (there's no data source behind it), just the interaction pattern in place
# so it's easy to wire up later. The button opens a modalDialog (see
# app.R's observeEvent(input$brief_preview_open, ...)).
intelligence_brief_roadmap_ui <- function(objective) {
  tagList(
    roadmap_ui(objective),
    bs4Dash::bs4Card(
      title = "Export", width = 12, status = "secondary",
      p(
        style = "color: #6c757d;",
        "Once the Intelligence Brief has real content, exporting will open an",
        "editable preview so the brief can be reviewed and adjusted before it's",
        "printed or shared. The button below previews that interaction pattern."
      ),
      actionButton("brief_preview_open", "Preview & Edit Brief", icon = icon("file-lines"))
    )
  )
}

# The one "Live" tile in the Intelligence Overview's roadmap strip --
# Biological Threat, the only objective with a wired data source today.
# Styled distinctly (institutional blue, "Live" badge) from the "Planned"
# tiles. Clicking it jumps to the current species' Biological Threat tab
# (handled server-side via observeEvent(input$overview_live_tile_click, ...)
# in app.R, mirroring the "Explore Evidence" button).
live_tile_ui <- function() {
  div(
    style = paste(
      "cursor: pointer; border: 1px solid #4A6C8C; border-radius: 6px;",
      "padding: 10px 12px; text-align: center; color: #4A6C8C;",
      "background: #fff; height: 100%;"
    ),
    onclick = "Shiny.setInputValue('overview_live_tile_click', 'click', {priority: 'event'});",
    icon("dna"),
    div(style = "font-size: 0.78rem; margin-top: 4px; font-weight: 600;", "Biological Threat"),
    span(class = "badge badge-primary", style = "font-size: 0.65rem; margin-top: 2px;", "Live")
  )
}

# Placeholder skeleton shown inside the preview/edit modal -- editable, but
# not backed by real data yet.
intelligence_brief_preview_modal <- function() {
  modalDialog(
    title = "Preview & Edit Brief",
    size = "l",
    easyClose = TRUE,
    p(
      style = "color: #6c757d; font-size: 0.85rem;",
      "This is a placeholder preview \u2014 content will be generated from the",
      "assessment once the Intelligence Brief objective is wired in."
    ),
    textAreaInput(
      "brief_preview_text", label = NULL, width = "100%", height = "320px",
      value = paste(
        "SPECIES / ASSESSMENT\n[to be filled in]\n\n",
        "CONFIDENCE\n[to be filled in]\n\n",
        "KEY FINDINGS\n[to be filled in]\n\n",
        "RECOMMENDED ACTIONS\n[to be filled in]",
        sep = ""
      )
    ),
    footer = tagList(
      modalButton("Close"),
      actionButton(
        "brief_export_pdf", "Print / Save as PDF",
        icon = icon("file-pdf"), class = "btn-secondary disabled",
        title = "Available once Intelligence Brief content is wired in"
      )
    )
  )
}
