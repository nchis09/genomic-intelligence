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
       icon = "share-nodes", status = "Live",
       desc = "Estimated transmissibility and spread indicators from phylogenetic and epidemiological linkage."),
  list(id = "geographic_temporal", label = "Geographic & Temporal Context",
       icon = "earth-africa", status = "Live",
       desc = "Where and when related cases or sequences have been detected."),
  list(id = "health_impact", label = "Health Impact",
       icon = "heart-pulse",
       desc = "Case counts, severity and burden metrics associated with this pathogen or lineage."),
  list(id = "countermeasure_readiness", label = "Countermeasure Readiness",
       icon = "syringe",
       desc = "Diagnostic, therapeutic and vaccine relevance/readiness for this pathogen or lineage."),
  list(id = "evidence_knowledge_gaps", label = "Evidence & Knowledge Gaps",
       icon = "magnifying-glass",
       desc = "Supporting literature, evidence strength, and open knowledge gaps behind the current assessment."),
  list(id = "intelligence_brief", label = "Intelligence Brief",
       icon = "file-lines", status = "Live",
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
  is_live <- !is.null(objective$status) && objective$status == "Live"
  color <- if (is_live) "#4A6C8C" else "#6c757d"
  border <- if (is_live) "#4A6C8C" else "#e9ecef"
  badge_class <- if (is_live) "badge badge-primary" else "badge badge-light"
  badge_text <- if (is_live) "Live" else "Planned"
  div(
    style = paste0(
      "cursor: pointer; border: 1px solid ", border, "; border-radius: 6px;",
      "padding: 10px 12px; text-align: center; color: ", color, ";",
      "background: #fff; height: 100%;"
    ),
    onclick = sprintf(
      "Shiny.setInputValue('overview_nav_click', '%s', {priority: 'event'});",
      roadmap_tab_name(objective$id)
    ),
    icon(objective$icon),
    div(style = paste0("font-size: 0.78rem; margin-top: 4px;", if (is_live) " font-weight: 600;" else ""), objective$label),
    span(class = badge_class, style = "font-size: 0.65rem; margin-top: 2px;", badge_text)
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

# "Live" tile for Pathogen Genomics on the Intelligence Overview page.
pathogen_genomics_tile_ui <- function() {
  div(
    style = paste(
      "cursor: pointer; border: 1px solid #4A6C8C; border-radius: 6px;",
      "padding: 10px 12px; text-align: center; color: #4A6C8C;",
      "background: #fff; height: 100%;"
    ),
    onclick = "Shiny.setInputValue('overview_pg_tile_click', 'click', {priority: 'event'});",
    icon("microscope"),
    div(style = "font-size: 0.78rem; margin-top: 4px; font-weight: 600;", "Pathogen Genomics"),
    span(class = "badge badge-primary", style = "font-size: 0.65rem; margin-top: 2px;", "Live")
  )
}

# Preview/edit modal for the Intelligence Brief. The textarea is pre-filled
# with the generated brief text (brief_markdown()); edits are what get
# exported — the pipeline JSON stays untouched. "Download .md" saves the
# current textarea content; "Print / Save as PDF" opens a minimal window
# containing just the brief text and triggers the browser print dialog.
intelligence_brief_preview_modal <- function(content = NULL) {
  modalDialog(
    title = "Preview & Edit Brief",
    size = "l",
    easyClose = TRUE,
    p(
      style = "color: #6c757d; font-size: 0.85rem;",
      "Generated from the pipeline's intelligence_brief.json — edit freely",
      "before exporting; the underlying file is not modified."
    ),
    textAreaInput(
      "brief_preview_text", label = NULL, width = "100%", height = "320px",
      value = content %||% "No intelligence brief content available."
    ),
    footer = tagList(
      modalButton("Close"),
      downloadButton("brief_download", "Download .md",
                     class = "btn-outline-secondary"),
      actionButton(
        "brief_export_pdf", "Print / Save as PDF",
        icon = icon("file-pdf"), class = "btn-secondary",
        onclick = paste(
          "var t=document.getElementById('brief_preview_text').value;",
          "var w=window.open('','_blank');",
          "w.document.write('<pre style=\"font-family:monospace;white-space:pre-wrap\">'",
          "+ t.replace(/&/g,'&amp;').replace(/</g,'&lt;') +'</pre>');",
          "w.document.close(); w.print();"
        )
      )
    )
  )
}
