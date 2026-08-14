# Genomic Intelligence Overview: the Home page. A single card whose header
# reads "Intelligence Brief — <species>" with the species dropdown merged
# inline into the header bar itself (no separate "Species" label/box), plus
# two CTAs and a de-emphasized roadmap strip for the other 9 not-yet-wired
# objectives. Relies on ROADMAP_OBJECTIVES / roadmap_tile_ui()
# (modules/placeholder.R).

# Reads the GIF summary text already written for MultiQC's intro_text so the
# copy only has to be maintained in one place.
read_gif_intro_text <- function(config_path = "../assets/multiqc_config.yml") {
  default_text <- paste(
    "The Genomic Intelligence Framework (GIF) is a public-health-support",
    "that transforms pathogen genomic data into actionable public-health",
    "intelligence. It systematically assesses pathogen identity, biological",
    "threat, and transmission to inform potential health impact."
  )
  if (!file.exists(config_path)) return(default_text)
  cfg <- tryCatch(yaml::read_yaml(config_path), error = function(e) NULL)
  if (is.null(cfg) || is.null(cfg$intro_text)) return(default_text)
  trimws(cfg$intro_text)
}

# Static shell for the Overview tab: the card header/body are filled in
# server-side via renderUI (uiOutput("overview_header") /
# uiOutput("overview_body")), since they depend on the discovered species
# list and the selected species. CTAs and the roadmap strip are static --
# they don't depend on species/outdir.
overview_ui <- function() {
  tagList(
    h2("Genomic Intelligence", style = "margin-bottom: 12px; font-weight: bold;"),
    bs4Dash::bs4Card(
      width = 12,
      status = "secondary",
      solidHeader = TRUE,
      title = uiOutput("overview_header", inline = TRUE),
      uiOutput("overview_body")
    ),
    fluidRow(
      style = "margin-top: 8px;",
      column(
        width = 6,
        actionButton("overview_explore_evidence", "Explore Evidence",
                     icon = icon("magnifying-glass"), class = "btn-primary btn-block")
      ),
      column(
        width = 6,
        actionButton("overview_generate_brief", "Generate Genomic Intelligence Brief",
                     icon = icon("file-lines"), class = "btn-outline-secondary btn-block")
      )
    ),
    h5("Roadmap", style = "margin-top: 28px; color: #6c757d;"),
    p(
      style = "color: #adb5bd; font-size: 0.8rem; margin-top: -6px;",
      "These intelligence objectives don't have a wired data source yet."
    ),
    fluidRow(
      lapply(ROADMAP_OBJECTIVES, function(obj) {
        column(width = 2, style = "margin-bottom: 10px;", roadmap_tile_ui(obj))
      })
    )
  )
}

# Card header content: "Intelligence Brief" plus an inline species dropdown,
# borderless/blended into the dark header bar so it reads as one continuous
# line ("Intelligence Brief — [bdbv ▾]"). Falls back to plain text when
# there's 0 or 1 species (nothing to choose between).
overview_header_ui <- function(species) {
  if (length(species) == 0) return(span("Intelligence Brief"))
  if (length(species) == 1) {
    return(span(paste0("Intelligence Brief \u2014 ", toupper(species[1]))))
  }
  div(
    class = "header-inline-select",
    span("Intelligence Brief \u2014"),
    selectInput("overview_species", label = NULL, choices = species,
                selected = species[1], selectize = FALSE)
  )
}

# Card body: currently a neutral placeholder -- no risk verdict or evidence
# facts are surfaced here yet.
overview_body_ui <- function(has_species) {
  if (!has_species) {
    return(p("No Biological Threat data found under the configured --outdir yet."))
  }
  tags$ul(tags$li("In planning."))
}
