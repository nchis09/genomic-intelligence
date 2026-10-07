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

# Tile order for the objectives grid: Intelligence Brief first (it's the
# product's front page), then the remaining objectives in declared order.
# ROADMAP_OBJECTIVES itself is untouched — sidebar menu slicing relies on it.
brief_objective <- function() {
  Filter(function(o) identical(o$id, "intelligence_brief"), ROADMAP_OBJECTIVES)[[1]]
}

ordered_objectives <- function() {
  Filter(function(o) !identical(o$id, "intelligence_brief"), ROADMAP_OBJECTIVES)
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
    h5("Intelligence Objectives", style = "margin-top: 28px; color: #6c757d; font-weight: bold;"),
    fluidRow(
      column(width = 2, style = "margin-bottom: 10px;", roadmap_tile_ui(brief_objective())),
      column(width = 2, style = "margin-bottom: 10px;", live_tile_ui()),
      column(width = 2, style = "margin-bottom: 10px;", pathogen_genomics_tile_ui()),
      lapply(ordered_objectives(), function(obj) {
        column(width = 2, style = "margin-bottom: 10px;", roadmap_tile_ui(obj))
      })
    )
  )
}

# Card header content: "Intelligence Brief" plus an inline species dropdown,
# borderless/blended into the dark header bar so it reads as one continuous
# line ("Intelligence Brief — [bdbv ▾]"). Falls back to plain text when
# there's 0 or 1 species (nothing to choose between). `family`/`n_families`
# add a family badge + a back-to-cards link for multi-pathogen runs.
overview_header_ui <- function(species, family = NULL, n_families = 1) {
  switcher <- if (!is.null(family)) {
    span(
      style = "font-size: 0.75rem; font-weight: 400; margin-left: 10px;",
      span(class = "badge badge-info", style = "margin-right: 4px;",
           family_display(family)),
      actionLink("fam_reset", "change", style = "color: #cfe2f3;")
    )
  } else NULL

  if (length(species) == 0) return(span("Intelligence Brief", switcher))
  if (length(species) == 1) {
    return(span(paste0("Intelligence Brief \u2014 ", toupper(species[1])), switcher))
  }
  div(
    class = "header-inline-select",
    span("Intelligence Brief \u2014"),
    selectInput("overview_species", label = NULL, choices = species,
                selected = species[1], selectize = FALSE),
    switcher
  )
}

# Card body: the Intelligence Brief front page — situation strip, verdict
# badges, bottom-line narrative, and one compact block per intelligence
# objective with deep links into the detail tabs. Rendered from the
# pre-generated intelligence_brief.json (GENERATE_INTELLIGENCE_BRIEF).
overview_body_ui <- function(species, outdir = NULL) {
  if (is.null(species)) {
    return(p("No Biological Threat data found in the pipeline output yet."))
  }
  brief_overview_ui(brief_read(outdir, species), species, outdir)
}

# ---------------------------------------------------------------------------
# Pathogen-first landing: when the run detected more than one pathogen
# family, the home page opens on family cards (Ebola / Influenza / other)
# instead of jumping straight into one species' brief — picking a card gates
# the sidebar to that family's tailored views.
# `families` is a data.frame(species, pathogen, family); family_display()
# gives the user-facing label.
# ---------------------------------------------------------------------------
family_display <- function(family) {
  switch(family,
    ebola = "Ebola",
    influenza = "Influenza",
    tools::toTitleCase(gsub("_", " ", family))
  )
}

family_icon <- function(family) {
  switch(family,
    ebola = "virus",
    influenza = "viruses",
    "dna"
  )
}

# Short per-family blurb shown on the detection cards.
family_blurb <- function(family, n_species) {
  sp_txt <- if (n_species == 1) "1 species group" else paste(n_species, "species groups")
  switch(family,
    ebola = paste0(
      "Non-segmented filovirus genome. Views: biological threat, phylotree, ",
      "mutation landscape, transmission context. Detected: ", sp_txt, "."),
    influenza = paste0(
      "Segmented genome (up to 8 segments). Segment-aware views: per-segment trees, ",
      "mutation landscape, genotype constellations & reassortment, phenotype markers, ",
      "genome integrity and group diversity. Detected: ", sp_txt, "."),
    paste("Detected:", sp_txt, "."))
}

# Landing page body: intro copy + one card per detected pathogen family
# (empty-state card when nothing was detected). Selection is always
# explicit — the family card click is the entry point to that family's
# tailored sidebar.
pathogen_select_ui <- function(families) {
  tagList(
    h2("Genomic Intelligence",
       style = "margin-bottom: 8px; font-weight: bold;"),
    p(class = "home-intro",
      "This dashboard is the readout of the Genomic Intelligence Framework — a ",
      "pipeline that takes raw pathogen genomes and turns them into intelligence ",
      "you can act on."
    ),
    p(class = "home-intro-sub",
      "The pipeline screened your sequences against the full reference panel and ",
      "grouped them by detected pathogen. Below are the pathogen families found ",
      "in this run — each opens a dashboard tailored to that pathogen's biology."
    ),

    if (is.null(families) || nrow(families) == 0) {
      bs4Dash::bs4Card(
        width = 12, status = "secondary", solidHeader = TRUE,
        title = "No pathogen detected",
        div(
          class = "fam-empty-state",
          icon("magnifying-glass", class = "fa-3x",
               style = "margin-bottom: 14px;"),
          h4("No pathogen was detected in this output",
             style = "margin-bottom: 8px;"),
          p(style = "max-width: 560px; margin: 0 auto;",
            "No species groups were found under this results directory. ",
            "Check that the pipeline ran to completion (classification stage) ",
            "and that the dashboard is pointed at the right --outdir.")
        )
      )
    } else {
      fams <- unique(families$family)
      fluidRow(
        lapply(fams, function(fam) {
          fam_sp <- families$species[families$family == fam]
          column(
            width = 5,
            style = "margin-bottom: 14px;",
            bs4Dash::bs4Card(
              title = NULL,
              width = 12,
              status = if (fam == "influenza") "info" else "primary",
              solidHeader = FALSE,
              div(
                style = "padding: 18px 18px 14px 18px;",
                div(
                  style = "display: flex; align-items: center; gap: 14px; margin-bottom: 12px;",
                  icon(family_icon(fam), class = "fa-3x",
                       style = "color: #4A6C8C;"),
                  div(
                    h4(family_display(fam),
                       style = "margin: 0; font-weight: 700;"),
                    div(style = "margin-top: 4px;",
                      lapply(fam_sp, function(sp) {
                        span(class = "badge badge-secondary",
                             style = "margin-right: 4px; font-size: 0.8rem;",
                             toupper(sp))
                      })
                    )
                  )
                ),
                p(class = "fam-card-blurb",
                  family_blurb(fam, length(fam_sp))),
                actionButton(
                  inputId = paste0("fam_select_", fam),
                  label = paste("Open", family_display(fam), "dashboard"),
                  icon = icon("arrow-right"),
                  class = "btn-primary btn-block"
                )
              )
            )
          )
        })
      )
    }
  )
}

# Pathogen Genomics landing page: choose a species and whether to view the
# phylotree or the mutation profile.
pg_home_ui <- function(species) {
  if (length(species) == 0) {
    return(bs4Dash::bs4Card(
      title = "Pathogen Genomics", width = 12, status = "secondary",
      div(
        style = "text-align: center; padding: 40px 0; color: #888;",
        icon("microscope", class = "fa-3x"),
        h4("No species found", style = "margin-top: 16px;"),
        p("No species were found in the pipeline output. Run the pipeline, or check that results/ exists.")
      )
    ))
  }

  bs4Dash::bs4Card(
    title = "Pathogen Genomics",
    width = 12,
    status = "secondary",
    lapply(species, function(sp) {
      div(
        style = "display: inline-block; width: 260px; margin: 12px; vertical-align: top;",
        bs4Dash::bs4Card(
          title = toupper(sp),
          width = 12,
          status = "primary",
          solidHeader = FALSE,
          actionButton(
            inputId = paste0("pg_home_phylotree_", sp),
            label = "Phylotree",
            icon = icon("share-nodes"),
            class = "btn-primary btn-block",
            style = "margin-bottom: 8px;"
          ),
          actionButton(
            inputId = paste0("pg_home_mutation_", sp),
            label = "Mutation",
            icon = icon("chart-bar"),
            class = "btn-outline-primary btn-block"
          )
        )
      )
    })
  )
}
