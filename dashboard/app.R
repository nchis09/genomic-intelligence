#!/usr/bin/env Rscript

# Genomic Intelligence Framework -- Dashboard (v3, Intelligence Overview)
#
# bs4Dash-based dashboard: an assessment-first "Intelligence Overview" home
# page (MONITOR/INVESTIGATE/ESCALATE + confidence + Key Intelligence Signals
# for Biological Threat, the only objective with real pipeline data today
# -- it carries the former "Pathogen Identification" content), a sidebar
# organized around the 8 intelligence objectives plus Evidence & Knowledge
# Gaps and Intelligence Brief, and roadmap placeholder pages for every
# objective that doesn't have a wired data source yet.
# Reads the flat TSV outputs already produced by bin/pathogen_identification.R
# (results/pathogen_identification/<species>/species_identification/*.tsv).
# Not wired into the Nextflow pipeline -- launch manually:
#
#   Rscript -e "shiny::runApp('dashboard')"
#
# See dashboard/README.md for details.

suppressPackageStartupMessages({
  library(shiny)
  library(bs4Dash)
  library(fresh)
  library(DT)
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(plotly)
  library(yaml)
})

source("modules/pathogen_identification.R")
source("modules/assessment.R")
source("modules/placeholder.R")
source("modules/home.R")

# Serve repo-root assets (GIF logo, institution logos) without copying them
# into dashboard/www.
addResourcePath("gif_assets", normalizePath(file.path("..", "assets")))

list_species <- function(outdir) {
  base <- file.path(outdir, "pathogen_identification")
  if (!dir.exists(base)) return(character())
  basename(list.dirs(base, recursive = FALSE))
}

# Per-logo pixel heights, tuned per source image (WHO_Hub.png and UVRI.jpg
# have a lot of built-in whitespace/lettering next to the mark, so they need
# a taller box than the others to read at the same visual size).
INSTITUTION_LOGOS <- c(
  "GOARN.png"   = 20,
  "MRC.jpeg"    = 20,
  "RKI.png"     = 20,
  "UVRI.jpg"    = 28,
  "WHO_Hub.png" = 40
)

gif_theme <- create_theme(
  bs4dash_status(primary = "#4A6C8C", info = "#3498DB", success = "#2ECC71",
                  warning = "#F39C12", danger = "#E74C3C"),
  bs4dash_color(gray_900 = "#1F2D3D")
)

footer_ui <- function() {
  tagList(
    p(
      style = "font-size: 0.72rem; color: #888; text-align: center; margin: 0 0 4px 0; width: 100%;",
      read_gif_intro_text()
    ),
    div(
      style = paste(
        "display: flex; align-items: center; justify-content: center;",
        "gap: 24px; flex-wrap: wrap; padding: 6px 0; width: 100%;"
      ),
      lapply(names(INSTITUTION_LOGOS), function(f) {
        tags$img(
          src = file.path("gif_assets", "institution_logo", f),
          height = paste0(INSTITUTION_LOGOS[[f]], "px")
        )
      })
    )
  )
}

# Custom CSS: center the brand logo in the navbar/sidebar brand box and put
# the framework name on its own line below the image instead of squeezing
# both onto one row (which caused the title text to overflow/wrap badly).
brand_css <- "
  .brand-link {
    display: flex !important;
    flex-direction: column !important;
    align-items: center !important;
    justify-content: center !important;
    height: auto !important;
    min-height: 140px;
    padding: 8px 4px !important;
    text-align: center;
    white-space: normal !important;
    cursor: pointer;
  }
  .brand-image {
    margin-right: 0 !important;
    margin-bottom: 4px !important;
    max-height: 92px;
    max-width: 92px;
  }
  .brand-text {
    font-size: 0.8rem !important;
    line-height: 1.15;
    white-space: normal !important;
    overflow: visible !important;
  }
  /* Leave room at the bottom so page content isn't hidden behind the
     fixed footer (which is taller than the AdminLTE default: caption +
     logo row). */
  .content-wrapper {
    padding-top: 16px !important;
    padding-bottom: 100px !important;
  }
  /* Only force the light-gray background in light mode -- leave dark mode's
     own (dark) content-wrapper background alone, otherwise dark-mode text
     (e.g. the 'Genomic Intelligence' heading, switched to white by bs4Dash)
     becomes invisible against this forced-light background. */
  body:not(.dark-mode) .content-wrapper {
    background-color: #f4f6f8 !important;
  }
  .badge-secondary { background-color: #adb5bd; color: #fff; }
  .badge-light { background-color: #e9ecef; color: #495057; }
  .badge-dark { background-color: #343a40; color: #fff; }
  /* Species dropdown merged inline into the Intelligence Brief card header
     bar -- borderless/transparent so it reads as one continuous line. */
  .header-inline-select {
    display: inline-flex;
    align-items: center;
    gap: 6px;
  }
  .header-inline-select .form-group {
    margin-bottom: 0 !important;
    display: inline-block;
  }
  .header-inline-select select.form-control {
    background-color: transparent !important;
    border: none !important;
    color: #fff !important;
    box-shadow: none !important;
    padding: 0 18px 0 2px !important;
    height: auto !important;
    width: auto !important;
    font-size: inherit;
    font-weight: 600;
  }
  .header-inline-select select.form-control:focus {
    outline: none !important;
    box-shadow: none !important;
  }

  /* ---- Pathogen Identification: comparison workspace ---- */
  .pi-summary-strip {
    display: flex;
    gap: 16px;
    flex-wrap: wrap;
    margin-bottom: 16px;
  }
  .pi-summary-card {
    flex: 1 1 180px;
    background: #fff;
    border: 1px solid #dee2e6;
    border-radius: 8px;
    padding: 16px 20px;
    text-align: center;
    min-width: 160px;
  }
  .pi-summary-card .pi-sc-label {
    font-size: 0.75rem;
    color: #6c757d;
    text-transform: uppercase;
    letter-spacing: 0.5px;
    margin-bottom: 4px;
  }
  .pi-summary-card .pi-sc-value {
    font-size: 1.6rem;
    font-weight: 700;
    color: #4A6C8C;
    line-height: 1.2;
  }
  .pi-summary-card .pi-sc-detail {
    font-size: 0.78rem;
    color: #888;
    margin-top: 2px;
  }
  .pi-info-card {
    flex: 1 1 220px;
    background: #eef3f8;
    border: 1px solid #c5d4e3;
    border-radius: 8px;
    padding: 16px 20px;
    display: flex;
    align-items: center;
    gap: 12px;
    min-width: 200px;
  }
  .pi-info-card .fa, .pi-info-card .fas {
    font-size: 1.2rem;
    color: #4A6C8C;
  }
  .pi-info-card span {
    font-size: 0.82rem;
    color: #4A6C8C;
  }
  .pi-chart-section {
    background: #f8f9fa;
    border: 1px solid #e9ecef;
    border-radius: 8px;
    padding: 20px 16px 12px;
    margin-bottom: 16px;
  }
  .pi-chart-section h5 {
    font-weight: 600;
    margin-bottom: 16px;
    color: #333;
  }
  .pi-legend {
    display: flex;
    justify-content: center;
    gap: 24px;
    flex-wrap: wrap;
    margin-top: 8px;
    margin-bottom: 8px;
  }
  .pi-legend-item {
    display: flex;
    align-items: center;
    gap: 6px;
    font-size: 0.82rem;
    color: #333;
  }
  .pi-legend-swatch {
    width: 14px;
    height: 14px;
    border-radius: 50%;
    display: inline-block;
  }
  .pi-sample-select-bar {
    display: flex;
    align-items: center;
    gap: 12px;
    flex-wrap: wrap;
    margin-bottom: 16px;
  }
  .pi-sample-select-bar .form-group {
    flex: 1;
    margin-bottom: 0 !important;
  }
  /* Smaller sidebar menu font */
  .sidebar-menu .nav-link {
    font-size: 0.82rem !important;
    padding: 8px 12px !important;
  }
  .sidebar-menu .nav-header {
    font-size: 0.7rem !important;
  }
  /* Smaller header & cell font for the identification results table */
  .dataTables_wrapper th {
    font-size: 0.78rem !important;
  }
  .dataTables_wrapper td {
    font-size: 0.82rem !important;
  }

  /* ---- Dark-mode overrides for Pathogen Identification section ---- */
  body.dark-mode .pi-summary-card {
    background: #2b3e50;
    border-color: #3d5469;
  }
  body.dark-mode .pi-summary-card .pi-sc-label {
    color: #adb5bd;
  }
  body.dark-mode .pi-summary-card .pi-sc-value {
    color: #e0e0e0;
  }
  body.dark-mode .pi-summary-card .pi-sc-detail {
    color: #adb5bd;
  }
  body.dark-mode .pi-info-card {
    background: #2b3e50;
    border-color: #3d5469;
  }
  body.dark-mode .pi-info-card .fa,
  body.dark-mode .pi-info-card .fas {
    color: #8ab4d6;
  }
  body.dark-mode .pi-info-card span {
    color: #cdd9e5;
  }
  body.dark-mode .pi-chart-section {
    background: #1e2d3a;
    border-color: #3d5469;
  }
  body.dark-mode .pi-chart-section h5 {
    color: #e0e0e0;
  }
  body.dark-mode .pi-legend-item {
    color: #e0e0e0;
  }
  body.dark-mode h3, body.dark-mode h4, body.dark-mode h5 {
    color: #e0e0e0;
  }
  body.dark-mode p {
    color: #cdd9e5;
  }
  /* DataTable text in dark mode */
  body.dark-mode .dataTables_wrapper th {
    color: #e0e0e0 !important;
  }
  body.dark-mode .dataTables_wrapper td {
    color: #cdd9e5 !important;
  }
  body.dark-mode .dataTables_wrapper .dataTables_info,
  body.dark-mode .dataTables_wrapper .dataTables_length label,
  body.dark-mode .dataTables_wrapper .dataTables_filter label {
    color: #adb5bd !important;
  }
  body.dark-mode .dataTables_wrapper .dataTables_paginate .paginate_button {
    color: #cdd9e5 !important;
  }
  /* Card text and captions */
  body.dark-mode .card-title {
    color: #e0e0e0 !important;
  }
  body.dark-mode caption {
    color: #adb5bd !important;
  }
  /* Selectize input in dark mode */
  body.dark-mode .selectize-input {
    background: #2b3e50 !important;
    border-color: #3d5469 !important;
    color: #e0e0e0 !important;
  }
  body.dark-mode .selectize-input .item {
    color: #e0e0e0 !important;
  }
  body.dark-mode .selectize-dropdown {
    background: #2b3e50 !important;
    border-color: #3d5469 !important;
    color: #e0e0e0 !important;
  }
  body.dark-mode .selectize-dropdown .option {
    color: #e0e0e0 !important;
  }
  body.dark-mode .selectize-dropdown .option.active {
    background: #3d5469 !important;
  }
"

# -----------------------------------------------------------------------------
# UI
# -----------------------------------------------------------------------------
ui <- bs4DashPage(
  freshTheme = gif_theme,
  title = "Genomic Intelligence Framework",
  header = bs4DashNavbar(
    title = bs4DashBrand(
      title = "Genomic Intelligence Framework",
      image = file.path("gif_assets", "logo.png"),
      color = "primary"
    )
  ),
  sidebar = bs4DashSidebar(
    status = "primary",
    sidebarMenuOutput("sidebarmenu")
  ),
  body = bs4DashBody(
    tags$head(
      tags$style(HTML(brand_css)),
      tags$script(HTML(
        "$(document).on('click', '.brand-link', function(e) {
           e.preventDefault();
           Shiny.setInputValue('brand_home_click', 'click', {priority: 'event'});
         });
         // Track dark-mode toggle and expose to Shiny
         $(function() {
           var sendDarkMode = function() {
             var isDark = $('body').hasClass('dark-mode');
             Shiny.setInputValue('is_dark_mode', isDark);
           };
           // Observe class changes on body
           var observer = new MutationObserver(function(mutations) {
             sendDarkMode();
           });
           observer.observe(document.body, {attributes: true, attributeFilter: ['class']});
           // Initial state once Shiny is ready
           $(document).on('shiny:connected', sendDarkMode);
         });"
      ))
    ),
    uiOutput("body_ui")
  ),
  footer = bs4DashFooter(
    left = footer_ui(),
    right = NULL,
    fixed = TRUE
  )
)

# -----------------------------------------------------------------------------
# Server
# -----------------------------------------------------------------------------
server <- function(input, output, session) {

  outdir_r <- reactive({ "../results" })
  species_rv <- reactive({ list_species(outdir_r()) })

  # -- Sidebar: Intelligence Overview + Biological Threat (species submenu,
  # formerly "Pathogen Identification") + the 6 remaining objective roadmap
  # items + a divider + Evidence & Knowledge Gaps / Intelligence Brief.
  # Rebuilt whenever the discovered species list changes. Compare is
  # intentionally not included (see plan).
  output$sidebarmenu <- renderMenu({
    species <- species_rv()

    pi_item <- if (length(species) > 0) {
      do.call(menuItem, c(
        list(text = "Biological Threat", icon = icon("dna"), startExpanded = TRUE),
        lapply(species, function(sp) {
          menuSubItem(text = toupper(sp), tabName = paste0("pi_", sp))
        })
      ))
    } else {
      menuItem("Biological Threat", tabName = "pi_home", icon = icon("dna"))
    }

    roadmap_items <- lapply(ROADMAP_OBJECTIVES, function(obj) {
      menuItem(obj$label, tabName = roadmap_tab_name(obj$id), icon = icon(obj$icon))
    })

    sidebarMenu(
      id = "sidebarmenu",
      menuItem("Intelligence Overview", tabName = "home", icon = icon("house")),
      sidebarHeader("INTELLIGENCE OBJECTIVES"),
      pi_item,
      roadmap_items[1:6],
      sidebarHeader("EVIDENCE & REPORTING"),
      roadmap_items[7:8]
    )
  })

  # -- Body: one tabItem per sidebar entry, rebuilt in lockstep with the menu.
  output$body_ui <- renderUI({
    species <- species_rv()

    pi_tabs <- if (length(species) > 0) {
      lapply(species, function(sp) {
        tabItem(tabName = paste0("pi_", sp), pathogen_identification_ui(sp))
      })
    } else {
      list(tabItem(
        tabName = "pi_home",
        bs4Dash::bs4Card(
          title = "Biological Threat", width = 12, status = "secondary",
          div(
            style = "text-align: center; padding: 40px 0; color: #888;",
            icon("dna", class = "fa-3x"),
            h4("No species found", style = "margin-top: 16px;"),
            p("No species were found in the pipeline output. Run the pipeline, or check that results/ exists.")
          )
        )
      ))
    }

    roadmap_tabs <- lapply(ROADMAP_OBJECTIVES, function(obj) {
      tab_ui <- if (identical(obj$id, "intelligence_brief")) {
        intelligence_brief_roadmap_ui(obj)
      } else {
        roadmap_ui(obj)
      }
      tabItem(tabName = roadmap_tab_name(obj$id), tab_ui)
    })

    do.call(tabItems, c(
      list(tabItem(tabName = "home", overview_ui())),
      pi_tabs,
      roadmap_tabs
    ))
  })

  # -- Register renderDT outputs for every discovered species.
  observeEvent(species_rv(), {
    for (sp in species_rv()) {
      pathogen_identification_register(input, output, session, sp, outdir_r)
    }
  }, ignoreNULL = FALSE)

  # -- Intelligence Overview: header (species dropdown, merged inline into
  # the card header bar) + body. Header depends only on species_rv() (not on
  # input$overview_species) so it doesn't re-render -- and lose the user's
  # in-progress selection -- every time the dropdown changes.
  output$overview_header <- renderUI({
    overview_header_ui(species_rv())
  })

  current_species <- reactive({
    species <- species_rv()
    if (length(species) == 0) return(NULL)
    sel <- input$overview_species
    if (is.null(sel) || !(sel %in% species)) species[1] else sel
  })

  output$overview_body <- renderUI({
    overview_body_ui(!is.null(current_species()))
  })

  # -- Overview CTAs and signal/roadmap-tile clicks jump the sidebar tab.
  observeEvent(input$overview_explore_evidence, {
    sp <- current_species()
    if (!is.null(sp)) updateTabItems(session, "sidebarmenu", selected = paste0("pi_", sp))
  })

  observeEvent(input$overview_generate_brief, {
    updateTabItems(session, "sidebarmenu", selected = roadmap_tab_name("intelligence_brief"))
  })

  observeEvent(input$brand_home_click, {
    updateTabItems(session, "sidebarmenu", selected = "home")
  })

  observeEvent(input$overview_live_tile_click, {
    sp <- current_species()
    if (!is.null(sp)) updateTabItems(session, "sidebarmenu", selected = paste0("pi_", sp))
  })

  observeEvent(input$overview_nav_click, {
    updateTabItems(session, "sidebarmenu", selected = input$overview_nav_click)
  })

  # -- Intelligence Brief "preview & edit before export" scaffold (UI only).
  observeEvent(input$brief_preview_open, {
    showModal(intelligence_brief_preview_modal())
  })
}

shinyApp(ui, server)
