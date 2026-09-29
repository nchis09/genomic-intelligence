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

# -- Auto-install missing packages on first run --
local({
  cran_pkgs <- c("shiny", "bs4Dash", "fresh", "DT", "readr", "dplyr", "tidyr",
                 "ggplot2", "ggrepel", "gridExtra", "scales", "plotly", "yaml",
                 "DBI", "duckdb", "leaflet", "httr", "jsonlite",
                 "RColorBrewer", "ape", "ggnewscale")
  bioc_pkgs <- c("treeio", "ggtree")

  missing_cran <- cran_pkgs[!sapply(cran_pkgs, requireNamespace, quietly = TRUE)]
  missing_bioc <- bioc_pkgs[!sapply(bioc_pkgs, requireNamespace, quietly = TRUE)]

  if (length(missing_cran) > 0) {
    message("[dashboard] Installing CRAN packages: ", paste(missing_cran, collapse = ", "))
    install.packages(missing_cran, repos = "https://cloud.r-project.org", quiet = TRUE)
  }

  if (length(missing_bioc) > 0) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "https://cloud.r-project.org", quiet = TRUE)
    }
    message("[dashboard] Installing Bioconductor packages: ", paste(missing_bioc, collapse = ", "))
    BiocManager::install(missing_bioc, ask = FALSE, update = FALSE, quiet = TRUE)
  }
})

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
source("modules/llm_note.R")
source("modules/pathogen_genomics.R")
source("modules/pathogen_mutation_profile.R")
source("modules/assessment.R")
source("modules/placeholder.R")
source("modules/pathogen_transmission.R")
source("modules/pathogen_geographic.R")
source("modules/countermeasure_metaanalysis.R")
source("modules/evidence_metaanalysis.R")
source("modules/domain_views.R")
source("modules/countermeasure_readiness.R")
source("modules/evidence_knowledge_gaps.R")
source("modules/intelligence_brief.R")
source("modules/home.R")

# Serve repo-root assets (GIF logo, institution logos) without copying them
# into dashboard/www.
addResourcePath("gif_assets", normalizePath(file.path("..", "assets")))



list_species <- function(outdir) {
  base <- file.path(outdir, "pathogen_identification")
  if (!dir.exists(base)) return(character())
  basename(list.dirs(base, recursive = FALSE))
}

# Resolve the pipeline --outdir the user typed into a real path. Accepts an
# absolute path, a path relative to dashboard/, or a bare directory name that
# lives at the repo root (e.g. "results" -> "../results").
resolve_outdir <- function(x) {
  if (is.null(x) || !nzchar(trimws(x))) return("../results")
  x <- trimws(x)
  if (grepl("^(/|~)", x)) return(path.expand(x))
  if (dir.exists(x)) return(x)
  cand <- file.path("..", x)
  if (dir.exists(cand)) return(cand)
  x
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
  /* Indent Phylotree/Mutation under Pathogen Genomics */
  .pg-child-items > .nav-item {
    padding-left: 12px;
    background-color: rgba(0, 0, 0, 0.03);
    border-left: 3px solid #4A6C8C;
  }
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
  /* ---- Dark-mode overrides for Pathogen Genomics section ---- */
  body.dark-mode .radio label,
  body.dark-mode .checkbox label {
    color: #cdd9e5 !important;
  }
  body.dark-mode .control-label {
    color: #adb5bd !important;
  }
  body.dark-mode .shiny-input-container .radio-inline label {
    color: #cdd9e5 !important;
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
    div(
      style = "padding: 10px 12px 0 12px;",
      textInput(
        inputId = "outdir",
        label = "Pipeline --outdir",
        value = "results",
        placeholder = "e.g. results"
      )
    ),
    sidebarMenuOutput("sidebarmenu")
  ),
  body = bs4DashBody(
    tags$head(
      tags$style(HTML(brand_css)),
      tags$script(HTML(
        "// Brief export: server sends a complete HTML doc to print
        Shiny.addCustomMessageHandler('open_print_window', function(m) {
          var w = window.open('', '_blank');
          if (!w) return;
          w.document.open(); w.document.write(m.html); w.document.close();
          w.focus(); w.print();
        });
        $(document).on('click', '.brand-link', function(e) {
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
         });
         // Auto-collapse other sidebar tree items when one is expanded
         $(document).on('click', '.sidebar .nav-item.has-treeview > a', function() {
           var $this = $(this).parent();
           $('.sidebar .nav-item.has-treeview.menu-open').not($this).each(function() {
             $(this).removeClass('menu-open');
             $(this).children('.nav-treeview').css('display', 'none');
           });
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

  outdir_r <- reactive({ resolve_outdir(input$outdir) })
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
        list(text = "Biological Threat", icon = icon("dna"), startExpanded = FALSE),
        lapply(species, function(sp) {
          menuSubItem(text = toupper(sp), tabName = paste0("pi_", sp))
        })
      ))
    } else {
      menuItem("Biological Threat", tabName = "pi_home", icon = icon("dna"))
    }

    pg_items <- if (length(species) > 0) {
      list(
        menuItem(text = "Pathogen Genomics", tabName = "pg_home", icon = icon("microscope")),
        tags$div(
          class = "pg-child-items",
          do.call(menuItem, c(
            list(text = "Phylotree", icon = icon("share-nodes"), startExpanded = FALSE),
            lapply(species, function(sp) {
              menuSubItem(text = toupper(sp), tabName = paste0("pg_", sp))
            })
          )),
          do.call(menuItem, c(
            list(text = "Mutation", icon = icon("chart-bar"), startExpanded = FALSE),
            lapply(species, function(sp) {
              menuSubItem(text = toupper(sp), tabName = paste0("mp_", sp))
            })
          ))
        )
      )
    } else {
      list(menuItem("Pathogen Genomics", tabName = "pg_home", icon = icon("microscope")))
    }

    roadmap_items <- lapply(ROADMAP_OBJECTIVES, function(obj) {
      if (identical(obj$id, "transmission_spread") && length(species) > 0) {
        do.call(menuItem, c(
          list(text = obj$label, icon = icon(obj$icon), startExpanded = FALSE),
          lapply(species, function(sp) menuSubItem(text = toupper(sp), tabName = paste0("ts_", sp)))
        ))
      } else if (identical(obj$id, "geographic_temporal") && length(species) > 0) {
        do.call(menuItem, c(
          list(text = obj$label, icon = icon(obj$icon), startExpanded = FALSE),
          lapply(species, function(sp) menuSubItem(text = toupper(sp), tabName = paste0("gt_", sp)))
        ))
      } else if (identical(obj$id, "countermeasure_readiness") && length(species) > 0) {
        do.call(menuItem, c(
          list(text = obj$label, icon = icon(obj$icon), startExpanded = FALSE),
          lapply(species, function(sp) menuSubItem(text = toupper(sp), tabName = paste0("cm_", sp)))
        ))
      } else if (identical(obj$id, "evidence_knowledge_gaps") && length(species) > 0) {
        do.call(menuItem, c(
          list(text = obj$label, icon = icon(obj$icon), startExpanded = FALSE),
          lapply(species, function(sp) menuSubItem(text = toupper(sp), tabName = paste0("ekg_", sp)))
        ))
      } else {
        menuItem(obj$label, tabName = roadmap_tab_name(obj$id), icon = icon(obj$icon))
      }
    })

    do.call(sidebarMenu, c(
      list(id = "sidebarmenu"),
      list(menuItem("Intelligence Overview", tabName = "home", icon = icon("house"))),
      list(sidebarHeader("INTELLIGENCE OBJECTIVES")),
      list(pi_item),
      pg_items,
      roadmap_items[1:4],
      list(sidebarHeader("EVIDENCE & REPORTING")),
      roadmap_items[5:6]
    ))
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

    pg_tabs <- if (length(species) > 0) {
      c(
        list(tabItem(tabName = "pg_home", pg_home_ui(species))),
        unlist(
          lapply(species, function(sp) {
            list(
              tabItem(tabName = paste0("pg_", sp), pathogen_genomics_ui(sp)),
              tabItem(tabName = paste0("mp_", sp), pathogen_mutation_profile_ui(sp))
            )
          }),
          recursive = FALSE
        )
      )
    } else {
      list(tabItem(
        tabName = "pg_home",
        bs4Dash::bs4Card(
          title = "Pathogen Genomics", width = 12, status = "secondary",
          div(
            style = "text-align: center; padding: 40px 0; color: #888;",
            icon("microscope", class = "fa-3x"),
            h4("No species found", style = "margin-top: 16px;"),
            p("No species were found in the pipeline output. Run the pipeline, or check that results/ exists.")
          )
        )
      ))
    }

    roadmap_tabs <- lapply(ROADMAP_OBJECTIVES, function(obj) {
      tab_ui <- if (identical(obj$id, "intelligence_brief")) {
        uiOutput("brief_tab_body")
      } else if (identical(obj$id, "transmission_spread")) {
        uiOutput("ts_roadmap_body")
      } else if (identical(obj$id, "geographic_temporal")) {
        uiOutput("gt_roadmap_body")
      } else if (identical(obj$id, "countermeasure_readiness")) {
        uiOutput("cm_roadmap_body")
      } else if (identical(obj$id, "evidence_knowledge_gaps")) {
        uiOutput("ekg_roadmap_body")
      } else {
        roadmap_ui(obj)
      }
      tabItem(tabName = roadmap_tab_name(obj$id), tab_ui)
    })

    ts_tabs <- lapply(species, function(sp) {
      tabItem(tabName = paste0("ts_", sp), transmission_ui(sp))
    })
    gt_tabs <- lapply(species, function(sp) {
      tabItem(tabName = paste0("gt_", sp), geographic_temporal_ui(sp, species))
    })
    cm_tabs <- lapply(species, function(sp) {
      tabItem(tabName = paste0("cm_", sp), countermeasure_readiness_ui(sp))
    })
    ekg_tabs <- lapply(species, function(sp) {
      tabItem(tabName = paste0("ekg_", sp), evidence_knowledge_gaps_ui(sp))
    })

    do.call(tabItems, c(
      list(tabItem(tabName = "home", overview_ui())),
      pi_tabs,
      pg_tabs,
      ts_tabs,
      gt_tabs,
      cm_tabs,
      ekg_tabs,
      roadmap_tabs
    ))
  })

  # -- Register renderDT outputs for every discovered species.
  observeEvent(species_rv(), {
    for (sp in species_rv()) {
      pathogen_identification_register(input, output, session, sp, outdir_r)
      pathogen_genomics_register(input, output, session, sp, outdir_r)
      pathogen_mutation_profile_register(input, output, session, sp, outdir_r)
      transmission_register(input, output, session, sp, outdir_r)
      geographic_temporal_register(input, output, session, sp, outdir_r)
      countermeasure_readiness_register(input, output, session, sp, outdir_r, species_rv)
      evidence_knowledge_gaps_register(input, output, session, sp, outdir_r, species_rv)

      local({
        s <- sp
        observeEvent(input[[paste0("pg_home_phylotree_", s)]], {
          updateTabItems(session, "sidebarmenu", selected = paste0("pg_", s))
        })
        observeEvent(input[[paste0("pg_home_mutation_", s)]], {
          updateTabItems(session, "sidebarmenu", selected = paste0("mp_", s))
        })
      })
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

  # Roadmap "Transmission & Spread" tile tab -> shows the currently selected
  # species' transmission module (per-species instances are registered above).
  output$ts_roadmap_body <- renderUI({
    sp <- current_species(); if (is.null(sp)) return(NULL)
    transmission_ui(sp)
  })
  # Roadmap "Geographic & Temporal Context" tile tab -> shows the currently
  # selected species' geographic module (per-species instances registered above).
  output$gt_roadmap_body <- renderUI({
    sp <- current_species(); if (is.null(sp)) return(NULL)
    geographic_temporal_ui(sp, species_rv())
  })
  output$cm_roadmap_body <- renderUI({
    sp <- current_species(); if (is.null(sp)) return(NULL)
    countermeasure_readiness_ui(sp)
  })
  output$ekg_roadmap_body <- renderUI({
    sp <- current_species(); if (is.null(sp)) return(NULL)
    evidence_knowledge_gaps_ui(sp)
  })

  output$overview_body <- renderUI({
    overview_body_ui(current_species(), outdir_r())
  })

  # -- Home-page genomics mini-plots (divergence + VIP) wired into the
  # "What is different about it" card. Divergence is read from the tree_tips
  # table of the nextstrain tree; it is a per-sample evolution metric.
  output$genomics_burden_plot <- renderPlot({
    sp <- current_species(); outdir <- outdir_r(); if (is.null(sp) || is.null(outdir)) return(NULL)
    db <- file.path(outdir, "knowledge_warehouse", "knowledge_warehouse.duckdb")
    if (!file.exists(db)) return(NULL)
    con <- DBI::dbConnect(duckdb::duckdb(), dbdir = db, read_only = TRUE)
    on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
    tree_id <- DBI::dbGetQuery(con,
      "SELECT tree_id FROM phylogenetic_trees WHERE species = ? AND tree_method = 'nextstrain' ORDER BY tree_id LIMIT 1",
      params = list(sp))
    if (!nrow(tree_id))
      tree_id <- DBI::dbGetQuery(con,
        "SELECT tree_id FROM phylogenetic_trees WHERE species = ? LIMIT 1",
        params = list(sp))
    if (!nrow(tree_id)) return(NULL)
    tips <- DBI::dbGetQuery(con,
      "SELECT is_query, div FROM tree_tips WHERE tree_id = ? AND div IS NOT NULL",
      params = list(tree_id$tree_id))
    if (!nrow(tips)) return(NULL)
    tips$group <- ifelse(tips$is_query, "Query", "Background")
    tips$group <- factor(tips$group, levels = c("Background", "Query"))
    query <- tips[tips$group == "Query", ]
    ggplot2::ggplot(tips, aes(group, div, fill = group)) +
      ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.45, show.legend = FALSE) +
      ggplot2::geom_jitter(data = query, color = "#C0392B", width = 0.1, height = 0, size = 2.2) +
      ggplot2::scale_fill_manual(values = c(Query = "#C0392B", Background = "#4A6C8C")) +
      ggplot2::labs(x = NULL, y = "Divergence (substs/site)", title = "Query vs background divergence") +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(legend.position = "none")
  })

  output$genomics_vip_plot <- renderPlot({
    sp <- current_species(); outdir <- outdir_r(); if (is.null(sp) || is.null(outdir)) return(NULL)
    df <- .brief_tsv(file.path(outdir, "pathogen_mutation_profile", sp, "mutation_profile"),
                     "01_plsda_vip.tsv")
    if (is.null(df) || !nrow(df)) return(NULL)
    df <- df[order(df$vip, decreasing = TRUE), ]
    if (nrow(df) > 10) df <- utils::head(df, 10)
    df$protein_name <- factor(df$protein_name, levels = rev(df$protein_name))
    ggplot2::ggplot(df, aes(protein_name, vip)) +
      ggplot2::geom_col(fill = "#4A6C8C", width = 0.7, show.legend = FALSE) +
      ggplot2::coord_flip() +
      ggplot2::labs(x = NULL, y = "VIP score", title = "Top proteins") +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(legend.position = "none",
                     plot.title = ggplot2::element_text(size = 10, face = "bold"))
  })

  # -- Home-page projection fan chart (Where it could go card).
  output$projection_fan_plot <- renderPlot({
    sp <- current_species(); outdir <- outdir_r(); if (is.null(sp) || is.null(outdir)) return(NULL)
    pj <- .brief_tsv(file.path(outdir, "transmission_context", sp), "query_projection.tsv")
    if (is.null(pj) || !nrow(pj)) return(NULL)
    pj <- pj |>
      dplyr::mutate(
        wk = as.numeric(.data$week_ahead),
        med = as.numeric(.data$proj_median),
        lo = as.numeric(.data$proj_lower95),
        hi = as.numeric(.data$proj_upper95),
        scenario = as.character(.data$scenario)
      ) |>
      dplyr::filter(!is.na(.data$wk), !is.na(.data$med)) |>
      dplyr::group_by(.data$scenario, .data$wk) |>
      dplyr::summarise(med = median(.data$med), lo = min(.data$lo), hi = max(.data$hi), .groups = "drop")
    pal <- c(baseline = "#4A6C8C", contained = "#3A916E", expanded = "#C0392B")
    ggplot2::ggplot(pj, aes(wk, med, color = scenario, fill = scenario)) +
      ggplot2::geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, color = NA) +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::scale_color_manual(values = pal) +
      ggplot2::scale_fill_manual(values = pal) +
      ggplot2::labs(x = "Weeks ahead", y = "Weekly cases", color = "Scenario", title = "8-week projection") +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(legend.position = "top",
                     legend.title = ggplot2::element_blank(),
                     legend.text = ggplot2::element_text(size = 7))
  })

  # -- Home-page spread map (How it's spreading card): query locations,
  # likely origins, and reported historical cases from transmission_spatial.
  output$spread_map_plot <- leaflet::renderLeaflet({
    brief <- brief_data()
    outdir <- outdir_r()
    if (is.null(brief) || is.null(outdir)) return(NULL)
    q_countries <- unique(brief$situation$query_countries %||% character())
    origins <- unique(brief$transmission$likely_origins %||% character())
    all_countries <- unique(c(q_countries, origins))
    if (!length(all_countries)) return(NULL)
    cc <- .brief_country_centroid(all_countries)
    cc$country <- all_countries
    cc <- cc[!is.na(cc$lat) & !is.na(cc$lon), ]
    q <- cc[cc$country %in% q_countries, ]
    o <- cc[cc$country %in% origins, ]
    m <- leaflet::leaflet() |> leaflet::addTiles()
    if (nrow(q)) {
      m <- m |> leaflet::addCircleMarkers(
        data = q, lng = ~lon, lat = ~lat, color = "#C0392B",
        fillColor = "#C0392B", fillOpacity = 0.85, radius = 7,
        popup = ~paste0("Query location: ", country)
      )
    }
    if (nrow(o)) {
      m <- m |> leaflet::addCircleMarkers(
        data = o, lng = ~lon, lat = ~lat, color = "#4A6C8C",
        fillColor = "#4A6C8C", fillOpacity = 0.85, radius = 7,
        popup = ~paste0("Likely origin: ", country)
      )
      for (i in seq_len(nrow(o))) {
        for (j in seq_len(nrow(q))) {
          m <- m |> leaflet::addPolylines(
            lng = c(o$lon[i], q$lon[j]), lat = c(o$lat[i], q$lat[j]),
            color = "#4A6C8C", dashArray = "5,8", opacity = 0.6, weight = 1.5
          )
        }
      }
    }
    if (nrow(cc)) {
      pad <- 8
      m <- m |> leaflet::fitBounds(
        min(cc$lon) - pad, min(cc$lat) - pad,
        max(cc$lon) + pad, max(cc$lat) + pad
      )
    }
    m
  })

  # -- Intelligence Brief tab: renders the pre-generated per-species
  # intelligence_brief.json for whichever species the Overview dropdown has
  # selected, plus the export card.
  brief_data <- reactive({
    sp <- current_species()
    if (is.null(sp)) return(NULL)
    brief_read(outdir_r(), sp)
  })

  output$brief_tab_body <- renderUI({
    sp <- current_species()
    if (is.null(sp)) {
      return(roadmap_ui(ROADMAP_OBJECTIVES[[which(vapply(ROADMAP_OBJECTIVES,
        function(o) identical(o$id, "intelligence_brief"), logical(1)))]]))
    }
    tagList(
      brief_full_ui(brief_data(), sp, outdir_r()),
      bs4Dash::bs4Card(
        title = "Export", width = 12, status = "secondary",
        p(style = "color: #6c757d;",
          "Open an editable preview of the brief before printing or sharing."),
        actionButton("brief_preview_open", "Preview & Edit Brief",
                     icon = icon("file-lines"))
      )
    )
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

  observeEvent(input$overview_pg_tile_click, {
    updateTabItems(session, "sidebarmenu", selected = "pg_home")
  })

  observeEvent(input$overview_nav_click, {
    tab <- input$overview_nav_click
    sp <- current_species()
    # transmission/geographic are now per-species expandable menus -> route the
    # roadmap tile to the current species' subitem tab.
    if (identical(tab, roadmap_tab_name("transmission_spread")) && !is.null(sp)) {
      tab <- paste0("ts_", sp)
    } else if (identical(tab, roadmap_tab_name("geographic_temporal")) && !is.null(sp)) {
      tab <- paste0("gt_", sp)
    } else if (identical(tab, roadmap_tab_name("countermeasure_readiness")) && !is.null(sp)) {
      tab <- paste0("cm_", sp)
    } else if (identical(tab, roadmap_tab_name("evidence_knowledge_gaps")) && !is.null(sp)) {
      tab <- paste0("ekg_", sp)
    }
    updateTabItems(session, "sidebarmenu", selected = tab)
  })

  # -- Intelligence Brief "preview & edit before export": the modal is
  # pre-filled with the generated brief as editable plain text.
  observeEvent(input$brief_preview_open, {
    sp <- current_species()
    showModal(intelligence_brief_preview_modal(brief_markdown(brief_data(), sp)))
  })

  # Sign-off metadata entered in the preview modal, collected once here.
  brief_signoff_meta <- reactive({
    list(lab          = input$report_lab_name,
         run_by       = input$report_run_by,
         validated_by = input$report_validated_by,
         report_date  = as.character(input$report_date %||% Sys.Date()),
         species      = current_species())
  })

  # Download the (possibly edited) brief text as a .md file, sign-off appended.
  output$brief_download <- downloadHandler(
    filename = function() {
      sp <- current_species() %||% "brief"
      paste0("intelligence_brief_", sp, "_", Sys.Date(), ".md")
    },
    content = function(file) {
      m <- brief_signoff_meta()
      writeLines(c(input$brief_preview_text %||% "",
                   brief_signoff_md(m$lab, m$run_by, m$validated_by, m$report_date)),
                 file)
    }
  )

  # Download the formatted standalone HTML report (same doc the print path uses).
  output$brief_download_html <- downloadHandler(
    filename = function() {
      sp <- current_species() %||% "brief"
      paste0("intelligence_brief_", sp, "_", Sys.Date(), ".html")
    },
    content = function(file) {
      writeLines(brief_report_html(input$brief_preview_text, brief_signoff_meta()), file)
    }
  )

  # Print / Save as PDF: build the formatted document server-side and hand it
  # to the browser print dialog (head JS handler 'open_print_window').
  observeEvent(input$brief_export_pdf, {
    session$sendCustomMessage("open_print_window",
      list(html = brief_report_html(input$brief_preview_text, brief_signoff_meta())))
  })
}

shinyApp(ui, server)
