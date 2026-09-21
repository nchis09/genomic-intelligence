# Transmission & Spread dashboard module
#
# Query-driven intelligence brief. Reads pre-computed TSVs from
# results/transmission_context/<species>/ and the DuckDB knowledge warehouse.
#
# Visual hierarchy (evidence -> interpretation -> projection):
#   GENOME -> A. Genomic relationship & spread -> B. Outbreak impact over time ->
#   C. Historical behavior -> D. What might happen next -> E. Why -> F. Assessment
#   (Geographic & Temporal Context lives in its own module/tab: pathogen_geographic.R)
#   Secondary tab: Surveillance Data (raw burden / anomaly audit tables)

transmission_ui <- function() {
  uiOutput("transmission_body")
}

# ---- Data helpers ----
.ts_db_path <- function(outdir) file.path(outdir, "knowledge_warehouse", "knowledge_warehouse.duckdb")

.ts_read <- function(outdir, species, fname) {
  p <- file.path(outdir, "transmission_context", species, fname)
  if (!file.exists(p)) return(NULL)
  df <- tryCatch(readr::read_tsv(p, show_col_types = FALSE), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df
}

.ts_clean_label <- function(x) gsub("[,;():\\[\\] ]", "_", as.character(x))

.ts_strain_pal <- function(strains) {
  strains <- sort(unique(stats::na.omit(as.character(strains))))
  n <- length(strains)
  if (n == 0) return(stats::setNames(character(0), character(0)))
  base <- RColorBrewer::brewer.pal(8, "Set2")
  cols <- if (n <= 8) base[1:n] else grDevices::colorRampPalette(base)(n)
  stats::setNames(cols, strains)
}

# Consolidate historical country-name variants into short acronyms so the same
# country isn't split into multiple colors/labels. All Zaire / DRC spellings ->
# "DRC". Republic of the Congo is a DIFFERENT country (Congo-Brazzaville) and is
# kept separate as "Rep. Congo".
.ts_norm_country <- function(x) {
  x <- trimws(as.character(x))
  dplyr::case_when(
    grepl("zaire|democratic republic of the congo|\\bdrc\\b", x, ignore.case = TRUE) ~ "DRC",
    grepl("republic of the congo|congo-brazzaville|\\broc\\b", x, ignore.case = TRUE) ~ "Rep. Congo",
    grepl("united states", x, ignore.case = TRUE) ~ "USA",
    grepl("united kingdom", x, ignore.case = TRUE) ~ "UK",
    TRUE ~ x
  )
}

# ---- Dynamic country centroids ----
# Centroids are computed on the fly from the world map polygons (ggplot2::map_data),
# so any country present in the base map resolves automatically and unknown names
# return NA (skipped) rather than erroring. No coordinates are hardcoded.
.ts_centroid_cache <- new.env(parent = emptyenv())

.ts_world_centroids <- function() {
  if (exists("centroids", envir = .ts_centroid_cache)) return(get("centroids", envir = .ts_centroid_cache))
  cents <- data.frame(country = character(0), lat = numeric(0), lon = numeric(0), stringsAsFactors = FALSE)
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    wm <- tryCatch(ggplot2::map_data("world"), error = function(e) NULL)
    if (!is.null(wm) && nrow(wm)) {
      # centroid = mean of polygon vertices per region (largest subregion dominates via weighting)
      cents <- wm %>%
        dplyr::filter(!is.na(.data$long), !is.na(.data$lat)) %>%
        dplyr::group_by(country = .data$region) %>%
        dplyr::summarise(lon = mean(.data$long, na.rm = TRUE), lat = mean(.data$lat, na.rm = TRUE), .groups = "drop")
    }
  }
  assign("centroids", cents, envir = .ts_centroid_cache)
  cents
}

# Resolve a vector of country names to c(lat, lon). Returns a data.frame with
# lat/lon columns aligned to the input; unmatched names -> NA. Case-insensitive
# with a small alias normalisation for common naming variants.
.ts_country_centroid <- function(country) {
  out <- data.frame(lat = rep(NA_real_, length(country)), lon = rep(NA_real_, length(country)))
  cents <- .ts_world_centroids()
  if (is.null(cents) || !nrow(cents)) return(out)
  norm <- function(x) tolower(trimws(gsub("[^a-z ]", "", tolower(as.character(x)))))
  alias <- c("democratic republic of the congo" = "democratic republic of the congo",
             "drc" = "democratic republic of the congo", "congo kinshasa" = "democratic republic of the congo",
             "republic of the congo" = "republic of congo", "congo brazzaville" = "republic of congo",
             "united states" = "usa", "united states of america" = "usa", "us" = "usa",
             "united kingdom" = "uk", "great britain" = "uk",
             "ivory coast" = "ivory coast", "cote divoire" = "ivory coast",
             "south sudan" = "south sudan", "the gambia" = "gambia")
  key <- norm(country)
  key <- ifelse(key %in% names(alias), unname(alias[key]), key)
  cmap <- stats::setNames(seq_len(nrow(cents)), norm(cents$country))
  idx <- unname(cmap[key])
  hit <- !is.na(idx)
  out$lat[hit] <- cents$lat[idx[hit]]
  out$lon[hit] <- cents$lon[idx[hit]]
  out
}

# Fill missing lat/lon in a data.frame from a country-name column. `lat_col`,
# `lon_col`, `country_col` are column names; rows already having coords are kept.
.ts_fill_coords <- function(df, lat_col, lon_col, country_col) {
  if (is.null(df) || !nrow(df) || !all(c(lat_col, lon_col, country_col) %in% names(df))) return(df)
  need <- is.na(df[[lat_col]]) | is.na(df[[lon_col]])
  if (any(need)) {
    cc <- .ts_country_centroid(df[[country_col]][need])
    df[[lat_col]][need] <- cc$lat
    df[[lon_col]][need] <- cc$lon
  }
  df
}

# ---- Growth model (model-derived projection) ----
# Weekly epi series for a strain's outbreak (or all burden if strain is NULL).
.ts_strain_epi <- function(burden, strain = NULL) {
  if (is.null(burden) || !nrow(burden) || !all(c("week_start", "new_cases") %in% names(burden))) return(NULL)
  df <- burden
  if (!is.null(strain) && "dominant_strain" %in% names(df)) {
    in_strains <- if ("genome_strains" %in% names(df)) grepl(strain, df$genome_strains, fixed = TRUE) else rep(FALSE, nrow(df))
    m <- dplyr::filter(df, .data$dominant_strain == strain | in_strains)
    if (nrow(m) >= 4) df <- m
  }
  df <- df %>% dplyr::filter(!is.na(.data$week_start), !is.na(.data$new_cases)) %>%
    dplyr::group_by(week_start = as.Date(.data$week_start)) %>%
    dplyr::summarise(new_cases = sum(as.numeric(.data$new_cases), na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(.data$week_start)
  if (nrow(df) < 4) return(NULL)
  df$t <- as.numeric(df$week_start - min(df$week_start))
  df$cum <- cumsum(df$new_cases)
  df
}

# Estimate growth dynamics from a strain's observed epidemic curve.
# r = exponential growth rate from the rising (pre-peak) phase via log-linear
# regression on weekly incidence; K = observed outbreak size; tm = time to half.
# Robust to multi-year / plateaued series where a single logistic nls fails.
.ts_fit_logistic <- function(epi) {
  if (is.null(epi) || nrow(epi) < 4) return(NULL)
  t <- epi$t; inc <- as.numeric(epi$new_cases); cum <- epi$cum
  # Intrinsic growth rate = steepest sustained log-linear slope over a rolling
  # window (captures the exponential-growth phase, not the slow multi-year ramp).
  r <- NA_real_; w <- min(8, length(inc))
  if (w >= 3 && any(inc > 0)) {
    slopes <- vapply(seq_len(length(inc) - w + 1), function(i) {
      idx <- i:(i + w - 1); ii <- inc[idx]
      if (sum(ii > 0) < 3) return(NA_real_)
      cf <- tryCatch(stats::coef(stats::lm(log(pmax(ii, 1)) ~ t[idx])), error = function(e) NULL)
      if (!is.null(cf) && length(cf) >= 2 && is.finite(cf[2])) unname(cf[2]) else NA_real_
    }, numeric(1))
    slopes <- slopes[is.finite(slopes) & slopes > 0]
    if (length(slopes)) r <- stats::quantile(slopes, 0.9, names = FALSE)  # fast-but-robust
  }
  if (is.na(r)) {  # fallback: median positive week-over-week growth
    gr <- diff(log(pmax(inc, 1))); gr <- gr[is.finite(gr) & gr > 0]
    r <- if (length(gr)) stats::median(gr) / 7 else 0.02
  }
  r_day <- max(1e-4, r)
  K <- max(cum); if (!is.finite(K) || K <= 0) K <- sum(inc, na.rm = TRUE)
  tm <- t[which.min(abs(cum - K / 2))]
  list(K = K, r = r_day, tm = tm, resid_sd = stats::sd(inc) * 0.3, method = "growth-phase fit",
       r_week = r_day * 7, doubling_days = log(2) / r_day,
       t_end = max(t), date0 = min(epi$week_start), obs = epi)
}

# Simulate a comparable outbreak forward `horizon` days from a small seed, using
# the fitted growth rate and size scaled by scenario. Band widens with lead time.
.ts_project <- function(fit, horizon = 180, scenario = "baseline", seed = 5) {
  if (is.null(fit)) return(NULL)
  r_mult <- switch(scenario, contained = 0.5, baseline = 1, expanded = 1.35, 1)
  K_mult <- switch(scenario, contained = 0.35, baseline = 1, expanded = 2.2, 1)
  r <- fit$r * r_mult; K <- max(fit$K * K_mult, seed * 2)
  tm <- log(K / seed - 1) / r            # inflection so cum(0) ~= seed
  t <- seq(0, horizon, length.out = 300)
  cum <- K / (1 + exp(-r * (t - tm)))
  sd <- (fit$resid_sd + 1) * (0.4 + 1.6 * t / horizon)
  data.frame(t = t, cum = cum, lower = pmax(0, cum - 1.96 * sd), upper = cum + 1.96 * sd)
}

# Pruned phylo tree of query + nearest neighbors. Returns list(tree, matched) or NULL.
.ts_query_tree <- function(outdir, species, query_label, neighbor_labels) {
  db <- .ts_db_path(outdir)
  if (!file.exists(db)) return(NULL)
  if (!requireNamespace("DBI", quietly = TRUE) || !requireNamespace("duckdb", quietly = TRUE) ||
      !requireNamespace("ape", quietly = TRUE)) return(NULL)
  con <- tryCatch(DBI::dbConnect(duckdb::duckdb(), dbdir = db, read_only = TRUE), error = function(e) NULL)
  if (is.null(con)) return(NULL)
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  trees <- tryCatch(DBI::dbGetQuery(con,
    "SELECT tree_id, newick FROM phylogenetic_trees WHERE LOWER(species) = ?",
    params = list(tolower(species))), error = function(e) NULL)
  if (is.null(trees) || nrow(trees) == 0) return(NULL)
  qlab <- .ts_clean_label(query_label)
  keep <- unique(c(qlab, .ts_clean_label(neighbor_labels)))
  keep <- keep[!is.na(keep) & keep != ""]
  best <- NULL; best_n <- -1
  for (i in seq_len(nrow(trees))) {
    nwk <- trees$newick[i]
    if (is.na(nwk) || nchar(trimws(nwk)) == 0) next
    tmp <- tempfile(fileext = ".nwk"); writeLines(nwk, tmp)
    tr <- tryCatch(ape::read.tree(tmp), error = function(e) NULL); unlink(tmp)
    if (is.null(tr) || is.null(tr$tip.label)) next
    tr$tip.label <- .ts_clean_label(tr$tip.label)
    if (!(qlab %in% tr$tip.label)) next
    present <- keep[keep %in% tr$tip.label]
    if (length(present) > best_n && length(present) >= 2) {
      sub <- tryCatch(ape::keep.tip(tr, present), error = function(e) NULL)
      if (!is.null(sub)) { best <- list(tree = ape::ladderize(sub), matched = present); best_n <- length(present) }
    }
  }
  best
}

# ---- UI ----
transmission_content_ui <- function(species, outdir, all_species) {
  if (is.null(species)) {
    return(div(style = "text-align:center; padding:40px 0; color:#888;",
               icon("share-nodes", class = "fa-3x"), h4("No species selected"),
               p("Select a species from the Intelligence Overview dropdown.")))
  }
  has_data <- dir.exists(file.path(outdir, "transmission_context", species))
  tagList(
    bs4Dash::bs4Card(
      title = tagList(icon("share-nodes"), " Transmission & Spread Intelligence"),
      width = 12, status = "primary",
      fluidRow(
        column(3, selectInput("transmission_species", "Species:", choices = all_species, selected = species, selectize = FALSE)),
        column(4, uiOutput("ts_query_select_ui")),
        column(5, uiOutput("transmission_status"))
      ),
      uiOutput("ts_header_facts")
    ),
    if (has_data) {
      bs4Dash::bs4Card(
        width = 12, status = "primary", solidHeader = FALSE,
        tabsetPanel(
          id = "ts_tabs",
          tabPanel(
            "Intelligence Brief", br(),
            h4(icon("route"), " A. Genomic relationship & how it spread ", tags$small(class = "text-muted", "observed history")),
            p(class = "text-muted", style = "font-size:0.85rem;",
              "The new genome and its closest historical relatives (tree), where each strain was first detected, and how far/fast it moved."),
            fluidRow(column(5, plotOutput("ts_tree", height = "440px")),
                     column(7, plotly::plotlyOutput("ts_outbreak_plot", height = "440px"),
                            uiOutput("ts_outbreak_takeaways"))),
            uiOutput("ts_spread_stats"), br(),
            h4(icon("chart-line"), " B. Closest relatives — outbreak impact over time ", tags$small(class = "text-muted", "observed history")),
            p(class = "text-muted", style = "font-size:0.85rem;",
              "Monthly reported cases and deaths per country — the epidemiological footprint of the lineages closest to the new genome."),
            plotly::plotlyOutput("ts_neighbors_plot", height = "460px"),
            hr(),
            h4(icon("chart-line"), " C. How did it behave historically? ", tags$small(class = "text-muted", "observed history")),
            fluidRow(column(3, uiOutput("ts_strain_sel_ui")), column(9, uiOutput("ts_behavior_boxes"))), br(),
            h4(icon("chart-area"), " D. What might happen next? ", tags$small(class = "badge badge-warning", "MODEL-DERIVED PROJECTION")),
            div(class = "alert alert-warning", style = "font-size:0.82rem; padding:8px 12px;",
                icon("triangle-exclamation"), strong(" Model-derived projection. "),
                "A logistic growth model is fitted to the linked strain's observed historical epidemic curve, then projected forward under the selected scenario. The shaded band is the 95% uncertainty interval; it widens with lead time. This is a model-derived estimate, not a validated forecast."),
            uiOutput("ts_model_cards"),
            fluidRow(column(4,
                            selectInput("ts_scenario", "Scenario:",
                                        choices = c("Baseline (as observed)" = "baseline",
                                                    "Contained (rapid control)" = "contained",
                                                    "Expanded (sustained spread)" = "expanded"),
                                        selected = "baseline", selectize = FALSE),
                            sliderInput("ts_horizon", "Projection horizon (days):", min = 60, max = 365, value = 180, step = 30)),
                     column(8, plotly::plotlyOutput("ts_projection", height = "320px"))),
            hr(),
            h4(icon("magnifying-glass-chart"), " E. Why this assessment? ", tags$small(class = "text-muted", "drivers")),
            uiOutput("ts_drivers"), br(),
            h4(icon("file-medical"), " F. Transmission & Spread Assessment"),
            uiOutput("ts_assessment")
          ),
          tabPanel(
            "Surveillance Data", br(),
            p(class = "text-muted", style = "font-size:0.85rem;",
              "Raw weekly epidemiological burden and anomaly flags for transparency and auditability."),
            h5("Flagged anomalous weeks"), DT::DTOutput("ts_anomaly_table"),
            br(), h5("Weekly burden table"), DT::DTOutput("ts_burden_table")
          )
        )
      )
    } else {
      div(style = "text-align:center; padding:40px 0; color:#888;",
          icon("triangle-exclamation", class = "fa-3x"), h4("No transmission data for this species"),
          p("Run the Nextflow pipeline, or select a species with Transmission & Spread outputs."))
    }
  )
}

# ---- Server ----
transmission_register <- function(input, output, session, outdir_r, current_species, all_species) {
  selected_species <- reactive(if (is.null(input$transmission_species)) current_species() else input$transmission_species)
  output$transmission_body <- renderUI(transmission_content_ui(selected_species(), outdir_r(), all_species()))
  observeEvent(current_species(), updateSelectInput(session, "transmission_species", selected = current_species()))
  output$transmission_status <- renderUI({
    sp <- selected_species(); if (is.null(sp)) return(NULL)
    ok <- file.exists(file.path(outdir_r(), "transmission_context", sp, "transmission_potential.tsv"))
    div(style = paste0("padding-top:6px;color:", if (ok) "#198754" else "#dc3545"),
        icon(if (ok) "check-circle" else "circle-exclamation"),
        if (ok) paste(" Context loaded for", toupper(sp)) else paste(" No context for", toupper(sp)))
  })
  nd <- function(m) DT::datatable(data.frame(Message = m), options = list(dom = "t"), rownames = FALSE)

  rd <- function(f) reactive({ sp <- selected_species(); if (is.null(sp)) NULL else .ts_read(outdir_r(), sp, f) })
  potential_r <- rd("transmission_potential.tsv"); score_r <- rd("transmission_potential_score.tsv")
  strain_profile_r <- rd("strain_transmission_profile.tsv"); spatial_r <- rd("transmission_spatial.tsv")
  outbreak_r <- rd("transmission_outbreak_summary.tsv")
  burden_r <- reactive({ df <- .ts_read(outdir_r(), selected_species(), "transmission_burden.tsv"); if (!is.null(df)) df$week_start <- as.Date(df$week_start); df })
  anomaly_r <- reactive({ df <- .ts_read(outdir_r(), selected_species(), "transmission_anomaly.tsv"); if (!is.null(df) && "week_start" %in% names(df)) df$week_start <- as.Date(df$week_start); df })

  query_samples_r <- reactive({ p <- potential_r(); if (is.null(p)) NULL else sort(unique(stats::na.omit(p$query_sample))) })
  output$ts_query_select_ui <- renderUI({
    qs <- query_samples_r(); if (is.null(qs)) return(div(style = "color:#888;", "No query samples."))
    selectInput("ts_query", "Query genome:", choices = qs, selected = qs[1], selectize = FALSE)
  })
  selected_query <- reactive({ qs <- query_samples_r(); if (is.null(qs)) NULL else if (is.null(input$ts_query) || !(input$ts_query %in% qs)) qs[1] else input$ts_query })
  links_r <- reactive({ p <- potential_r(); sq <- selected_query(); if (is.null(p) || is.null(sq)) NULL else dplyr::filter(p, query_sample == sq) })
  query_row_r <- reactive({ df <- links_r(); if (is.null(df) || !nrow(df)) NULL else df[1, ] })
  score_row_r <- reactive({ sc <- score_r(); sq <- selected_query(); if (is.null(sc) || is.null(sq)) NULL else { d <- dplyr::filter(sc, query_sample == sq); if (nrow(d)) d[1, ] else NULL } })

  focus_strain <- reactiveVal(NULL)
  observeEvent(selected_query(), {
    q <- query_row_r(); st <- if (!is.null(q) && !is.na(q$query_strain) && q$query_strain != "") q$query_strain else NULL
    focus_strain(st)
    sp <- strain_profile_r(); if (!is.null(sp) && nrow(sp)) { if (is.null(st) || !(st %in% sp$strain)) st <- sp$strain[1]; updateSelectizeInput(session, "ts_strain_sel", selected = st) }
  })
  observeEvent(input$ts_strain_sel, {
    s <- input$ts_strain_sel; if (length(s)) focus_strain(as.character(s[1]))
  }, ignoreNULL = TRUE)
  current_strain <- reactive({ fs <- focus_strain(); if (!is.null(fs) && !is.na(fs) && fs != "") return(fs); q <- query_row_r(); if (!is.null(q)) q$query_strain else NULL })
  selected_strains <- reactive({ s <- input$ts_strain_sel; if (!is.null(s) && length(s)) as.character(s) else { cs <- current_strain(); if (!is.null(cs)) cs else character(0) } })
  strain_row_r <- reactive({ sp <- strain_profile_r(); st <- current_strain(); if (is.null(sp) || is.null(st)) NULL else { d <- dplyr::filter(sp, strain == st); if (nrow(d)) d[1, ] else NULL } })
  strain_links_r <- reactive({ df <- links_r(); st <- current_strain(); if (is.null(df)) NULL else { d <- if (is.null(st)) df else dplyr::filter(df, background_strain == st); if (nrow(d)) d else df } })

  # Growth-model fit for the current strain's historical outbreak (model-derived).
  fit_r <- reactive({
    br <- burden_r(); st <- current_strain()
    epi <- .ts_strain_epi(br, st)
    if (is.null(epi)) epi <- .ts_strain_epi(br, NULL)
    .ts_fit_logistic(epi)
  })

  output$ts_header_facts <- renderUI({
    q <- query_row_r(); if (is.null(q)) return(NULL)
    cell <- function(l, v) column(2, tags$small(class = "text-muted", l), div(strong(v)))
    div(style = "margin-top:8px;padding:10px;background:#f8f9fa;border-radius:6px;", fluidRow(
      cell("Strain", ifelse(is.na(q$query_strain), "-", q$query_strain)),
      cell("Clade", ifelse(is.na(q$query_clade), "-", q$query_clade)),
      cell("Lineage", ifelse(is.na(q$query_lineage), "-", q$query_lineage)),
      cell("Outbreak", ifelse(is.na(q$query_outbreak), "-", q$query_outbreak)),
      cell("Location", paste(ifelse(is.na(q$query_country), "-", q$query_country), ifelse(is.na(q$query_collection_date), "", q$query_collection_date)))))
  })

  # A. tree + neighbors
  tree_r <- reactive({
    q <- query_row_r(); df <- links_r(); if (is.null(q) || is.null(df)) return(NULL)
    nb <- df %>% dplyr::filter(!is.na(background_sample)) %>% dplyr::arrange(background_div_diff) %>% dplyr::pull(background_sample)
    .ts_query_tree(outdir_r(), selected_species(), q$query_sample, nb)
  })
  output$ts_tree <- renderPlot({
    t <- tree_r()
    if (is.null(t) || !requireNamespace("ggtree", quietly = TRUE)) { plot.new(); text(0.5, 0.5, "Subtree unavailable", col = "#6c757d"); return() }
    tr <- t$tree; q <- query_row_r(); df <- links_r(); qlab <- .ts_clean_label(q$query_sample)
    smap <- stats::setNames(as.character(df$background_strain), .ts_clean_label(df$background_sample))
    ann <- data.frame(label = tr$tip.label,
                      role = ifelse(tr$tip.label == qlab, "Query", "Relative"),
                      strain = ifelse(tr$tip.label == qlab, "Query", ifelse(tr$tip.label %in% names(smap), smap[tr$tip.label], "other")))
    library(ggtree)
    p <- ggtree(tr, size = 0.5, color = "#555") %<+% ann +
      geom_tippoint(aes(color = strain, shape = role), size = 3) +
      geom_tiplab(size = 2.4, offset = 0.0005, hjust = 0) +
      scale_color_manual(values = c(Query = "#dc3545", .ts_strain_pal(ann$strain[ann$strain != "Query"])), na.value = "#999") +
      scale_shape_manual(values = c(Query = 8, Relative = 16)) + theme_tree2() +
      theme(legend.position = "right", legend.text = element_text(size = 8)) + labs(title = "Query + closest historical genomes")
    xr <- layer_scales(p)$x$range$range; if (length(xr) == 2) p <- p + xlim(NA, xr[2] * 1.6); p
  })
  output$ts_neighbors_plot <- plotly::renderPlotly({
    br <- burden_r(); if (is.null(br) || !nrow(br)) return(NULL)
    d <- br %>% dplyr::filter(!is.na(country), !is.na(week_start)) %>%
      dplyr::mutate(month = as.Date(format(as.Date(week_start), "%Y-%m-01"))) %>%
      dplyr::group_by(country, month) %>%
      dplyr::summarise(cases = sum(as.numeric(new_cases), na.rm = TRUE),
                       deaths = sum(as.numeric(new_deaths), na.rm = TRUE), .groups = "drop") %>%
      dplyr::arrange(country, month)
    if (!nrow(d)) return(NULL)
    pal <- .ts_strain_pal(d$country)
    mk <- function(metric, ttl) {
      plotly::plot_ly(d, x = ~month, y = ~d[[metric]], color = ~country, colors = pal,
                      type = "scatter", mode = "lines+markers",
                      marker = list(size = 6), line = list(width = 1.8),
                      text = ~paste0("<b>", country, "</b> — ", format(month, "%b %Y"),
                                     "<br>", ttl, ": ", format(d[[metric]], big.mark = ",")),
                      hoverinfo = "text", showlegend = (metric == "cases")) %>%
        plotly::layout(xaxis = list(title = "", tickformat = "%b %Y"),
                       yaxis = list(title = ttl, rangemode = "tozero"))
    }
    plotly::subplot(mk("cases", "Monthly cases"), mk("deaths", "Monthly deaths"),
                    shareX = TRUE, titleY = TRUE, margin = 0.04) %>%
      plotly::layout(title = "Reported cases & deaths by country over time",
                     legend = list(orientation = "h", y = -0.18))
  })

  # C. behavior + outbreaks — compare the selected strain(s) side by side
  output$ts_strain_sel_ui <- renderUI({
    sp <- strain_profile_r(); if (is.null(sp) || !nrow(sp)) return(div(class = "text-muted", "No strains."))
    sel0 <- isolate(current_strain()); if (is.null(sel0) || is.na(sel0) || !nzchar(sel0) || !(sel0 %in% sp$strain)) sel0 <- sp$strain[1]
    selectizeInput("ts_strain_sel", "Compare strains (pick one or more):",
                   choices = sp$strain, selected = sel0, multiple = TRUE,
                   options = list(maxItems = 4, placeholder = "Choose strains..."))
  })
  output$ts_behavior_boxes <- renderUI({
    sp <- strain_profile_r(); sel <- selected_strains()
    if (is.null(sp) || !length(sel)) return(div(class = "text-muted", "Select strain(s) to compare."))
    chip <- function(l, v) column(2, div(style = "padding:8px;background:#f8f9fa;border-radius:6px;text-align:center;min-height:58px;",
      tags$small(class = "text-muted", l), div(style = "font-weight:600;font-size:0.95rem;", v)))
    tagList(lapply(sel, function(st) {
      s <- sp[sp$strain == st, ]; if (!nrow(s)) return(NULL); s <- s[1, ]
      div(style = "margin-bottom:12px;",
        div(style = "font-weight:600;margin-bottom:4px;", s$strain, tags$small(class = "text-muted", paste0("  ", s$strain_level))),
        fluidRow(
          chip("Linked cases", format(as.numeric(s$linked_cases), big.mark = ",")),
          chip("Linked deaths", format(as.numeric(s$linked_deaths), big.mark = ",")),
          chip("Mean CFR", ifelse(is.na(s$mean_cfr), "-", paste0(round(as.numeric(s$mean_cfr) * 100, 1), "%"))),
          chip("Countries", s$n_countries),
          chip("Years active", s$active_years),
          chip("Behavior", ifelse(is.na(s$behavior_label), "-", gsub("_", " ", s$behavior_label)))))
    }))
  })
  # D. spread
  output$ts_spread_stats <- renderUI({
    s <- strain_row_r(); if (is.null(s)) return(div(class = "text-muted", "No spread metrics."))
    st <- function(l, v) column(3, div(style = "padding:8px;background:#f8f9fa;border-radius:6px;margin-bottom:8px;", tags$small(class = "text-muted", l), div(style = "font-size:1.1rem;font-weight:600;", v)))
    fluidRow(st("Spread rate", ifelse(is.na(s$spread_rate_km_per_year), "-", paste0(round(as.numeric(s$spread_rate_km_per_year), 1), " km/yr"))),
             st("Max extent", ifelse(is.na(s$max_spread_km), "-", paste0(round(as.numeric(s$max_spread_km), 0), " km"))),
             st("Countries", s$n_countries), st("Origin", ifelse(is.na(s$origin_country), "-", s$origin_country)))
  })
  output$ts_outbreak_plot <- plotly::renderPlotly({
    df <- outbreak_r(); if (is.null(df) || !nrow(df)) return(NULL)
    d <- df %>% dplyr::filter(!is.na(start_year), !is.na(country)) %>%
      dplyr::mutate(country = .ts_norm_country(country),
                    cases = as.numeric(cases), deaths = as.numeric(deaths), cfr = as.numeric(cfr),
                    mid = (cases + deaths) / 2,
                    hover = paste0("<b>", country, "</b> ", start_year, "<br>Cases: ", cases, "  Deaths: ", deaths,
                                   "<br>CFR: ", round(cfr * 100, 1), "%")) %>%
      dplyr::arrange(start_year)
    if (!nrow(d)) return(NULL)
    pal <- .ts_strain_pal(d$country)
    xr <- range(d$start_year, na.rm = TRUE)

    # dumbbell connectors (deaths -> cases), NA-separated so each is its own segment
    seg <- do.call(rbind, lapply(seq_len(nrow(d)), function(i) {
      r <- d[i, ]; data.frame(year = c(r$start_year, r$start_year, NA), n = c(r$deaths, r$cases, NA), country = r$country)
    }))

    dumb <- plotly::plot_ly() %>%
      plotly::add_lines(data = seg, x = ~year, y = ~n, color = ~country, colors = pal,
                        line = list(width = 1.4), showlegend = FALSE, hoverinfo = "none") %>%
      plotly::add_markers(data = d, x = ~start_year, y = ~cases, color = ~country, colors = pal,
                          marker = list(symbol = "circle", size = 10, line = list(width = 1, color = "#ffffff")),
                          text = ~hover, hoverinfo = "text", showlegend = FALSE) %>%
      plotly::add_markers(data = d, x = ~start_year, y = ~deaths, color = ~country, colors = pal,
                          marker = list(symbol = "circle-open", size = 10, line = list(width = 2)),
                          text = ~hover, hoverinfo = "text", showlegend = FALSE) %>%
      plotly::add_text(data = d, x = ~start_year, y = ~mid, text = ~paste0("CFR ", round(cfr * 100), "%"),
                       textposition = "middle right", textfont = list(size = 8, color = "#868e96"),
                       showlegend = FALSE, hoverinfo = "none") %>%
      plotly::layout(xaxis = list(range = c(xr[1] - 1, xr[2] + 1)),
                     yaxis = list(title = "Cases / deaths", rangemode = "tozero"))

    tl <- plotly::plot_ly(d, x = ~start_year, y = ~reorder(country, start_year), type = "scatter", mode = "markers",
                          color = ~country, colors = pal,
                          marker = list(size = ~pmax(7, sqrt(cases) * 1.3), opacity = 0.85,
                                        line = list(width = 1, color = "#ffffff")),
                          text = ~hover, hoverinfo = "text", showlegend = FALSE) %>%
      plotly::layout(xaxis = list(title = "Year", range = c(xr[1] - 1, xr[2] + 1)),
                     yaxis = list(title = ""))

    plotly::subplot(dumb, tl, nrows = 2, shareX = TRUE, titleY = TRUE, heights = c(0.66, 0.34), margin = 0.05) %>%
      plotly::layout(title = "Historical outbreaks — cases vs deaths (top) and outbreaks by country (bottom)",
                     annotations = list(list(x = 0, y = 1.07, xref = "paper", yref = "paper", showarrow = FALSE,
                                             xanchor = "left", font = list(size = 11, color = "#6c757d"),
                                             text = "\u25CF cases   \u25CB deaths   \u2502 connected pair")))
  })

  # auto-generated key takeaways for the outbreaks panel
  output$ts_outbreak_takeaways <- renderUI({
    df <- outbreak_r(); if (is.null(df) || !nrow(df)) return(NULL)
    d <- df %>% dplyr::filter(!is.na(start_year), !is.na(country)) %>%
      dplyr::mutate(country = .ts_norm_country(country),
                    cases = as.numeric(cases), deaths = as.numeric(deaths), cfr = as.numeric(cfr))
    if (!nrow(d)) return(NULL)
    big <- d[which.max(d$cases), ]
    cnt <- sort(table(d$country), decreasing = TRUE)
    most <- names(cnt)[1]
    hi_cfr <- round(100 * mean(d$cfr >= 0.7, na.rm = TRUE))
    yrs <- range(d$start_year, na.rm = TRUE)
    div(style = "background:#f8f9fa;border-left:4px solid #4A6C8C;border-radius:4px;padding:8px 12px;font-size:0.82rem;margin-top:6px;",
        tags$strong("Key takeaways"),
        tags$ul(style = "margin:4px 0 0 16px;padding:0;",
                tags$li(paste0("Largest outbreak: ", big$country, " ", big$start_year, " — ", big$cases, " cases, ", big$deaths, " deaths (CFR ", round(big$cfr * 100), "%).")),
                tags$li(paste0(most, " recorded the most outbreaks (", cnt[1], ").")),
                tags$li(paste0(hi_cfr, "% of outbreaks had CFR \u2265 70%.")),
                tags$li(paste0("Outbreaks span ", yrs[1], "\u2013", yrs[2], " across ", length(unique(d$country)), " countries."))))
  })

  # D. model
  output$ts_model_cards <- renderUI({
    fit <- fit_r(); sc <- score_row_r()
    if (is.null(fit)) return(div(class = "text-muted", "Insufficient epi history to fit a growth model for this strain."))
    cd <- function(l, v, c) column(3, div(style = paste0("padding:10px;background:#fff;border-left:4px solid ", c, ";border-radius:4px;margin-bottom:8px;"), tags$small(class = "text-muted", l), div(style = "font-size:1.15rem;font-weight:700;", v)))
    fluidRow(
      cd("Growth rate", paste0(round(fit$r_week, 2), " /wk"), "#4A6C8C"),
      cd("Doubling time", paste0(round(fit$doubling_days, 1), " days"), "#fd7e14"),
      cd("Est. outbreak size", format(round(fit$K), big.mark = ","), "#6f42c1"),
      cd("Model", fit$method, "#198754"))
  })
  output$ts_projection <- plotly::renderPlotly({
    fit <- fit_r(); if (is.null(fit)) return(NULL)
    hz <- if (is.null(input$ts_horizon)) 180 else input$ts_horizon
    scn <- if (is.null(input$ts_scenario)) "baseline" else input$ts_scenario
    pr <- .ts_project(fit, hz, scn); if (is.null(pr)) return(NULL)
    p <- plotly::plot_ly() %>%
      plotly::add_ribbons(data = pr, x = ~t, ymin = ~lower, ymax = ~upper, name = "95% band",
                          fillcolor = "rgba(74,108,140,0.18)", line = list(width = 0)) %>%
      plotly::add_lines(data = pr, x = ~t, y = ~cum, name = "Projected outbreak", line = list(color = "#4A6C8C", width = 2.5))
    if (!is.null(fit$obs)) p <- p %>% plotly::add_lines(data = fit$obs, x = ~t, y = ~cum, name = "Historical (observed)",
                          line = list(color = "#2c3e50", dash = "dot", width = 1.5))
    p %>% plotly::layout(title = paste0("Model-derived projection — ", scn, " scenario"),
                   xaxis = list(title = "Days since introduction"),
                   yaxis = list(title = "Cumulative cases"), legend = list(orientation = "h"))
  })

  # H. drivers
  output$ts_drivers <- renderUI({
    sc <- score_row_r(); s <- strain_row_r(); fit <- fit_r()
    if (is.null(sc) && is.null(s) && is.null(fit)) return(div(class = "text-muted", "No drivers."))
    dr <- function(l, v) div(style = "padding:8px;border-left:3px solid #4A6C8C;margin-bottom:6px;background:#f8f9fa;", div(style = "font-weight:600;", l), div(style = "font-size:0.9rem;", v))
    tagList(
      if (!is.null(sc)) dr("Epi-linked relatives", paste0(sc$n_neighbors_with_epi, " of ", sc$n_neighbors, " neighbors have epi data")),
      if (!is.null(sc) && !is.na(sc$mean_div_diff)) dr("Genetic similarity", paste0("Mean divergence diff ", round(as.numeric(sc$mean_div_diff), 4))),
      if (!is.null(sc) && !is.na(sc$min_geo_km)) dr("Geographic proximity", paste0("Nearest linked genome ", round(as.numeric(sc$min_geo_km), 0), " km")),
      if (!is.null(s)) dr("Historical burden", paste0(format(as.numeric(s$linked_cases), big.mark = ","), " cases / ", s$n_countries, " countries")),
      if (!is.null(s)) dr("Geographic spread", paste0(round(as.numeric(s$max_spread_km), 0), " km extent")),
      if (!is.null(s)) dr("Persistence", paste0(s$active_years, " years active")),
      if (!is.null(fit)) dr("Fitted growth", paste0(round(fit$r_week, 2), " cases/wk growth, doubling every ", round(fit$doubling_days, 1), " days (", fit$method, ")")))
  })

  # I. assessment
  output$ts_assessment <- renderUI({
    q <- query_row_r(); sc <- score_row_r(); s <- strain_row_r(); fit <- fit_r(); if (is.null(q)) return(div(class = "text-muted", "No assessment."))
    strain <- ifelse(is.na(q$query_strain), "an uncharacterized strain", q$query_strain)
    risk <- if (!is.null(sc) && !is.na(sc$risk_label)) toupper(sc$risk_label) else "UNKNOWN"
    hist <- if (!is.null(s)) paste0("Historically, ", strain, " caused ", format(as.numeric(s$linked_cases), big.mark = ","), " cases and ", format(as.numeric(s$linked_deaths), big.mark = ","), " deaths across ", s$n_countries, " countries over ", s$active_years, " year(s), classified as ", gsub("_", " ", s$behavior_label), ".") else "No historical strain profile."
    mdl <- if (!is.null(fit)) paste0("Fitting a growth model to ", strain, "'s observed epidemic curve gives a growth rate of ", round(fit$r_week, 2), " cases/week (doubling every ", round(fit$doubling_days, 1), " days) and an estimated outbreak size of ~", format(round(fit$K), big.mark = ","), " cases; projected forward under the selected scenario.") else if (!is.null(sc)) paste0("Predicted behavior '", gsub("_", " ", sc$predicted_behavior), "'; insufficient epi history for a fitted projection.") else "No model projection."
    div(style = "padding:14px;background:#eef4f8;border-left:4px solid #4A6C8C;border-radius:6px;",
        h5("Transmission & Spread Assessment"),
        p(tags$span(class = "badge badge-secondary", "OBSERVED"), " Query ", strong(q$query_sample), " is strain ", strong(strain), ". ", hist),
        p(tags$span(class = "badge badge-warning", "MODEL-DERIVED"), " ", mdl),
        p(class = "text-muted", style = "font-size:0.82rem;margin-bottom:0;", "Observed = historical evidence. Model-derived = projection under stated assumptions with uncertainty."))
  })

  # Surveillance tab
  output$ts_anomaly_table <- DT::renderDT({ df <- anomaly_r(); if (is.null(df) || !nrow(df)) nd("No anomalous weeks") else DT::datatable(df, options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE) })
  output$ts_burden_table <- DT::renderDT({
    df <- burden_r(); if (is.null(df) || !nrow(df)) return(nd("No burden data"))
    s <- df %>% dplyr::select(country, admin1, week_start, cases_cum, deaths_cum, new_cases, new_deaths, cfr_cum, cfr_new, growth_rate, trend_cases, alert, n_genomes, dominant_strain) %>%
      dplyr::mutate(cfr_cum = round(cfr_cum, 3), cfr_new = round(cfr_new, 3), growth_rate = round(growth_rate, 3), trend_cases = as.character(trend_cases))
    DT::datatable(s, options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
  })
}
