# Geographic & Temporal Context module
#
# Synchronized genomic-epidemiological reconstruction. A single Play/slider
# time control drives four coordinated views:
#   MAP       - where lineages appeared and spread (genomes reveal by date,
#               chronological spread paths, outbreak circles sized by burden)
#   PHYLO     - how the query genome relates to its closest historical genomes
#   TIMELINE  - when transmission intensified (monthly cases/deaths + cursor)
#   NARRATIVE - a chronological genomic-epidemiological story
# Links are tiered (confirmed / inferred / probable) - never claims direct
# person-to-person transmission. Reuses .ts_* helpers from pathogen_transmission.R.

# ---- Data helpers ----

.gt_decimal_to_date <- function(num) {
  yr <- floor(num); mo <- floor((num - yr) * 12) + 1
  mo <- pmax(1, pmin(12, mo))
  as.Date(paste0(yr, "-", sprintf("%02d", mo), "-15"))
}

.gt_parse_date <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (inherits(x, c("POSIXct", "POSIXlt"))) return(as.Date(x))
  d <- suppressWarnings(as.Date(as.character(x)))
  num <- suppressWarnings(as.numeric(x))
  need <- is.na(d) & !is.na(num)
  if (any(need)) d[need] <- .gt_decimal_to_date(num[need])
  d
}

# All genomes for a species from the DuckDB warehouse (samples + tree_tips + geo).
.gt_samples <- function(outdir, species) {
  db <- .ts_db_path(outdir)
  if (!file.exists(db)) return(NULL)
  if (!requireNamespace("DBI", quietly = TRUE) || !requireNamespace("duckdb", quietly = TRUE)) return(NULL)
  con <- tryCatch(DBI::dbConnect(duckdb::duckdb(), dbdir = db, read_only = TRUE), error = function(e) NULL)
  if (is.null(con)) return(NULL)
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  df <- tryCatch(DBI::dbGetQuery(con, "
    SELECT s.sample_name, s.collection_date, s.country, s.admin1,
           COALESCE(s.clade, t.clade) AS clade, s.lineage, s.outbreak, s.is_query,
           gl.latitude, gl.longitude, t.tip_date, t.div
    FROM samples s
    LEFT JOIN sample_geo_location sgl ON s.sample_id = sgl.sample_id
    LEFT JOIN geographic_locations gl ON sgl.location_id = gl.location_id
    LEFT JOIN tree_tips t ON s.sample_id = t.sample_id
    WHERE LOWER(s.species) = ?",
    params = list(tolower(species))), error = function(e) NULL)
  if (is.null(df) || !nrow(df)) return(NULL)
  df <- df %>% dplyr::distinct(sample_name, .keep_all = TRUE)
  df$is_query <- tolower(as.character(df$is_query)) %in% c("true", "1", "t", "yes")
  df$date <- dplyr::coalesce(.gt_parse_date(df$collection_date), .gt_parse_date(df$tip_date))
  for (cc in c("clade", "lineage", "outbreak")) {
    v <- as.character(df[[cc]]); v[v %in% c("", "NA", "nan", "None", "unassigned")] <- NA; df[[cc]] <- v
  }
  df$strain <- dplyr::coalesce(df$clade, df$lineage, df$outbreak)
  df <- .ts_fill_coords(df, "latitude", "longitude", "country")
  df
}

# Build the chronological narrative event list. Each event: date, type, tier, text(html).
.gt_events <- function(samples, outbreaks, links, qname) {
  ev <- list()
  add <- function(date, type, tier, text) {
    if (is.null(date) || length(date) == 0 || is.na(date) || !is.finite(as.numeric(date))) return(invisible(NULL))
    ev[[length(ev) + 1]] <<- list(date = as.Date(date), type = type, tier = tier, text = text)
  }
  bg <- samples[!samples$is_query & !is.na(samples$date), ]
  bg_s <- if ("strain" %in% names(bg)) bg[!is.na(bg$strain), ] else bg[0, ]
  if (nrow(bg_s)) {
    fs <- bg_s %>% dplyr::group_by(strain) %>%
      dplyr::summarise(date = min(date), n = dplyr::n(),
                       country = dplyr::first(country[order(date)]), .groups = "drop")
    for (i in seq_len(nrow(fs)))
      add(fs$date[i], "origin", "observed",
          paste0("Strain <b>", fs$strain[i], "</b> first sampled in <b>", .ts_norm_country(fs$country[i]),
                 "</b> (", fs$n[i], " genome(s) in total)."))
    for (st in unique(stats::na.omit(bg_s$strain))) {
      g <- bg[bg$strain == st & !is.na(bg$country), ]; g <- g[order(g$date), ]
      seen <- character(0)
      for (j in seq_len(nrow(g))) {
        ctry <- g$country[j]
        if (!is.na(ctry) && !(ctry %in% seen)) {
          if (length(seen) > 0)
            add(g$date[j], "spread", "inferred",
                paste0("<b>", st, "</b> detected in <b>", .ts_norm_country(ctry),
                       "</b> — first record outside ", paste(.ts_norm_country(seen), collapse = ", "), "."))
          seen <- c(seen, ctry)
        }
      }
    }
  }
  if (!is.null(outbreaks) && nrow(outbreaks)) {
    ob <- outbreaks[!is.na(outbreaks$start_year) & !is.na(outbreaks$country), ]
    for (i in seq_len(nrow(ob))) {
      cfr_txt <- if (!is.na(ob$cfr[i])) paste0(" (CFR ", round(as.numeric(ob$cfr[i]) * 100), "%)") else ""
      add(as.Date(paste0(ob$start_year[i], "-06-01")), "outbreak", "contextual",
          paste0("Outbreak in <b>", .ts_norm_country(ob$country[i]), "</b> (", ob$start_year[i], "): ",
                 format(as.numeric(ob$cases[i]), big.mark = ","), " cases, ",
                 format(as.numeric(ob$deaths[i]), big.mark = ","), " deaths", cfr_txt, "."))
    }
  }
  q <- samples[samples$is_query, ]
  if (!is.null(qname) && nrow(q)) q <- q[q$sample_name == qname, ]
  if (nrow(q)) {
    qd <- q$date[1]; if (is.na(qd)) qd <- Sys.Date()
    add(qd, "query", "observed",
        paste0("Query genome <b>", q$sample_name[1], "</b> collected in <b>", .ts_norm_country(q$country[1]), "</b>."))
  }
  if (!is.null(links) && nrow(links)) {
    lk <- links %>% dplyr::filter(!is.na(background_sample))
    # Order by genetic proximity when available, else temporal proximity.
    if ("background_div_diff" %in% names(lk) && any(!is.na(lk$background_div_diff)))
      lk <- lk %>% dplyr::arrange(background_div_diff)
    else if ("temporal_distance_days" %in% names(lk))
      lk <- lk %>% dplyr::arrange(temporal_distance_days)
    lk <- lk %>% dplyr::distinct(background_sample, .keep_all = TRUE)
    if (nrow(lk)) {
      for (i in seq_len(min(5, nrow(lk)))) {
        r <- lk[i, ]
        same_c <- !is.na(r$background_country) && !is.na(r$query_country) &&
          identical(.ts_norm_country(r$background_country), .ts_norm_country(r$query_country))
        same_s <- !is.na(r$background_strain) && !is.na(r$query_strain) &&
          identical(as.character(r$background_strain), as.character(r$query_strain))
        td <- suppressWarnings(as.numeric(r$temporal_distance_days))
        # Tiers use strain/country/temporal signals (div_diff is NA in this data).
        tier <- if (same_s && same_c && (is.na(td) || td <= 730)) "confirmed"
                else if (same_s || (same_c && (is.na(td) || td <= 1825))) "inferred"
                else "probable"
        bd <- .gt_parse_date(r$background_collection_date); if (is.na(bd)) bd <- Sys.Date()
        dd <- if (!is.na(r$background_div_diff)) paste0("div diff ", signif(r$background_div_diff, 3))
              else if (!is.na(td)) paste0(round(td / 30.4), " mo apart") else "related"
        add(bd, "link", tier,
            paste0("<b>", r$background_sample, "</b> sampled in <b>", .ts_norm_country(r$background_country),
                   "</b> — ", tier, " link to query (", dd, ")."))
      }
    }
  }
  if (length(ev)) ev <- ev[order(vapply(ev, function(e) as.numeric(e$date), numeric(1)))]
  ev
}

# Great-circle distance in km (vectorised).
.gt_haversine <- function(lon1, lat1, lon2, lat2) {
  r <- 6371; tr <- pi / 180
  dlat <- (lat2 - lat1) * tr; dlon <- (lon2 - lon1) * tr
  a <- sin(dlat / 2)^2 + cos(lat1 * tr) * cos(lat2 * tr) * sin(dlon / 2)^2
  2 * r * asin(pmin(1, sqrt(a)))
}

# Build a plausible transmission network. Each genome's "parent" is the earlier
# genome (within a 5-year lookback) minimising a combined score of genetic
# divergence difference (tree_tips.div, when present), great-circle distance and
# temporal gap, with bonuses for shared country/strain. This is a heuristic
# reconstruction of "who plausibly seeded whom" — NOT a proven person-to-person
# chain. Returns one row per genome (its incoming edge), revealed at the child date.
.gt_edges <- function(samples) {
  s <- samples[!is.na(samples$date) & !is.na(samples$latitude) & !is.na(samples$longitude), ]
  if (nrow(s) < 2) return(NULL)
  s <- s[order(s$date), ]
  n <- nrow(s); has_div <- any(!is.na(s$div))
  out <- vector("list", n)
  for (i in 2:n) {
    cand <- which(s$date[seq_len(i - 1)] >= s$date[i] - 1825)  # earlier, within ~5y
    if (!length(cand)) cand <- i - 1                            # dormant lineage: immediate predecessor
    geo  <- .gt_haversine(s$longitude[i], s$latitude[i], s$longitude[cand], s$latitude[cand])
    tgap <- as.numeric(s$date[i] - s$date[cand])
    score <- geo / 1000 + tgap / 365
    if (has_div) {
      dd <- abs(s$div[i] - s$div[cand])
      med <- stats::median(dd, na.rm = TRUE); if (is.na(med)) med <- 0
      dd[is.na(dd)] <- med
      score <- score + dd * 50
    }
    same_c <- !is.na(s$country[i]) & !is.na(s$country[cand]) & s$country[cand] == s$country[i]
    same_s <- !is.na(s$strain[i])  & !is.na(s$strain[cand])  & s$strain[cand]  == s$strain[i]
    score <- score - same_c * 0.5 - same_s * 0.5
    score[!is.finite(score)] <- Inf
    pos <- which.min(score); if (!length(pos) || is.na(pos)) pos <- 1
    j <- cand[pos]
    tier <- if (isTRUE(same_s[pos]) && isTRUE(same_c[pos]) && tgap[pos] <= 730) "confirmed"
            else if (isTRUE(same_s[pos]) || isTRUE(same_c[pos])) "inferred" else "probable"
    out[[i]] <- data.frame(
      from = s$sample_name[j], to = s$sample_name[i],
      from_lon = s$longitude[j], from_lat = s$latitude[j],
      to_lon = s$longitude[i], to_lat = s$latitude[i],
      date = s$date[i], tier = tier,
      geo_km = round(geo[pos]), tgap_days = round(tgap[pos]),
      div_diff = if (has_div) signif(abs(s$div[i] - s$div[j]), 3) else NA_real_,
      from_country = s$country[j], to_country = s$country[i], stringsAsFactors = FALSE)
  }
  do.call(rbind, out)
}

# ---- UI ----
geographic_temporal_ui <- function(sp, all_species = sp) {
  ns <- function(id) paste0(sp, "_", id)
  tagList(
    br(),
    fluidRow(
      column(8, h4(icon("earth-africa"), " Geographic & Temporal Context ",
                   tags$small(class = "text-muted", "genomic-epidemiological reconstruction")),
             p(class = "text-muted", style = "font-size:0.85rem;",
               "Press play or drag the timeline to watch the lineage's history unfold: where genomes appeared, how they spread, and which outbreaks they coincided with.")),
      column(4, selectInput(ns("gt_species_jump"), "Species:", choices = all_species, selected = sp, selectize = FALSE))
    ),
    fluidRow(
      column(3, uiOutput(ns("gt_query_select_ui"))),
      column(6, uiOutput(ns("gt_time_ui"))),
      column(3, uiOutput(ns("gt_now_label")))
    ),
    fluidRow(
      column(12, div(style = "display:flex;gap:22px;align-items:center;flex-wrap:wrap;padding:2px 0;",
        checkboxInput(ns("gt_show_network"), "Transmission network", TRUE, width = "auto"),
        checkboxInput(ns("gt_show_paths"), "Spread paths", TRUE, width = "auto"),
        tags$span(class = "text-muted", style = "font-size:0.78rem;",
                  "Network edges = inferred transmission (genetic divergence + geography + time); line points to the newer genome.")))
    ),
    fluidRow(
      column(7, leaflet::leafletOutput(ns("gt_map"), height = "540px")),
      column(5, plotOutput(ns("gt_tree"), height = "540px"))
    ),
    fluidRow(
      column(8, plotly::plotlyOutput(ns("gt_timeline"), height = "200px")),
      column(4, plotly::plotlyOutput(ns("gt_pie"), height = "200px"))
    ),
    br(),
    fluidRow(column(12, uiOutput(ns("gt_legend")))),
    fluidRow(column(12, h5(icon("clock"), " Reconstruction narrative"), uiOutput(ns("gt_narrative"))))
  )
}

geographic_temporal_register <- function(input, output, session, sp, outdir_r) {
  ns <- function(id) paste0(sp, "_", id)

  # corner species dropdown -> jump to that species' geographic tab
  observeEvent(input[[ns("gt_species_jump")]], {
    tgt <- input[[ns("gt_species_jump")]]
    if (!is.null(tgt) && tgt != sp) updateTabItems(session, "sidebarmenu", selected = paste0("gt_", tgt))
  }, ignoreInit = TRUE)

  potential_r <- reactive({ if (is.null(sp)) NULL else .ts_read(outdir_r(), sp, "transmission_potential.tsv") })
  burden_r <- reactive({ df <- .ts_read(outdir_r(), sp, "transmission_burden.tsv"); if (!is.null(df) && "week_start" %in% names(df)) df$week_start <- as.Date(df$week_start); df })
  outbreak_r <- reactive({ .ts_read(outdir_r(), sp, "transmission_outbreak_summary.tsv") })

  # All genomes for the species (DuckDB); fall back to potential.tsv neighbors.
  samples_r <- reactive({
    s <- .gt_samples(outdir_r(), sp)
    if (!is.null(s) && nrow(s)) return(s)
    p <- potential_r(); if (is.null(p)) return(NULL)
    bg <- p %>% dplyr::filter(!is.na(background_sample)) %>% dplyr::transmute(
      sample_name = background_sample, date = .gt_parse_date(background_collection_date),
      country = background_country, admin1 = background_admin1, clade = background_clade,
      lineage = background_lineage, outbreak = NA_character_, strain = background_strain,
      is_query = FALSE, latitude = background_latitude, longitude = background_longitude, div = background_div)
    q <- p %>% dplyr::filter(!is.na(query_sample)) %>% dplyr::transmute(
      sample_name = query_sample, date = .gt_parse_date(query_collection_date),
      country = query_country, admin1 = query_admin1, clade = query_clade,
      lineage = query_lineage, outbreak = query_outbreak, strain = query_strain,
      is_query = TRUE, latitude = query_latitude, longitude = query_longitude, div = query_div)
    s <- dplyr::bind_rows(bg, q) %>% dplyr::distinct(sample_name, .keep_all = TRUE)
    .ts_fill_coords(s, "latitude", "longitude", "country")
  })

  query_samples_r <- reactive({ p <- potential_r(); if (is.null(p)) NULL else sort(unique(stats::na.omit(p$query_sample))) })
  output[[ns("gt_query_select_ui")]] <- renderUI({
    qs <- query_samples_r(); if (is.null(qs)) return(div(style = "color:#888;", "No query samples."))
    selectInput(ns("gt_query"), "Query genome:", choices = qs, selected = qs[1], selectize = FALSE)
  })
  selected_query <- reactive({ qs <- query_samples_r(); if (is.null(qs)) NULL else if (is.null(input[[ns("gt_query")]]) || !(input[[ns("gt_query")]] %in% qs)) qs[1] else input[[ns("gt_query")]] })
  links_r <- reactive({ p <- potential_r(); sq <- selected_query(); if (is.null(p) || is.null(sq)) NULL else dplyr::filter(p, query_sample == sq) })
  query_row_r <- reactive({ df <- links_r(); if (is.null(df) || !nrow(df)) NULL else df[1, ] })

  # ---- Time control ----
  output[[ns("gt_time_ui")]] <- renderUI({
    s <- samples_r(); if (is.null(s)) return(NULL)
    d <- s$date[!is.na(s$date)]; if (!length(d)) return(NULL)
    mn <- min(d); mx <- max(d); if (mx <= mn) mx <- mn + 1
    step <- max(10, round(as.numeric(mx - mn) / 90))
    tagList(
      tags$small(class = "text-muted", "Reconstruction timeline — press play"),
      sliderInput(ns("gt_time"), NULL, min = mn, max = mx, value = mn, step = step,
                  timeFormat = "%Y-%m", width = "100%",
                  animate = animationOptions(interval = 140, loop = FALSE))
    )
  })
  current_time <- reactive({
    t <- input[[ns("gt_time")]]
    if (!is.null(t)) return(as.Date(t))
    s <- samples_r(); if (is.null(s)) return(Sys.Date())
    d <- s$date[!is.na(s$date)]; if (!length(d)) Sys.Date() else min(d)
  })
  output[[ns("gt_now_label")]] <- renderUI({
    ct <- current_time()
    div(style = "padding-top:22px;text-align:right;",
        tags$small(class = "text-muted", "Reconstruction date"),
        div(style = "font-size:1.25rem;font-weight:700;", format(ct, "%b %Y")))
  })

  # Narrative events (precomputed, sorted by date)
  events_r <- reactive({
    s <- samples_r(); if (is.null(s)) return(list())
    .gt_events(s, outbreak_r(), links_r(), selected_query())
  })

  # Transmission network edges (computed once per species; cached).
  edges_r <- reactive({ s <- samples_r(); if (is.null(s)) NULL else .gt_edges(s) })

  # ---- Map ----
  # Stable strain palette across the whole species so colors don't shift as time advances.
  strain_pal <- reactive({ s <- samples_r(); if (is.null(s)) NULL else leaflet::colorFactor(.ts_strain_pal(s$strain), domain = s$strain) })

  output[[ns("gt_map")]] <- leaflet::renderLeaflet({
    m <- leaflet::leaflet() %>% leaflet::addTiles()
    isolate({
      s <- samples_r()
      if (!is.null(s) && nrow(s)) {
        bb <- s[!is.na(s$latitude) & !is.na(s$longitude), ]
        if (nrow(bb)) m <- m %>% leaflet::fitBounds(min(bb$longitude), min(bb$latitude), max(bb$longitude), max(bb$latitude))
      }
    })
    m
  })

  # Incremental updates: genomes carry layerId = sample_name and network edges
  # carry layerId = child sample_name, so each tick only adds newly-revealed and
  # removes newly-hidden layers (fast even for ebov).
  tier_col <- c(confirmed = "#198754", inferred = "#fd7e14", probable = "#6c757d")
  shown_genomes <- reactiveVal(character(0))
  shown_edges   <- reactiveVal(character(0))
  last_bb       <- reactiveVal(NULL)
  observe({
    ct <- current_time(); s <- samples_r(); ob <- outbreak_r(); sq <- selected_query()
    ed <- edges_r()
    show_net   <- isTRUE(input[[ns("gt_show_network")]])
    show_paths <- isTRUE(input[[ns("gt_show_paths")]])
    proxy <- leaflet::leafletProxy(ns("gt_map"))
    if (is.null(s) || !nrow(s)) { shown_genomes(character(0)); shown_edges(character(0)); return() }
    pal <- strain_pal()

    rev <- s[!s$is_query & !is.na(s$date) & s$date <= ct & !is.na(s$latitude) & !is.na(s$longitude), ]
    cur <- shown_genomes(); tgt <- rev$sample_name
    to_del <- setdiff(cur, tgt); to_add <- setdiff(tgt, cur)
    if (length(to_del)) proxy <- proxy %>% leaflet::removeMarker(layerId = to_del)
    if (length(to_add) && !is.null(pal)) {
      ad <- rev[rev$sample_name %in% to_add, ]
      proxy <- proxy %>% leaflet::addCircleMarkers(data = ad, ~longitude, ~latitude, layerId = ~sample_name,
        radius = 5, color = ~pal(strain), fillOpacity = 0.75, stroke = FALSE, group = "genomes",
        popup = ~paste0("<b>", sample_name, "</b><br>", .ts_norm_country(country),
                        ifelse(is.na(admin1) | admin1 == "National", "", paste0(", ", admin1)),
                        "<br>", format(date, "%Y-%m-%d"), "<br>Strain: ", ifelse(is.na(strain), "?", strain)))
    }
    shown_genomes(tgt)

    # --- Transmission network edges (incremental, layerId = child name) ---
    if (show_net && !is.null(ed) && nrow(ed)) {
      etgt <- ed$to[ed$date <= ct]
      ecur <- shown_edges()
      e_del <- setdiff(ecur, etgt); e_add <- setdiff(etgt, ecur)
      if (length(e_del)) proxy <- proxy %>% leaflet::removeShape(layerId = e_del)
      if (length(e_add)) {
        ae <- ed[ed$to %in% e_add, ]
        for (k in seq_len(nrow(ae))) {
          r <- ae[k, ]
          proxy <- proxy %>% leaflet::addPolylines(
            lng = c(r$from_lon, r$to_lon), lat = c(r$from_lat, r$to_lat),
            layerId = r$to, group = "network",
            color = unname(tier_col[r$tier]), weight = 2, opacity = 0.7,
            popup = paste0("<b>", r$from, " \u2192 ", r$to, "</b><br>",
                           .ts_norm_country(r$from_country), " \u2192 ", .ts_norm_country(r$to_country),
                           "<br>", ifelse(is.na(r$div_diff), "", paste0("\u0394div ", r$div_diff, " \u00b7 ")),
                           r$geo_km, " km \u00b7 ", round(r$tgap_days / 30.4), " mo",
                           "<br><i>", r$tier, " transmission (inferred)</i>"))
        }
      }
      shown_edges(etgt)
    } else if (length(shown_edges())) {
      proxy <- proxy %>% leaflet::removeShape(layerId = shown_edges())
      shown_edges(character(0))
    }

    # Small layers (paths / outbreaks / query) redrawn each tick.
    proxy <- proxy %>% leaflet::clearGroup("paths") %>% leaflet::clearGroup("outbreaks") %>% leaflet::clearGroup("query")
    if (show_paths && nrow(rev) && !is.null(pal)) {
      for (st in unique(stats::na.omit(rev$strain))) {
        g <- rev[rev$strain == st, ]; g <- g[order(g$date), ]
        if (nrow(g) >= 2)
          proxy <- proxy %>% leaflet::addPolylines(data = g, ~longitude, ~latitude, color = pal(st),
                                                 weight = 1.5, opacity = 0.55, dashArray = "4", group = "paths")
      }
    }
    if (!is.null(ob) && nrow(ob)) {
      o <- ob[!is.na(ob$start_year) & !is.na(ob$country), ]
      if (nrow(o)) {
        o <- o[o$start_year <= as.numeric(format(ct, "%Y")), ]
        if (nrow(o)) {
          cc <- .ts_country_centroid(o$country); o$lat <- cc$lat; o$lon <- cc$lon
          o <- o[!is.na(o$lat) & !is.na(o$lon), ]
          if (nrow(o))
            proxy <- proxy %>% leaflet::addCircleMarkers(data = o, ~lon, ~lat,
              radius = ~pmax(6, sqrt(as.numeric(cases)) / 6), color = "#b03a2e", weight = 2,
              fillColor = "#e74c3c", fillOpacity = 0.25, group = "outbreaks",
              popup = ~paste0("<b>", .ts_norm_country(country), "</b> ", start_year,
                              "<br>", format(as.numeric(cases), big.mark = ","), " cases / ",
                              format(as.numeric(deaths), big.mark = ","), " deaths"))
        }
      }
    }
    q <- if (!is.null(sq)) s[s$is_query & s$sample_name == sq, ] else s[s$is_query, ]
    if (!nrow(q)) q <- s[s$is_query, ]
    if (nrow(q)) {
      q <- q[1, ]
      if (!is.na(q$date) && q$date <= ct && !is.na(q$latitude) && !is.na(q$longitude))
        proxy <- proxy %>% leaflet::addCircleMarkers(lng = q$longitude, lat = q$latitude, radius = 10,
          color = "#dc3545", weight = 3, fillColor = "#dc3545", fillOpacity = 0.9, group = "query",
          popup = paste0("<b>QUERY: ", q$sample_name, "</b><br>", .ts_norm_country(q$country), "<br>", q$date))
    }

    # --- Follow the action: zoom toward the revealed genomes (country level) ---
    if ((length(to_add) || length(to_del)) && nrow(rev)) {
      bb <- c(min(rev$longitude), min(rev$latitude), max(rev$longitude), max(rev$latitude))
      lb <- last_bb()
      if (is.null(lb) || any(abs(bb - lb) > 0.5)) {
        last_bb(bb)
        proxy <- proxy %>% leaflet::flyToBounds(bb[1], bb[2], bb[3], bb[4],
                                                options = list(maxZoom = 5, animate = TRUE, duration = 0.6))
      }
    }
  })

  # ---- Phylo panel (tips reveal as they were sampled) ----
  tree_r <- reactive({
    q <- query_row_r(); df <- links_r(); if (is.null(q) || is.null(df)) return(NULL)
    nb <- df %>% dplyr::filter(!is.na(background_sample)) %>% dplyr::arrange(background_div_diff) %>% dplyr::pull(background_sample)
    .ts_query_tree(outdir_r(), sp, q$query_sample, nb)
  })
  output[[ns("gt_tree")]] <- renderPlot({
    t <- tree_r(); ct <- current_time()
    if (is.null(t) || !requireNamespace("ggtree", quietly = TRUE)) {
      plot.new(); text(0.5, 0.5, "Subtree unavailable", col = "#6c757d"); return()
    }
    tr <- t$tree; q <- query_row_r(); df <- links_r(); qlab <- .ts_clean_label(q$query_sample)
    smap <- stats::setNames(as.character(df$background_strain), .ts_clean_label(df$background_sample))
    dmap <- stats::setNames(.gt_parse_date(df$background_collection_date), .ts_clean_label(df$background_sample))
    ann <- data.frame(
      label = tr$tip.label,
      role = ifelse(tr$tip.label == qlab, "Query", "Relative"),
      strain = ifelse(tr$tip.label == qlab, "Query", ifelse(tr$tip.label %in% names(smap), smap[tr$tip.label], "other")),
      revealed = ifelse(tr$tip.label == qlab, TRUE,
                        ifelse(tr$tip.label %in% names(dmap), !is.na(dmap[tr$tip.label]) & dmap[tr$tip.label] <= ct, TRUE)))
    library(ggtree)
    p <- ggtree(tr, size = 0.5, color = "#555") %<+% ann +
      geom_tippoint(aes(color = strain, shape = role, alpha = revealed), size = 3) +
      geom_tiplab(size = 2.4, offset = 0.0005, hjust = 0, aes(alpha = revealed)) +
      scale_color_manual(values = c(Query = "#dc3545", .ts_strain_pal(ann$strain[ann$strain != "Query"])), na.value = "#999") +
      scale_shape_manual(values = c(Query = 8, Relative = 16)) +
      scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.22), guide = "none") +
      theme_tree2() + theme(legend.position = "right", legend.text = element_text(size = 8)) +
      labs(title = "Query + closest historical genomes")
    xr <- layer_scales(p)$x$range$range; if (length(xr) == 2) p <- p + xlim(NA, xr[2] * 1.6); p
  })

  # ---- Epi timeline (monthly burden + cursor at current time) ----
  output[[ns("gt_timeline")]] <- plotly::renderPlotly({
    br <- burden_r(); ct <- current_time()
    if (is.null(br) || !nrow(br)) return(NULL)
    d <- br %>% dplyr::filter(!is.na(week_start)) %>%
      dplyr::mutate(month = as.Date(format(week_start, "%Y-%m-01"))) %>%
      dplyr::group_by(month) %>%
      dplyr::summarise(cases = sum(as.numeric(new_cases), na.rm = TRUE),
                       deaths = sum(as.numeric(new_deaths), na.rm = TRUE), .groups = "drop") %>%
      dplyr::arrange(month)
    if (!nrow(d)) return(NULL)
    plotly::plot_ly(d, x = ~month) %>%
      plotly::add_bars(y = ~cases, name = "Cases", marker = list(color = "rgba(74,108,140,0.5)")) %>%
      plotly::add_lines(y = ~deaths, name = "Deaths", line = list(color = "#dc3545", width = 1.6)) %>%
      plotly::layout(
        title = list(text = "Reported burden over time — cursor marks reconstruction date", font = list(size = 12)),
        xaxis = list(title = ""), yaxis = list(title = "Monthly cases / deaths", rangemode = "tozero"),
        legend = list(orientation = "h", y = -0.25),
        shapes = list(list(type = "line", x0 = ct, x1 = ct, y0 = 0, y1 = 1, yref = "paper",
                           line = list(color = "#2c3e50", width = 2, dash = "dot"))))
  })

  # ---- Cumulative burden pie (grows as the reconstruction date advances) ----
  output[[ns("gt_pie")]] <- plotly::renderPlotly({
    br <- burden_r(); ct <- current_time()
    if (is.null(br) || !nrow(br) || !"week_start" %in% names(br)) return(NULL)
    d <- br[!is.na(br$week_start) & br$week_start <= ct, ]
    cases  <- sum(as.numeric(d$new_cases),  na.rm = TRUE)
    deaths <- sum(as.numeric(d$new_deaths), na.rm = TRUE)
    if (!is.finite(cases) || cases <= 0) return(NULL)
    nonfatal <- max(0, cases - deaths)
    cfr <- round(100 * deaths / cases, 1)
    plotly::plot_ly(labels = c("Deaths", "Non-fatal cases"), values = c(deaths, nonfatal),
                    type = "pie", hole = 0.55, sort = FALSE,
                    marker = list(colors = c("#dc3545", "#4A6C8C")),
                    textinfo = "label+percent", hoverinfo = "label+value",
                    insidetextorientation = "radial") %>%
      plotly::layout(
        title = list(text = paste0("Cumulative to ", format(ct, "%b %Y"), "<br>",
                                   format(cases, big.mark = ","), " cases \u00b7 CFR ", cfr, "%"),
                     font = list(size = 12)),
        showlegend = FALSE, margin = list(t = 58, b = 8, l = 8, r = 8))
  })

  # ---- Confidence-tier legend ----
  output[[ns("gt_legend")]] <- renderUI({
    div(style = "font-size:0.8rem;color:#555;padding:4px 0;",
        tags$strong("Link confidence: "),
        tags$span(class = "badge", style = "background:#198754;color:#fff;", "CONFIRMED"), " same strain + same country + recent  ",
        tags$span(class = "badge", style = "background:#fd7e14;color:#fff;", "INFERRED"), " same strain or co-located  ",
        tags$span(class = "badge", style = "background:#6c757d;color:#fff;", "PROBABLE"), " broader geo/genetic proximity  ",
        tags$span(class = "badge", style = "background:#4A6C8C;color:#fff;", "CONTEXTUAL"), " epidemiological record  ",
        tags$span(class = "badge", style = "background:#2c3e50;color:#fff;", "OBSERVED"), " sampled genome")
  })

  # ---- Narrative feed ----
  output[[ns("gt_narrative")]] <- renderUI({
    ev <- events_r(); ct <- current_time()
    if (!length(ev)) return(div(class = "text-muted", "No events to narrate."))
    shown <- Filter(function(e) !is.na(e$date) && e$date <= ct, ev)
    if (!length(shown))
      return(div(class = "text-muted", style = "padding:10px;", "Press play or advance the timeline to begin the reconstruction."))
    badge <- function(tier) switch(tier,
      confirmed  = tags$span(class = "badge", style = "background:#198754;color:#fff;", "CONFIRMED"),
      inferred   = tags$span(class = "badge", style = "background:#fd7e14;color:#fff;", "INFERRED"),
      probable   = tags$span(class = "badge", style = "background:#6c757d;color:#fff;", "PROBABLE"),
      contextual = tags$span(class = "badge", style = "background:#4A6C8C;color:#fff;", "CONTEXTUAL"),
      observed   = tags$span(class = "badge", style = "background:#2c3e50;color:#fff;", "OBSERVED"),
      tags$span(class = "badge badge-secondary", toupper(tier)))
    ic <- function(type) switch(type,
      origin = "circle-dot", spread = "arrow-right", outbreak = "triangle-exclamation",
      query = "star", link = "link", "circle")
    div(style = "max-height:320px;overflow-y:auto;border:1px solid #e9ecef;border-radius:6px;padding:8px;",
        lapply(shown, function(e)
          div(style = "padding:6px 8px;border-left:3px solid #dee2e6;margin-bottom:5px;background:#f8f9fa;border-radius:3px;font-size:0.88rem;",
              tags$small(class = "text-muted", style = "font-weight:600;", format(e$date, "%Y-%m")), "  ",
              icon(ic(e$type)), " ", badge(e$tier), " ", HTML(e$text))))
  })
}
