# Geographic & Temporal Context module
#
# Where and when the new genome's closest historical relatives were detected.
# Reads transmission_potential.tsv from results/transmission_context/<species>/.
# Reuses the .ts_* helpers defined in pathogen_transmission.R (sourced first).

geographic_temporal_ui <- function() {
  tagList(
    br(),
    h4(icon("earth-africa"), " Geographic & Temporal Context ", tags$small(class = "text-muted", "observed history")),
    p(class = "text-muted", style = "font-size:0.85rem;",
      "Where and when the new genome's closest historical relatives were detected. Drag or animate the year slider to watch the spread unfold."),
    fluidRow(
      column(3, uiOutput("gt_query_select_ui")),
      column(9, uiOutput("gt_expansion_slider_ui"))
    ),
    leaflet::leafletOutput("gt_expansion_map", height = "540px")
  )
}

geographic_temporal_register <- function(input, output, session, outdir_r, current_species) {
  potential_r <- reactive({ sp <- current_species(); if (is.null(sp)) NULL else .ts_read(outdir_r(), sp, "transmission_potential.tsv") })

  query_samples_r <- reactive({ p <- potential_r(); if (is.null(p)) NULL else sort(unique(stats::na.omit(p$query_sample))) })
  output$gt_query_select_ui <- renderUI({
    qs <- query_samples_r(); if (is.null(qs)) return(div(style = "color:#888;", "No query samples."))
    selectInput("gt_query", "Query genome:", choices = qs, selected = qs[1], selectize = FALSE)
  })
  selected_query <- reactive({ qs <- query_samples_r(); if (is.null(qs)) NULL else if (is.null(input$gt_query) || !(input$gt_query %in% qs)) qs[1] else input$gt_query })
  links_r <- reactive({ p <- potential_r(); sq <- selected_query(); if (is.null(p) || is.null(sq)) NULL else dplyr::filter(p, query_sample == sq) })

  output$gt_expansion_slider_ui <- renderUI({
    df <- links_r(); if (is.null(df) || !nrow(df)) return(NULL)
    y <- suppressWarnings(as.numeric(format(as.Date(df$background_collection_date), "%Y"))); y <- y[!is.na(y)]
    if (!length(y)) return(NULL)
    sliderInput("gt_expansion_year", "Detections up to year:", min = min(y), max = max(y), value = max(y), step = 1, sep = "", animate = TRUE)
  })

  output$gt_expansion_map <- leaflet::renderLeaflet({
    df <- links_r(); m <- leaflet::leaflet() %>% leaflet::addTiles(); if (is.null(df) || !nrow(df)) return(m)
    df <- .ts_fill_coords(df, "background_latitude", "background_longitude", "background_country")
    d <- df %>% dplyr::filter(!is.na(background_latitude), !is.na(background_collection_date)) %>%
      dplyr::mutate(yr = as.numeric(format(as.Date(background_collection_date), "%Y")))
    if (!nrow(d)) return(m)
    if (!is.null(input$gt_expansion_year)) d <- dplyr::filter(d, yr <= input$gt_expansion_year)
    if (!nrow(d)) return(m)
    pal <- leaflet::colorFactor(.ts_strain_pal(d$background_strain), domain = d$background_strain)
    m %>% leaflet::addCircleMarkers(data = d, ~background_longitude, ~background_latitude, radius = 6,
        color = ~pal(background_strain), fillOpacity = 0.7,
        popup = ~paste0("<b>", background_sample, "</b><br>", background_country, "<br>", background_collection_date)) %>%
      leaflet::addLegend(position = "bottomright", pal = pal, values = d$background_strain, title = "Strain")
  })
}
