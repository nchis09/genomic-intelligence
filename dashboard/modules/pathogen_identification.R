# Biological Threat module (formerly "Pathogen Identification")
#
# Interactive comparison workspace: a multi-select sample picker drives a
# summary strip + plotly comparison charts above the existing DT results table.
# Produces one tab per species discovered under --outdir; output IDs are
# namespaced manually with the species name (no NS()).

# ---------------------------------------------------------------------------
# Table configs (unchanged data contract with bin/pathogen_identification.R)
# ---------------------------------------------------------------------------
PI_SUMMARY_TABLE <- list(
  title    = "Species Identification Summary",
  caption  = "Location in 'Closest reference' may read 'Not Provided' when the source record didn't include it.",
  filename = "identification_summary.tsv",
  colnames = c(
    "", "Sample", "Species identified", "Coverage", "QC overall score",
    "Closest reference (accession|location|date)",
    "% Identity to closest reference", "Genetic distance to closest reference",
    "Distance to assigned clade", "MSA identity",
    "Closest reference outbreak"
  ),
  pct_cols   = c("coverage", "pct_identity_closest_ref", "msa_identity"),
  round_cols = c("qc_overall_score", "genetic_distance_closest_ref", "distance_to_assigned_clade")
)

PI_UNRESOLVED_TABLE <- list(
  title    = "Unresolved Samples (Ambiguous / Rejected)",
  filename = "unresolved_samples.tsv",
  colnames = c("Sample", "Assignment", "Coverage", "QC Status"),
  pct_cols   = c("coverage"),
  round_cols = c()
)

# Consistent sample colors (max 5)
PI_COLORS <- c("#2ECC71", "#3498DB", "#F39C12", "#E74C3C", "#9B59B6")

# Metrics to chart: list(col, label, subtitle, higher_is_better, is_pct)
PI_METRICS <- list(
  list(col = "coverage",                    label = "Coverage (%)",                      subtitle = "Higher is better",             higher = TRUE,  pct = TRUE),
  list(col = "pct_identity_closest_ref",    label = "% Identity to\nclosest reference",  subtitle = "Higher is better",             higher = TRUE,  pct = TRUE),
  list(col = "genetic_distance_closest_ref",label = "Genetic distance to\nclosest reference", subtitle = "Lower is better",         higher = FALSE, pct = FALSE),
  list(col = "qc_overall_score",            label = "QC overall\nscore",                 subtitle = "Lower is better",              higher = FALSE, pct = FALSE),
  list(col = "distance_to_assigned_clade",  label = "Distance to\nassigned clade",       subtitle = "Lower is better",              higher = FALSE, pct = FALSE),
  list(col = "msa_identity",                label = "MSA Identity (%)",                  subtitle = "Higher is better\n(mean pairwise)", higher = TRUE,  pct = TRUE)
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
pi_table_path <- function(outdir, species, filename) {
  file.path(outdir, "pathogen_identification", species, "species_identification", filename)
}

pi_read_table_safe <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(readr::read_tsv(path, show_col_types = FALSE), error = function(e) NULL)
}

pi_id <- function(species, suffix) {
  paste0("pi_", species, "_", suffix)
}

# Build one plotly dot chart for a single metric
pi_metric_chart <- function(sel_df, metric, color_map) {
  col  <- metric$col
  vals <- sel_df[[col]]
  if (all(is.na(vals))) return(NULL)

  display_vals <- if (metric$pct) vals * 100 else vals
  labels <- sprintf(if (metric$pct) "%.1f%%" else "%.4f", display_vals)

  colors <- color_map[sel_df$sample]

  p <- plot_ly(
    x = sel_df$sample, y = display_vals,
    type = "scatter", mode = "markers+lines+text",
    marker = list(size = 12, color = colors),
    line = list(color = "#ccc", width = 1.5),
    text = labels, textposition = "top center",
    textfont = list(size = 11, color = colors),
    hoverinfo = "x+y",
    showlegend = FALSE
  ) %>%
    layout(
      title = list(text = paste0("<b>", gsub("\n", " ", metric$label), "</b>",
                                 "<br><sup>", metric$subtitle, "</sup>"),
                   font = list(size = 13), x = 0.5, y = 0.95),
      xaxis = list(title = "", tickangle = -30, tickfont = list(size = 10)),
      yaxis = list(
        title = "",
        tickformat = if (metric$pct) ".0f" else ".4f",
        range = list(0, max(display_vals, na.rm = TRUE) * 1.25),
        automargin = TRUE
      ),
      margin = list(t = 85, b = 50, l = 55, r = 20, pad = 8),
      paper_bgcolor = "transparent", plot_bgcolor = "transparent",
      height = 310
    ) %>%
    config(displayModeBar = FALSE)
  p
}

# ---------------------------------------------------------------------------
# UI for one species
# ---------------------------------------------------------------------------
pathogen_identification_ui <- function(species) {
  tagList(
    h3("Pathogen Identification"),
    p(style = "color: #6c757d; margin-bottom: 16px;",
      "Compare genomic evidence across selected samples"),

    # Sample selector bar
    div(class = "pi-sample-select-bar",
      selectizeInput(
        pi_id(species, "sample_select"),
        label = NULL,
        choices = NULL,
        multiple = TRUE,
        options = list(
          placeholder = "Select samples to compare",
          maxItems = 5,
          plugins = list("remove_button")
        )
      ),
      span(style = "white-space: nowrap;",
        uiOutput(pi_id(species, "select_count"), inline = TRUE)
      ),
      actionButton(pi_id(species, "reset_select"), "Reset selection",
                   icon = icon("rotate-left"),
                   class = "btn-outline-secondary btn-sm",
                   style = "white-space: nowrap;")
    ),

    # Summary strip (dynamic)
    uiOutput(pi_id(species, "summary_strip")),

    # Comparison charts (dynamic)
    uiOutput(pi_id(species, "chart_section")),

    # Results table
    bs4Dash::bs4Card(
      title = "Identification Results (All samples)",
      width = 12,
      collapsible = TRUE,
      DT::DTOutput(pi_id(species, "results_table"))
    ),

    # Unresolved samples
    bs4Dash::bs4Card(
      title = PI_UNRESOLVED_TABLE$title,
      width = 12,
      collapsible = TRUE,
      collapsed = TRUE,
      DT::DTOutput(pi_id(species, "unresolved_table"))
    )
  )
}

# ---------------------------------------------------------------------------
# Server registration for one species
# ---------------------------------------------------------------------------
pathogen_identification_register <- function(input, output, session, species, outdir) {
  sp <- species

  # -- Data reactive
  summary_data <- reactive({
    path <- pi_table_path(outdir(), sp, PI_SUMMARY_TABLE$filename)
    df <- pi_read_table_safe(path)
    validate(need(!is.null(df), paste0("No data found at: ", path)))
    df
  })

  # -- Populate selectize choices from data
  observe({
    df <- summary_data()
    updateSelectizeInput(session, pi_id(sp, "sample_select"),
                         choices = df$sample, selected = character(0),
                         server = FALSE)
  })

  # -- Reset button
  observeEvent(input[[pi_id(sp, "reset_select")]], {
    updateSelectizeInput(session, pi_id(sp, "sample_select"),
                         selected = character(0))
  })

  # -- Selected samples reactive
  selected_samples <- reactive({
    sel <- input[[pi_id(sp, "sample_select")]]
    if (is.null(sel)) character(0) else sel
  })

  # -- Color map reactive
  color_map <- reactive({
    sel <- selected_samples()
    if (length(sel) == 0) return(setNames(character(0), character(0)))
    setNames(PI_COLORS[seq_along(sel)], sel)
  })

  # -- Select count badge
  output[[pi_id(sp, "select_count")]] <- renderUI({
    n <- length(selected_samples())
    if (n == 0) return(NULL)
    span(class = "badge badge-primary", style = "font-size: 0.85rem; padding: 5px 10px;",
         paste0(n, " sample", ifelse(n != 1, "s", ""), " selected"))
  })

  # -- Summary strip
  output[[pi_id(sp, "summary_strip")]] <- renderUI({
    sel <- selected_samples()
    if (length(sel) == 0) return(NULL)
    df <- summary_data()
    sel_df <- df[df$sample %in% sel, , drop = FALSE]
    n_total <- nrow(df)

    species_name <- if (nrow(sel_df) > 0 && "species_identified" %in% names(sel_df)) {
      unique(sel_df$species_identified)[1]
    } else toupper(sp)

    outbreaks <- if ("closest_outbreak" %in% names(sel_df)) {
      unique(na.omit(sel_df$closest_outbreak))
    } else character(0)
    outbreak_text <- if (length(outbreaks) == 0) "N/A" else paste(outbreaks, collapse = ", ")

    div(class = "pi-summary-strip",
      div(class = "pi-summary-card",
        div(class = "pi-sc-label", "Pathogen (all selected samples)"),
        div(class = "pi-sc-value", toupper(species_name))
      ),
      div(class = "pi-summary-card",
        div(class = "pi-sc-label", "Samples compared"),
        div(class = "pi-sc-value", length(sel)),
        div(class = "pi-sc-detail", paste0("of ", n_total, " total"))
      ),
      div(class = "pi-summary-card",
        div(class = "pi-sc-label", "Closest outbreak(s)"),
        div(class = "pi-sc-value", style = "font-size: 1.1rem;", outbreak_text)
      ),
      div(class = "pi-info-card",
        icon("circle-info"),
        span("The charts below compare key genomic metrics across the selected samples.",
             "Select up to 5 samples to compare.")
      )
    )
  })

  # -- Comparison charts section
  output[[pi_id(sp, "chart_section")]] <- renderUI({
    sel <- selected_samples()
    if (length(sel) < 2) return(NULL)

    chart_outputs <- lapply(seq_along(PI_METRICS), function(i) {
      column(width = 4, style = "margin-bottom: 12px;",
        plotly::plotlyOutput(pi_id(sp, paste0("chart_", i)), height = "310px")
      )
    })

    cmap <- color_map()
    legend_items <- lapply(names(cmap), function(s) {
      div(class = "pi-legend-item",
        span(class = "pi-legend-swatch", style = paste0("background:", cmap[[s]], ";")),
        s
      )
    })

    div(class = "pi-chart-section",
      div(style = "display: flex; justify-content: space-between; align-items: center; margin-bottom: 8px;",
        h5("Comparison of Genomic Metrics", style = "margin: 0;"),
        actionButton(pi_id(sp, "reset_select2"), "Reset selection",
                     icon = icon("rotate-left"),
                     class = "btn-outline-secondary btn-sm")
      ),
      fluidRow(chart_outputs),
      div(class = "pi-legend", legend_items)
    )
  })

  # -- Render each metric chart
  lapply(seq_along(PI_METRICS), function(i) {
    output[[pi_id(sp, paste0("chart_", i))]] <- plotly::renderPlotly({
      sel <- selected_samples()
      req(length(sel) >= 2)
      df <- summary_data()
      sel_df <- df[df$sample %in% sel, , drop = FALSE]
      sel_df <- sel_df[match(sel, sel_df$sample), , drop = FALSE]
      pi_metric_chart(sel_df, PI_METRICS[[i]], color_map())
    })
  })

  # -- Second reset button (inside chart section)
  observeEvent(input[[pi_id(sp, "reset_select2")]], {
    updateSelectizeInput(session, pi_id(sp, "sample_select"),
                         selected = character(0))
  })

  # -- Results table with checkbox column (renders only when data changes,
  # not when selection changes -- checkbox state is synced via JS observer)
  output[[pi_id(sp, "results_table")]] <- DT::renderDT({
    df <- summary_data()

    # Add checkbox column (all unchecked initially; JS observer syncs state)
    cb <- sapply(seq_len(nrow(df)), function(i) {
      s <- df$sample[i]
      sprintf(
        '<input type="checkbox" class="pi-row-cb" data-sample="%s" data-species="%s"/>',
        htmltools::htmlEscape(s), htmltools::htmlEscape(sp)
      )
    })
    display_df <- cbind(data.frame(` ` = cb, check.names = FALSE), df)

    dt <- DT::datatable(
      display_df,
      colnames = PI_SUMMARY_TABLE$colnames,
      caption = PI_SUMMARY_TABLE$caption,
      escape = FALSE,
      selection = "none",
      options = list(
        pageLength = 15, scrollX = TRUE,
        columnDefs = list(
          list(orderable = FALSE, className = "dt-center", targets = 0)
        ),
        drawCallback = DT::JS("
          function(settings) {
            var table = this;
            table.find('.pi-row-cb').off('change').on('change', function() {
              var sample = $(this).data('sample');
              var species = $(this).data('species');
              var checked = this.checked;
              var inputId = 'pi_' + species + '_table_toggle';
              Shiny.setInputValue(inputId,
                {sample: sample, checked: checked, ts: Date.now()},
                {priority: 'event'});
            });
          }
        ")
      )
    )

    pct_cols <- intersect(PI_SUMMARY_TABLE$pct_cols, names(df))
    if (length(pct_cols) > 0) {
      pct_idx <- which(names(display_df) %in% pct_cols)
      dt <- DT::formatPercentage(dt, pct_idx, digits = 1)
    }
    round_cols <- intersect(PI_SUMMARY_TABLE$round_cols, names(df))
    if (length(round_cols) > 0) {
      round_idx <- which(names(display_df) %in% round_cols)
      dt <- DT::formatRound(dt, round_idx, digits = 4)
    }
    dt
  })

  # -- Table checkbox -> selectize sync
  observeEvent(input[[pi_id(sp, "table_toggle")]], {
    evt <- input[[pi_id(sp, "table_toggle")]]
    current <- selected_samples()
    s <- evt$sample
    if (isTRUE(evt$checked)) {
      if (!(s %in% current) && length(current) < 5) {
        updateSelectizeInput(session, pi_id(sp, "sample_select"),
                             selected = c(current, s))
      }
    } else {
      updateSelectizeInput(session, pi_id(sp, "sample_select"),
                           selected = setdiff(current, s))
    }
  })

  # -- Unresolved table
  unresolved_data <- reactive({
    path <- pi_table_path(outdir(), sp, PI_UNRESOLVED_TABLE$filename)
    df <- pi_read_table_safe(path)
    validate(need(!is.null(df), paste0("No data found at: ", path)))
    df
  })

  output[[pi_id(sp, "unresolved_table")]] <- DT::renderDT({
    df <- unresolved_data()
    dt <- DT::datatable(df, colnames = PI_UNRESOLVED_TABLE$colnames,
                        options = list(pageLength = 15, scrollX = TRUE))
    pct_cols <- intersect(PI_UNRESOLVED_TABLE$pct_cols, names(df))
    if (length(pct_cols) > 0) dt <- DT::formatPercentage(dt, pct_cols, digits = 1)
    dt
  })
}
