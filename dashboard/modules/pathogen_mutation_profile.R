# Pathogen Mutation Profile module
#
# Dashboard tab for analysis #1: overall mutation profile (query vs background).
# Reads published TSVs from outdir/pathogen_mutation_profile/<species>/
# and renders an interactive boxplot + summary tables.

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
mp_id <- function(species, suffix) {
  paste0("mp_", species, "_", suffix)
}

mp_table_path <- function(outdir, species, filename) {
  file.path(outdir, "pathogen_mutation_profile", species, "mutation_profile", filename)
}

mp_read_table <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(readr::read_tsv(path, show_col_types = FALSE), error = function(e) NULL)
}

# ---------------------------------------------------------------------------
# UI for one species
# ---------------------------------------------------------------------------
pathogen_mutation_profile_ui <- function(species) {
  tagList(
    h3("Pathogen Mutation Profile"),
    p(style = "color: #6c757d; margin-bottom: 16px;",
      "Query-vs-background mutation burden and protein-level MANOVA."),

    fluidRow(
      column(9,
        bs4Dash::bs4Card(
          title = "PLS-DA Score Plot (PC1 vs PC2)",
          width = 12,
          status = "primary",
          solidHeader = TRUE,
          radioButtons(
            inputId = mp_id(species, "plsda_colour"),
            label = "Colour samples by:",
            choices = c("country", "outbreak"),
            selected = "country",
            inline = TRUE
          ),
          div(style = "position: relative;",
              plotOutput(mp_id(species, "plsda_score_plot"), click = mp_id(species, "plsda_click"), height = "400px"),
              uiOutput(mp_id(species, "plsda_click_info"))
          )
        )
      ),
      column(3,
        bs4Dash::bs4Card(
          title = "Find a sample",
          width = 12,
          status = "info",
          solidHeader = TRUE,
          selectizeInput(
            inputId = mp_id(species, "plsda_sample_search"),
            label = NULL,
            choices = NULL,
            selected = NULL,
            options = list(placeholder = "Search sample..."),
            width = "100%"
          ),
          uiOutput(mp_id(species, "plsda_search_info"))
        )
      )
    ),

    fluidRow(
      column(12,
        bs4Dash::bs4Card(
          title = "MANOVA",
          width = 12,
          status = "info",
          solidHeader = TRUE,
          DT::DTOutput(mp_id(species, "manova_table"))
        )
      )
    ),

    fluidRow(
      column(9,
        bs4Dash::bs4Card(
          title = "Per-Protein Mutation Positions (Background vs Query)",
          width = 12,
          status = "warning",
          solidHeader = TRUE,
          plotOutput(mp_id(species, "protein_burden_summary_plot"),
                    click = mp_id(species, "protein_burden_click"),
                    height = "400px")
        )
      ),
      column(3,
        bs4Dash::bs4Card(
          title = "Select query sample",
          width = 12,
          status = "info",
          solidHeader = TRUE,
          selectizeInput(
            inputId = mp_id(species, "protein_sample_search"),
            label = NULL,
            choices = NULL,
            selected = NULL,
            options = list(placeholder = "Search query sample..."),
            width = "100%"
          ),
          uiOutput(mp_id(species, "protein_sample_info"))
        )
      )
    ),

    fluidRow(
      column(12,
        bs4Dash::bs4Card(
          title = "Mutation Landscape",
          width = 12,
          status = "warning",
          solidHeader = TRUE,
          div(
            style = "display: flex; gap: 16px; align-items: center; margin-bottom: 12px; flex-wrap: wrap;",
            div(style = "width: 220px;",
              selectizeInput(
                inputId = mp_id(species, "landscape_protein"),
                label = NULL,
                choices = NULL,
                selected = NULL,
                multiple = FALSE,
                options = list(placeholder = "Select protein..."),
                width = "100%"
              )
            ),
            checkboxInput(
              mp_id(species, "landscape_highlight_query"),
              label = "Highlight query mutations",
              value = TRUE
            )
          ),
          plotOutput(mp_id(species, "landscape_plot"),
                     click = mp_id(species, "landscape_click"),
                     brush = brushOpts(id = mp_id(species, "landscape_brush"), resetOnNew = TRUE),
                     height = "520px"),
          actionButton(mp_id(species, "landscape_reset_zoom"), "Reset zoom", class = "btn-sm")
        )
      )
    ),

    fluidRow(
      column(12,
        uiOutput(mp_id(species, "position_detail_panel"))
      )
    ),

    fluidRow(
      column(12,
        uiOutput(mp_id(species, "mutation_intelligence_card"))
      )
    ),

    fluidRow(
      column(12,
        bs4Dash::bs4Card(
          title = "Data Table / Export",
          width = 12,
          status = "secondary",
          solidHeader = TRUE,
          collapsed = TRUE,
          collapsible = TRUE,
          DT::DTOutput(mp_id(species, "mutation_data_table"))
        )
      )
    )

  )
}

# ---------------------------------------------------------------------------
# Server registration for one species
# ---------------------------------------------------------------------------
pathogen_mutation_profile_register <- function(input, output, session, species, outdir) {
  sp <- species

  manova_data <- reactive({
    mp_read_table(mp_table_path(outdir(), sp, "01_mutation_profile_manova.tsv"))
  })

  protein_samples_data <- reactive({
    mp_read_table(mp_table_path(outdir(), sp, "01_protein_burden_samples.tsv"))
  })

  protein_summary_data <- reactive({
    mp_read_table(mp_table_path(outdir(), sp, "01_protein_burden_summary.tsv"))
  })

  plsda_scores_data <- reactive({
    mp_read_table(mp_table_path(outdir(), sp, "01_plsda_scores.tsv"))
  })



  mutation_detail_data <- reactive({
    mp_read_table(mp_table_path(outdir(), sp, "01_mutation_detail.tsv"))
  })

  mutation_phenotypes_data <- reactive({
    mp_read_table(mp_table_path(outdir(), sp, "01_mutation_phenotypes.tsv"))
  })

  aa_frequencies_data <- reactive({
    mp_read_table(mp_table_path(outdir(), sp, "01_position_aa_frequencies.tsv"))
  })

  click_info <- reactiveVal(NULL)
  protein_selected_sample <- reactiveVal(NULL)
  selected_mutation_id <- reactiveVal(NULL)

  observe({
    df <- plsda_scores_data()
    if (is.null(df) || nrow(df) == 0) return()
    ids <- as.character(df$sample_id)
    labs <- paste0(df$sample_name, " (", df$sample_id, ")")
    if (length(ids) != length(labs)) return()
    choices <- setNames(ids, labs)
    updateSelectizeInput(
      session,
      inputId = mp_id(sp, "plsda_sample_search"),
      choices = choices,
      selected = NULL,
      server = FALSE
    )
  })

  observe({
    df <- protein_samples_data()
    if (is.null(df) || nrow(df) == 0) return()
    query_df <- df |> filter(is_query)
    if (nrow(query_df) == 0) return()
    choices <- unique(query_df$sample_name)
    updateSelectizeInput(
      session,
      inputId = mp_id(sp, "protein_sample_search"),
      choices = choices,
      selected = character(0),
      server = FALSE
    )
  })

  observeEvent(input[[mp_id(sp, "protein_sample_search")]], {
    val <- input[[mp_id(sp, "protein_sample_search")]]
    if (!is.null(val) && val != "") {
      protein_selected_sample(val)
    } else {
      protein_selected_sample(NULL)
    }
  }, ignoreNULL = FALSE)

  observeEvent(input[[mp_id(sp, "protein_burden_click")]], {
    click <- input[[mp_id(sp, "protein_burden_click")]]
    if (is.null(click)) return()
    df <- protein_samples_data()
    if (is.null(df) || nrow(df) == 0) return()
    if (!"mutation_positions" %in% names(df)) {
      df <- df |> rename(mutation_positions = mutation_count)
    }
    query_df <- df |> filter(is_query)
    if (nrow(query_df) == 0) return()
    selected <- nearPoints(query_df, click,
                           xvar = "protein_name", yvar = "mutation_positions",
                           threshold = 10, maxpoints = 1)
    if (nrow(selected) > 0) {
      sname <- selected$sample_name[1]
      protein_selected_sample(sname)
      updateSelectizeInput(session, mp_id(sp, "protein_sample_search"),
                           selected = sname)
    }
  })

  output[[mp_id(sp, "manova_table")]] <- DT::renderDT({
    req(manova_data())
    DT::datatable(
      manova_data(),
      options = list(pageLength = 10, scrollX = TRUE),
      rownames = FALSE
    ) |>
      DT::formatRound(columns = c("Wilks", "statistic", "num_Df", "den_Df"), digits = 3) |>
      DT::formatSignif(columns = c("p_value"), digits = 3)
  })

  output[[mp_id(sp, "protein_burden_summary_plot")]] <- renderPlot({
    df <- protein_samples_data()
    if (is.null(df) || nrow(df) == 0) {
      plot.new()
      text(0.5, 0.5, "No protein position data available.", cex = 1.1, col = "#6c757d")
      return()
    }

    if (!"mutation_positions" %in% names(df)) {
      df <- df |> rename(mutation_positions = mutation_count)
    }

    pval <- protein_summary_data()
    if (!is.null(pval) && nrow(pval) > 0) {
      pval <- pval |>
        filter(is_query) |>
        distinct(protein_name, .keep_all = TRUE) |>
        dplyr::select(protein_name, p_value)
      max_y <- df |>
        group_by(protein_name) |>
        summarise(y = max(mutation_positions, na.rm = TRUE), .groups = "drop")
      pval <- pval |>
        left_join(max_y, by = "protein_name") |>
        mutate(
          label = ifelse(p_value < 0.001, "p < 0.001", paste0("p = ", signif(p_value, 2)))
        )
    } else {
      pval <- NULL
    }

    p <- ggplot(df, aes(x = .data[["protein_name"]], y = .data[["mutation_positions"]])) +
      geom_boxplot(
        data = df |> filter(!is_query),
        aes(x = .data[["protein_name"]], y = .data[["mutation_positions"]]),
        fill = "#BDC3C7",
        alpha = 0.7,
        outlier.shape = NA
      ) +
      geom_jitter(
        data = df |> filter(is_query),
        aes(x = .data[["protein_name"]], y = .data[["mutation_positions"]]),
        colour = "#E74C3C",
        size = 3,
        width = 0.2,
        height = 0
      ) +
      labs(x = "Protein", y = "Amino acid substitution positions per sample") +
      theme_bw(base_size = 12) +
      theme(legend.position = "none")

    sel <- protein_selected_sample()
    if (!is.null(sel) && sel != "") {
      sel_df <- df |> filter(is_query, .data[["sample_name"]] == sel)
      if (nrow(sel_df) > 0) {
        p <- p + geom_point(
          data = sel_df,
          aes(x = .data[["protein_name"]], y = .data[["mutation_positions"]]),
          colour = "#8E44AD",
          size = 4
        )
      }
    }

    if (!is.null(pval) && nrow(pval) > 0) {
      p <- p + geom_text(
        data = pval,
        aes(x = .data[["protein_name"]], y = .data[["y"]], label = .data[["label"]]),
        vjust = -0.5,
        size = 3,
        colour = "#2C3E50"
      )
    }

    p
  })

  output[[mp_id(sp, "protein_sample_info")]] <- renderUI({
    sel <- protein_selected_sample()
    if (is.null(sel) || sel == "") {
      return(div(style = "margin-top: 10px; color: #6c757d; font-size: 0.85rem;",
                 "Click a dot or search to select a sample."))
    }
    div(
      style = "margin-top: 10px; padding: 10px; background-color: #f8f9fa; border-radius: 4px;",
      h5(sel, style = "margin: 0; font-weight: 600; color: #8E44AD;")
    )
  })

  output[[mp_id(sp, "protein_burden_samples_table")]] <- DT::renderDT({
    req(protein_samples_data())
    DT::datatable(
      protein_samples_data(),
      options = list(pageLength = 15, scrollX = TRUE),
      rownames = FALSE
    ) |>
      DT::formatRound(columns = c("mutation_positions"), digits = 0)
  })

  # =========================================================================
  # Mutation Landscape + Position Detail + Intelligence Card
  # =========================================================================

  clicked_position <- reactiveVal(NULL)
  landscape_zoom <- reactiveVal(NULL)

  # Populate protein selector for the landscape
  observe({
    df <- mutation_detail_data()
    if (is.null(df) || nrow(df) == 0) return()
    proteins <- sort(unique(df$protein_name))
    updateSelectizeInput(session, mp_id(sp, "landscape_protein"),
                         choices = proteins, selected = proteins[1])
  })

  # Aggregated mutation data: one row per unique mutation, all proteins
  .catalogue_data <- reactive({
    df <- mutation_detail_data()
    if (is.null(df) || nrow(df) == 0) return(NULL)

    n_total <- length(unique(df$sample_id))
    n_query <- length(unique(df$sample_id[df$is_query]))
    n_bg <- n_total - n_query

    agg <- df |>
      group_by(mutation_id, mutation_label, protein_name, position, ref_aa, alt_aa, mutation_type) |>
      summarise(
        n_samples = dplyr::n_distinct(sample_id),
        n_query_mut = dplyr::n_distinct(sample_id[is_query]),
        n_bg_mut = dplyr::n_distinct(sample_id[!is_query]),
        query_samples = paste(unique(sample_name[is_query]), collapse = ", "),
        .groups = "drop"
      ) |>
      mutate(
        n_query_total = n_query,
        n_bg_total = n_bg,
        query_prev_pct = if (n_query > 0) round(100 * n_query_mut / n_query, 1) else 0,
        bg_prev_pct = if (n_bg > 0) round(100 * n_bg_mut / n_bg, 1) else 0,
        enrichment = case_when(
          n_bg_mut == 0 & n_query_mut > 0 ~ Inf,
          n_query_mut == 0 ~ 0,
          TRUE ~ round((n_query_mut / n_query) / (n_bg_mut / n_bg), 2)
        ),
        in_query = n_query_mut > 0
      )

    agg
  })

  # --- Mutation landscape lollipop (ProteinPaint style) ---
  output[[mp_id(sp, "landscape_plot")]] <- renderPlot({
    prot <- input[[mp_id(sp, "landscape_protein")]]
    show_query <- isTRUE(input[[mp_id(sp, "landscape_highlight_query")]])
    agg <- .catalogue_data()

    if (is.null(agg) || is.null(prot) || prot == "") {
      plot.new()
      text(0.5, 0.5, "Select a protein to view its mutation landscape.",
           cex = 1.1, col = "#6c757d")
      return()
    }

    prot_data <- agg |> filter(protein_name == prot)
    if (nrow(prot_data) == 0) {
      plot.new()
      text(0.5, 0.5, paste("No mutations found in", prot),
           cex = 1.1, col = "#6c757d")
      return()
    }

    pheno <- mutation_phenotypes_data()
    pheno_ids <- if (!is.null(pheno) && nrow(pheno) > 0) unique(pheno$mutation_id) else character(0)

    plot_df <- prot_data |>
      mutate(
        n_total = n_query_mut + n_bg_mut,
        mutation_class = coalesce(mutation_type, "Unknown")
      ) |>
      filter(n_total > 0) |>
      group_by(position) |>
      arrange(desc(n_total), .by_group = TRUE) |>
      mutate(y = -row_number()) |>
      ungroup() |>
      mutate(
        stroke_col = if_else(
          mutation_id %in% pheno_ids,
          "#27AE60",
          if_else(show_query & in_query, "#E74C3C", "transparent")
        ),
        count_label = if_else(n_total >= 1, as.character(n_total), "")
      )

    if (nrow(plot_df) == 0) {
      plot.new()
      text(0.5, 0.5, "No mutations to display.",
           cex = 1.1, col = "#6c757d")
      return()
    }

    z <- landscape_zoom()
    if (!is.null(z)) {
      zmin <- min(z, na.rm = TRUE)
      zmax <- max(z, na.rm = TRUE)
      plot_df <- plot_df |>
        filter(position >= zmin, position <= zmax)
    }

    if (nrow(plot_df) == 0) {
      plot.new()
      text(0.5, 0.5, "No mutations in selected zoom region.",
           cex = 1.1, col = "#6c757d")
      return()
    }

    x_limits <- range(plot_df$position) + c(-1, 1)
    y_min <- min(plot_df$y, -1, na.rm = TRUE) - 1.5

    p <- ggplot(plot_df, aes(x = position, y = y)) +
      geom_segment(aes(xend = position, y = 0, yend = y),
                   colour = "#BDC3C7", linewidth = 0.5) +
      geom_point(aes(size = n_total, fill = mutation_class, colour = stroke_col),
                 shape = 21, stroke = 1.2) +
      geom_text(aes(label = count_label),
                colour = "white", size = 2.2, fontface = "bold",
                vjust = 0.5, hjust = 0.5) +
      ggrepel::geom_text_repel(
        aes(label = mutation_label),
        colour = "#2C3E50", size = 2.4,
        nudge_x = 3, nudge_y = 0,
        hjust = 0, vjust = 0.5,
        force = 2, force_pull = 0.5,
        box.padding = 0.3, point.padding = 0.3,
        max.overlaps = Inf,
        segment.size = 0.25, min.segment.length = 0
      ) +
      scale_size_continuous(range = c(4, 16), guide = "none") +
      scale_fill_brewer(palette = "Set1", na.value = "#95A5A6",
                        name = "Mutation class") +
      scale_colour_identity(guide = "none") +
      scale_x_continuous(limits = x_limits) +
      scale_y_continuous(expand = expansion(add = c(0.5, 0.5))) +
      coord_cartesian(ylim = c(y_min, 0.5)) +
      labs(
        x = "Amino acid position",
        y = NULL,
        title = paste(prot, "\u2014 Mutation Landscape"),
        caption = if (show_query) "Red outline = also in query sample(s)" else NULL
      ) +
      theme_bw(base_size = 11) +
      theme(
        plot.title = element_text(size = 13, face = "bold"),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = "right",
        legend.key.size = unit(0.4, "cm"),
        plot.margin = margin(t = 5, r = 5, b = 0, l = 5)
      )

    p
  })

  # --- Landscape zoom controls ---
  observeEvent(input[[mp_id(sp, "landscape_brush")]], {
    br <- input[[mp_id(sp, "landscape_brush")]]
    if (is.null(br) || is.null(br$xmin) || is.null(br$xmax)) return()
    landscape_zoom(c(br$xmin, br$xmax))
  })

  observeEvent(input[[mp_id(sp, "landscape_reset_zoom")]], {
    landscape_zoom(NULL)
  })

  # --- Landscape click handler ---
  observeEvent(input[[mp_id(sp, "landscape_click")]], {
    click <- input[[mp_id(sp, "landscape_click")]]
    if (is.null(click)) return()

    prot <- input[[mp_id(sp, "landscape_protein")]]
    agg <- .catalogue_data()
    if (is.null(agg) || is.null(prot) || prot == "") return()

    prot_data <- agg |> filter(protein_name == prot)
    if (nrow(prot_data) == 0) return()

    plot_df <- prot_data |>
      mutate(n_total = n_query_mut + n_bg_mut) |>
      filter(n_total > 0) |>
      group_by(position) |>
      arrange(desc(n_total), .by_group = TRUE) |>
      mutate(y = -row_number()) |>
      ungroup()

    z <- landscape_zoom()
    if (!is.null(z)) {
      zmin <- min(z, na.rm = TRUE)
      zmax <- max(z, na.rm = TRUE)
      plot_df <- plot_df |>
        filter(position >= zmin, position <= zmax)
    }

    if (nrow(plot_df) == 0) {
      selected_mutation_id(NULL)
      clicked_position(NULL)
      return()
    }

    selected <- nearPoints(plot_df, click,
                           xvar = "position", yvar = "y",
                           threshold = 25, maxpoints = 1)
    if (nrow(selected) == 0) {
      selected_mutation_id(NULL)
      clicked_position(NULL)
      return()
    }

    selected_mutation_id(selected$mutation_id[1])
    clicked_position(selected$position[1])
  })

  # --- Position detail panel (shown when multiple substitutions at one position) ---
  output[[mp_id(sp, "position_detail_panel")]] <- renderUI({
    pos <- clicked_position()
    if (is.null(pos)) return(NULL)

    prot <- input[[mp_id(sp, "landscape_protein")]]
    agg <- .catalogue_data()
    if (is.null(agg) || is.null(prot)) return(NULL)

    subs <- agg |> filter(protein_name == prot, position == pos) |>
      arrange(desc(bg_prev_pct))
    if (nrow(subs) == 0) return(NULL)

    pheno <- mutation_phenotypes_data()
    pheno_ids <- if (!is.null(pheno) && nrow(pheno) > 0) unique(pheno$mutation_id) else character(0)

    bs4Dash::bs4Card(
      title = paste0(prot, " position ", pos, " \u2014 ", nrow(subs), " substitutions"),
      width = 12,
      status = "info",
      solidHeader = TRUE,
      collapsible = FALSE,
      div(
        style = "display: flex; flex-wrap: wrap; gap: 12px;",
        lapply(seq_len(nrow(subs)), function(i) {
          row <- subs[i, ]
          is_pheno <- row$mutation_id %in% pheno_ids
          border_col <- if (row$in_query) "#E74C3C"
                        else if (is_pheno) "#27AE60"
                        else "#dee2e6"
          bg_col <- if (row$in_query) "#FDEDEC" else "#f8f9fa"

          actionLink(
            inputId = mp_id(sp, paste0("pos_sub_", row$mutation_id)),
            label = div(
              style = paste0(
                "padding: 10px 14px; border: 2px solid ", border_col, ";",
                " border-radius: 6px; background: ", bg_col, ";",
                " cursor: pointer; min-width: 140px; text-align: center;"
              ),
              strong(row$mutation_label, style = "font-size: 1rem; display: block;"),
              if (!is.na(row$ref_aa) && !is.na(row$alt_aa)) {
                span(paste0(row$ref_aa, " → ", row$alt_aa),
                     style = "font-size: 0.85rem; color: #495057; display: block;")
              } else NULL,
              span(paste0("Query: ", row$query_prev_pct, "%"),
                   style = "font-size: 0.8rem; color: #E74C3C; display: block;"),
              span(paste0("Background: ", row$bg_prev_pct, "%"),
                   style = "font-size: 0.8rem; color: #6c757d; display: block;"),
              if (is_pheno) span("Has biological annotation",
                                 style = "font-size: 0.7rem; color: #27AE60; display: block; margin-top: 2px;")
            )
          )
        })
      )
    )
  })

  # Wire position detail clicks to intelligence card
  observe({
    prot <- input[[mp_id(sp, "landscape_protein")]]
    agg <- .catalogue_data()
    pos <- clicked_position()
    if (is.null(agg) || is.null(prot) || is.null(pos)) return()

    subs <- agg |> filter(protein_name == prot, position == pos)
    if (nrow(subs) == 0) return()

    lapply(seq_len(nrow(subs)), function(i) {
      mid <- subs$mutation_id[i]
      observeEvent(input[[mp_id(sp, paste0("pos_sub_", mid))]], {
        clicked_position(NULL)
        selected_mutation_id(mid)
      }, ignoreInit = TRUE, once = TRUE)
    })
  })

  # --- Data Table / Export (collapsed view) ---
  output[[mp_id(sp, "mutation_data_table")]] <- DT::renderDT({
    agg <- .catalogue_data()
    if (is.null(agg) || nrow(agg) == 0) {
      return(DT::datatable(data.frame(Message = "No mutations to display."),
                           rownames = FALSE, options = list(dom = "t")))
    }

    display_df <- agg |>
      mutate(
        enrichment_label = case_when(
          is.infinite(enrichment) ~ "Query-only",
          enrichment == 0 ~ "Absent in query",
          enrichment > 2 ~ paste0(enrichment, "x elevated"),
          enrichment < 0.5 ~ paste0(enrichment, "x lower"),
          TRUE ~ paste0(enrichment, "x (similar)")
        )
      ) |>
      dplyr::select(
        Protein = protein_name,
        Mutation = mutation_label,
        Position = position,
        Ref = ref_aa,
        Alt = alt_aa,
        `Query (%)` = query_prev_pct,
        `Background (%)` = bg_prev_pct,
        Enrichment = enrichment_label,
        `Query Samples` = query_samples
      ) |>
      arrange(Protein, Position)

    DT::datatable(
      display_df,
      selection = "single",
      extensions = "Buttons",
      options = list(
        pageLength = 20, scrollX = TRUE,
        dom = "Bfrtip",
        buttons = list("csv", "excel")
      ),
      rownames = FALSE
    )
  })

  # Data table row selection -> intelligence card
  observeEvent(input[[mp_id(sp, "mutation_data_table_rows_selected")]], {
    sel_row <- input[[mp_id(sp, "mutation_data_table_rows_selected")]]
    if (is.null(sel_row) || length(sel_row) == 0) {
      selected_mutation_id(NULL)
      return()
    }
    agg <- .catalogue_data()
    if (is.null(agg)) return()
    # Data table is sorted by Protein, Position — rebuild same order
    sorted_agg <- agg |> arrange(protein_name, position)
    if (sel_row > nrow(sorted_agg)) return()
    selected_mutation_id(sorted_agg$mutation_id[sel_row])
  })

  # --- Mutation Intelligence Card ---
  output[[mp_id(sp, "mutation_intelligence_card")]] <- renderUI({
    mid <- selected_mutation_id()
    df <- mutation_detail_data()
    if (is.null(mid) || is.null(df) || nrow(df) == 0) return(NULL)

    mut_rows <- df |> filter(mutation_id == mid)
    mut_label <- mut_rows$mutation_label[1]
    protein <- mut_rows$protein_name[1]
    pos <- mut_rows$position[1]
    ref_aa_val <- mut_rows$ref_aa[1]
    alt_aa_val <- mut_rows$alt_aa[1]

    n_query_total <- length(unique(df$sample_id[df$is_query]))
    n_bg_total <- length(unique(df$sample_id[!df$is_query]))
    n_query_mut <- length(unique(mut_rows$sample_id[mut_rows$is_query]))
    n_bg_mut <- length(unique(mut_rows$sample_id[!mut_rows$is_query]))

    query_prev <- if (n_query_total > 0) round(100 * n_query_mut / n_query_total, 1) else 0
    bg_prev <- if (n_bg_total > 0) round(100 * n_bg_mut / n_bg_total, 1) else 0
    query_samples_list <- unique(mut_rows$sample_name[mut_rows$is_query])

    status <- if (n_query_mut > 0 & n_bg_mut == 0) "Query-only (novel)"
              else if (query_prev > bg_prev + 20) "Elevated in query"
              else if (query_prev < bg_prev - 20) "Lower in query"
              else if (bg_prev >= 80) "Common (fixed/near-fixed)"
              else "Similar to background"

    status_color <- switch(status,
      "Query-only (novel)" = "#9B59B6",
      "Elevated in query" = "#E74C3C",
      "Lower in query" = "#3498DB",
      "Common (fixed/near-fixed)" = "#95A5A6",
      "#27AE60"
    )

    # Phenotype section
    pheno <- mutation_phenotypes_data()
    pheno_ui <- NULL
    has_pheno <- FALSE
    if (!is.null(pheno) && nrow(pheno) > 0) {
      matched <- pheno |> filter(mutation_id == mid)
      if (nrow(matched) > 0) {
        has_pheno <- TRUE
        pheno_ui <- div(
          style = "margin-top: 12px; padding: 12px; background: #f0fff0; border-left: 4px solid #27AE60; border-radius: 4px;",
          h6("Phenotype Characteristics", style = "margin: 0 0 8px 0; font-weight: 700; color: #27AE60;"),
          lapply(seq_len(nrow(matched)), function(i) {
            row <- matched[i, ]
            div(style = "margin-bottom: 6px;",
              p(style = "margin: 2px 0; font-size: 0.85rem;",
                strong("Phenotype: "), row$phenotype,
                if (!is.na(row$effect)) paste0(" | Effect: ", row$effect) else "",
                if (!is.na(row$source)) paste0(" | Source: ", row$source) else "")
            )
          })
        )
      }
    }
    if (!has_pheno) {
      pheno_ui <- div(
        style = "margin-top: 12px; padding: 10px 12px; background: #f8f9fa; border-left: 4px solid #95A5A6; border-radius: 4px;",
        p(style = "margin: 0; font-size: 0.85rem; color: #6c757d;",
          "No phenotype characteristics recorded for this position.")
      )
    }

    aa_change <- if (!is.na(ref_aa_val) && !is.na(alt_aa_val)) {
      paste0(ref_aa_val, " \u2192 ", alt_aa_val)
    } else {
      ""
    }

    bs4Dash::bs4Card(
      title = div(
        style = "display: flex; align-items: center; gap: 12px;",
        span(mut_label, style = "font-weight: 700; font-size: 1.1rem;"),
        span(paste0(protein, " position ", pos,
                    if (nchar(aa_change) > 0) paste0(" (", aa_change, ")") else ""),
             style = "color: #6c757d; font-size: 0.9rem;")
      ),
      width = 12,
      status = "primary",
      solidHeader = TRUE,
      collapsible = FALSE,

      # Row 1: Stats + Temporal + Map
      fluidRow(
        column(3,
          div(style = "padding: 8px;",
            div(style = "display: flex; gap: 16px; margin-bottom: 12px;",
              div(style = "text-align: center; flex: 1; padding: 10px; background: #f8f9fa; border-radius: 6px;",
                h4(paste0(query_prev, "%"), style = "margin: 0; color: #E74C3C; font-weight: 700;"),
                p("Query", style = "margin: 2px 0 0 0; font-size: 0.75rem; color: #6c757d;"),
                p(paste(n_query_mut, "/", n_query_total), style = "margin: 0; font-size: 0.7rem;")
              ),
              div(style = "text-align: center; flex: 1; padding: 10px; background: #f8f9fa; border-radius: 6px;",
                h4(paste0(bg_prev, "%"), style = "margin: 0; color: #6c757d; font-weight: 700;"),
                p("Background", style = "margin: 2px 0 0 0; font-size: 0.75rem; color: #6c757d;"),
                p(paste(n_bg_mut, "/", n_bg_total), style = "margin: 0; font-size: 0.7rem;")
              )
            ),
            div(style = paste0("padding: 6px 10px; border-radius: 4px; border-left: 4px solid ",
                                status_color, "; background: ", status_color, "15; margin-bottom: 8px;"),
              strong(style = paste0("color: ", status_color, "; font-size: 0.85rem;"), status)
            ),
            if (length(query_samples_list) > 0) {
              div(style = "margin-top: 8px;",
                p(strong("In query:"), style = "margin: 0 0 2px 0; font-size: 0.8rem;"),
                p(paste(query_samples_list, collapse = ", "),
                  style = "font-size: 0.8rem; color: #495057; margin: 0;")
              )
            }
          )
        ),
        column(5,
          plotOutput(mp_id(sp, "intelligence_trajectory_plot"), height = "200px")
        ),
        column(4,
          leaflet::leafletOutput(mp_id(sp, "intelligence_geo_map"), height = "200px")
        )
      ),

      # Row 2: Sequence Logo + Property Annotation
      fluidRow(
        column(4,
          plotOutput(mp_id(sp, "intelligence_seqlogo_plot"), height = "180px")
        ),
        column(8,
          uiOutput(mp_id(sp, "intelligence_aa_properties"))
        )
      ),

      # Row 3: Phenotype
      pheno_ui
    )
  })

  # --- Temporal prevalence line plot ---
  output[[mp_id(sp, "intelligence_trajectory_plot")]] <- renderPlot({
    mid <- selected_mutation_id()
    df <- mutation_detail_data()
    if (is.null(mid) || is.null(df) || nrow(df) == 0) return(NULL)

    year_totals <- df |>
      filter(!is.na(collection_year)) |>
      distinct(sample_id, collection_year) |>
      count(collection_year, name = "total")

    mut_counts <- df |>
      filter(mutation_id == mid, !is.na(collection_year)) |>
      distinct(sample_id, collection_year) |>
      count(collection_year, name = "n_mutant")

    plot_df <- year_totals |>
      left_join(mut_counts, by = "collection_year") |>
      mutate(
        n_mutant = tidyr::replace_na(n_mutant, 0L),
        prevalence = round(100 * n_mutant / total, 1)
      )

    if (nrow(plot_df) == 0) {
      plot.new()
      text(0.5, 0.5, "No temporal data.", cex = 0.9, col = "#6c757d")
      return()
    }

    mut_label <- unique(df$mutation_label[df$mutation_id == mid])[1]
    overall_prev <- round(100 * sum(plot_df$n_mutant) / sum(plot_df$total), 1)

    p <- ggplot(plot_df, aes(x = collection_year, y = prevalence))

    if (nrow(plot_df) > 1) {
      p <- p + geom_line(colour = "#E74C3C", linewidth = 1)
    }

    p <- p +
      geom_point(colour = "#E74C3C", size = 3) +
      geom_hline(yintercept = overall_prev, linetype = "dashed", colour = "#95A5A6", linewidth = 0.5) +
      scale_y_continuous(limits = c(0, max(plot_df$prevalence * 1.1, 5))) +
      labs(x = NULL, y = "Prevalence (%)",
           title = paste(mut_label, "\u2014 Temporal Trend")) +
      theme_bw(base_size = 11) +
      theme(plot.title = element_text(size = 10, face = "bold"),
            plot.margin = margin(t = 5, r = 10, b = 5, l = 5))
    p
  })

  # --- Geographic map (leaflet) ---
  output[[mp_id(sp, "intelligence_geo_map")]] <- leaflet::renderLeaflet({
    mid <- selected_mutation_id()
    df <- mutation_detail_data()
    if (is.null(mid) || is.null(df) || nrow(df) == 0) return(leaflet::leaflet() |> leaflet::addTiles())

    # Country centroid lookup (common pathogen-relevant countries)
    country_centroids <- data.frame(
      country = c("Uganda", "Democratic Republic of the Congo", "Sudan", "South Sudan",
                  "Gabon", "Republic of the Congo", "Ivory Coast", "Guinea",
                  "Liberia", "Sierra Leone", "Nigeria", "Mali", "Senegal",
                  "United States", "United Kingdom", "Spain", "Italy",
                  "China", "India", "Brazil", "South Africa", "Kenya",
                  "Tanzania", "Ethiopia", "Cameroon", "Ghana", "Angola",
                  "Mozambique", "Zimbabwe", "Zambia", "Rwanda", "Burundi"),
      lat = c(1.37, -4.04, 15.50, 6.88,
              -0.80, -0.23, 7.54, 9.95,
              6.43, 8.46, 9.08, 17.57, 14.50,
              37.09, 55.38, 40.46, 41.87,
              35.86, 20.59, -14.24, -30.56, -0.02,
              -6.37, 9.14, 7.37, 7.95, -11.20,
              -18.67, -19.02, -13.13, -1.94, -3.37),
      lng = c(32.29, 21.76, 32.53, 31.31,
              11.61, 15.83, -5.55, -9.70,
              -9.43, -11.78, 8.68, -4.00, -14.45,
              -95.71, -3.44, -3.75, 12.57,
              104.20, 78.96, -51.93, 22.94, 37.91,
              34.89, 40.49, 12.35, -1.02, 17.87,
              35.53, 29.15, 28.32, 29.87, 29.92),
      stringsAsFactors = FALSE
    )

    # Compute per-country prevalence
    country_totals <- df |>
      filter(!is.na(country), country != "") |>
      distinct(sample_id, country) |>
      count(country, name = "total")

    mut_country <- df |>
      filter(mutation_id == mid, !is.na(country), country != "") |>
      distinct(sample_id, country) |>
      count(country, name = "n_mutant")

    if (nrow(country_totals) == 0 || nrow(mut_country) == 0) {
      return(leaflet::leaflet() |> leaflet::addTiles() |>
        leaflet::setView(lng = 20, lat = 0, zoom = 2))
    }

    prev <- country_totals |>
      left_join(mut_country, by = "country") |>
      mutate(
        n_mutant = tidyr::replace_na(n_mutant, 0L),
        prevalence = round(100 * n_mutant / total, 1)
      ) |>
      filter(n_mutant > 0) |>
      inner_join(country_centroids, by = "country")

    if (nrow(prev) == 0) {
      return(leaflet::leaflet() |> leaflet::addTiles() |>
        leaflet::setView(lng = 20, lat = 0, zoom = 2))
    }

    mut_label <- unique(df$mutation_label[df$mutation_id == mid])[1]

    leaflet::leaflet(prev) |>
      leaflet::addTiles() |>
      leaflet::addCircleMarkers(
        lng = ~lng, lat = ~lat,
        radius = ~sqrt(prevalence) * 3 + 5,
        color = "#E74C3C", fillColor = "#E74C3C",
        fillOpacity = 0.6, weight = 1,
        popup = ~paste0(
          "<strong>", country, "</strong><br>",
          "Prevalence: ", prevalence, "%<br>",
          "Samples: ", n_mutant, " / ", total
        )
      ) |>
      leaflet::fitBounds(
        lng1 = min(prev$lng) - 5, lat1 = min(prev$lat) - 5,
        lng2 = max(prev$lng) + 5, lat2 = max(prev$lat) + 5
      )
  })

  # --- Amino acid dictionary: properties and biochemistry color codes ---
  .aa_dict <- list(
    G = list(code3 = "Gly", name = "Glycine",       class = "Nonpolar", subclass = "Aliphatic (smallest)", color = "#000000", charge = "Neutral", size = "Tiny", special = "Backbone flexibility; no side chain"),
    A = list(code3 = "Ala", name = "Alanine",       class = "Nonpolar", subclass = "Aliphatic",           color = "#000000", charge = "Neutral", size = "Small", special = "Helix-forming; minimal steric effects"),
    V = list(code3 = "Val", name = "Valine",        class = "Nonpolar", subclass = "Aliphatic, branched", color = "#000000", charge = "Neutral", size = "Medium", special = "Beta-sheet preference; hydrophobic core"),
    L = list(code3 = "Leu", name = "Leucine",       class = "Nonpolar", subclass = "Aliphatic, branched", color = "#000000", charge = "Neutral", size = "Large", special = "Hydrophobic core; helix-forming"),
    I = list(code3 = "Ile", name = "Isoleucine",    class = "Nonpolar", subclass = "Aliphatic, branched", color = "#000000", charge = "Neutral", size = "Large", special = "Beta-branched; restricted conformations"),
    P = list(code3 = "Pro", name = "Proline",       class = "Nonpolar", subclass = "Cyclic imino",        color = "#000000", charge = "Neutral", size = "Medium", special = "Helix breaker; introduces backbone rigidity"),
    F = list(code3 = "Phe", name = "Phenylalanine", class = "Nonpolar", subclass = "Aromatic",            color = "#000000", charge = "Neutral", size = "Large", special = "Aromatic stacking; hydrophobic core"),
    W = list(code3 = "Trp", name = "Tryptophan",    class = "Nonpolar", subclass = "Aromatic (indole)",   color = "#000000", charge = "Neutral", size = "Largest", special = "Membrane anchoring; aromatic interactions"),
    M = list(code3 = "Met", name = "Methionine",    class = "Nonpolar", subclass = "Sulfur-containing",   color = "#000000", charge = "Neutral", size = "Large", special = "Flexible hydrophobic; oxidation-sensitive"),
    S = list(code3 = "Ser", name = "Serine",        class = "Polar",    subclass = "Hydroxyl",            color = "#009900", charge = "Neutral", size = "Small", special = "H-bond donor/acceptor; phosphorylation site"),
    T = list(code3 = "Thr", name = "Threonine",     class = "Polar",    subclass = "Hydroxyl, branched",  color = "#009900", charge = "Neutral", size = "Medium", special = "H-bond donor/acceptor; glycosylation site"),
    C = list(code3 = "Cys", name = "Cysteine",      class = "Polar",    subclass = "Sulfhydryl (thiol)",  color = "#009900", charge = "Neutral", size = "Small", special = "Disulfide bonds; redox-sensitive; metal coordination"),
    Y = list(code3 = "Tyr", name = "Tyrosine",      class = "Polar",    subclass = "Aromatic, hydroxyl",  color = "#009900", charge = "Neutral", size = "Large", special = "Phosphorylation site; aromatic H-bonding"),
    N = list(code3 = "Asn", name = "Asparagine",    class = "Polar",    subclass = "Amide",               color = "#009900", charge = "Neutral", size = "Medium", special = "N-glycosylation site; H-bond donor/acceptor"),
    Q = list(code3 = "Gln", name = "Glutamine",     class = "Polar",    subclass = "Amide",               color = "#009900", charge = "Neutral", size = "Large", special = "H-bond donor/acceptor; deamidation-prone"),
    D = list(code3 = "Asp", name = "Aspartate",     class = "Acidic",   subclass = "Negatively charged",  color = "#CC0000", charge = "Negative (-1)", size = "Medium", special = "Salt bridges; metal coordination; catalytic residue"),
    E = list(code3 = "Glu", name = "Glutamate",     class = "Acidic",   subclass = "Negatively charged",  color = "#CC0000", charge = "Negative (-1)", size = "Large", special = "Salt bridges; catalytic residue; longer reach than Asp"),
    K = list(code3 = "Lys", name = "Lysine",        class = "Basic",    subclass = "Positively charged",  color = "#0000CC", charge = "Positive (+1)", size = "Large", special = "Salt bridges; ubiquitination; acetylation site"),
    R = list(code3 = "Arg", name = "Arginine",      class = "Basic",    subclass = "Positively charged (guanidinium)", color = "#0000CC", charge = "Positive (+1)", size = "Largest", special = "Strong salt bridges; RNA binding; delocalized charge"),
    H = list(code3 = "His", name = "Histidine",     class = "Basic",    subclass = "Aromatic, imidazole", color = "#0000CC", charge = "Positive (pH-dependent)", size = "Large", special = "Catalytic residue; metal coordination; pH-sensing")
  )

  .aa_get <- function(code) {
    if (is.null(code) || is.na(code) || !code %in% names(.aa_dict)) {
      return(list(code3 = "???", name = "Unknown", class = "Unknown", subclass = "Unknown",
                  color = "#6c757d", charge = "Unknown", size = "Unknown", special = ""))
    }
    .aa_dict[[code]]
  }

  # --- Sequence logo: single position, WT at bottom, alts stacked above ---
  output[[mp_id(sp, "intelligence_seqlogo_plot")]] <- renderPlot({
    mid <- selected_mutation_id()
    df <- mutation_detail_data()
    aa_freq <- aa_frequencies_data()
    if (is.null(mid) || is.null(df) || nrow(df) == 0) return(NULL)

    mut_rows <- df |> filter(mutation_id == mid)
    protein <- mut_rows$protein_name[1]
    pos <- mut_rows$position[1]
    ref_aa_val <- mut_rows$ref_aa[1]
    alt_aa_val <- mut_rows$alt_aa[1]

    if (is.null(aa_freq) || nrow(aa_freq) == 0) {
      plot.new()
      text(0.5, 0.5, "No sequence data.", cex = 0.9, col = "#6c757d")
      return()
    }

    # Get frequencies at this position
    pos_data <- aa_freq |>
      filter(protein_name == protein, position == pos) |>
      arrange(desc(frequency))

    if (nrow(pos_data) == 0) {
      plot.new()
      text(0.5, 0.5, "No data at this position.", cex = 0.9, col = "#6c757d")
      return()
    }

    # Order: WT at bottom, then alternatives sorted by frequency (smallest alt closest to WT)
    wt_row <- pos_data |> filter(amino_acid == ref_aa_val)
    alt_rows <- pos_data |> filter(amino_acid != ref_aa_val) |> arrange(frequency)

    # Stack: WT at bottom, then alts ascending frequency
    logo_df <- bind_rows(wt_row, alt_rows) |>
      mutate(
        aa_color = sapply(amino_acid, function(a) .aa_get(a)$color),
        is_wt = (amino_acid == ref_aa_val),
        ymax = cumsum(frequency),
        ymin = lag(ymax, default = 0),
        ymid = (ymin + ymax) / 2,
        height = ymax - ymin
      )

    # Use grid graphics for proper letter stretching (logo style)
    # Each letter fills full width, height = frequency proportion
    par(mar = c(2, 3, 2, 0.5))
    plot(NULL, xlim = c(0, 1), ylim = c(0, 1),
         xlab = "", ylab = "", xaxt = "n", yaxt = "n",
         main = paste0("Position ", pos), cex.main = 1.1, font.main = 2)
    axis(2, at = seq(0, 1, 0.25), labels = paste0(seq(0, 100, 25), "%"),
         las = 1, cex.axis = 0.9)

    # Draw each letter as a full-width stretched character
    for (i in seq_len(nrow(logo_df))) {
      row <- logo_df[i, ]
      if (row$height < 0.001) next

      # Letter height in user coords determines cex
      # Full plot height = 1.0 user units ~= plot height in inches
      # We want each letter to exactly fill its vertical band
      letter_cex <- row$height * 12  # scale factor: larger freq = taller letter

      text(x = 0.5, y = row$ymid, labels = row$amino_acid,
           col = row$aa_color, cex = max(letter_cex, 0.8),
           font = 2, family = "mono")
    }

    # Add frequency labels on right
    for (i in seq_len(nrow(logo_df))) {
      row <- logo_df[i, ]
      if (row$height < 0.005) next
      pct_label <- paste0(round(row$frequency * 100, 1), "%")
      text(x = 0.92, y = row$ymid, labels = pct_label,
           col = "#6c757d", cex = 0.7, adj = 0)
    }
  })

  # --- Amino acid property annotation panel ---
  output[[mp_id(sp, "intelligence_aa_properties")]] <- renderUI({
    mid <- selected_mutation_id()
    df <- mutation_detail_data()
    if (is.null(mid) || is.null(df) || nrow(df) == 0) return(NULL)

    mut_rows <- df |> filter(mutation_id == mid)
    ref_aa_val <- mut_rows$ref_aa[1]
    alt_aa_val <- mut_rows$alt_aa[1]
    mut_label <- mut_rows$mutation_label[1]

    if (is.na(ref_aa_val) || is.na(alt_aa_val)) return(NULL)

    ref_props <- .aa_get(ref_aa_val)
    alt_props <- .aa_get(alt_aa_val)

    # Determine if conservative or non-conservative
    same_class <- ref_props$class == alt_props$class
    same_charge <- ref_props$charge == alt_props$charge
    is_conservative <- same_class && same_charge

    # Generate functional implication
    implication <- if (is_conservative) {
      paste0("Conservative substitution (",  ref_props$class, " \u2192 ", alt_props$class,
             "). Side-chain properties are similar; likely minimal structural impact.")
    } else {
      # Describe the specific change
      changes <- c()
      if (ref_props$class != alt_props$class) {
        changes <- c(changes, paste0(ref_props$class, " \u2192 ", alt_props$class))
      }
      if (ref_props$charge != alt_props$charge) {
        changes <- c(changes, paste0("charge: ", ref_props$charge, " \u2192 ", alt_props$charge))
      }
      if (ref_props$size != alt_props$size) {
        changes <- c(changes, paste0("size: ", ref_props$size, " \u2192 ", alt_props$size))
      }
      paste0("Non-conservative substitution (", paste(changes, collapse = "; "),
             "). May affect protein folding, stability, or function.")
    }

    # Specific structural notes
    structural_note <- ""
    if (alt_aa_val == "P") {
      structural_note <- "Proline introduces backbone rigidity and eliminates the amide H-bond donor; can break helices and disrupt turns."
    } else if (ref_aa_val == "P") {
      structural_note <- "Loss of proline removes backbone constraint; may increase local flexibility."
    } else if (ref_aa_val == "G") {
      structural_note <- "Loss of glycine reduces backbone flexibility; the new side chain may cause steric clashes."
    } else if (alt_aa_val == "G") {
      structural_note <- "Glycine substitution removes side chain; increases flexibility but may destabilize the fold."
    } else if (ref_aa_val == "C" || alt_aa_val == "C") {
      structural_note <- "Cysteine change may affect disulfide bond formation or redox sensitivity."
    } else if (!same_charge && (ref_props$charge != "Neutral" || alt_props$charge != "Neutral")) {
      structural_note <- "Charge change may disrupt salt bridges, protein-protein interactions, or ligand binding."
    }

    # Severity indicator
    severity_color <- if (is_conservative) "#27AE60" else "#E74C3C"
    severity_icon <- if (is_conservative) "\u2713" else "\u26A0"
    severity_label <- if (is_conservative) "Conservative" else "Non-conservative"

    div(style = "padding: 8px 12px;",
      # Header
      h6(paste0(mut_label, " \u2014 Amino Acid Property Change"),
         style = "margin: 0 0 10px 0; font-weight: 700; font-size: 0.95rem;"),

      # Reference AA
      div(style = "display: flex; align-items: center; gap: 8px; margin-bottom: 8px;",
        span(ref_aa_val, style = paste0(
          "display: inline-block; width: 28px; height: 28px; line-height: 28px;",
          " text-align: center; font-weight: 700; font-size: 1rem; font-family: monospace;",
          " border-radius: 4px; color: white; background: ", ref_props$color, ";"
        )),
        div(
          p(style = "margin: 0; font-size: 0.85rem; font-weight: 600;",
            paste0(ref_props$name, " (", ref_props$code3, ") \u2014 Wild-type")),
          p(style = "margin: 0; font-size: 0.78rem; color: #6c757d;",
            paste0(ref_props$class, ", ", ref_props$subclass, " | ", ref_props$charge, " | Size: ", ref_props$size))
        )
      ),

      # Arrow
      div(style = "text-align: left; margin: 4px 0 4px 10px; font-size: 1.1rem; color: #6c757d;", "\u2193"),

      # Mutant AA
      div(style = "display: flex; align-items: center; gap: 8px; margin-bottom: 12px;",
        span(alt_aa_val, style = paste0(
          "display: inline-block; width: 28px; height: 28px; line-height: 28px;",
          " text-align: center; font-weight: 700; font-size: 1rem; font-family: monospace;",
          " border-radius: 4px; color: white; background: ", alt_props$color, ";"
        )),
        div(
          p(style = "margin: 0; font-size: 0.85rem; font-weight: 600;",
            paste0(alt_props$name, " (", alt_props$code3, ") \u2014 Mutant")),
          p(style = "margin: 0; font-size: 0.78rem; color: #6c757d;",
            paste0(alt_props$class, ", ", alt_props$subclass, " | ", alt_props$charge, " | Size: ", alt_props$size))
        )
      ),

      # Classification badge
      div(style = paste0("padding: 6px 10px; border-radius: 4px; border-left: 4px solid ",
                          severity_color, "; background: ", severity_color, "15; margin-bottom: 8px;"),
        strong(style = paste0("color: ", severity_color, "; font-size: 0.85rem;"),
               paste0(severity_icon, " ", severity_label, " substitution"))
      ),

      # Functional implication
      p(style = "margin: 4px 0; font-size: 0.82rem; color: #2C3E50;", implication),

      # Structural note (if applicable)
      if (nchar(structural_note) > 0) {
        p(style = "margin: 4px 0; font-size: 0.80rem; color: #6c757d; font-style: italic;",
          structural_note)
      }
    )
  })

  observeEvent(input[[mp_id(sp, "plsda_click")]], {
    click <- input[[mp_id(sp, "plsda_click")]]
    if (is.null(click) || is.null(click$coords_css)) {
      return()
    }
    df <- plsda_scores_data()
    if (is.null(df) || nrow(df) == 0) {
      click_info(NULL)
      return()
    }
    selected <- nearPoints(df, click, xvar = "PC1", yvar = "PC2",
                           threshold = 15, maxpoints = 1)
    if (nrow(selected) == 0) {
      click_info(NULL)
    } else {
      click_info(list(coords = click$coords_css, row = selected[1, , drop = FALSE]))
    }
  })

  output[[mp_id(sp, "plsda_score_plot")]] <- renderPlot({
    df <- plsda_scores_data()
    if (is.null(df) || nrow(df) == 0) {
      plot.new()
      text(0.5, 0.5, "No PLS-DA scores available.", cex = 1.1, col = "#6c757d")
      return()
    }

    colour_by <- input[[mp_id(sp, "plsda_colour")]]
    if (is.null(colour_by) || !colour_by %in% names(df)) {
      colour_by <- "country"
    }

    df <- df |>
      mutate(is_query = factor(is_query, levels = c("FALSE", "TRUE"),
                               labels = c("Background", "Query")))

    if (colour_by == "country") {
      df <- df |> mutate(country = ifelse(is.na(country), "Missing", as.character(country)))
    } else if (colour_by == "outbreak") {
      df <- df |> mutate(outbreak = ifelse(is.na(outbreak), "Missing", as.character(outbreak)))
    }

    df_ellipse <- df |> filter(.data[[colour_by]] != "Missing")

    p <- ggplot(df, aes(x = .data[["PC1"]], y = .data[["PC2"]])) +
      geom_point(aes(colour = .data[[colour_by]], shape = is_query),
                 size = 2.5, alpha = 0.8) +
      labs(title = paste("PLS-DA scores —", toupper(sp)),
           x = "PC1", y = "PC2",
           colour = tools::toTitleCase(colour_by),
           shape = "Group") +
      theme_bw(base_size = 12) +
      scale_shape_manual(values = c("Background" = 19, "Query" = 17)) +
      scale_colour_brewer(palette = "Set1", na.value = "grey50")

    if (nrow(df_ellipse) >= 2) {
      p <- p + stat_ellipse(
        aes(x = .data[["PC1"]], y = .data[["PC2"]],
            colour = .data[[colour_by]], group = .data[[colour_by]]),
        data = df_ellipse,
        alpha = 0.7,
        linewidth = 1.2,
        level = 0.95,
        inherit.aes = FALSE
      )
    }

    selected_id <- input[[mp_id(sp, "plsda_sample_search")]]
    if (!is.null(selected_id) && selected_id != "") {
      df_selected <- df |> filter(.data[["sample_id"]] == as.integer(selected_id))
      if (nrow(df_selected) > 0) {
        p <- p + geom_point(
          data = df_selected,
          aes(x = .data[["PC1"]], y = .data[["PC2"]]),
          colour = "black", size = 5, shape = 1, stroke = 1.2
        )
      }
    }

    p
  })

  output[[mp_id(sp, "plsda_click_info")]] <- renderUI({
    info <- click_info()
    if (is.null(info)) {
      return(div(style = "display: none;"))
    }

    s <- info$row[1, ]
    is_query_label <- if (isTRUE(s$is_query)) "Query" else "Background"
    year <- if (is.na(s$collection_year)) "Missing" else s$collection_year
    country <- if (is.na(s$country)) "Missing" else s$country
    outbreak <- if (is.na(s$outbreak)) "Missing" else s$outbreak
    px_x <- info$coords$x
    px_y <- info$coords$y

    div(
      style = paste0(
        "position: absolute; left: ", px_x, "px; top: ", px_y, "px;",
        "z-index: 1000; min-width: 180px; max-width: 260px;",
        "padding: 10px; background-color: #ffffff; border: 1px solid #dee2e6;",
        "border-radius: 4px; box-shadow: 0 4px 8px rgba(0,0,0,0.15);"
      ),
      h5(s$sample_name, style = "margin: 0 0 6px 0; font-weight: 600;"),
      p(style = "margin: 2px 0; font-size: 0.85rem;", strong("Group: "), is_query_label),
      p(style = "margin: 2px 0; font-size: 0.85rem;", strong("Year: "), year),
      p(style = "margin: 2px 0; font-size: 0.85rem;", strong("Country: "), country),
      p(style = "margin: 2px 0; font-size: 0.85rem;", strong("Outbreak: "), outbreak),
      p(style = "margin: 2px 0; font-size: 0.85rem;", strong("Mutation burden: "), s$mutation_burden)
    )
  })

  output[[mp_id(sp, "plsda_search_info")]] <- renderUI({
    selected_id <- input[[mp_id(sp, "plsda_sample_search")]]
    if (is.null(selected_id) || selected_id == "") {
      return(div(style = "margin-top: 10px; color: #6c757d; font-size: 0.85rem;", "Select a sample to see details."))
    }

    df <- plsda_scores_data()
    if (is.null(df) || nrow(df) == 0) return(NULL)

    selected <- df |> filter(.data[["sample_id"]] == as.integer(selected_id))
    if (nrow(selected) == 0) {
      return(div(style = "margin-top: 10px; color: #6c757d;", "Select a sample to see details."))
    }

    s <- selected[1, ]
    is_query_label <- if (isTRUE(s$is_query)) "Query" else "Background"
    year <- if (is.na(s$collection_year)) "Missing" else s$collection_year
    country <- if (is.na(s$country)) "Missing" else s$country
    outbreak <- if (is.na(s$outbreak)) "Missing" else s$outbreak

    div(
      style = "margin-top: 10px; padding: 10px; background-color: #f8f9fa; border-radius: 4px;",
      h5(s$sample_name, style = "margin: 0 0 6px 0; font-weight: 600;"),
      p(style = "margin: 2px 0; font-size: 0.85rem;", strong("Group: "), is_query_label),
      p(style = "margin: 2px 0; font-size: 0.85rem;", strong("Year: "), year),
      p(style = "margin: 2px 0; font-size: 0.85rem;", strong("Country: "), country),
      p(style = "margin: 2px 0; font-size: 0.85rem;", strong("Outbreak: "), outbreak),
      p(style = "margin: 2px 0; font-size: 0.85rem;", strong("Mutation burden: "), s$mutation_burden)
    )
  })

}
