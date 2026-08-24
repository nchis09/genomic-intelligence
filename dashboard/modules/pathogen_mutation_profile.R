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
      column(6,
        bs4Dash::bs4Card(
          title = "Per-Protein Mutation Positions (Background vs Query)",
          width = 12,
          status = "warning",
          solidHeader = TRUE,
          plotOutput(mp_id(species, "protein_burden_summary_plot"), height = "400px")
        )
      ),
      column(6,
        bs4Dash::bs4Card(
          title = "Per-Protein Mutation Positions by Sample",
          width = 12,
          status = "warning",
          solidHeader = TRUE,
          DT::DTOutput(mp_id(species, "protein_burden_samples_table"))
        )
      )
    ),

    fluidRow(
      column(6,
        bs4Dash::bs4Card(
          title = "PLS-DA Loadings (PC1 vs PC2)",
          width = 12,
          status = "success",
          solidHeader = TRUE,
          plotOutput(mp_id(species, "plsda_loading_plot"))
        )
      ),
      column(6,
        bs4Dash::bs4Card(
          title = "PLS-DA VIP",
          width = 12,
          status = "success",
          solidHeader = TRUE,
          plotOutput(mp_id(species, "plsda_vip_plot"))
        )
      )
    ),

    fluidRow(
      column(12,
        bs4Dash::bs4Card(
          title = "PLS-DA Scores Table",
          width = 12,
          status = "success",
          solidHeader = TRUE,
          DT::DTOutput(mp_id(species, "plsda_scores_table"))
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

  plsda_loadings_data <- reactive({
    mp_read_table(mp_table_path(outdir(), sp, "01_plsda_loadings.tsv"))
  })

  plsda_vip_data <- reactive({
    mp_read_table(mp_table_path(outdir(), sp, "01_plsda_vip.tsv"))
  })

  click_info <- reactiveVal(NULL)

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
      labs(x = "Protein", y = "Mutation positions per sample") +
      theme_bw(base_size = 12) +
      theme(legend.position = "none")

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

  output[[mp_id(sp, "protein_burden_samples_table")]] <- DT::renderDT({
    req(protein_samples_data())
    DT::datatable(
      protein_samples_data(),
      options = list(pageLength = 15, scrollX = TRUE),
      rownames = FALSE
    ) |>
      DT::formatRound(columns = c("mutation_positions"), digits = 0)
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

  output[[mp_id(sp, "plsda_loading_plot")]] <- renderPlot({
    df <- plsda_loadings_data()
    if (is.null(df) || nrow(df) == 0) {
      plot.new()
      text(0.5, 0.5, "No PLS-DA loadings available.", cex = 1.1, col = "#6c757d")
      return()
    }

    ggplot(df, aes(x = .data[["PC1"]], y = .data[["PC2"]], label = .data[["protein_name"]])) +
      geom_point(alpha = 0.7, size = 3, colour = "#2C3E50") +
      geom_text(vjust = -0.5, hjust = 0.5, size = 3, colour = "#2C3E50") +
      labs(title = paste("PLS-DA loadings —", toupper(sp)),
           x = "PC1", y = "PC2") +
      theme_bw(base_size = 12)
  })

  output[[mp_id(sp, "plsda_vip_plot")]] <- renderPlot({
    df <- plsda_vip_data()
    if (is.null(df) || nrow(df) == 0) {
      plot.new()
      text(0.5, 0.5, "No PLS-DA VIP available.", cex = 1.1, col = "#6c757d")
      return()
    }

    df <- df |> arrange(desc(.data[["vip"]])) |> head(20)

    ggplot(df, aes(x = reorder(.data[["protein_name"]], .data[["vip"]]), y = .data[["vip"]])) +
      geom_col(fill = "#27AE60") +
      coord_flip() +
      labs(title = paste("Top 20 PLS-DA VIP —", toupper(sp)),
           x = "Protein", y = "VIP (|PC1| + |PC2|)") +
      theme_bw(base_size = 12)
  })

  output[[mp_id(sp, "plsda_scores_table")]] <- DT::renderDT({
    req(plsda_scores_data())
    DT::datatable(
      plsda_scores_data(),
      options = list(pageLength = 10, scrollX = TRUE),
      rownames = FALSE
    ) |>
      DT::formatRound(columns = c("PC1", "PC2"), digits = 3)
  })
}
