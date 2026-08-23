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
      column(6,
        bs4Dash::bs4Card(
          title = "Mutation Burden",
          width = 12,
          status = "primary",
          solidHeader = TRUE,
          plotOutput(mp_id(species, "burden_boxplot"))
        )
      ),
      column(6,
        bs4Dash::bs4Card(
          title = "Group Summary",
          width = 12,
          status = "primary",
          solidHeader = TRUE,
          DT::DTOutput(mp_id(species, "burden_summary_table"))
        )
      )
    ),

    fluidRow(
      column(12,
        bs4Dash::bs4Card(
          title = "Per-Sample Burden",
          width = 12,
          status = "info",
          solidHeader = TRUE,
          DT::DTOutput(mp_id(species, "burden_samples_table"))
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
          title = "Per-Protein Burden Summary",
          width = 12,
          status = "warning",
          solidHeader = TRUE,
          DT::DTOutput(mp_id(species, "protein_burden_summary_table"))
        )
      ),
      column(6,
        bs4Dash::bs4Card(
          title = "Per-Protein Burden by Sample",
          width = 12,
          status = "warning",
          solidHeader = TRUE,
          DT::DTOutput(mp_id(species, "protein_burden_samples_table"))
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

  samples_data <- reactive({
    mp_read_table(mp_table_path(outdir(), sp, "01_mutation_burden_samples.tsv"))
  })

  summary_data <- reactive({
    mp_read_table(mp_table_path(outdir(), sp, "01_mutation_burden_summary.tsv"))
  })

  manova_data <- reactive({
    mp_read_table(mp_table_path(outdir(), sp, "01_mutation_profile_manova.tsv"))
  })

  protein_samples_data <- reactive({
    mp_read_table(mp_table_path(outdir(), sp, "01_protein_burden_samples.tsv"))
  })

  protein_summary_data <- reactive({
    mp_read_table(mp_table_path(outdir(), sp, "01_protein_burden_summary.tsv"))
  })

  output[[mp_id(sp, "burden_boxplot")]] <- renderPlot({
    df <- samples_data()
    if (is.null(df) || nrow(df) == 0) {
      plot.new()
      text(0.5, 0.5, "No mutation burden data available.", cex = 1.1, col = "#6c757d")
      return()
    }

    df <- df |>
      mutate(is_query = factor(is_query, levels = c("FALSE", "TRUE"),
                               labels = c("Background", "Query")))

    ggplot(df, aes(x = is_query, y = mutation_burden, fill = is_query)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.7) +
      geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
      scale_fill_manual(values = c("Background" = "#3498DB", "Query" = "#E74C3C")) +
      labs(title = paste("Mutation burden per sample —", toupper(sp)),
           x = "Group", y = "Mutation count", fill = NULL) +
      theme_bw(base_size = 12) +
      theme(legend.position = "none")
  })

  output[[mp_id(sp, "burden_samples_table")]] <- DT::renderDT({
    req(samples_data())
    DT::datatable(
      samples_data(),
      options = list(pageLength = 15, scrollX = TRUE),
      rownames = FALSE,
      filter = "top"
    ) |>
      DT::formatRound(columns = c("mutation_burden"), digits = 0)
  })

  output[[mp_id(sp, "burden_summary_table")]] <- DT::renderDT({
    req(summary_data())
    DT::datatable(
      summary_data(),
      options = list(pageLength = 5, scrollX = TRUE),
      rownames = FALSE
    ) |>
      DT::formatRound(columns = c("mean_burden", "median_burden", "sd_burden"), digits = 2) |>
      DT::formatSignif(columns = c("statistic", "p_value", "estimate"), digits = 3)
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

  output[[mp_id(sp, "protein_burden_summary_table")]] <- DT::renderDT({
    req(protein_summary_data())
    DT::datatable(
      protein_summary_data(),
      options = list(pageLength = 10, scrollX = TRUE),
      rownames = FALSE
    ) |>
      DT::formatRound(columns = c("mean_burden", "median_burden", "sd_burden"), digits = 2) |>
      DT::formatSignif(columns = c("statistic", "p_value", "estimate"), digits = 3)
  })

  output[[mp_id(sp, "protein_burden_samples_table")]] <- DT::renderDT({
    req(protein_samples_data())
    DT::datatable(
      protein_samples_data(),
      options = list(pageLength = 15, scrollX = TRUE),
      rownames = FALSE
    ) |>
      DT::formatRound(columns = c("mutation_count"), digits = 0)
  })
}
