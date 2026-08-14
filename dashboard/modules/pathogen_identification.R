# Biological Threat module (formerly "Pathogen Identification")
#
# Renders the 4 curated tables/plots produced by bin/pathogen_identification.R
# for a single species, and wires up the corresponding server-side renderers.
# Not a true Shiny module (no NS()) because tabs are generated dynamically per
# species discovered under --outdir; output IDs are namespaced manually with
# the species name instead. Internal `pi_`/`pi_<species>` prefixes are kept
# as-is (purely internal ids) even though the user-facing label is now
# "Biological Threat".

PI_TABLES <- list(
  species_assignment = list(
    title    = "Species Assignment",
    filename = "00_species_assignment.tsv",
    plot_col = "coverage",
    plot_lab = "Coverage"
  ),
  sequence_similarity = list(
    title    = "Sequence Similarity",
    filename = "02_sequence_similarity.tsv",
    plot_col = "percent_nucleotide_identity",
    plot_lab = "Percent Nucleotide Identity"
  ),
  msa_profile = list(
    title    = "MSA Profile",
    filename = "03_msa_profile.tsv",
    plot_col = "mean_identity",
    plot_lab = "Mean Identity"
  ),
  phylogenetic_placement = list(
    title    = "Phylogenetic Placement",
    filename = "04_phylogenetic_placement.tsv",
    plot_col = "distance_to_assigned_clade",
    plot_lab = "Distance to Assigned Clade"
  )
)

pi_table_path <- function(outdir, species, filename) {
  file.path(outdir, "pathogen_identification", species, "species_identification", filename)
}

pi_read_table_safe <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(readr::read_tsv(path, show_col_types = FALSE), error = function(e) NULL)
}

pi_output_id <- function(species, key, suffix) {
  paste0("pi_", species, "_", key, "_", suffix)
}

# UI for one species' Pathogen Identification content: one bs4Card per table.
pathogen_identification_ui <- function(species) {
  cards <- lapply(names(PI_TABLES), function(key) {
    tbl <- PI_TABLES[[key]]
    bs4Dash::bs4Card(
      title = tbl$title,
      width = 12,
      collapsible = TRUE,
      DT::DTOutput(pi_output_id(species, key, "table")),
      plotOutput(pi_output_id(species, key, "plot"), height = "280px")
    )
  })
  tagList(
    h3(paste0("Biological Threat \u2014 ", toupper(species))),
    cards
  )
}

# Registers renderDT/renderPlot outputs for one species. `outdir` is a
# reactive expression (e.g. reactive(input$outdir)).
pathogen_identification_register <- function(output, species, outdir) {
  for (key in names(PI_TABLES)) {
    local({
      key_local <- key
      species_local <- species
      tbl <- PI_TABLES[[key_local]]

      data_reactive <- reactive({
        path <- pi_table_path(outdir(), species_local, tbl$filename)
        df <- pi_read_table_safe(path)
        validate(need(!is.null(df), paste0("No data found at: ", path)))
        df
      })

      output[[pi_output_id(species_local, key_local, "table")]] <- DT::renderDT({
        DT::datatable(data_reactive(), options = list(pageLength = 15, scrollX = TRUE))
      })

      output[[pi_output_id(species_local, key_local, "plot")]] <- renderPlot({
        df <- data_reactive()
        if (!tbl$plot_col %in% names(df) || !"sample" %in% names(df)) return(NULL)
        df <- df %>% filter(!is.na(.data[[tbl$plot_col]]))
        validate(need(nrow(df) > 0, "No numeric data to plot."))
        ggplot(df, aes(x = sample, y = .data[[tbl$plot_col]])) +
          geom_col(fill = "#4A6C8C") +
          labs(x = "Sample", y = tbl$plot_lab) +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
      })
    })
  }
}
