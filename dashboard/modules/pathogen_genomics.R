# Pathogen Genomics module
#
# Species-tabbed section for genomic characterization data.
# Currently an empty shell; content will be populated incrementally.

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
pg_id <- function(species, suffix) {
  paste0("pg_", species, "_", suffix)
}

# ---------------------------------------------------------------------------
# UI for one species
# ---------------------------------------------------------------------------
pathogen_genomics_ui <- function(species) {
  tagList(
    h3("Pathogen Genomics"),
    p(style = "color: #6c757d; margin-bottom: 16px;",
      "Genomic characterization and molecular epidemiology for the detected pathogen."),

    bs4Dash::bs4Card(
      title = paste0("Genomic Characterization \u2014 ", toupper(species)),
      width = 12,
      status = "primary",
      solidHeader = TRUE,
      div(
        style = "text-align: center; padding: 48px 0; color: #6c757d;",
        icon("microscope", class = "fa-3x"),
        h4("Content coming soon", style = "margin-top: 16px;"),
        p(style = "max-width: 560px; margin: 8px auto;",
          "This section will display genomic characterization data including ",
          "mutation profiles, clade assignments, sequence quality metrics, ",
          "and phylogenetic context for the detected pathogen.")
      )
    )
  )
}

# ---------------------------------------------------------------------------
# Server registration for one species (placeholder — no reactives yet)
# ---------------------------------------------------------------------------
pathogen_genomics_register <- function(input, output, session, species, outdir) {
  # Future: register reactive data sources and outputs here
}
