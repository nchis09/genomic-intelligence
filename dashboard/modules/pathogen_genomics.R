# Pathogen Genomics module
#
# Species-tabbed section for phylogenetic tree visualization.
# Provides an interactive tree explorer using treeio + ggtree
# when the knowledge warehouse DuckDB export is available.

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
pg_id <- function(species, suffix) {
  paste0("pg_", species, "_", suffix)
}

# Path to the DuckDB knowledge warehouse file
pg_duckdb_path <- function(outdir) {
  file.path(outdir, "knowledge_warehouse", "knowledge_warehouse.duckdb")
}

# Read tree data from DuckDB for interactive visualization.
# Returns list(trees=, tips=) on success, or list(error="...") on failure.
pg_read_tree_data <- function(duckdb_path, species) {
  pg_err <- function(msg) { message("[pg] ", msg); list(error = msg) }

  if (!file.exists(duckdb_path)) {
    return(pg_err(paste0("DuckDB file not found:\n", duckdb_path)))
  }
  if (!requireNamespace("DBI", quietly = TRUE)) {
    return(pg_err("R package 'DBI' is not installed.\nActivate the pgirl_dashboard conda env."))
  }
  if (!requireNamespace("duckdb", quietly = TRUE)) {
    return(pg_err("R package 'duckdb' is not installed.\nActivate the pgirl_dashboard conda env."))
  }

  con <- tryCatch(
    DBI::dbConnect(duckdb::duckdb(), dbdir = duckdb_path, read_only = TRUE),
    error = function(e) e
  )
  if (inherits(con, "error")) {
    return(pg_err(paste0("Cannot open DuckDB:\n", con$message)))
  }
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)

  # Verify required tables exist
  tables <- DBI::dbListTables(con)
  missing <- setdiff(c("phylogenetic_trees", "tree_tips"), tables)
  if (length(missing) > 0) {
    return(pg_err(paste0(
      "DuckDB is missing table(s): ", paste(missing, collapse = ", "),
      ".\nRe-run the pipeline so EXPORT_KNOWLEDGE_DB writes them."
    )))
  }

  # Parameterized queries (avoid SQL injection / quoting issues)
  sp_lower <- tolower(species)

  trees <- tryCatch(
    DBI::dbGetQuery(con,
      "SELECT tree_id, tree_method, newick FROM phylogenetic_trees WHERE LOWER(species) = ?",
      params = list(sp_lower)
    ),
    error = function(e) e
  )
  if (inherits(trees, "error")) {
    return(pg_err(paste0("Query phylogenetic_trees failed:\n", trees$message)))
  }
  if (nrow(trees) == 0) {
    return(pg_err(paste0(
      "No trees found for species '", species, "' in the database.\n",
      "Available species: ",
      paste(tryCatch(
        DBI::dbGetQuery(con, "SELECT DISTINCT species FROM phylogenetic_trees")$species,
        error = function(e) "(query failed)"
      ), collapse = ", ")
    )))
  }

  tips <- tryCatch(
    DBI::dbGetQuery(con,
      "SELECT t.* FROM tree_tips t
       JOIN phylogenetic_trees pt ON t.tree_id = pt.tree_id
       WHERE LOWER(pt.species) = ?",
      params = list(sp_lower)
    ),
    error = function(e) e
  )
  if (inherits(tips, "error")) {
    return(pg_err(paste0("Query tree_tips failed:\n", tips$message)))
  }

  message("[pg] Loaded ", nrow(trees), " trees and ", nrow(tips), " tips for ", species)
  list(trees = trees, tips = tips)
}

# ---------------------------------------------------------------------------
# UI for one species
# ---------------------------------------------------------------------------
pathogen_genomics_ui <- function(species) {
  tagList(
    h3("Pathogen Genomics"),
    p(style = "color: #6c757d; margin-bottom: 16px;",
      "Interactive phylogenetic tree explorer for the detected pathogen."),

    fluidRow(
      column(3,
        bs4Dash::bs4Card(
          title = "Annotation Layers",
          width = 12,
          status = "secondary",
          solidHeader = FALSE,
          checkboxGroupInput(
            inputId = pg_id(species, "annotations"),
            label = "Show layers:",
            choices = c(
              "Clade coloring" = "clade",
              "Outbreak labels" = "outbreak",
              "Country" = "country",
              "Genome coverage" = "coverage",
              "Mutation count" = "mutations",
              "Query highlight" = "query",
              "Bootstrap support" = "bootstrap"
            ),
            selected = c("clade", "query", "outbreak")
          ),
          radioButtons(
            inputId = pg_id(species, "tree_select"),
            label = "Tree type:",
            choices = c("Augur (evolutionary)" = "augur", "IQ-TREE (bootstrap)" = "iqtree"),
            selected = "augur"
          ),
          radioButtons(
            inputId = pg_id(species, "layout_select"),
            label = "Layout:",
            choices = c("Rectangular" = "rectangular", "Circular" = "circular"),
            selected = "rectangular"
          ),
          hr(),
          sliderInput(
            inputId = pg_id(species, "tip_height"),
            label = "Tip spacing (px):",
            min = 6, max = 30, value = 12, step = 2
          ),
          sliderInput(
            inputId = pg_id(species, "zoom"),
            label = "Horizontal zoom:",
            min = 0.5, max = 3, value = 1, step = 0.1
          )
        )
      ),
      column(9,
        bs4Dash::bs4Card(
          title = paste0("Phylogenetic Tree \u2014 ", toupper(species)),
          width = 12,
          status = "primary",
          solidHeader = TRUE,
          div(style = "min-height: 500px; overflow-y: auto;",
            uiOutput(pg_id(species, "interactive_tree_ui"))
          )
        )
      )
    ),
    fluidRow(
      column(12,
        bs4Dash::bs4Card(
          title = "Tip Details",
          width = 12,
          status = "secondary",
          solidHeader = FALSE,
          collapsed = TRUE,
          collapsible = TRUE,
          DT::DTOutput(pg_id(species, "tip_table"))
        )
      )
    )
  )
}

# ---------------------------------------------------------------------------
# Server registration for one species
# ---------------------------------------------------------------------------
pathogen_genomics_register <- function(input, output, session, species, outdir) {
  sp <- species

  # -- Interactive tree (requires treeio + ggtree + ape) --
  tree_data_rv <- reactive({
    duckdb_path <- pg_duckdb_path(outdir())
    tryCatch(
      pg_read_tree_data(duckdb_path, sp),
      error = function(e) list(error = paste0("Reactive error: ", e$message))
    )
  })

  # Dynamic plot height: scale with number of tips so annotations stay readable
  output[[pg_id(sp, "interactive_tree_ui")]] <- renderUI({
    data <- tree_data_rv()
    n_tips <- 30  # default
    if (!is.null(data$trees) && is.null(data$error)) {
      newick <- data$trees$newick[1]
      if (!is.null(newick) && !is.na(newick) && nchar(newick) > 0) {
        tmp <- tempfile(fileext = ".nwk")
        writeLines(newick, tmp)
        tr <- tryCatch(ape::read.tree(tmp), error = function(e) NULL)
        unlink(tmp)
        if (!is.null(tr)) n_tips <- length(tr$tip.label)
      }
    }
    tip_h <- input[[pg_id(sp, "tip_height")]]
    if (is.null(tip_h)) tip_h <- 12
    plot_h <- max(500, n_tips * tip_h)
    plotOutput(pg_id(sp, "interactive_tree"), height = paste0(plot_h, "px"))
  })

  output[[pg_id(sp, "interactive_tree")]] <- renderPlot({
    tryCatch({
    data <- tree_data_rv()
    if (!is.null(data$error)) {
      plot.new()
      text(0.5, 0.5, paste0("No tree data available.\n\n", data$error),
           cex = 1.1, col = "#6c757d")
      return()
    }
    if (is.null(data$trees)) {
      plot.new()
      text(0.5, 0.5, "No tree data returned (unexpected NULL).",
           cex = 1.1, col = "#6c757d")
      return()
    }

    # Check for required packages
    pkgs_available <- all(sapply(c("ape", "treeio", "ggtree"), requireNamespace, quietly = TRUE))
    if (!pkgs_available) {
      plot.new()
      text(0.5, 0.5, "Interactive tree requires: ape, treeio, ggtree\nInstall via BiocManager or conda.",
           cex = 1.1, col = "#6c757d")
      return()
    }

    library(ape)
    library(ggtree)

    # Select tree based on user choice
    tree_choice <- input[[pg_id(sp, "tree_select")]]
    trees_df <- data$trees
    tips_df <- data$tips

    selected_tree <- if (!is.null(tree_choice) && tree_choice == "iqtree") {
      trees_df[grepl("iqtree", tolower(trees_df$tree_method)), ]
    } else {
      trees_df[!grepl("iqtree", tolower(trees_df$tree_method)), ]
    }

    if (nrow(selected_tree) == 0) {
      selected_tree <- trees_df[1, , drop = FALSE]
    } else {
      selected_tree <- selected_tree[1, , drop = FALSE]
    }

    newick <- selected_tree$newick
    tree_id <- selected_tree$tree_id

    if (is.na(newick) || nchar(trimws(newick)) == 0) {
      plot.new()
      text(0.5, 0.5, "Tree newick is empty.", cex = 1.2, col = "#6c757d")
      return()
    }

    # Parse newick
    tmp <- tempfile(fileext = ".nwk")
    writeLines(newick, tmp)
    tr <- tryCatch(ape::read.tree(tmp), error = function(e) NULL)
    unlink(tmp)

    if (is.null(tr)) {
      plot.new()
      text(0.5, 0.5, "Failed to parse tree newick.", cex = 1.2, col = "#6c757d")
      return()
    }

    # Normalize tip labels: replace Newick metacharacters and spaces with _
    # so they match the _clean_newick_name() format used by the database.
    tr$tip.label <- gsub("[,;():\\[\\] ]", "_", tr$tip.label)

    # Ladderize for cleaner layout, and truncate extreme branch lengths
    # so one long branch doesn't compress all other tips into a sliver.
    tr <- ape::ladderize(tr)
    if (!is.null(tr$edge.length) && length(tr$edge.length) > 0) {
      q95 <- quantile(tr$edge.length[tr$edge.length > 0], 0.95, na.rm = TRUE)
      tr$edge.length[tr$edge.length > 3 * q95] <- 3 * q95
    }

    # Get annotations selection
    annotations <- input[[pg_id(sp, "annotations")]]
    layout_choice <- input[[pg_id(sp, "layout_select")]]
    if (is.null(layout_choice)) layout_choice <- "rectangular"

    # Dark mode
    dark <- isTRUE(input$is_dark_mode)
    bg_col <- if (dark) "#1e1e2e" else "white"
    text_col <- if (dark) "#e0e0e0" else "#333333"
    branch_col <- if (dark) "#8899aa" else "#555555"

    # Base tree plot
    p <- ggtree(tr, layout = layout_choice, size = 0.4, color = branch_col) +
      theme(
        plot.background = element_rect(fill = bg_col, color = NA),
        panel.background = element_rect(fill = bg_col, color = NA),
        legend.background = element_rect(fill = bg_col, color = NA),
        legend.text = element_text(color = text_col, size = 8),
        legend.title = element_text(color = text_col, size = 9, face = "bold")
      )

    # Expand x-axis based on zoom slider so tip labels/annotations have room
    zoom_factor <- input[[pg_id(sp, "zoom")]]
    if (is.null(zoom_factor)) zoom_factor <- 1
    x_range <- layer_scales(p)$x$range$range
    if (!is.null(x_range) && length(x_range) == 2) {
      expand <- 1 + (0.4 * zoom_factor)
      p <- p + xlim(NA, x_range[2] * expand)
    }

    # Join tip metadata — prefix columns to avoid ggtree name collisions.
    # If the selected tree has no tips (e.g. IQ-TREE trees don't store tips
    # separately), fall back to tips from another tree of the same species.
    tip_meta <- tips_df[tips_df$tree_id == tree_id, ]
    if (nrow(tip_meta) == 0 && nrow(tips_df) > 0) {
      tip_meta <- tips_df[!duplicated(tips_df$label), ]
    }
    keep_cols <- intersect(
      c("label", "is_query", "clade", "outbreak", "country",
        "genome_coverage", "nuc_mutation_count", "aa_mutation_count"),
      names(tip_meta)
    )
    tip_meta <- tip_meta[, keep_cols, drop = FALSE]
    # Also normalize database labels to match the tree tip cleaning
    tip_meta$label <- gsub("[,;():\\[\\] ]", "_", tip_meta$label)
    tip_meta <- tip_meta[!duplicated(tip_meta$label), ]

    # Log match rate for debugging
    tree_labels <- tr$tip.label
    matched <- sum(tree_labels %in% tip_meta$label)
    message(sprintf("[pg] Tip join: %d/%d tree labels matched (%s, tree_id=%s)",
                    matched, length(tree_labels), sp, tree_id))

    tip_meta$clade[is.na(tip_meta$clade) | tip_meta$clade == ""] <- "Unknown"
    tip_meta$outbreak[is.na(tip_meta$outbreak) | tip_meta$outbreak == ""] <- "Unknown"
    tip_meta$country[is.na(tip_meta$country) | tip_meta$country == ""] <- "Unknown"
    tip_meta$is_query[is.na(tip_meta$is_query)] <- FALSE
    tip_meta$genome_coverage[is.na(tip_meta$genome_coverage)] <- 0
    tip_meta$nuc_mutation_count[is.na(tip_meta$nuc_mutation_count)] <- 0

    if (nrow(tip_meta) > 0) {
      # Prefix non-label columns to avoid collision with ggtree internal data
      ann <- tip_meta
      data_cols <- setdiff(names(ann), "label")
      names(ann)[names(ann) %in% data_cols] <- paste0("pg_", data_cols)
      p <- p %<+% ann

      # Clade coloring
      if ("clade" %in% annotations) {
        n_clades <- length(unique(tip_meta$clade))
        clade_cols <- if (n_clades <= 8) {
          setNames(RColorBrewer::brewer.pal(max(3, n_clades), "Set2")[1:n_clades],
                   sort(unique(tip_meta$clade)))
        } else {
          setNames(colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(n_clades),
                   sort(unique(tip_meta$clade)))
        }
        p <- p + geom_tippoint(aes(color = pg_clade), size = 2.5, alpha = 0.8) +
          scale_color_manual(values = clade_cols, name = "Clade", na.value = "#999999")
      }

      # Outbreak labels
      if ("outbreak" %in% annotations && "pg_outbreak" %in% names(p$data)) {
        outbreak_vals <- unique(tip_meta$outbreak)
        n_ob <- length(outbreak_vals)
        ob_cols <- if (n_ob <= 8) {
          setNames(RColorBrewer::brewer.pal(max(3, n_ob), "Dark2")[1:n_ob],
                   sort(outbreak_vals))
        } else {
          setNames(colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(n_ob),
                   sort(outbreak_vals))
        }
        p <- p + ggnewscale::new_scale_fill() +
          geom_tippoint(aes(fill = pg_outbreak), shape = 21, size = 2.5,
                        color = "transparent", alpha = 0.8) +
          scale_fill_manual(values = ob_cols, name = "Outbreak", na.value = "#999999")
      }

      # Country
      if ("country" %in% annotations && "pg_country" %in% names(p$data)) {
        p <- p + geom_tiplab(
          aes(label = pg_country),
          size = 2, offset = 0.001, align = FALSE,
          color = if (dark) "#90caf9" else "#1565c0",
          fontface = "italic"
        )
      }

      # Genome coverage bar at tips
      if ("coverage" %in% annotations && "pg_genome_coverage" %in% names(p$data)) {
        p <- p + geom_tippoint(
          aes(size = pg_genome_coverage),
          shape = 16, alpha = 0.5,
          color = if (dark) "#4dd0e1" else "#00838f"
        ) +
          scale_size_continuous(name = "Coverage", range = c(1, 5))
      }

      # Mutation count
      if ("mutations" %in% annotations && "pg_nuc_mutation_count" %in% names(p$data)) {
        p <- p + ggnewscale::new_scale_color() +
          geom_tippoint(
            aes(color = pg_nuc_mutation_count),
            shape = 18, size = 2.5, alpha = 0.8
          ) +
          scale_color_gradient(low = if (dark) "#a5d6a7" else "#c8e6c9",
                               high = if (dark) "#e53935" else "#b71c1c",
                               name = "Nuc. mutations")
      }

      # Query highlight — use the prefixed pg_is_query column directly
      # instead of matching labels with !!, which fails when a label
      # like "SAMPLE" collides with an R symbol name.
      if ("query" %in% annotations && "pg_is_query" %in% names(p$data)) {
        p <- p + geom_tiplab(
          aes(subset = (pg_is_query == TRUE)),
          size = 2.5, offset = 0.0005,
          color = if (dark) "#ffab40" else "#e65100",
          fontface = "bold"
        )
      }

      # Bootstrap support (IQ-TREE only)
      if ("bootstrap" %in% annotations && grepl("iqtree", tolower(selected_tree$tree_method))) {
        n_tips <- length(tr$tip.label)
        if (!is.null(tr$node.label) && length(tr$node.label) > 0) {
          bs_vals <- suppressWarnings(as.numeric(tr$node.label))
          bs_vals[is.na(bs_vals)] <- 0
          high_support <- (n_tips + 1):(n_tips + tr$Nnode)
          high_support <- high_support[bs_vals >= 70]
          if (length(high_support) > 0) {
            p <- p + geom_point2(
              aes(subset = (node %in% !!high_support)),
              color = if (dark) "#66bb6a" else "#2e7d32",
              size = 2, alpha = 0.7, shape = 16
            )
          }
        }
      }
    }

    # Force print inside tryCatch so ggplot rendering errors are caught here
    print(p)
    }, error = function(e) {
      plot.new()
      text(0.5, 0.5, paste0("Render error:\n\n", conditionMessage(e)),
           cex = 1.0, col = "#c0392b")
    })
  })

  # -- Tip metadata table --
  output[[pg_id(sp, "tip_table")]] <- DT::renderDT({
    data <- tree_data_rv()
    if (!is.null(data$error) || is.null(data$tips)) return(DT::datatable(data.frame()))

    tips_df <- data$tips
    display_cols <- intersect(
      c("label", "is_query", "clade", "outbreak", "country", "div",
        "genome_coverage", "nextclade_qc", "nuc_mutation_count", "aa_mutation_count"),
      names(tips_df)
    )
    df <- tips_df[, display_cols, drop = FALSE]
    df <- df[!duplicated(df$label), ]

    DT::datatable(
      df,
      rownames = FALSE,
      options = list(
        pageLength = 15,
        scrollX = TRUE,
        dom = "frtip"
      ),
      class = "compact stripe"
    )
  })
}
