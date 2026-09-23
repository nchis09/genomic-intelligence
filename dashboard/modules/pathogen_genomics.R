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
              "Query highlight" = "query"
            ),
            selected = c("clade", "query", "outbreak")
          ),
          radioButtons(
            inputId = pg_id(species, "tree_select"),
            label = "Tree type:",
            choices = c("Augur (evolutionary)" = "augur"),
            selected = "augur"
          ),
          radioButtons(
            inputId = pg_id(species, "layout_select"),
            label = "Layout:",
            choices = c("Rectangular" = "rectangular", "Circular" = "circular"),
            selected = "rectangular"
          ),
          hr(),
          checkboxInput(
            inputId = pg_id(species, "fit_view"),
            label = "Fit tree to view",
            value = TRUE
          ),
          sliderInput(
            inputId = pg_id(species, "tip_height"),
            label = "Tip spacing / vertical zoom (px):",
            min = 6, max = 30, value = 12, step = 2
          ),
          tags$small(class = "text-muted", "Tip spacing applies when 'Fit tree to view' is off."),
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
          uiOutput(pg_id(species, "tree_note")),
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

  # Parsed tree + joined tip metadata, shared by the plot and the LLM note.
  # Only invalidates when the underlying DuckDB data changes.
  tree_obj_r <- reactive({
    data <- tree_data_rv()
    if (!is.null(data$error)) return(list(error = data$error))
    if (is.null(data$trees)) return(list(error = "No tree data returned (unexpected NULL)."))
    if (!requireNamespace("ape", quietly = TRUE))
      return(list(error = "R package 'ape' is not installed."))

    trees_df <- data$trees
    tips_df <- data$tips

    # Prefer the Nextstrain/Augur tree (full metadata); fall back to nextclade
    non_iq <- trees_df[!grepl("iqtree", tolower(trees_df$tree_method)), , drop = FALSE]
    priority <- c("nextstrain", "augur", "nextclade")
    ranks <- match(tolower(non_iq$tree_method), priority, nomatch = 99)
    selected_tree <- if (nrow(non_iq) == 0) {
      trees_df[1, , drop = FALSE]
    } else {
      non_iq[order(ranks), ][1, , drop = FALSE]
    }

    newick <- selected_tree$newick
    if (is.na(newick) || nchar(trimws(newick)) == 0)
      return(list(error = "Tree newick is empty."))

    tmp <- tempfile(fileext = ".nwk")
    writeLines(newick, tmp)
    tr <- tryCatch(ape::read.tree(tmp), error = function(e) NULL)
    unlink(tmp)
    if (is.null(tr)) return(list(error = "Failed to parse tree newick."))

    # Normalize tip labels to match the database's _clean_newick_name() format
    tr$tip.label <- gsub("[,;():\\[\\] ]", "_", tr$tip.label)

    # Ladderize + truncate extreme branch lengths so one long branch
    # doesn't compress all other tips into a sliver
    tr <- ape::ladderize(tr)
    if (!is.null(tr$edge.length) && length(tr$edge.length) > 0) {
      q95 <- quantile(tr$edge.length[tr$edge.length > 0], 0.95, na.rm = TRUE)
      tr$edge.length[tr$edge.length > 3 * q95] <- 3 * q95
    }

    # Tip metadata for the selected tree; fall back to any tips of this species
    tip_meta <- tips_df[tips_df$tree_id == selected_tree$tree_id, ]
    if (nrow(tip_meta) == 0 && nrow(tips_df) > 0) {
      tip_meta <- tips_df[!duplicated(tips_df$label), ]
    }
    keep_cols <- intersect(
      c("label", "is_query", "clade", "outbreak", "country", "tip_date",
        "genome_coverage", "nuc_mutation_count", "aa_mutation_count"),
      names(tip_meta)
    )
    tip_meta <- tip_meta[, keep_cols, drop = FALSE]
    tip_meta$label <- gsub("[,;():\\[\\] ]", "_", tip_meta$label)
    tip_meta <- tip_meta[!duplicated(tip_meta$label), ]

    message(sprintf("[pg] Tip join: %d/%d tree labels matched (%s, tree_id=%s)",
                    sum(tr$tip.label %in% tip_meta$label), length(tr$tip.label),
                    sp, selected_tree$tree_id))

    tip_meta$clade[is.na(tip_meta$clade) | tip_meta$clade == ""] <- "Unknown"
    tip_meta$outbreak[is.na(tip_meta$outbreak) | tip_meta$outbreak == ""] <- "Unknown"
    tip_meta$country[is.na(tip_meta$country) | tip_meta$country == ""] <- "Unknown"
    tip_meta$is_query[is.na(tip_meta$is_query)] <- FALSE
    tip_meta$genome_coverage[is.na(tip_meta$genome_coverage)] <- 0
    tip_meta$nuc_mutation_count[is.na(tip_meta$nuc_mutation_count)] <- 0

    list(tr = tr, tree_id = selected_tree$tree_id,
         tree_method = selected_tree$tree_method, tip_meta = tip_meta)
  })

  # -- LLM interpretation note (local Ollama, deterministic template fallback) --
  evo_summary_r <- reactive({
    obj <- tree_obj_r()
    if (!is.null(obj$error) || is.null(obj$tr)) return(NULL)
    tryCatch(.pg_evo_summary(obj$tr, obj$tip_meta), error = function(e) {
      message("[pg] evo summary failed: ", e$message)
      NULL
    })
  })

  note_rv <- reactiveVal(NULL)
  note_busy <- reactiveVal(FALSE)

  # Pre-generated note written by the pipeline's GENERATE_TREE_NOTES step
  # (results/pathogen_genomics/<species>/tree_note.json). Read once — the
  # dashboard prefers it so no live Ollama call is needed at view time.
  pregen_note <- reactive({
    f <- file.path(outdir(), "pathogen_genomics", sp, "tree_note.json")
    if (!file.exists(f) || !requireNamespace("jsonlite", quietly = TRUE))
      return(NULL)
    tryCatch({
      j <- jsonlite::fromJSON(readLines(f, warn = FALSE), simplifyVector = FALSE)
      if (is.null(j$text) || !nzchar(j$text)) return(NULL)
      list(text = j$text, source = j$source %||% "template",
           model = j$model, error = j$error)
    }, error = function(e) NULL)
  })

  gen_note <- function(live_only = FALSE) {
    s <- evo_summary_r()
    if (is.null(s)) return()
    if (!live_only) {
      pg <- pregen_note()
      if (!is.null(pg)) { note_rv(pg); return() }
    }
    note_busy(TRUE)
    on.exit(note_busy(FALSE), add = TRUE)
    note_rv(tryCatch(.pg_tree_note(s), error = function(e)
      list(text = .pg_note_template(s), source = "template",
           model = NULL, error = e$message)))
  }

  # Auto-generate once per tree load (summary only invalidates on data change)
  observeEvent(evo_summary_r(), gen_note())
  # Regenerate always calls the LLM live, bypassing the pre-generated file
  observeEvent(input[[pg_id(sp, "regen_note")]], gen_note(live_only = TRUE),
               ignoreInit = TRUE)

  output[[pg_id(sp, "tree_note")]] <- renderUI({
    if (is.null(evo_summary_r())) return(NULL)
    n <- note_rv()
    regen <- actionLink(pg_id(sp, "regen_note"), "Regenerate",
                        style = "font-size:0.75rem;")
    if (is.null(n)) {
      return(div(style = "background:#f8f9fa;border-left:4px solid #4A6C8C;border-radius:4px;padding:8px 12px;font-size:0.85rem;margin-bottom:8px;",
        tags$em(class = "text-muted",
                if (isTRUE(note_busy())) "Generating interpretation\u2026" else "Interpretation pending\u2026"),
        " ", regen))
    }
    badge <- if (identical(n$source, "ollama")) {
      tags$span(class = "badge badge-info", style = "margin-right:6px;",
                paste0("AI \u00b7 ", n$model))
    } else {
      tags$span(class = "badge badge-secondary", style = "margin-right:6px;",
                "Auto-summary")
    }
    div(style = "background:#eef4f8;border-left:4px solid #4A6C8C;border-radius:4px;padding:10px 12px;font-size:0.85rem;margin-bottom:8px;",
      div(style = "margin-bottom:4px;", badge, regen),
      div(n$text),
      if (!is.null(n$error)) {
        tags$small(class = "text-muted", style = "display:block;margin-top:4px;",
                   paste0("LLM note: ", n$error))
      })
  })

  # Dynamic plot height: fit-to-view compresses the whole tree into the
  # visible panel; otherwise scale with tip count so annotations stay readable
  output[[pg_id(sp, "interactive_tree_ui")]] <- renderUI({
    obj <- tree_obj_r()
    n_tips <- if (!is.null(obj$tr)) length(obj$tr$tip.label) else 30
    tip_h <- input[[pg_id(sp, "tip_height")]]
    if (is.null(tip_h)) tip_h <- 12
    fit <- isTRUE(input[[pg_id(sp, "fit_view")]])
    plot_h <- if (fit) 700 else min(20000, max(500, n_tips * tip_h))
    plotOutput(pg_id(sp, "interactive_tree"), height = paste0(plot_h, "px"))
  })

  output[[pg_id(sp, "interactive_tree")]] <- renderPlot({
    tryCatch({
    obj <- tree_obj_r()
    if (!is.null(obj$error)) {
      plot.new()
      text(0.5, 0.5, paste0("No tree data available.\n\n", obj$error),
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

    tr <- obj$tr
    tip_meta <- obj$tip_meta

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
      if ("bootstrap" %in% annotations && grepl("iqtree", tolower(obj$tree_method))) {
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
      c("label", "is_query", "clade", "outbreak", "country", "tip_date", "div",
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
