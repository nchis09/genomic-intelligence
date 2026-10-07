# Influenza Segment Signatures module
#
# Per-species page for the SEGMENT_SIGNATURES pipeline outputs
# (results/segment_signatures/<meta.id>/*.tsv):
#   - Constellation & Reassortment: isolate x segment label matrix + flags
#   - Phenotype Markers: curated resistance/virulence catalogue hits
#   - Genome Integrity: defective-interfering screen (per-segment QC)
#   - Group Diversity: shared/private mutations per group_id x segment
#
# Every segment-bearing card gets its own segment dropdown — per the
# dashboard convention of per-card (not global) controls.

fs_id <- function(species, suffix) {
  paste0("fs_", species, "_", suffix)
}

fs_sig_dir <- function(outdir, species) {
  base <- file.path(outdir, "segment_signatures")
  if (!dir.exists(base)) return(NULL)
  hits <- list.dirs(base, recursive = FALSE, full.names = TRUE)
  hits <- hits[grepl(paste0("_", species, "$"), basename(hits), ignore.case = TRUE)]
  if (length(hits) == 0) NULL else hits[1]
}

fs_read <- function(outdir, species, filename) {
  dir <- fs_sig_dir(outdir, species)
  if (is.null(dir)) return(NULL)
  f <- file.path(dir, filename)
  if (!file.exists(f)) return(NULL)
  tryCatch(readr::read_tsv(f, show_col_types = FALSE), error = function(e) NULL)
}

SEGMENT_ORDER <- c("pb2", "pb1", "pa", "ha", "np", "na", "mp", "ns")

# Small inline segment dropdown for card headers ("drop menu in the corner")
fs_segment_select <- function(species, key, segments) {
  if (length(segments) <= 1 && !("genome" %in% segments)) return(NULL)
  div(
    style = "width: 130px; margin-left: auto;",
    selectizeInput(
      inputId = fs_id(species, paste0(key, "_segment")),
      label = NULL,
      choices = c("All" = "", setNames(segments, toupper(segments))),
      selected = "",
      multiple = FALSE,
      width = "100%"
    )
  )
}

# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------
flu_signatures_ui <- function(species) {
  tagList(
    h3("Influenza Segment Signatures"),
    p(class = "subtle-text", style = "margin-bottom: 16px;",
      "Per-segment genotype constellations, reassortment flags, genome-integrity (DI) screen, curated phenotype-marker hits and group-level diversity."),

    # -- Constellation & Reassortment --------------------------------------
    bs4Dash::bs4Card(
      title = "Genomic Constellation & Reassortment",
      width = 12,
      status = "primary",
      solidHeader = TRUE,
      collapsible = TRUE,
      uiOutput(fs_id(species, "constellation_note")),
      plotOutput(fs_id(species, "constellation_plot"), height = "auto"),
      hr(),
      h5("Reassortment flags", style = "font-weight: 600;"),
      DT::DTOutput(fs_id(species, "reassortment_table"))
    ),

    # -- Phenotype markers --------------------------------------------------
    bs4Dash::bs4Card(
      title = "Phenotype Markers (resistance / virulence / host adaptation)",
      width = 12,
      status = "warning",
      solidHeader = TRUE,
      collapsible = TRUE,
      div(
        style = "display: flex; gap: 16px; align-items: center; margin-bottom: 10px; flex-wrap: wrap;",
        div(style = "width: 220px;",
          selectizeInput(
            inputId = fs_id(species, "marker_category"),
            label = NULL,
            choices = NULL,
            multiple = FALSE,
            options = list(placeholder = "Category..."),
            width = "100%"
          )
        ),
        uiOutput(fs_id(species, "marker_segment_ui"))
      ),
      DT::DTOutput(fs_id(species, "markers_table"))
    ),

    # -- Genome integrity (DI screen) ---------------------------------------
    bs4Dash::bs4Card(
      title = "Genome Integrity — DI / truncation screen",
      width = 12,
      status = "info",
      solidHeader = TRUE,
      collapsible = TRUE,
      collapsed = FALSE,
      div(
        style = "display: flex; gap: 16px; align-items: center; margin-bottom: 10px;",
        p(class = "subtle-text", style = "font-size: 0.82rem; margin: 0; flex: 1;",
          "Consensus-level screen from each record's winning Nextclade QC row — true DI confirmation needs read-level coverage."),
        uiOutput(fs_id(species, "di_segment_ui"))
      ),
      DT::DTOutput(fs_id(species, "di_table"))
    ),

    # -- Group diversity ----------------------------------------------------
    bs4Dash::bs4Card(
      title = "Group Diversity (group_id × segment)",
      width = 12,
      status = "success",
      solidHeader = TRUE,
      collapsible = TRUE,
      div(
        style = "display: flex; gap: 16px; align-items: center; margin-bottom: 10px;",
        p(class = "subtle-text", style = "font-size: 0.82rem; margin: 0; flex: 1;",
          "Shared vs private mutations among isolates sharing a metadata group_id — intra-host / same-site diversity from consensus genomes."),
        uiOutput(fs_id(species, "gd_segment_ui"))
      ),
      DT::DTOutput(fs_id(species, "gd_table"))
    )
  )
}

# ---------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------
flu_signatures_register <- function(input, output, session, species, outdir) {
  sp <- species

  constellation_data <- reactive({
    fs_read(outdir(), sp, "constellation.tsv")
  })
  reassortment_data <- reactive({
    fs_read(outdir(), sp, "reassortment_flags.tsv")
  })
  markers_data <- reactive({
    fs_read(outdir(), sp, "markers.tsv")
  })
  di_data <- reactive({
    fs_read(outdir(), sp, "di_candidates.tsv")
  })
  gd_data <- reactive({
    fs_read(outdir(), sp, "group_diversity.tsv")
  })

  fs_seg <- function(key) {
    sel <- input[[fs_id(sp, paste0(key, "_segment"))]]
    if (is.null(sel) || sel == "") NULL else tolower(sel)
  }

  seg_values <- function(df) {
    if (is.null(df) || !"segment" %in% names(df)) return(character())
    segs <- unique(na.omit(df$segment))
    segs[order(match(segs, SEGMENT_ORDER, nomatch = 99))]
  }

  # -- Segment dropdowns (per-card corner controls) --
  output[[fs_id(sp, "marker_segment_ui")]] <- renderUI({
    fs_segment_select(sp, "marker", seg_values(markers_data()))
  })
  output[[fs_id(sp, "di_segment_ui")]] <- renderUI({
    fs_segment_select(sp, "di", seg_values(di_data()))
  })
  output[[fs_id(sp, "gd_segment_ui")]] <- renderUI({
    fs_segment_select(sp, "gd", seg_values(gd_data()))
  })

  # -- Marker category dropdown --
  observe({
    df <- markers_data()
    if (is.null(df) || nrow(df) == 0) return()
    cats <- sort(unique(na.omit(df$category)))
    updateSelectizeInput(session, fs_id(sp, "marker_category"),
                         choices = c("All categories" = "", cats),
                         selected = "", server = FALSE)
  })

  # -- Constellation summary note --
  output[[fs_id(sp, "constellation_note")]] <- renderUI({
    df <- constellation_data()
    if (is.null(df) || nrow(df) == 0) {
      return(div(
        class = "subtle-text",
        style = "padding: 24px; text-align: center;",
        "No constellation data — re-run with segment signatures enabled."))
    }
    iso <- unique(df$isolate)
    n_novel <- length(unique(df$isolate[df$is_novel %in% c(TRUE, "true", "True")]))
    flags <- reassortment_data()
    n_flag <- if (!is.null(flags)) sum(flags$reassortment_suspected %in% c(TRUE, "true", "True")) else 0
    div(
      style = "display: flex; gap: 16px; flex-wrap: wrap; margin-bottom: 12px;",
      div(class = "pi-summary-card",
        div(class = "pi-sc-label", "Isolates"),
        div(class = "pi-sc-value", length(iso))),
      div(class = "pi-summary-card",
        div(class = "pi-sc-label", "Reassortant"),
        div(class = "pi-sc-value", n_flag,
            style = if (n_flag > 0) "color: #E74C3C;" else NULL)),
      div(class = "pi-summary-card",
        div(class = "pi-sc-label", "Novel constellations"),
        div(class = "pi-sc-value", n_novel,
            style = if (n_novel > 0) "color: #E67E22;" else NULL))
    )
  })

  # -- Constellation tile matrix --
  constellation_plot_obj <- reactive({
    df <- constellation_data()
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df$segment <- factor(df$segment, levels = SEGMENT_ORDER)
    labels <- sort(unique(na.omit(df$label)))
    pal <- setNames(
      colorRampPalette(RColorBrewer::brewer.pal(8, "Set3"))(max(3, length(labels)))[seq_along(labels)],
      labels)
    df$label <- factor(df$label, levels = labels)
    novel <- unique(df$isolate[df$is_novel %in% c(TRUE, "true", "True")])

    p <- ggplot2::ggplot(df, ggplot2::aes(x = segment, y = isolate, fill = label)) +
      ggplot2::geom_tile(color = "white", linewidth = 0.4) +
      ggplot2::geom_text(ggplot2::aes(label = label), size = 3, color = "#333") +
      ggplot2::scale_fill_manual(values = pal, name = "Label", na.value = "#eeeeee") +
      ggplot2::labs(x = NULL, y = NULL,
                    title = "Genotype / lineage label per segment") +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 0, face = "bold"),
        legend.position = "none",
        panel.grid = ggplot2::element_blank()
      )
    if (length(novel) > 0) {
      p <- p + ggplot2::geom_tile(
        data = df[df$isolate %in% novel, ],
        ggplot2::aes(x = segment, y = isolate),
        fill = NA, color = "#E67E22", linewidth = 0.8
      )
    }
    p
  })

  output[[fs_id(sp, "constellation_plot")]] <- renderPlot({
    p <- constellation_plot_obj()
    if (is.null(p)) {
      plot.new()
      text(0.5, 0.5, "No constellation data.", cex = 1.1, col = "#6c757d")
      return()
    }
    print(p)
  }, height = function() {
    df <- constellation_data()
    n <- if (is.null(df)) 2 else max(2, length(unique(df$isolate)))
    max(160, n * 44 + 90)
  })

  # -- Reassortment table --
  output[[fs_id(sp, "reassortment_table")]] <- DT::renderDT({
    df <- reassortment_data()
    if (is.null(df) || nrow(df) == 0) {
      return(DT::datatable(data.frame(Message = "No reassortment flag data."),
                           rownames = FALSE, options = list(dom = "t")))
    }
    DT::datatable(
      df,
      rownames = FALSE,
      options = list(pageLength = 10, scrollX = TRUE, dom = "frtip"),
      class = "compact stripe"
    ) |>
      DT::formatStyle(
        "reassortment_suspected",
        target = "row",
        backgroundColor = DT::styleEqual(c(TRUE, "true", "True"), "#FDEDEC")
      )
  })

  # -- Phenotype markers table --
  output[[fs_id(sp, "markers_table")]] <- DT::renderDT({
    df <- markers_data()
    if (is.null(df) || nrow(df) == 0) {
      return(DT::datatable(
        data.frame(Message = "No curated marker hits detected in query isolates."),
        rownames = FALSE, options = list(dom = "t")))
    }
    cat <- input[[fs_id(sp, "marker_category")]]
    if (!is.null(cat) && cat != "") df <- df[df$category == cat, ]
    seg <- fs_seg("marker")
    if (!is.null(seg)) df <- df[df$segment == seg, ]

    keep <- intersect(c("sample", "segment", "mutation", "category", "phenotype",
                        "drug", "confidence", "coordinate_caveat", "source"),
                      names(df))
    d <- df[, keep, drop = FALSE]
    DT::datatable(
      d, rownames = FALSE,
      options = list(pageLength = 15, scrollX = TRUE, dom = "frtip"),
      class = "compact stripe",
      caption = "Marker hits are phenotype ASSOCIATIONS from the curated catalogue — not clinical calls. Rows with a coordinate caveat matched across numbering systems and need manual verification."
    ) |>
      DT::formatStyle(
        "coordinate_caveat",
        target = "row",
        backgroundColor = DT::styleEqual(c(1, "1", TRUE), "#FFF6E5")
      )
  })

  # -- Genome integrity table --
  output[[fs_id(sp, "di_table")]] <- DT::renderDT({
    df <- di_data()
    if (is.null(df) || nrow(df) == 0) {
      return(DT::datatable(data.frame(Message = "No genome-integrity data."),
                           rownames = FALSE, options = list(dom = "t")))
    }
    seg <- fs_seg("di")
    if (!is.null(seg)) df <- df[df$segment == seg, ]
    keep <- intersect(c("isolate", "record", "segment", "coverage",
                        "total_deletions", "total_frameshifts",
                        "total_stop_codons", "failed_cdses", "total_missing",
                        "di_candidate", "reasons"),
                      names(df))
    DT::datatable(
      df[, keep, drop = FALSE], rownames = FALSE,
      options = list(pageLength = 15, scrollX = TRUE, dom = "frtip"),
      class = "compact stripe"
    ) |>
      DT::formatStyle(
        "di_candidate",
        target = "row",
        backgroundColor = DT::styleEqual(c(TRUE, "true", "True"), "#FDEDEC")
      )
  })

  # -- Group diversity table --
  output[[fs_id(sp, "gd_table")]] <- DT::renderDT({
    df <- gd_data()
    if (is.null(df) || nrow(df) == 0) {
      return(DT::datatable(
        data.frame(Message = "No group_id metadata or only singleton groups."),
        rownames = FALSE, options = list(dom = "t")))
    }
    seg <- fs_seg("gd")
    if (!is.null(seg)) df <- df[df$segment == seg, ]
    DT::datatable(
      df, rownames = FALSE,
      options = list(pageLength = 15, scrollX = TRUE, dom = "frtip"),
      class = "compact stripe"
    )
  })
}
