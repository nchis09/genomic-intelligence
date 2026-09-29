# Countermeasure Readiness dashboard module
#
# Per-species view of countermeasure evidence: which diagnostics, vaccines,
# therapeutics, interventions and policies were reported and how each
# product/test performed, pooled with standard meta-analysis. No readiness
# status is assigned — the module reports products and performance only.

countermeasure_readiness_ui <- function(sp) {
  ns <- function(id) paste0(sp, "_", id)
  uiOutput(ns("cm_body"))
}

.cm_fmt_prop <- function(x) {
  ifelse(is.finite(x), paste0(formatC(x * 100, digits = 1, format = "f"), "%"), "")
}

# Policy summary: which policies were recommended, by whom, and their outcomes.
.cm_policy_summary <- function(ev) {
  if (is.null(ev) || !nrow(ev)) return(NULL)
  df <- ev[tolower(ev$domain) == "policy", , drop = FALSE]
  if (!nrow(df)) return(NULL)
  tbl <- data.frame(
    policy = .dv_get(df, "policy_name"),
    body = .dv_get(df, "issuing_body"),
    type = .dv_get(df, "policy_type"),
    recommendation = .dv_get(df, "recommendation"),
    outcome = .dv_get(df, "policy_outcome"),
    evidence = .dv_get(df, "evidence_basis"),
    travel = .dv_get(df, "travel_measure"),
    finding = .dv_get(df, "finding"),
    pmid = as.character(df$pmid),
    stringsAsFactors = FALSE
  )
  # Fall back to the finding text as the row identity when the structured
  # policy_name was not extracted (common in this domain).
  key <- ifelse(nzchar(tbl$policy), tbl$policy, substr(tbl$finding, 1, 80))
  tbl$key <- key
  tbl <- tbl[nzchar(key), , drop = FALSE]
  if (!nrow(tbl)) return(NULL)
  tbl |>
    dplyr::group_by(key) |>
    dplyr::summarise(
      policy = if (any(nzchar(policy))) dplyr::first(policy[nzchar(policy)]) else NA_character_,
      issuing_body = paste(unique(body[nzchar(body)]), collapse = "; "),
      type = paste(unique(type[nzchar(type)]), collapse = "; "),
      recommendation = paste(unique(recommendation[nzchar(recommendation)]), collapse = " | "),
      outcome = paste(unique(outcome[nzchar(outcome)]), collapse = " | "),
      travel_measure = paste(unique(travel[nzchar(travel)]), collapse = "; "),
      n_claims = dplyr::n(),
      pmids = paste(unique(pmid), collapse = ", "),
      .groups = "drop"
    ) |>
    dplyr::mutate(policy = ifelse(is.na(policy), "(unstructured claim)", policy))
}

countermeasure_content_ui <- function(sp, outdir, species) {
  ns <- function(id) paste0(sp, "_", id)

  bs4Dash::bs4Card(
    title = "Countermeasure Readiness",
    width = 12,
    status = "primary",
    fluidRow(
      column(3, selectInput(ns("cm_species"), "Species", choices = species, selected = sp, width = "100%")),
      column(3, selectInput(
        ns("cm_domain"), "Countermeasure",
        choices = setNames(names(DOM_COUNTERMEASURE), gsub("_", " ", names(DOM_COUNTERMEASURE))),
        width = "100%"
      ))
    ),
    uiOutput(ns("cm_inner"))
  )
}

countermeasure_readiness_register <- function(input, output, session, sp, outdir_r, species_rv) {
  ns <- function(id) paste0(sp, "_", id)

  sp_r <- reactive({
    sel <- input[[ns("cm_species")]]
    if (!is.null(sel) && nzchar(sel)) sel else sp
  })

  domain_r <- reactive({
    sel <- input[[ns("cm_domain")]]
    if (!is.null(sel) && nzchar(sel) && sel %in% names(DOM_COUNTERMEASURE)) sel else "vaccine_therapeutic"
  })

  data_r <- reactive({
    cur <- sp_r(); if (is.null(cur)) return(NULL)
    cm_gap_data(outdir_r(), cur)
  })

  ev_r <- reactive({
    d <- data_r(); if (is.null(d)) return(NULL); d$evidence
  })

  output[[ns("cm_body")]] <- renderUI({
    countermeasure_content_ui(sp, outdir_r(), species_rv())
  })

  output[[ns("cm_inner")]] <- renderUI({
    ev <- ev_r()
    if (is.null(ev) || !nrow(ev)) {
      return(div(
        style = "text-align: center; padding: 48px 0; color: #6c757d;",
        icon("syringe", class = "fa-3x"),
        h4("No countermeasure data", style = "margin-top: 16px;"),
        p("No literature evidence is available for this species.")
      ))
    }
    domain <- domain_r()
    spec <- DOM_COUNTERMEASURE[[domain]]
    narrative <- isTRUE(spec$narrative)
    has_metrics <- length(spec$metrics %||% list()) > 0
    tagList(
      uiOutput(ns("cm_coverage")),
      if (narrative) bs4Dash::bs4Card(
        title = "Narrative summary", width = 12, status = "info",
        collapsible = FALSE,
        uiOutput(ns("cm_narrative"))
      ),
      if (domain == "policy") {
        DT::DTOutput(ns("cm_policy_table"))
      } else {
        tagList(
          if (has_metrics) fluidRow(
            column(4, uiOutput(ns("cm_metric_ui"))),
            column(8, plotOutput(ns("cm_forest"), height = "420px"))
          ),
          DT::DTOutput(ns("cm_product_table"))
        )
      }
    )
  })

  # ---- Paper coverage strip: screened vs with-claims -------------------------

  output[[ns("cm_coverage")]] <- renderUI({
    d <- data_r(); if (is.null(d)) return(NULL)
    cov <- cm_paper_coverage(d$papers)
    if (is.null(cov)) return(NULL)
    row <- cov[cov$domain == domain_r(), , drop = FALSE]
    if (!nrow(row)) return(NULL)
    total <- colSums(cov[, c("screened", "with_claims")], na.rm = TRUE)
    div(
      style = "color:#555;margin:4px 0 10px 0;font-size:12px;",
      sprintf("%s: %d of %d screened papers yielded claims  |  overall: %d of %d",
              gsub("_", " ", domain_r()), row$with_claims, row$screened,
              total["with_claims"], total["screened"])
    )
  })

  # ---- Metric selector (within the chosen countermeasure domain) -------------

  domain_specs_r <- reactive({
    ev <- ev_r(); if (is.null(ev)) return(list())
    specs <- ma_domain_specs(domain_r())
    avail <- list()
    for (spec in specs) {
      ms <- ma_metric_studies(ev, domain_r(), spec$field, spec$type, spec$group)
      if (!is.null(ms) && nrow(ms$studies) >= 1) avail[[spec$field]] <- spec
    }
    avail
  })

  output[[ns("cm_metric_ui")]] <- renderUI({
    specs <- domain_specs_r()
    if (!length(specs)) return(div(style = "padding:24px;color:#888;", "No performance metrics extracted for this domain/species."))
    choices <- setNames(names(specs), vapply(names(specs), ma_metric_label, character(1)))
    selectInput(ns("cm_metric"), "Metric", choices = choices, width = "100%")
  })

  metric_r <- reactive({
    specs <- domain_specs_r()
    if (!length(specs)) return(NULL)
    sel <- input[[ns("cm_metric")]]
    if (is.null(sel) || !(sel %in% names(specs))) sel <- names(specs)[1]
    specs[[sel]]
  })

  output[[ns("cm_forest")]] <- renderPlot({
    ev <- ev_r(); if (is.null(ev)) return(NULL)
    spec <- metric_r(); if (is.null(spec)) return(NULL)
    ms <- ma_metric_studies(ev, domain_r(), spec$field, spec$type, spec$group)
    pr <- ma_pool_metric(ms, ev)
    ma_forest_gg(ms, pr$pooled,
                 title = paste0(gsub("_", " ", domain_r()), ": ", ma_metric_label(spec$field)))
  })

  # ---- Narrative (text-heavy domains) ----------------------------------------

  output[[ns("cm_narrative")]] <- renderUI({
    ev <- ev_r(); if (is.null(ev) || !nrow(ev)) return(NULL)
    note <- .dv_note(domain_r(), ev, sp_r(), data_r()$summaries)
    tagList(
      p(htmltools::htmlEscape(note$text)),
      p(style = "color:#888;font-size:11px;",
        paste0("Summary generated at pipeline time (", note$source,
               if (nzchar(note$model %||% "")) paste0(" · ", note$model) else "", ")."))
    )
  })

  # ---- Per-product pooled performance table (wide: one row per product) ------

  output[[ns("cm_product_table")]] <- DT::renderDT({
    ev <- ev_r(); if (is.null(ev) || !nrow(ev)) return(NULL)
    domain <- domain_r()
    spec <- DOM_COUNTERMEASURE[[domain]]
    metrics <- spec$metrics
    if (is.null(metrics) || !length(metrics)) return(NULL)

    tbl <- .dv_group_metric_table(ev, domain, spec$group, metrics)
    if (is.null(tbl) || !nrow(tbl)) {
      return(DT::datatable(data.frame(Note = "No extractable performance values for this domain."),
                           rownames = FALSE, options = list(dom = "t")))
    }

    # Append domain context columns (what the product was tested on/where).
    df <- ev[tolower(ev$domain) == domain, , drop = FALSE]
    df$.grp <- .dv_group_col(df, spec$group)
    for (col in (spec$context_cols %||% character(0))) {
      if (!col %in% names(df)) next
      vals <- vapply(tbl$Group, function(g)
        .dv_context(df[df$.grp == g, , drop = FALSE], col), character(1))
      if (any(nzchar(vals))) tbl[[ma_metric_label(col)]] <- vals
    }

    names(tbl)[1] <- switch(domain,
      diagnostic = "Test",
      vaccine_therapeutic = "Product",
      intervention = "Intervention",
      "Item")
    DT::datatable(tbl, rownames = FALSE, selection = "none",
                  options = list(pageLength = 15, dom = "ftip", scrollX = TRUE))
  })

  # ---- Policy ----------------------------------------------------------------

  output[[ns("cm_policy_table")]] <- DT::renderDT({
    ev <- ev_r(); if (is.null(ev)) return(NULL)
    tbl <- .cm_policy_summary(ev)
    if (is.null(tbl) || !nrow(tbl)) {
      return(DT::datatable(
        data.frame(Note = "No policy evidence extracted for this species."),
        rownames = FALSE, options = list(dom = "t")
      ))
    }
    DT::datatable(tbl, rownames = FALSE, selection = "none",
                  options = list(pageLength = 10, dom = "ftip", scrollX = TRUE))
  })
}
