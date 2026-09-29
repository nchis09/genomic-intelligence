# Evidence & Knowledge Gaps dashboard module
#
# Per-species view of what the literature says about the pathogen itself:
# clinical severity, epidemiology, transmission, seroprevalence, reservoirs,
# genomic surveillance and variant phenotypes — with pooled estimates,
# forest plots, field completeness and knowledge-gap cards. Countermeasure
# domains (diagnostic / vaccine / intervention / policy) live in the
# Countermeasure Readiness module.

evidence_knowledge_gaps_ui <- function(sp) {
  ns <- function(id) paste0(sp, "_", id)
  uiOutput(ns("ekg_body"))
}

.ekg_format_details <- function(details) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) return("")
  if (is.na(details) || details == "" || details == "{}") return("")
  tryCatch({
    d <- jsonlite::fromJSON(details, simplifyVector = TRUE)
    if (length(d) == 0) return("")
    d <- d[!vapply(d, function(x) is.null(x) || (length(x) == 1 && is.na(x)) || (length(x) == 1 && x == "") || (length(x) == 1 && x == "null"), logical(1))]
    paste(names(d), vapply(d, function(x) paste(as.character(x), collapse = ", "), character(1)), sep = ": ", collapse = "; ")
  }, error = function(e) as.character(details))
}

.ekg_data <- function(outdir, sp) {
  cm_gap_data(outdir, sp)
}

.ekg_priority_color <- function(p) {
  switch(as.character(p),
         "high" = "danger",
         "medium" = "warning",
         "low" = "info",
         "secondary")
}

.ekg_fmt_est <- function(x, type) {
  if (!is.finite(x)) return("")
  if (type == "proportion") return(paste0(formatC(x * 100, digits = 1, format = "f"), "%"))
  formatC(x, digits = 2, format = "f")
}

# Pathogen-domain gaps: missing domains and fields never populated.
.ekg_domain_gaps <- function(ev, species) {
  expected <- MA_EVIDENCE_DOMAINS
  present <- unique(tolower(ev$domain))
  gaps <- list()

  for (d in expected) {
    if (!(d %in% present)) {
      gaps[[length(gaps) + 1]] <- list(
        topic = gsub("_", " ", d),
        gap = sprintf("No literature evidence was extracted for the '%s' domain in %s.", gsub("_", " ", d), toupper(species)),
        priority = "high",
        n_papers_touching_topic = 0L,
        notes = "Either no papers were retrieved or none yielded extractable claims."
      )
    }
  }

  comp <- ma_field_completeness(ev, MA_EVIDENCE_DOMAINS)
  if (!is.null(comp)) {
    comp <- comp[comp$completeness == 0 & comp$claims > 0, ]
    for (i in seq_len(nrow(comp))) {
      gaps[[length(gaps) + 1]] <- list(
        topic = paste0(gsub("_", " ", comp$domain[i]), ": ", gsub("_", " ", comp$field[i])),
        gap = sprintf("Field '%s' was never populated across %d '%s' claims.",
                      gsub("_", " ", comp$field[i]), comp$claims[i], gsub("_", " ", comp$domain[i])),
        priority = "medium",
        n_papers_touching_topic = comp$claims[i],
        notes = "Papers mention the domain but not this metric, or extraction could not locate it."
      )
    }
  }
  gaps
}

evidence_knowledge_gaps_content <- function(sp, outdir, species) {
  ns <- function(id) paste0(sp, "_", id)

  bs4Dash::bs4Card(
    title = "Evidence & Knowledge Gaps",
    width = 12,
    status = "primary",
    fluidRow(
      column(3, selectInput(ns("ekg_species"), "Species", choices = species, selected = sp, width = "100%"))
    ),
    uiOutput(ns("ekg_inner"))
  )
}

ekg_inner_ui <- function(sp) {
  ns <- function(id) paste0(sp, "_", id)
  bs4Dash::bs4TabCard(
    title = NULL,
    width = 12,
    side = "right",
    bs4Dash::bs4TabItem(
      tabName = ns("tab_meta"),
      active = TRUE,
      fluidRow(
        column(3, selectInput(ns("ekg_domain"), "Domain",
                              choices = c("All domains" = "all",
                                          setNames(MA_EVIDENCE_DOMAINS,
                                                   gsub("_", " ", MA_EVIDENCE_DOMAINS))),
                              width = "100%"))
      ),
      uiOutput(ns("ekg_domain_ui"))
    ),
    bs4Dash::bs4TabItem(
      tabName = ns("tab_gaps"),
      fluidRow(
        column(7, plotOutput(ns("coverage_plot"), height = "420px"),
                  plotOutput(ns("paper_cov_plot"), height = "220px")),
        column(5, uiOutput(ns("ekg_overview")), uiOutput(ns("gaps_ui")))
      )
    ),
    bs4Dash::bs4TabItem(
      tabName = ns("tab_evidence"),
      DT::DTOutput(ns("evidence_table"))
    )
  )
}

evidence_knowledge_gaps_register <- function(input, output, session, sp, outdir_r, species_rv) {
  ns <- function(id) paste0(sp, "_", id)

  sp_r <- reactive({
    sel <- input[[ns("ekg_species")]]
    if (!is.null(sel) && nzchar(sel)) sel else sp
  })

  data_r <- reactive({
    cur <- sp_r(); if (is.null(cur)) return(NULL)
    .ekg_data(outdir_r(), cur)
  })

  ev_r <- reactive({
    d <- data_r(); if (is.null(d)) return(NULL); d$evidence
  })

  gaps_r <- reactive({
    d <- data_r(); if (is.null(d)) return(NULL); d$gaps
  })

  output[[ns("ekg_body")]] <- renderUI({
    evidence_knowledge_gaps_content(sp, outdir_r(), species_rv())
  })

  output[[ns("ekg_inner")]] <- renderUI({
    ev <- ev_r()
    if (is.null(ev) || !nrow(ev)) {
      return(div(
        style = "text-align: center; padding: 48px 0; color: #6c757d;",
        icon("magnifying-glass", class = "fa-3x"),
        h4("No evidence data", style = "margin-top: 16px;"),
        p("No extracted evidence is available for this species.")
      ))
    }
    ekg_inner_ui(sp)
  })

  # ---- Domain-specific evidence view -----------------------------------------

  ekg_domain_r <- reactive({
    sel <- input[[ns("ekg_domain")]]
    if (!is.null(sel) && nzchar(sel)) sel else "all"
  })

  output[[ns("ekg_domain_ui")]] <- renderUI({
    domain <- ekg_domain_r()
    if (domain == "all") {
      return(tagList(
        DT::DTOutput(ns("ma_summary")),
        br(),
        fluidRow(
          column(4, uiOutput(ns("ma_sel_ui"))),
          column(8, plotOutput(ns("ma_plot"), height = "380px"))
        )
      ))
    }
    spec <- DOM_EVIDENCE[[domain]]
    tagList(
      uiOutput(ns("ekg_coverage")),
      if (isTRUE(spec$narrative)) bs4Dash::bs4Card(
        title = "Narrative summary", width = 12, status = "info",
        collapsible = FALSE,
        uiOutput(ns("ekg_dom_note"))
      ),
      fluidRow(
        column(6, DT::DTOutput(ns("ekg_dom_table"))),
        column(6, plotOutput(ns("ekg_dom_plot"), height = "400px"))
      )
    )
  })

  output[[ns("ekg_coverage")]] <- renderUI({
    d <- data_r(); if (is.null(d)) return(NULL)
    cov <- cm_paper_coverage(d$papers)
    if (is.null(cov)) return(NULL)
    row <- cov[cov$domain == ekg_domain_r(), , drop = FALSE]
    if (!nrow(row)) return(NULL)
    div(
      style = "color:#555;margin:4px 0 10px 0;font-size:12px;",
      sprintf("%s: %d of %d screened papers yielded claims",
              gsub("_", " ", ekg_domain_r()), row$with_claims, row$screened)
    )
  })

  output[[ns("ekg_dom_table")]] <- DT::renderDT({
    ev <- ev_r(); if (is.null(ev) || !nrow(ev)) return(NULL)
    domain <- ekg_domain_r()
    if (domain == "all") return(NULL)
    tbl <- .dv_domain_table(ev, domain)
    if (is.null(tbl) || !nrow(tbl)) {
      return(DT::datatable(data.frame(Note = "No structured values extracted for this domain."),
                           rownames = FALSE, options = list(dom = "t")))
    }
    DT::datatable(tbl, rownames = FALSE, selection = "none",
                  options = list(pageLength = 12, dom = "ftip", scrollX = TRUE))
  })

  output[[ns("ekg_dom_plot")]] <- renderPlot({
    ev <- ev_r(); if (is.null(ev) || !nrow(ev)) return(NULL)
    domain <- ekg_domain_r()
    if (domain == "all") return(NULL)
    .dv_domain_plot(ev, domain)
  })

  output[[ns("ekg_dom_note")]] <- renderUI({
    ev <- ev_r(); if (is.null(ev) || !nrow(ev)) return(NULL)
    note <- .dv_note(ekg_domain_r(), ev, sp_r(), data_r()$summaries)
    tagList(
      p(htmltools::htmlEscape(note$text)),
      p(style = "color:#888;font-size:11px;",
        paste0("Summary generated at pipeline time (", note$source,
               if (nzchar(note$model %||% "")) paste0(" · ", note$model) else "", ")."))
    )
  })

  # ---- Pooled estimates ----------------------------------------------------

  output[[ns("ma_summary")]] <- DT::renderDT({
    ev <- ev_r(); if (is.null(ev) || !nrow(ev)) return(NULL)
    s <- ma_all_metric_summary(ev, MA_EVIDENCE_DOMAINS)
    if (is.null(s) || !nrow(s)) return(NULL)
    disp <- data.frame(
      Domain = s$domain,
      Metric = s$metric,
      `Studies (k)` = s$k,
      `Pooled` = mapply(function(v, t) .ekg_fmt_est(v, t), s$pooled, s$type),
      `CI / IQR` = mapply(function(lo, hi, t) {
        if (!is.finite(lo) || !is.finite(hi)) return("")
        paste0(.ekg_fmt_est(lo, t), " - ", .ekg_fmt_est(hi, t))
      }, s$lo, s$hi, s$type),
      `I2` = ifelse(is.finite(s$i2), paste0(round(s$i2 * 100), "%"), "-"),
      Method = ifelse(s$method == "meta", "random-effects", s$method),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    DT::datatable(
      disp, rownames = FALSE, selection = "none",
      options = list(pageLength = 10, dom = "ftip", scrollX = TRUE)
    )
  })

  metric_specs_r <- reactive({
    ev <- ev_r(); if (is.null(ev) || !nrow(ev)) return(list())
    ma_available_specs(ev, MA_EVIDENCE_DOMAINS)
  })

  output[[ns("ma_sel_ui")]] <- renderUI({
    specs <- metric_specs_r()
    if (!length(specs)) return(div(style = "padding:24px;color:#888;", "No metrics with extractable values."))
    choices <- setNames(
      names(specs),
      vapply(names(specs), function(k) paste0(gsub("_", " ", strsplit(k, "|", fixed = TRUE)[[1]][1]), " | ",
                                             ma_metric_label(strsplit(k, "|", fixed = TRUE)[[1]][2])), character(1))
    )
    selectInput(ns("ma_metric"), "Metric", choices = choices, width = "100%")
  })

  output[[ns("ma_plot")]] <- renderPlot({
    specs <- metric_specs_r()
    sel <- input[[ns("ma_metric")]]
    if (!length(specs)) return(NULL)
    if (is.null(sel) || !(sel %in% names(specs))) sel <- names(specs)[1]
    parts <- strsplit(sel, "|", fixed = TRUE)[[1]]
    spec <- specs[[sel]]
    ev <- ev_r()
    ms <- ma_metric_studies(ev, parts[1], spec$field, spec$type, spec$group)
    pr <- ma_pool_metric(ms, ev)
    ma_forest_gg(ms, pr$pooled, title = paste0(parts[1], ": ", ma_metric_label(spec$field)))
  })

  # ---- Coverage & gaps -----------------------------------------------------

  output[[ns("coverage_plot")]] <- renderPlot({
    ev <- ev_r(); if (is.null(ev) || !nrow(ev)) return(NULL)
    comp <- ma_field_completeness(ev, MA_EVIDENCE_DOMAINS)
    if (is.null(comp) || !nrow(comp)) return(NULL)
    comp$domain <- gsub("_", " ", comp$domain)
    comp$field <- gsub("_", " ", comp$field)
    ggplot(comp, aes(x = field, y = domain, fill = completeness)) +
      geom_tile(color = "white") +
      geom_text(aes(label = ifelse(completeness > 0, paste0(round(completeness * 100), "%"), "0")),
                size = 3, colour = "#333") +
      scale_fill_gradient(low = "#f7fbff", high = "#08519c", limits = c(0, 1),
                          labels = scales::percent_format()) +
      labs(x = NULL, y = NULL, fill = "% claims\nwith value",
           title = "Field completeness by domain") +
      theme_minimal(base_size = 11) +
      theme(axis.text.x = element_text(angle = 40, hjust = 1),
            panel.grid = element_blank())
  })

  output[[ns("paper_cov_plot")]] <- renderPlot({
    d <- data_r(); if (is.null(d)) return(NULL)
    cov <- cm_paper_coverage(d$papers)
    if (is.null(cov) || !nrow(cov)) return(NULL)
    dd <- data.frame(
      domain = rep(cov$domain_label, 2),
      n = c(cov$screened - cov$with_claims, cov$with_claims),
      kind = rep(c("screened, no claims", "yielded claims"), each = nrow(cov))
    )
    dd$kind <- factor(dd$kind, levels = c("yielded claims", "screened, no claims"))
    ggplot(dd, aes(x = n, y = stats::reorder(domain, n), fill = kind)) +
      geom_col() +
      scale_fill_manual(values = c("yielded claims" = "#2c7fb8",
                                   "screened, no claims" = "#d9d9d9")) +
      labs(x = "papers", y = NULL, fill = NULL,
           title = "Papers screened vs yielding claims") +
      theme_minimal(base_size = 11) +
      theme(legend.position = "top")
  })

  output[[ns("ekg_overview")]] <- renderUI({
    d <- data_r(); if (is.null(d)) return(NULL)
    ov <- .dv_overview(sp_r(), d$summaries)
    if (is.null(ov)) {
      return(bs4Dash::bs4Card(
        title = "Species summary", width = 12, status = "info", collapsible = FALSE,
        p("No pre-generated summary — rerun the pipeline so EVIDENCE_SUMMARY writes domain_summaries.tsv.")
      ))
    }
    bs4Dash::bs4Card(
      title = paste0(toupper(sp_r()), " — what we know / what we don't"),
      width = 12, status = "info", collapsible = FALSE,
      p(htmltools::htmlEscape(ov$text)),
      p(style = "color:#888;font-size:11px;",
        paste0("Generated at pipeline time (", ov$source,
               if (nzchar(ov$model %||% "")) paste0(" · ", ov$model) else "", ")."))
    )
  })

  output[[ns("gaps_ui")]] <- renderUI({
    ev <- ev_r()
    if (is.null(ev)) return(div(style = "color:#888;padding:24px;", "No evidence loaded."))
    cm_gaps <- gaps_r()
    domain_gaps <- .ekg_domain_gaps(ev, sp_r())

    cards <- list()
    if (!is.null(cm_gaps) && nrow(cm_gaps)) {
      for (i in seq_len(nrow(cm_gaps))) {
        cards[[length(cards) + 1]] <- list(
          topic = cm_gaps$topic[i], gap = cm_gaps$gap[i],
          priority = cm_gaps$priority[i], n = cm_gaps$n_papers_touching_topic[i],
          notes = cm_gaps$notes[i]
        )
      }
    }
    cards <- c(cards, domain_gaps)
    if (!length(cards)) {
      return(div(style = "color:#888;padding:24px;", "No knowledge gaps identified."))
    }
    ord <- order(match(vapply(cards, function(g) g$priority %||% "low", character(1)),
                       c("high", "medium", "low")))
    items <- lapply(cards[ord], function(g) {
      color <- .ekg_priority_color(g$priority)
      bs4Dash::bs4Card(
        title = g$topic,
        status = color,
        solidHeader = TRUE,
        width = 12,
        p(strong("Gap: "), g$gap),
        p(strong("Priority: "), tolower(g$priority)),
        p(strong("Related claims: "), as.character(g$n %||% g$n_papers_touching_topic %||% 0)),
        if (!is.null(g$notes) && !is.na(g$notes) && g$notes != "") p(strong("Notes: "), g$notes) else NULL
      )
    })
    do.call(tagList, items)
  })

  # ---- Underlying claims (drill-down) ---------------------------------------

  output[[ns("evidence_table")]] <- DT::renderDT({
    df <- ev_r(); if (is.null(df) || !nrow(df)) return(NULL)
    df <- df[tolower(df$domain) %in% MA_EVIDENCE_DOMAINS, , drop = FALSE]
    if (!nrow(df)) return(NULL)
    display <- df |>
      dplyr::mutate(details_summary = vapply(details, .ekg_format_details, character(1))) |>
      dplyr::select(pmid, domain, topic, claim_type,
                    finding, evidence_level, geography, details_summary)
    DT::datatable(
      display,
      options = list(pageLength = 10, dom = 'ftip', scrollX = TRUE),
      rownames = FALSE,
      selection = 'none'
    )
  })
}
