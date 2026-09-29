# Domain view registry
#
# Each literature domain gets its own display contract instead of one generic
# table + forest: which columns the summary table has, which plot to draw and
# whether a narrative (Ollama) note is produced. Everything is per-species —
# callers pass an already species-filtered evidence data frame (claim-level,
# domain fields already widened into columns).

if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a)) b else a

.dv_get <- function(df, col) {
  if (!col %in% names(df)) return(rep("", nrow(df)))
  v <- as.character(df[[col]])
  v[is.na(v)] <- ""
  trimws(v)
}

.dv_nonempty <- function(df, col) {
  v <- .dv_get(df, col)
  v[!v %in% c("", "null", "NA")]
}

# Frequency table of a categorical field (bar chart for textual domains).
.dv_freq_bar <- function(df, col, title) {
  if (is.null(df) || !nrow(df)) return(NULL)
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(NULL)
  vals <- .dv_multi_values(df, col)
  if (!length(vals)) return(NULL)
  tab <- sort(table(vals), decreasing = TRUE)
  tab <- head(tab, 15)
  dd <- data.frame(value = names(tab), n = as.integer(tab), stringsAsFactors = FALSE)
  ggplot2::ggplot(dd, ggplot2::aes(x = stats::reorder(value, n), y = n)) +
    ggplot2::geom_col(fill = "#2c7fb8") +
    ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = "claims", title = title) +
    ggplot2::theme_minimal(base_size = 11)
}

# Split a possibly-JSON list field ("[\"fever\",\"headache\"]" or "fever; headache")
# into atomic values across rows.
.dv_multi_values <- function(df, col) {
  vals <- .dv_nonempty(df, col)
  out <- character(0)
  for (v in vals) {
    if (grepl("^\\[", v)) {
      parsed <- tryCatch(jsonlite::fromJSON(v), error = function(e) NULL)
      if (!is.null(parsed)) {
        out <- c(out, trimws(as.character(parsed)))
        next
      }
    }
    out <- c(out, trimws(unlist(strsplit(v, ";|,", perl = TRUE))))
  }
  out[nzchar(out)]
}

.dv_fmt_ci <- function(est, lo, hi, prop = TRUE) {
  if (!is.finite(est)) return("")
  f <- function(x) if (prop) paste0(formatC(x * 100, digits = 1, format = "f"), "%") else formatC(x, digits = 2, format = "f")
  if (is.finite(lo) && is.finite(hi)) paste0(f(est), " (", f(lo), " - ", f(hi), ")") else f(est)
}

# Resolved group vector: primary field with the fallback chain applied.
.dv_group_col <- function(df, field) {
  grp <- .dv_get(df, field)
  for (fb in (.ma_group_fallbacks[[field]] %||% character(0))) {
    need <- !nzchar(grp)
    if (!any(need)) break
    alt <- .dv_get(df, fb)
    grp[need & nzchar(alt)] <- alt[need & nzchar(alt)]
  }
  grp[!nzchar(grp)] <- "unspecified"
  grp
}

.dv_group_metric_table <- function(ev, domain, group_field, metrics) {
  df <- ev[tolower(ev$domain) == domain, , drop = FALSE]
  if (!nrow(df)) return(NULL)
  grp <- .dv_group_col(df, group_field)
  df$.grp <- grp
  # One row per product/test — even when no metric parses, so safety/phase
  # context still shows instead of the product silently disappearing.
  out <- data.frame(group = unique(df$.grp), stringsAsFactors = FALSE)
  for (field in names(metrics)) {
    gp <- ma_group_pooled(df, domain, field, ".grp", metrics[[field]])
    if (is.null(gp) || !nrow(gp)) next
    prop <- metrics[[field]] == "proportion"
    cell <- data.frame(
      group = gp$group,
      v = mapply(.dv_fmt_ci, gp$est, gp$lo, gp$hi,
                 MoreArgs = list(prop = prop)),
      stringsAsFactors = FALSE
    )
    names(cell)[2] <- ma_metric_label(field)
    out <- merge(out, cell, by = "group", all.x = TRUE)
  }
  if (!nrow(out)) return(NULL)
  out[is.na(out)] <- ""
  # studies per group
  k <- as.integer(table(df$.grp)[out$group])
  data.frame(Group = out$group, `Studies (k)` = k, out[, -1, drop = FALSE],
             check.names = FALSE, stringsAsFactors = FALSE)
}

# Unique joined context values for a group of rows (e.g. study phases seen).
.dv_context <- function(df, col, sep = "; ") {
  v <- unique(.dv_nonempty(df, col))
  paste(head(v, 4), collapse = sep)
}

# ---------------------------------------------------------------------------
# LLM narrative for text-heavy domains (reuses .ollama_generate from llm_note.R)
# ---------------------------------------------------------------------------

.dv_note_cache <- new.env(parent = emptyenv())

# Deterministic fallback when Ollama is offline.
.dv_note_template <- function(domain, facts) {
  dom <- gsub("_", " ", domain)
  top <- function(x) {
    x <- unlist(x)
    if (!length(x)) return(NULL)
    paste(names(head(sort(x, decreasing = TRUE), 3)), collapse = "; ")
  }
  parts <- c(sprintf("The %s literature for %s yields %d claims across %d papers.",
                     dom, toupper(facts$species %||% "this species"),
                     facts$n_claims %||% 0, facts$n_papers %||% 0))
  for (k in names(facts$frequencies %||% list())) {
    t <- top(facts$frequencies[[k]])
    if (!is.null(t) && nzchar(t))
      parts <- c(parts, sprintf("Most-mentioned %s: %s.", gsub("_", " ", k), t))
  }
  paste(parts, collapse = " ")
}

# Build a compact facts object for the narrative prompt.
.dv_facts <- function(domain, ev, species) {
  df <- ev[tolower(ev$domain) == domain, , drop = FALSE]
  facts <- list(
    species = species, domain = domain,
    n_claims = nrow(df),
    n_papers = length(unique(df$pmid)),
    frequencies = list()
  )
  freq_fields <- switch(domain,
    policy = c("issuing_body", "policy_type", "intervention", "travel_measure", "policy_outcome"),
    intervention = c("intervention_name", "intervention_type", "implementation_setting", "target_population"),
    reservoir = c("host_species", "detection_method", "geographic_scope"),
    clinical = c("clinical_severity", "treatment_setting"),
    transmission = c("transmission_route", "transmission_setting", "risk_factors"),
    genomic_surveillance = c("sequencing_platform", "lineage", "clade", "genotype"),
    variant_phenotype = c("gene", "phenotype"),
    NULL)
  for (f in freq_fields %||% character(0)) {
    vals <- .dv_nonempty(df, f)
    if (length(vals)) facts$frequencies[[f]] <- as.list(table(vals))
  }
  # Sampled findings give the LLM the substance (cap to keep prompts small).
  finds <- .dv_nonempty(df, "finding")
  facts$findings <- head(unique(finds), 40)
  facts
}

# Narrative for a domain: pre-generated pipeline summary first (the dashboard
# never calls the LLM itself), template fallback when no summary exists.
.dv_note <- function(domain, ev, species, summaries = NULL) {
  if (!is.null(summaries) && nrow(summaries)) {
    row <- summaries[summaries$domain == domain, , drop = FALSE]
    if (nrow(row) && !is.na(row$summary[1]) && nzchar(row$summary[1]))
      return(list(text = row$summary[1],
                  source = row$source[1] %||% "pipeline",
                  model = row$model[1] %||% ""))
  }
  facts <- .dv_facts(domain, ev, species)
  list(text = .dv_note_template(domain, facts), source = "template", model = "")
}

.dv_overview <- function(species, summaries) {
  if (!is.null(summaries) && nrow(summaries)) {
    row <- summaries[summaries$domain == "species_overview", , drop = FALSE]
    if (nrow(row) && nzchar(row$summary[1] %||% ""))
      return(list(text = row$summary[1],
                  source = row$source[1] %||% "pipeline",
                  model = row$model[1] %||% ""))
  }
  NULL
}

# ---------------------------------------------------------------------------
# Domain registry — what each domain shows
# ---------------------------------------------------------------------------

DOM_COUNTERMEASURE <- list(
  diagnostic = list(
    group = "test_name",
    metrics = c(sensitivity = "proportion", specificity = "proportion",
                positive_predictive_value = "proportion",
                negative_predictive_value = "proportion",
                limit_of_detection = "value"),
    context_cols = c("test_type", "sample_type", "target_gene_or_antigen", "study_design")
  ),
  vaccine_therapeutic = list(
    group = "product_name",
    metrics = c(efficacy = "proportion", effectiveness = "proportion"),
    context_cols = c("product_type", "study_phase", "study_design", "endpoint",
                     "neutralization_result", "dose", "treatment_regimen", "adverse_events")
  ),
  intervention = list(
    group = "intervention_name",
    metrics = c(effectiveness = "proportion", coverage = "proportion"),
    context_cols = c("intervention_type", "implementation_setting", "target_population",
                     "outcome", "barriers", "facilitators"),
    narrative = TRUE
  ),
  policy = list(
    group = "policy_name",
    metrics = NULL,
    narrative = TRUE
  )
)

DOM_EVIDENCE <- list(
  clinical = list(
    metrics = c(case_fatality_rate = "proportion", hospitalization_rate = "proportion",
                incubation_period = "value"),
    freq_field = "symptoms"
  ),
  epidemiological_context = list(
    metrics = c(attack_rate = "proportion", incidence = "proportion", prevalence = "proportion"),
    table_fn = "outbreaks"
  ),
  transmission = list(
    metrics = c(r0 = "value", serial_interval = "value", generation_time = "value",
                secondary_attack_rate = "proportion"),
    freq_field = "transmission_route"
  ),
  seroprevalence = list(
    metrics = c(seroprevalence_value = "proportion"),
    group = "host_species"
  ),
  reservoir = list(
    metrics = c(prevalence = "proportion"),
    group = "host_species",
    narrative = TRUE
  ),
  genomic_surveillance = list(
    metrics = c(sequence_count = "total"),
    freq_field = "lineage",
    narrative = TRUE
  ),
  variant_phenotype = list(
    table_fn = "variants"
  )
)

# ---------------------------------------------------------------------------
# Per-domain summary tables and plots (Evidence & Knowledge Gaps module)
# ---------------------------------------------------------------------------

# Metric summary table: one row per populated numeric metric.
.dv_metric_summary_table <- function(ev, domain, metrics) {
  rows <- list()
  for (field in names(metrics)) {
    ms <- ma_metric_studies(ev, domain, field, metrics[[field]])
    if (is.null(ms)) next
    pr <- ma_pool_metric(ms, ev)
    p <- pr$pooled
    if (is.null(p)) next
    prop <- metrics[[field]] == "proportion"
    rows[[length(rows) + 1]] <- data.frame(
      Metric = ma_metric_label(field),
      `Studies (k)` = p$k,
      Pooled = .dv_fmt_ci(p$estimate, p$ci_lo, p$ci_hi, prop),
      `I2` = if (!is.null(p$i2) && is.finite(p$i2)) paste0(round(p$i2 * 100), "%") else "-",
      check.names = FALSE, stringsAsFactors = FALSE
    )
  }
  if (!length(rows)) return(NULL)
  do.call(rbind, rows)
}

# Per-outbreak table for epidemiological_context — totals are shown per
# outbreak, never pooled across them.
.dv_outbreak_table <- function(ev) {
  df <- ev[tolower(ev$domain) == "epidemiological_context", , drop = FALSE]
  if (!nrow(df)) return(NULL)
  key <- .dv_get(df, "outbreak_name")
  key[!nzchar(key)] <- paste0("claim ", df$pmid[!nzchar(key)])
  df$.key <- key
  out <- dplyr::group_by(df, .data$.key) |>
    dplyr::summarise(
      period = paste(unique(c(.dv_nonempty(dplyr::pick(dplyr::everything()), "outbreak_start"),
                              .dv_nonempty(dplyr::pick(dplyr::everything()), "outbreak_end"))), collapse = " - "),
      cases = .dv_context(dplyr::pick(dplyr::everything()), "case_count"),
      deaths = .dv_context(dplyr::pick(dplyr::everything()), "death_count"),
      attack_rate = .dv_context(dplyr::pick(dplyr::everything()), "attack_rate"),
      geography = .dv_context(dplyr::pick(dplyr::everything()), "geographic_scope"),
      affected = .dv_context(dplyr::pick(dplyr::everything()), "affected_population"),
      claims = dplyr::n(),
      pmids = paste(unique(.data$pmid), collapse = ", "),
      .groups = "drop"
    )
  names(out)[1] <- "Outbreak"
  out
}

# Mutation / variant phenotype table.
.dv_variant_table <- function(ev) {
  df <- ev[tolower(ev$domain) == "variant_phenotype", , drop = FALSE]
  if (!nrow(df)) return(NULL)
  out <- data.frame(
    `Variant / mutation` = ifelse(nzchar(.dv_get(df, "variant_name")),
                                  .dv_get(df, "variant_name"),
                                  ifelse(nzchar(.dv_get(df, "mutation")),
                                         .dv_get(df, "mutation"), substr(.dv_get(df, "finding"), 1, 60))),
    Gene = .dv_get(df, "gene"),
    `AA change` = .dv_get(df, "aa_change"),
    Phenotype = .dv_get(df, "phenotype"),
    `Effect` = .dv_get(df, "effect_size"),
    Assay = .dv_get(df, "assay"),
    `Evidence` = .dv_get(df, "evidence_strength"),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  keep <- nzchar(out$`Variant / mutation`) | nzchar(out$Phenotype)
  out[keep, , drop = FALSE]
}

# Frequency table of a field (used for genomic surveillance, routes, hosts).
.dv_freq_table <- function(df, cols) {
  rows <- list()
  for (col in cols) {
    vals <- .dv_multi_values(df, col)
    if (!length(vals)) next
    tab <- sort(table(vals), decreasing = TRUE)
    for (v in names(tab))
      rows[[length(rows) + 1]] <- data.frame(
        Field = ma_metric_label(col), Value = v, Claims = as.integer(tab[[v]]),
        stringsAsFactors = FALSE)
  }
  if (!length(rows)) return(NULL)
  do.call(rbind, rows)
}

# Domain summary table dispatcher.
.dv_domain_table <- function(ev, domain) {
  spec <- DOM_EVIDENCE[[domain]]
  df <- ev[tolower(ev$domain) == domain, , drop = FALSE]
  if (!nrow(df)) return(NULL)
  if (identical(spec$table_fn, "outbreaks")) return(.dv_outbreak_table(ev))
  if (identical(spec$table_fn, "variants")) return(.dv_variant_table(ev))
  if (domain == "genomic_surveillance") {
    t <- .dv_freq_table(df, c("lineage", "clade", "genotype", "sequencing_platform",
                              "sequencing_method", "geographic_distribution",
                              "temporal_distribution"))
    if (!is.null(t)) return(t)
  }
  if (!is.null(spec$group)) {
    t <- .dv_group_metric_table(ev, domain, spec$group, spec$metrics)
    if (!is.null(t) && nrow(t)) return(t)
  }
  t <- .dv_metric_summary_table(ev, domain, spec$metrics %||% character(0))
  if (!is.null(t)) return(t)
  # Fallback: frequency table over the domain's categorical/text fields so
  # metric-light domains (transmission routes, settings…) still have content.
  candidates <- setdiff(names(df), c("claim_id", "run_id", "pmid", "pmcid", "species",
                                     "domain", "topic", "claim_type", "finding", "quote",
                                     "evidence_level", "population", "geography",
                                     "source_file", "details", "loaded_at", "n", ".grp",
                                     "table", names(which(sapply(df, is.numeric)))))
  populated <- candidates[vapply(df[candidates], function(x)
    sum(nzchar(trimws(as.character(x))) & !is.na(x)) > 0, logical(1))]
  .dv_freq_table(df, head(populated, 8))
}

# Domain plot dispatcher.
.dv_domain_plot <- function(ev, domain) {
  spec <- DOM_EVIDENCE[[domain]]
  df <- ev[tolower(ev$domain) == domain, , drop = FALSE]
  if (!nrow(df)) return(NULL)
  if (!is.null(spec$freq_field))
    return(.dv_freq_bar(df, spec$freq_field, paste0(gsub("_", " ", domain), " — by ", gsub("_", " ", spec$freq_field))))
  if (identical(spec$table_fn, "outbreaks")) {
    # Bar of reported case counts per outbreak (descriptive, not pooled).
    vals <- df[nzchar(.dv_get(df, "case_count")), , drop = FALSE]
    if (!nrow(vals)) return(NULL)
    est <- ma_parse_num(.dv_get(vals, "case_count"))
    ok <- is.finite(est)
    if (!sum(ok)) return(NULL)
    lab <- .dv_get(vals, "outbreak_name")
    lab[!nzchar(lab)] <- paste0("pmid:", vals$pmid[!nzchar(lab)])
    dd <- data.frame(label = lab[ok], est = est[ok], stringsAsFactors = FALSE)
    dd <- dd[!duplicated(dd), , drop = FALSE]
    if (!requireNamespace("ggplot2", quietly = TRUE)) return(NULL)
    return(ggplot2::ggplot(dd, ggplot2::aes(x = stats::reorder(label, est), y = est)) +
      ggplot2::geom_col(fill = "#7570b3") +
      ggplot2::coord_flip() +
      ggplot2::labs(x = NULL, y = "reported cases", title = "Cases by outbreak") +
      ggplot2::theme_minimal(base_size = 11))
  }
  if (identical(spec$table_fn, "variants")) return(NULL)
  # default: forest for the first metric with data
  for (field in names(spec$metrics %||% character(0))) {
    ms <- ma_metric_studies(ev, domain, field, spec$metrics[[field]], spec$group)
    if (is.null(ms)) next
    pr <- ma_pool_metric(ms, ev)
    return(ma_forest_gg(ms, pr$pooled, title = paste0(gsub("_", " ", domain), ": ", ma_metric_label(field))))
  }
  NULL
}
