# Evidence meta-analysis helpers for the dashboard.
#
# Turns per-domain extracted fields (columns produced by
# extract_literature_evidence.py / evidence_fields) into pooled estimates,
# forest plots and field-completeness matrices. Everything is per-species:
# callers pass the already-filtered evidence data frame.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a)) b else a

# ---------------------------------------------------------------------------
# Value parsing
# ---------------------------------------------------------------------------

# First (or midpoint of a range) numeric in a string. "84.0%" -> 84, "25-36%" -> 30.5,
# "17 days" -> 17, "482" -> 482.
.ma_first_num <- function(v) {
  if (is.na(v) || !nzchar(trimws(v))) return(NA_real_)
  hits <- regmatches(v, gregexpr("[0-9]+\\.?[0-9]*", v))[[1]]
  if (!length(hits)) return(NA_real_)
  nums <- suppressWarnings(as.numeric(hits))
  nums <- nums[is.finite(nums)]
  if (!length(nums)) return(NA_real_)
  if (length(nums) > 1 && grepl("-|–|\\bto\\b", v)) mean(c(nums[1], nums[length(nums)])) else nums[1]
}

ma_parse_num <- function(x) unname(vapply(as.character(x), .ma_first_num, numeric(1), USE.NAMES = FALSE))

# Interpret a string as a proportion in (0,1): "%" -> /100; bare 84 -> 0.84; 0.84 -> 0.84.
ma_parse_prop <- function(x) {
  s <- as.character(x)
  v <- ma_parse_num(s)
  has_pct <- grepl("%", s)
  out <- ifelse(has_pct, v / 100, ifelse(v > 1 & v <= 100, v / 100, v))
  out[!is.finite(out) | out < 0 | out > 1] <- NA_real_
  out
}

# Wilson score interval for a proportion. Returns c(lo, hi).
.ma_wilson <- function(p, n, z = 1.959964) {
  if (is.na(p) || is.na(n) || n <= 0) return(c(NA_real_, NA_real_))
  den <- 1 + z^2 / n
  centre <- (p + z^2 / (2 * n)) / den
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / den
  c(max(0, centre - half), min(1, centre + half))
}

# ---------------------------------------------------------------------------
# DerSimonian–Laird random-effects pooling
# ---------------------------------------------------------------------------

.ma_dl <- function(est, var) {
  ok <- is.finite(est) & is.finite(var) & var > 0
  est <- est[ok]; var <- var[ok]
  k <- length(est)
  if (k == 0L) return(NULL)
  w <- 1 / var
  fixed <- sum(w * est) / sum(w)
  q <- sum(w * (est - fixed)^2)
  cstat <- sum(w) - sum(w^2) / sum(w)
  tau2 <- if (cstat > 0) max(0, (q - (k - 1)) / cstat) else 0
  wr <- 1 / (var + tau2)
  mu <- sum(wr * est) / sum(wr)
  se <- sqrt(1 / sum(wr))
  i2 <- if (q > k - 1) max(0, (q - (k - 1)) / q) else 0
  list(
    k = k, estimate = mu, se = se,
    ci_lo = mu - 1.959964 * se, ci_hi = mu + 1.959964 * se,
    i2 = i2, tau2 = tau2, q = q,
    p_het = stats::pchisq(q, df = max(k - 1, 1), lower.tail = FALSE)
  )
}

# Pool proportions via logit transform; se_i = 1/(p(1-p)n).
# Boundary estimates (p = 0 or 1) get a continuity correction at 1/(2n).
ma_pool_props <- function(p, n) {
  ok <- is.finite(p) & is.finite(n) & n > 0 & p >= 0 & p <= 1
  p <- p[ok]; n <- n[ok]
  if (!length(p)) return(NULL)
  lo_bound <- 1 / (2 * n)
  p <- pmin(pmax(p, lo_bound), 1 - lo_bound)
  res <- .ma_dl(stats::qlogis(p), 1 / (p * (1 - p) * n))
  if (is.null(res)) return(NULL)
  res$estimate <- stats::plogis(res$estimate)
  res$ci_lo <- stats::plogis(res$ci_lo)
  res$ci_hi <- stats::plogis(res$ci_hi)
  res
}

# Descriptive summary when denominators/SEs are unavailable.
ma_describe <- function(x, weight = NULL) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NULL)
  w <- rep(1, length(x))
  if (!is.null(weight) && length(weight) == length(x)) {
    w <- suppressWarnings(as.numeric(weight))
    w[!is.finite(w) | w <= 0] <- 0.2
  }
  list(
    k = length(x),
    estimate = sum(x * w) / sum(w),
    ci_lo = unname(stats::quantile(x, 0.25, na.rm = TRUE)),
    ci_hi = unname(stats::quantile(x, 0.75, na.rm = TRUE)),
    median = stats::median(x), lo = min(x), hi = max(x),
    i2 = NA_real_
  )
}

# ---------------------------------------------------------------------------
# Domain metric registry
# ---------------------------------------------------------------------------

# Domains describing countermeasures (what was used and how it performed)
MA_COUNTERMEASURE_DOMAINS <- c("vaccine_therapeutic", "diagnostic", "intervention", "policy")

# Domains describing the pathogen itself (what we learned about the strain)
MA_EVIDENCE_DOMAINS <- c(
  "clinical", "epidemiological_context", "transmission", "seroprevalence",
  "reservoir", "genomic_surveillance", "variant_phenotype"
)

MA_DOMAIN_METRICS <- list(
  diagnostic = list(
    list(field = "sensitivity", type = "proportion", group = "test_name"),
    list(field = "specificity", type = "proportion", group = "test_name"),
    list(field = "positive_predictive_value", type = "proportion", group = "test_name"),
    list(field = "negative_predictive_value", type = "proportion", group = "test_name"),
    list(field = "limit_of_detection", type = "value", group = "test_name")
  ),
  clinical = list(
    list(field = "case_fatality_rate", type = "proportion"),
    list(field = "hospitalization_rate", type = "proportion"),
    list(field = "incubation_period", type = "value")
  ),
  epidemiological_context = list(
    list(field = "attack_rate", type = "proportion"),
    list(field = "incidence", type = "proportion"),
    list(field = "prevalence", type = "proportion"),
    list(field = "case_count", type = "total"),
    list(field = "death_count", type = "total")
  ),
  transmission = list(
    list(field = "r0", type = "value"),
    list(field = "serial_interval", type = "value"),
    list(field = "generation_time", type = "value"),
    list(field = "secondary_attack_rate", type = "proportion")
  ),
  seroprevalence = list(
    list(field = "seroprevalence_value", type = "proportion", group = "host_species")
  ),
  reservoir = list(
    list(field = "prevalence", type = "proportion", group = "host_species")
  ),
  vaccine_therapeutic = list(
    list(field = "efficacy", type = "proportion", group = "product_name"),
    list(field = "effectiveness", type = "proportion", group = "product_name")
  ),
  intervention = list(
    list(field = "effectiveness", type = "proportion", group = "intervention_name"),
    list(field = "coverage", type = "proportion", group = "intervention_name")
  ),
  genomic_surveillance = list(
    list(field = "sequence_count", type = "total")
  ),
  policy = list(),
  variant_phenotype = list()
)

ma_metric_label <- function(field) {
  paste0(toupper(substring(field, 1, 1)), gsub("_", " ", substring(field, 2)))
}

# Identity-field fallbacks: product_name -> product_type, etc.
.ma_group_fallbacks <- list(
  product_name = c("product_type", "target"),
  test_name = c("test_type", "target_gene_or_antigen", "target"),
  intervention_name = c("intervention_type"),
  policy_name = c("issuing_body", "policy_type"),
  host_species = c("host_type", "host_population")
)

# ---------------------------------------------------------------------------
# Study-level extraction for one metric
# ---------------------------------------------------------------------------

# Returns list(studies = df, kind = "proportion"|"value"|"total").
# studies: label, est, lo, hi, n, pmid, group.
ma_metric_studies <- function(ev, domain, field, type, group = NULL) {
  if (is.null(ev) || !nrow(ev)) return(NULL)
  df <- ev[tolower(ev$domain) == domain, , drop = FALSE]
  if (!nrow(df) || !(field %in% names(df))) return(NULL)

  raw <- df[[field]]
  keep <- !is.na(raw) & nzchar(trimws(as.character(raw)))
  df <- df[keep, , drop = FALSE]
  raw <- raw[keep]
  if (!nrow(df)) return(NULL)

  grp <- if (!is.null(group) && group %in% names(df)) as.character(df[[group]]) else ""
  # Fallback chain when the primary identity field is empty (extractor often
  # populates product_type but not a named product).
  grp_fb <- if (!is.null(group)) .ma_group_fallbacks[[group]] else NULL
  for (fb in grp_fb %||% character(0)) {
    if (!(fb %in% names(df))) next
    need <- is.na(grp) | !nzchar(grp)
    if (!any(need)) break
    alt <- as.character(df[[fb]])
    grp[need & !is.na(alt) & nzchar(alt)] <- alt[need & !is.na(alt) & nzchar(alt)]
  }
  grp[is.na(grp) | !nzchar(grp)] <- "unspecified"

  # denominator: n (all domains), then sample_size (seroprevalence / genomic)
  n <- suppressWarnings(ma_parse_num(if ("n" %in% names(df)) df$n else NA))
  if (all(is.na(n)) && "sample_size" %in% names(df)) {
    n <- suppressWarnings(ma_parse_num(df$sample_size))
  }

  if (type == "proportion") {
    est <- ma_parse_prop(raw)
  } else {
    est <- ma_parse_num(raw)
  }
  ok <- is.finite(est)
  df <- df[ok, , drop = FALSE]; est <- est[ok]; n <- n[ok]; grp <- grp[ok]
  if (!length(est)) return(NULL)

  pmid <- as.character(df$pmid)
  pmid[is.na(pmid) | !nzchar(pmid)] <- "unknown"
  label <- paste0(pmid, ifelse(nzchar(grp), paste0(" | ", grp), ""))
  # Disambiguate repeats of the same pmid+group (e.g. multiple sample types)
  if (anyDuplicated(label)) {
    extra <- if ("sample_type" %in% names(df)) as.character(df$sample_type)
             else if ("finding" %in% names(df)) substr(as.character(df$finding), 1, 30)
             else ""
    extra[is.na(extra)] <- ""
    label <- make.unique(paste0(label, ifelse(nzchar(extra), paste0(" | ", extra), "")), sep = " ")
  }
  studies <- data.frame(
    label = label, pmid = pmid, group = grp,
    est = est, n = n, stringsAsFactors = FALSE
  )

  if (type == "proportion") {
    ci <- t(vapply(seq_along(est), function(i) .ma_wilson(est[i], n[i]), numeric(2)))
    studies$lo <- ci[, 1]; studies$hi <- ci[, 2]
  } else {
    studies$lo <- NA_real_; studies$hi <- NA_real_
  }
  list(studies = studies, kind = type, field = field, domain = domain, group = group)
}

# ---------------------------------------------------------------------------
# Pooled estimate for one metric spec
# ---------------------------------------------------------------------------

ma_pool_metric <- function(ms, ev) {
  if (is.null(ms) || !nrow(ms$studies)) return(NULL)
  s <- ms$studies
  weights <- vapply(ev$evidence_level[match(s$pmid, ev$pmid)] %||% "unknown",
                    function(l) cm_score_evidence(l), numeric(1))
  if (ms$kind == "proportion") {
    pooled <- ma_pool_props(s$est, s$n)
    if (is.null(pooled)) pooled <- ma_describe(s$est, weights)
    return(list(pooled = pooled, pooled_kind = if (!is.null(pooled$i2) && is.finite(pooled$i2)) "meta" else "descriptive"))
  }
  if (ms$kind == "total") {
    return(list(pooled = list(k = nrow(s), estimate = sum(s$est, na.rm = TRUE),
                              ci_lo = min(s$est, na.rm = TRUE), ci_hi = max(s$est, na.rm = TRUE),
                              i2 = NA_real_), pooled_kind = "total"))
  }
  list(pooled = ma_describe(s$est, weights), pooled_kind = "descriptive")
}

# ---------------------------------------------------------------------------
# Summary of every metric across every domain
# ---------------------------------------------------------------------------

ma_all_metric_summary <- function(ev, domains = names(MA_DOMAIN_METRICS)) {
  if (is.null(ev) || !nrow(ev)) return(NULL)
  rows <- list()
  for (domain in domains) {
    for (spec in MA_DOMAIN_METRICS[[domain]] %||% list()) {
      ms <- ma_metric_studies(ev, domain, spec$field, spec$type, spec$group)
      if (is.null(ms)) next
      pr <- ma_pool_metric(ms, ev)
      p <- pr$pooled
      if (is.null(p)) next
      rows[[length(rows) + 1]] <- data.frame(
        domain = gsub("_", " ", domain),
        metric = ma_metric_label(spec$field),
        k = p$k,
        pooled = p$estimate,
        lo = p$ci_lo, hi = p$ci_hi,
        i2 = if (!is.null(p$i2)) p$i2 else NA_real_,
        method = pr$pooled_kind,
        type = spec$type,
        stringsAsFactors = FALSE
      )
    }
  }
  if (!length(rows)) return(NULL)
  do.call(rbind, rows)
}

# ---------------------------------------------------------------------------
# Forest plot
# ---------------------------------------------------------------------------

ma_forest_gg <- function(ms, pooled, title = NULL) {
  if (is.null(ms) || !nrow(ms$studies)) return(NULL)
  s <- ms$studies[order(ms$studies$est), ]
  s$y <- seq_len(nrow(s))
  is_prop <- ms$kind == "proportion"

  p <- ggplot(s, aes(x = est, y = y))
  if (is_prop && any(is.finite(s$lo) & is.finite(s$hi))) {
    p <- p + geom_segment(aes(x = lo, xend = hi, y = y, yend = y), colour = "grey60", na.rm = TRUE)
  }
  p <- p +
    geom_point(size = 2.4, colour = "#2c7fb8") +
    scale_y_continuous(breaks = s$y, labels = s$label) +
    labs(x = ma_metric_label(ms$field), y = NULL,
         title = title %||% paste0(ms$domain, ": ", ms$field)) +
    theme_minimal(base_size = 12)

  if (is_prop) p <- p + scale_x_continuous(labels = scales::percent_format(accuracy = 1))

  if (!is.null(pooled) && is.finite(pooled$estimate)) {
    p <- p +
      geom_vline(xintercept = pooled$estimate, linetype = 2, colour = "#d95f02") +
      annotate("point", x = pooled$estimate, y = max(s$y) + 0.8,
               shape = 18, size = 4, colour = "#d95f02")
    if (all(is.finite(c(pooled$ci_lo, pooled$ci_hi)))) {
      p <- p +
        annotate("segment", x = pooled$ci_lo, xend = pooled$ci_hi,
                 y = max(s$y) + 0.8, yend = max(s$y) + 0.8, colour = "#d95f02") +
        annotate("segment", x = pooled$ci_lo, xend = pooled$ci_lo,
                 y = max(s$y) + 0.6, yend = max(s$y) + 1.0, colour = "#d95f02") +
        annotate("segment", x = pooled$ci_hi, xend = pooled$ci_hi,
                 y = max(s$y) + 0.6, yend = max(s$y) + 1.0, colour = "#d95f02") +
        annotate("text", x = pooled$estimate, y = max(s$y) + 0.8,
                 label = " pooled", hjust = -0.15, vjust = -0.8, colour = "#d95f02", size = 3.5)
    }
  }
  p
}

# ---------------------------------------------------------------------------
# Field completeness matrix (for gap scoring)
# ---------------------------------------------------------------------------

ma_field_completeness <- function(ev, domains = names(MA_DOMAIN_METRICS)) {
  if (is.null(ev) || !nrow(ev)) return(NULL)
  out <- list()
  for (domain in domains) {
    specs <- MA_DOMAIN_METRICS[[domain]]
    if (!length(specs)) next
    df <- ev[tolower(ev$domain) == domain, , drop = FALSE]
    total <- nrow(df)
    for (spec in specs) {
      f <- spec$field
      if (!(f %in% names(df))) {
        frac <- 0
      } else {
        vals <- as.character(df[[f]])
        frac <- if (total) sum(!is.na(vals) & nzchar(trimws(vals))) / total else 0
      }
      out[[length(out) + 1]] <- data.frame(
        domain = domain, field = f, claims = total, completeness = frac,
        stringsAsFactors = FALSE
      )
    }
  }
  if (!length(out)) return(NULL)
  do.call(rbind, out)
}

# ---------------------------------------------------------------------------
# Per-group pooled estimates (e.g. pooled sensitivity per diagnostic test)
# ---------------------------------------------------------------------------

ma_domain_specs <- function(domain) MA_DOMAIN_METRICS[[domain]] %||% list()

# Pool a single metric within each value of `group` (test_name, product_name…).
# Returns df: group, k, est, lo, hi, i2, method.
ma_group_pooled <- function(ev, domain, field, group = NULL, type = "proportion") {
  ms <- ma_metric_studies(ev, domain, field, type, group)
  if (is.null(ms) || !nrow(ms$studies)) return(NULL)
  s <- ms$studies
  out <- lapply(unique(s$group), function(g) {
    sg <- s[s$group == g, , drop = FALSE]
    pooled <- if (type == "proportion") {
      p <- ma_pool_props(sg$est, sg$n)
      if (is.null(p)) ma_describe(sg$est) else p
    } else {
      ma_describe(sg$est)
    }
    if (is.null(pooled)) return(NULL)
    data.frame(
      group = g, k = nrow(sg), est = pooled$estimate,
      lo = pooled$ci_lo, hi = pooled$ci_hi,
      i2 = if (!is.null(pooled$i2)) pooled$i2 else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  out <- out[!vapply(out, is.null, logical(1))]
  if (!length(out)) return(NULL)
  do.call(rbind, out)
}

# All metric specs that actually have data, for one or more domains.
ma_available_specs <- function(ev, domains = names(MA_DOMAIN_METRICS)) {
  specs <- list()
  if (is.null(ev) || !nrow(ev)) return(specs)
  for (domain in domains) {
    for (spec in ma_domain_specs(domain)) {
      ms <- ma_metric_studies(ev, domain, spec$field, spec$type, spec$group)
      if (!is.null(ms) && nrow(ms$studies) >= 1) {
        specs[[paste0(domain, "|", spec$field)]] <- spec
      }
    }
  }
  specs
}

# ---------------------------------------------------------------------------
# Pivot evidence_fields (long) to wide columns for the DuckDB path
# ---------------------------------------------------------------------------

# Fallback: expand the `details` JSON column of evidence_extracted into columns.
ma_details_wide <- function(ev) {
  if (is.null(ev) || !("details" %in% names(ev)) ||
      !requireNamespace("jsonlite", quietly = TRUE)) return(NULL)
  rows <- lapply(seq_len(nrow(ev)), function(i) {
    det <- ev$details[i]
    if (is.na(det) || !nzchar(det) || det == "{}") return(NULL)
    d <- tryCatch(jsonlite::fromJSON(det, simplifyVector = TRUE), error = function(e) NULL)
    if (is.null(d) || !length(d)) return(NULL)
    data.frame(
      claim_id = if ("claim_id" %in% names(ev)) ev$claim_id[i] else i,
      field_name = names(d),
      field_value = vapply(d, function(x) paste(as.character(x), collapse = "; "), character(1)),
      stringsAsFactors = FALSE
    )
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) return(NULL)
  ma_fields_wide(do.call(rbind, rows))
}

ma_fields_wide <- function(fld) {
  if (is.null(fld) || !nrow(fld)) return(NULL)
  ids <- unique(fld$claim_id)
  fns <- unique(fld$field_name)
  m <- matrix(NA_character_, nrow = length(ids), ncol = length(fns),
              dimnames = list(NULL, fns))
  idx <- match(fld$claim_id, ids)
  for (i in seq_len(nrow(fld))) {
    v <- fld$field_value[i]
    m[idx[i], fld$field_name[i]] <- if (is.na(v) || !nzchar(v) || v == "null") NA_character_ else v
  }
  data.frame(claim_id = ids, m, check.names = FALSE, stringsAsFactors = FALSE)
}
