# Intelligence Brief module — shared between pipeline and dashboard
#
# Pipeline side (driven by bin/generate_intelligence_brief.R inside
# modules/local/generate_intelligence_brief): assembles a structured fact
# list from the species' pipeline outputs, calls the local Ollama model via
# dashboard/modules/llm_note.R for the bottom-line narrative, and falls back
# to a deterministic template — intelligence_brief.json is always produced.
#
# Dashboard side (sourced by app.R): reads
# results/intelligence_brief/<species>/intelligence_brief.json and renders
# the compact front-page summary (Intelligence Overview card) and the full
# Intelligence Brief tab.
#
# The JSON is the versioned artifact; the dashboard only renders it. The
# MONITOR/INVESTIGATE/ESCALATE verdict stays dashboard-side in
# modules/assessment.R (compute_pi_assessment).

# ---------------------------------------------------------------------------
# Tolerant readers (NULL on missing/empty — every section degrades gracefully)
# ---------------------------------------------------------------------------

if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a)) b else a

.brief_tsv <- function(dir, fname) {
  if (is.null(dir) || !nzchar(dir)) return(NULL)
  p <- if (grepl("\\.tsv$", dir)) dir else file.path(dir, fname)
  if (!file.exists(p)) return(NULL)
  df <- tryCatch(readr::read_tsv(p, show_col_types = FALSE), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df
}

.brief_json <- function(p) {
  if (is.null(p) || !nzchar(p) || !file.exists(p) ||
      !requireNamespace("jsonlite", quietly = TRUE)) return(NULL)
  tryCatch(jsonlite::fromJSON(p, simplifyVector = FALSE), error = function(e) NULL)
}

.brief_num  <- function(x) { v <- suppressWarnings(as.numeric(x)); v[!is.na(v)] }
.brief_med  <- function(x) { v <- .brief_num(x); if (!length(v)) return(NULL); stats::median(v) }
.brief_rng  <- function(x) { v <- .brief_num(x); if (!length(v)) return(NULL); c(min = min(v), max = max(v)) }
.brief_uniq <- function(x) { x <- as.character(x); unique(x[!is.na(x) & nzchar(x)]) }
.brief_mode <- function(x) { u <- .brief_uniq(x); if (!length(u)) return(NULL); names(which.max(table(x[!is.na(x) & nzchar(x)]))) }
.brief_dates_rng <- function(x) {
  d <- suppressWarnings(as.Date(x)); d <- d[!is.na(d)]
  if (!length(d)) return(NULL)
  list(min = as.character(min(d)), max = as.character(max(d)))
}
# First non-null/non-NA scalar of a column
.brief_gv <- function(df, col) {
  if (is.null(df) || !(col %in% names(df)) || !nrow(df)) return(NULL)
  v <- df[[col]][1]
  if (length(v) == 0 || (length(v) == 1 && is.na(v))) NULL else v
}
.brief_drop_null <- function(x) x[!vapply(x, function(e) is.null(e) || (length(e) == 1 && is.na(e)), logical(1))]

# Recursively prune NULLs and empty lists so absent values serialize as
# missing keys (jsonlite turns list(NULL) elements into {} otherwise, which
# then reads back as empty list rather than NULL).
.brief_prune <- function(x) {
  if (!is.list(x)) return(x)
  x <- lapply(x, .brief_prune)
  x[!vapply(x, function(e) is.null(e) || (is.list(e) && length(e) == 0), logical(1))]
}

# ---------------------------------------------------------------------------
# Fact assembly — the "what the workflow found" layer
# ---------------------------------------------------------------------------
# All inputs optional; missing pieces are recorded in data_gaps.
.pg_brief_facts <- function(species,
                            pi_dir = NULL, tc_dir = NULL, mp_dir = NULL,
                            tree_note_path = NULL, spread_json_path = NULL,
                            duckdb_path = NULL, lit_dir = NULL) {
  gaps <- character()

  # ---- Pathogen identification ----
  id_sum      <- .brief_tsv(pi_dir, "identification_summary.tsv")
  unresolved  <- .brief_tsv(pi_dir, "unresolved_samples.tsv")
  if (is.null(id_sum)) gaps <- c(gaps, "pathogen_identification")

  situation <- list(
    n_query_genomes      = if (!is.null(id_sum)) nrow(id_sum) else NULL,
    sample_ids           = if (!is.null(id_sum)) .brief_uniq(id_sum$sample) else NULL,
    species_identified   = if (!is.null(id_sum)) .brief_mode(id_sum$species_identified) else NULL,
    median_coverage      = if (!is.null(id_sum)) .brief_med(id_sum$coverage) else NULL
  )

  identification <- list(
    closest_reference        = if (!is.null(id_sum)) .brief_mode(id_sum$closest_reference) else NULL,
    pct_identity_median      = if (!is.null(id_sum)) .brief_med(id_sum$pct_identity_closest_ref) else NULL,
    pct_identity_range       = if (!is.null(id_sum)) .brief_rng(id_sum$pct_identity_closest_ref) else NULL,
    genetic_distance_median  = if (!is.null(id_sum)) .brief_med(id_sum$genetic_distance_closest_ref) else NULL,
    genetic_distance_max     = if (!is.null(id_sum)) { v <- .brief_num(id_sum$genetic_distance_closest_ref); if (length(v)) max(v) else NULL } else NULL,
    distance_to_clade_median = if (!is.null(id_sum)) .brief_med(id_sum$distance_to_assigned_clade) else NULL,
    closest_outbreak_match   = if (!is.null(id_sum)) .brief_uniq(id_sum$closest_outbreak) else NULL,
    n_unresolved             = if (!is.null(unresolved)) nrow(unresolved) else NULL,
    unresolved_ids           = if (!is.null(unresolved)) .brief_uniq(unresolved$sample) else NULL,
    unresolved_assignments   = if (!is.null(unresolved) && "assignment" %in% names(unresolved)) .brief_uniq(unresolved$assignment) else NULL
  )

  # ---- Transmission context ----
  potential <- .brief_tsv(tc_dir, "transmission_potential.tsv")
  score     <- .brief_tsv(tc_dir, "transmission_potential_score.tsv")
  src       <- .brief_tsv(tc_dir, "query_source_inference.tsv")
  cprof     <- .brief_tsv(tc_dir, "query_cluster_profile.tsv")
  clusters  <- .brief_tsv(tc_dir, "genetic_clusters.tsv")
  qctx      <- .brief_tsv(tc_dir, "query_epi_context.tsv")
  proj      <- .brief_tsv(tc_dir, "query_projection.tsv")
  strain_p  <- .brief_tsv(tc_dir, "strain_transmission_profile.tsv")
  obsum     <- .brief_tsv(tc_dir, "transmission_outbreak_summary.tsv")
  spatial   <- .brief_tsv(tc_dir, "transmission_spatial.tsv")
  if (is.null(potential)) gaps <- c(gaps, "transmission_context")

  # situation geography/dates — transmission tables carry them when PI doesn't
  q_countries <- unique(c(
    if (!is.null(potential)) .brief_uniq(potential$query_country) else character(),
    if (!is.null(qctx))      .brief_uniq(qctx$query_country) else character()
  ))
  q_admin1 <- unique(c(
    if (!is.null(potential)) .brief_uniq(potential$query_admin1) else character(),
    if (!is.null(qctx))      .brief_uniq(qctx$query_admin1) else character()
  ))
  q_dates <- c(
    if (!is.null(potential)) potential$query_collection_date else character(),
    if (!is.null(qctx))      qctx$query_collection_date else character()
  )
  situation$query_countries <- if (length(q_countries)) q_countries else NULL
  situation$query_admin1    <- if (length(q_admin1)) q_admin1 else NULL
  situation$collection_date_range <- .brief_dates_rng(q_dates)
  situation$query_strains   <- if (!is.null(potential)) .brief_uniq(potential$query_strain) else NULL

  transmission <- list(
    n_queries              = if (!is.null(potential)) length(.brief_uniq(potential$query_sample)) else NULL,
    risk_labels            = if (!is.null(score)) .brief_uniq(score$risk_label) else NULL,
    predicted_behaviors    = if (!is.null(score)) .brief_uniq(score$predicted_behavior) else NULL,
    likely_origins         = unique(c(
      if (!is.null(src))   .brief_uniq(src$likely_origin_country) else character(),
      if (!is.null(score)) .brief_uniq(score$likely_source_countries) else character()
    )),
    n_introductions_range  = if (!is.null(src)) .brief_rng(src$n_introductions_into_query_country) else NULL,
    n_query_clusters       = if (!is.null(cprof)) length(.brief_uniq(cprof$cluster_id)) else NULL,
    largest_cluster_size   = if (!is.null(cprof)) { v <- .brief_num(cprof$cluster_size); if (length(v)) max(v) else NULL } else NULL,
    cluster_countries      = if (!is.null(clusters)) .brief_uniq(clusters$country) else NULL,
    cluster_date_range     = if (!is.null(clusters)) .brief_dates_rng(clusters$tip_date) else NULL,
    n_neighbors_with_epi_min = if (!is.null(score)) { v <- .brief_num(score$n_neighbors_with_epi); if (length(v)) min(v) else NULL } else NULL,
    r_at_sampling_median   = if (!is.null(qctx)) .brief_med(qctx$r_at_sampling) else NULL,
    r_at_sampling_range    = if (!is.null(qctx)) .brief_rng(qctx$r_at_sampling) else NULL,
    growth_phases          = if (!is.null(qctx)) .brief_uniq(qctx$growth_phase) else NULL,
    weekly_cases_at_sampling = if (!is.null(qctx)) .brief_med(qctx$weekly_cases_at_sampling) else NULL,
    cfr_at_sampling        = if (!is.null(qctx)) .brief_med(qctx$cfr_at_sampling) else NULL,
    sampling_fraction      = if (!is.null(qctx)) .brief_med(qctx$sampling_fraction) else NULL
  )
  if (length(transmission$likely_origins) == 0) transmission$likely_origins <- NULL

  # ---- Distribution / lineage history ----
  query_strains <- situation$query_strains %||% character()
  strain_rows <- NULL
  if (!is.null(strain_p)) {
    keep <- rep(TRUE, nrow(strain_p))
    if (length(query_strains) && "strain" %in% names(strain_p))
      keep <- strain_p$strain %in% query_strains
    # Query-matching strains first, then the most impactful background strains
    sr <- strain_p[keep, , drop = FALSE]
    other <- strain_p[!keep, , drop = FALSE]
    if (nrow(other) && "linked_cases" %in% names(other)) {
      other <- other[order(-suppressWarnings(as.numeric(other$linked_cases))), , drop = FALSE]
      sr <- rbind(sr, utils::head(other, 3))
    }
    strain_rows <- sr
  }
  distribution <- list(
    query_countries = situation$query_countries,
    strain_profiles = if (!is.null(strain_rows) && nrow(strain_rows)) {
      lapply(seq_len(nrow(strain_rows)), function(i) .brief_drop_null(list(
        strain              = .brief_gv(strain_rows[i, ], "strain"),
        first_seen          = .brief_gv(strain_rows[i, ], "first_seen"),
        last_seen           = .brief_gv(strain_rows[i, ], "last_seen"),
        active_years        = .brief_gv(strain_rows[i, ], "active_years"),
        n_countries         = .brief_gv(strain_rows[i, ], "n_countries"),
        countries           = .brief_gv(strain_rows[i, ], "countries"),
        origin_country      = .brief_gv(strain_rows[i, ], "origin_country"),
        spread_rate_km_per_year = .brief_gv(strain_rows[i, ], "spread_rate_km_per_year"),
        linked_cases        = .brief_gv(strain_rows[i, ], "linked_cases"),
        linked_deaths       = .brief_gv(strain_rows[i, ], "linked_deaths"),
        mean_cfr            = .brief_gv(strain_rows[i, ], "mean_cfr"),
        behavior_label      = .brief_gv(strain_rows[i, ], "behavior_label")
      )))
    } else NULL,
    epi_countries = if (!is.null(spatial)) .brief_uniq(spatial$country) else NULL
  )

  # ---- Historical outbreak context ----
  historical_outbreaks <- if (!is.null(obsum)) {
    lapply(seq_len(nrow(obsum)), function(i) .brief_drop_null(list(
      start_year      = .brief_gv(obsum[i, ], "start_year"),
      country         = .brief_gv(obsum[i, ], "country"),
      admin1          = .brief_gv(obsum[i, ], "admin1"),
      cases           = .brief_gv(obsum[i, ], "cases"),
      deaths          = .brief_gv(obsum[i, ], "deaths"),
      cfr             = .brief_gv(obsum[i, ], "cfr"),
      species_subtype = .brief_gv(obsum[i, ], "species_subtype")
    )))
  } else NULL

  # ---- Outlook: final-horizon week per query, median across queries per scenario ----
  outlook <- NULL
  if (!is.null(proj) && all(c("scenario", "week_ahead", "proj_median") %in% names(proj))) {
    by_scenario <- split(proj, proj$scenario)
    outlook <- list(scenarios = lapply(by_scenario, function(d) {
      d <- d[d$week_ahead == max(suppressWarnings(as.numeric(d$week_ahead)), na.rm = TRUE), , drop = FALSE]
      .brief_drop_null(list(
        week_ahead = .brief_gv(d, "week_ahead"),
        median     = .brief_med(d$proj_median),
        lower95    = .brief_med(d$proj_lower95),
        upper95    = .brief_med(d$proj_upper95),
        r_used     = .brief_med(d$r_used)
      ))
    }))
  }

  # ---- Genomics / mutation profile ----
  burden_sum <- .brief_tsv(mp_dir, "01_mutation_burden_summary.tsv")
  pheno      <- .brief_tsv(mp_dir, "01_mutation_phenotypes.tsv")
  vip        <- .brief_tsv(mp_dir, "01_plsda_vip.tsv")
  pos_sum    <- .brief_tsv(mp_dir, "01_position_summary.tsv")
  if (is.null(burden_sum) && is.null(pheno) && is.null(vip)) gaps <- c(gaps, "mutation_profile")

  burden_q  <- if (!is.null(burden_sum) && "is_query" %in% names(burden_sum)) burden_sum[burden_sum$is_query %in% c(TRUE, "TRUE", "true", 1), , drop = FALSE] else NULL
  burden_bg <- if (!is.null(burden_sum) && "is_query" %in% names(burden_sum)) burden_sum[!(burden_sum$is_query %in% c(TRUE, "TRUE", "true", 1)), , drop = FALSE] else NULL

  vip_top <- NULL
  if (!is.null(vip) && all(c("protein_name", "vip") %in% names(vip))) {
    vip <- vip[order(-suppressWarnings(as.numeric(vip$vip))), , drop = FALSE]
    vip_top <- utils::head(paste0(vip$protein_name), 3)
  }
  pos_query <- NULL
  if (!is.null(pos_sum) && "has_query_mutation" %in% names(pos_sum)) {
    qp <- pos_sum[pos_sum$has_query_mutation %in% c(TRUE, "TRUE", "true", 1), , drop = FALSE]
    if (nrow(qp) && "shannon_entropy" %in% names(qp))
      qp <- qp[order(-suppressWarnings(as.numeric(qp$shannon_entropy))), , drop = FALSE]
    pos_query <- list(
      n_positions = nrow(qp),
      top = if (nrow(qp)) utils::head(paste0(qp$protein_name, ":", qp$position), 5) else NULL
    )
  }

  tn <- .brief_json(tree_note_path)
  genomics <- list(
    query_burden_median      = if (!is.null(burden_q) && nrow(burden_q)) .brief_gv(burden_q, "median_burden") else NULL,
    background_burden_median = if (!is.null(burden_bg) && nrow(burden_bg)) .brief_gv(burden_bg, "median_burden") else NULL,
    burden_test_p            = if (!is.null(burden_q) && nrow(burden_q)) .brief_gv(burden_q, "p_value") else NULL,
    n_phenotype_mutations    = if (!is.null(pheno)) nrow(pheno) else 0,
    phenotype_mutations      = if (!is.null(pheno)) {
      lbl <- if ("mutation_label" %in% names(pheno)) pheno$mutation_label else pheno$mutation_id
      ph  <- if ("phenotype" %in% names(pheno)) pheno$phenotype else rep("", length(lbl))
      utils::head(paste0(lbl, ifelse(nzchar(ph), paste0(" (", ph, ")"), "")), 5)
    } else NULL,
    top_vip_proteins         = vip_top,
    query_mutated_positions  = pos_query,
    tree_note_text           = if (!is.null(tn)) tn$text else NULL,
    tree_note_source         = if (!is.null(tn)) tn$source else NULL
  )

  # ---- Per-query spread assessments (pre-generated narratives) ----
  sa <- .brief_json(spread_json_path)
  spread_assessments <- if (!is.null(sa) && !is.null(sa$assessments)) sa$assessments else NULL
  if (is.null(spread_assessments)) gaps <- c(gaps, "spread_assessment")

  # ---- Evidence base: literature domain rollups. Preferred source is the
  # DuckDB warehouse (literature_domains); when those tables are empty — e.g.
  # literature evidence isn't ingested into the warehouse in this run — fall
  # back to counting the published per-domain metadata JSONs under
  # literature_retrieval/literature_metadata/<species>/<domain>/.
  evidence_base <- NULL
  if (!is.null(duckdb_path) && file.exists(duckdb_path) &&
      requireNamespace("DBI", quietly = TRUE) && requireNamespace("duckdb", quietly = TRUE)) {
    con <- tryCatch(DBI::dbConnect(duckdb::duckdb(), dbdir = duckdb_path, read_only = TRUE),
                    error = function(e) NULL)
    if (!is.null(con)) {
      dom <- tryCatch(DBI::dbGetQuery(con,
        "SELECT domain, SUM(total_papers) AS n_papers, SUM(clean_count) AS n_clean
         FROM literature_domains WHERE lower(species) = lower(?) GROUP BY domain",
        params = list(species)), error = function(e) NULL)
      tryCatch(DBI::dbDisconnect(con, shutdown = TRUE), error = function(e) NULL)
      if (!is.null(dom) && nrow(dom)) {
        evidence_base <- list(
          domains     = stats::setNames(as.list(dom$n_papers), dom$domain),
          clean       = stats::setNames(as.list(dom$n_clean), dom$domain),
          total_papers = sum(dom$n_papers, na.rm = TRUE)
        )
      }
    }
  }
  if (is.null(evidence_base) && !is.null(lit_dir) && dir.exists(lit_dir)) {
    doms <- basename(list.dirs(lit_dir, recursive = FALSE))
    counts <- vapply(doms, function(d)
      length(list.files(file.path(lit_dir, d), pattern = "\\.json$")),
      integer(1))
    if (length(doms) && sum(counts) > 0)
      evidence_base <- list(
        domains      = as.list(stats::setNames(counts, doms)),
        total_papers = sum(counts)
      )
  }
  if (is.null(evidence_base)) gaps <- c(gaps, "literature_evidence")

  .brief_drop_null(list(
    situation = situation, identification = identification,
    genomics = genomics, transmission = transmission,
    distribution = distribution, historical_outbreaks = historical_outbreaks,
    outlook = outlook, evidence_base = evidence_base,
    spread_assessments = spread_assessments,
    data_gaps = if (length(gaps)) unique(gaps) else NULL
  ))
}

# ---------------------------------------------------------------------------
# Caveats — deterministic, transparent, shown verbatim in the brief
# ---------------------------------------------------------------------------
.pg_brief_caveats <- function(f) {
  cv <- character()
  idn <- f$identification; tr <- f$transmission; gn <- f$genomics
  n_q <- f$situation$n_query_genomes %||% 0
  if (!is.null(idn$n_unresolved) && idn$n_unresolved > 0)
    cv <- c(cv, sprintf("%d of %d query genomes did not pass identification QC (%s) — species identification is provisional.",
                        idn$n_unresolved, n_q,
                        paste(idn$unresolved_assignments %||% "REJECTED", collapse = "/")))
  cov <- f$situation$median_coverage
  if (!is.null(cov) && cov < 0.9)
    cv <- c(cv, sprintf("Median genome coverage is low (%.0f%%).", cov * 100))
  if (!is.null(tr$n_neighbors_with_epi_min) && tr$n_neighbors_with_epi_min == 0)
    cv <- c(cv, "Some queries have no neighbours with linked epidemiological data — projections for them are unreliable.")
  if (length(tr$growth_phases) && all(tr$growth_phases == "unknown"))
    cv <- c(cv, "Growth phase at sampling is unknown for all queries — projections assume baseline behaviour.")
  if (is.null(gn) || (is.null(gn$n_phenotype_mutations) || gn$n_phenotype_mutations == 0))
    cv <- c(cv, "No phenotype-linked mutation evidence found for query genomes.")
  gap_msgs <- c(
    pathogen_identification = "Pathogen identification outputs were not generated this run.",
    transmission_context    = "Transmission-context tables were not generated this run — spread and outlook sections are empty.",
    mutation_profile        = "Mutation-profile outputs were not generated this run.",
    spread_assessment       = "Per-query spread assessments were not generated this run.",
    literature_evidence     = "No literature-evidence counts are available for this species."
  )
  for (g in f$data_gaps %||% character())
    cv <- c(cv, gap_msgs[[g]] %||% sprintf("Section not generated this run: %s.", gsub("_", " ", g)))
  cv
}

# ---------------------------------------------------------------------------
# Bottom-line narrative — Ollama first, deterministic template fallback
# ---------------------------------------------------------------------------
.pg_brief_template <- function(f) {
  s <- f$situation; idn <- f$identification; tr <- f$transmission; out <- f$outlook
  parts <- character()
  sp_name <- s$species_identified %||% "an uncharacterized pathogen"

  loc <- c(s$query_admin1, s$query_countries)
  loc <- loc[nzchar(loc)]
  dr <- s$collection_date_range
  n_q <- s$n_query_genomes %||% 0
  parts <- c(parts, sprintf(
    "%d query genome%s were identified as %s%s%s.",
    n_q, if (n_q == 1) "" else "s", sp_name,
    if (length(loc)) paste0(" collected in ", paste(loc, collapse = ", ")) else "",
    if (!is.null(dr)) paste0(" (", dr$min, " to ", dr$max, ")") else ""
  ))
  if (!is.null(idn$closest_reference)) {
    parts <- c(parts, sprintf(
      "The closest known reference is %s (median identity %.1f%%, median genetic distance %.4f)%s.",
      idn$closest_reference,
      (idn$pct_identity_median %||% NA_real_) * 100,
      idn$genetic_distance_median %||% NA_real_,
      if (length(idn$closest_outbreak_match))
        paste0("; nearest historical outbreak match: ", paste(idn$closest_outbreak_match, collapse = ", "))
      else ""
    ))
  }
  if (!is.null(tr$n_query_clusters)) {
    parts <- c(parts, sprintf(
      "The queries fall into %d genetic cluster%s%s%s.",
      tr$n_query_clusters, if (tr$n_query_clusters == 1) "" else "s",
      if (length(tr$cluster_countries)) paste0(" spanning ", paste(utils::head(tr$cluster_countries, 4), collapse = ", ")) else "",
      if (!is.null(tr$n_introductions_range))
        sprintf("; %s introduction%s into the query country inferred (likely origin: %s)",
                if (tr$n_introductions_range["min"] == tr$n_introductions_range["max"])
                  as.character(tr$n_introductions_range["min"])
                else sprintf("%d–%d", tr$n_introductions_range["min"], tr$n_introductions_range["max"]),
                if (tr$n_introductions_range["max"] == 1) "" else "s",
                paste(tr$likely_origins %||% "unclear", collapse = ", "))
      else ""
    ))
  }
  if (!is.null(tr$r_at_sampling_median))
    parts <- c(parts, sprintf(
      "Estimated reproduction number at sampling is %.2f%s.",
      tr$r_at_sampling_median,
      if (!is.null(tr$r_at_sampling_range))
        sprintf(" (range %.2f–%.2f)", tr$r_at_sampling_range["min"], tr$r_at_sampling_range["max"])
      else ""
    ))
  if (!is.null(out$scenarios$baseline$median))
    parts <- c(parts, sprintf(
      "Under the baseline scenario the model projects ~%s weekly cases by week %d (95%% PI %s–%s).",
      format(round(as.numeric(out$scenarios$baseline$median)), big.mark = ","),
      as.integer(out$scenarios$baseline$week_ahead %||% 8),
      format(round(as.numeric(out$scenarios$baseline$lower95 %||% 0)), big.mark = ","),
      format(round(as.numeric(out$scenarios$baseline$upper95 %||% 0)), big.mark = ",")
    ))
  strain <- if (length(f$distribution$strain_profiles %||% list()))
    f$distribution$strain_profiles[[1]] else NULL
  if (!is.null(strain))
    parts <- c(parts, sprintf(
      "The closest historical strain (%s) was active %s–%s across %s countr%s%s.",
      strain$strain, strain$first_seen %||% "?", strain$last_seen %||% "?",
      strain$n_countries %||% "an unknown number of",
      if ((strain$n_countries %||% 2) == 1) "y" else "ies",
      if (!is.null(strain$linked_cases))
        sprintf(", linked to %s cases and %s deaths",
                format(as.numeric(strain$linked_cases), big.mark = ","),
                format(as.numeric(strain$linked_deaths %||% 0), big.mark = ","))
      else ""
    ))
  cv <- .pg_brief_caveats(f)
  if (length(cv))
    parts <- c(parts, paste0("Key caveats: ", paste(utils::head(cv, 3), collapse = " ")))
  paste(parts, collapse = " ")
}

.ollama_brief <- function(facts, timeout = 90) {
  fail <- function(e) list(text = NULL, source = "none", model = NULL, error = e)
  # Compact facts for the prompt — drop bulky/derivative structures.
  llm_facts <- facts
  llm_facts$spread_assessments <- NULL
  llm_facts$situation$sample_ids <- NULL
  facts_json <- tryCatch(
    jsonlite::toJSON(llm_facts, auto_unbox = TRUE, na = "null", pretty = FALSE,
                     digits = NA),
    error = function(e) NULL
  )
  if (is.null(facts_json)) return(fail("Could not serialize facts"))

  prompt <- paste0(
    "You are a genomic epidemiologist writing the bottom-line summary of a ",
    "Genomic Intelligence Brief for field epidemiologists deciding how to ",
    "respond to newly sequenced pathogen genomes. Below is a JSON object of ",
    "computed facts about the query genomes of one species.\n\n",
    "Write one paragraph of 4-6 plain-language sentences covering, in order: ",
    "(1) what was detected and where/when it was sampled; (2) how confident the ",
    "identification is — closest known reference, % identity, nearest historical ",
    "outbreak; (3) the spread evidence — genetic clustering, inferred origin and ",
    "introductions, R at sampling; (4) the baseline projection; and (5) the single ",
    "most important caveat. Write for someone who will act on this — concrete, ",
    "no jargon, no markdown, no bullet points.\n",
    "Rules: use ONLY the numbers given — never invent figures. If a field is ",
    "absent, do not mention it. If the historical record is sparse or queries ",
    "failed QC, say the assessment is provisional rather than stating it ",
    "confidently.\n\n",
    "FACTS:\n", facts_json
  )
  .ollama_generate(prompt, num_predict = 450, timeout = timeout)
}

# Orchestrator: LLM first, deterministic template fallback.
.pg_intelligence_brief <- function(facts) {
  if (is.null(facts) || is.null(facts$situation) || is.null(facts$situation$n_query_genomes))
    return(list(text = "No data available to compose an intelligence brief.",
                source = "none", model = NULL, error = NULL))
  r <- .ollama_brief(facts)
  if (!is.null(r$text) && nzchar(r$text)) return(r)
  list(text = .pg_brief_template(facts), source = "template",
       model = r$model, error = r$error)
}

# Full payload written to intelligence_brief.json.
pg_build_intelligence_brief <- function(species, pi_dir = NULL, tc_dir = NULL,
                                        mp_dir = NULL, tree_note_path = NULL,
                                        spread_json_path = NULL, duckdb_path = NULL,
                                        lit_dir = NULL) {
  facts <- .pg_brief_facts(species, pi_dir, tc_dir, mp_dir,
                           tree_note_path, spread_json_path, duckdb_path, lit_dir)
  narrative <- .pg_intelligence_brief(facts)
  .brief_prune(list(
    species     = tolower(species),
    generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    situation    = facts$situation,
    identification = facts$identification,
    genomics     = facts$genomics,
    transmission = facts$transmission,
    distribution = facts$distribution,
    outlook      = facts$outlook,
    historical_outbreaks = facts$historical_outbreaks,
    evidence_base = facts$evidence_base,
    spread_assessments = facts$spread_assessments,
    narrative    = list(text = narrative$text, source = narrative$source,
                        model = narrative$model %||% "", error = narrative$error %||% ""),
    caveats      = .pg_brief_caveats(facts),
    data_gaps    = facts$data_gaps
  ))
}

# ---------------------------------------------------------------------------
# Dashboard side — JSON reading + renderers
# ---------------------------------------------------------------------------

brief_json_path <- function(outdir, species) {
  file.path(outdir, "intelligence_brief", species, "intelligence_brief.json")
}

brief_read <- function(outdir, species) {
  if (is.null(species)) return(NULL)
  .brief_json(brief_json_path(outdir, species))
}

.brief_pct <- function(x, digits = 1) if (is.null(x) || is.na(x)) "—" else paste0(formatC(as.numeric(x) * 100, digits = digits, format = "f"), "%")
.brief_fmt <- function(x, digits = 2, comma = FALSE) {
  if (is.null(x) || (length(x) == 1 && is.na(x))) return("—")
  v <- as.numeric(x)
  if (comma) format(round(v), big.mark = ",") else formatC(v, digits = digits, format = "f")
}
.brief_rng_txt <- function(rng, fmt = "%.2f") {
  if (is.null(rng)) return("—")
  sprintf(paste0(fmt, "–", fmt), rng$min, rng$max)
}

# Small stat block used across the brief sections.
.brief_stat <- function(label, value, detail = NULL) {
  div(
    class = "pi-summary-card", style = "min-width: 130px; padding: 12px 14px;",
    div(class = "pi-sc-label", label),
    div(class = "pi-sc-value", style = "font-size: 1.15rem;", value),
    if (!is.null(detail)) div(class = "pi-sc-detail", detail)
  )
}

# Section block: title + stat strip + optional extra content + deep-link.
.brief_section <- function(icon_name, title, subtitle, stats, extra = NULL,
                           nav_tab = NULL) {
  link <- if (!is.null(nav_tab)) {
    tags$a(
      href = "#", style = "font-size: 0.78rem; margin-left: auto;",
      onclick = sprintf(
        "Shiny.setInputValue('overview_nav_click', '%s', {priority: 'event'}); return false;",
        nav_tab),
      "Details ", icon("arrow-right")
    )
  }
  div(
    style = "border: 1px solid #e9ecef; border-radius: 8px; padding: 14px 16px; margin-bottom: 12px; background: #fff;",
    div(
      style = "display: flex; align-items: baseline; gap: 8px; margin-bottom: 10px;",
      icon(icon_name, style = "color: #4A6C8C;"),
      span(style = "font-weight: 700; font-size: 0.95rem;", title),
      tags$small(class = "text-muted", subtitle),
      link
    ),
    if (length(stats)) div(class = "pi-summary-strip", style = "margin-bottom: 4px;", stats),
    extra
  )
}

# Compact front-page body for the Intelligence Overview card.
brief_overview_ui <- function(brief, species, outdir = NULL) {
  if (is.null(brief)) {
    return(div(
      style = "text-align: center; padding: 28px 0; color: #888;",
      icon("file-lines", class = "fa-2x"),
      p(style = "margin-top: 10px;",
        "No intelligence brief generated yet for ", strong(toupper(species)),
        ". Run the pipeline (GENERATE_INTELLIGENCE_BRIEF) or check ",
        tags$code(brief_json_path(outdir %||% "results", species)), ".")
    ))
  }
  s <- brief$situation; idn <- brief$identification
  tr <- brief$transmission; out <- brief$outlook

  # Verdict badges (dashboard-side heuristic over the identification table)
  verdict <- NULL
  if (!is.null(outdir) && exists("compute_pi_assessment")) {
    a <- tryCatch(compute_pi_assessment(species, outdir), error = function(e) NULL)
    if (!is.null(a)) verdict <- a
  }

  badge_row <- div(
    style = "display: flex; gap: 8px; flex-wrap: wrap; align-items: center; margin-bottom: 12px;",
    if (!is.null(verdict))
      span(class = paste("badge", paste0("badge-", assessment_status(verdict$assessment))),
           style = "font-size: 0.8rem; padding: 5px 10px;", verdict$assessment),
    if (!is.null(verdict))
      span(class = paste("badge", confidence_badge_class(verdict$confidence)),
           style = "font-size: 0.75rem; padding: 4px 9px;",
           paste("Confidence:", verdict$confidence)),
    if (length(tr$risk_labels))
      span(class = "badge badge-info", style = "font-size: 0.75rem; padding: 4px 9px;",
           paste("Transmission risk:", paste(tr$risk_labels, collapse = "/"))),
    span(class = "badge badge-light", style = "font-size: 0.7rem; padding: 4px 9px;",
         paste0("Narrative: ", brief$narrative$source %||% "unknown")),
    if (!is.null(brief$generated_at))
      tags$small(class = "text-muted", paste("Generated", brief$generated_at))
  )

  situation_strip <- div(
    class = "pi-summary-strip",
    .brief_stat("Query genomes", s$n_query_genomes %||% "—"),
    .brief_stat("Countries", length(s$query_countries %||% character()),
                paste(utils::head(s$query_countries %||% character(), 3), collapse = ", ")),
    .brief_stat("Collection dates",
                if (!is.null(s$collection_date_range))
                  paste0(format(as.Date(s$collection_date_range$min), "%b %d"), " – ",
                         format(as.Date(s$collection_date_range$max), "%b %d, %Y"))
                else "—"),
    .brief_stat("Median coverage", .brief_pct(s$median_coverage))
  )

  narrative_block <- div(
    style = "border-left: 4px solid #4A6C8C; background: #eef3f8; padding: 12px 16px; margin-bottom: 14px; font-size: 0.92rem; line-height: 1.5;",
    brief$narrative$text
  )

  id_stats <- list(
    .brief_stat("Identified as", toupper(s$species_identified %||% "—")),
    .brief_stat("Closest reference", idn$closest_reference %||% "—",
                paste0(.brief_pct(idn$pct_identity_median), " identity")),
    .brief_stat("Genetic distance", .brief_fmt(idn$genetic_distance_median, 4),
                paste0("max ", .brief_fmt(idn$genetic_distance_max, 4))),
    .brief_stat("Outbreak match", paste(idn$closest_outbreak_match %||% "—", collapse = ", "))
  )

  gn <- brief$genomics
  gen_stats <- list(
    .brief_stat("Query mutation burden", .brief_fmt(gn$query_burden_median, 0),
                paste0("background ", .brief_fmt(gn$background_burden_median, 0))),
    .brief_stat("Burden p-value", .brief_fmt(gn$burden_test_p, 3)),
    .brief_stat("Phenotype-linked muts", gn$n_phenotype_mutations %||% 0),
    .brief_stat("Top discriminating proteins",
                paste(gn$top_vip_proteins %||% "—", collapse = ", "))
  )

  tr_stats <- list(
    .brief_stat("Genetic clusters", tr$n_query_clusters %||% "—",
                paste0("largest ", tr$largest_cluster_size %||% "—", " genomes")),
    .brief_stat("Likely origin", paste(tr$likely_origins %||% "—", collapse = ", ")),
    .brief_stat("Introductions", .brief_rng_txt(tr$n_introductions_range, "%d")),
    .brief_stat("R at sampling", .brief_fmt(tr$r_at_sampling_median, 2),
                if (!is.null(tr$r_at_sampling_range))
                  .brief_rng_txt(tr$r_at_sampling_range, "%.2f"))
  )

  baseline <- out$scenarios$baseline
  strain1 <- if (length(brief$distribution$strain_profiles %||% list()))
    brief$distribution$strain_profiles[[1]] else NULL
  out_stats <- list(
    .brief_stat("Baseline projection", 
                if (!is.null(baseline$median)) paste0("~", .brief_fmt(baseline$median, 0, comma = TRUE), " wkly cases") else "—",
                if (!is.null(baseline$week_ahead)) paste0("week ", baseline$week_ahead) else NULL),
    .brief_stat("Closest strain history",
                if (!is.null(strain1)) paste0(strain1$n_countries %||% "—", " countries") else "—",
                if (!is.null(strain1)) paste0(.brief_fmt(strain1$linked_cases, 0, comma = TRUE), " cases, ",
                                              .brief_fmt(strain1$linked_deaths, 0, comma = TRUE), " deaths") else NULL),
    .brief_stat("Historical outbreaks", length(brief$historical_outbreaks %||% list()))
  )

  caveats_block <- if (length(brief$caveats)) {
    div(
      class = "alert alert-warning", style = "font-size: 0.8rem; padding: 8px 12px; margin-bottom: 4px;",
      icon("triangle-exclamation"), strong(" Caveats: "),
      paste(brief$caveats, collapse = " ")
    )
  }

  tagList(
    badge_row,
    situation_strip,
    narrative_block,
    .brief_section("dna", "What is it", "— identification",
                   id_stats, nav_tab = paste0("pi_", species)),
    .brief_section("microscope", "What's different about it", "— genomics",
                   gen_stats, nav_tab = paste0("pg_", species)),
    .brief_section("share-nodes", "How it's spreading", "— transmission evidence",
                   tr_stats, nav_tab = paste0("ts_", species)),
    .brief_section("earth-africa", "Where it could go", "— distribution & outlook",
                   out_stats, nav_tab = paste0("gt_", species)),
    caveats_block
  )
}

# Full-page brief for the dedicated Intelligence Brief tab.
brief_full_ui <- function(brief, species, outdir = NULL) {
  if (is.null(brief)) {
    return(bs4Dash::bs4Card(
      title = "Intelligence Brief", width = 12, status = "secondary",
      div(style = "text-align: center; padding: 40px 0; color: #888;",
          icon("file-lines", class = "fa-3x"),
          h4("No intelligence brief", style = "margin-top: 16px;"),
          p("No intelligence_brief.json found for this species. Run the pipeline with the brief stage enabled."))
    ))
  }
  s <- brief$situation; tr <- brief$transmission
  strain_tbl <- NULL
  sp_list <- brief$distribution$strain_profiles %||% list()
  if (length(sp_list)) {
    strain_tbl <- tags$table(
      class = "table table-sm table-striped", style = "font-size: 0.82rem;",
      tags$thead(tags$tr(lapply(c("Strain", "First seen", "Last seen", "Countries",
                                  "Origin", "Cases", "Deaths", "CFR", "Behaviour"),
                                tags$th))),
      tags$tbody(lapply(sp_list, function(sp) {
        tags$tr(
          tags$td(sp$strain %||% "—"), tags$td(sp$first_seen %||% "—"),
          tags$td(sp$last_seen %||% "—"), tags$td(sp$n_countries %||% "—"),
          tags$td(sp$origin_country %||% "—"),
          tags$td(.brief_fmt(sp$linked_cases, 0, comma = TRUE)),
          tags$td(.brief_fmt(sp$linked_deaths, 0, comma = TRUE)),
          tags$td(.brief_pct(sp$mean_cfr)),
          tags$td(gsub("_", " ", sp$behavior_label %||% "—"))
        )
      }))
    )
  }

  proj_rows <- brief$outlook$scenarios %||% list()
  proj_tbl <- if (length(proj_rows)) {
    tags$table(
      class = "table table-sm", style = "font-size: 0.82rem; max-width: 620px;",
      tags$thead(tags$tr(lapply(c("Scenario", "Week", "Median weekly cases", "95% PI", "R used"), tags$th))),
      tags$tbody(lapply(names(proj_rows), function(sc) {
        d <- proj_rows[[sc]]
        tags$tr(
          tags$td(tools::toTitleCase(sc)), tags$td(d$week_ahead %||% "—"),
          tags$td(.brief_fmt(d$median, 0, comma = TRUE)),
          tags$td(paste0(.brief_fmt(d$lower95, 0, comma = TRUE), " – ",
                         .brief_fmt(d$upper95, 0, comma = TRUE))),
          tags$td(.brief_fmt(d$r_used, 2))
        )
      }))
    )
  }

  sa <- brief$spread_assessments %||% list()
  per_query_block <- if (length(sa)) {
    tagList(
      h5(icon("viruses"), " Per-query transmission assessments"),
      lapply(names(sa), function(q) {
        a <- sa[[q]]
        div(style = "border: 1px solid #e9ecef; border-radius: 6px; padding: 10px 14px; margin-bottom: 8px;",
            div(style = "display:flex; gap:8px; align-items:center; margin-bottom:4px;",
                strong(q),
                span(class = "badge badge-light", style = "font-size:0.65rem;", a$source %||% "")),
            p(style = "font-size: 0.82rem; margin: 0; white-space: pre-wrap;", a$text))
      })
    )
  }

  hist_rows <- brief$historical_outbreaks %||% list()
  hist_tbl <- if (length(hist_rows)) {
    tags$table(
      class = "table table-sm table-striped", style = "font-size: 0.82rem; max-width: 720px;",
      tags$thead(tags$tr(lapply(c("Year", "Country", "Admin1", "Cases", "Deaths", "CFR"), tags$th))),
      tags$tbody(lapply(hist_rows, function(o) {
        tags$tr(tags$td(o$start_year %||% "—"), tags$td(o$country %||% "—"),
                tags$td(o$admin1 %||% "—"), tags$td(.brief_fmt(o$cases, 0, comma = TRUE)),
                tags$td(.brief_fmt(o$deaths, 0, comma = TRUE)), tags$td(.brief_pct(o$cfr)))
      }))
    )
  }

  eb <- brief$evidence_base
  evidence_line <- if (!is.null(eb) && !is.null(eb$total_papers)) {
    doms <- names(eb$domains %||% list())
    p(class = "text-muted", style = "font-size: 0.8rem;",
      sprintf("Evidence base: %s papers screened across %d literature domains (%s).",
              format(eb$total_papers, big.mark = ","), length(doms),
              paste(doms, collapse = ", ")))
  }

  tagList(
    bs4Dash::bs4Card(
      title = tagList(icon("file-lines"), sprintf(" Intelligence Brief — %s", toupper(species))),
      width = 12, status = "primary", solidHeader = TRUE,
      brief_overview_ui(brief, species, outdir)
    ),
    bs4Dash::bs4Card(
      title = tagList(icon("earth-africa"), " Lineage distribution & historical burden"),
      width = 12, status = "secondary",
      p(class = "text-muted", style = "font-size:0.85rem;",
        "Strain-level history for the strains closest to the query genomes, plus the historical outbreak record."),
      strain_tbl,
      if (!is.null(hist_tbl)) tagList(h6("Historical outbreak record"), hist_tbl)
    ),
    bs4Dash::bs4Card(
      title = tagList(icon("chart-area"), " Scenario outlook"),
      width = 12, status = "secondary",
      div(class = "alert alert-warning", style = "font-size:0.82rem; padding:8px 12px;",
          icon("triangle-exclamation"), strong(" Model-derived. "),
          "Branching-process projections — estimates, not validated forecasts."),
      proj_tbl,
      per_query_block,
      evidence_line
    )
  )
}

# Plain-text serialization used to pre-fill the preview/edit modal.
brief_markdown <- function(brief, species) {
  if (is.null(brief)) return("No intelligence brief content available.")
  s <- brief$situation; idn <- brief$identification; tr <- brief$transmission
  out <- brief$outlook
  lines <- c(
    sprintf("GENOMIC INTELLIGENCE BRIEF — %s", toupper(species)),
    sprintf("Generated: %s", brief$generated_at %||% "unknown"),
    "",
    "BOTTOM LINE",
    brief$narrative$text %||% "",
    "",
    "SITUATION",
    sprintf("  Query genomes: %s (%s)", s$n_query_genomes %||% "—",
            paste(s$sample_ids %||% character(), collapse = ", ")),
    sprintf("  Identified as: %s", s$species_identified %||% "—"),
    sprintf("  Countries: %s", paste(s$query_countries %||% "—", collapse = ", ")),
    sprintf("  Collection dates: %s – %s",
            s$collection_date_range$min %||% "—", s$collection_date_range$max %||% "—"),
    "",
    "IDENTIFICATION",
    sprintf("  Closest reference: %s", idn$closest_reference %||% "—"),
    sprintf("  Median identity to closest ref: %s", .brief_pct(idn$pct_identity_median)),
    sprintf("  Median genetic distance: %s", .brief_fmt(idn$genetic_distance_median, 4)),
    sprintf("  Closest outbreak match: %s", paste(idn$closest_outbreak_match %||% "—", collapse = ", ")),
    "",
    "TRANSMISSION & SPREAD",
    sprintf("  Genetic clusters: %s (largest %s genomes)", tr$n_query_clusters %||% "—",
            tr$largest_cluster_size %||% "—"),
    sprintf("  Likely origin: %s", paste(tr$likely_origins %||% "—", collapse = ", ")),
    sprintf("  Introductions into query country: %s", .brief_rng_txt(tr$n_introductions_range, "%d")),
    sprintf("  R at sampling (median): %s", .brief_fmt(tr$r_at_sampling_median, 2)),
    ""
  )
  base <- out$scenarios$baseline
  if (!is.null(base))
    lines <- c(lines, "OUTLOOK (baseline)",
               sprintf("  ~%s weekly cases by week %s (95%% PI %s–%s)",
                       .brief_fmt(base$median, 0, comma = TRUE), base$week_ahead %||% "—",
                       .brief_fmt(base$lower95, 0, comma = TRUE), .brief_fmt(base$upper95, 0, comma = TRUE)), "")
  if (length(brief$caveats))
    lines <- c(lines, "CAVEATS", paste0("  - ", brief$caveats), "")
  if (!is.null(brief$evidence_base$total_papers))
    lines <- c(lines, sprintf("Evidence base: %s papers screened across %d domains.",
                              format(brief$evidence_base$total_papers, big.mark = ","),
                              length(brief$evidence_base$domains %||% list())))
  paste(lines, collapse = "\n")
}
