#!/usr/bin/env Rscript
# generate_spread_assessment.R
#
# Pre-generate the per-query "Transmission & Spread Assessment" narrative for
# one species, so the dashboard renders it instantly instead of calling Ollama
# at view time.
#
# Reads the transmission-context TSVs produced by compute_transmission_context.py
# (+ model_transmission_potential.py), assembles a structured fact list per
# query genome, calls the local Ollama model via dashboard/modules/llm_note.R,
# and writes spread_assessment.json keyed by query_sample. If Ollama is
# unreachable the deterministic template text is written instead — the file is
# always produced.
#
# Usage:
#   Rscript generate_spread_assessment.R --tc_dir transmission_context/ebov \
#       --species ebov --outdir . \
#       --llm_note_r llm_note.R --module_r pathogen_transmission.R

suppressWarnings(suppressMessages({
  library(readr)
  library(dplyr)
}))

if (requireNamespace("optparse", quietly = TRUE)) {
  option_list <- list(
    optparse::make_option("--tc_dir",      type = "character", help = "transmission_context/<species> directory"),
    optparse::make_option("--species",     type = "character", help = "Species key (e.g. ebov, bdbv, sudv)"),
    optparse::make_option("--outdir",      type = "character", default = ".", help = "Output directory"),
    optparse::make_option("--llm_note_r",  type = "character", default = "llm_note.R",
                help = "Path to dashboard/modules/llm_note.R (staged by Nextflow)"),
    optparse::make_option("--module_r",    type = "character", default = "pathogen_transmission.R",
                help = "Path to dashboard/modules/pathogen_transmission.R (staged; provides .ts_strain_epi/.ts_fit_logistic)")
  )
  opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))
} else {
  # Minimal fallback: parse "--key value" pairs (optparse not installed).
  args <- commandArgs(trailingOnly = TRUE)
  opt <- list(tc_dir = NULL, species = NULL, outdir = ".",
              llm_note_r = "llm_note.R", module_r = "pathogen_transmission.R")
  i <- 1
  while (i <= length(args)) {
    if (startsWith(args[i], "--") && i < length(args)) {
      opt[[sub("^--", "", args[i])]] <- args[i + 1]
      i <- i + 2
    } else i <- i + 1
  }
}

if (is.null(opt$tc_dir) || is.null(opt$species))
  stop("--tc_dir and --species are required")
if (!file.exists(opt$llm_note_r))
  stop("llm_note.R not found at: ", opt$llm_note_r)
source(opt$llm_note_r)
if (file.exists(opt$module_r)) source(opt$module_r)   # growth-fit helpers

species <- tolower(opt$species)
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
out_json <- file.path(opt$outdir, "spread_assessment.json")

write_assessments <- function(assessments) {
  payload <- list(
    species = species,
    generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    assessments = assessments
  )
  writeLines(jsonlite::toJSON(payload, auto_unbox = TRUE, pretty = TRUE), out_json)
  message("[spread_assessment] wrote ", out_json, " (", length(assessments), " queries)")
}

fail <- function(msg) {
  message("[spread_assessment] ", msg)
  write_assessments(list())
  quit(save = "no", status = 0)
}

# -- Load transmission-context TSVs -------------------------------------------
read_tc <- function(f) {
  p <- file.path(opt$tc_dir, f)
  if (!file.exists(p)) return(NULL)
  df <- tryCatch(readr::read_tsv(p, show_col_types = FALSE), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df
}

potential <- read_tc("transmission_potential.tsv")
if (is.null(potential)) fail("transmission_potential.tsv missing or empty.")

score    <- read_tc("transmission_potential_score.tsv")
strain_p <- read_tc("strain_transmission_profile.tsv")
burden   <- read_tc("transmission_burden.tsv")
qctx     <- read_tc("query_epi_context.tsv")
cprof    <- read_tc("query_cluster_profile.tsv")
proj     <- read_tc("query_projection.tsv")
src      <- read_tc("query_source_inference.tsv")

if (!is.null(burden) && "week_start" %in% names(burden))
  burden$week_start <- as.Date(burden$week_start)

row1 <- function(df, col, val) {
  if (is.null(df) || !(col %in% names(df))) return(NULL)
  d <- df[df[[col]] == val, , drop = FALSE]
  if (nrow(d)) d[1, ] else NULL
}
gv <- function(row, col) {           # scalar-or-NULL getter
  if (is.null(row) || !(col %in% names(row))) return(NULL)
  v <- row[[col]][1]
  if (length(v) == 0 || (length(v) == 1 && is.na(v))) NULL else v
}

# Growth-model fit for a strain's historical outbreak (same as dashboard fit_r).
growth_fit <- function(strain) {
  if (!exists(".ts_strain_epi") || !exists(".ts_fit_logistic") || is.null(burden)) return(NULL)
  epi <- .ts_strain_epi(burden, strain)
  if (is.null(epi)) epi <- .ts_strain_epi(burden, NULL)
  tryCatch(.ts_fit_logistic(epi), error = function(e) NULL)
}

# -- Per-query facts + assessment ----------------------------------------------
queries <- sort(unique(stats::na.omit(potential$query_sample)))
if (!length(queries)) fail("No query samples in transmission_potential.tsv.")

assessments <- list()
for (qs in queries) {
  q   <- row1(potential, "query_sample", qs)
  sc  <- row1(score, "query_sample", qs)
  qc  <- row1(qctx, "query_sample", qs)
  cp  <- row1(cprof, "query_sample", qs)
  sr  <- row1(src, "query_sample", qs)
  st  <- gv(q, "query_strain")
  sp  <- if (!is.null(st) && !is.null(strain_p)) row1(strain_p, "strain", st) else NULL
  pj  <- if (!is.null(proj)) proj[proj$query_sample == qs & proj$scenario == "baseline", , drop = FALSE] else NULL
  pj8 <- if (!is.null(pj) && nrow(pj)) pj[which.max(pj$week_ahead), , drop = FALSE] else NULL
  fit <- growth_fit(st)

  facts <- list(
    query_sample           = qs,
    query_strain           = st,
    query_clade            = gv(q, "query_clade"),
    query_lineage          = gv(q, "query_lineage"),
    query_outbreak         = gv(q, "query_outbreak"),
    query_country          = gv(q, "query_country"),
    query_admin1           = gv(q, "query_admin1"),
    query_collection_date  = gv(q, "query_collection_date"),
    # historical strain profile
    strain_linked_cases    = gv(sp, "linked_cases"),
    strain_linked_deaths   = gv(sp, "linked_deaths"),
    strain_n_countries     = gv(sp, "n_countries"),
    strain_active_years    = gv(sp, "active_years"),
    strain_behavior_label  = gv(sp, "behavior_label"),
    strain_max_spread_km   = gv(sp, "max_spread_km"),
    strain_mean_cfr        = gv(sp, "mean_cfr"),
    # epi context at sampling
    r_at_sampling          = gv(qc, "r_at_sampling"),
    r_lower                = gv(qc, "r_lower"),
    r_upper                = gv(qc, "r_upper"),
    growth_phase           = gv(qc, "growth_phase"),
    sampling_fraction      = gv(qc, "sampling_fraction"),
    weekly_cases_at_sampling = gv(qc, "weekly_cases_at_sampling"),
    # genetic cluster
    cluster_id             = gv(cp, "cluster_id"),
    cluster_size           = gv(cp, "cluster_size"),
    cluster_countries      = gv(cp, "countries"),
    cluster_linked_cases   = gv(cp, "linked_cases"),
    # inferred origin
    likely_origin_country  = gv(sr, "likely_origin_country"),
    n_introductions        = gv(sr, "n_introductions_into_query_country"),
    # transmission-potential score
    risk_label             = gv(sc, "risk_label"),
    predicted_behavior     = gv(sc, "predicted_behavior"),
    n_neighbors            = gv(sc, "n_neighbors"),
    n_neighbors_with_epi   = gv(sc, "n_neighbors_with_epi"),
    # branching-process projection (baseline, final horizon week)
    proj_country           = gv(pj8, "country"),
    proj_r_used            = gv(pj8, "r_used"),
    proj_week8_median      = gv(pj8, "proj_median"),
    proj_week8_lower95     = gv(pj8, "proj_lower95"),
    proj_week8_upper95     = gv(pj8, "proj_upper95"),
    # growth-model fit
    fit_r_week             = if (!is.null(fit)) fit$r_week else NULL,
    fit_doubling_days      = if (!is.null(fit)) fit$doubling_days else NULL,
    fit_K                  = if (!is.null(fit)) fit$K else NULL,
    fit_method             = if (!is.null(fit)) fit$method else NULL
  )
  facts <- facts[!vapply(facts, is.null, logical(1))]

  res <- tryCatch(.pg_spread_assessment(facts), error = function(e)
    list(text = .pg_spread_facts_template(facts), source = "template",
         model = NULL, error = e$message))
  assessments[[qs]] <- list(text = res$text, source = res$source,
                            model = if (is.null(res$model)) "" else as.character(res$model),
                            error = if (is.null(res$error)) "" else as.character(res$error))
  message("[spread_assessment] ", qs, " -> ", res$source)
}

write_assessments(assessments)
