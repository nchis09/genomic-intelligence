#!/usr/bin/env Rscript
# summarize_geographic.R
#
# Pre-generate the geographic & temporal spread narrative for one species, so
# the dashboard renders it instantly instead of calling Ollama at view time.
#
# Reads the transmission-context TSVs produced by compute_transmission_context.py
# (+ model_transmission_potential.py), assembles species-level facts — strain
# spread histories, outbreak timeline, and where the query genomes most likely
# came from — calls the local Ollama model via dashboard/modules/llm_note.R,
# and writes geo_summary.json. Falls back to a deterministic template when
# Ollama is unreachable — the file is always produced.
#
# Usage:
#   Rscript summarize_geographic.R --tc_dir transmission_context/ebov \
#       --species ebov --outdir . --llm_note_r llm_note.R

suppressWarnings(suppressMessages({
  library(readr)
  library(dplyr)
}))

if (requireNamespace("optparse", quietly = TRUE)) {
  option_list <- list(
    optparse::make_option("--tc_dir", type = "character", help = "transmission_context/<species> directory"),
    optparse::make_option("--species", type = "character", help = "Species key (e.g. ebov, bdbv, sudv)"),
    optparse::make_option("--outdir", type = "character", default = ".", help = "Output directory"),
    optparse::make_option("--llm_note_r", type = "character", default = "llm_note.R",
              help = "Path to dashboard/modules/llm_note.R (staged by Nextflow)")
  )
  opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))
} else {
  args <- commandArgs(trailingOnly = TRUE)
  opt <- list(tc_dir = NULL, species = NULL, outdir = ".", llm_note_r = "llm_note.R")
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

species <- tolower(opt$species)
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
out_json <- file.path(opt$outdir, "geo_summary.json")

write_summary <- function(text, source = "none", model = NULL, error = NULL) {
  payload <- list(species = species, text = text, source = source,
                  model = model, error = error,
                  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
  writeLines(jsonlite::toJSON(payload, auto_unbox = TRUE, pretty = TRUE), out_json)
  message("[geo_summary] wrote ", out_json, " (source=", source, ")")
}

fail <- function(msg) {
  message("[geo_summary] ", msg)
  write_summary(paste0("No geographic reconstruction available. ", msg),
                source = "none", error = msg)
  quit(save = "no", status = 0)
}

read_tc <- function(f) {
  p <- file.path(opt$tc_dir, f)
  if (!file.exists(p)) return(NULL)
  df <- tryCatch(readr::read_tsv(p, show_col_types = FALSE), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df
}

potential <- read_tc("transmission_potential.tsv")
if (is.null(potential)) fail("transmission_potential.tsv missing or empty.")

strain_p  <- read_tc("strain_transmission_profile.tsv")
outbreaks <- read_tc("transmission_outbreak_summary.tsv")
burden    <- read_tc("transmission_burden.tsv")
src       <- read_tc("query_source_inference.tsv")
score     <- read_tc("transmission_potential_score.tsv")

# -- Strain spread histories -------------------------------------------------
# Chronological country path per background strain.
bg <- potential[!is.na(potential$background_sample), , drop = FALSE]
spread <- list()
if (nrow(bg) && "background_strain" %in% names(bg)) {
  bs <- bg[!is.na(bg$background_strain) & bg$background_strain != "", , drop = FALSE]
  for (st in unique(bs$background_strain)) {
    g <- bs[bs$background_strain == st, , drop = FALSE]
    g$yr <- suppressWarnings(as.integer(substr(g$background_collection_date, 1, 4)))
    g <- g[!is.na(g$background_country), , drop = FALSE]
    g <- g[order(g$yr), ]
    path <- unique(g$background_country)
    spread[[st]] <- list(
      first_country = path[1],
      first_year = min(g$yr, na.rm = TRUE),
      countries = path,
      n_genomes = nrow(bs[bs$background_strain == st, , drop = FALSE])
    )
  }
}
# Merge precomputed strain profiles when present
strain_profiles <- list()
if (!is.null(strain_p)) {
  for (i in seq_len(nrow(strain_p))) {
    r <- strain_p[i, ]
    strain_profiles[[as.character(r$strain)]] <- list(
      linked_cases = r$linked_cases, linked_deaths = r$linked_deaths,
      n_countries = r$n_countries, active_years = r$active_years,
      max_spread_km = r$max_spread_km, behavior_label = r$behavior_label,
      mean_cfr = r$mean_cfr)
  }
}

# -- Outbreak timeline ---------------------------------------------------------
ob <- NULL
if (!is.null(outbreaks)) {
  o <- outbreaks[!is.na(outbreaks$start_year) & !is.na(outbreaks$country), , drop = FALSE]
  o <- o[order(o$start_year), ]
  ob <- lapply(seq_len(nrow(o)), function(i) list(
    country = o$country[i], year = o$start_year[i],
    cases = o$cases[i], deaths = o$deaths[i], cfr = o$cfr[i]))
}

# -- Burden peak ---------------------------------------------------------------
peak <- NULL
if (!is.null(burden) && "week_start" %in% names(burden)) {
  b <- burden %>% dplyr::mutate(month = substr(as.character(week_start), 1, 7)) %>%
    dplyr::group_by(month) %>%
    dplyr::summarise(cases = sum(as.numeric(new_cases), na.rm = TRUE), .groups = "drop")
  b <- b[!is.na(b$cases), ]
  if (nrow(b)) {
    pk <- b[which.max(b$cases), ]
    peak <- list(month = pk$month, monthly_cases = pk$cases,
                 total_cases = sum(b$cases))
  }
}

# -- Query origins -------------------------------------------------------------
queries <- sort(unique(stats::na.omit(potential$query_sample)))
q_facts <- list()
for (qs in queries) {
  q <- potential[potential$query_sample == qs, , drop = FALSE]
  if (!nrow(q)) next
  links <- q[!is.na(q$background_sample), , drop = FALSE]
  # strongest links: confirmed-tier first, then by genetic/temporal proximity
  if ("background_div_diff" %in% names(links) && any(!is.na(links$background_div_diff)))
    links <- links[order(links$background_div_diff), , drop = FALSE]
  else if ("temporal_distance_days" %in% names(links))
    links <- links[order(links$temporal_distance_days), , drop = FALSE]
  top <- head(unique(links[!duplicated(links$background_sample), ,
                           drop = FALSE]$background_sample), 3)
  top_rows <- links[links$background_sample %in% top, , drop = FALSE]
  sr <- if (!is.null(src)) src[src$query_sample == qs, , drop = FALSE] else NULL
  q_facts[[qs]] <- list(
    query_country = q$query_country[1],
    query_strain = q$query_strain[1],
    query_collection_date = q$query_collection_date[1],
    likely_origin_country = if (nrow(sr)) sr$likely_origin_country[1] else NULL,
    n_introductions = if (nrow(sr)) sr$n_introductions_into_query_country[1] else NULL,
    closest_backgrounds = lapply(seq_len(nrow(top_rows)), function(i) list(
      sample = top_rows$background_sample[i],
      strain = top_rows$background_strain[i],
      country = top_rows$background_country[i],
      date = top_rows$background_collection_date[i],
      geo_distance_km = top_rows$geo_distance_km[i],
      temporal_distance_days = top_rows$temporal_distance_days[i]))
  )
}

facts <- list(species = toupper(species),
              strain_spread = spread,
              strain_profiles = strain_profiles,
              outbreak_timeline = ob,
              burden_peak = peak,
              query_origins = q_facts)
facts <- facts[lengths(facts) > 0]

# -- Template fallback ----------------------------------------------------------
template_text <- function() {
  parts <- c()
  if (length(spread)) {
    paths <- vapply(names(spread), function(st) {
      s <- spread[[st]]
      paste0(st, " (first ", s$first_country, " ", s$first_year,
             if (length(s$countries) > 1) paste0(" -> ", paste(s$countries[-1], collapse = ", ")) else "", ")")
    }, character(1))
    parts <- c(parts, paste("Spread histories:", paste(paths, collapse = "; ")))
  }
  if (length(q_facts)) {
    orgs <- vapply(names(q_facts), function(qs) {
      q <- q_facts[[qs]]
      org <- q$likely_origin_country %||% NULL
      nb <- if (length(q$closest_backgrounds)) q$closest_backgrounds[[1]]$country else NULL
      tgt <- org %||% nb %||% "unknown"
      paste0(qs, " -> most likely from ", tgt)
    }, character(1))
    parts <- c(parts, paste("Query origins:", paste(orgs, collapse = "; ")))
  }
  if (!length(parts)) return("No spread reconstruction available for this species.")
  paste(parts, collapse = ". ")
}

# -- LLM call -------------------------------------------------------------------
prompt <- paste0(
  "You are a genomic epidemiologist writing for a public-health reader. Below ",
  "is a JSON fact sheet about the geographic spread of ", toupper(species),
  " built only from this analysis run's data: strain spread histories, ",
  "historical outbreaks, burden peaks, and where each query genome most ",
  "likely originated.\n\n",
  "Write one paragraph (4-6 sentences): how the strain(s) spread geographically ",
  "over time and from where; where each query genome most plausibly originated ",
  "(name the closest background country/strain and the confidence implied by ",
  "shared strain vs broader proximity); and one caveat that these are inferred ",
  "links, not proven transmission. Use only the numbers and names given — no ",
  "outside knowledge, no markdown.\n\nFACTS:\n",
  jsonlite::toJSON(facts, auto_unbox = TRUE, na = "null")
)

res <- tryCatch(.ollama_generate(prompt, num_predict = 700, timeout = 90),
                error = function(e) list(text = NULL, source = "none",
                                         model = NULL, error = e$message))
if (is.null(res$text) || !nzchar(res$text)) {
  write_summary(template_text(), source = "template",
                model = NULL, error = res$error)
} else {
  write_summary(res$text, source = "ollama", model = res$model, error = NULL)
}
