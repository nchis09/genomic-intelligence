#!/usr/bin/env Rscript
# Thin wrapper around the asreview CLI.
# Converts deduplicated paper JSONs into an ASReview-compatible CSV, runs
# asreview simulate, and exports the top-ranked records.

suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(yaml))

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

parse_arg <- function(flag, argv) {
  idx <- which(argv == flag)
  if (length(idx) == 0 || idx >= length(argv)) return(NA_character_)
  argv[idx + 1]
}

argv <- commandArgs(trailingOnly = TRUE)

input_dir <- parse_arg("--input-dir", argv)
outdir <- parse_arg("--outdir", argv)
species <- parse_arg("--species", argv)
domain <- parse_arg("--domain", argv)
terms_file <- parse_arg("--terms-file", argv)
min_year <- parse_arg("--min-year", argv)
n_prior_included <- as.integer(parse_arg("--n-prior-included", argv) %||% "5")
n_prior_excluded <- as.integer(parse_arg("--n-prior-excluded", argv) %||% "5")
n_stop <- as.integer(parse_arg("--n-stop", argv) %||% "10")
top_n <- as.integer(parse_arg("--top-n", argv) %||% "50")

if (is.na(input_dir) || is.na(outdir) || is.na(species) || is.na(domain) || is.na(terms_file)) {
  stop("Usage: screen_literature.R --input-dir <dir> --outdir <dir> --species <species> --domain <domain> --terms-file <yml> [options]")
}

current_year <- as.integer(format(Sys.Date(), "%Y"))
if (is.na(min_year) || nchar(min_year) == 0) {
  min_year <- current_year - 10
} else {
  min_year <- as.integer(min_year)
}

terms <- yaml::read_yaml(terms_file)
domain_terms <- terms$domains[[domain]]
if (is.null(domain_terms) || is.null(domain_terms$required)) {
  stop("Domain '", domain, "' not found or has no required terms in ", terms_file)
}
species_syn <- terms$species_synonyms[[species]]
if (is.null(species_syn) || is.null(species_syn$exact)) {
  stop("Species '", species, "' not found or has no exact synonyms in ", terms_file)
}
excluded_terms <- terms$excluded_terms[[species]] %||% character(0)
domain_excluded_terms <- terms$domain_excluded_terms[[domain]] %||% character(0)

clean_in <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub("[^a-z0-9 -]", " ", x)
  x <- gsub("[[:space:]]+", " ", x)
  trimws(x)
}

in_files <- list.files(input_dir, pattern = "\\.json$", full.names = TRUE)
# Ignore summary / helper files and any JSONs we may have produced previously.
in_files <- in_files[!grepl("(^screening_summary\\.json$|_summary\\.json$)", basename(in_files))]

records <- lapply(in_files, fromJSON, simplifyVector = TRUE)
input_count <- length(records)

if (input_count == 0) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  summary <- list(
    species = species,
    domain = domain,
    input_count = 0,
    year_removed = 0,
    seed_included = 0,
    seed_excluded = 0,
    output_count = 0,
    method = "asreview"
  )
  write_json(summary, file.path(outdir, "screening_summary.json"), pretty = TRUE, auto_unbox = TRUE)
  quit(save = "no", status = 0)
}

df <- tibble(
  record_id = seq_len(input_count),
  pmid = vapply(records, function(r) as.character(r$pmid %||% ""), character(1)),
  title = vapply(records, function(r) as.character(r$title %||% ""), character(1)),
  abstract = vapply(records, function(r) as.character(r$abstract %||% ""), character(1)),
  year = vapply(records, function(r) as.integer(r$year %||% NA), integer(1))
)

# Hard 10-year (or user-specified) cut-off
df <- df %>% filter(!is.na(year) & year >= min_year)
year_removed <- input_count - nrow(df)

if (nrow(df) == 0) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  summary <- list(
    species = species,
    domain = domain,
    input_count = input_count,
    year_removed = year_removed,
    seed_included = 0,
    seed_excluded = 0,
    output_count = 0,
    method = "asreview"
  )
  write_json(summary, file.path(outdir, "screening_summary.json"), pretty = TRUE, auto_unbox = TRUE)
  quit(save = "no", status = 0)
}

# Text used for keyword cross-checks
df$text <- clean_in(paste(df$title, df$abstract, sep = " "))

has_match <- function(text, terms) {
  if (length(terms) == 0) return(FALSE)
  vapply(text, function(t) any(vapply(terms, function(term) grepl(term, t, fixed = TRUE), logical(1))), logical(1))
}

req_terms <- tolower(as.character(domain_terms$required))
exact_terms <- tolower(as.character(species_syn$exact))
broad_terms <- tolower(as.character(species_syn$broad %||% character(0)))
exc_terms <- tolower(as.character(excluded_terms))
domain_exc_terms <- tolower(as.character(domain_excluded_terms))

df <- df %>%
  mutate(
    has_required = has_match(text, req_terms),
    has_exact = has_match(text, exact_terms),
    has_broad = if (length(broad_terms) > 0) has_match(text, broad_terms) else FALSE,
    has_excluded = if (length(exc_terms) > 0) has_match(text, exc_terms) else FALSE,
    has_domain_excluded = if (length(domain_exc_terms) > 0) has_match(text, domain_exc_terms) else FALSE,
    # Positive seed: matches domain term, has exact species term, no excluded species/domain terms
    included = case_when(
      has_excluded ~ 0L,
      has_domain_excluded ~ 0L,
      has_required & has_exact ~ 1L,
      TRUE ~ 0L
    )
  )

seed_included <- sum(df$included == 1)
seed_excluded <- sum(df$included == 0)

# Prepare ASReview input: it only needs title, abstract, and an included label.
screen_in <- df %>%
  select(record_id, title, abstract, included) %>%
  mutate(title = as.character(title), abstract = as.character(abstract), included = as.integer(included))

write.csv(screen_in, file.path(outdir, "screen_input.csv"), row.names = FALSE, na = "")

# Edge case: if there is no variation, just keep the keyword-matched records.
if (seed_included == 0 || seed_excluded == 0) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  keep_idx <- df$record_id[df$included == 1]
  for (i in keep_idx) {
    r <- records[[i]]
    if (is.character(r$authors) && !is.list(r$authors)) r$authors <- as.list(r$authors)
    if (is.character(r$keywords) && !is.list(r$keywords)) r$keywords <- as.list(r$keywords)
    out_name <- if (nchar(r$pmid %||% "") > 0) paste0(r$pmid, ".json") else paste0(i, ".json")
    write_json(r, file.path(outdir, out_name), pretty = TRUE, auto_unbox = TRUE)
  }
  summary <- list(
    species = species,
    domain = domain,
    input_count = input_count,
    year_removed = year_removed,
    seed_included = seed_included,
    seed_excluded = seed_excluded,
    output_count = length(keep_idx),
    method = "keyword_only",
    note = if (seed_included == 0) "No positive seed records found; keeping none." else "Only one label class found; using keyword filter."
  )
  write_json(summary, file.path(outdir, "screening_summary.json"), pretty = TRUE, auto_unbox = TRUE)
  quit(save = "no", status = 0)
}

n_prior_inc <- min(n_prior_included, seed_included)
n_prior_exc <- min(n_prior_excluded, seed_excluded)
# n_stop is the total number of label actions for asreview simulate.
# It must be at least the priors + 1 so the model runs and a last ranking
# table is generated.
min_stop <- n_prior_inc + n_prior_exc + 1
n_stop_val <- min(nrow(screen_in), max(min_stop, n_stop))

# Run asreview simulate. The model is seeded with n_prior included/excluded
# records and then performs n_stop label actions before saving the project.
sim_cmd <- c(
  "asreview",
  "simulate",
  "screen_input.csv",
  "-o", "screen.asreview",
  "--n-prior-included", as.character(n_prior_inc),
  "--n-prior-excluded", as.character(n_prior_exc),
  "--n-stop", as.character(n_stop_val),
  "--seed", "42"
)

sim_log <- file.path(outdir, "asreview_simulate.log")
status <- system2("asreview", args = sim_cmd[-1], stdout = sim_log, stderr = sim_log)

if (status != 0 || !file.exists("screen.asreview")) {
  stop("asreview simulate failed (status ", status, "). Log: ", sim_log)
}

# Export the last ranking from the asreview state file using the asreview
# Python API (no custom screening logic, just project I/O).
export_code <- "import asreview
import pandas as pd

p = asreview.Project.load('screen.asreview', project_path='.')
last = p.db.get_last_ranking_table()

# If the model never ran, fall back to the query order from the results table.
if last.empty:
    res = p.db.get_results_table().sort_values('time')
    last = pd.DataFrame({'record_id': res['record_id'], 'ranking': range(1, len(res) + 1)})

# ASReview stores record indices as 0-based, but the R side uses 1-based IDs.
last['record_id'] = last['record_id'].astype(int) + 1

p.close()
last.to_csv('screen_ranking.csv', index=False)"
writeLines(export_code, "export_ranking.py")
status2 <- system2("python", args = "export_ranking.py", stdout = "asreview_export.log", stderr = "asreview_export.log")

if (status2 != 0 || !file.exists("screen_ranking.csv")) {
  stop("Failed to export asreview ranking (status ", status2, ").")
}

ranking <- read.csv("screen_ranking.csv", stringsAsFactors = FALSE)
# Keep the top-N ranked records. The table is already sorted, but sorting by
# the ranking column is defensive.
keep_ids <- ranking %>%
  arrange(ranking) %>%
  head(top_n) %>%
  pull(record_id)

keep_ids <- keep_ids[keep_ids %in% df$record_id]

# Write the surviving original records.
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
for (i in keep_ids) {
  r <- records[[i]]
  if (is.character(r$authors) && !is.list(r$authors)) r$authors <- as.list(r$authors)
  if (is.character(r$keywords) && !is.list(r$keywords)) r$keywords <- as.list(r$keywords)
  out_name <- if (nchar(r$pmid %||% "") > 0) paste0(r$pmid, ".json") else paste0(i, ".json")
  write_json(r, file.path(outdir, out_name), pretty = TRUE, auto_unbox = TRUE)
}

summary <- list(
  species = species,
  domain = domain,
  input_count = input_count,
  year_removed = year_removed,
  seed_included = as.integer(seed_included),
  seed_excluded = as.integer(seed_excluded),
  output_count = length(keep_ids),
  method = "asreview",
  top_n = as.integer(top_n)
)
write_json(summary, file.path(outdir, "screening_summary.json"), pretty = TRUE, auto_unbox = TRUE)
