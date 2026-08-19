#!/usr/bin/env Rscript
# Thin wrapper around revtools: drop empty abstracts and deduplicate per-domain JSON files.

suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(dplyr))

# r-revtools is not on conda-forge; install revtools from CRAN the first time this runs.
if (!requireNamespace("revtools", quietly = TRUE)) {
  install.packages("revtools", repos = "https://cloud.r-project.org", dependencies = c("Depends", "Imports", "LinkingTo"))
}
suppressPackageStartupMessages(library(revtools))

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

if (is.na(input_dir) || is.na(outdir) || is.na(species) || is.na(domain)) {
  stop("Usage: deduplicate_literature.R --input-dir <dir> --outdir <dir> --species <species> --domain <domain>")
}

in_files <- list.files(input_dir, pattern = "\\.json$", full.names = TRUE)
# Avoid reprocessing a summary file if the module is ever rerun on its own output,
# and skip the placeholder written when PubMed found no publications.
in_files <- in_files[!grepl("^(deduplication_summary|no_results)\\.json$", basename(in_files))]

records <- lapply(in_files, fromJSON, simplifyVector = TRUE)
input_count <- length(records)

# Drop records with missing/empty abstracts.
has_abstract <- function(r) {
  a <- r$abstract
  !is.null(a) && length(a) > 0 && !is.na(a[1]) && nchar(trimws(as.character(a[1]))) > 0
}
records <- records[vapply(records, has_abstract, logical(1))]
no_abstract_removed <- input_count - length(records)

if (length(records) == 0) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  summary <- list(
    species = species,
    domain = domain,
    method = "revtools",
    input_count = input_count,
    no_abstract_removed = no_abstract_removed,
    duplicate_removed = 0,
    output_count = 0
  )
  write_json(summary, file.path(outdir, "deduplication_summary.json"), pretty = TRUE, auto_unbox = TRUE)
  quit(save = "no", status = 0)
}

# Build a revtools-compatible data frame.
df <- tibble(
  id = seq_along(records),
  pmid = vapply(records, function(r) as.character(r$pmid %||% ""), character(1)),
  title = vapply(records, function(r) as.character(r$title %||% ""), character(1)),
  author = vapply(records, function(r) paste(unlist(r$authors %||% character(0)), collapse = " "), character(1)),
  year = vapply(records, function(r) as.character(r$year %||% ""), character(1)),
  journal = vapply(records, function(r) as.character(r$journal %||% ""), character(1)),
  doi = vapply(records, function(r) as.character(r$doi %||% ""), character(1)),
  abstract = vapply(records, function(r) as.character(r$abstract %||% ""), character(1))
)

# Deduplicate with revtools (title-only to avoid multi-column data.frame == bug).
dup <- find_duplicates(df, match_variable = "title")
uniq <- extract_unique_references(df, dup)

if ("id" %in% names(uniq)) {
  keep_idx <- as.integer(uniq$id)
} else if ("record_id" %in% names(uniq)) {
  keep_idx <- as.integer(uniq$record_id)
} else {
  keep_idx <- as.integer(rownames(uniq))
}

# Write selected records back out.
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
for (i in keep_idx) {
  r <- records[[i]]
  if (is.character(r$authors) && !is.list(r$authors)) r$authors <- as.list(r$authors)
  if (is.character(r$keywords) && !is.list(r$keywords)) r$keywords <- as.list(r$keywords)
  out_path <- file.path(outdir, paste0(as.character(r$pmid %||% i), ".json"))
  write_json(r, out_path, pretty = TRUE, auto_unbox = TRUE)
}

summary <- list(
  species = species,
  domain = domain,
  method = "revtools",
  input_count = input_count,
  no_abstract_removed = no_abstract_removed,
  duplicate_removed = length(records) - length(keep_idx),
  output_count = length(keep_idx)
)
write_json(summary, file.path(outdir, "deduplication_summary.json"), pretty = TRUE, auto_unbox = TRUE)
