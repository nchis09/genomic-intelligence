#!/usr/bin/env Rscript
# Remove records with no abstract and deduplicate per-paper PubMed JSON files.

suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(dplyr))

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
# Avoid reprocessing a summary file if the module is ever rerun on its own output.
in_files <- in_files[!grepl("^deduplication_summary\\.json$", basename(in_files))]

records <- lapply(in_files, function(f) fromJSON(f, simplifyVector = TRUE))

# Build a searchable metadata table while preserving the original index into `records`.
meta <- bind_rows(lapply(seq_along(records), function(i) {
  r <- records[[i]]
  tibble(
    idx = i,
    pmid = as.character(r$pmid %||% ""),
    title = as.character(r$title %||% ""),
    doi = as.character(r$doi %||% ""),
    abstract = as.character(r$abstract %||% ""),
    year = as.character(r$year %||% ""),
    journal = as.character(r$journal %||% ""),
    pmcid = as.character(r$pmcid %||% ""),
    publication_date = as.character(r$publication_date %||% ""),
    cited_by_count = as.integer(r$cited_by_count %||% 0),
    is_oa = list(r$is_oa),
    authors = list(r$authors %||% character(0)),
    keywords = list(r$keywords %||% character(0))
  )
}))

input_count <- nrow(meta)

# Drop records with missing/empty abstracts.
meta <- meta %>% filter(!is.na(abstract) & trimws(abstract) != "")
no_abstract_removed <- input_count - nrow(meta)

# --- Deduplication backends --------------------------------------------

dedup_revtools <- function(df) {
  if (!requireNamespace("revtools", quietly = TRUE)) stop("revtools not installed")
  d <- df
  d$record_id <- d$idx
  d$author <- sapply(d$authors, function(x) paste(unlist(x), collapse = " "))
  dup <- revtools::find_duplicates(d, match_type = "exact")
  uniq <- revtools::extract_unique_references(d, dup)
  as.integer(uniq$record_id)
}

dedup_metagear <- function(df) {
  if (!requireNamespace("metagear", quietly = TRUE)) stop("metagear not installed")
  # metagear has no stable deduplication API we can rely on; any error falls back to custom.
  result <- metagear::find_duplicates(df, match_on = c("title", "doi"))
  as.integer(df$idx[result$unique])
}

# Fallback: exact DOI match, exact normalized title match, then fuzzy title similarity with stringdist.
dedup_custom <- function(df) {
  if (nrow(df) == 0) return(integer(0))
  df$doi_norm <- tolower(trimws(df$doi))
  df$title_norm <- tolower(trimws(gsub("[^[:alnum:]]+", " ", df$title)))
  df$title_norm <- trimws(gsub("  +", " ", df$title_norm))

  # Prefer the richest record when deciding which duplicate to keep.
  df$score <- (
    (nchar(df$abstract) > 0) +
    (nchar(df$doi) > 0) +
    (nchar(df$year) > 0) +
    (nchar(df$journal) > 0) +
    (nchar(df$pmcid) > 0) +
    (sapply(df$authors, length) > 0)
  )

  # Deterministic tie-breaker by pmid.
  df <- df[order(df$score, decreasing = TRUE, df$pmid, decreasing = FALSE), ]

  has_stringdist <- requireNamespace("stringdist", quietly = TRUE)
  keep <- rep(TRUE, nrow(df))
  seen_doi <- character(0)
  seen_title <- character(0)

  for (i in seq_len(nrow(df))) {
    d <- df$doi_norm[i]
    t <- df$title_norm[i]

    if (nchar(d) > 0 && d %in% seen_doi) {
      keep[i] <- FALSE
      next
    }

    if (nchar(t) > 0 && t %in% seen_title) {
      keep[i] <- FALSE
      next
    }

    if (has_stringdist && nchar(t) > 0 && length(seen_title) > 0) {
      sim <- tryCatch(
        stringdist::stringsim(t, seen_title, method = "lv"),
        error = function(e) numeric(0)
      )
      if (length(sim) > 0 && any(sim >= 0.85, na.rm = TRUE)) {
        keep[i] <- FALSE
        next
      }
    }

    if (nchar(d) > 0) seen_doi <- c(seen_doi, d)
    if (nchar(t) > 0) seen_title <- c(seen_title, t)
  }

  df$idx[keep]
}

# --- Run deduplication -------------------------------------------------

keep_idx <- NULL
method <- "custom"

if (nrow(meta) > 0) {
  for (m in c("revtools", "metagear")) {
    tryCatch({
      if (m == "revtools") {
        keep_idx <- dedup_revtools(meta)
      } else {
        keep_idx <- dedup_metagear(meta)
      }
      method <- m
      break
    }, error = function(e) {
      message(m, " failed or not available: ", conditionMessage(e))
    })
  }
  if (is.null(keep_idx)) {
    keep_idx <- dedup_custom(meta)
    method <- "custom_stringdist"
  }
} else {
  keep_idx <- integer(0)
  method <- "none"
}

# --- Write outputs -----------------------------------------------------

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

out_meta <- meta[meta$idx %in% keep_idx, ]
out_meta <- out_meta[order(out_meta$idx), ]

for (i in out_meta$idx) {
  r <- records[[i]]
  if (is.character(r$authors) && !is.list(r$authors)) r$authors <- as.list(r$authors)
  if (is.character(r$keywords) && !is.list(r$keywords)) r$keywords <- as.list(r$keywords)

  out_path <- file.path(outdir, paste0(r$pmid, ".json"))
  write_json(r, out_path, pretty = TRUE, auto_unbox = TRUE)
}

summary <- list(
  species = species,
  domain = domain,
  method = method,
  input_count = input_count,
  no_abstract_removed = no_abstract_removed,
  duplicate_removed = input_count - no_abstract_removed - nrow(out_meta),
  output_count = nrow(out_meta)
)

write_json(summary, file.path(outdir, "deduplication_summary.json"), pretty = TRUE, auto_unbox = TRUE)
