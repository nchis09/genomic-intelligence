#!/usr/bin/env Rscript
#
# fetch_rhdx.R
#
# Search HDX for the provided disease term using the rhdx R package,
# select the first downloadable CSV/XLSX/XLS/JSON resource, and save it
# as a flattened CSV plus a search-summary TSV.
#
# Usage:
#   Rscript fetch_rhdx.R --disease ebola --rhdx_dir tools/rhdx \
#       [--species bdbv] [--mapping database/hdx_ebola_datasets.yml] \
#       [--outdir .] [--rows 20]
#

# ---- Parse arguments ----
args <- commandArgs(trailingOnly = TRUE)

parse_arg <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) > 0 && idx < length(args)) return(args[idx + 1])
  return(default)
}

disease      <- parse_arg("--disease")
outdir       <- parse_arg("--outdir", ".")
rows         <- as.integer(parse_arg("--rows", "20"))
rhdx_dir     <- parse_arg("--rhdx_dir", "tools/rhdx")
species      <- parse_arg("--species", NULL)
mapping_file <- parse_arg("--mapping", NULL)

if (is.null(disease) || is.null(rhdx_dir)) {
  stop("Usage: Rscript fetch_rhdx.R --disease <term> --rhdx_dir <dir> [--species <virus>] [--mapping <yml>] [--outdir <dir>] [--rows <n>]")
}

if (!dir.exists(rhdx_dir)) {
  stop(paste("rhdx source directory not found:", rhdx_dir))
}

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---- Ensure rhdx is available ----
if (!require("rhdx", character.only = TRUE, quietly = TRUE)) {
  if (!require("remotes", character.only = TRUE, quietly = TRUE)) {
    install.packages("remotes", repos = "https://cloud.r-project.org/")
  }
  message(paste("Installing rhdx from local directory:", rhdx_dir))
  remotes::install_local(rhdx_dir, dependencies = c("Depends", "Imports", "LinkingTo"), upgrade = "never", force = FALSE)
}

suppressPackageStartupMessages({
  library(rhdx)
  library(readr)
  library(readxl)
  library(dplyr)
  library(tibble)
  library(yaml)
  library(jsonlite)
})

# ---- Load optional species -> dataset mapping ----
dataset_to_species <- list()
if (!is.null(mapping_file) && file.exists(mapping_file)) {
  dataset_to_species <- yaml.load_file(mapping_file)
}

# ---- Configure HDX ----
set_rhdx_config(hdx_site = "prod")

message(paste("Searching HDX for:", disease))

datasets <- search_datasets(disease, rows = rows)

if (length(datasets) == 0) {
  warning("No HDX datasets found for term: ", disease)
  search_summary <- tibble(
    rank = integer(0),
    title = character(0),
    name = character(0),
    organization = character(0),
    source = character(0)
  )
  write.table(search_summary,
              file = file.path(outdir, "rhdx_search_results.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  file.create(file.path(outdir, "epi_data.csv"))
  quit(save = "no", status = 0)
}

# ---- Build search summary ----
get_dataset_title <- function(ds) {
  out <- ds$data$title
  if (is.null(out) || length(out) == 0) out <- NA_character_
  as.character(out)[1]
}

get_dataset_name <- function(ds) {
  out <- ds$data$name
  if (is.null(out) || length(out) == 0) out <- NA_character_
  as.character(out)[1]
}

get_dataset_org <- function(ds) {
  org <- ds$data$organization
  if (is.null(org) || length(org) == 0) return(NA_character_)
  out <- org$title
  if (is.null(out) || length(out) == 0) out <- org$name
  if (is.null(out) || length(out) == 0) out <- NA_character_
  as.character(out)[1]
}

species_match_for <- function(ds) {
  if (is.null(species) || length(dataset_to_species) == 0) return(NA)
  mapped <- dataset_to_species[[get_dataset_name(ds)]]
  if (is.null(mapped)) return(FALSE)
  species %in% mapped
}

search_summary <- tibble(
  rank = seq_along(datasets),
  title = vapply(datasets, get_dataset_title, character(1)),
  name = vapply(datasets, get_dataset_name, character(1)),
  organization = vapply(datasets, get_dataset_org, character(1)),
  source = rep("HDX", length(datasets)),
  species_match = vapply(datasets, species_match_for, logical(1)),
  selected = rep(FALSE, length(datasets))
)

# ---- Find and download the best usable resource ----
find_usable_resource <- function(datasets) {
  formats <- c("csv", "xlsx", "xls", "json", "txt")
  candidates <- list()

  for (i in seq_along(datasets)) {
    ds <- datasets[[i]]

    # Species filter: only consider datasets mapped to the requested species
    if (!is.null(species) && length(dataset_to_species) > 0) {
      ds_name <- get_dataset_name(ds)
      mapped_species <- dataset_to_species[[ds_name]]
      if (is.null(mapped_species) || !(species %in% mapped_species)) next
    }

    resources <- tryCatch(get_resources(ds, format = formats), error = function(e) list())
    if (length(resources) == 0) next

    # Prefer CSV over XLSX over XLS over JSON over TXT
    format_score <- function(fmt) {
      switch(fmt,
             csv = 1L,
             xlsx = 2L,
             xls = 3L,
             json = 4L,
             txt = 5L,
             6L)
    }
    fmt_scores <- vapply(resources, function(r) format_score(tolower(r$get_format())), integer(1))
    resources <- resources[order(fmt_scores)]

    # Use the best tabular resource in the dataset
    rs <- resources[[1]]
    fmt <- tolower(rs$get_format())

    # Parse dataset date for ranking (prefer the most recent)
    date_val <- tryCatch({
      date_str <- ds$data$dataset_date
      if (is.null(date_str) || length(date_str) == 0) stop("no date")
      # HDX returns ranges like "[2021-01-28T00:00:00 TO 2021-03-17T23:59:59]"
      first_date <- regmatches(as.character(date_str),
                               regexpr("\\d{4}-\\d{2}-\\d{2}", as.character(date_str)))
      as.Date(first_date, format = "%Y-%m-%d")
    }, error = function(e) as.Date(NA))

    candidates[[length(candidates) + 1]] <- list(
      ds = ds,
      rs = rs,
      rank = i,
      format = fmt,
      date = date_val
    )
  }

  if (length(candidates) == 0) return(NULL)

  # Rank candidates by date (descending), then by original rank (ascending)
  ranks <- vapply(candidates, function(x) x$rank, integer(1))
  dates <- vapply(candidates, function(x) as.numeric(x$date), numeric(1))
  dates[is.na(dates)] <- 0

  ord <- order(-dates, ranks)
  candidates[[ord[1]]]
}

hit <- find_usable_resource(datasets)

if (!is.null(hit)) {
  search_summary$selected[hit$rank] <- TRUE
}

write.table(search_summary,
            file = file.path(outdir, "rhdx_search_results.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

if (is.null(hit)) {
  msg <- if (!is.null(species)) {
    paste0("No mapped CSV/XLSX/XLS/JSON/TXT resource found for species '", species, "' under term: ", disease)
  } else {
    paste0("No CSV/XLSX/XLS/JSON/TXT resource found for term: ", disease)
  }
  warning(msg)
  file.create(file.path(outdir, "epi_data.csv"))
  quit(save = "no", status = 0)
}

message(paste("Selected dataset:", hit$ds$data$title))
message(paste("Selected resource:", hit$rs$data$name, "(", hit$format, ")"))

# ---- Download/read the resource ----
# rhdx's read_resource() assumes comma-separated CSVs. Some HDX resources
# (e.g. Guinea 2021 Ebola data) are semicolon-delimited, so we download the raw
# file and detect the delimiter ourselves.
download_and_read <- function(rs, fmt, folder) {
  raw_path <- download_resource(rs, folder = folder)
  message(paste("Downloaded resource to:", raw_path))

  if (fmt %in% c("csv", "txt")) {
    first_line <- readLines(raw_path, n = 1, warn = FALSE)
    n_semicolon <- lengths(regmatches(first_line, gregexpr(";", first_line)))
    n_comma     <- lengths(regmatches(first_line, gregexpr(",", first_line)))
    if (n_semicolon > n_comma) {
      message("Reading semicolon-delimited file")
      df <- read_csv2(raw_path, show_col_types = FALSE)
    } else {
      message("Reading comma-delimited file")
      df <- read_csv(raw_path, show_col_types = FALSE)
    }
  } else if (fmt %in% c("xlsx", "xls")) {
    message("Reading Excel file")
    df <- read_excel(raw_path)
  } else if (fmt == "json") {
    message("Reading JSON file")
    df <- as_tibble(jsonlite::fromJSON(raw_path, flatten = TRUE, simplifyDataFrame = TRUE))
  } else {
    stop(paste("Unsupported resource format:", fmt))
  }
  df
}

df <- tryCatch({
  download_and_read(hit$rs, hit$format, outdir)
}, error = function(e) {
  stop(paste("Unable to read HDX resource:", e$message))
})

if (is.null(df) || nrow(df) == 0) {
  warning("Resource returned empty data frame")
  file.create(file.path(outdir, "epi_data.csv"))
  quit(save = "no", status = 0)
}

out_csv <- file.path(outdir, "epi_data.csv")
write_csv(df, out_csv)
message(paste("Wrote:", out_csv, "(", nrow(df), "rows x", ncol(df), "cols)"))
