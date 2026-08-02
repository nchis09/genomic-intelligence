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

epi_dir <- file.path(outdir, "epi_data")
dir.create(epi_dir, showWarnings = FALSE, recursive = TRUE)

if (length(datasets) == 0) {
  warning("No HDX datasets found for term: ", disease)
  search_summary <- tibble(
    rank = integer(0),
    title = character(0),
    name = character(0),
    organization = character(0),
    source = character(0),
    species_match = character(0),
    downloaded = logical(0)
  )
  write.table(search_summary,
              file = file.path(outdir, "rhdx_search_results.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
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

species_list_for <- function(ds) {
  if (length(dataset_to_species) == 0) return(NA_character_)
  mapped <- dataset_to_species[[get_dataset_name(ds)]]
  if (is.null(mapped)) return(NA_character_)
  paste(mapped, collapse = ",")
}

search_summary <- tibble(
  rank = seq_along(datasets),
  title = vapply(datasets, get_dataset_title, character(1)),
  name = vapply(datasets, get_dataset_name, character(1)),
  organization = vapply(datasets, get_dataset_org, character(1)),
  source = rep("HDX", length(datasets)),
  species_match = vapply(datasets, species_list_for, character(1)),
  downloaded = rep(FALSE, length(datasets))
)

# ---- Helper: find best resource for a dataset ----
format_score <- function(fmt) {
  switch(fmt,
         csv = 1L,
         xlsx = 2L,
         xls = 3L,
         json = 4L,
         txt = 5L,
         6L)
}

find_best_resource <- function(ds) {
  formats <- c("csv", "xlsx", "xls", "json", "txt")
  resources <- tryCatch(get_resources(ds, format = formats), error = function(e) list())
  if (length(resources) == 0) return(NULL)

  fmt_scores <- vapply(resources, function(r) format_score(tolower(r$get_format())), integer(1))
  resources <- resources[order(fmt_scores)]

  rs <- resources[[1]]
  fmt <- tolower(rs$get_format())
  list(ds = ds, rs = rs, format = fmt)
}

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

# ---- Download all mapped datasets ----
raw_dir <- file.path(outdir, "raw_downloads")
dir.create(raw_dir, showWarnings = FALSE, recursive = TRUE)

mapped_names <- names(dataset_to_species)
if (!is.null(species) && species != "all" && length(dataset_to_species) > 0) {
  mapped_names <- mapped_names[vapply(mapped_names, function(n) species %in% dataset_to_species[[n]], logical(1))]
}

downloaded <- 0
for (ds_name in mapped_names) {
  idx <- which(search_summary$name == ds_name)
  if (length(idx) == 0) {
    message(paste("Mapped dataset not found in HDX search results:", ds_name))
    next
  }

  ds <- datasets[[idx[1]]]
  hit <- find_best_resource(ds)
  if (is.null(hit)) {
    message(paste("No usable resource for dataset:", ds_name))
    next
  }

  message(paste("Downloading:", ds_name, "-", hit$rs$data$name, "(", hit$format, ")"))

  df <- tryCatch({
    download_and_read(hit$rs, hit$format, raw_dir)
  }, error = function(e) {
    message(paste("ERROR reading resource for", ds_name, ":", e$message))
    NULL
  })

  if (is.null(df) || nrow(df) == 0) {
    message(paste("Skipping empty resource for dataset:", ds_name))
    next
  }

  safe_name <- gsub("[^A-Za-z0-9_-]", "_", ds_name)
  out_csv <- file.path(epi_dir, paste0(safe_name, ".csv"))
  write_csv(df, out_csv)
  message(paste("Wrote:", out_csv, "(", nrow(df), "rows x", ncol(df), "cols)"))

  search_summary$downloaded[idx[1]] <- TRUE
  downloaded <- downloaded + 1
}

write.table(search_summary,
            file = file.path(outdir, "rhdx_search_results.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

if (downloaded == 0) {
  warning(paste("No mapped datasets could be downloaded for term:", disease))
  quit(save = "no", status = 0)
}

message(paste("Downloaded", downloaded, "dataset(s) to", epi_dir))
