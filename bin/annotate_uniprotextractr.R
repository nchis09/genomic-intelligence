#!/usr/bin/env Rscript
#
# annotate_uniprotextractr.R
#
# Wrapper script that sources the cloned UniProtExtractR repo function
# and parses the UniProtKB TSV download (produced by extract_query_proteins.py)
# to extract clean fields for 9 categories:
#   DNA binding, Pathway, Transmembrane, Signal peptide, Protein families,
#   Domain [FT], Motif, Involvement in disease, Subcellular location [CC]
#
# Usage:
#   Rscript annotate_uniprotextractr.R \
#     --input query_uniprot_download.tsv \
#     --extractr_dir tools/UniProtExtractR \
#     --outdir uniprotextractr_results \
#     --prefix bdbv
#

suppressPackageStartupMessages({
  library(stringr)
  library(stringi)
  library(tibble)
})

# ---- Parse arguments ----
args <- commandArgs(trailingOnly = TRUE)

parse_arg <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) > 0 && idx < length(args)) return(args[idx + 1])
  return(default)
}

input_tsv    <- parse_arg("--input")
extractr_dir <- parse_arg("--extractr_dir", "tools/UniProtExtractR")
outdir       <- parse_arg("--outdir", "uniprotextractr_results")
prefix       <- parse_arg("--prefix", "query")
map_file     <- parse_arg("--mapping", NULL)

if (is.null(input_tsv)) {
  stop("Usage: Rscript annotate_uniprotextractr.R --input <tsv> [--extractr_dir <dir>] [--outdir <dir>] [--prefix <str>] [--mapping <tsv>]")
}

# ---- Source UniProtExtractR function from cloned repo ----
extractr_script <- file.path(extractr_dir, "R", "UniProtExtractR_uniprotextract.R")
if (!file.exists(extractr_script)) {
  stop(paste("UniProtExtractR script not found at:", extractr_script))
}
source(extractr_script)
message("Loaded UniProtExtractR::uniprotextract function")

# ---- Read UniProtKB TSV ----
message(paste("Reading UniProtKB TSV:", input_tsv))
up_df <- read.table(input_tsv, header = TRUE, sep = "\t",
                     quote = "", fill = TRUE, comment.char = "",
                     stringsAsFactors = FALSE, na.strings = c("", "NA", "NULL"))

message(paste("  Loaded", nrow(up_df), "protein entries with", ncol(up_df), "columns"))
message(paste("  Columns:", paste(colnames(up_df), collapse = ", ")))

if (nrow(up_df) == 0) {
  message("No data to process. Writing empty output.")
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  write.csv(data.frame(), file.path(outdir, paste0(prefix, "_uniprotextractr.csv")), row.names = FALSE)
  quit(save = "no", status = 0)
}

# ---- Load optional mapping file ----
map_up <- NULL
if (!is.null(map_file) && file.exists(map_file)) {
  map_up <- read.table(map_file, header = TRUE, sep = "\t",
                        quote = "", fill = TRUE, stringsAsFactors = FALSE)
  message(paste("  Loaded mapping file with", nrow(map_up), "entries"))
}

# ---- Run UniProtExtractR ----
message("\nRunning uniprotextract()...")
result <- tryCatch({
  uniprotextract(my.uniprot.df = up_df, map.up = map_up, write.local = FALSE)
}, error = function(e) {
  message(paste("ERROR in uniprotextract:", e$message))
  return(up_df)  # return original on failure
})

message(paste("  Result:", nrow(result), "rows x", ncol(result), "columns"))

# ---- Show which new columns were added ----
new_cols <- setdiff(colnames(result), colnames(up_df))
if (length(new_cols) > 0) {
  message(paste("  New extracted columns:", paste(new_cols, collapse = ", ")))
} else {
  message("  No new columns extracted (input may lack the 9 supported categories)")
}

# ---- Write outputs ----
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Full result
full_out <- file.path(outdir, paste0(prefix, "_uniprotextractr.csv"))
write.csv(result, full_out, row.names = FALSE)
message(paste("Wrote full result to:", full_out))

# Also write a clean TSV with just key extracted fields
extracted_cols <- c("Entry", "Entry.Name", "Protein.names", "Gene.Names", "Organism", "Length")
# Add any new extracted columns
extracted_cols <- c(extracted_cols, new_cols)
# Keep only columns that exist
extracted_cols <- extracted_cols[extracted_cols %in% colnames(result)]

if (length(extracted_cols) > 0) {
  clean_out <- file.path(outdir, paste0(prefix, "_uniprotextractr_clean.tsv"))
  write.table(result[, extracted_cols, drop = FALSE], clean_out,
              sep = "\t", row.names = FALSE, quote = FALSE)
  message(paste("Wrote clean extracted fields to:", clean_out))
}

message("\nDone.")
