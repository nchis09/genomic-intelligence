#!/usr/bin/env Rscript
#
# annotate_uniprotr.R
#
# Wrapper script that sources the cloned UniprotR repo functions
# and retrieves all available functional annotation for the
# species-specific UniProt accessions produced by extract_query_proteins.py.
#
# Usage:
#   Rscript annotate_uniprotr.R \
#     --accessions query_accessions.txt \
#     --uniprotr_dir tools/UniprotR \
#     --outdir uniprotr_results \
#     --prefix bdbv
#

suppressPackageStartupMessages({
  library(httr)
  library(curl)
  library(progress)
  library(jsonlite)
})

# ---- Parse arguments ----
args <- commandArgs(trailingOnly = TRUE)

parse_arg <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) > 0 && idx < length(args)) return(args[idx + 1])
  return(default)
}

acc_file     <- parse_arg("--accessions")
uniprotr_dir <- parse_arg("--uniprotr_dir", "tools/UniprotR")
outdir       <- parse_arg("--outdir", "uniprotr_results")
prefix       <- parse_arg("--prefix", "query")

if (is.null(acc_file)) {
  stop("Usage: Rscript annotate_uniprotr.R --accessions <file> [--uniprotr_dir <dir>] [--outdir <dir>] [--prefix <str>]")
}

# ---- Source UniprotR functions from cloned repo ----
r_files <- list.files(file.path(uniprotr_dir, "R"), pattern = "\\.R$|\\.r$", full.names = TRUE)
for (rf in r_files) {
  tryCatch(source(rf, local = TRUE), error = function(e) {
    message(paste("  Note: could not source", basename(rf), "-", e$message))
  })
}

# ---- Read accessions ----
accessions <- readLines(acc_file)
accessions <- trimws(accessions[nchar(trimws(accessions)) > 0])
message(paste("Querying UniProt for", length(accessions), "accessions:", paste(accessions, collapse = ", ")))

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---- Helper: safe call with error handling ----
safe_call <- function(func_name, acc_list, outdir_path) {
  tryCatch({
    func <- match.fun(func_name)
    result <- func(acc_list, directorypath = outdir_path)
    if (!is.null(result) && nrow(result) > 0) {
      outfile <- file.path(outdir_path, paste0(prefix, "_", func_name, ".csv"))
      write.csv(result, outfile, row.names = TRUE)
      message(paste("  ", func_name, "->", nrow(result), "rows ->", outfile))
      return(result)
    } else {
      message(paste("  ", func_name, "-> no data returned"))
      return(NULL)
    }
  }, error = function(e) {
    message(paste("  ", func_name, "-> ERROR:", e$message))
    return(NULL)
  })
}

# ---- Call UniprotR retrieval functions ----
message("\n=== Retrieving protein function ===")
func_result <- safe_call("GetProteinFunction", accessions, outdir)

message("\n=== Retrieving GO terms ===")
go_result <- safe_call("GetProteinGOInfo", accessions, outdir)

message("\n=== Retrieving subcellular location ===")
loc_result <- safe_call("GetSubcellular_location", accessions, outdir)

message("\n=== Retrieving pathology/biotech ===")
path_result <- safe_call("GetPathology_Biotech", accessions, outdir)

message("\n=== Retrieving family & domains ===")
fam_result <- safe_call("GetFamily_Domains", accessions, outdir)

message("\n=== Retrieving PTM & processing ===")
ptm_result <- safe_call("GetPTM_Processing", accessions, outdir)

message("\n=== Retrieving names & taxonomy ===")
taxa_result <- safe_call("GetNamesTaxa", accessions, outdir)

message("\n=== Retrieving structure info ===")
struct_result <- safe_call("GetStructureInfo", accessions, outdir)

message("\n=== Retrieving protein interactions ===")
interact_result <- safe_call("GetProteinInteractions", accessions, outdir)

message("\n=== Retrieving publications ===")
pub_result <- safe_call("GetPublication", accessions, outdir)

message("\n=== Retrieving miscellaneous ===")
misc_result <- safe_call("GetMiscellaneous", accessions, outdir)

# ---- Combine key results into a single summary TSV ----
message("\n=== Building combined summary ===")

# Start with accession list
combined <- data.frame(accession = accessions, stringsAsFactors = FALSE)

# Merge each result by row name (accession)
merge_by_rowname <- function(combined, result, label) {
  if (!is.null(result) && nrow(result) > 0) {
    result$accession <- rownames(result)
    # Prefix columns to avoid clashes
    orig_cols <- setdiff(colnames(result), "accession")
    colnames(result)[colnames(result) %in% orig_cols] <- paste0(label, ".", orig_cols)
    combined <- merge(combined, result, by = "accession", all.x = TRUE)
  }
  return(combined)
}

combined <- merge_by_rowname(combined, taxa_result, "taxa")
combined <- merge_by_rowname(combined, func_result, "func")
combined <- merge_by_rowname(combined, go_result, "go")
combined <- merge_by_rowname(combined, loc_result, "loc")
combined <- merge_by_rowname(combined, path_result, "pathology")
combined <- merge_by_rowname(combined, fam_result, "family")
combined <- merge_by_rowname(combined, ptm_result, "ptm")
combined <- merge_by_rowname(combined, struct_result, "struct")

summary_file <- file.path(outdir, paste0(prefix, "_uniprotr_combined.tsv"))
write.table(combined, summary_file, sep = "\t", row.names = FALSE, quote = FALSE)
message(paste("Wrote combined summary:", summary_file, "->", ncol(combined), "columns"))

message("\nDone.")
