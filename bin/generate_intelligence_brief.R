#!/usr/bin/env Rscript
# generate_intelligence_brief.R
#
# Pre-generate the per-species "Intelligence Brief" (intelligence_brief.json)
# so the dashboard renders it instantly instead of composing it at view time.
#
# Assembles a structured fact list from the species' pipeline outputs
# (pathogen identification TSVs, transmission-context TSVs, mutation-profile
# TSVs, tree_note.json, spread_assessment.json, and the DuckDB warehouse for
# literature-domain counts), then asks dashboard/modules/intelligence_brief.R
# for the bottom-line narrative — local Ollama first, deterministic template
# fallback. The JSON is always written.
#
# Usage:
#   Rscript generate_intelligence_brief.R --species bdbv \
#       --pi_dir species_identification --tc_dir bdbv \
#       --mp_dir . --tree_note tree_note.json --spread_json spread_assessment.json \
#       --duckdb knowledge_warehouse.duckdb --outdir . \
#       --llm_note_r llm_note.R --module_r intelligence_brief.R

suppressWarnings(suppressMessages({
  library(readr)
  library(dplyr)
}))

if (requireNamespace("optparse", quietly = TRUE)) {
  option_list <- list(
    optparse::make_option("--species",     type = "character", help = "Species key (e.g. ebov, bdbv, sudv)"),
    optparse::make_option("--pi_dir",      type = "character", default = NULL, help = "Dir with species_identification TSVs"),
    optparse::make_option("--tc_dir",      type = "character", default = NULL, help = "transmission_context/<species> dir"),
    optparse::make_option("--mp_dir",      type = "character", default = NULL, help = "Dir with mutation_profile TSVs"),
    optparse::make_option("--tree_note",   type = "character", default = NULL, help = "tree_note.json path"),
    optparse::make_option("--spread_json", type = "character", default = NULL, help = "spread_assessment.json path"),
    optparse::make_option("--duckdb",      type = "character", default = NULL, help = "knowledge_warehouse.duckdb path"),
    optparse::make_option("--lit_dir",     type = "character", default = NULL, help = "literature_metadata/<species> dir (fallback evidence-base counts)"),
    optparse::make_option("--outdir",      type = "character", default = ".", help = "Output directory"),
    optparse::make_option("--llm_note_r",  type = "character", default = "llm_note.R",
                help = "Path to dashboard/modules/llm_note.R (staged by Nextflow)"),
    optparse::make_option("--module_r",    type = "character", default = "intelligence_brief.R",
                help = "Path to dashboard/modules/intelligence_brief.R (staged by Nextflow)")
  )
  opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))
} else {
  # Minimal fallback: parse "--key value" pairs (optparse not installed).
  args <- commandArgs(trailingOnly = TRUE)
  opt <- list(species = NULL, pi_dir = NULL, tc_dir = NULL, mp_dir = NULL,
              tree_note = NULL, spread_json = NULL, duckdb = NULL, lit_dir = NULL,
              outdir = ".",
              llm_note_r = "llm_note.R", module_r = "intelligence_brief.R")
  i <- 1
  while (i <= length(args)) {
    if (startsWith(args[i], "--") && i < length(args)) {
      opt[[sub("^--", "", args[i])]] <- args[i + 1]
      i <- i + 2
    } else i <- i + 1
  }
}

if (is.null(opt$species)) stop("--species is required")
if (!file.exists(opt$llm_note_r))  stop("llm_note.R not found at: ", opt$llm_note_r)
if (!file.exists(opt$module_r))    stop("intelligence_brief.R not found at: ", opt$module_r)
source(opt$llm_note_r)
source(opt$module_r)

species <- tolower(opt$species)
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
out_json <- file.path(opt$outdir, "intelligence_brief.json")

write_brief <- function(payload) {
  # digits = NA keeps full precision — the default (4) rounds small p-values
  # and genetic distances to 0.
  writeLines(jsonlite::toJSON(payload, auto_unbox = TRUE, pretty = TRUE,
                              digits = NA), out_json)
  message("[intelligence_brief] wrote ", out_json)
}

empty_brief <- function(species, msg) {
  list(species = species,
       generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
       narrative = list(text = msg, source = "none", model = "", error = msg),
       caveats = list(msg), data_gaps = list("all"))
}

payload <- tryCatch(
  pg_build_intelligence_brief(
    species          = species,
    pi_dir           = opt$pi_dir,
    tc_dir           = opt$tc_dir,
    mp_dir           = opt$mp_dir,
    tree_note_path   = opt$tree_note,
    spread_json_path = opt$spread_json,
    duckdb_path      = opt$duckdb,
    lit_dir          = opt$lit_dir
  ),
  error = function(e) {
    message("[intelligence_brief] assembly failed: ", e$message)
    empty_brief(species, paste("Brief assembly failed:", e$message))
  }
)

write_brief(payload)
message("[intelligence_brief] ", species, " narrative source: ",
        payload$narrative$source %||% "unknown")
