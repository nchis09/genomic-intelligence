#!/usr/bin/env Rscript
# generate_tree_notes.R
#
# Pre-generate the plain-language phylogenetic-tree interpretation note for
# one species, so the dashboard can render it instantly instead of calling
# Ollama at view time.
#
# Reads the exported DuckDB knowledge warehouse, picks the best tree
# (nextstrain > augur > nextclade, skipping iqtree), builds the evolutionary
# summary via dashboard/modules/llm_note.R, calls the local Ollama model, and
# writes tree_note.json. If Ollama is unreachable the deterministic template
# note is written instead — the file is always produced.
#
# Usage:
#   Rscript generate_tree_notes.R --duckdb knowledge_warehouse.duckdb \
#       --species ebov --outdir . [--llm-note-r path/to/llm_note.R]

suppressWarnings(suppressMessages({
  library(optparse)
  library(DBI)
}))

option_list <- list(
  make_option("--duckdb",      type = "character", help = "Path to knowledge_warehouse.duckdb"),
  make_option("--species",     type = "character", help = "Species key (e.g. ebov, bdbv, sudv)"),
  make_option("--outdir",      type = "character", default = ".", help = "Output directory"),
  make_option("--llm-note-r",  type = "character", default = "llm_note.R",
              help = "Path to dashboard/modules/llm_note.R (staged by Nextflow)")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$duckdb) || is.null(opt$species))
  stop("--duckdb and --species are required")

if (!file.exists(opt$`llm-note-r`))
  stop("llm_note.R not found at: ", opt$`llm-note-r`)
source(opt$`llm-note-r`)

species <- tolower(opt$species)
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
out_json <- file.path(opt$outdir, "tree_note.json")

write_note <- function(text, source, model = NULL, error = NULL) {
  payload <- list(
    species = species, text = text, source = source,
    model = model, error = error,
    generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  )
  writeLines(jsonlite::toJSON(payload, auto_unbox = TRUE, pretty = TRUE), out_json)
  message("[tree_notes] wrote ", out_json, " (source=", source, ")")
}

fail <- function(msg) {
  message("[tree_notes] ", msg)
  write_note(paste0("No tree data available to interpret. ", msg),
             source = "none", error = msg)
  quit(save = "no", status = 0)
}

# -- Load tree + tips ---------------------------------------------------------
con <- tryCatch(
  DBI::dbConnect(duckdb::duckdb(), dbdir = opt$duckdb, read_only = TRUE),
  error = function(e) NULL
)
if (is.null(con)) fail("Cannot open DuckDB.")
on.exit(try(DBI::dbDisconnect(con, shutdown = TRUE), silent = TRUE), add = TRUE)

trees <- DBI::dbGetQuery(con,
  "SELECT tree_id, tree_method, newick FROM phylogenetic_trees WHERE LOWER(species) = ?",
  params = list(species))
tips <- DBI::dbGetQuery(con,
  "SELECT t.* FROM tree_tips t JOIN phylogenetic_trees pt ON t.tree_id = pt.tree_id WHERE LOWER(pt.species) = ?",
  params = list(species))

if (nrow(trees) == 0) fail("No trees for species ", species)

# Prefer nextstrain/augur (full metadata); skip iqtree trees
non_iq <- trees[!grepl("iqtree", tolower(trees$tree_method)), , drop = FALSE]
ranks <- match(tolower(non_iq$tree_method), c("nextstrain", "augur", "nextclade"), nomatch = 99)
sel <- if (nrow(non_iq) == 0) trees[1, ] else non_iq[order(ranks), ][1, ]

if (is.na(sel$newick) || nchar(trimws(sel$newick)) == 0) fail("Tree newick is empty.")

tmp <- tempfile(fileext = ".nwk")
writeLines(sel$newick, tmp)
tr <- tryCatch(ape::read.tree(tmp), error = function(e) NULL)
unlink(tmp)
if (is.null(tr)) fail("Failed to parse tree newick.")

# Same normalization as the dashboard: clean labels, ladderize, cap long branches
tr$tip.label <- gsub("[,;():\\[\\] ]", "_", tr$tip.label)
tr <- ape::ladderize(tr)
if (!is.null(tr$edge.length) && length(tr$edge.length) > 0) {
  q95 <- quantile(tr$edge.length[tr$edge.length > 0], 0.95, na.rm = TRUE)
  tr$edge.length[tr$edge.length > 3 * q95] <- 3 * q95
}

tip_meta <- tips[tips$tree_id == sel$tree_id, ]
if (nrow(tip_meta) == 0 && nrow(tips) > 0)
  tip_meta <- tips[!duplicated(tips$label), ]
keep <- intersect(c("label", "is_query", "clade", "outbreak", "country", "tip_date",
                    "genome_coverage", "nuc_mutation_count", "aa_mutation_count"),
                  names(tip_meta))
tip_meta <- tip_meta[, keep, drop = FALSE]
tip_meta$label <- gsub("[,;():\\[\\] ]", "_", tip_meta$label)
tip_meta <- tip_meta[!duplicated(tip_meta$label), ]
for (cc in intersect(c("clade", "outbreak", "country"), names(tip_meta)))
  tip_meta[[cc]][is.na(tip_meta[[cc]]) | tip_meta[[cc]] == ""] <- "Unknown"
tip_meta$is_query[is.na(tip_meta$is_query)] <- FALSE

# -- Summary + note -----------------------------------------------------------
s <- tryCatch(.pg_evo_summary(tr, tip_meta), error = function(e) NULL)
if (is.null(s)) fail("Could not build evolutionary summary.")

note <- tryCatch(.pg_tree_note(s), error = function(e)
  list(text = .pg_note_template(s), source = "template", model = NULL, error = e$message))

write_note(note$text, note$source, note$model, note$error)
