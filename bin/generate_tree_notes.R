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

write_note <- function(path, text, source, model = NULL, error = NULL) {
  payload <- list(
    species = species, text = text, source = source,
    model = model, error = error,
    generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  )
  writeLines(jsonlite::toJSON(payload, auto_unbox = TRUE, pretty = TRUE), path)
  message("[tree_notes] wrote ", path, " (source=", source, ")")
}

fail <- function(msg) {
  message("[tree_notes] ", msg)
  write_note(out_json,
             paste0("No tree data available to interpret. ", msg),
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
  "SELECT tree_id, tree_method, newick, segment FROM phylogenetic_trees WHERE LOWER(species) = ?",
  params = list(species))
tips <- DBI::dbGetQuery(con,
  "SELECT t.* FROM tree_tips t JOIN phylogenetic_trees pt ON t.tree_id = pt.tree_id WHERE LOWER(pt.species) = ?",
  params = list(species))

if (nrow(trees) == 0) fail("No trees for species ", species)
if (!("segment" %in% names(trees))) trees$segment <- NA_character_

# Prefer nextstrain/augur (full metadata); skip iqtree trees
pick_tree <- function(cand) {
  non_iq <- cand[!grepl("iqtree", tolower(cand$tree_method)), , drop = FALSE]
  ranks <- match(tolower(non_iq$tree_method), c("nextstrain", "augur", "nextclade"), nomatch = 99)
  if (nrow(non_iq) == 0) cand[1, ] else non_iq[order(ranks), ][1, ]
}

# Builds the note for one selected tree; NULL when the tree is unusable.
note_for <- function(sel) {
  if (is.na(sel$newick) || nchar(trimws(sel$newick)) == 0) return(NULL)

  tmp <- tempfile(fileext = ".nwk")
  writeLines(sel$newick, tmp)
  tr <- tryCatch(ape::read.tree(tmp), error = function(e) NULL)
  unlink(tmp)
  if (is.null(tr)) return(NULL)

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

  s <- tryCatch(.pg_evo_summary(tr, tip_meta), error = function(e) NULL)
  if (is.null(s)) return(NULL)

  tryCatch(.pg_tree_note(s), error = function(e)
    list(text = .pg_note_template(s), source = "template", model = NULL, error = e$message))
}

# -- Per-segment notes ---------------------------------------------------------
# Segmented pathogens (influenza) get one note per segment —
# tree_note_<seg>.json — plus tree_note.json mirroring the preferred segment
# (ha when present) so species-level consumers keep working.
seg_key <- ifelse(is.na(trees$segment) | trees$segment == "", "genome",
                  tolower(trees$segment))
segs <- unique(seg_key)
pref_seg <- if ("ha" %in% segs) "ha" else segs[1]

wrote_any <- FALSE
for (seg in segs) {
  sel <- pick_tree(trees[seg_key == seg, , drop = FALSE])
  note <- note_for(sel)
  if (is.null(note)) {
    message("[tree_notes] no usable tree for segment ", seg)
    next
  }
  fname <- if (seg == pref_seg) "tree_note.json" else paste0("tree_note_", seg, ".json")
  write_note(file.path(opt$outdir, fname), note$text, note$source,
             note$model, note$error)
  wrote_any <- TRUE
}

if (!wrote_any) fail("Could not build evolutionary summary for any tree.")
