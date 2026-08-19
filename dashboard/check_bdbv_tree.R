#!/usr/bin/env Rscript
# Quick diagnostic: check BDBV tree data in DuckDB

duckdb_path <- file.path("..", "results", "knowledge_warehouse", "knowledge_warehouse.duckdb")
cat("DuckDB path:", normalizePath(duckdb_path, mustWork = FALSE), "\n")
cat("File exists:", file.exists(duckdb_path), "\n\n")

if (!file.exists(duckdb_path)) {
  stop("DuckDB file not found!")
}

library(DBI)
library(duckdb)

con <- dbConnect(duckdb(), dbdir = duckdb_path, read_only = TRUE)
on.exit(dbDisconnect(con, shutdown = TRUE))

cat("Tables:", paste(dbListTables(con), collapse = ", "), "\n\n")

# All trees
trees <- dbGetQuery(con, "SELECT tree_id, species, tree_method, LENGTH(newick) AS nwk_len FROM phylogenetic_trees")
cat("=== All trees ===\n")
print(trees)

# BDBV non-iqtree trees (what Augur selection picks)
bdbv_augur <- dbGetQuery(con,
  "SELECT tree_id, tree_method, LENGTH(newick) AS nwk_len, LEFT(newick, 200) AS nwk_preview
   FROM phylogenetic_trees
   WHERE LOWER(species) = 'bdbv' AND LOWER(tree_method) NOT LIKE '%iqtree%'"
)
cat("\n=== BDBV Augur trees ===\n")
print(bdbv_augur)

# BDBV tips count
bdbv_tips <- dbGetQuery(con,
  "SELECT pt.tree_id, pt.tree_method, COUNT(*) AS n_tips
   FROM tree_tips t JOIN phylogenetic_trees pt ON t.tree_id = pt.tree_id
   WHERE LOWER(pt.species) = 'bdbv'
   GROUP BY pt.tree_id, pt.tree_method"
)
cat("\n=== BDBV tip counts ===\n")
print(bdbv_tips)

# Try parsing the BDBV Augur newick
if (nrow(bdbv_augur) > 0) {
  newick <- dbGetQuery(con, sprintf(
    "SELECT newick FROM phylogenetic_trees WHERE tree_id = %d", bdbv_augur$tree_id[1]
  ))$newick
  cat("\n=== Attempting to parse BDBV Augur newick (", nchar(newick), " chars) ===\n")
  tmp <- tempfile(fileext = ".nwk")
  writeLines(newick, tmp)
  tr <- tryCatch(ape::read.tree(tmp), error = function(e) {
    cat("PARSE ERROR:", e$message, "\n")
    NULL
  })
  unlink(tmp)
  if (!is.null(tr)) {
    cat("Parsed OK:", length(tr$tip.label), "tips,", tr$Nnode, "internal nodes\n")
  }
}

cat("\nDone.\n")
