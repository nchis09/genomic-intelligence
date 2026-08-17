#!/usr/bin/env Rscript
#
# annotate_phylo_trees.R
#
# Produce fully-annotated phylogenetic tree figures (PNG + SVG) for a given
# species using data from the knowledge warehouse DuckDB export.
#
# Two trees are generated:
#   1. Augur evolutionary tree (time-resolved, from Nextstrain auspice)
#   2. IQ-TREE bootstrap tree (ML with bootstrap support)
#
# Annotation layers (applied via ggtree + ggtreeExtra):
#   - Clade coloring (tip + branch)
#   - Outbreak labels (outer ring)
#   - Country/geography (outer ring)
#   - Collection date (outer ring)
#   - Genome coverage (bar ring)
#   - Mutation counts (bar ring)
#   - Query sample highlighting (tip shape)
#   - Bootstrap support (IQ-TREE only, node labels)
#
# Usage:
#   Rscript bin/annotate_phylo_trees.R \
#     --duckdb results/knowledge_warehouse/knowledge_warehouse.duckdb \
#     --species bdbv \
#     --outdir results/pathogen_genomics/bdbv \
#     [--augur-results-dir <path>]   # optional: augur intermediate files
#     [--dark-mode]                   # produce dark-background variant
#
# Requires: pgirl_ggtree conda environment

suppressPackageStartupMessages({
  library(optparse)
  library(ape)
  library(treeio)
  library(ggtree)
  library(ggtreeExtra)
  library(ggplot2)
  library(ggnewscale)
  library(dplyr)
  library(tidyr)
  library(RColorBrewer)
})

# ---------------------------------------------------------------------------
# CLI arguments
# ---------------------------------------------------------------------------
option_list <- list(
  make_option("--duckdb", type = "character", help = "Path to DuckDB file"),
  make_option("--species", type = "character", help = "Species code (e.g. bdbv, sudv)"),
  make_option("--outdir", type = "character", help = "Output directory for tree figures"),
  make_option("--augur-results-dir", type = "character", default = NULL,
              help = "Optional path to augur results dir with branch_lengths.json, muts.json"),
  make_option("--dark-mode", action = "store_true", default = FALSE,
              help = "Produce dark-background variant")
)
opt <- parse_args(OptionParser(option_list = option_list))

stopifnot(!is.null(opt$duckdb), !is.null(opt$species), !is.null(opt$outdir))

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# Connect to DuckDB
# ---------------------------------------------------------------------------
if (!requireNamespace("duckdb", quietly = TRUE)) {
  # Fall back: try DBI + duckdb
  stop("Package 'duckdb' is required. Add r-duckdb to conda env.")
}
library(DBI)
con <- DBI::dbConnect(duckdb::duckdb(), dbdir = opt$duckdb, read_only = TRUE)
on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)

# ---------------------------------------------------------------------------
# Load tree data
# ---------------------------------------------------------------------------
species_filter <- tolower(opt$species)

trees_df <- DBI::dbGetQuery(con, sprintf(
  "SELECT tree_id, tree_method, newick FROM phylogenetic_trees WHERE LOWER(species) = '%s'",
  species_filter
))

tips_all <- DBI::dbGetQuery(con, sprintf(
  "SELECT t.* FROM tree_tips t
   JOIN phylogenetic_trees pt ON t.tree_id = pt.tree_id
   WHERE LOWER(pt.species) = '%s'",
  species_filter
))

if (nrow(trees_df) == 0) {
  message("No trees found for species: ", opt$species)
  quit(status = 0)
}

# ---------------------------------------------------------------------------
# Color palettes
# ---------------------------------------------------------------------------
make_clade_palette <- function(clades) {
  n <- length(unique(clades))
  if (n <= 8) {
    cols <- brewer.pal(max(3, n), "Set2")[1:n]
  } else if (n <= 12) {
    cols <- brewer.pal(n, "Set3")
  } else {
    cols <- colorRampPalette(brewer.pal(12, "Set3"))(n)
  }
  setNames(cols, sort(unique(clades)))
}

make_outbreak_palette <- function(outbreaks) {
  n <- length(unique(outbreaks))
  if (n <= 9) {
    cols <- brewer.pal(max(3, n), "Pastel1")[1:n]
  } else {
    cols <- colorRampPalette(brewer.pal(9, "Pastel1"))(n)
  }
  setNames(cols, sort(unique(outbreaks)))
}

# ---------------------------------------------------------------------------
# Theme helpers
# ---------------------------------------------------------------------------
dark_mode <- isTRUE(opt$`dark-mode`)

bg_color   <- if (dark_mode) "#1e1e2e" else "white"
text_color <- if (dark_mode) "#e0e0e0" else "#333333"
grid_color <- if (dark_mode) "#3a3a4a" else "#e0e0e0"

tree_theme <- theme(

  plot.background  = element_rect(fill = bg_color, color = NA),
  panel.background = element_rect(fill = bg_color, color = NA),
  legend.background = element_rect(fill = bg_color, color = NA),
  legend.text  = element_text(color = text_color, size = 7),
  legend.title = element_text(color = text_color, size = 8, face = "bold"),
  plot.title   = element_text(color = text_color, size = 12, face = "bold"),
  plot.subtitle = element_text(color = text_color, size = 9)
)

# ---------------------------------------------------------------------------
# Helper: parse newick + join tip metadata → treedata
# ---------------------------------------------------------------------------
build_treedata <- function(newick_str, tips_df) {
  # Write newick to temp file for ape::read.tree

  tmp <- tempfile(fileext = ".nwk")
  writeLines(newick_str, tmp)
  tr <- tryCatch(read.tree(tmp), error = function(e) NULL)
  unlink(tmp)
  if (is.null(tr)) return(NULL)

  # Join tip metadata
  tip_labels <- tr$tip.label
  tips_matched <- tips_df %>%
    filter(label %in% tip_labels) %>%
    distinct(label, .keep_all = TRUE)

  # Create treedata with metadata
  td <- as.treedata(tr)
  td
}

# ---------------------------------------------------------------------------
# Helper: produce annotated tree plot
# ---------------------------------------------------------------------------
render_annotated_tree <- function(newick_str, tips_df, tree_method, title_prefix) {
  tmp <- tempfile(fileext = ".nwk")
  writeLines(newick_str, tmp)
  tr <- tryCatch(read.tree(tmp), error = function(e) NULL)
  unlink(tmp)
  if (is.null(tr)) {
    message("  Could not parse newick for method: ", tree_method)
    return(NULL)
  }

  # Match tip metadata — select only needed columns with `label` first
  # (ggtree's %<+% operator expects the join key as the first column and
  # will fail if `label` appears more than once in the data.frame)
  keep_cols <- intersect(
    c("label", "is_query", "clade", "outbreak", "country",
      "genome_coverage", "nuc_mutation_count", "aa_mutation_count", "div",
      "nextclade_qc", "tip_date"),
    names(tips_df)
  )
  tip_meta <- tips_df %>%
    select(all_of(keep_cols)) %>%
    filter(label %in% tr$tip.label) %>%
    distinct(label, .keep_all = TRUE) %>%
    mutate(
      clade = ifelse(is.na(clade) | clade == "", "Unknown", clade),
      outbreak = ifelse(is.na(outbreak) | outbreak == "", "Unknown", outbreak),
      country = ifelse(is.na(country) | country == "", "Unknown", country),
      is_query = ifelse(is.na(is_query), FALSE, is_query),
      genome_coverage = ifelse(is.na(genome_coverage), 0, genome_coverage),
      nuc_mutation_count = ifelse(is.na(nuc_mutation_count), 0, as.numeric(nuc_mutation_count)),
      aa_mutation_count = ifelse(is.na(aa_mutation_count), 0, as.numeric(aa_mutation_count)),
      total_mutations = nuc_mutation_count + aa_mutation_count
    )

  # Determine if we have enough metadata
  has_metadata <- nrow(tip_meta) > 0

  # Palettes
  clade_pal <- if (has_metadata) make_clade_palette(tip_meta$clade) else NULL
  outbreak_pal <- if (has_metadata) make_outbreak_palette(tip_meta$outbreak) else NULL

  # Determine number of tips for sizing

  n_tips <- length(tr$tip.label)
  fig_height <- max(8, min(30, n_tips * 0.15))
  fig_width <- max(12, min(20, 12 + n_tips * 0.02))

  # Base tree
  p <- ggtree(tr, layout = "rectangular", size = 0.4,
              color = if (dark_mode) "#8899aa" else "#555555") +
    labs(
      title = paste0(title_prefix, " \u2014 ", toupper(opt$species)),
      subtitle = paste0(n_tips, " tips | Method: ", tree_method)
    ) +
    tree_theme

  # Add bootstrap support for IQ-TREE trees (values in node labels)
  if (grepl("iqtree", tolower(tree_method))) {
    bs_vals <- suppressWarnings(as.numeric(tr$node.label))
    bs_vals[is.na(bs_vals)] <- 0
    # Identify internal node IDs with high bootstrap support
    high_bs_nodes <- (n_tips + 1):(n_tips + tr$Nnode)
    high_bs_nodes <- high_bs_nodes[bs_vals >= 70]
    if (length(high_bs_nodes) > 0) {
      p <- p + geom_point2(
        aes(subset = (node %in% !!high_bs_nodes)),
        color = if (dark_mode) "#66bb6a" else "#2e7d32",
        size = 2, alpha = 0.7, shape = 16
      )
    }
  }

  if (!has_metadata) {
    # No metadata — just return the base tree with tip labels
    p <- p + geom_tiplab(size = 2, color = text_color)
    return(list(plot = p, height = fig_height, width = fig_width))
  }

  # Merge tip metadata into the tree's internal data via %<+%.
  # To avoid column-name collisions (ggtree warns and can drop columns),
  # we join only `label` + prefixed columns, then use the prefixed names
  # in all subsequent aes() calls.
  ann <- tip_meta
  data_cols <- setdiff(names(ann), "label")
  prefixed  <- paste0("pg_", data_cols)
  names(ann)[names(ann) %in% data_cols] <- prefixed
  p <- p %<+% ann

  # Tip points colored by clade, shaped by query status
  p <- p +
    geom_tippoint(
      aes(color = pg_clade, shape = pg_is_query),
      size = 2, alpha = 0.85
    ) +
    scale_color_manual(values = clade_pal, name = "Clade") +
    scale_shape_manual(values = c("TRUE" = 18, "FALSE" = 16),
                       labels = c("TRUE" = "Query", "FALSE" = "Background"),
                       name = "Sample Type")

  # Tip labels for query samples only (to avoid clutter)
  p <- p + geom_tiplab(
    aes(subset = pg_is_query, label = label),
    size = 2.5, offset = 0.001,
    color = if (dark_mode) "#ffab40" else "#e65100",
    fontface = "bold"
  )

  # --- Outer annotation rings via ggtreeExtra ---
  # Use the original (unprefixed) tip_meta for geom_fruit data argument

  # Ring 1: Outbreak (heatmap tile)
  p <- p + new_scale_fill() +
    geom_fruit(
      data = tip_meta,
      geom = geom_tile,
      mapping = aes(y = label, fill = outbreak),
      width = 0.02, offset = 0.05
    ) +
    scale_fill_manual(values = outbreak_pal, name = "Outbreak")

  # Ring 2: Country (heatmap tile)
  country_vals <- unique(tip_meta$country)
  n_countries <- length(country_vals)
  country_cols <- if (n_countries <= 9) {
    setNames(brewer.pal(max(3, n_countries), "Set1")[1:n_countries], sort(country_vals))
  } else {
    setNames(colorRampPalette(brewer.pal(9, "Set1"))(n_countries), sort(country_vals))
  }

  p <- p + new_scale_fill() +
    geom_fruit(
      data = tip_meta,
      geom = geom_tile,
      mapping = aes(y = label, fill = country),
      width = 0.02, offset = 0.02
    ) +
    scale_fill_manual(values = country_cols, name = "Country")

  # Ring 3: Genome coverage (bar)
  p <- p + new_scale_fill() +
    geom_fruit(
      data = tip_meta,
      geom = geom_bar,
      mapping = aes(y = label, x = genome_coverage),
      stat = "identity", width = 0.6, offset = 0.02,
      fill = if (dark_mode) "#4fc3f7" else "#0277bd",
      alpha = 0.7
    )

  # Ring 4: Mutation count (bar)
  p <- p + new_scale_fill() +
    geom_fruit(
      data = tip_meta,
      geom = geom_bar,
      mapping = aes(y = label, x = total_mutations),
      stat = "identity", width = 0.6, offset = 0.02,
      fill = if (dark_mode) "#ef5350" else "#b71c1c",
      alpha = 0.7
    )

  return(list(plot = p, height = fig_height, width = fig_width))
}

# ---------------------------------------------------------------------------
# Generate trees
# ---------------------------------------------------------------------------
message("Generating annotated phylogenetic trees for species: ", opt$species)
message("  Found ", nrow(trees_df), " tree(s) in database")

for (i in seq_len(nrow(trees_df))) {
  row <- trees_df[i, ]
  method <- row$tree_method
  newick <- row$newick
  tree_id <- row$tree_id

  if (is.na(newick) || nchar(trimws(newick)) == 0) {
    message("  Skipping tree_id=", tree_id, " (empty newick)")
    next
  }

  # Get tips for this specific tree
  tips_for_tree <- tips_all %>% filter(tree_id == !!tree_id)

  # Determine title prefix
  title_prefix <- if (grepl("iqtree", tolower(method))) {
    "IQ-TREE Bootstrap Tree"
  } else if (grepl("augur|nextstrain|auspice", tolower(method))) {
    "Augur Evolutionary Tree"
  } else {
    paste0("Phylogenetic Tree (", method, ")")
  }

  message("  Rendering: ", title_prefix, " (", nrow(tips_for_tree), " tips in DB)")

  result <- render_annotated_tree(newick, tips_for_tree, method, title_prefix)

  if (is.null(result)) next

  # File naming
  method_slug <- gsub("[^a-z0-9]", "_", tolower(method))
  suffix <- if (dark_mode) "_dark" else ""
  base_name <- paste0("annotated_tree_", method_slug, suffix)

  # Save PNG

  png_path <- file.path(opt$outdir, paste0(base_name, ".png"))
  ggsave(png_path, plot = result$plot,
         width = result$width, height = result$height,
         dpi = 300, bg = bg_color)
  message("  Saved: ", png_path)

  # Save SVG (requires svglite)
  if (requireNamespace("svglite", quietly = TRUE)) {
    svg_path <- file.path(opt$outdir, paste0(base_name, ".svg"))
    ggsave(svg_path, plot = result$plot,
           width = result$width, height = result$height,
           bg = bg_color)
    message("  Saved: ", svg_path)
  } else {
    message("  SKIP SVG: package 'svglite' not installed")
  }
}

message("Done.")
