#!/usr/bin/env Rscript
# Plot an annotated rectangular phylogeny from a Newick tree and tip metadata.
#
# Branches are colored by the most common country of their descendant tips
# (inferred by majority rule). A binary heatmap of the top 5 amino-acid
# mutations per protein can be attached to the right of the tree.
# Query tips are highlighted with red diamonds.
#
# Example:
#   Rscript plot_ggtreeExtra.R
#     --tree=output/nextstrain_ebola/results/bdbv/all-outbreaks/tree.nwk
#     --metadata=output/nextstrain_ebola/annotations/bdbv_tip_metadata.tsv
#     --mutation-matrix=output/nextstrain_ebola/annotations/bdbv_mutation_matrix.tsv
#     --output=output/nextstrain_ebola/figures/bdbv_ggtreeExtra.png
#     --title="BDBV annotated tree"
#     --layout=rectangular

suppressPackageStartupMessages({
  library(ggtree)
  library(ggtreeExtra)
  library(treeio)
  library(ggplot2)
  library(ggnewscale)
  library(dplyr)
  library(tidyr)
  library(RColorBrewer)
  library(rlang)
})

# --- CLI helpers -----------------------------------------------------------
parse_arg <- function(args, key, default = NULL) {
  pattern <- paste0("^--", key, "=")
  hit <- grep(pattern, args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^--", key, "="), "", hit[length(hit)])
}

has_flag <- function(args, key) {
  any(grepl(paste0("^--", key, "$"), args))
}

args <- commandArgs(trailingOnly = TRUE)

tree_file       <- parse_arg(args, "tree")
meta_file       <- parse_arg(args, "metadata")
out_file        <- parse_arg(args, "output", "ggtreeExtra_plot.png")
title_text      <- parse_arg(args, "title", "Annotated phylogeny")
mut_matrix_file <- parse_arg(args, "mutation-matrix", NULL)
query_var       <- parse_arg(args, "query-var", "is_query")
layout          <- parse_arg(args, "layout", "rectangular")
branch_length   <- parse_arg(args, "branch-length", "branch.length")
open_angle      <- as.numeric(parse_arg(args, "open-angle", "30"))
tip_size        <- as.numeric(parse_arg(args, "tip-size", "1.0"))
fig_width       <- as.numeric(parse_arg(args, "width", "18"))
fig_height      <- as.numeric(parse_arg(args, "height", "14"))
fig_dpi         <- as.integer(parse_arg(args, "dpi", "200"))

if (is.null(tree_file) || is.null(meta_file)) {
  stop("Usage: Rscript plot_ggtreeExtra.R --tree=<nwk> --metadata=<tsv> [--mutation-matrix=...] [--output=...] [--title=...] [--query-var=is_query] [--layout=rectangular] [--branch-length=branch.length] [--open-angle=30] [--tip-size=1.0] [--width=18] [--height=14] [--dpi=200]",
       call. = FALSE)
}

# --- Load data -------------------------------------------------------------
tree <- read.tree(tree_file)
meta <- read.delim(
  meta_file,
  stringsAsFactors = FALSE,
  check.names      = FALSE,
  na.strings       = c("", "NA", "None", "null")
)

if (!"label" %in% names(meta)) {
  stop("metadata TSV must contain a 'label' column matching tree tip labels")
}

# Convert is_query to logical if present
if (query_var %in% names(meta)) {
  meta[[query_var]] <- as.logical(meta[[query_var]])
  meta[[query_var]][is.na(meta[[query_var]])] <- FALSE
}

# Clean country column (used for branch coloring)
if ("country" %in% names(meta)) {
  meta$country <- as.character(meta$country)
  meta$country[meta$country %in% c("", "NA", "None", "null")] <- NA
}

# Keep only tips present in the tree
meta <- meta %>% dplyr::filter(label %in% tree$tip.label)
missing_tips <- setdiff(tree$tip.label, meta$label)
if (length(missing_tips) > 0) {
  warning(
    length(missing_tips), " tree tip(s) not found in metadata; they will be greyed out"
  )
}

# --- Infer ancestral country for internal nodes ----------------------------
infer_node_country <- function(tree, meta, country_col = "country") {
  n_tips <- length(tree$tip.label)
  n_nodes <- tree$Nnode
  labels <- c(tree$tip.label, tree$node.label)
  tip_country <- meta[[country_col]]
  names(tip_country) <- meta$label

  all_desc <- tidytree::offspring(tree, .node = 1:(n_tips + n_nodes), type = "tips")
  country <- vapply(1:(n_tips + n_nodes), function(i) {
    if (i <= n_tips) {
      return(tip_country[labels[i]])
    }
    desc <- all_desc[[i]]
    vals <- tip_country[tree$tip.label[desc]]
    vals <- vals[!is.na(vals)]
    if (length(vals) == 0) return(NA_character_)
    tbl <- sort(table(vals), decreasing = TRUE)
    names(tbl)[1]
  }, character(1))

  country[country == "NULL"] <- NA_character_
  data.frame(
    node    = 1:(n_tips + n_nodes),
    label   = labels,
    country = country,
    stringsAsFactors = FALSE
  )
}

node_country <- infer_node_country(tree, meta)

# Carry the query flag into the node data so it survives full_join
if (query_var %in% names(meta)) {
  node_country[[query_var]] <- meta[[query_var]][match(node_country$label, meta$label)]
  node_country[[query_var]][is.na(node_country[[query_var]])] <- FALSE
}

# --- Merge country onto the tree object ------------------------------------
tree2 <- treeio::full_join(tree, node_country, by = "label")

# --- Base plot -------------------------------------------------------------
tree_args <- list(
  mapping       = aes(color = country),
  layout        = layout,
  branch.length = branch_length,
  size          = 0.3
)
if (layout %in% c("fan", "circular", "radial")) {
  tree_args$open.angle <- open_angle
}
p <- do.call(ggtree, c(list(tree2), tree_args))

# --- Tip points and query highlight ----------------------------------------
p <- p +
  geom_tippoint(size = tip_size, shape = 16) +
  geom_tippoint(
    aes(subset = !!sym(query_var)),
    color = "red",
    fill  = "red",
    size  = 4,
    shape = 23,
    stroke = 0.8
  )

# --- Country color scale ---------------------------------------------------
if ("country" %in% names(node_country)) {
  n_countries <- length(unique(na.omit(node_country$country)))
  if (n_countries <= 8) {
    p <- p + scale_color_brewer(name = "Country", palette = "Set1", na.value = "grey70")
  } else {
    p <- p + scale_color_discrete(name = "Country", na.value = "grey70")
  }
}

# --- Top-AA-mutation binary heatmap ----------------------------------------
if (!is.null(mut_matrix_file) && file.exists(mut_matrix_file)) {
  mut_mat <- read.delim(
    mut_matrix_file,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    row.names = 1
  )
  # keep only tips present in the tree
  mut_mat <- mut_mat[rownames(mut_mat) %in% tree$tip.label, , drop = FALSE]

  if (ncol(mut_mat) > 0) {
    # convert to presence/absence labels for a discrete fill scale
    mut_mat[mut_mat == "0"] <- "absent"
    mut_mat[mut_mat == "1"] <- "present"
    mut_mat[is.na(mut_mat)] <- "absent"

    p <- gheatmap(
      p,
      mut_mat,
      width             = 0.25,
      offset            = 0.01,
      color             = "grey80",
      colnames          = TRUE,
      colnames_position = "bottom",
      colnames_angle    = 90,
      colnames_offset_y = -1,
      font.size         = 3
    ) +
    scale_fill_manual(
      name   = "AA mutation",
      values = c(absent = "white", present = "steelblue"),
      breaks = c("present", "absent")
    )
  }
}

# --- Finalize and save -----------------------------------------------------
p <- p +
  ggtitle(title_text) +
  theme(
    plot.title        = element_text(hjust = 0.5, size = 16, face = "bold"),
    legend.position   = "right",
    legend.box        = "vertical",
    legend.key.size   = unit(0.4, "cm"),
    legend.text       = element_text(size = 8),
    legend.title      = element_text(size = 9)
  )

ggsave(
  out_file,
  plot     = p,
  width    = fig_width,
  height   = fig_height,
  units    = "in",
  dpi      = fig_dpi,
  limitsize = FALSE
)

message("Saved annotated tree to: ", normalizePath(out_file, mustWork = FALSE))
