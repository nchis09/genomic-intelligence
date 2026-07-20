#!/usr/bin/env Rscript
# Decision-oriented figures for the Genomic Intelligence Engine.
#
# Reads output/intelligence_object.json (produced by the assessment pipeline)
# and writes 8 PNGs into output/figures/:
#   1. time_scaled_phylogeny.png
#   2. evidence_weighted_assessment.png  (transparent, evidence-traced scores)
#   3. source_attribution.png            (origin probability distribution)
#   4. genetic_relatedness.png           (SNP/temporal distances to nearest genomes)
#   5. mutation_phenotype_heatmap.png
#   6. evidence_consistency.png          (agreement / conflict across evidence streams)
#   7. knowledge_gaps.png                (prioritised missing-data gaps)
#   8. lineage_persistence_timeline.png
#
# Uses only CRAN packages available on most R installs: ape, ggplot2, dplyr,
# tidyr, jsonlite, RColorBrewer.  ggtree/ggraph fallbacks are avoided to keep
# the script portable. Generic geographic maps are replaced by quantitative
# source-attribution and relatedness figures.

suppressPackageStartupMessages({
  library(jsonlite)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ape)
  library(RColorBrewer)
})

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

# ---------------------------------------------------------------------------
# Config / CLI
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
INTEL_PATH <- if (length(args) >= 1) args[1] else "output/genomic_intelligence/intelligence_object.json"
OUT_DIR    <- if (length(args) >= 2) args[2] else "output/genomic_intelligence/figures"

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
safe <- function(x) {
  x[is.null(x)] <- NA
  x
}

read_intel <- function(path) {
  if (!file.exists(path)) stop("Intelligence object not found: ", path)
  fromJSON(path, simplifyDataFrame = FALSE, flatten = TRUE)
}

# Convert a list of named vectors/lists to a data.frame with id column.
list_to_df <- function(x, col_name = "value") {
  if (is.null(x) || length(x) == 0) return(data.frame(name = character(), value = numeric(), stringsAsFactors = FALSE))
  if (is.data.frame(x)) return(x)
  data.frame(
    name = names(x),
    value = unlist(x),
    stringsAsFactors = FALSE
  ) %>% filter(!is.na(name))
}

# ---------------------------------------------------------------------------
# 1. Time-scaled phylogeny
# ---------------------------------------------------------------------------
plot_time_scaled_phylogeny <- function(obj, out_path) {
  tree_file <- obj$tree_file
  if (is.null(tree_file) || !file.exists(tree_file)) {
    png(out_path, width = 1200, height = 800, res = 120)
    plot.new()
    title("No phylogenetic tree available")
    dev.off()
    return(invisible(NULL))
  }

  tree <- read.tree(tree_file)
  tips <- obj$tree_tips
  tip_df <- if (!is.null(tips) && length(tips) > 0) {
    bind_rows(lapply(tips, function(t) data.frame(
      name = safe(t$name),
      country = safe(t$country),
      date = safe(t$date),
      is_sample = isTRUE(t$is_sample),
      stringsAsFactors = FALSE
    )))
  } else data.frame(name = tree$tip.label, country = NA, is_sample = FALSE, stringsAsFactors = FALSE)

  # match order of tree tips
  tip_df <- tip_df[match(tree$tip.label, tip_df$name), ]
  tip_df$country[is.na(tip_df$country)] <- "unknown"
  countries <- unique(tip_df$country)
  pal <- colorRampPalette(brewer.pal(min(9, max(3, length(countries))), "Set1"))(length(countries))
  country_col <- setNames(pal, countries)

  png(out_path, width = 1400, height = 900, res = 120)
  par(mar = c(5, 1, 4, 1))
  max_height <- max(node.depth.edgelength(tree), na.rm = TRUE)
  plot(tree, show.tip.label = FALSE, edge.width = 1.2, x.lim = c(0, max_height * 1.15))
  tiplabels(pch = 21, col = "black", bg = country_col[tip_df$country], cex = 1.2)
  sample_idx <- which(tip_df$is_sample)
  if (length(sample_idx) > 0) {
    tiplabels(tip = sample_idx, pch = 23, col = "black", bg = "red", cex = 2)
  }
  axisPhylo(1, backward = FALSE)
  title("Time-scaled phylogeny (branch length ~ years from root)")
  legend("topright", legend = names(country_col), fill = country_col, bty = "n", cex = 0.8)
  dev.off()
}

# ---------------------------------------------------------------------------
# 2. Evidence-weighted assessment dashboard (transparent, evidence-traced)
# ---------------------------------------------------------------------------
plot_evidence_weighted_assessment <- function(obj, out_path) {
  threat <- obj$analyses$evidence_weighted_threat$metrics$scores
  if (is.null(threat)) {
    png(out_path, width = 1000, height = 600, res = 120)
    plot.new()
    title("No assessment scores available")
    dev.off()
    return(invisible(NULL))
  }

  df <- data.frame(
    criterion = names(threat),
    score = as.numeric(unlist(threat)),
    stringsAsFactors = FALSE
  ) %>%
    mutate(criterion = gsub("_", " ", criterion),
           criterion = tools::toTitleCase(criterion),
           criterion = reorder(criterion, score))

  p <- ggplot(df, aes(x = criterion, y = score, fill = score)) +
    geom_bar(stat = "identity", width = 0.7, color = "grey30") +
    scale_fill_gradientn(colors = c("#2c7bb6", "#ffffbf", "#d7191c"), limits = c(0, 5)) +
    geom_text(aes(label = score), hjust = -0.2, size = 4) +
    coord_flip() +
    scale_y_continuous(limits = c(0, 5.5), breaks = 0:5) +
    labs(
      title = "Evidence-weighted assessment",
      subtitle = paste(
        "Overall rank:",
        obj$analyses$evidence_weighted_threat$metrics$overall_label,
        "| Scores traceable to analyzer findings (score_basis)"
      ),
      x = NULL,
      y = "Score (0-5)",
      fill = "Score"
    ) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "none")

  ggsave(out_path, p, width = 10, height = 6, dpi = 150)
}

# ---------------------------------------------------------------------------
# 3. Source attribution: origin probability distribution and introduction scenarios
# ---------------------------------------------------------------------------
plot_source_attribution <- function(obj, out_path) {
  geo <- obj$analyses$phylogeographic_analysis$metrics
  if (is.null(geo) || is.null(geo$origin_support)) {
    png(out_path, width = 1000, height = 600, res = 120)
    plot.new()
    title("No source-attribution data available")
    dev.off()
    return(invisible(NULL))
  }

  support <- geo$origin_support
  prob <- support$origin_probability %||% list()

  if (length(prob) == 0 && !is.null(support$evidence)) {
    # Curated fallback: single point
    prob <- setNames(1, paste(support$evidence, "(curated)"))
  }

  if (length(prob) == 0) {
    png(out_path, width = 1000, height = 600, res = 120)
    plot.new()
    title("No origin probabilities available")
    dev.off()
    return(invisible(NULL))
  }

  df <- list_to_df(prob) %>%
    mutate(name = gsub("_", " ", name),
           name = tools::toTitleCase(name),
           name = reorder(name, value))

  df$label <- paste0(formatC(df$value * 100, format = "f", digits = 1), "%")

  p <- ggplot(df, aes(x = name, y = value, fill = value)) +
    geom_bar(stat = "identity", color = "grey30", width = 0.7) +
    scale_fill_gradientn(colors = c("#2c7bb6", "#ffffbf", "#d7191c"), limits = c(0, 1)) +
    geom_text(aes(label = label), hjust = -0.2, size = 4) +
    coord_flip() +
    scale_y_continuous(limits = c(0, 1.05), labels = function(x) paste0(x * 100, "%")) +
    labs(
      title = "Source attribution: inferred origin probabilities",
      subtitle = paste("Method:", support$method %||% "unknown"),
      x = NULL,
      y = "Probability"
    ) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "none")

  ggsave(out_path, p, width = 10, height = max(4, 1 + nrow(df) * 0.5), dpi = 150)
}

# ---------------------------------------------------------------------------
# 4. Genetic relatedness: SNP distance and temporal gap to nearest genomes
# ---------------------------------------------------------------------------
plot_genetic_relatedness <- function(obj, out_path) {
  gr <- obj$analyses$genetic_relatedness_analysis
  if (is.null(gr) || is.null(gr$metrics)) {
    png(out_path, width = 1000, height = 600, res = 120)
    plot.new()
    title("No genetic relatedness data available")
    dev.off()
    return(invisible(NULL))
  }

  snps <- gr$metrics$snps_to_closest_genome %||% NA
  gap <- gr$metrics$min_temporal_gap_days %||% NA

  df <- data.frame(
    metric = c("SNPs to closest genome", "Minimum temporal gap (days)"),
    value = c(as.numeric(snps), as.numeric(gap)),
    stringsAsFactors = FALSE
  ) %>% filter(!is.na(value))

  if (nrow(df) == 0) {
    png(out_path, width = 1000, height = 600, res = 120)
    plot.new()
    title("Genetic relatedness: no quantified distances available")
    dev.off()
    return(invisible(NULL))
  }

  df$metric <- factor(df$metric, levels = c("SNPs to closest genome", "Minimum temporal gap (days)"))
  df$label <- prettyNum(round(df$value, 1), big.mark = ",")

  # Each metric is on a different scale, so use free y-axis facets.
  p <- ggplot(df, aes(x = metric, y = value, fill = metric)) +
    geom_bar(stat = "identity", width = 0.6, color = "grey30") +
    geom_text(aes(label = label), vjust = -0.5, size = 4) +
    scale_fill_brewer(palette = "Set2") +
    facet_wrap(~metric, scales = "free_y") +
    labs(
      title = "Genetic relatedness to nearest curated genome",
      subtitle = "Transmission networks cannot be inferred from sequence data alone",
      x = NULL,
      y = "Value"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "none",
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )

  ggsave(out_path, p, width = 9, height = 6, dpi = 150)
}

# ---------------------------------------------------------------------------
# 5. Mutation - phenotype heatmap
# ---------------------------------------------------------------------------
plot_mutation_phenotype_heatmap <- function(obj, out_path) {
  variants <- obj$variants
  phenotypes <- obj$matched_phenotypes

  if (length(phenotypes) == 0) {
    p <- ggplot() + theme_void() +
      labs(title = "Mutation-phenotype heatmap", subtitle = "No phenotype associations detected")
  } else {
    df <- bind_rows(lapply(phenotypes, function(p) data.frame(
      variant = safe(p$genotype_description) %||% safe(p$variant_key) %||% "unknown",
      category = safe(p$phenotype_category) %||% "unspecified",
      impact = safe(p$phenotype_effect) %||% "known",
      stringsAsFactors = FALSE
    )))
    df$variant <- ifelse(df$variant == "" | is.na(df$variant), "unknown", df$variant)
    df$category <- ifelse(df$category == "" | is.na(df$category), "unspecified", df$category)

    p <- ggplot(df, aes(x = category, y = reorder(variant, category), fill = category)) +
      geom_tile(color = "white") +
      scale_fill_brewer(palette = "Set2") +
      labs(title = "Mutation-phenotype heatmap", x = "Phenotype category", y = "Variant") +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  }

  ggsave(out_path, p, width = 10, height = max(4, min(12, 1 + length(variants) * 0.4)), dpi = 150)
}

# ---------------------------------------------------------------------------
# 6. Evidence consistency summary
# ---------------------------------------------------------------------------
plot_evidence_consistency <- function(obj, out_path) {
  cons <- obj$analyses$evidence_consistency
  checks <- cons$metrics$consistency_checks
  if (is.null(checks) || length(checks) == 0) {
    png(out_path, width = 1000, height = 600, res = 120)
    plot.new()
    title("No consistency checks available")
    dev.off()
    return(invisible(NULL))
  }

  df <- bind_rows(lapply(checks, function(c) data.frame(
    check = safe(c$check) %||% "unknown",
    consistent = if (isTRUE(c$consistent)) "consistent" else if (isFALSE(c$consistent)) "conflicting" else "indeterminate",
    stringsAsFactors = FALSE
  )))

  df$consistent <- factor(df$consistent, levels = c("consistent", "indeterminate", "conflicting"))
  counts <- df %>% group_by(consistent) %>% summarise(n = n(), .groups = "drop")

  p <- ggplot(counts, aes(x = consistent, y = n, fill = consistent)) +
    geom_bar(stat = "identity", width = 0.6, color = "grey30") +
    geom_text(aes(label = n), vjust = -0.5, size = 5) +
    scale_fill_manual(values = c("consistent" = "#1a9641", "indeterminate" = "#fdae61", "conflicting" = "#d7191c")) +
    labs(
      title = "Evidence consistency across genomic, epidemiological and literature streams",
      subtitle = "Consistency of phylogeographic, temporal and phenotypic signals",
      x = NULL,
      y = "Number of checks"
    ) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "none")

  ggsave(out_path, p, width = 9, height = 6, dpi = 150)
}

# ---------------------------------------------------------------------------
# 7. Knowledge gaps
# ---------------------------------------------------------------------------
plot_knowledge_gaps <- function(obj, out_path) {
  gaps <- obj$analyses$knowledge_gaps$metrics$gaps
  if (is.null(gaps) || length(gaps) == 0) {
    png(out_path, width = 1000, height = 600, res = 120)
    plot.new()
    title("No knowledge gaps recorded")
    dev.off()
    return(invisible(NULL))
  }

  # Categorise gaps by keyword
  categories <- sapply(gaps, function(g) {
    if (grepl("tree|phylogen|clock|molecular-clock|genetic", g, ignore.case = TRUE)) "phylogenetic"
    else if (grepl("epi|outbreak|R0|serial|transmission param|exposure|travel", g, ignore.case = TRUE)) "epidemiological"
    else if (grepl("lineage|metadata|host|reservoir|country|date", g, ignore.case = TRUE)) "metadata"
    else if (grepl("phenotype|literature|reference|variant|catalogue", g, ignore.case = TRUE)) "phenotypic/literature"
    else "other"
  })

  df <- data.frame(category = categories, stringsAsFactors = FALSE) %>%
    group_by(category) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(category = reorder(category, n))

  p <- ggplot(df, aes(x = category, y = n, fill = category)) +
    geom_bar(stat = "identity", width = 0.6, color = "grey30") +
    geom_text(aes(label = n), hjust = -0.2, size = 4) +
    scale_fill_brewer(palette = "Set2") +
    coord_flip() +
    labs(
      title = "Knowledge gaps prioritised for additional data collection",
      subtitle = paste("Total gaps:", length(gaps)),
      x = NULL,
      y = "Count"
    ) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "none")

  ggsave(out_path, p, width = 10, height = max(4, 1 + nrow(df) * 0.5), dpi = 150)
}

# ---------------------------------------------------------------------------
# 8. Lineage persistence timeline
# ---------------------------------------------------------------------------
plot_lineage_persistence_timeline <- function(obj, out_path) {
  sample <- obj$sample
  lineage_id <- sample$lineage_metadata$lineage_id %||% NA

  # Use contextual tree tips only — no need for a whole-genome metadata export.
  tips <- obj$tree_tips
  tip_df <- bind_rows(lapply(tips, function(t) data.frame(
    collection_country = safe(t$country),
    collection_year = suppressWarnings(as.numeric(substr(safe(t$date), 1, 4))),
    stringsAsFactors = FALSE
  ))) %>% filter(!is.na(collection_year) & collection_country != "")
  df <- tip_df %>% group_by(collection_country, collection_year) %>% summarise(n = n(), .groups = "drop")

  if (nrow(df) == 0) {
    # No contextual dates available: just place the sample country in the current year.
    df <- data.frame(
      collection_country = sample$country,
      collection_year = suppressWarnings(as.numeric(substr(sample$collection_date, 1, 4))),
      n = 1,
      stringsAsFactors = FALSE
    ) %>% filter(!is.na(collection_year))
  }

  df$collection_country <- trimws(df$collection_country)
  df$sample <- ifelse(df$collection_country == sample$country, "sample country", "other")

  p <- ggplot(df, aes(x = collection_year, y = reorder(collection_country, collection_year))) +
    geom_line(aes(group = collection_country), color = "grey70", linewidth = 0.5) +
    geom_point(aes(size = n, color = sample), alpha = 0.8) +
    scale_color_manual(values = c("sample country" = "red", "other" = "steelblue")) +
    scale_size_continuous(range = c(2, 8)) +
    labs(title = "Lineage persistence timeline",
         subtitle = paste("First / last detection by country for lineage", lineage_id),
         x = "Year", y = "Country", size = "Genomes") +
    theme_minimal(base_size = 12)

  ggsave(out_path, p, width = 10, height = max(5, min(14, 0.5 * length(unique(df$collection_country)) + 2)), dpi = 150)
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
obj <- read_intel(INTEL_PATH)

message("Generating decision-oriented figures into ", OUT_DIR, " ...")

plot_time_scaled_phylogeny(obj, file.path(OUT_DIR, "01_time_scaled_phylogeny.png"))
plot_evidence_weighted_assessment(obj, file.path(OUT_DIR, "02_evidence_weighted_assessment.png"))
plot_source_attribution(obj, file.path(OUT_DIR, "03_source_attribution.png"))
plot_genetic_relatedness(obj, file.path(OUT_DIR, "04_genetic_relatedness.png"))
plot_mutation_phenotype_heatmap(obj, file.path(OUT_DIR, "05_mutation_phenotype_heatmap.png"))
plot_evidence_consistency(obj, file.path(OUT_DIR, "06_evidence_consistency.png"))
plot_knowledge_gaps(obj, file.path(OUT_DIR, "07_knowledge_gaps.png"))
plot_lineage_persistence_timeline(obj, file.path(OUT_DIR, "08_lineage_persistence_timeline.png"))

message("Done. Figures written to ", OUT_DIR)
