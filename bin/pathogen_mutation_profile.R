#!/usr/bin/env Rscript
#
# pathogen_mutation_profile.R
#
# Analysis #1: overall mutation profile (query vs background).
# Reads the exported knowledge warehouse DuckDB and writes TSVs only:
#   01_mutation_burden_samples.tsv    per-sample burden
#   01_mutation_burden_summary.tsv    group summary + statistical test
#   01_mutation_profile_manova.tsv    MANOVA on protein-level burden
#
# Usage:
#   Rscript bin/pathogen_mutation_profile.R \
#       --duckdb results/knowledge_warehouse/knowledge_warehouse.duckdb \
#       --species bdbv \
#       --run-id <run_id> \
#       --outdir . \
#       [--parametric]

suppressPackageStartupMessages({
  library(optparse)
  library(DBI)
  library(duckdb)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(readr)
  library(broom)
})

# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------
write_tsv <- function(df, filename, subdir = NULL) {
  if (!is.data.frame(df)) df <- as.data.frame(df)
  if (nrow(df) == 0) df <- df[0, , drop = FALSE]
  out_path <- outdir
  if (!is.null(subdir)) {
    out_path <- file.path(outdir, subdir)
    dir.create(out_path, showWarnings = FALSE, recursive = TRUE)
  }
  write.table(df, file.path(out_path, filename), sep = "\t", row.names = FALSE, quote = FALSE, na = "")
}

write_mqc_tsv <- function(df, filename, id, section_name, description, subdir = NULL) {
  if (!is.data.frame(df)) df <- as.data.frame(df)
  if (nrow(df) == 0) df <- df[0, , drop = FALSE]
  out_path <- outdir
  if (!is.null(subdir)) {
    out_path <- file.path(outdir, subdir)
    dir.create(out_path, showWarnings = FALSE, recursive = TRUE)
  }
  full_path <- file.path(out_path, filename)
  con_out <- file(full_path, open = "w")
  writeLines(c(
    paste0("# id: '", id, "'"),
    paste0("# section_name: '", section_name, "'"),
    paste0("# description: '", description, "'"),
    "# plot_type: 'table'",
    "# pconfig:",
    paste0("#     id: '", id, "_table'"),
    paste0("#     namespace: '", section_name, "'")
  ), con_out)
  write.table(df, con_out, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
  close(con_out)
}

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
option_list <- list(
  make_option(c("--duckdb"), type = "character", default = NULL, help = "Path to exported knowledge warehouse .duckdb"),
  make_option(c("--species"), type = "character", default = NULL, help = "Species code (e.g. bdbv, sudv)"),
  make_option(c("--run-id"), type = "character", default = NULL, help = "Run prefix for output filenames"),
  make_option(c("--outdir"), type = "character", default = ".", help = "Output directory"),
  make_option(c("--parametric"), action = "store_true", default = FALSE, help = "Use parametric tests where appropriate")
)

opts <- parse_args(OptionParser(option_list = option_list))
stopifnot(!is.null(opts$duckdb), !is.null(opts$species))

outdir <- opts$outdir
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
run_id <- if (!is.null(opts$`run-id`)) opts$`run-id` else opts$species
run_prefix <- run_id

# ---------------------------------------------------------------------------
# Connect to DuckDB
# ---------------------------------------------------------------------------
con <- dbConnect(duckdb::duckdb(), dbdir = opts$duckdb, read_only = TRUE)
on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)

sp_lower <- tolower(opts$species)

# ---------------------------------------------------------------------------
# Load samples
# ---------------------------------------------------------------------------
samples <- dbGetQuery(con,
  "SELECT sample_id, sample_name, is_query, collection_date, country, host, outbreak, clade, lineage,
          EXTRACT(YEAR FROM collection_date) AS collection_year
   FROM samples
   WHERE LOWER(species) = ? AND run_id = ?",
  params = list(sp_lower, run_id)
)

if (nrow(samples) == 0) {
  message("No samples for species: ", opts$species)
  quit(status = 0)
}

samples$is_query <- as.logical(samples$is_query)

# ---------------------------------------------------------------------------
# Load mutations (may be empty for some species/samples)
# ---------------------------------------------------------------------------
mutations_long <- dbGetQuery(con,
  "SELECT s.sample_id, s.sample_name, s.is_query,
          m.mutation_id, m.mutation_label, m.position, m.mutation_type,
          COALESCE(p.protein_name, g.gene_name, 'unknown') AS protein_name,
          COALESCE(g.gene_name, 'unknown') AS gene_name
   FROM samples s
   JOIN sample_mutation sm ON s.sample_id = sm.sample_id
   JOIN mutations m ON sm.mutation_id = m.mutation_id
   LEFT JOIN proteins p ON m.protein_id = p.protein_id
   LEFT JOIN genes g ON m.gene_id = g.gene_id
   WHERE LOWER(s.species) = ? AND s.run_id = ?",
  params = list(sp_lower, run_id)
)

if (nrow(mutations_long) > 0) {
  mutations_long$is_query <- as.logical(mutations_long$is_query)
  mutations_long <- mutations_long |>
    # Some annotation sources (e.g. nextclade) create sample-specific protein
    # names like "SAMPLE-0001_GP". Use the stable gene name when available
    # so the same protein is grouped across samples.
    mutate(protein_name = ifelse(gene_name != "unknown", gene_name, protein_name)) |>
    mutate(mutation_uid = paste0(protein_name, "_", mutation_label)) |>
    distinct(sample_id, mutation_uid, .keep_all = TRUE)
}

# ---------------------------------------------------------------------------
# Build sample metadata and full sample list
# ---------------------------------------------------------------------------
sample_meta <- samples |>
  select(sample_id, sample_name, is_query, collection_year, country, host, outbreak, clade, lineage)

sample_order <- sample_meta$sample_id
all_samples <- tibble(sample_id = sample_order)

# ---------------------------------------------------------------------------
# Build binary sample x mutation matrix, including all samples
# ---------------------------------------------------------------------------
mutation_mat <- if (nrow(mutations_long) > 0) {
  mutations_long |>
    select(sample_id, mutation_uid) |>
    mutate(present = 1L) |>
    pivot_wider(
      id_cols = sample_id,
      names_from = mutation_uid,
      values_from = present,
      values_fill = 0L
    ) |>
    right_join(all_samples, by = "sample_id") |>
    mutate(across(-sample_id, ~ replace_na(., 0))) |>
    arrange(match(sample_id, sample_order)) |>
    column_to_rownames("sample_id") |>
    as.matrix()
} else {
  matrix(0L, nrow = length(sample_order), ncol = 0,
         dimnames = list(as.character(sample_order), NULL))
}

# ---------------------------------------------------------------------------
# Build protein burden matrix, including all samples
# ---------------------------------------------------------------------------
protein_burden <- if (nrow(mutations_long) > 0) {
  mutations_long |>
    count(sample_id, protein_name, name = "burden") |>
    pivot_wider(
      id_cols = sample_id,
      names_from = protein_name,
      values_from = burden,
      values_fill = 0L
    ) |>
    right_join(all_samples, by = "sample_id") |>
    mutate(across(-sample_id, ~ replace_na(., 0))) |>
    arrange(match(sample_id, sample_order)) |>
    column_to_rownames("sample_id") |>
    as.matrix()
} else {
  matrix(0L, nrow = length(sample_order), ncol = 0,
         dimnames = list(as.character(sample_order), NULL))
}

# ---------------------------------------------------------------------------
# Per-sample mutation burden
# ---------------------------------------------------------------------------
log_info <- function(...) message("[pg_mutation_profile] ", ...)
log_info("Computing mutation burden")

sample_meta <- sample_meta |>
  mutate(
    mutation_burden = as.integer(rowSums(mutation_mat[as.character(sample_id), , drop = TRUE]))
  )

n_query <- sum(sample_meta$is_query)
n_bg <- sum(!sample_meta$is_query)

# ---------------------------------------------------------------------------
# Per-sample burden TSV
# ---------------------------------------------------------------------------
write_tsv(
  sample_meta |> select(sample_id, sample_name, is_query, collection_year, country, mutation_burden),
  "01_mutation_burden_samples.tsv",
  subdir = "mutation_profile"
)

# ---------------------------------------------------------------------------
# Burden comparison: Wilcoxon (default) or t-test (--parametric)
# ---------------------------------------------------------------------------
burden_test <- if (n_query > 0 && n_bg > 0) {
  use_t <- opts$parametric && n_query >= 2 && n_bg >= 2
  if (use_t) {
    broom::tidy(t.test(mutation_burden ~ is_query, data = sample_meta))
  } else {
    broom::tidy(wilcox.test(mutation_burden ~ is_query, data = sample_meta, exact = FALSE))
  }
} else {
  tibble(statistic = NA_real_, p.value = NA_real_, estimate = NA_real_)
}

test_name <- if (n_query == 0 || n_bg == 0) {
  "single_group"
} else if (opts$parametric && n_query >= 2 && n_bg >= 2) {
  "t-test"
} else {
  "Wilcoxon rank-sum"
}

burden_summary <- sample_meta |>
  group_by(is_query) |>
  summarise(
    n = n(),
    mean_burden = mean(mutation_burden, na.rm = TRUE),
    median_burden = median(mutation_burden, na.rm = TRUE),
    sd_burden = sd(mutation_burden, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    test = test_name,
    statistic = if ("statistic" %in% names(burden_test)) burden_test$statistic else NA_real_,
    p_value = if ("p.value" %in% names(burden_test)) burden_test$p.value else NA_real_,
    estimate = if ("estimate" %in% names(burden_test)) burden_test$estimate else NA_real_
  )

write_tsv(burden_summary, "01_mutation_burden_summary.tsv", subdir = "mutation_profile")
write_mqc_tsv(
  burden_summary,
  paste0(run_prefix, "_01_mutation_burden_summary_mqc.tsv"),
  id = "mutation_burden_summary",
  section_name = "Mutation Burden Summary",
  description = "Overall mutation burden comparison between query and background samples.",
  subdir = "mutation_profile/mqc"
)

# ---------------------------------------------------------------------------
# MANOVA on protein burden matrix
# ---------------------------------------------------------------------------
log_info("Running MANOVA on protein burden")

manova_df <- sample_meta |>
  select(sample_id, is_query) |>
  inner_join(
    as.data.frame(protein_burden) |>
      tibble::rownames_to_column("sample_id") |>
      mutate(sample_id = as.integer(sample_id)),
    by = "sample_id"
  )

protein_names <- setdiff(colnames(manova_df), c("sample_id", "is_query"))

manova_out <- if (n_query >= 2 && n_bg >= 1 && length(protein_names) > 1) {
  f <- as.formula(paste("cbind(", paste(protein_names, collapse = ","), ") ~ is_query"))
  mfit <- tryCatch(manova(f, data = manova_df), error = function(e) NULL)
  if (!is.null(mfit)) {
    s <- summary(mfit, test = "Wilks")
    as.data.frame(s$stats) |>
      tibble::rownames_to_column("term") |>
      filter(term != "Residuals") |>
      select(term, Wilks = `Wilks`, statistic = `approx F`, num_Df = `num Df`, den_Df = `den Df`, p_value = `Pr(>F)`) |>
      mutate(note = "MANOVA on protein-level burden")
  } else {
    tibble(term = "is_query", note = "MANOVA model failed")
  }
} else {
  note <- if (length(protein_names) <= 1) {
    "Only one or zero proteins; MANOVA skipped"
  } else if (n_query < 2) {
    "Not enough query samples for MANOVA (need >= 2)"
  } else {
    "No background samples for MANOVA"
  }
  tibble(term = "is_query", note = note)
}

write_tsv(manova_out, "01_mutation_profile_manova.tsv", subdir = "mutation_profile")
write_mqc_tsv(
  manova_out,
  paste0(run_prefix, "_01_mutation_profile_manova_mqc.tsv"),
  id = "mutation_profile_manova",
  section_name = "MANOVA Protein Burden",
  description = "MANOVA test of protein-level mutation burden between query and background.",
  subdir = "mutation_profile/mqc"
)

# ---------------------------------------------------------------------------
# Per-sample per-protein burden (long format) and per-protein summary
# ---------------------------------------------------------------------------
log_info("Computing per-sample per-protein burden")

protein_burden_samples_long <- if (length(protein_names) > 0) {
  as.data.frame(protein_burden) |>
    tibble::rownames_to_column("sample_id") |>
    mutate(sample_id = as.integer(sample_id)) |>
    inner_join(
      sample_meta |> select(sample_id, sample_name, is_query),
      by = "sample_id"
    ) |>
    tidyr::pivot_longer(
      cols = -c(sample_id, sample_name, is_query),
      names_to = "protein_name",
      values_to = "mutation_count"
    )
} else {
  tibble(
    sample_id = integer(),
    sample_name = character(),
    is_query = logical(),
    protein_name = character(),
    mutation_count = integer()
  )
}

write_tsv(
  protein_burden_samples_long |>
    select(sample_id, sample_name, is_query, protein_name, mutation_count),
  "01_protein_burden_samples.tsv",
  subdir = "mutation_profile"
)
write_mqc_tsv(
  protein_burden_samples_long,
  paste0(run_prefix, "_01_protein_burden_samples_mqc.tsv"),
  id = "protein_burden_samples",
  section_name = "Protein Burden (per sample)",
  description = "Long-format per-sample per-protein mutation counts.",
  subdir = "mutation_profile/mqc"
)

log_info("Running per-protein query-vs-background tests")

protein_test <- if (length(protein_names) > 0) {
  lapply(protein_names, function(prot) {
    d <- protein_burden_samples_long |> filter(protein_name == prot)
    nq <- sum(d$is_query)
    nb <- sum(!d$is_query)
    if (nq > 0 && nb > 0) {
      use_t <- opts$parametric && nq >= 2 && nb >= 2
      res <- if (use_t) {
        broom::tidy(t.test(mutation_count ~ is_query, data = d))
      } else {
        broom::tidy(wilcox.test(mutation_count ~ is_query, data = d, exact = FALSE))
      }
      tibble(
        protein_name = prot,
        test = if (use_t) "t-test" else "Wilcoxon rank-sum",
        statistic = if ("statistic" %in% names(res)) res$statistic else NA_real_,
        p_value = if ("p.value" %in% names(res)) res$p.value else NA_real_,
        estimate = if ("estimate" %in% names(res)) res$estimate else NA_real_
      )
    } else {
      tibble(
        protein_name = prot,
        test = "single_group",
        statistic = NA_real_,
        p_value = NA_real_,
        estimate = NA_real_
      )
    }
  }) |>
    bind_rows()
} else {
  tibble(
    protein_name = character(),
    test = character(),
    statistic = numeric(),
    p_value = numeric(),
    estimate = numeric()
  )
}

protein_burden_summary <- if (nrow(protein_burden_samples_long) > 0) {
  protein_burden_samples_long |>
    group_by(protein_name, is_query) |>
    summarise(
      n = n(),
      mean_burden = mean(mutation_count, na.rm = TRUE),
      median_burden = median(mutation_count, na.rm = TRUE),
      sd_burden = sd(mutation_count, na.rm = TRUE),
      .groups = "drop"
    ) |>
    left_join(protein_test, by = "protein_name")
} else {
  tibble(
    protein_name = character(),
    is_query = logical(),
    n = integer(),
    mean_burden = numeric(),
    median_burden = numeric(),
    sd_burden = numeric(),
    test = character(),
    statistic = numeric(),
    p_value = numeric(),
    estimate = numeric()
  )
}

write_tsv(protein_burden_summary, "01_protein_burden_summary.tsv", subdir = "mutation_profile")
write_mqc_tsv(
  protein_burden_summary,
  paste0(run_prefix, "_01_protein_burden_summary_mqc.tsv"),
  id = "protein_burden_summary",
  section_name = "Protein Burden Summary",
  description = "Per-protein mutation burden comparison between query and background samples.",
  subdir = "mutation_profile/mqc"
)

log_info("Done. Outputs in ", file.path(outdir, "mutation_profile"))
