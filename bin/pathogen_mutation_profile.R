#!/usr/bin/env Rscript
#
# pathogen_mutation_profile.R
#
# Analysis #1: overall mutation profile (query vs background).
# Reads the exported knowledge warehouse DuckDB and writes TSVs only:
#   01_mutation_burden_samples.tsv    per-sample burden
#   01_mutation_burden_summary.tsv    group summary + statistical test
#   01_mutation_profile_manova.tsv    MANOVA on protein-level burden
#   01_plsda_scores.tsv               PLS-DA sample scores
#   01_plsda_loadings.tsv             PLS-DA protein loadings
#   01_plsda_vip.tsv                  PLS-DA protein VIP
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
  library(mixOmics)
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
  make_option(c("--translations-dir"), type = "character", default = NULL, help = "Directory containing per-protein aligned FASTA translations (from nextclade)"),
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
          m.ref_aa, m.alt_aa,
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
  dplyr::select(sample_id, sample_name, is_query, collection_year, country, host, outbreak, clade, lineage)

sample_order <- sample_meta$sample_id
all_samples <- tibble(sample_id = sample_order)

# ---------------------------------------------------------------------------
# Build binary sample x mutation matrix, including all samples
# ---------------------------------------------------------------------------
mutation_mat <- if (nrow(mutations_long) > 0) {
  mutations_long |>
    dplyr::select(sample_id, mutation_uid) |>
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
  sample_meta |> dplyr::select(sample_id, sample_name, is_query, collection_year, country, mutation_burden),
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
  dplyr::select(sample_id, is_query) |>
  inner_join(
    as.data.frame(protein_burden) |>
      tibble::rownames_to_column("sample_id") |>
      mutate(sample_id = as.integer(sample_id)),
    by = "sample_id"
  )

protein_names <- setdiff(colnames(manova_df), c("sample_id", "is_query"))

# Fallback row used whenever MANOVA cannot be computed. Keeps the same
# columns as the success path so downstream readers (dashboard formatRound
# on Wilks/statistic/num_Df/den_Df) never hit missing-column errors.
manova_skipped <- function(note) {
  tibble(
    term = "is_query", Wilks = NA_real_, statistic = NA_real_,
    num_Df = NA_real_, den_Df = NA_real_, p_value = NA_real_,
    note = note
  )
}

manova_out <- if (n_query >= 2 && n_bg >= 1 && length(protein_names) > 1) {
  f <- as.formula(paste("cbind(", paste(protein_names, collapse = ","), ") ~ is_query"))
  mfit <- tryCatch(manova(f, data = manova_df), error = function(e) NULL)
  # summary() can still fail when residual rank < number of responses
  # (e.g. very few samples relative to proteins), so guard it as well.
  s <- if (!is.null(mfit)) {
    tryCatch(summary(mfit, test = "Wilks"), error = function(e) NULL)
  }
  if (!is.null(s)) {
    as.data.frame(s$stats) |>
      tibble::rownames_to_column("term") |>
      filter(term != "Residuals") |>
      dplyr::select(term, Wilks = `Wilks`, statistic = `approx F`, num_Df = `num Df`, den_Df = `den Df`, p_value = `Pr(>F)`) |>
      mutate(note = "MANOVA on protein-level burden")
  } else {
    manova_skipped("MANOVA model failed")
  }
} else {
  note <- if (length(protein_names) <= 1) {
    "Only one or zero proteins; MANOVA skipped"
  } else if (n_query < 2) {
    "Not enough query samples for MANOVA (need >= 2)"
  } else {
    "No background samples for MANOVA"
  }
  manova_skipped(note)
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
      sample_meta |> dplyr::select(sample_id, sample_name, is_query),
      by = "sample_id"
    ) |>
    tidyr::pivot_longer(
      cols = -c(sample_id, sample_name, is_query),
      names_to = "protein_name",
      values_to = "mutation_positions"
    )
} else {
  tibble(
    sample_id = integer(),
    sample_name = character(),
    is_query = logical(),
    protein_name = character(),
    mutation_positions = integer()
  )
}

write_tsv(
  protein_burden_samples_long |>
    dplyr::select(sample_id, sample_name, is_query, protein_name, mutation_positions),
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
        broom::tidy(t.test(mutation_positions ~ is_query, data = d))
      } else {
        broom::tidy(wilcox.test(mutation_positions ~ is_query, data = d, exact = FALSE))
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
      mean_mutation_positions = mean(mutation_positions, na.rm = TRUE),
      median_mutation_positions = median(mutation_positions, na.rm = TRUE),
      sd_mutation_positions = sd(mutation_positions, na.rm = TRUE),
      .groups = "drop"
    ) |>
    left_join(protein_test, by = "protein_name")
} else {
  tibble(
    protein_name = character(),
    is_query = logical(),
    n = integer(),
    mean_mutation_positions = numeric(),
    median_mutation_positions = numeric(),
    sd_mutation_positions = numeric(),
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

# ---------------------------------------------------------------------------
# PLS-DA on protein burden to separate query vs background
# ---------------------------------------------------------------------------
log_info("Running PLS-DA on protein burden")

plsda_out <- if (n_query >= 2 && n_bg >= 2 && length(protein_names) >= 2) {
  X <- as.matrix(protein_burden)
  rownames(X) <- as.character(sample_meta$sample_id)
  Y <- factor(
    sample_meta$is_query,
    levels = c(FALSE, TRUE),
    labels = c("Background", "Query")
  )
  pls <- tryCatch(
    plsda(X, Y, ncomp = 2, scale = TRUE),
    error = function(e) NULL
  )

  if (!is.null(pls)) {
    scores <- as.data.frame(pls$variates$X) |>
      tibble::rownames_to_column("sample_id") |>
      mutate(sample_id = as.integer(sample_id)) |>
      inner_join(
        sample_meta |> dplyr::select(sample_id, sample_name, is_query, collection_year, country, outbreak, mutation_burden),
        by = "sample_id"
      ) |>
      dplyr::select(sample_id, sample_name, is_query, collection_year, country, outbreak, mutation_burden, PC1 = comp1, PC2 = comp2)

    loadings <- as.data.frame(pls$loadings$X) |>
      tibble::rownames_to_column("protein_name") |>
      dplyr::select(protein_name, PC1 = comp1, PC2 = comp2)

    vip <- loadings |>
      mutate(vip = abs(PC1) + abs(PC2)) |>
      dplyr::select(protein_name, vip)

    list(scores = scores, loadings = loadings, vip = vip)
  } else {
    list(
      scores = tibble(sample_id = integer(), sample_name = character(), is_query = logical(), collection_year = integer(), country = character(), outbreak = character(), mutation_burden = integer(), PC1 = numeric(), PC2 = numeric()),
      loadings = tibble(protein_name = character(), PC1 = numeric(), PC2 = numeric()),
      vip = tibble(protein_name = character(), vip = numeric())
    )
  }
} else {
  list(
    scores = tibble(sample_id = integer(), sample_name = character(), is_query = logical(), collection_year = integer(), country = character(), outbreak = character(), mutation_burden = integer(), PC1 = numeric(), PC2 = numeric()),
    loadings = tibble(protein_name = character(), PC1 = numeric(), PC2 = numeric()),
    vip = tibble(protein_name = character(), vip = numeric())
  )
}

write_tsv(plsda_out$scores, "01_plsda_scores.tsv", subdir = "mutation_profile")
write_mqc_tsv(
  plsda_out$scores,
  paste0(run_prefix, "_01_plsda_scores_mqc.tsv"),
  id = "plsda_scores",
  section_name = "PLS-DA Scores",
  description = "PLS-DA sample scores on protein-level mutation burden.",
  subdir = "mutation_profile/mqc"
)

write_tsv(plsda_out$loadings, "01_plsda_loadings.tsv", subdir = "mutation_profile")
write_mqc_tsv(
  plsda_out$loadings,
  paste0(run_prefix, "_01_plsda_loadings_mqc.tsv"),
  id = "plsda_loadings",
  section_name = "PLS-DA Loadings",
  description = "PLS-DA protein loadings on protein-level mutation burden.",
  subdir = "mutation_profile/mqc"
)

write_tsv(plsda_out$vip, "01_plsda_vip.tsv", subdir = "mutation_profile")
write_mqc_tsv(
  plsda_out$vip,
  paste0(run_prefix, "_01_plsda_vip_mqc.tsv"),
  id = "plsda_vip",
  section_name = "PLS-DA VIP",
  description = "PLS-DA variable importance (sum of absolute PC1 and PC2 loadings).",
  subdir = "mutation_profile/mqc"
)

# ---------------------------------------------------------------------------
# Position-specific amino acid frequencies from aligned translations
# ---------------------------------------------------------------------------
translations_dir <- opts$`translations-dir`
aa_frequencies <- NULL

if (!is.null(translations_dir) && dir.exists(translations_dir)) {
  log_info("Computing position-specific AA frequencies from translations")

  fasta_files <- list.files(translations_dir, pattern = "\\.fasta$", full.names = TRUE)

  if (length(fasta_files) > 0) {
    aa_freq_list <- lapply(fasta_files, function(fpath) {
      protein <- tools::file_path_sans_ext(basename(fpath))
      lines <- readLines(fpath, warn = FALSE)

      # Parse FASTA: extract sequences (skip headers)
      header_idx <- grep("^>", lines)
      if (length(header_idx) == 0) return(NULL)

      seqs <- vapply(seq_along(header_idx), function(i) {
        start <- header_idx[i] + 1L
        end <- if (i < length(header_idx)) header_idx[i + 1L] - 1L else length(lines)
        paste(lines[start:end], collapse = "")
      }, character(1))

      # All sequences should be the same length (aligned)
      seq_len_max <- max(nchar(seqs))
      if (seq_len_max == 0) return(NULL)

      # Split into character matrix
      char_mat <- do.call(rbind, strsplit(seqs, ""))
      n_seqs <- nrow(char_mat)
      n_pos <- ncol(char_mat)

      # Compute frequencies at each position
      pos_freq <- do.call(rbind, lapply(seq_len(n_pos), function(pos_i) {
        col <- char_mat[, pos_i]
        # Exclude gaps, unknown (X, -, *)
        valid <- col[!col %in% c("X", "-", "*", "x")]
        if (length(valid) == 0) return(NULL)
        tbl <- table(valid)
        tibble(
          protein_name = protein,
          position = pos_i,
          amino_acid = names(tbl),
          count = as.integer(tbl),
          total_valid = length(valid),
          frequency = round(as.numeric(tbl) / length(valid), 4)
        )
      }))

      pos_freq
    })

    aa_frequencies <- do.call(rbind, aa_freq_list)

    if (!is.null(aa_frequencies) && nrow(aa_frequencies) > 0) {
      # Mark the reference (most frequent) AA at each position
      aa_frequencies <- aa_frequencies |>
        group_by(protein_name, position) |>
        mutate(is_reference = frequency == max(frequency)) |>
        ungroup()

      write_tsv(aa_frequencies, "01_position_aa_frequencies.tsv", subdir = "mutation_profile")
      log_info("Wrote position AA frequencies: ", nrow(aa_frequencies), " rows across ",
               length(unique(aa_frequencies$protein_name)), " proteins")
    } else {
      log_info("No valid AA frequency data computed from translations")
      write_tsv(
        tibble(protein_name = character(), position = integer(), amino_acid = character(),
               count = integer(), total_valid = integer(), frequency = numeric(), is_reference = logical()),
        "01_position_aa_frequencies.tsv", subdir = "mutation_profile"
      )
    }
  } else {
    log_info("No FASTA files found in translations directory: ", translations_dir)
    write_tsv(
      tibble(protein_name = character(), position = integer(), amino_acid = character(),
             count = integer(), total_valid = integer(), frequency = numeric(), is_reference = logical()),
      "01_position_aa_frequencies.tsv", subdir = "mutation_profile"
    )
  }
} else {
  log_info("No translations directory provided or it does not exist; skipping AA frequency computation")
  write_tsv(
    tibble(protein_name = character(), position = integer(), amino_acid = character(),
           count = integer(), total_valid = integer(), frequency = numeric(), is_reference = logical()),
    "01_position_aa_frequencies.tsv", subdir = "mutation_profile"
  )
}

# ---------------------------------------------------------------------------
# Per-position mutation summary — drives the cleaner Mutation Landscape
# ---------------------------------------------------------------------------
compute_position_summary <- function(aa_freq, mut_long, n_total_samples) {
  # Count mutations per position from long mutation table
  mut_counts <- if (nrow(mut_long) > 0) {
    mut_long |>
      group_by(protein_name, position) |>
      summarise(
        reference_aa = dplyr::first(ref_aa, default = NA_character_),
        n_distinct_alt = dplyr::n_distinct(alt_aa[!is.na(alt_aa)]),
        n_mutated = dplyr::n_distinct(sample_id),
        n_query_mutated = dplyr::n_distinct(sample_id[is_query]),
        n_bg_mutated = dplyr::n_distinct(sample_id[!is_query]),
        .groups = "drop"
      )
  } else {
    tibble(
      protein_name = character(), position = integer(), reference_aa = character(),
      n_distinct_alt = integer(), n_mutated = integer(),
      n_query_mutated = integer(), n_bg_mutated = integer()
    )
  }

  if (!is.null(aa_freq) && nrow(aa_freq) > 0) {
    aa_valid <- aa_freq |>
      filter(!amino_acid %in% c("X", "-", "*", "x"))

    aa_summary <- aa_valid |>
      group_by(protein_name, position) |>
      summarise(
        reference_aa = amino_acid[is_reference][1],
        total_valid = max(total_valid, na.rm = TRUE),
        n_distinct_aa = sum(count > 0, na.rm = TRUE),
        reference_count = count[amino_acid == amino_acid[is_reference][1]][1],
        shannon_entropy = -sum(frequency[count > 0] * log2(frequency[count > 0]), na.rm = TRUE),
        .groups = "drop"
      ) |>
      mutate(
        reference_count = tidyr::replace_na(reference_count, 0L),
        n_mutated = total_valid - reference_count,
        mutated_fraction = if_else(total_valid > 0, n_mutated / total_valid, 0)
      )

    out <- aa_summary |>
      left_join(
        mut_counts |>
          dplyr::select(protein_name, position, n_query_mutated, n_bg_mutated),
        by = c("protein_name", "position")
      ) |>
      mutate(
        n_query_mutated = tidyr::replace_na(n_query_mutated, 0L),
        n_bg_mutated = tidyr::replace_na(n_bg_mutated, 0L)
      )
  } else {
    out <- mut_counts |>
      mutate(
        total_valid = n_total_samples,
        n_distinct_aa = n_distinct_alt + 1L,
        mutated_fraction = if_else(total_valid > 0, n_mutated / total_valid, 0),
        shannon_entropy = NA_real_
      )
  }

  out |>
    mutate(
      has_query_mutation = n_query_mutated > 0,
      shannon_entropy = tidyr::replace_na(shannon_entropy, 0)
    ) |>
    dplyr::select(
      protein_name, position, reference_aa, total_valid, n_distinct_aa,
      n_mutated, n_query_mutated, n_bg_mutated, mutated_fraction,
      shannon_entropy, has_query_mutation
    )
}

log_info("Computing per-position mutation summary")
position_summary <- compute_position_summary(aa_frequencies, mutations_long, nrow(samples))
write_tsv(position_summary, "01_position_summary.tsv", subdir = "mutation_profile")
log_info("Wrote position summary: ", nrow(position_summary), " rows")

# ---------------------------------------------------------------------------
# Per-sample mutation detail (long format) — drives the dashboard Mutation
# Landscape, catalogue, trajectory and related views.
# ---------------------------------------------------------------------------
mutation_detail <- if (nrow(mutations_long) > 0) {
  mutations_long |>
    dplyr::select(
      sample_id, sample_name, is_query,
      mutation_id, mutation_uid, mutation_label, mutation_type,
      protein_name, gene_name, position, ref_aa, alt_aa
    ) |>
    # Attach per-sample metadata used by the dashboard trajectory
    # (collection_year) and geo-map (country) views.
    left_join(
      sample_meta |> dplyr::select(sample_id, collection_year, country, host, outbreak, clade, lineage),
      by = "sample_id"
    ) |>
    arrange(protein_name, position, sample_id)
} else {
  tibble(
    sample_id = integer(), sample_name = character(), is_query = logical(),
    mutation_id = integer(), mutation_uid = character(), mutation_label = character(),
    mutation_type = character(), protein_name = character(), gene_name = character(),
    position = integer(), ref_aa = character(), alt_aa = character(),
    collection_year = numeric(), country = character(), host = character(),
    outbreak = character(), clade = character(), lineage = character()
  )
}
write_tsv(mutation_detail, "01_mutation_detail.tsv", subdir = "mutation_profile")

# ---------------------------------------------------------------------------
# Mutation phenotype annotations — drives the dashboard phenotype filtering
# and Mutation Intelligence Card phenotype panel.
# ---------------------------------------------------------------------------
log_info("Querying phenotype annotations")
mutation_phenotypes <- dbGetQuery(con,
  "SELECT DISTINCT
          m.mutation_id,
          m.mutation_label,
          m.ref_aa, m.position, m.alt_aa,
          COALESCE(g.gene_name, p.protein_name, 'unknown') AS protein_name,
          mp.phenotype, mp.effect, mp.evidence, mp.source
   FROM mutations m
   JOIN mutation_phenotypes mp ON m.mutation_id = mp.mutation_id
   LEFT JOIN proteins p ON m.protein_id = p.protein_id
   LEFT JOIN genes g ON m.gene_id = g.gene_id
   WHERE m.mutation_id IN (
         SELECT sm.mutation_id
         FROM sample_mutation sm
         JOIN samples s ON sm.sample_id = s.sample_id
         WHERE LOWER(s.species) = ? AND s.run_id = ?)",
  params = list(sp_lower, run_id)
)

if (nrow(mutation_phenotypes) > 0) {
  mutation_phenotypes <- mutation_phenotypes |>
    mutate(mutation_uid = paste0(protein_name, "_", mutation_label))
} else {
  mutation_phenotypes <- tibble(
    mutation_id = integer(), mutation_label = character(),
    ref_aa = character(), position = integer(), alt_aa = character(),
    protein_name = character(), phenotype = character(),
    effect = character(), evidence = character(), source = character(),
    mutation_uid = character()
  )
}
write_tsv(mutation_phenotypes, "01_mutation_phenotypes.tsv", subdir = "mutation_profile")

# ---------------------------------------------------------------------------
# Mutation functional / domain context — HMMER/Pfam domains, UniProt
# functions, GO terms and other protein-level annotations overlapping the
# mutation position or attached to the protein.
# ---------------------------------------------------------------------------
log_info("Querying mutation context")
mutation_context <- dbGetQuery(con,
  "WITH current_mutations AS (
     SELECT m.mutation_id, m.mutation_label, m.ref_aa, m.position, m.alt_aa,
            m.protein_id, m.gene_id
     FROM mutations m
     WHERE m.mutation_id IN (
       SELECT sm.mutation_id
       FROM sample_mutation sm
       JOIN samples s ON sm.sample_id = s.sample_id
       WHERE LOWER(s.species) = ? AND s.run_id = ?
     )
   )
   SELECT DISTINCT cm.mutation_id,
          cm.mutation_label,
          cm.ref_aa, cm.position, cm.alt_aa,
          COALESCE(g.gene_name, p.protein_name, 'unknown') AS protein_name,
          'domain' AS context_type,
          pd.domain_name AS context_label,
          CAST(NULL AS TEXT) AS context_description,
          pd.source,
          pd.env_start AS start_pos,
          pd.env_end AS end_pos,
          pd.evalue
   FROM current_mutations cm
   LEFT JOIN proteins p ON cm.protein_id = p.protein_id
   LEFT JOIN genes g ON cm.gene_id = g.gene_id
   JOIN protein_domains pd ON cm.protein_id = pd.protein_id
   WHERE cm.position BETWEEN pd.env_start AND pd.env_end
   UNION ALL
   SELECT DISTINCT cm.mutation_id, cm.mutation_label, cm.ref_aa, cm.position, cm.alt_aa,
          COALESCE(g.gene_name, p.protein_name, 'unknown'),
          'function', pf.source, pf.function_text, pf.source,
          CAST(NULL AS INTEGER) AS start_pos,
          CAST(NULL AS INTEGER) AS end_pos,
          CAST(NULL AS REAL) AS evalue
   FROM current_mutations cm
   LEFT JOIN proteins p ON cm.protein_id = p.protein_id
   LEFT JOIN genes g ON cm.gene_id = g.gene_id
   JOIN protein_functions pf ON cm.protein_id = pf.protein_id",
  params = list(sp_lower, run_id)
)

if (nrow(mutation_context) > 0) {
  mutation_context <- mutation_context |>
    mutate(mutation_uid = paste0(protein_name, "_", mutation_label))
} else {
  mutation_context <- tibble(
    mutation_id = integer(), mutation_label = character(),
    ref_aa = character(), position = integer(), alt_aa = character(),
    protein_name = character(), context_type = character(),
    context_label = character(), context_description = character(),
    source = character(), start_pos = integer(), end_pos = integer(),
    evalue = numeric(), mutation_uid = character()
  )
}
write_tsv(mutation_context, "01_mutation_context.tsv", subdir = "mutation_profile")

log_info("Done. Outputs in ", file.path(outdir, "mutation_profile"))
