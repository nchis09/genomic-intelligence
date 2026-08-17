#!/usr/bin/env Rscript

# PATHOGEN_IDENTIFICATION — per-species statistical tables from a portable
# DuckDB export of the genomic-intelligence knowledge warehouse (see
# bin/export_knowledge_db.py).

suppressPackageStartupMessages({
  library(optparse)
  library(DBI)
  library(duckdb)
  library(dplyr)
  library(dbplyr)
  library(tibble)
  library(tidyr)
  library(stringr)
  library(jsonlite)
  library(ape)
  library(vegan)
  library(betapart)
  library(Biostrings)
})

# -----------------------------------------------------------------------------
# Output / statistical helpers
# -----------------------------------------------------------------------------
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

# Writes a MultiQC custom-content table: same tab-separated data as write_tsv(),
# prefixed with a '# id:'/'# plot_type:' comment-header block so MultiQC picks
# it up automatically. Using the same `id` across all species' files lets
# MultiQC merge per-species rows into a single combined section (see
# bin/build_knowledge_db.py::write_mqc_summary for the established convention).
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

safe_div <- function(x, y) ifelse(y == 0 | is.na(y), NA_real_, as.numeric(x) / as.numeric(y))

percentile <- function(x, ref) {
  if (length(ref) == 0) return(rep(NA_real_, length(x)))
  sapply(x, function(v) {
    if (is.na(v)) return(NA_real_)
    mean(ref <= v, na.rm = TRUE) * 100
  })
}

safe_fisher <- function(mat, min_expected = 5) {
  if (any(dim(mat) < 2) || any(is.na(mat))) {
    return(list(estimate = NA_real_, conf.int = c(NA_real_, NA_real_), p.value = NA_real_, method = "skipped"))
  }
  expected <- suppressWarnings(chisq.test(mat)$expected)
  if (any(expected < min_expected, na.rm = TRUE)) {
    return(list(estimate = NA_real_, conf.int = c(NA_real_, NA_real_), p.value = NA_real_, method = "expected_counts_too_low"))
  }
  ft <- suppressWarnings(fisher.test(mat))
  list(estimate = as.numeric(ft$estimate), conf.int = as.numeric(ft$conf.int), p.value = ft$p.value, method = ft$method)
}

safe_adonis2 <- function(dist_mat, group, n_query) {
  if (n_query < 3) {
    return(tibble(term = character(), df = numeric(), SumOfSqs = numeric(), MeanSqs = numeric(), R2 = numeric(), F = numeric(), Pr = numeric()))
  }
  df <- tibble(sample = rownames(dist_mat), group = group)
  dist_obj <- vegan::as.dist(dist_mat)
  mod <- suppressMessages(vegan::adonis2(dist_obj ~ group, data = df, permutations = 999, by = "margin"))
  out <- as_tibble(mod, rownames = "term")
  dplyr::rename(out, sum_sq = SumOfSqs, mean_sq = MeanSqs, f_statistic = F, p_value = Pr)
}

cophenetic_matrix <- function(newick) {
  phy <- newick_to_phylo(newick)
  if (is.null(phy)) return(NULL)
  as.matrix(ape::cophenetic.phylo(phy))
}

# -----------------------------------------------------------------------------
# Argument parsing
# -----------------------------------------------------------------------------
option_list <- list(
  make_option(c("--species"), type = "character", default = NULL,
              help = "Species to characterize (e.g. BDBV, EBOV)", metavar = "SPECIES"),
  make_option(c("--run-id"), type = "character", default = NULL,
              help = "analysis_runs.run_id", metavar = "RUN_ID"),
  make_option(c("--duckdb-file"), type = "character", default = NULL,
              help = "Path to the exported knowledge_warehouse.duckdb file", metavar = "FILE"),
  make_option(c("--outdir"), type = "character", default = ".",
              help = "Output directory", metavar = "OUTDIR"),
  make_option(c("--skip-plots"), action = "store_true", default = FALSE,
              help = "Skip PNG generation")
)

parser <- OptionParser(option_list = option_list)
opts <- parse_args(parser)

if (is.null(opts$species) || is.null(opts$`run-id`)) {
  stop("--species and --run-id are required")
}
if (is.null(opts$`duckdb-file`)) {
  stop("--duckdb-file is required")
}

outdir <- opts$outdir
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# Logging helpers
# -----------------------------------------------------------------------------
log_info <- function(...) message(format(Sys.time(), "%H:%M:%S"), " [INFO] ", paste(...))
log_warn <- function(...) message(format(Sys.time(), "%H:%M:%S"), " [WARN] ", paste(...))

# -----------------------------------------------------------------------------
# Database connection
# -----------------------------------------------------------------------------
con <- tryCatch(
  dbConnect(duckdb::duckdb(), dbdir = opts$`duckdb-file`, read_only = TRUE),
  error = function(e) {
    stop("Could not connect to knowledge warehouse: ", conditionMessage(e))
  }
)

on.exit(dbDisconnect(con), add = TRUE)

# -----------------------------------------------------------------------------
# Safe SQL helpers
# -----------------------------------------------------------------------------
get_rows <- function(con, sql, ..., species = opts$species, run_id = opts$`run-id`) {
  args <- list(...)
  if (grepl("\\?species\\b", sql)) args$species <- species
  if (grepl("\\?run_id\\b", sql)) args$run_id <- run_id
  q <- do.call(DBI::sqlInterpolate, c(list(con, sql), args))
  DBI::dbGetQuery(con, q)
}

has_rows <- function(df) is.data.frame(df) && nrow(df) > 0

# -----------------------------------------------------------------------------
# Data extraction layer
# -----------------------------------------------------------------------------
extract_inputs <- function() {
  inputs <- list()
  inputs$completeness <- list()

  # Tree(s)
  trees <- as.data.frame(get_rows(con, "
    SELECT tree_id, tree_method, newick
    FROM phylogenetic_trees
    WHERE species = ?species AND run_id = ?run_id
  "))
  if (nrow(trees) > 0) {
    trees$tree_id <- as.integer(trees$tree_id)
    trees$tree_method <- as.character(trees$tree_method)
    trees$newick <- as.character(trees$newick)
  }
  inputs$completeness$phylogenetic_trees <- list(nrow = nrow(trees), available_methods = unique(trees$tree_method))

  inputs$tree_iq <- trees %>%
    filter(tolower(tree_method) == "iqtree2") %>%
    dplyr::slice(1)

  inputs$tree_augur <- trees %>%
    filter(tolower(tree_method) %in% c("augur", "nextstrain", "auspice")) %>%
    dplyr::slice(1)

  inputs$tree_nextclade <- trees %>%
    filter(tolower(tree_method) == "nextclade") %>%
    dplyr::slice(1)

  # Tree tips
  if (has_rows(inputs$tree_iq)) {
    inputs$tip_meta <- get_rows(con, "
      SELECT tt.sample_id, tt.label, tt.is_query, tt.div,
             s.sample_name as sample, s.clade, s.lineage,
             s.country, s.host, s.collection_date, s.outbreak,
             s.nextclade_qc as sample_qc, s.genome_coverage,
             s.aa_mutation_count, s.nuc_substitution_count
      FROM tree_tips tt
      JOIN samples s ON s.sample_id = tt.sample_id
      WHERE tt.tree_id = ?tree_id
    ", tree_id = inputs$tree_iq$tree_id[[1]])
  } else {
    inputs$tip_meta <- tibble(
      sample_id = integer(), label = character(), is_query = logical(),
      div = numeric(), clade = character(), lineage = character(),
      sample = character(), country = character(), host = character(),
      collection_date = as.Date(character()), outbreak = character(),
      sample_qc = character(), genome_coverage = numeric(),
      aa_mutation_count = integer(), nuc_substitution_count = integer()
    )
  }
  inputs$completeness$tip_meta <- list(nrow = nrow(inputs$tip_meta))

  # Nextclade tree tips (best-match screening tree for species-assignment)
  if (has_rows(inputs$tree_nextclade)) {
    inputs$tip_meta_nextclade <- get_rows(con, "
      SELECT tt.sample_id, tt.label, tt.is_query, tt.div,
             s.sample_name as sample,
             COALESCE(s.clade, tt.clade) as clade,
             s.lineage,
             s.country, s.host, s.collection_date,
             COALESCE(s.outbreak, tt.outbreak) as outbreak,
             COALESCE(s.nextclade_qc, tt.nextclade_qc) as sample_qc,
             COALESCE(s.genome_coverage, tt.genome_coverage) as genome_coverage,
             COALESCE(s.aa_mutation_count, tt.aa_mutation_count) as aa_mutation_count,
             COALESCE(s.nuc_substitution_count, tt.nuc_mutation_count) as nuc_substitution_count
      FROM tree_tips tt
      LEFT JOIN samples s ON s.sample_id = tt.sample_id
      WHERE tt.tree_id = ?tree_id
    ", tree_id = inputs$tree_nextclade$tree_id[[1]])
  } else {
    inputs$tip_meta_nextclade <- tibble(
      sample_id = integer(), label = character(), is_query = logical(),
      div = numeric(), clade = character(), lineage = character(),
      sample = character(), country = character(), host = character(),
      collection_date = as.Date(character()), outbreak = character(),
      sample_qc = character(), genome_coverage = numeric(),
      aa_mutation_count = integer(), nuc_substitution_count = integer()
    )
  }
  inputs$completeness$tip_meta_nextclade <- list(nrow = nrow(inputs$tip_meta_nextclade))

  # Samples
  inputs$samples <- get_rows(con, "
    SELECT sample_id, sample_name as sample, pathogen, species, is_query,
           nextclade_qc, qc_score, genome_coverage, aa_mutation_count, nuc_substitution_count,
           alignment_score, divergence, best_dataset_file,
           nextclade_json::text as nextclade_json,
           country, host, collection_date, outbreak, clade, lineage
    FROM samples
    WHERE species = ?species AND run_id = ?run_id
  ")
  inputs$completeness$samples <- list(nrow = nrow(inputs$samples))

  # Nextstrain MSA fallback (used when no Nextclade screening tree exists for
  # this species, e.g. SUDV -- see register_nextstrain_msa in build_knowledge_db.py)
  msa_fallback_info <- get_rows(con, "
    SELECT file_path
    FROM pipeline_outputs
    WHERE run_id = ?run_id AND process_name = ?process_name
  ", process_name = paste0("nextstrain_msa_", tolower(opts$species)))
  inputs$nextstrain_msa_path <- if (nrow(msa_fallback_info) > 0) msa_fallback_info$file_path[1] else NA_character_

  msa_meta_info <- get_rows(con, "
    SELECT file_path
    FROM pipeline_outputs
    WHERE run_id = ?run_id AND process_name = ?process_name
  ", process_name = paste0("nextstrain_metadata_extended_", tolower(opts$species)))
  inputs$nextstrain_metadata_extended_path <- if (nrow(msa_meta_info) > 0) msa_meta_info$file_path[1] else NA_character_

  # Genes
  inputs$genes <- get_rows(con, "
    SELECT g.gene_id, g.gene_name, g.start_pos, g.end_pos, g.strand
    FROM genes g
    JOIN reference_genomes rg ON rg.ref_id = g.ref_id
    WHERE rg.species = ?species
  ")
  inputs$completeness$genes <- list(nrow = nrow(inputs$genes))

  # Mutations
  inputs$mutations <- get_rows(con, "
    SELECT m.mutation_id, m.mutation_label, m.protein_id, m.gene_id,
           g.gene_name as gene, m.ref_aa, m.position, m.alt_aa, m.mutation_type
    FROM mutations m
    LEFT JOIN genes g ON g.gene_id = m.gene_id
    LEFT JOIN reference_genomes rg ON rg.ref_id = g.ref_id
    WHERE rg.species = ?species
  ")
  inputs$completeness$mutations <- list(nrow = nrow(inputs$mutations))

  # Sample-mutation links
  if (has_rows(inputs$samples) && has_rows(inputs$mutations)) {
    inputs$sample_mutation <- get_rows(con, "
      SELECT sm.sample_id, sm.mutation_id
      FROM sample_mutation sm
      JOIN samples s ON s.sample_id = sm.sample_id
      WHERE s.species = ?species AND s.run_id = ?run_id
    ")
  } else {
    inputs$sample_mutation <- tibble(sample_id = integer(), mutation_id = integer())
  }
  inputs$completeness$sample_mutation <- list(nrow = nrow(inputs$sample_mutation))

  # Mutation phenotypes
  inputs$mutation_phenotypes <- get_rows(con, "
    SELECT mp.phenotype_id, mp.mutation_id, mp.phenotype, mp.effect, mp.evidence, mp.source
    FROM mutation_phenotypes mp
    JOIN mutations m ON m.mutation_id = mp.mutation_id
    LEFT JOIN genes g ON g.gene_id = m.gene_id
    LEFT JOIN reference_genomes rg ON rg.ref_id = g.ref_id
    WHERE rg.species = ?species
  ")
  inputs$completeness$mutation_phenotypes <- list(nrow = nrow(inputs$mutation_phenotypes))

  # Mutation community (tip x mutation) for quick checks
  if (has_rows(inputs$tip_meta)) {
    inputs$mutation_community <- get_rows(con, "
      SELECT tt.label, m.mutation_label
      FROM sample_mutation sm
      JOIN mutations m ON m.mutation_id = sm.mutation_id
      JOIN samples s ON s.sample_id = sm.sample_id
      JOIN tree_tips tt ON tt.sample_id = s.sample_id
      WHERE tt.tree_id = ?tree_id
        AND s.species = ?species
    ", tree_id = inputs$tree_iq$tree_id[[1]])
  } else {
    inputs$mutation_community <- tibble(label = character(), mutation_label = character())
  }
  inputs$completeness$mutation_community <- list(nrow = nrow(inputs$mutation_community))

  # Literature / evidence
  inputs$literature <- get_rows(con, "
    SELECT le.field as extraction_field, le.value as extraction_value, lp.pmid as pmid
    FROM literature_extractions le
    JOIN literature_papers lp ON lp.paper_id = le.paper_id
    WHERE lp.species = ?species
  ")
  inputs$completeness$literature <- list(nrow = nrow(inputs$literature))

  inputs
}

# -----------------------------------------------------------------------------
# Convert Newick to ape::phylo
# -----------------------------------------------------------------------------
newick_to_phylo <- function(newick) {
  if (is.null(newick) || is.na(newick) || nchar(newick) == 0) return(NULL)
  tryCatch(
    ape::read.tree(text = newick),
    error = function(e) {
      log_warn("Could not parse Newick: ", conditionMessage(e))
      NULL
    }
  )
}

# -----------------------------------------------------------------------------
# Helpers for Nextclade JSON extraction
# -----------------------------------------------------------------------------
get_json_field <- function(json_text, field) {
  if (is.na(json_text) || is.null(json_text) || nchar(as.character(json_text)) == 0) return(NA)
  x <- tryCatch(jsonlite::fromJSON(json_text), error = function(e) NULL)
  if (is.null(x)) return(NA)
  ux <- unlist(x, use.names = TRUE)
  if (field %in% names(ux)) return(ux[[field]])
  NA
}

flatten_nextclade_json <- function(json_text, sample_name) {
  if (is.na(json_text) || is.null(json_text) || nchar(as.character(json_text)) == 0) {
    return(tibble(sample = sample_name))
  }
  x <- tryCatch(jsonlite::fromJSON(json_text), error = function(e) NULL)
  if (is.null(x)) return(tibble(sample = sample_name))
  flat <- as_tibble(as.list(unlist(x, use.names = TRUE)))
  flat$sample <- sample_name
  flat
}

# -----------------------------------------------------------------------------
# MSA-based fallback nearest-reference (used when no Nextclade screening tree
# exists for this species, e.g. SUDV). Computes percent identity / distance
# directly from the pre-aligned nextstrain_ebola alignment.fasta, since a
# proper alignment makes a phylogenetic tree unnecessary for a simple nearest-
# neighbor lookup: columns are already in the same coordinate space.
# -----------------------------------------------------------------------------
compute_msa_fallback_nearest <- function(inputs, query_samples) {
  empty <- tibble(
    sample = character(), nearest_reference = character(),
    nearest_reference_taxon = character(), nearest_reference_outbreak = character(),
    msa_pct_identity = numeric(), msa_genetic_distance = numeric(),
    msa_mean_distance_to_refs = numeric()
  )
  if (is.na(inputs$nextstrain_msa_path) || length(query_samples) == 0) return(empty)
  if (!file.exists(inputs$nextstrain_msa_path)) return(empty)

  seqs <- tryCatch(Biostrings::readBStringSet(inputs$nextstrain_msa_path), error = function(e) NULL)
  if (is.null(seqs)) return(empty)
  seq_chars <- as.character(seqs)

  # Exclude ALL of this species' query samples (not just the ones we're
  # computing a fallback for) so another query sample in the same alignment
  # is never mistaken for a reference isolate.
  all_query_samples <- if (has_rows(inputs$samples)) {
    inputs$samples$sample[as.logical(inputs$samples$is_query)]
  } else {
    query_samples
  }
  ref_names <- setdiff(names(seq_chars), all_query_samples)
  if (length(ref_names) == 0) return(empty)

  meta <- NULL
  if (!is.na(inputs$nextstrain_metadata_extended_path) && file.exists(inputs$nextstrain_metadata_extended_path)) {
    meta <- tryCatch(
      as_tibble(read.delim(inputs$nextstrain_metadata_extended_path, sep = "\t", colClasses = "character", check.names = FALSE)),
      error = function(e) NULL
    )
  }

  ref_label <- function(accession) {
    label <- accession
    if (!is.null(meta) && accession %in% meta$accession) {
      row <- meta[meta$accession == accession, ][1, ]
      loc <- if ("location" %in% names(row) && !is.na(row$location) && nchar(row$location) > 0) row$location else "Not Provided"
      country <- if ("country" %in% names(row) && !is.na(row$country)) row$country else NA_character_
      date <- if ("date" %in% names(row) && !is.na(row$date)) row$date else NA_character_
      cdate <- paste0(coalesce(country, ""), ifelse(!is.na(date), paste0("/", date), ""))
      label <- paste0(accession, "|", loc, "|", cdate)
    }
    label
  }
  ref_field <- function(accession, field) {
    if (is.null(meta) || !(field %in% names(meta)) || !(accession %in% meta$accession)) return(NA_character_)
    val <- meta[[field]][meta$accession == accession][1]
    if (is.na(val) || nchar(as.character(val)) == 0) NA_character_ else as.character(val)
  }

  out <- empty
  for (q in query_samples) {
    if (!q %in% names(seq_chars)) next
    q_seq <- unlist(strsplit(seq_chars[[q]], ""))
    if (length(q_seq) == 0) next

    dists <- numeric(0)
    idents <- numeric(0)
    for (r in ref_names) {
      r_seq <- unlist(strsplit(seq_chars[[r]], ""))
      if (length(r_seq) != length(q_seq)) next
      valid <- q_seq != "-" & r_seq != "-"
      n_valid <- sum(valid)
      if (n_valid == 0) next
      ident <- sum(q_seq[valid] == r_seq[valid]) / n_valid
      idents[r] <- ident
      dists[r] <- 1 - ident
    }
    if (length(dists) == 0) next

    best <- names(dists)[which.min(dists)]
    out <- bind_rows(out, tibble(
      sample = q,
      nearest_reference = ref_label(best),
      nearest_reference_taxon = ref_field(best, "clade"),
      nearest_reference_outbreak = ref_field(best, "outbreak"),
      msa_pct_identity = as.numeric(idents[best]),
      msa_genetic_distance = as.numeric(dists[best]),
      msa_mean_distance_to_refs = as.numeric(mean(dists, na.rm = TRUE))
    ))
  }
  if (!has_rows(out)) return(empty)
  out
}

# -----------------------------------------------------------------------------
# 00 Species assignment
# -----------------------------------------------------------------------------
make_00_species_assignment <- function(inputs) {
  empty <- tibble(
    sample = character(), assigned_taxon = character(),
    genome_length = integer(), coverage = numeric(),
    qc_status = character(), qc_overall_score = numeric(), nt_identity_best_ref = numeric(),
    nt_identity_assigned_refs = numeric(), genetic_distance = numeric(),
    nearest_reference = character(), nearest_reference_taxon = character(),
    nearest_reference_outbreak = character(), msa_mean_distance_to_refs = numeric(),
    phylo_clade = character(), bootstrap = logical(),
    diagnostic_sites_supported = integer(), conflicting_sites = integer(),
    assignment = character()
  )
  if (!has_rows(inputs$samples)) return(empty)

  samples <- inputs$samples %>% filter(as.logical(is_query))
  if (!has_rows(samples)) return(empty)

  # Nearest reference from the best-match Nextclade tree
  nearest <- tibble(sample = character(),
                    nearest_reference = character(),
                    nearest_reference_taxon = character(),
                    nearest_reference_outbreak = character(),
                    msa_pct_identity = numeric(),
                    msa_genetic_distance = numeric(),
                    msa_mean_distance_to_refs = numeric())
  if (has_rows(inputs$tree_nextclade) && has_rows(inputs$tip_meta_nextclade)) {
    phy <- newick_to_phylo(inputs$tree_nextclade$newick[[1]])
    if (!is.null(phy)) {
      coph <- tryCatch(as.matrix(ape::cophenetic.phylo(phy)), error = function(e) NULL)
      if (!is.null(coph)) {
        refs <- inputs$tip_meta_nextclade %>% filter(!as.logical(is_query))
        query_tips <- inputs$tip_meta_nextclade %>% filter(as.logical(is_query))
        for (i in seq_len(nrow(query_tips))) {
          q_label <- as.character(query_tips$label[i])
          q_sample <- as.character(query_tips$sample[i])
          if (q_label %in% rownames(coph)) {
            sub <- coph[q_label, refs$label, drop = FALSE]
            vals <- sub[1, ]
            vals <- vals[!is.na(vals) & names(vals) != q_label]
            if (length(vals) > 0) {
              best <- names(vals)[which.min(vals)]
              taxon <- refs$clade[match(best, refs$label)]
              if (is.na(taxon) || is.null(taxon)) taxon <- NA_character_
              ref_outbreak <- refs$outbreak[match(best, refs$label)]
              if (is.na(ref_outbreak) || is.null(ref_outbreak)) ref_outbreak <- NA_character_
              nearest <- bind_rows(nearest, tibble(
                sample = q_sample,
                nearest_reference = as.character(best),
                nearest_reference_taxon = as.character(taxon),
                nearest_reference_outbreak = as.character(ref_outbreak),
                msa_pct_identity = NA_real_,
                msa_genetic_distance = NA_real_,
                msa_mean_distance_to_refs = NA_real_
              ))
            }
          }
        }
      }
    }
  }

  # Fall back to a direct MSA-based nearest-reference for any query sample
  # that didn't get a match from the Nextclade tree above (e.g. SUDV, which
  # has no usable Nextclade screening tree in this run).
  missing_samples <- setdiff(samples$sample, nearest$sample)
  if (length(missing_samples) > 0) {
    msa_fallback <- compute_msa_fallback_nearest(inputs, missing_samples)
    if (has_rows(msa_fallback)) nearest <- bind_rows(nearest, msa_fallback)
  }

  out <- samples %>%
    mutate(
      assigned_taxon = species,
      genome_length_raw = as.integer(sapply(nextclade_json, get_json_field, "genomeLength")),
      alignment_end = as.integer(sapply(nextclade_json, get_json_field, "alignmentEnd")),
      alignment_start = as.integer(sapply(nextclade_json, get_json_field, "alignmentStart")),
      total_substitutions = as.integer(sapply(nextclade_json, get_json_field, "totalSubstitutions")),
      genome_length = coalesce(genome_length_raw, alignment_end - alignment_start + 1),
      coverage = as.numeric(genome_coverage),
      qc_status = nextclade_qc,
      qc_overall_score = as.numeric(qc_score),
      nt_identity_best_ref = as.numeric(sapply(nextclade_json, get_json_field, "identity")),
      nt_identity_assigned_refs = NA_real_,
      derived_divergence = ifelse(
        !is.na(total_substitutions) & !is.na(alignment_end) & !is.na(alignment_start) & (alignment_end - alignment_start + 1) > 0,
        total_substitutions / (alignment_end - alignment_start + 1),
        NA_real_
      ),
      genetic_distance = coalesce(as.numeric(divergence), derived_divergence),
      phylo_clade = clade,
      bootstrap = NA,
      diagnostic_sites_supported = as.integer(sapply(nextclade_json, get_json_field, "privateNucMutations.totalPrivateSubstitutions")),
      conflicting_sites = NA_integer_
    ) %>%
    left_join(nearest, by = "sample") %>%
    mutate(
      nt_identity_best_ref = ifelse(
        is.na(nt_identity_best_ref) & !is.na(genetic_distance) & genetic_distance <= 1,
        1 - genetic_distance,
        nt_identity_best_ref
      ),
      # Prefer the MSA-fallback identity/distance (specific to the actual
      # nearest_reference chosen) over Nextclade's generic best-reference
      # estimate, when the fallback was used (e.g. SUDV).
      nt_identity_best_ref = coalesce(msa_pct_identity, nt_identity_best_ref),
      genetic_distance = coalesce(msa_genetic_distance, genetic_distance),
      assignment = case_when(
        is.na(qc_status) | is.na(coverage) ~ "AMBIGUOUS",
        qc_status == "good" & coverage >= 0.95 & !is.na(genetic_distance) & genetic_distance <= 50 ~ "CONFIRMED",
        qc_status == "good" & coverage >= 0.90 ~ "AMBIGUOUS",
        TRUE ~ "REJECTED"
      )
    ) %>%
    select(sample, assigned_taxon, genome_length, coverage, qc_status, qc_overall_score,
           nt_identity_best_ref, nt_identity_assigned_refs, genetic_distance,
           nearest_reference, nearest_reference_taxon, nearest_reference_outbreak,
           msa_mean_distance_to_refs,
           phylo_clade, bootstrap,
           diagnostic_sites_supported, conflicting_sites, assignment)

  # NOTE: unlike earlier versions, this returns ALL assignments (CONFIRMED,
  # AMBIGUOUS, REJECTED) rather than pre-filtering to CONFIRMED here. The
  # CONFIRMED-only filter for the user-facing summary table now happens in
  # make_identification_summary(); AMBIGUOUS/REJECTED rows are captured
  # separately by make_unresolved_samples() instead of being discarded.
  if (!has_rows(out)) return(empty)
  out
}

# -----------------------------------------------------------------------------
# 01 QC profile
# -----------------------------------------------------------------------------
make_01_qc_profile <- function(inputs) {
  empty <- tibble(
    sample = character(), coverage = numeric(), cds_coverage = numeric(),
    substitutions = integer(), deletions = integer(), insertions = integer(),
    frame_shifts = integer(), missing_bases = integer(), non_acgtns = integer(),
    amino_acid_substitutions = integer(), alignment_score = numeric(),
    private_nucleotide_mutations = integer(), snp_clusters = integer(),
    stop_codons = integer(), qc_status = character()
  )
  if (!has_rows(inputs$samples)) return(empty)

  samples <- inputs$samples %>% filter(as.logical(is_query))
  if (!has_rows(samples)) return(empty)

  flat <- bind_rows(lapply(seq_len(nrow(samples)), function(i) {
    flatten_nextclade_json(samples$nextclade_json[[i]], samples$sample[[i]])
  }))

  get_col <- function(df, col, default = NA) {
    if (col %in% names(df)) df[[col]] else rep(default, nrow(df))
  }

  out <- tibble(
    sample = as.character(get_col(flat, "sample", "")),
    coverage = as.numeric(get_col(flat, "coverage")),
    cds_coverage = NA_real_,
    substitutions = as.integer(get_col(flat, "totalSubstitutions")),
    deletions = as.integer(get_col(flat, "totalDeletions")),
    insertions = as.integer(get_col(flat, "totalInsertions")),
    frame_shifts = as.integer(get_col(flat, "totalFrameShifts")),
    missing_bases = as.integer(get_col(flat, "totalMissing")),
    non_acgtns = as.integer(get_col(flat, "totalNonACGTNs")),
    amino_acid_substitutions = as.integer(get_col(flat, "totalAminoacidSubstitutions")),
    alignment_score = as.numeric(get_col(flat, "alignmentScore")),
    private_nucleotide_mutations = as.integer(get_col(flat, "privateNucMutations.totalPrivateSubstitutions")),
    snp_clusters = NA_integer_,
    stop_codons = NA_integer_,
    qc_status = as.character(get_col(flat, "qc.overallStatus"))
  )
  out
}

# -----------------------------------------------------------------------------
# 02 Sequence similarity
# -----------------------------------------------------------------------------
make_02_sequence_similarity <- function(inputs) {
  empty <- tibble(
    sample = character(), reference = character(),
    percent_nucleotide_identity = numeric(), number_of_mismatches = integer(),
    number_of_gaps = integer(), alignment_coverage = numeric(),
    p_distance = numeric()
  )
  if (!has_rows(inputs$samples)) return(empty)

  # Query the registered aligned FASTA for the best-match species
  fasta_info <- get_rows(con, "
    SELECT file_path
    FROM pipeline_outputs
    WHERE run_id = ?run_id
      AND file_name = ?file_name
      AND file_type = 'fasta'
  ", file_name = paste0("input_FASTA_", tolower(opts$species), ".aligned.fasta"))

  if (nrow(fasta_info) == 0 || is.na(fasta_info$file_path[1])) return(empty)

  fasta_path <- fasta_info$file_path[1]
  if (!file.exists(fasta_path)) {
    log_warn("Aligned FASTA not found: ", fasta_path)
    return(empty)
  }

  seqs <- tryCatch(Biostrings::readBStringSet(fasta_path), error = function(e) {
    log_warn("Could not read aligned FASTA: ", conditionMessage(e))
    NULL
  })
  if (is.null(seqs)) return(empty)

  seq_chars <- as.character(seqs)
  names(seq_chars) <- names(seqs)

  queries <- inputs$samples$sample[as.logical(inputs$samples$is_query)]
  refs <- inputs$tip_meta_nextclade %>%
    filter(!as.logical(is_query), label %in% names(seq_chars))

  if (length(queries) == 0) return(empty)

  # If Nextclade's aligned FASTA has no reference sequences, report one summary row per query from the TSV JSON
  if (nrow(refs) == 0) {
    out <- tibble(
      sample = character(), reference = character(),
      percent_nucleotide_identity = numeric(), number_of_mismatches = integer(),
      number_of_gaps = integer(), alignment_coverage = numeric(),
      p_distance = numeric()
    )
    for (q in queries) {
      srow <- inputs$samples[inputs$samples$sample == q, ]
      if (nrow(srow) != 1) next
      js <- srow$nextclade_json[1]
      identity <- as.numeric(get_json_field(js, "identity"))
      total_subs <- as.integer(get_json_field(js, "totalSubstitutions"))
      total_del <- as.integer(get_json_field(js, "totalDeletions"))
      total_ins <- as.integer(get_json_field(js, "totalInsertions"))
      total_miss <- as.integer(get_json_field(js, "totalMissing"))
      aln_start <- as.integer(get_json_field(js, "alignmentStart"))
      aln_end <- as.integer(get_json_field(js, "alignmentEnd"))
      genome_len <- as.integer(get_json_field(js, "genomeLength"))
      if (!is.na(aln_end) && !is.na(aln_start) && aln_end >= aln_start) {
        aln_len <- aln_end - aln_start + 1
      } else {
        aln_len <- genome_len
      }
      if (is.na(aln_len) || aln_len <= 0) next
      p_dist <- if (!is.na(total_subs)) total_subs / aln_len else NA_real_
      if (is.na(identity) && !is.na(p_dist) && p_dist <= 1) identity <- 1 - p_dist
      gaps <- as.integer(sum(c(total_del, total_ins, total_miss), na.rm = TRUE))
      coverage <- if (!is.na(total_miss) && aln_len > 0) (aln_len - total_miss) / aln_len else NA_real_
      out <- bind_rows(out, tibble(
        sample = q, reference = "best_reference",
        percent_nucleotide_identity = as.numeric(identity),
        number_of_mismatches = as.integer(total_subs),
        number_of_gaps = gaps,
        alignment_coverage = as.numeric(coverage),
        p_distance = as.numeric(p_dist)
      ))
    }
    if (!has_rows(out)) return(empty)
    return(out)
  }

  out <- tibble(
    sample = character(), reference = character(),
    percent_nucleotide_identity = numeric(), number_of_mismatches = integer(),
    number_of_gaps = integer(), alignment_coverage = numeric(),
    p_distance = numeric()
  )

  for (q in queries) {
    if (!q %in% names(seq_chars)) next
    q_ungapped <- gsub("-", "", seq_chars[q])
    if (nchar(q_ungapped) == 0) next
    q_dna <- tryCatch(Biostrings::DNAString(q_ungapped), error = function(e) NULL)
    if (is.null(q_dna)) next

    for (r in refs$label) {
      r_ungapped <- gsub("-", "", seq_chars[r])
      if (nchar(r_ungapped) == 0) next
      r_dna <- tryCatch(Biostrings::DNAString(r_ungapped), error = function(e) NULL)
      if (is.null(r_dna)) next

      aln <- tryCatch(Biostrings::pairwiseAlignment(
        q_dna, r_dna,
        type = "global",
        substitutionMatrix = Biostrings::nucleotideSubstitutionMatrix(match = 1, mismatch = -1, baseOnly = FALSE),
        gapOpening = 5, gapExtension = 2
      ), error = function(e) NULL)
      if (is.null(aln)) next

      ap <- Biostrings::alignedPattern(aln)
      nmatch <- Biostrings::nmatch(aln)
      nmismatch <- Biostrings::nmismatch(aln)
      aln_len <- nchar(as.character(ap))
      n_gaps <- aln_len - nmatch - nmismatch
      p_dist <- if (aln_len > 0) Biostrings::nedit(aln) / aln_len else NA_real_
      identity <- if (aln_len > 0) nmatch / aln_len else NA_real_
      coverage <- if (nchar(as.character(r_dna)) > 0) nchar(as.character(q_dna)) / nchar(as.character(r_dna)) else NA_real_

      out <- bind_rows(out, tibble(
        sample = q,
        reference = r,
        percent_nucleotide_identity = as.numeric(identity),
        number_of_mismatches = as.integer(nmismatch),
        number_of_gaps = as.integer(n_gaps),
        alignment_coverage = as.numeric(coverage),
        p_distance = as.numeric(p_dist)
      ))
    }
  }

  if (!has_rows(out)) return(empty)
  out
}

# -----------------------------------------------------------------------------
# 03 MSA profile
# -----------------------------------------------------------------------------
make_03_msa_profile <- function(inputs) {
  empty <- tibble(
    sample = character(), msa_method = character(),
    n_seqs = integer(), alignment_length = integer(),
    mean_identity = numeric()
  )
  if (!has_rows(inputs$samples)) return(empty)

  fasta_info <- get_rows(con, "
    SELECT file_path
    FROM pipeline_outputs
    WHERE run_id = ?run_id
      AND file_name = ?file_name
      AND file_type = 'fasta'
  ", file_name = paste0("input_FASTA_", tolower(opts$species), ".aligned.fasta"))

  if (nrow(fasta_info) == 0 || is.na(fasta_info$file_path[1])) return(empty)
  fasta_path <- fasta_info$file_path[1]
  if (!file.exists(fasta_path)) {
    log_warn("Aligned FASTA not found: ", fasta_path)
    return(empty)
  }

  seqs <- tryCatch(Biostrings::readBStringSet(fasta_path), error = function(e) {
    log_warn("Could not read aligned FASTA: ", conditionMessage(e))
    NULL
  })
  if (is.null(seqs)) return(empty)

  seq_chars <- as.character(seqs)
  n_seqs <- length(seq_chars)
  if (n_seqs == 0) return(empty)
  aln_len <- nchar(seq_chars[1])

  queries <- inputs$samples$sample[as.logical(inputs$samples$is_query)]
  if (length(queries) == 0) return(empty)

  out <- empty
  for (q in queries) {
    q_seq_val <- seq_chars[q]
    if (is.na(q_seq_val)) next
    q_seq <- unlist(strsplit(q_seq_val, ""))
    if (length(q_seq) == 0) next
    ids <- numeric()
    for (nm in names(seq_chars)) {
      if (nm == q) next
      o_seq_val <- seq_chars[nm]
      if (is.na(o_seq_val)) next
      o_seq <- unlist(strsplit(o_seq_val, ""))
      if (length(o_seq) != length(q_seq)) next
      valid <- q_seq != "-" & o_seq != "-"
      if (sum(valid) == 0) next
      ids <- c(ids, sum(q_seq[valid] == o_seq[valid]) / sum(valid))
    }
    out <- bind_rows(out, tibble(
      sample = q,
      msa_method = "nextclade",
      n_seqs = as.integer(n_seqs),
      alignment_length = as.integer(aln_len),
      mean_identity = if (length(ids) > 0) as.numeric(mean(ids, na.rm = TRUE)) else NA_real_
    ))
  }
  if (!has_rows(out)) return(empty)
  out
}

# -----------------------------------------------------------------------------
# 04 Phylogenetic placement
# -----------------------------------------------------------------------------
make_04_phylogenetic_placement <- function(inputs) {
  empty <- tibble(
    sample = character(),
    nearest_1 = character(), distance_to_nearest_1 = numeric(),
    nearest_2 = character(), distance_to_nearest_2 = numeric(),
    nearest_3 = character(), distance_to_nearest_3 = numeric(),
    nearest_4 = character(), distance_to_nearest_4 = numeric(),
    nearest_5 = character(), distance_to_nearest_5 = numeric(),
    distance_to_assigned_clade = numeric(),
    distance_to_other_species = numeric()
  )
  if (!has_rows(inputs$samples) || !has_rows(inputs$tree_nextclade) || !has_rows(inputs$tip_meta_nextclade)) {
    return(empty)
  }

  phy <- newick_to_phylo(inputs$tree_nextclade$newick[[1]])
  if (is.null(phy)) return(empty)

  coph <- tryCatch(as.matrix(ape::cophenetic.phylo(phy)), error = function(e) NULL)
  if (is.null(coph)) return(empty)

  # Only keep query tips that belong to the current species (best-match)
  query_tips <- inputs$tip_meta_nextclade %>%
    filter(as.logical(is_query), sample %in% inputs$samples$sample)
  queries_label <- as.character(query_tips$label)
  queries_sample <- as.character(query_tips$sample)
  refs <- inputs$tip_meta_nextclade %>% filter(!as.logical(is_query))
  ref_labels <- as.character(refs$label[refs$label %in% rownames(coph)])

  if (length(queries_label) == 0 || length(ref_labels) == 0) return(empty)

  out <- empty
  for (i in seq_along(queries_label)) {
    q_label <- queries_label[i]
    q_sample <- queries_sample[i]
    if (!q_label %in% rownames(coph)) next
    dists <- coph[q_label, ref_labels]
    dists <- dists[!is.na(dists)]
    if (length(dists) == 0) next

    ord <- order(dists, na.last = TRUE)
    top5 <- ref_labels[ord][1:min(5, length(ord))]
    top5d <- dists[ord][1:min(5, length(ord))]

    row <- tibble(
      sample = q_sample,
      nearest_1 = as.character(top5[1]),
      distance_to_nearest_1 = as.numeric(top5d[1]),
      nearest_2 = as.character(ifelse(length(top5) >= 2, top5[2], NA)),
      distance_to_nearest_2 = as.numeric(ifelse(length(top5) >= 2, top5d[2], NA)),
      nearest_3 = as.character(ifelse(length(top5) >= 3, top5[3], NA)),
      distance_to_nearest_3 = as.numeric(ifelse(length(top5) >= 3, top5d[3], NA)),
      nearest_4 = as.character(ifelse(length(top5) >= 4, top5[4], NA)),
      distance_to_nearest_4 = as.numeric(ifelse(length(top5) >= 4, top5d[4], NA)),
      nearest_5 = as.character(ifelse(length(top5) >= 5, top5[5], NA)),
      distance_to_nearest_5 = as.numeric(ifelse(length(top5) >= 5, top5d[5], NA)),
      distance_to_assigned_clade = as.numeric(mean(dists, na.rm = TRUE)),
      distance_to_other_species = NA_real_
    )
    out <- bind_rows(out, row)
  }

  if (!has_rows(out)) return(empty)
  out
}

# -----------------------------------------------------------------------------
# Unified identification summary
#
# One row per CONFIRMED sample, combining the species assignment, the
# similarity/genetic-distance to that sample's own closest reference
# (joined from 02_sequence_similarity on the *specific* nearest_reference
# tip, not the generic Nextclade dataset identity), MSA identity, and
# distance to the assigned phylogenetic clade. Replaces the previous 4
# separate tables with a single dashboard-facing table.
# -----------------------------------------------------------------------------
make_identification_summary <- function(species_df, seq_sim_df, msa_df, phylo_df) {
  empty <- tibble(
    sample = character(), species_identified = character(), coverage = numeric(),
    qc_overall_score = numeric(),
    closest_reference = character(), pct_identity_closest_ref = numeric(),
    genetic_distance_closest_ref = numeric(),
    distance_to_assigned_clade = numeric(), msa_identity = numeric(),
    closest_outbreak = character()
  )
  if (!has_rows(species_df)) return(empty)

  confirmed <- species_df %>% filter(assignment == "CONFIRMED")
  if (!has_rows(confirmed)) return(empty)

  seq_sim_closest <- tibble(sample = character(), percent_nucleotide_identity = numeric(), p_distance = numeric())
  if (has_rows(seq_sim_df)) {
    seq_sim_closest <- seq_sim_df %>%
      inner_join(confirmed %>% select(sample, nearest_reference), by = "sample") %>%
      filter(reference == nearest_reference) %>%
      select(sample, percent_nucleotide_identity, p_distance) %>%
      distinct(sample, .keep_all = TRUE)
  }

  msa_identity <- tibble(sample = character(), mean_identity = numeric())
  if (has_rows(msa_df)) {
    msa_identity <- msa_df %>% select(sample, mean_identity) %>% distinct(sample, .keep_all = TRUE)
  }

  clade_distance <- tibble(sample = character(), distance_to_assigned_clade = numeric())
  if (has_rows(phylo_df)) {
    clade_distance <- phylo_df %>% select(sample, distance_to_assigned_clade) %>% distinct(sample, .keep_all = TRUE)
  }

  out <- confirmed %>%
    left_join(seq_sim_closest, by = "sample") %>%
    left_join(msa_identity, by = "sample") %>%
    left_join(clade_distance, by = "sample") %>%
    mutate(
      pct_identity_closest_ref = coalesce(percent_nucleotide_identity, nt_identity_best_ref),
      genetic_distance_closest_ref = coalesce(p_distance, genetic_distance),
      # MSA-fallback mean distance to all references (used when there's no
      # Nextclade tree to compute a proper clade-distance from, e.g. SUDV).
      distance_to_assigned_clade = coalesce(distance_to_assigned_clade, msa_mean_distance_to_refs),
      mean_identity = coalesce(mean_identity, nt_identity_best_ref)
    ) %>%
    select(
      sample,
      species_identified = assigned_taxon,
      coverage,
      qc_overall_score,
      closest_reference = nearest_reference,
      pct_identity_closest_ref,
      genetic_distance_closest_ref,
      distance_to_assigned_clade,
      msa_identity = mean_identity,
      closest_outbreak = nearest_reference_outbreak
    )

  if (!has_rows(out)) return(empty)
  out
}

# -----------------------------------------------------------------------------
# Unresolved samples log
#
# Samples that did NOT pass identification (AMBIGUOUS/REJECTED) for this
# species, kept for traceability instead of silently disappearing from the
# per-species output.
# -----------------------------------------------------------------------------
make_unresolved_samples <- function(species_df) {
  empty <- tibble(
    sample = character(), assignment = character(),
    coverage = numeric(), qc_status = character()
  )
  if (!has_rows(species_df)) return(empty)

  out <- species_df %>%
    filter(assignment %in% c("AMBIGUOUS", "REJECTED")) %>%
    select(sample, assignment, coverage, qc_status)

  if (!has_rows(out)) return(empty)
  out
}

# -----------------------------------------------------------------------------
# 01 Genetic distinctiveness
# -----------------------------------------------------------------------------
make_01_genetic_distinctiveness <- function(inputs) {
  if (!has_rows(inputs$tip_meta)) {
    return(tibble(sample = character(), is_query = logical(), div = numeric(),
                  mean_pairwise_distance = numeric(), nearest_background_distance = numeric(),
                  nearest_background_sample = character(), nearest_background_percentile = numeric(),
                  long_branch_score = numeric(), long_branch_z = numeric()))
  }
  tip_meta <- inputs$tip_meta %>%
    mutate(label = as.character(label), is_query = as.logical(is_query))
  cop <- cophenetic_matrix(if (has_rows(inputs$tree_iq)) inputs$tree_iq$newick[[1]] else NULL)
  if (is.null(cop)) {
    return(tip_meta %>% select(sample, is_query, div) %>%
             mutate(mean_pairwise_distance = NA_real_, nearest_background_distance = NA_real_,
                    nearest_background_sample = NA_character_, nearest_background_percentile = NA_real_,
                    long_branch_score = NA_real_, long_branch_z = NA_real_))
  }
  labels <- rownames(cop)
  tip_meta <- tip_meta %>% filter(label %in% labels)
  mean_dist <- rowMeans(cop, na.rm = TRUE)
  bg_labels <- tip_meta %>% filter(!is_query) %>% pull(label)

  nearest_bg <- sapply(tip_meta$label, function(lbl) {
    if (length(bg_labels) == 0) return(NA_real_)
    min(cop[lbl, bg_labels], na.rm = TRUE)
  })
  nearest_bg_sample <- sapply(tip_meta$label, function(lbl) {
    if (length(bg_labels) == 0) return(NA_character_)
    bg_dists <- cop[lbl, bg_labels, drop = TRUE]
    names(bg_dists)[which.min(bg_dists)]
  })
  bg_bg_dists <- if (length(bg_labels) >= 2) as.vector(cop[bg_labels, bg_labels]) else NA_real_
  nearest_bg_pct <- percentile(nearest_bg, bg_bg_dists)

  n <- ncol(cop)
  lb <- sapply(tip_meta$label, function(lbl) {
    sqrt(sum((cop[lbl, labels] - mean_dist)^2, na.rm = TRUE) / max(1, n - 1))
  })
  lb_z <- (lb - mean(lb, na.rm = TRUE)) / sd(lb, na.rm = TRUE)

  tibble(
    sample = tip_meta$sample,
    is_query = tip_meta$is_query,
    div = tip_meta$div,
    mean_pairwise_distance = mean_dist[tip_meta$label],
    nearest_background_distance = nearest_bg[tip_meta$label],
    nearest_background_sample = nearest_bg_sample[tip_meta$label],
    nearest_background_percentile = nearest_bg_pct[tip_meta$label],
    long_branch_score = lb[tip_meta$label],
    long_branch_z = lb_z[tip_meta$label]
  )
}

# -----------------------------------------------------------------------------
# 02 Divergence distribution (branch lengths)
# -----------------------------------------------------------------------------
make_02_divergence_distribution <- function(inputs) {
  empty <- tibble(sample = character(), is_query = logical(), terminal_branch_length = numeric(),
                  terminal_branch_length_percentile = numeric(), distance_to_nearest_background = numeric(),
                  distance_to_nearest_background_percentile = numeric(),
                  distance_to_nearest_background_sample = character(),
                  distance_to_mrca_with_nearest_background = numeric())
  if (!has_rows(inputs$tip_meta)) return(empty)
  tip_meta <- inputs$tip_meta %>%
    mutate(label = as.character(label), is_query = as.logical(is_query))
  phy <- newick_to_phylo(if (has_rows(inputs$tree_iq)) inputs$tree_iq$newick[[1]] else NULL)
  if (is.null(phy)) {
    return(tip_meta %>% select(sample, is_query) %>%
             mutate(terminal_branch_length = NA_real_, terminal_branch_length_percentile = NA_real_,
                    distance_to_nearest_background = NA_real_, distance_to_nearest_background_percentile = NA_real_,
                    distance_to_nearest_background_sample = NA_character_,
                    distance_to_mrca_with_nearest_background = NA_real_))
  }
  cop <- cophenetic_matrix(if (has_rows(inputs$tree_iq)) inputs$tree_iq$newick[[1]] else NULL)
  labels <- rownames(cop)
  tip_meta <- tip_meta %>% filter(label %in% labels)
  if (nrow(tip_meta) == 0) return(empty)

  node_depth <- node.depth.edgelength(phy)
  tip_nodes <- match(tip_meta$label, phy$tip.label)
  parent_nodes <- sapply(tip_nodes, function(tip) {
    row <- which(phy$edge[, 2] == tip)
    if (length(row) == 0) return(NA_integer_)
    phy$edge[row, 1]
  })
  terminal_bl <- node_depth[tip_nodes] - node_depth[parent_nodes]
  bg_labels <- tip_meta %>% filter(!is_query) %>% pull(label)

  nearest_bg_dist <- sapply(tip_meta$label, function(lbl) {
    if (length(bg_labels) == 0) return(NA_real_)
    min(cop[lbl, bg_labels], na.rm = TRUE)
  })
  nearest_bg_sample <- sapply(tip_meta$label, function(lbl) {
    if (length(bg_labels) == 0) return(NA_character_)
    bg_dists <- cop[lbl, bg_labels, drop = TRUE]
    names(bg_dists)[which.min(bg_dists)]
  })
  dist_to_mrca <- sapply(seq_along(tip_meta$label), function(i) {
    if (length(bg_labels) == 0) return(NA_real_)
    bg_dists <- cop[tip_meta$label[i], bg_labels, drop = TRUE]
    nb <- names(bg_dists)[which.min(bg_dists)]
    mrca <- getMRCA(phy, c(tip_meta$label[i], nb))
    if (is.null(mrca)) return(NA_real_)
    node_depth[tip_nodes[i]] - node_depth[mrca]
  })

  bl_pct <- percentile(terminal_bl, terminal_bl[!tip_meta$is_query])
  dist_pct <- percentile(nearest_bg_dist, nearest_bg_dist[!tip_meta$is_query])

  tibble(
    sample = tip_meta$sample,
    is_query = tip_meta$is_query,
    terminal_branch_length = terminal_bl,
    terminal_branch_length_percentile = bl_pct,
    distance_to_nearest_background = nearest_bg_dist,
    distance_to_nearest_background_percentile = dist_pct,
    distance_to_nearest_background_sample = nearest_bg_sample,
    distance_to_mrca_with_nearest_background = dist_to_mrca
  )
}

# -----------------------------------------------------------------------------
# 03 Clade placement
# -----------------------------------------------------------------------------
make_03_clade_placement <- function(inputs) {
  empty <- tibble(sample = character(), is_query = logical(), assigned_clade = character(),
                  assigned_lineage = character(), nearest_background_sample = character(),
                  nearest_background_distance = numeric(), min_distance_to_assigned_clade = numeric(),
                  mean_distance_to_assigned_clade = numeric(), n_background_in_assigned_clade = integer())
  if (!has_rows(inputs$tip_meta)) return(empty)
  tip_meta <- inputs$tip_meta %>%
    mutate(label = as.character(label), is_query = as.logical(is_query))
  cop <- cophenetic_matrix(if (has_rows(inputs$tree_iq)) inputs$tree_iq$newick[[1]] else NULL)
  if (is.null(cop)) {
    return(tip_meta %>% select(sample, is_query, assigned_clade = clade, assigned_lineage = lineage) %>%
             mutate(nearest_background_sample = NA_character_, nearest_background_distance = NA_real_,
                    min_distance_to_assigned_clade = NA_real_, mean_distance_to_assigned_clade = NA_real_,
                    n_background_in_assigned_clade = NA_integer_))
  }
  labels <- rownames(cop)
  tip_meta <- tip_meta %>% filter(label %in% labels)
  if (nrow(tip_meta) == 0) return(empty)
  bg_labels <- tip_meta %>% filter(!is_query) %>% pull(label)

  nearest_bg_dist <- sapply(tip_meta$label, function(lbl) {
    if (length(bg_labels) == 0) return(NA_real_)
    min(cop[lbl, bg_labels], na.rm = TRUE)
  })
  nearest_bg_sample <- sapply(tip_meta$label, function(lbl) {
    if (length(bg_labels) == 0) return(NA_character_)
    bg_dists <- cop[lbl, bg_labels, drop = TRUE]
    names(bg_dists)[which.min(bg_dists)]
  })

  clade_mat <- sapply(tip_meta$label, function(lbl) {
    my_clade <- tip_meta$clade[tip_meta$label == lbl]
    clade_bg <- tip_meta %>% filter(!is_query, clade == my_clade) %>% pull(label)
    c(min = if (length(clade_bg) > 0) min(cop[lbl, clade_bg], na.rm = TRUE) else NA_real_,
      mean = if (length(clade_bg) > 0) mean(cop[lbl, clade_bg], na.rm = TRUE) else NA_real_,
      n = length(clade_bg))
  })

  tibble(
    sample = tip_meta$sample,
    is_query = tip_meta$is_query,
    assigned_clade = tip_meta$clade,
    assigned_lineage = tip_meta$lineage,
    nearest_background_sample = nearest_bg_sample,
    nearest_background_distance = nearest_bg_dist,
    min_distance_to_assigned_clade = clade_mat["min", ],
    mean_distance_to_assigned_clade = clade_mat["mean", ],
    n_background_in_assigned_clade = as.integer(clade_mat["n", ])
  )
}

# -----------------------------------------------------------------------------
# 04 Distance distributions (query vs background)
# -----------------------------------------------------------------------------
make_04_distance_distributions <- function(inputs) {
  empty <- tibble(comparison = character(), n_query = integer(), n_background = integer(),
                  mean_distance = numeric(), median_distance = numeric(), min_distance = numeric(),
                  max_distance = numeric(), q25 = numeric(), q75 = numeric(), sd_distance = numeric())
  if (!has_rows(inputs$tip_meta)) return(empty)
  tip_meta <- inputs$tip_meta %>% mutate(label = as.character(label), is_query = as.logical(is_query))
  cop <- cophenetic_matrix(if (has_rows(inputs$tree_iq)) inputs$tree_iq$newick[[1]] else NULL)
  if (is.null(cop)) return(empty)
  labels <- rownames(cop)
  tip_meta <- tip_meta %>% filter(label %in% labels)
  q_labels <- tip_meta %>% filter(is_query) %>% pull(label)
  b_labels <- tip_meta %>% filter(!is_query) %>% pull(label)
  if (length(q_labels) == 0) return(empty)

  summarize_dist <- function(v) tibble(
    mean_distance = mean(v, na.rm = TRUE),
    median_distance = median(v, na.rm = TRUE),
    min_distance = min(v, na.rm = TRUE),
    max_distance = max(v, na.rm = TRUE),
    q25 = quantile(v, 0.25, na.rm = TRUE),
    q75 = quantile(v, 0.75, na.rm = TRUE),
    sd_distance = sd(v, na.rm = TRUE)
  )

  out <- tibble()
  if (length(b_labels) > 0) {
    out <- bind_rows(out, cbind(
      tibble(comparison = "query_vs_background", n_query = length(q_labels), n_background = length(b_labels)),
      summarize_dist(as.vector(cop[q_labels, b_labels]))))
  }
  if (length(b_labels) >= 2) {
    out <- bind_rows(out, cbind(
      tibble(comparison = "background_vs_background", n_query = 0L, n_background = length(b_labels)),
      summarize_dist(as.vector(cop[b_labels, b_labels]))))
  }
  if (length(q_labels) >= 2) {
    out <- bind_rows(out, cbind(
      tibble(comparison = "query_vs_query", n_query = length(q_labels), n_background = 0L),
      summarize_dist(as.vector(cop[q_labels, q_labels]))))
  }
  out
}

# -----------------------------------------------------------------------------
# 05 Temporal signal (root-to-tip regression)
# -----------------------------------------------------------------------------
make_05_temporal_signal <- function(inputs) {
  empty <- tibble(group = character(), n = integer(), n_with_date = integer(),
                  slope = numeric(), intercept = numeric(), r_squared = numeric(),
                  p_value = numeric(), rtt_r_squared = numeric(), rtt_slope = numeric(),
                  note = character())
  if (!has_rows(inputs$tip_meta)) return(empty)
  phy <- newick_to_phylo(if (has_rows(inputs$tree_iq)) inputs$tree_iq$newick[[1]] else NULL)
  if (is.null(phy)) return(empty)
  rtt <- node.depth.edgelength(phy)[1:length(phy$tip.label)]
  tip_meta <- inputs$tip_meta %>%
    mutate(label = as.character(label), is_query = as.logical(is_query), date = as.Date(collection_date))
  df <- tibble(label = phy$tip.label, rtt = rtt) %>%
    left_join(tip_meta, by = "label") %>%
    filter(!is.na(date))
  if (nrow(df) < 10) {
    return(tibble(group = "all", n = nrow(df), n_with_date = nrow(df), slope = NA_real_,
                  intercept = NA_real_, r_squared = NA_real_, p_value = NA_real_,
                  rtt_r_squared = NA_real_, rtt_slope = NA_real_, note = "insufficient_dated_samples"))
  }

  fit_group <- function(gdf) {
    if (nrow(gdf) < 5) {
      return(tibble(group = gdf$group[1], n = nrow(gdf), n_with_date = nrow(gdf),
                    slope = NA_real_, intercept = NA_real_, r_squared = NA_real_, p_value = NA_real_,
                    rtt_r_squared = NA_real_, rtt_slope = NA_real_, note = "n<5"))
    }
    gdf$days <- as.numeric(gdf$date - min(gdf$date, na.rm = TRUE))
    fit <- lm(rtt ~ days, data = gdf)
    s <- summary(fit)
    tibble(group = gdf$group[1], n = nrow(gdf), n_with_date = nrow(gdf),
           slope = coef(fit)["days"], intercept = coef(fit)["(Intercept)"],
           r_squared = s$r.squared, p_value = s$coefficients["days", "Pr(>|t|)"],
           rtt_r_squared = NA_real_, rtt_slope = NA_real_, note = "current_root_regression")
  }

  bind_rows(
    df %>% mutate(group = "all") %>% fit_group(),
    df %>% filter(!is_query) %>% mutate(group = "background") %>% fit_group(),
    df %>% filter(is_query) %>% mutate(group = "query") %>% fit_group()
  )
}

# -----------------------------------------------------------------------------
# 06 Differential mutations (query vs background, Fisher exact)
# -----------------------------------------------------------------------------
make_06_differential_mutations <- function(inputs) {
  empty <- tibble(mutation_id = integer(), mutation_label = character(), gene = character(),
                  ref_aa = character(), position = integer(), alt_aa = character(),
                  n_query = integer(), n_background = integer(), query_freq = numeric(),
                  background_freq = numeric(), private_query = logical(), n_clades = integer(),
                  odds_ratio = numeric(), lower_ci = numeric(), upper_ci = numeric(),
                  p_value = numeric(), fdr = numeric(), test_note = character())
  if (!has_rows(inputs$sample_mutation)) return(empty)
  sm <- inputs$sample_mutation %>%
    left_join(inputs$samples %>% select(sample_id, is_query), by = "sample_id") %>%
    mutate(is_query = as.logical(is_query))
  nq <- sum(sm$is_query == TRUE, na.rm = TRUE)
  nb <- sum(sm$is_query == FALSE, na.rm = TRUE)
  if (nq == 0) return(empty)

  mut_summary <- sm %>%
    group_by(mutation_id) %>%
    summarise(n_query = sum(is_query == TRUE, na.rm = TRUE),
              n_background = sum(is_query == FALSE, na.rm = TRUE), .groups = "drop") %>%
    mutate(a = n_query, b = nq - n_query, c = n_background, d = nb - n_background,
           odds_ratio = NA_real_, lower_ci = NA_real_, upper_ci = NA_real_,
           p_value = NA_real_, test_note = NA_character_)

  for (i in seq_len(nrow(mut_summary))) {
    mat <- matrix(c(mut_summary$a[i], mut_summary$b[i], mut_summary$c[i], mut_summary$d[i]),
                  nrow = 2, byrow = TRUE)
    ft <- safe_fisher(mat)
    mut_summary$odds_ratio[i] <- ft$estimate
    mut_summary$lower_ci[i] <- ft$conf.int[1]
    mut_summary$upper_ci[i] <- ft$conf.int[2]
    mut_summary$p_value[i] <- ft$p.value
    mut_summary$test_note[i] <- ft$method
  }

  mut_summary <- mut_summary %>%
    left_join(inputs$mutations, by = "mutation_id") %>%
    mutate(private_query = n_query > 0 & n_background == 0,
           query_freq = safe_div(n_query, nq),
           background_freq = safe_div(n_background, nb))

  clades_per_mut <- sm %>%
    left_join(inputs$tip_meta %>% select(sample_id, clade), by = "sample_id") %>%
    group_by(mutation_id) %>%
    summarise(n_clades = n_distinct(clade, na.rm = TRUE), .groups = "drop")

  mut_summary %>%
    left_join(clades_per_mut, by = "mutation_id") %>%
    mutate(fdr = p.adjust(p_value, method = "BH")) %>%
    select(mutation_id, mutation_label, gene, ref_aa, position, alt_aa,
           n_query, n_background, query_freq, background_freq, private_query,
           n_clades, odds_ratio, lower_ci, upper_ci, p_value, fdr, test_note)
}

# -----------------------------------------------------------------------------
# 07 Mutation portfolio (Jaccard / Sorensen / PERMANOVA)
# -----------------------------------------------------------------------------
make_07_mutation_portfolio <- function(inputs) {
  empty <- tibble(
    comparison = character(), n_pairs = integer(),
    mean_distance = numeric(), median_distance = numeric(), sd_distance = numeric(),
    permanova_term = character(), permanova_df = numeric(), permanova_sum_sq = numeric(),
    permanova_r2 = numeric(), permanova_f = numeric(), permanova_p = numeric(),
    permanova_test_note = character()
  )
  if (!has_rows(inputs$sample_mutation) || !has_rows(inputs$mutations)) return(empty)
  sm <- inputs$sample_mutation %>%
    left_join(inputs$mutations %>% select(mutation_id, mutation_label), by = "mutation_id") %>%
    left_join(inputs$samples %>% select(sample_id, is_query), by = "sample_id") %>%
    filter(!is.na(mutation_label))
  if (!has_rows(sm)) return(empty)
  mat <- sm %>%
    distinct(sample_id, mutation_label, is_query) %>%
    mutate(present = 1L) %>%
    pivot_wider(id_cols = sample_id, names_from = mutation_label, values_from = present, values_fill = 0L) %>%
    column_to_rownames("sample_id") %>%
    as.matrix()
  if (nrow(mat) < 3 || ncol(mat) < 2) return(empty)
  group <- as.character(inputs$samples$is_query[match(rownames(mat), inputs$samples$sample_id)])
  group[is.na(group)] <- "unknown"
  n_query <- sum(group == "TRUE", na.rm = TRUE)

  jac <- as.matrix(vegan::vegdist(mat, method = "jaccard"))
  sor <- as.matrix(vegan::vegdist(mat, method = "bray"))
  qidx <- which(group == "TRUE")
  bidx <- which(group == "FALSE")

  add_row <- function(m, label, idx1, idx2) {
    if (length(idx1) == 0 || length(idx2) == 0) return(NULL)
    vals <- as.vector(m[idx1, idx2])
    tibble(comparison = label, n_pairs = length(vals), mean_distance = mean(vals), median_distance = median(vals), sd_distance = sd(vals))
  }
  desc <- bind_rows(
    add_row(jac, "query_vs_background_jaccard", qidx, bidx),
    add_row(sor, "query_vs_background_sorensen", qidx, bidx),
    if (length(qidx) >= 2) add_row(jac, "query_vs_query_jaccard", qidx, qidx) else NULL,
    if (length(bidx) >= 2) add_row(jac, "background_vs_background_jaccard", bidx, bidx) else NULL
  )
  if (!has_rows(desc)) return(empty)

  adonis_out <- safe_adonis2(jac, group, n_query)
  if (nrow(adonis_out) > 0) {
    perm <- adonis_out %>%
      filter(term != "Residual") %>%
      transmute(permanova_term = term, permanova_df = df, permanova_sum_sq = sum_sq,
                permanova_r2 = R2, permanova_f = f_statistic, permanova_p = p_value,
                permanova_test_note = "permanova_adonis2") %>%
      dplyr::slice(1)
    desc <- bind_cols(desc, perm[rep(1, nrow(desc)), ])
  }
  desc
}

# -----------------------------------------------------------------------------
# 08 Per-position evolution
# -----------------------------------------------------------------------------
make_08_per_position_evolution <- function(inputs) {
  empty <- tibble(gene = character(), position = integer(), ref_aa = character(),
                  n_total = integer(), n_query = integer(), n_background = integer(),
                  query_freq = numeric(), background_freq = numeric(),
                  n_clades = integer(), n_countries = integer(),
                  first_seen_date = as.Date(character()), last_seen_date = as.Date(character()),
                  substitutions_per_year = numeric(), poisson_trend_p = numeric())
  if (!has_rows(inputs$sample_mutation) || !has_rows(inputs$mutations)) return(empty)
  nq <- sum(inputs$samples$is_query == TRUE, na.rm = TRUE)
  nb <- sum(inputs$samples$is_query == FALSE, na.rm = TRUE)
  sm <- inputs$sample_mutation %>%
    left_join(inputs$mutations %>% select(mutation_id, gene, position, ref_aa), by = "mutation_id") %>%
    left_join(inputs$samples %>% select(sample_id, is_query, collection_date, country, clade), by = "sample_id") %>%
    filter(!is.na(gene), !is.na(position))
  if (!has_rows(sm)) return(empty)
  sm %>%
    group_by(gene, position, ref_aa) %>%
    summarise(
      n_total = n_distinct(sample_id),
      n_query = n_distinct(sample_id[is_query == TRUE]),
      n_background = n_distinct(sample_id[is_query == FALSE]),
      n_clades = n_distinct(clade, na.rm = TRUE),
      n_countries = n_distinct(country, na.rm = TRUE),
      first_seen_date = min(as.Date(collection_date), na.rm = TRUE),
      last_seen_date = max(as.Date(collection_date), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      query_freq = safe_div(n_query, nq),
      background_freq = safe_div(n_background, nb),
      substitutions_per_year = safe_div(n_total, as.numeric(last_seen_date - first_seen_date) / 365.25),
      poisson_trend_p = NA_real_
    )
}

# -----------------------------------------------------------------------------
# 09 Mutation spatiotemporal context (for query mutations)
# -----------------------------------------------------------------------------
make_09_mutation_spatiotemporal_context <- function(inputs) {
  empty <- tibble(mutation_id = integer(), mutation_label = character(), gene = character(),
                  n_background_with_mutation = integer(), first_seen_date = as.Date(character()),
                  last_seen_date = as.Date(character()), countries = character(),
                  hosts = character(), outbreaks = character(), recent_2yr_count = integer())
  if (!has_rows(inputs$sample_mutation) || !has_rows(inputs$samples)) return(empty)
  query_mut_ids <- inputs$sample_mutation %>%
    left_join(inputs$samples %>% select(sample_id, is_query), by = "sample_id") %>%
    filter(is_query == TRUE) %>%
    distinct(mutation_id) %>%
    pull(mutation_id)
  if (length(query_mut_ids) == 0) return(empty)
  inputs$sample_mutation %>%
    filter(mutation_id %in% query_mut_ids) %>%
    left_join(inputs$mutations %>% select(mutation_id, mutation_label, gene), by = "mutation_id") %>%
    left_join(inputs$samples %>% select(sample_id, is_query, collection_date, country, host, outbreak), by = "sample_id") %>%
    group_by(mutation_id, mutation_label, gene) %>%
    summarise(
      n_background_with_mutation = sum(is_query == FALSE, na.rm = TRUE),
      first_seen_date = min(as.Date(collection_date), na.rm = TRUE),
      last_seen_date = max(as.Date(collection_date), na.rm = TRUE),
      countries = paste(sort(unique(na.omit(country))), collapse = ";"),
      hosts = paste(sort(unique(na.omit(host))), collapse = ";"),
      outbreaks = paste(sort(unique(na.omit(outbreak))), collapse = ";"),
      recent_2yr_count = sum(is_query == FALSE & as.Date(collection_date) >= (max(as.Date(collection_date), na.rm = TRUE) - 365.25 * 2), na.rm = TRUE),
      .groups = "drop"
    )
}

# -----------------------------------------------------------------------------
# 10 Phenotype effects
# -----------------------------------------------------------------------------
make_10_phenotype_effects <- function(inputs) {
  empty <- tibble(sample = character(), mutation_label = character(), gene = character(),
                  position = integer(), ref_aa = character(), alt_aa = character(),
                  phenotype = character(), effect = character(), evidence = character(),
                  source = character())
  if (!has_rows(inputs$sample_mutation) || !has_rows(inputs$mutation_phenotypes)) return(empty)
  inputs$sample_mutation %>%
    left_join(inputs$samples %>% select(sample_id, sample, is_query), by = "sample_id") %>%
    filter(is_query == TRUE) %>%
    left_join(inputs$mutations %>% select(mutation_id, mutation_label, gene, position, ref_aa, alt_aa), by = "mutation_id") %>%
    left_join(inputs$mutation_phenotypes, by = "mutation_id") %>%
    filter(!is.na(phenotype)) %>%
    select(sample, mutation_label, gene, position, ref_aa, alt_aa, phenotype, effect, evidence, source)
}

# -----------------------------------------------------------------------------
# 10b Phenotype enrichment (query mutations vs background mutations)
# -----------------------------------------------------------------------------
make_10b_phenotype_enrichment <- function(inputs) {
  empty <- tibble(n_query_annotated = integer(), n_query_not_annotated = integer(),
                  n_background_annotated = integer(), n_background_not_annotated = integer(),
                  odds_ratio = numeric(), lower_ci = numeric(), upper_ci = numeric(),
                  p_value = numeric(), test_note = character())
  if (!has_rows(inputs$sample_mutation)) return(empty)
  annotated_ids <- inputs$mutation_phenotypes %>%
    filter(!is.na(phenotype)) %>%
    distinct(mutation_id) %>%
    pull(mutation_id)
  sm <- inputs$sample_mutation %>%
    left_join(inputs$samples %>% select(sample_id, is_query), by = "sample_id") %>%
    mutate(annotated = mutation_id %in% annotated_ids, is_query = as.logical(is_query))
  sample_counts <- sm %>%
    group_by(sample_id, is_query) %>%
    summarise(n_annotated = sum(annotated), n_total = n(), .groups = "drop")
  q <- sample_counts %>% filter(is_query == TRUE)
  b <- sample_counts %>% filter(is_query == FALSE)
  nq_a <- sum(q$n_annotated); nq_na <- sum(q$n_total) - nq_a
  nb_a <- sum(b$n_annotated); nb_na <- sum(b$n_total) - nb_a
  ft <- safe_fisher(matrix(c(nq_a, nq_na, nb_a, nb_na), nrow = 2, byrow = TRUE))
  tibble(n_query_annotated = nq_a, n_query_not_annotated = nq_na,
         n_background_annotated = nb_a, n_background_not_annotated = nb_na,
         odds_ratio = ft$estimate, lower_ci = ft$conf.int[1], upper_ci = ft$conf.int[2],
         p_value = ft$p.value, test_note = ft$method)
}

# -----------------------------------------------------------------------------
# 11 Metadata context (country / host / outbreak)
# -----------------------------------------------------------------------------
make_11_metadata_context <- function(inputs) {
  empty <- tibble(variable = character(), value = character(), query_n = integer(),
                  background_n = integer(), query_freq = numeric(), background_freq = numeric())
  if (!has_rows(inputs$samples)) return(empty)
  samples <- inputs$samples %>% mutate(is_query = as.logical(is_query))
  nq <- sum(samples$is_query == TRUE, na.rm = TRUE)
  nb <- sum(samples$is_query == FALSE, na.rm = TRUE)
  if (nq == 0) return(empty)
  summarize_var <- function(var_name) {
    samples %>%
      filter(!is.na(.data[[var_name]]), .data[[var_name]] != "") %>%
      group_by(value = .data[[var_name]]) %>%
      summarise(query_n = sum(is_query == TRUE, na.rm = TRUE),
                background_n = sum(is_query == FALSE, na.rm = TRUE), .groups = "drop") %>%
      mutate(variable = var_name,
             query_freq = safe_div(query_n, nq),
             background_freq = safe_div(background_n, nb))
  }
  bind_rows(
    summarize_var("country"),
    summarize_var("host"),
    summarize_var("outbreak")
  ) %>% select(variable, value, query_n, background_n, query_freq, background_freq)
}

# -----------------------------------------------------------------------------
# 12 Literature evidence
# -----------------------------------------------------------------------------
make_12_literature_evidence <- function(inputs) {
  empty <- tibble(field = character(), value = character(), n = integer(), pmids = character())
  if (!has_rows(inputs$literature)) return(empty)
  inputs$literature %>%
    filter(!is.na(extraction_field)) %>%
    group_by(field = extraction_field, value = extraction_value) %>%
    summarise(n = n(), pmids = paste(unique(na.omit(pmid)), collapse = ";"), .groups = "drop") %>%
    arrange(desc(n))
}

# -----------------------------------------------------------------------------
# 13 Evidence synthesis (no automated risk/novelty calls)
# -----------------------------------------------------------------------------
make_13_evidence_synthesis <- function(inputs, outputs) {
  questions <- c("Q2 Genetic isolation", "Q4 Divergence distribution", "Q7 Distinguishing mutations",
                 "Q10 Phenotype effects", "Q11 Metadata context")
  tsvs <- c("01_genetic_distinctiveness.tsv", "02_divergence_distribution.tsv", "06_differential_mutations.tsv",
            "10b_phenotype_enrichment.tsv", "11_metadata_context.tsv")
  max_pct <- if (has_rows(outputs$distinctiveness)) max(outputs$distinctiveness$nearest_background_percentile, na.rm = TRUE) else NA_real_
  max_bl_pct <- if (has_rows(outputs$divergence)) max(outputs$divergence$terminal_branch_length_percentile, na.rm = TRUE) else NA_real_
  n_diff <- if (has_rows(outputs$diffmut)) sum(outputs$diffmut$fdr < 0.05 & outputs$diffmut$odds_ratio > 1 & !is.na(outputs$diffmut$fdr), na.rm = TRUE) else 0L
  n_private <- if (has_rows(outputs$diffmut)) sum(outputs$diffmut$private_query, na.rm = TRUE) else 0L
  phen <- if (has_rows(outputs$phen_enrich)) outputs$phen_enrich$p_value[[1]] else NA_real_
  has_meta <- has_rows(outputs$metadata)
  max_meta_diff <- if (has_meta) max(abs(outputs$metadata$query_freq - outputs$metadata$background_freq), na.rm = TRUE) else NA_real_

  support <- c(
    if (is.na(max_pct)) "low" else if (max_pct > 95) "high" else if (max_pct > 80) "moderate" else "low",
    if (is.na(max_bl_pct)) "low" else if (max_bl_pct > 95) "high" else if (max_bl_pct > 80) "moderate" else "low",
    if (n_diff > 0) "high" else if (n_private > 0) "moderate" else "low",
    if (is.na(phen)) "low" else if (phen < 0.05) "high" else if (has_rows(outputs$phen_effects)) "moderate" else "low",
    if (!has_meta) "low" else if (max_meta_diff > 0.5) "moderate" else "low"
  )
  stats <- c(
    sprintf("max_nearest_background_percentile=%.1f", max_pct),
    sprintf("max_terminal_branch_length_percentile=%.1f", max_bl_pct),
    sprintf("n_fdr05_or_private=%d", n_diff + n_private),
    sprintf("phenotype_enrichment_p=%.3g", phen),
    sprintf("max_freq_diff=%.2f", max_meta_diff)
  )
  caveats <- c(
    "nearest_background_percentile requires comparable background sampling",
    "terminal_branch_length_percentile is sensitive to tree rooting and rate variation",
    "Fisher exact tests are skipped when expected counts < 5; interpret counts directly",
    "phenotype annotations may be incomplete or conditional",
    "metadata frequency is descriptive; no formal association test"
  )
  gaps <- c(
    if (is.na(max_pct)) "no_tree_or_distances" else "",
    if (is.na(max_bl_pct)) "no_tree_or_branch_lengths" else "",
    if (!has_rows(outputs$diffmut) || !any(grepl("Fisher", outputs$diffmut$test_note), na.rm = TRUE)) "mutation_counts_too_small_for_fisher" else "",
    if (is.na(phen)) "no_phenotype_annotations" else "",
    if (!has_meta) "no_metadata" else ""
  )
  tibble(question = questions, support_level = support, key_statistic = stats,
         supporting_tsv = tsvs, caveat = caveats, evidence_gap = gaps)
}

# -----------------------------------------------------------------------------
# 14 nf-metro one-line summary
# -----------------------------------------------------------------------------
make_14_nfmetro_summary <- function(inputs, outputs) {
  nq <- if (has_rows(inputs$samples)) sum(inputs$samples$is_query == TRUE, na.rm = TRUE) else 0L
  nb <- if (has_rows(inputs$samples)) sum(inputs$samples$is_query == FALSE, na.rm = TRUE) else 0L
  max_pct <- if (has_rows(outputs$distinctiveness)) max(outputs$distinctiveness$nearest_background_percentile, na.rm = TRUE) else NA_real_
  n_private <- if (has_rows(outputs$diffmut)) sum(outputs$diffmut$private_query, na.rm = TRUE) else 0L
  n_diff <- if (has_rows(outputs$diffmut)) sum(outputs$diffmut$fdr < 0.05 & outputs$diffmut$odds_ratio > 1 & !is.na(outputs$diffmut$fdr), na.rm = TRUE) else 0L
  slope <- if (has_rows(outputs$temporal)) outputs$temporal$slope[outputs$temporal$group == "all"] else NA_real_
  n_phen <- if (has_rows(outputs$phen_effects)) nrow(outputs$phen_effects) else 0L
  gaps <- if (has_rows(outputs$synthesis)) sum(outputs$synthesis$evidence_gap != "", na.rm = TRUE) else 0L

  tibble(
    run_id = opts$`run-id`,
    species = opts$species,
    n_query = nq,
    n_background = nb,
    max_nearest_background_percentile = max_pct,
    n_private_mutations = n_private,
    n_differential_mutations_fdr05 = n_diff,
    temporal_signal_slope = slope,
    n_phenotype_annotated_query_mutations = n_phen,
    evidence_gap_count = gaps
  )
}

# -----------------------------------------------------------------------------
# Main execution
# -----------------------------------------------------------------------------
log_info("Preparing inputs for species ", opts$species)
inputs <- extract_inputs()

log_info("Building species assignment outputs")
outputs <- list()
outputs$species <- make_00_species_assignment(inputs)
outputs$seq_sim <- make_02_sequence_similarity(inputs)
outputs$msa <- make_03_msa_profile(inputs)
outputs$phylo <- make_04_phylogenetic_placement(inputs)
outputs$summary <- make_identification_summary(outputs$species, outputs$seq_sim, outputs$msa, outputs$phylo)
outputs$unresolved <- make_unresolved_samples(outputs$species)

log_info("Writing species identification outputs")
write_tsv(outputs$summary, "identification_summary.tsv", subdir = "species_identification")
write_tsv(outputs$unresolved, "unresolved_samples.tsv", subdir = "species_identification")

log_info("Writing MultiQC custom-content tables")
run_prefix <- opts$`run-id`
write_mqc_tsv(outputs$summary, paste0(run_prefix, "_identification_summary_mqc.tsv"),
              id = "identification_summary", section_name = "Species Identification Summary",
              description = "Per-sample confirmed species assignment, coverage, QC status, closest reference (and its outbreak), identity/genetic distance to that reference, phylogenetic clade, and MSA identity.",
              subdir = "species_identification/mqc")

log_info("Done. Outputs in ", outdir)
