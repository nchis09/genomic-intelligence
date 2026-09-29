#!/usr/bin/env Rscript
# Aggregate per-domain evidence_extracted.tsv files into:
#   - countermeasure_readiness.tsv
#   - knowledge_gaps.tsv
#
# Usage:
#   Rscript aggregate_countermeasures_gaps.R --species ebov --evidence-dir literature_evidence/ebov --outdir .

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a)) b else a

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  opt <- list(species = NULL, evidence_dir = NULL, outdir = ".")
  i <- 1
  while (i <= length(args)) {
    if (startsWith(args[i], "--") && i < length(args)) {
      key <- sub("^--", "", args[i])
      if (key == "species") opt$species <- args[i + 1]
      if (key == "evidence-dir") opt$evidence_dir <- args[i + 1]
      if (key == "outdir") opt$outdir <- args[i + 1]
      i <- i + 2
    } else {
      i <- i + 1
    }
  }
  opt
}

score_evidence <- function(level) {
  level <- tolower(level %||% "unknown")
  c("rct" = 1.0, "observational" = 0.75, "case_report" = 0.5,
    "in_silico" = 0.35, "review" = 0.4, "expert_opinion" = 0.25,
    "abstract_only" = 0.15, "unknown" = 0.0)[level] %||% 0.0
}

status_from_score <- function(score) {
  if (score >= 0.75) return("available")
  if (score >= 0.45) return("under_study")
  if (score >= 0.20) return("limited")
  "unknown"
}

read_evidence_tsvs <- function(evidence_dir) {
  paths <- list.files(evidence_dir, pattern = "evidence_extracted\\.tsv$",
                      recursive = TRUE, full.names = TRUE)
  if (!length(paths)) return(NULL)
  rows <- lapply(paths, function(p) {
    tryCatch(
      read_tsv(p, show_col_types = FALSE, col_types = cols(.default = "c")),
      error = function(e) { message("[warn] failed to read ", p, ": ", e$message); NULL }
    )
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) return(NULL)
  bind_rows(rows)
}

build_countermeasure_readiness <- function(ev) {
  if (is.null(ev) || !nrow(ev)) {
    return(tibble(
      species = character(), countermeasure = character(), status = character(),
      best_evidence = character(), best_pmid = character(), n_supporting_papers = integer(),
      readiness_score = numeric()
    ))
  }

  # Map topic to a countermeasure category
  ev <- ev |>
    mutate(countermeasure = case_when(
      tolower(topic) == "vaccine" ~ "vaccine",
      tolower(topic) == "therapeutic" & grepl("(?i)monoclonal|antibody|mab", product_name) ~ "monoclonal_antibody",
      tolower(topic) == "therapeutic" ~ "antiviral_therapeutic",
      tolower(topic) == "diagnostic" ~ "diagnostic",
      tolower(topic) == "surveillance" ~ "surveillance_package",
      TRUE ~ NA_character_
    )) |>
    filter(!is.na(countermeasure)) |>
    mutate(score = sapply(evidence_level, score_evidence))

  if (!nrow(ev)) {
    return(tibble(
      species = character(), countermeasure = character(), status = character(),
      best_evidence = character(), best_pmid = character(), n_supporting_papers = integer(),
      readiness_score = numeric()
    ))
  }

  ev |>
    group_by(species, countermeasure) |>
    arrange(desc(score)) |>
    summarise(
      best_evidence = first(finding),
      best_pmid = first(pmid),
      n_supporting_papers = length(unique(pmid)),
      readiness_score = max(score, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(status = sapply(readiness_score, status_from_score))
}

build_knowledge_gaps <- function(ev, species) {
  priority <- c("vaccine", "monoclonal_antibody", "antiviral_therapeutic", "diagnostic", "surveillance_package")
  present <- character()
  if (!is.null(ev) && nrow(ev)) {
    present <- tolower(ev$topic)
  }

  gaps <- lapply(priority, function(cm) {
    matched <- switch(cm,
      vaccine = "vaccine",
      monoclonal_antibody = "therapeutic",
      antiviral_therapeutic = "therapeutic",
      diagnostic = "diagnostic",
      surveillance_package = "surveillance"
    )
    topic_label <- cm
    topic_ebov <- ev |> filter(tolower(topic) == matched)
    n_total <- length(unique(topic_ebov$pmid))
    n_strong <- length(unique(topic_ebov$pmid[topic_ebov$evidence_level %in% c("RCT", "observational")]))

    if (n_total == 0) {
      return(tibble(
        species = species,
        topic = topic_label,
        gap = sprintf("No literature evidence found for %s efficacy/performance for %s.", topic_label, toupper(species)),
        priority = "high",
        n_papers_touching_topic = 0L,
        notes = "Search returned no directly relevant evidence."
      ))
    }
    if (n_strong == 0) {
      return(tibble(
        species = species,
        topic = topic_label,
        gap = sprintf("Only weak or abstract-only evidence for %s in %s; no solid efficacy/performance data.", topic_label, toupper(species)),
        priority = "medium",
        n_papers_touching_topic = as.integer(n_total),
        notes = "Relevant papers are reviews, expert opinion, or abstract-only."
      ))
    }
    NULL
  })
  gaps <- gaps[!vapply(gaps, is.null, logical(1))]
  if (!length(gaps)) {
    return(tibble(
      species = character(), topic = character(), gap = character(),
      priority = character(), n_papers_touching_topic = integer(), notes = character()
    ))
  }
  bind_rows(gaps)
}

main <- function() {
  opt <- parse_args()
  if (is.null(opt$species) || is.null(opt$evidence_dir)) {
    stop("Usage: aggregate_countermeasures_gaps.R --species <sp> --evidence-dir <dir> [--outdir <dir>]")
  }
  dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

  ev <- read_evidence_tsvs(opt$evidence_dir)
  if (is.null(ev) || !nrow(ev)) {
    message("[aggregate] No evidence TSVs found in ", opt$evidence_dir)
  }

  cm <- build_countermeasure_readiness(ev)
  kg <- build_knowledge_gaps(ev, opt$species)

  out_cm <- file.path(opt$outdir, "countermeasure_readiness.tsv")
  out_kg <- file.path(opt$outdir, "knowledge_gaps.tsv")

  write_tsv(cm, out_cm)
  write_tsv(kg, out_kg)

  message("[aggregate] wrote ", out_cm, " (", nrow(cm), " rows)")
  message("[aggregate] wrote ", out_kg, " (", nrow(kg), " rows)")
}

main()
