# Dashboard meta-analysis helpers for countermeasure readiness and knowledge gaps.
# Reads the raw evidence_extracted table (DuckDB or TSV) and derives countermeasure
# readiness and gap tables on the fly, in line with the pipeline knowledge warehouse.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a)) b else a

# ---------------------------------------------------------------------------
# Scoring helpers (mirror bin/aggregate_countermeasures_gaps.R)
# ---------------------------------------------------------------------------

cm_score_evidence <- function(level) {
  if (is.null(level) || length(level) == 0 || is.na(level) || !nzchar(level)) level <- "unknown"
  level <- tolower(level)
  c("rct" = 1.0, "observational" = 0.75, "case_report" = 0.5,
    "in_silico" = 0.35, "review" = 0.4, "expert_opinion" = 0.25,
    "abstract_only" = 0.15, "unknown" = 0.0)[level] %||% 0.0
}

cm_status_from_score <- function(score) {
  if (score >= 0.75) return("available")
  if (score >= 0.45) return("under_study")
  if (score >= 0.20) return("limited")
  "unknown"
}

# ---------------------------------------------------------------------------
# Load evidence for a species from DuckDB or TSVs
# ---------------------------------------------------------------------------

cm_load_evidence <- function(outdir, species) {
  db <- file.path(outdir, "knowledge_warehouse", "knowledge_warehouse.duckdb")

  if (file.exists(db) && requireNamespace("DBI", quietly = TRUE) &&
      requireNamespace("duckdb", quietly = TRUE)) {
    con <- DBI::dbConnect(duckdb::duckdb(), dbdir = db, read_only = TRUE)
    on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
    ev <- tryCatch(
      DBI::dbGetQuery(con, "SELECT * FROM evidence_extracted WHERE species = ?", params = list(species)),
      error = function(e) NULL
    )
    if (!is.null(ev) && nrow(ev)) {
      # Materialise per-domain fields from evidence_fields as columns, so
      # downstream meta-analysis works the same as the per-domain TSVs.
      fld <- tryCatch(
        DBI::dbGetQuery(con,
          "SELECT f.claim_id, f.field_name, f.field_value
           FROM evidence_fields f JOIN evidence_extracted e ON f.claim_id = e.claim_id
           WHERE e.species = ?", params = list(species)),
        error = function(e) NULL
      )
      wide <- ma_fields_wide(fld)
      if (is.null(wide)) wide <- ma_details_wide(ev)
      if (!is.null(wide) && "claim_id" %in% names(ev)) {
        # The base table has legacy product_name/sensitivity/specificity columns
        # that collide with the widened domain fields — drop them so the join
        # yields clean names instead of .x/.y suffixes.
        drop <- intersect(setdiff(names(ev), "claim_id"), names(wide))
        if (length(drop)) ev <- dplyr::select(ev, -dplyr::all_of(drop))
        ev <- dplyr::left_join(ev, wide, by = "claim_id")
      }
      ev <- .cm_fill_paper_identity(ev)
      return(ev)
    }
  }

  # Fallback: read published TSVs
  base <- file.path(outdir, "literature_retrieval", "literature_evidence", species)
  if (!dir.exists(base)) return(NULL)
  paths <- list.files(base, pattern = "evidence_extracted\\.tsv$",
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
  .cm_fill_paper_identity(bind_rows(rows))
}

# Within a (domain, pmid) group, when exactly one distinct non-empty value of
# an identity field exists, fill it onto sibling claims that left it empty —
# the same paper's safety/efficacy claims refer to the same product. Never
# merges across distinct values.
.cm_fill_paper_identity <- function(ev) {
  if (is.null(ev) || !nrow(ev)) return(ev)
  id_fields <- intersect(c("product_name", "test_name", "intervention_name",
                           "policy_name"), names(ev))
  for (f in id_fields) {
    v <- as.character(ev[[f]])
    grp <- interaction(ev$domain, ev$pmid, drop = TRUE)
    mode <- tapply(seq_along(v), grp, function(i) {
      vals <- unique(v[i][nzchar(v[i]) & !is.na(v[i])])
      if (length(vals) == 1) vals else NA_character_
    })
    fill <- mode[as.character(grp)]
    need <- !nzchar(v) | is.na(v)
    v[need & !is.na(fill)] <- fill[need & !is.na(fill)]
    ev[[f]] <- v
  }
  ev
}

# ---------------------------------------------------------------------------
# Derive countermeasure readiness
# ---------------------------------------------------------------------------

cm_build_readiness <- function(ev) {
  if (is.null(ev) || !nrow(ev)) {
    return(tibble(
      species = character(), countermeasure = character(), status = character(),
      best_evidence = character(), best_pmid = character(),
      n_supporting_papers = integer(), readiness_score = numeric()
    ))
  }

  if (!"product_name" %in% names(ev)) ev <- ev |> mutate(product_name = "")

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
    mutate(score = sapply(evidence_level, cm_score_evidence))

  if (!nrow(ev)) return(tibble(
    species = character(), countermeasure = character(), status = character(),
    best_evidence = character(), best_pmid = character(),
    n_supporting_papers = integer(), readiness_score = numeric()
  ))

  ev |>
    group_by(species, countermeasure) |>
    arrange(desc(score), .by_group = TRUE) |>
    summarise(
      best_evidence = first(finding),
      best_pmid = first(pmid),
      n_supporting_papers = length(unique(pmid)),
      readiness_score = max(score, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(status = sapply(readiness_score, cm_status_from_score))
}

# ---------------------------------------------------------------------------
# Derive knowledge gaps
# ---------------------------------------------------------------------------

cm_build_gaps <- function(ev, species) {
  priority <- c("vaccine", "monoclonal_antibody", "antiviral_therapeutic", "diagnostic", "surveillance_package")

  if (is.null(ev) || !nrow(ev)) {
    return(tibble(
      species = species,
      topic = priority,
      gap = sprintf("No literature evidence found for %s efficacy/performance for %s.", priority, toupper(species)),
      priority = "high",
      n_papers_touching_topic = 0L,
      notes = "Search returned no directly relevant evidence."
    ))
  }

  present <- tolower(ev$topic)

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
    n_strong <- length(unique(topic_ebov$pmid[tolower(topic_ebov$evidence_level) %in% c("rct", "observational")]))

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

# ---------------------------------------------------------------------------
# Main retrieval function used by dashboard modules
# ---------------------------------------------------------------------------

cm_gap_data <- function(outdir, species) {
  ev <- cm_load_evidence(outdir, species)
  list(
    evidence = ev,
    countermeasures = cm_build_readiness(ev),
    gaps = cm_build_gaps(ev, species),
    papers = cm_load_papers(outdir, species),
    summaries = cm_load_summaries(outdir, species)
  )
}

# ---------------------------------------------------------------------------
# Pre-generated domain narratives (domain_summaries table or TSV)
# ---------------------------------------------------------------------------

# summaries <- df: domain, summary, source, model, n_claims, n_papers
cm_load_summaries <- function(outdir, species) {
  db <- file.path(outdir, "knowledge_warehouse", "knowledge_warehouse.duckdb")
  if (file.exists(db) && requireNamespace("DBI", quietly = TRUE) &&
      requireNamespace("duckdb", quietly = TRUE)) {
    con <- DBI::dbConnect(duckdb::duckdb(), dbdir = db, read_only = TRUE)
    on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
    rows <- tryCatch(
      DBI::dbGetQuery(con,
        "SELECT domain, summary, source, model, n_claims, n_papers
         FROM domain_summaries WHERE species = ?", params = list(species)),
      error = function(e) NULL
    )
    if (!is.null(rows) && nrow(rows)) return(rows)
  }

  tsv <- file.path(outdir, "literature_retrieval", "literature_summaries",
                   species, "domain_summaries.tsv")
  if (!file.exists(tsv)) return(NULL)
  tryCatch(read_tsv(tsv, show_col_types = FALSE,
                    col_types = cols(.default = "c")),
           error = function(e) NULL)
}

# ---------------------------------------------------------------------------
# Paper universe per domain (literature_papers or search TSVs)
# ---------------------------------------------------------------------------

# papers  <- df: domain, screened, with_claims  (per species)
cm_load_papers <- function(outdir, species) {
  db <- file.path(outdir, "knowledge_warehouse", "knowledge_warehouse.duckdb")

  if (file.exists(db) && requireNamespace("DBI", quietly = TRUE) &&
      requireNamespace("duckdb", quietly = TRUE)) {
    con <- DBI::dbConnect(duckdb::duckdb(), dbdir = db, read_only = TRUE)
    on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
    papers <- tryCatch(
      DBI::dbGetQuery(con,
        "SELECT domain, pmid, status, title, year, doi, journal
         FROM literature_papers WHERE species = ?", params = list(species)),
      error = function(e) NULL
    )
    if (!is.null(papers) && nrow(papers)) {
      papers$with_claims <- tolower(papers$status) == "clean"
      return(papers)
    }
  }

  # Fallback: search TSVs (universe) + evidence pmids (with claims)
  base <- file.path(outdir, "literature_retrieval", "literature_search", species)
  evbase <- file.path(outdir, "literature_retrieval", "literature_evidence", species)
  if (!dir.exists(base)) return(NULL)
  search_paths <- list.files(base, pattern = "results\\.tsv$",
                             recursive = TRUE, full.names = TRUE)
  if (!length(search_paths)) return(NULL)
  rows <- lapply(search_paths, function(p) {
    d <- tryCatch(read_tsv(p, show_col_types = FALSE, col_types = cols(.default = "c")),
                  error = function(e) NULL)
    if (is.null(d) || !nrow(d)) return(NULL)
    dom <- basename(dirname(p))
    ev_tsv <- file.path(evbase, dom, "evidence_extracted.tsv")
    pmids_claims <- character(0)
    if (file.exists(ev_tsv)) {
      ev <- tryCatch(read_tsv(ev_tsv, show_col_types = FALSE, col_types = cols(.default = "c")),
                     error = function(e) NULL)
      if (!is.null(ev) && "pmid" %in% names(ev))
        pmids_claims <- unique(ev$pmid[!is.na(ev$pmid) & nzchar(ev$pmid)])
    }
    data.frame(domain = dom, pmid = as.character(d$id),
               status = ifelse(d$id %in% pmids_claims, "clean", "searched"),
               title = d$title %||% "", year = d$year %||% "",
               doi = d$doi %||% "", journal = NA_character_,
               stringsAsFactors = FALSE)
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) return(NULL)
  papers <- bind_rows(rows)
  papers$with_claims <- papers$status == "clean"
  papers
}

# Per-domain screened/with-claims summary for the coverage strip.
cm_paper_coverage <- function(papers) {
  if (is.null(papers) || !nrow(papers)) return(NULL)
  papers |>
    dplyr::group_by(domain) |>
    dplyr::summarise(
      screened = dplyr::n_distinct(pmid),
      with_claims = sum(with_claims, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(domain_label = gsub("_", " ", domain)) |>
    dplyr::arrange(dplyr::desc(screened))
}
