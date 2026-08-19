# Assessment heuristic (Pathogen Identification objective only)
#
# First-pass, rule-based, fully transparent heuristic. It is NOT a
# statistical/ML model and is meant to be replaced/tuned later.
#
# Assessment (biological/public-health severity: MONITOR / INVESTIGATE /
# ESCALATE) and Confidence (data-quality: High / Medium / Low / Insufficient
# evidence) are two SEPARATE axes. Poor sequencing quality or missing tables
# never by themselves push the assessment toward ESCALATE -- they only lower
# Confidence. Genuine novelty signals (large genetic distance, low sequence
# identity, distant phylogenetic placement) can pull the assessment toward
# INVESTIGATE, but only when Confidence is adequate (Medium/High). ESCALATE
# requires multiple corroborating novelty signals together, so with only
# Pathogen Identification data available this heuristic will rarely (if
# ever) reach ESCALATE in practice -- intentionally conservative until
# richer evidence (epi/outbreak/phenotype data) is wired in later.
#
# Relies on PI_SUMMARY_TABLE / pi_table_path() / pi_read_table_safe(),
# defined in modules/pathogen_identification.R (must be sourced first).
#
# Since PATHOGEN_IDENTIFICATION now emits a single identification_summary.tsv
# (CONFIRMED samples only; qc_status == "good" is implied by that filter, so
# it's no longer tracked as a separate confidence signal here), n_total is 1
# and evidence_availability is simply whether that table has rows.

compute_pi_assessment <- function(species, outdir) {
  sa <- pi_read_table_safe(pi_table_path(outdir, species, PI_SUMMARY_TABLE$filename))

  n_total <- 1
  n_available <- if (!is.null(sa) && nrow(sa) > 0) 1 else 0
  evidence_availability <- n_available / n_total

  first_num <- function(df, col) {
    if (is.null(df) || !(col %in% names(df)) || nrow(df) == 0) return(NA_real_)
    suppressWarnings(as.numeric(df[[col]][1]))
  }

  coverage         <- first_num(sa, "coverage")
  identity         <- first_num(sa, "pct_identity_closest_ref")
  genetic_distance <- first_num(sa, "genetic_distance_closest_ref")
  clade_dist       <- first_num(sa, "distance_to_assigned_clade")

  signals <- character()

  # ---- Confidence / data-quality axis ------------------------------------
  confidence_penalty <- 0

  if (n_available == 0) {
    return(list(
      species = species,
      assessment = "MONITOR",
      confidence = "Insufficient evidence",
      evidence_availability = 0,
      n_available = 0,
      n_total = n_total,
      signals = "No confirmed Pathogen Identification data available yet for this species."
    ))
  }

  if (!is.na(coverage)) {
    if (coverage < 0.90) {
      confidence_penalty <- confidence_penalty + 1
      signals <- c(signals, sprintf("Genome coverage is low (%.0f%%) \u2014 lowers confidence.", coverage * 100))
    } else {
      signals <- c(signals, sprintf("Genome coverage: %.0f%%.", coverage * 100))
    }
  }

  confidence <- if (confidence_penalty >= 2) {
    "Low"
  } else if (confidence_penalty >= 1) {
    "Medium"
  } else {
    "High"
  }
  confidence_adequate <- confidence %in% c("Medium", "High")

  # ---- Assessment axis (novelty signals only) ----------------------------
  novelty_signals <- 0

  if (confidence_adequate) {
    # genetic_distance_closest_ref is normally a 0-1 fraction (Biostrings
    # edit-distance / alignment length), but falls back to the legacy
    # (possibly raw-count) genetic_distance field when no aligned-FASTA row
    # is available -- guard against both scales.
    dist_is_fraction <- !is.na(genetic_distance) && genetic_distance <= 1
    dist_elevated <- if (dist_is_fraction) {
      !is.na(genetic_distance) && genetic_distance > 0.05
    } else {
      !is.na(genetic_distance) && genetic_distance > 50
    }
    if (dist_elevated) {
      novelty_signals <- novelty_signals + 1
      signals <- c(signals, sprintf(
        "Genetic distance to closest reference is elevated (%.4f) \u2014 possible novelty signal.",
        genetic_distance
      ))
    }
    if (!is.na(identity) && identity < 0.95) {
      novelty_signals <- novelty_signals + 1
      signals <- c(signals, sprintf(
        "Sequence identity to closest reference is %.1f%% \u2014 below the typical range.",
        identity * 100
      ))
    }
    if (!is.na(clade_dist) && clade_dist > 0.05) {
      novelty_signals <- novelty_signals + 1
      signals <- c(signals, sprintf(
        "Phylogenetic distance to assigned clade is elevated (%.3f) \u2014 possible novel lineage.",
        clade_dist
      ))
    }
  } else {
    signals <- c(signals, "Confidence is too low to reliably assess novelty signals.")
  }

  assessment <- if (novelty_signals >= 2) {
    "ESCALATE"
  } else if (novelty_signals >= 1) {
    "INVESTIGATE"
  } else {
    "MONITOR"
  }

  list(
    species = species,
    assessment = assessment,
    confidence = confidence,
    evidence_availability = evidence_availability,
    n_available = n_available,
    n_total = n_total,
    signals = utils::head(signals, 5)
  )
}

# Plain factual evidence summary for the Overview header card -- no
# assessment/confidence verdict, just what data is available. Used instead
# of compute_pi_assessment() where GIF should present/contextualize evidence
# rather than make a risk call.
pi_evidence_facts <- function(species, outdir) {
  sa <- pi_read_table_safe(pi_table_path(outdir, species, PI_SUMMARY_TABLE$filename))

  n_total <- 1
  n_available <- if (!is.null(sa) && nrow(sa) > 0) 1 else 0

  first_num <- function(df, col) {
    if (is.null(df) || !(col %in% names(df)) || nrow(df) == 0) return(NA_real_)
    suppressWarnings(as.numeric(df[[col]][1]))
  }

  list(
    species = species,
    n_available = n_available,
    n_total = n_total,
    qc_status = if (n_available > 0) "good" else NA_character_,
    coverage = first_num(sa, "coverage")
  )
}

# MONITOR -> institutional blue/neutral (steady-state, not "all clear").
# INVESTIGATE -> amber. ESCALATE -> red. `success` (green) is intentionally
# not used here -- it's reserved for genuinely reassuring states elsewhere.
assessment_status <- function(assessment) {
  switch(assessment,
    MONITOR = "primary",
    INVESTIGATE = "warning",
    ESCALATE = "danger",
    "secondary"
  )
}

# Confidence uses a separate neutral gray scale so it's never visually
# confused with the assessment's semantic color.
confidence_badge_class <- function(confidence) {
  switch(confidence,
    "High" = "badge-dark",
    "Medium" = "badge-secondary",
    "Low" = "badge-light",
    "badge-light"
  )
}
