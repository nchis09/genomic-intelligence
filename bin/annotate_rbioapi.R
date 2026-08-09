#!/usr/bin/env Rscript
# ============================================================
# annotate_rbioapi.R
#
# Query UniProt (mutagenesis + variation), STRING (virus-host
# interactions), and Reactome (pathways) using rbioapi for
# query-sample-specific proteins and mutations.
#
# Inputs:
#   --mutations   TSV with columns: sample, gene, position, ref_aa, alt_aa, mutation_label
#   --accessions  TXT file with one UniProt accession per line
#   --rbioapi_dir Path to cloned rbioapi repo (tools/rbioapi)
#   --species     Species ID (bdbv, ebov, sudv, tafv, restv)
#   --outdir      Output directory
#   --prefix      Output file prefix
#
# Outputs:
#   {prefix}_mutagenesis.tsv        - Curated mutagenesis at query positions
#   {prefix}_variation.tsv          - Natural variants with phenotype/disease
#   {prefix}_string_interactions.tsv - Virus-host protein interactions
#   {prefix}_reactome_pathways.tsv   - Host pathways for interactors
# ============================================================

# --- Parse CLI args ---
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 0) return(default)
  if (idx + 1 > length(args)) return(default)
  args[idx + 1]
}

mutations_file <- get_arg("--mutations")
accessions_file <- get_arg("--accessions")
rbioapi_dir     <- get_arg("--rbioapi_dir", "tools/rbioapi")
species         <- get_arg("--species", "bdbv")
outdir          <- get_arg("--outdir", ".")
prefix          <- get_arg("--prefix", "query")
uniprot_tsv_file <- get_arg("--uniprot_tsv")

if (is.null(mutations_file) || is.null(accessions_file)) {
  stop("Usage: Rscript annotate_rbioapi.R --mutations <tsv> --accessions <txt> --rbioapi_dir <dir> --species <id> --outdir <dir> --prefix <str>")
}

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# --- Null coalescing operator ---
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

# --- Source rbioapi R files ---
cat("Loading rbioapi from:", rbioapi_dir, "\n")
r_files <- list.files(file.path(rbioapi_dir, "R"), pattern = "\\.R$", full.names = TRUE)
for (f in r_files) {
  tryCatch(source(f), error = function(e) cat("  WARNING: failed to source", f, ":", e$message, "\n"))
}

# Manually call .onLoad to set rbioapi options
if (exists(".onLoad", mode = "function")) {
  .onLoad(NULL, "rbioapi")
}

# Set sensible options for batch mode
options(rba_verbose = TRUE, rba_skip_error = TRUE, rba_timeout = 120, rba_retry_max = 1)

cat("rbioapi loaded successfully.\n\n")

# --- Read inputs ---
mutations <- read.table(mutations_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
accessions <- readLines(accessions_file)
accessions <- accessions[nchar(accessions) > 0 & !grepl("^#", accessions)]
accessions <- unique(accessions)  # defensive: avoid duplicate API calls

cat("Mutations:", nrow(mutations), "rows\n")
cat("Accessions:", length(accessions), "entries\n\n")

# --- Accession -> gene / organism lookup (from UniProtKB download TSV, if provided) ---
# gene: labels mutagenesis/variation rows with the gene they belong to, so
#   they can be directly cross-referenced against query_mutations.tsv
#   (which is keyed by gene + position).
# organism: identifies which species an accession's curated data actually
#   belongs to. Strategy 6 in extract_query_proteins.py (reviewed-canonical
#   discovery) falls back to a related species' reviewed entry when the
#   exact species has no reviewed Swiss-Prot coverage (e.g. Bundibugyo/
#   Tai Forest/Reston ebolavirus) - that provenance is otherwise lost by
#   the time accessions reach this script, making cross-species rows
#   indistinguishable from genuinely species-specific ones.
accession_gene <- new.env()
accession_organism <- new.env()
if (!is.null(uniprot_tsv_file) && file.exists(uniprot_tsv_file)) {
  up_tsv <- tryCatch(
    read.table(uniprot_tsv_file, header = TRUE, sep = "\t", quote = "",
               fill = TRUE, comment.char = "", stringsAsFactors = FALSE,
               na.strings = c("", "NA", "NULL")),
    error = function(e) NULL
  )
  if (!is.null(up_tsv) && "Entry" %in% colnames(up_tsv)) {
    gene_col <- if ("Gene.Names" %in% colnames(up_tsv)) "Gene.Names" else
                if ("Gene Names" %in% colnames(up_tsv)) "Gene Names" else NA
    organism_col <- if ("Organism" %in% colnames(up_tsv)) "Organism" else NA
    for (i in seq_len(nrow(up_tsv))) {
      acc <- up_tsv$Entry[i]
      gene <- if (!is.na(gene_col)) up_tsv[[gene_col]][i] else ""
      assign(acc, gene %||% "", envir = accession_gene)
      organism <- if (!is.na(organism_col)) up_tsv[[organism_col]][i] else ""
      assign(acc, organism %||% "", envir = accession_organism)
    }
  }
}
lookup_gene <- function(acc) {
  if (exists(acc, envir = accession_gene)) get(acc, envir = accession_gene) %||% "" else ""
}
lookup_organism <- function(acc) {
  if (exists(acc, envir = accession_organism)) get(acc, envir = accession_organism) %||% "" else ""
}

# --- Species -> precise UniProt organism name, for cross-species detection ---
# Mirrors ORGANISM_MAP in bin/extract_query_proteins.py.
SPECIES_ORGANISM_MAP <- list(
  bdbv  = "Bundibugyo ebolavirus",
  ebov  = "Zaire ebolavirus",
  sudv  = "Sudan ebolavirus",
  tafv  = "Tai Forest ebolavirus",
  restv = "Reston ebolavirus"
)
expected_organism <- SPECIES_ORGANISM_MAP[[species]] %||% ""

# Returns TRUE/FALSE/NA (NA if organism is unknown for this accession).
is_cross_species_accession <- function(acc) {
  organism <- lookup_organism(acc)
  if (organism == "" || expected_organism == "") return(NA)
  !grepl(expected_organism, organism, fixed = TRUE)
}

# --- Species taxon ID mapping for STRING ---
species_taxon <- list(
  bdbv  = 186538,   # Bundibugyo virus
  ebov  = 186539,   # Zaire ebolavirus
  sudv  = 186540,   # Sudan virus
  tafv  = 186537,   # Tai Forest ebolavirus
  restv = 186536    # Reston virus
)
taxon_id <- species_taxon[[species]]
if (is.null(taxon_id)) {
  warning("Unknown species '", species, "', STRING queries will be skipped")
  taxon_id <- NA
}

# ============================================================
# 1. UniProt Mutagenesis — curated mutation effects at positions
# ============================================================
cat("=== UniProt Mutagenesis ===\n")

mutagenesis_rows <- list()

for (acc in accessions) {
  cat("  Querying mutagenesis for", acc, "...\n")
  result <- tryCatch({
    rba_uniprot_mutagenesis(accession = acc)
  }, error = function(e) {
    cat("    ERROR:", e$message, "\n")
    NULL
  })

  if (!is.null(result) && is.list(result) && length(result) > 0) {
    # result is a list with "features" sub-list
    features <- if ("features" %in% names(result)) result$features else result
    full_seq <- result$sequence %||% ""
    gene     <- lookup_gene(acc)
    source_species <- lookup_organism(acc)
    is_xspecies <- is_cross_species_accession(acc)
    for (feat in features) {
      if (is.list(feat)) {
        # EBI Proteins API returns flat "begin"/"end" string fields
        # (NOT nested under location$begin$position) — see
        # https://www.ebi.ac.uk/proteins/api/mutagenesis/{accession}
        pos_start <- suppressWarnings(as.integer(feat$begin %||% NA))
        pos_end   <- suppressWarnings(as.integer(feat$end %||% NA))
        mut_desc  <- feat$description %||% ""
        # Mutagenesis features don't carry the original (wild-type) AA;
        # derive it from the full protein sequence at position_start.
        orig_aa   <- if (!is.na(pos_start) && nchar(full_seq) >= pos_start) {
          substr(full_seq, pos_start, pos_start)
        } else ""
        alt_aas   <- if (!is.null(feat$alternativeSequence)) paste(unlist(feat$alternativeSequence), collapse = ",") else ""

        # Check if any query mutation falls within this position range
        matching_muts <- NULL
        if (!is.na(pos_start) && !is.na(pos_end)) {
          for (i in seq_len(nrow(mutations))) {
            qpos <- mutations$position[i]
            qgene <- mutations$gene[i]
            if (!is.na(qpos) && qpos >= pos_start && qpos <= pos_end &&
                (gene == "" || is.na(qgene) || grepl(qgene, gene, fixed = TRUE))) {
              matching_muts <- c(matching_muts, mutations$mutation_label[i])
            }
          }
        }

        mutagenesis_rows[[length(mutagenesis_rows) + 1]] <- data.frame(
          uniprot_accession = acc,
          gene = gene,
          source_species = source_species,
          is_cross_species = is_xspecies,
          position_start = pos_start,
          position_end = pos_end,
          original_aa = orig_aa,
          alternative_aa = alt_aas,
          description = mut_desc,
          query_mutation_match = if (!is.null(matching_muts)) paste(matching_muts, collapse = ";") else "",
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

mutagenesis_df <- if (length(mutagenesis_rows) > 0) do.call(rbind, mutagenesis_rows) else data.frame()
if (nrow(mutagenesis_df) > 0) {
  # Filter to rows that match query mutation positions + all others
  mut_out <- file.path(outdir, paste0(prefix, "_mutagenesis.tsv"))
  write.table(mutagenesis_df, mut_out, sep = "\t", quote = FALSE, row.names = FALSE)
  n_matched <- sum(mutagenesis_df$query_mutation_match != "")
  cat("  Mutagenesis entries:", nrow(mutagenesis_df), "(", n_matched, "match query positions)\n")
  cat("  Wrote:", mut_out, "\n\n")
} else {
  cat("  No mutagenesis data found.\n\n")
  mut_out <- file.path(outdir, paste0(prefix, "_mutagenesis.tsv"))
  write.table(data.frame(uniprot_accession=character(), gene=character(),
    source_species=character(), is_cross_species=logical(),
    position_start=integer(), position_end=integer(), original_aa=character(),
    alternative_aa=character(), description=character(), query_mutation_match=character()),
    mut_out, sep="\t", quote=FALSE, row.names=FALSE)
}

# ============================================================
# 2. UniProt Variation — natural variants with disease annotations
# ============================================================
cat("=== UniProt Variation ===\n")

variation_rows <- list()

for (acc in accessions) {
  cat("  Querying variation for", acc, "...\n")
  result <- tryCatch({
    rba_uniprot_variation(id = acc, id_type = "uniprot")
  }, error = function(e) {
    cat("    ERROR:", e$message, "\n")
    NULL
  })

  if (!is.null(result) && is.list(result) && length(result) > 0) {
    features <- if ("features" %in% names(result)) result$features else result
    gene <- lookup_gene(acc)
    source_species <- lookup_organism(acc)
    is_xspecies <- is_cross_species_accession(acc)
    for (feat in features) {
      if (is.list(feat)) {
        # EBI Proteins API returns flat "begin"/"end" string fields, plus
        # "wildType"/"mutatedType" (not "originalSequence"/"consequence")
        # and free-text "descriptions" (not a "disease" field) — see
        # https://www.ebi.ac.uk/proteins/api/variation/{accession}
        pos_start <- suppressWarnings(as.integer(feat$begin %||% NA))
        pos_end   <- suppressWarnings(as.integer(feat$end %||% NA))
        orig_aa   <- feat$wildType %||% ""
        alt_aas   <- feat$alternativeSequence %||% feat$mutatedType %||% ""
        consequence <- feat$consequenceType %||% ""
        notes      <- if (!is.null(feat$descriptions)) {
          paste(sapply(feat$descriptions, function(d) d$value %||% ""), collapse = "; ")
        } else ""
        source     <- feat$sourceType %||% ""
        evidence   <- if (!is.null(feat$evidences)) paste(sapply(feat$evidences, function(e) e$code), collapse = ";") else ""

        matching_muts <- NULL
        if (!is.na(pos_start) && !is.na(pos_end)) {
          for (i in seq_len(nrow(mutations))) {
            qpos <- mutations$position[i]
            qgene <- mutations$gene[i]
            if (!is.na(qpos) && qpos >= pos_start && qpos <= pos_end &&
                (gene == "" || is.na(qgene) || grepl(qgene, gene, fixed = TRUE))) {
              matching_muts <- c(matching_muts, mutations$mutation_label[i])
            }
          }
        }

        variation_rows[[length(variation_rows) + 1]] <- data.frame(
          uniprot_accession = acc,
          gene = gene,
          source_species = source_species,
          is_cross_species = is_xspecies,
          position_start = pos_start,
          position_end = pos_end,
          original_aa = orig_aa,
          alternative_aa = alt_aas,
          consequence = consequence,
          notes = notes,
          source_type = source,
          evidence = evidence,
          query_mutation_match = if (!is.null(matching_muts)) paste(matching_muts, collapse = ";") else "",
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

variation_df <- if (length(variation_rows) > 0) do.call(rbind, variation_rows) else data.frame()
if (nrow(variation_df) > 0) {
  var_out <- file.path(outdir, paste0(prefix, "_variation.tsv"))
  write.table(variation_df, var_out, sep = "\t", quote = FALSE, row.names = FALSE)
  n_matched <- sum(variation_df$query_mutation_match != "")
  cat("  Variation entries:", nrow(variation_df), "(", n_matched, "match query positions)\n")
  cat("  Wrote:", var_out, "\n\n")
} else {
  cat("  No variation data found.\n\n")
  var_out <- file.path(outdir, paste0(prefix, "_variation.tsv"))
  write.table(data.frame(uniprot_accession=character(), gene=character(),
    source_species=character(), is_cross_species=logical(),
    position_start=integer(), position_end=integer(), original_aa=character(),
    alternative_aa=character(), consequence=character(), notes=character(),
    source_type=character(), evidence=character(), query_mutation_match=character()),
    var_out, sep="\t", quote=FALSE, row.names=FALSE)
}

# ============================================================
# 3. STRING — virus-host protein interactions
# ============================================================
cat("=== STRING Interactions ===\n")

string_df <- data.frame()

if (!is.na(taxon_id)) {
  # Map UniProt accessions to STRING IDs for the viral species
  cat("  Mapping accessions to STRING IDs (taxon:", taxon_id, ")...\n")
  string_ids <- tryCatch({
    rba_string_map_ids(ids = accessions, species = taxon_id)
  }, error = function(e) {
    cat("    ERROR mapping IDs:", e$message, "\n")
    NULL
  })

  if (!is.null(string_ids) && is.data.frame(string_ids) && nrow(string_ids) > 0) {
    string_identifiers <- string_ids$stringId
    cat("  Mapped", nrow(string_ids), "of", length(accessions), "accessions to STRING\n")

    # Get interaction network — add_nodes to find host interactors
    cat("  Retrieving interaction network...\n")
    interactions <- tryCatch({
      rba_string_interactions_network(
        ids = string_identifiers,
        species = taxon_id,
        required_score = 400,
        add_nodes = 50
      )
    }, error = function(e) {
      cat("    ERROR getting interactions:", e$message, "\n")
      NULL
    })

    if (!is.null(interactions) && is.data.frame(interactions) && nrow(interactions) > 0) {
      string_df <- interactions
      str_out <- file.path(outdir, paste0(prefix, "_string_interactions.tsv"))
      write.table(string_df, str_out, sep = "\t", quote = FALSE, row.names = FALSE)
      cat("  Interaction entries:", nrow(string_df), "\n")
      cat("  Wrote:", str_out, "\n\n")
    } else {
      cat("  No interactions found.\n\n")
    }
  } else {
    cat("  No STRING IDs mapped for this species.\n\n")
  }
} else {
  cat("  Skipped (no taxon ID for species).\n\n")
}

if (nrow(string_df) == 0) {
  str_out <- file.path(outdir, paste0(prefix, "_string_interactions.tsv"))
  write.table(data.frame(stringId_A=character(), stringId_B=character(),
    preferredName_A=character(), preferredName_B=character(),
    score=numeric(), stringsAsFactors=FALSE),
    str_out, sep="\t", quote=FALSE, row.names=FALSE)
}

# ============================================================
# 4. Reactome — pathways for host interactors
# ============================================================
cat("=== Reactome Pathways ===\n")

reactome_rows <- list()

if (nrow(string_df) > 0 && "preferredName_B" %in% colnames(string_df)) {
  # Get unique interactor protein names (host side)
  host_proteins <- unique(string_df$preferredName_B)
  # Filter to non-viral (host) interactors — those not in our accessions
  viral_names <- if ("preferredName_A" %in% colnames(string_df)) unique(string_df$preferredName_A) else character()
  host_proteins <- setdiff(host_proteins, viral_names)

  cat("  Host interactor proteins:", length(host_proteins), "\n")

  for (protein in head(host_proteins, 20)) {  # Limit to 20 to avoid rate limits
    cat("  Querying Reactome for", protein, "...\n")
    pathways <- tryCatch({
      rba_reactome_mapping(
        id = protein,
        resource = "UniProt",
        map_to = "pathways",
        species = "Homo sapiens"
      )
    }, error = function(e) {
      cat("    ERROR:", e$message, "\n")
      NULL
    })

    if (!is.null(pathways) && length(pathways) > 0) {
      for (pw in pathways) {
        if (is.list(pw)) {
          reactome_rows[[length(reactome_rows) + 1]] <- data.frame(
            host_protein = protein,
            pathway_id = pw$dbId %||% pw$stId %||% "",
            pathway_name = pw$name %||% "",
            species = pw$species %||% "",
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
} else {
  cat("  No STRING interactions to query Reactome with.\n")
}

reactome_df <- if (length(reactome_rows) > 0) do.call(rbind, reactome_rows) else data.frame()
if (nrow(reactome_df) > 0) {
  # Deduplicate
  reactome_df <- unique(reactome_df)
  rea_out <- file.path(outdir, paste0(prefix, "_reactome_pathways.tsv"))
  write.table(reactome_df, rea_out, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("  Pathway entries:", nrow(reactome_df), "\n")
  cat("  Wrote:", rea_out, "\n\n")
} else {
  cat("  No pathway data found.\n\n")
  rea_out <- file.path(outdir, paste0(prefix, "_reactome_pathways.tsv"))
  write.table(data.frame(host_protein=character(), pathway_id=character(),
    pathway_name=character(), species=character()),
    rea_out, sep="\t", quote=FALSE, row.names=FALSE)
}

cat("Done.\n")
