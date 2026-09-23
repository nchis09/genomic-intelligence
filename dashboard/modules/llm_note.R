# LLM tree-interpretation helper
#
# Turns a phylogenetic tree + tip metadata into a short plain-language note
# explaining how the query genome(s) relate to the background isolates.
#
# Pipeline: .pg_evo_summary() builds a deterministic fact list (the
# "evolutionary summary" JSON layer) -> .ollama_note() sends it to a local
# Ollama model -> .pg_tree_note() orchestrates with a template fallback so
# the UI always has something to show.
#
# Config (env vars):
#   PG_OLLAMA_HOST   default http://localhost:11434
#   PG_OLLAMA_MODEL  default: first installed model (via /api/tags)

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
.pg_ollama_host <- function() {
  h <- Sys.getenv("PG_OLLAMA_HOST", "http://localhost:11434")
  sub("/+$", "", h)
}

.pg_ollama_model <- function() Sys.getenv("PG_OLLAMA_MODEL", "")

# rlang-style default (defined locally so this file works outside Shiny too)
if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a)) b else a

# ---------------------------------------------------------------------------
# Tree helpers (pure ape, no extra deps)
# ---------------------------------------------------------------------------

# Parent lookup vector: parent[child] = parent node (0 at the root)
.pg_parent_vec <- function(tr) {
  parent <- integer(max(tr$edge))
  parent[tr$edge[, 2]] <- tr$edge[, 1]
  parent
}

# Patristic distances from one query tip to ALL tips in a single DFS.
# Tracks the deepest qtip-ancestor on the path to each node (the MRCA), then
# dist(q, t) = dr[q] + dr[t] - 2*dr[mrca]. O(n) per query — no keep.tip, no
# cophenetic matrix — so it stays fast on 4000+ tip trees.
.pg_dists_from <- function(tr, qtip, dr, kids, parent) {
  n_tip <- length(tr$tip.label)
  in_anc <- logical(n_tip + tr$Nnode)
  nd <- qtip
  while (TRUE) {
    in_anc[nd] <- TRUE
    if (parent[nd] == 0) break
    nd <- parent[nd]
  }
  root <- which(parent == 0)[1]
  mrca <- integer(n_tip + tr$Nnode)
  stack_n <- root
  stack_c <- if (in_anc[root]) root else 0L
  while (length(stack_n)) {
    nd <- stack_n[length(stack_n)]
    cur <- stack_c[length(stack_c)]
    stack_n <- stack_n[-length(stack_n)]
    stack_c <- stack_c[-length(stack_c)]
    if (in_anc[nd]) cur <- nd
    if (nd <= n_tip) { mrca[nd] <- cur; next }
    ch <- kids[[as.character(nd)]]
    if (!is.null(ch)) {
      stack_n <- c(stack_n, ch)
      stack_c <- c(stack_c, rep(cur, length(ch)))
    }
  }
  dr[qtip] + dr[seq_len(n_tip)] - 2 * dr[mrca[seq_len(n_tip)]]
}

# k nearest background tips to a query tip, ranked by true patristic distance.
.pg_neighbors <- function(tr, qtip, bg_tips, k = 8,
                          kids = NULL, parent = NULL, dr = NULL) {
  if (length(bg_tips) == 0) return(stats::setNames(numeric(0), character(0)))
  if (is.null(parent)) parent <- .pg_parent_vec(tr)
  if (is.null(kids))   kids   <- split(tr$edge[, 2], tr$edge[, 1])
  if (is.null(dr))     dr     <- ape::node.depth.edgelength(tr)
  d <- .pg_dists_from(tr, qtip, dr, kids, parent)
  vals <- d[bg_tips]
  names(vals) <- tr$tip.label[bg_tips]
  head(sort(vals), k)
}

# ---------------------------------------------------------------------------
# Evolutionary summary — deterministic facts about query vs background
# ---------------------------------------------------------------------------

# Empty / "Unknown" / "unassigned" all mean "no information" in the metadata.
.pg_known <- function(x) !is.na(x) & nzchar(x) & x != "Unknown" & tolower(x) != "unassigned"

.pg_evo_summary <- function(tr, tip_meta) {
  if (is.null(tr) || is.null(tip_meta) || nrow(tip_meta) == 0) return(NULL)

  meta <- tip_meta[match(tr$tip.label, tip_meta$label), , drop = FALSE]
  is_q <- meta$is_query %in% c(TRUE, 1, "TRUE", "true", "1")
  is_q[is.na(is_q)] <- FALSE
  q_idx <- which(is_q)
  bg_idx <- which(!is_q)

  bg <- meta[bg_idx, , drop = FALSE]
  n_top <- function(x, n = 5) {
    x <- x[.pg_known(x)]
    if (!length(x)) return(character(0))
    head(names(sort(table(x), decreasing = TRUE)), n)
  }
  rng_date <- function(x) {
    x <- suppressWarnings(as.Date(x))
    x <- x[!is.na(x)]
    if (!length(x)) return(NULL)
    c(min = as.character(min(x)), max = as.character(max(x)))
  }

  bg_div <- suppressWarnings(as.numeric(bg$div))
  bg_div <- bg_div[!is.na(bg_div)]
  bg_div_med <- if (length(bg_div)) stats::median(bg_div) else NA_real_

  # Shared across all queries: one split, one parent vector, one root-to-node
  # distance vector — per-query cost is a single O(n) DFS.
  kids <- split(tr$edge[, 2], tr$edge[, 1])
  parent <- .pg_parent_vec(tr)
  dr <- ape::node.depth.edgelength(tr)

  queries <- lapply(q_idx, function(qi) {
    nb <- .pg_neighbors(tr, qi, bg_idx, k = 8,
                        kids = kids, parent = parent, dr = dr)
    nb_meta <- meta[match(names(nb), tr$tip.label), , drop = FALSE]
    qdiv <- suppressWarnings(as.numeric(meta$div[qi]))
    list(
      label = tr$tip.label[qi],
      clade = meta$clade[qi], outbreak = meta$outbreak[qi],
      country = meta$country[qi], date = as.character(meta$tip_date[qi]),
      div = qdiv,
      nuc_mutations = suppressWarnings(as.numeric(meta$nuc_mutation_count[qi])),
      aa_mutations = suppressWarnings(as.numeric(meta$aa_mutation_count[qi])),
      genome_coverage = suppressWarnings(as.numeric(meta$genome_coverage[qi])),
      div_vs_bg_median = if (!is.na(qdiv) && !is.na(bg_div_med)) qdiv - bg_div_med else NA,
      div_rel = if (is.na(qdiv) || is.na(bg_div_med)) NA_character_ else
        if (abs(qdiv - bg_div_med) < 0.2 * abs(bg_div_med)) "typical of" else
        if (qdiv > bg_div_med) "more diverged than" else "less diverged than",
      neighbors = lapply(seq_len(nrow(nb_meta)), function(i) {
        list(label = nb_meta$label[i], clade = nb_meta$clade[i],
             outbreak = nb_meta$outbreak[i], country = nb_meta$country[i],
             date = as.character(nb_meta$tip_date[i]),
             patristic_dist = unname(nb[i]))
      })
    )
  })

  # Group queries that share the identical nearest-neighbour set — they sit in
  # the same spot on the tree, so one sentence covers all of them instead of
  # repeating "X clusters closest with the same 8 genomes" per query.
  nb_key <- vapply(queries, function(q)
    paste(sort(vapply(q$neighbors, function(x) x$label, character(1))),
          collapse = "|"), character(1))
  groups <- lapply(split(seq_along(queries), nb_key), function(ii) {
    qs <- queries[ii]
    uniq <- function(f) {
      v <- unique(vapply(qs, function(q) {
        x <- q[[f]]; if (is.null(x) || length(x) == 0) NA_character_ else as.character(x)
      }, character(1)))
      if (length(v) == 1) v else NA_character_
    }
    dr <- unique(vapply(qs, function(q) q$div_rel, character(1)))
    dr <- dr[!is.na(dr)]
    list(
      labels    = vapply(qs, function(q) q$label, character(1)),
      clade     = uniq("clade"),
      outbreak  = uniq("outbreak"),
      country   = uniq("country"),
      div_rel   = if (length(dr) == 1) dr else NA_character_,
      neighbors = qs[[1]]$neighbors
    )
  })
  groups <- unname(groups)

  # Cluster evidence per query: (a) known — the query's own outbreak/clade, or
  # a dominant outbreak/clade among the neighbours themselves; (b) geographic —
  # neighbours dominated by one country even when no outbreak/clade is known.
  cluster_info <- lapply(queries, function(q) {
    nb <- q$neighbors
    if (!length(nb)) return(list(known = FALSE, country = NA_character_))
    same_ob <- if (.pg_known(q$outbreak))
      mean(vapply(nb, function(x) identical(x$outbreak, q$outbreak), logical(1))) else 0
    same_cl <- if (.pg_known(q$clade))
      mean(vapply(nb, function(x) identical(x$clade, q$clade), logical(1))) else 0
    nb_field <- function(f) {
      v <- vapply(nb, function(x) {
        x <- x[[f]]; if (is.null(x) || length(x) == 0) NA_character_ else as.character(x)
      }, character(1))
      v[.pg_known(v)]
    }
    nb_ob <- nb_field("outbreak"); nb_cl <- nb_field("clade"); nb_co <- nb_field("country")
    dom <- max(same_ob, same_cl,
               if (length(nb_ob)) max(table(nb_ob)) / length(nb) else 0,
               if (length(nb_cl)) max(table(nb_cl)) / length(nb) else 0)
    top_co <- if (length(nb_co) && max(table(nb_co)) / length(nb) >= 0.5)
      names(which.max(table(nb_co))) else NA_character_
    list(known = dom >= 0.5, country = top_co)
  })
  q_known <- vapply(cluster_info, function(x) isTRUE(x$known), logical(1))
  geo <- vapply(cluster_info, function(x) x$country, character(1))
  geo <- geo[.pg_known(geo)]

  list(
    n_tips = length(tr$tip.label),
    n_query = length(q_idx),
    n_background = length(bg_idx),
    queries = queries,
    query_in_known_cluster = length(q_known) > 0 && all(q_known),
    query_geo_cluster = if (length(geo)) names(which.max(table(geo))) else NA_character_,
    groups = groups,
    background = list(
      n_clades = length(unique(bg$clade[.pg_known(bg$clade)])),
      top_clades = n_top(bg$clade),
      n_outbreaks = length(unique(bg$outbreak[.pg_known(bg$outbreak)])),
      top_outbreaks = n_top(bg$outbreak),
      n_countries = length(unique(bg$country[.pg_known(bg$country)])),
      top_countries = n_top(bg$country),
      date_range = rng_date(bg$tip_date),
      div_median = bg_div_med,
      nuc_mut_median = suppressWarnings(stats::median(as.numeric(bg$nuc_mutation_count), na.rm = TRUE))
    )
  )
}

# ---------------------------------------------------------------------------
# Deterministic fallback note (no LLM needed)
# ---------------------------------------------------------------------------

# Describe a neighbour set: clade/lineage, outbreak, countries, date range —
# e.g. " — mostly EBOV-Clade-3 isolates from the EBOV-2018 outbreak in
# Democratic Republic of the Congo (2018\u20132021)". Empty string if nothing known.
.pg_nb_desc <- function(nb) {
  pick <- function(f) {
    v <- unique(vapply(nb, function(x) {
      x <- x[[f]]; if (is.null(x) || length(x) == 0) NA_character_ else as.character(x)
    }, character(1)))
    v[.pg_known(v)]
  }
  nb_cl <- pick("clade")
  nb_ob <- pick("outbreak")
  nb_c  <- pick("country")
  nb_d  <- suppressWarnings(as.Date(vapply(nb, function(x) {
    x <- x$date; if (is.null(x) || length(x) == 0) NA_character_ else as.character(x)
  }, character(1))))
  nb_d  <- nb_d[!is.na(nb_d)]

  desc <- character(0)
  if (length(nb_cl))      desc <- c(desc, paste0(paste(head(nb_cl, 2), collapse = "/"), " isolates"))
  if (length(nb_ob) == 1) desc <- c(desc, paste0("from the ", nb_ob, " outbreak"))
  if (length(nb_c))       desc <- c(desc, paste0("in ", paste(head(nb_c, 3), collapse = ", ")))
  if (length(nb_d))       desc <- c(desc, sprintf("(%s\u2013%s)", format(min(nb_d), "%Y"), format(max(nb_d), "%Y")))
  if (!length(desc)) return("")
  paste0(" — mostly ", paste(desc, collapse = " "))
}

.pg_note_template <- function(s) {
  if (is.null(s)) return("No tree data available to interpret.")
  bg <- s$background
  parts <- character(0)

  parts <- c(parts, sprintf(
    "This tree places %d genome%s against %d background isolate%s spanning %d countr%s and %d known outbreak%s (%s).",
    s$n_query, if (s$n_query == 1) "" else "s",
    s$n_background, if (s$n_background == 1) "" else "s",
    bg$n_countries, if (bg$n_countries == 1) "y" else "ies",
    bg$n_outbreaks, if (bg$n_outbreaks == 1) "" else "s",
    if (!is.null(bg$date_range)) paste0(bg$date_range["min"], " to ", bg$date_range["max"]) else "dates unknown"
  ))

  # One sentence per neighbour-group: queries sharing the same closest
  # background genomes are combined instead of repeated one by one.
  for (g in s$groups) {
    nb <- g$neighbors
    labels <- g$labels
    single <- length(labels) == 1
    subj <- if (single) labels else
      paste0(paste(head(labels, -1), collapse = ", "), " and ", tail(labels, 1))
    clade_phrase <- if (.pg_known(g$clade))
      paste0("clade ", g$clade)
      else if (single) "an uncharacterized clade" else "the same part of the tree"
    if (length(nb)) {
      parts <- c(parts, sprintf(
        "%s %s in %s and cluster%s closest with %s%d background genome%s%s.",
        subj,
        if (single) "falls" else "all fall",
        clade_phrase,
        if (single) "s" else "",
        if (single) "" else "the same ",
        length(nb), if (length(nb) == 1) "" else "s",
        .pg_nb_desc(nb)
      ))
    } else {
      parts <- c(parts, sprintf(
        "%s branch%s on %s own with no close background relatives in this tree.",
        subj, if (single) "es" else "", if (single) "its" else "their"))
    }
    if (!is.na(g$div_rel)) {
      div_word <- c("more diverged than" = "above",
                    "less diverged than" = "below",
                    "typical of" = "close to")[g$div_rel]
      parts <- c(parts, sprintf("%s divergence is %s the background median.",
        if (single) "Its" else "Their", div_word))
    }
  }

  parts <- c(parts, if (isTRUE(s$query_in_known_cluster))
    "Overall the query sits inside a known transmission cluster, consistent with the outbreak it was sampled from."
    else if (.pg_known(s$query_geo_cluster)) sprintf(
      "Overall the query genomes cluster with background isolates from %s — a geographic grouping, though no shared outbreak or clade annotation confirms a known transmission chain.",
      s$query_geo_cluster)
    else "Overall the query does not nest cleanly inside a single known cluster — treat placement with care.")
  paste(parts, collapse = " ")
}

# ---------------------------------------------------------------------------
# Ollama call
# ---------------------------------------------------------------------------
.ollama_note <- function(summary, timeout = 60) {
  fail <- function(e) list(text = NULL, source = "none", model = NULL, error = e)
  if (!requireNamespace("httr", quietly = TRUE) || !requireNamespace("jsonlite", quietly = TRUE))
    return(fail("httr/jsonlite not installed"))

  host <- .pg_ollama_host()
  model <- .pg_ollama_model()

  if (!nzchar(model)) {
    tags <- tryCatch(
      httr::GET(paste0(host, "/api/tags"), httr::timeout(5)),
      error = function(e) NULL
    )
    if (is.null(tags) || httr::status_code(tags) != 200)
      return(fail(paste0("Ollama unreachable at ", host)))
    tl <- tryCatch(
      jsonlite::fromJSON(httr::content(tags, "text", encoding = "UTF-8"), simplifyVector = FALSE),
      error = function(e) NULL
    )
    models <- vapply(tl$models %||% list(), function(m) m$name %||% m$model %||% "", character(1))
    models <- models[nzchar(models)]
    if (!length(models)) return(fail("No Ollama models installed (ollama pull <model>)"))
    model <- models[1]
  }

  # Send the grouped view only — per-query neighbour lists duplicate the same
  # records and push the LLM towards one repetitive sentence per genome.
  llm_summary <- summary
  llm_summary$queries <- NULL
  facts <- tryCatch(
    jsonlite::toJSON(llm_summary, auto_unbox = TRUE, na = "null", pretty = FALSE),
    error = function(e) NULL
  )
  if (is.null(facts)) return(fail("Could not serialize summary"))

  prompt <- paste0(
    "You are a genomic epidemiologist writing for a public-health reader who ",
    "cannot interpret phylogenetic trees. Below is a JSON summary of a tree ",
    "comparing query genome(s) to background isolates.\n\n",
    "Write 3-5 plain-language sentences explaining: where the query genomes sit, ",
    "which background genomes they are closest to — name the clade or lineage, ",
    "the countries those neighbours were isolated in, and their sampling dates ",
    "when given — whether the queries belong to a known transmission cluster, ",
    "outbreak, or geographic grouping (query_geo_cluster), and one caveat.\n",
    "Queries sharing the same neighbours are already combined in 'groups': ",
    "describe each group once (e.g. 'UG_28, UG_27 and UG_25 all cluster with...'), ",
    "never write one sentence per genome. Be concrete, use the numbers given, ",
    "no jargon, no markdown.\n\n",
    "SUMMARY:\n", facts
  )

  body <- list(model = model, prompt = prompt, stream = FALSE,
               options = list(temperature = 0.3, num_predict = 350))
  resp <- tryCatch(
    httr::POST(paste0(host, "/api/generate"),
               body = jsonlite::toJSON(body, auto_unbox = TRUE),
               httr::content_type_json(), httr::timeout(timeout)),
    error = function(e) e
  )
  if (inherits(resp, "error")) return(fail(paste0("Ollama request failed: ", resp$message)))
  if (httr::status_code(resp) != 200)
    return(fail(paste0("Ollama HTTP ", httr::status_code(resp))))

  txt <- tryCatch(
    jsonlite::fromJSON(httr::content(resp, "text", encoding = "UTF-8"))$response,
    error = function(e) NULL
  )
  if (is.null(txt) || !nzchar(trimws(txt))) return(fail("Empty Ollama response"))
  list(text = trimws(txt), source = "ollama", model = model, error = NULL)
}

# ---------------------------------------------------------------------------
# Orchestrator: LLM first, template fallback
# ---------------------------------------------------------------------------
.pg_tree_note <- function(summary) {
  if (is.null(summary)) return(list(text = "No tree data available to interpret.",
                                    source = "none", model = NULL, error = NULL))
  r <- .ollama_note(summary)
  if (!is.null(r$text) && nzchar(r$text)) return(r)
  list(text = .pg_note_template(summary), source = "template",
       model = r$model, error = r$error)
}
