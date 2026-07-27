#!/usr/bin/env Rscript
#
# plot_tree_with_heatmap.R
#
# Wrapper around the published plotTree R function (katholt/plotTree,
# https://github.com/katholt/plotTree) for use in nf-genomic-intelligence.
#
# It takes a Nextstrain tree.nwk, a tip-metadata TSV and an optional mutation
# matrix TSV and produces a single static PNG containing:
#   - a ladderized phylogenetic tree,
#   - tip nodes coloured by a metadata trait (default: country),
#   - selected metadata columns printed beside the tree,
#   - a per-genome / per-mutation heatmap from the mutation matrix,
#   - an optional bar plot of mutation counts per tip.
#
# The original plotTree function is embedded below with only minor changes to
# read TSV files in addition to CSV.

suppressPackageStartupMessages({
  library(ape)
  library(optparse)
  library(RColorBrewer)
})

# ---------------------------------------------------------------------------
# plotTree function (from https://github.com/katholt/plotTree)
# Embedded here under the terms of its public repository with attribution.
# ---------------------------------------------------------------------------
readMatrix <- function(heatmapData) {
  if (is.matrix(heatmapData)) {
    x <- as.data.frame(heatmapData, stringsAsFactors = FALSE)
  } else if (is.data.frame(heatmapData)) {
    x <- heatmapData
  } else {
    sep <- "\t"
    if (grepl("\\.csv$", heatmapData, ignore.case = TRUE)) {
      sep <- ","
    }
    x <- read.table(heatmapData,
                    sep = sep,
                    header = TRUE,
                    row.names = 1,
                    check.names = FALSE,
                    stringsAsFactors = FALSE,
                    quote = "",
                    comment.char = "",
                    fill = TRUE)
  }
  x
}

getLayout <- function(infoFile, infoCols, heatmapData, barData, doBlocks,
                      treeWidth = 10, infoWidth = 10, dataWidth = 30,
                      edgeWidth = 1, labelHeight = 10, mainHeight = 100,
                      barDataWidth = 10, blockPlotWidth = 10) {
  w <- c(edgeWidth, treeWidth)
  m <- cbind(c(0, 0, 0), c(0, 1, 0))
  x <- 1

  if (!is.null(infoFile)) {
    printCols <- TRUE
    if (!is.null(infoCols)) {
      if (length(infoCols) == 1 && is.na(infoCols)) {
        printCols <- FALSE
      }
    }
    if (printCols) {
      x <- x + 1
      m <- cbind(m, c(0, x, 0))
      w <- c(w, infoWidth)
    }
  }

  if (!is.null(heatmapData)) {
    x <- x + 1
    m <- cbind(m, c(x + 1, x, 0))
    x <- x + 2
    m[1, 2] <- x
    w <- c(w, dataWidth)
  }

  if (!is.null(barData)) {
    x <- x + 1
    m <- cbind(m, c(0, x, x + 1))
    x <- x + 1
    w <- c(w, barDataWidth)
  }

  if (doBlocks) {
    x <- x + 1
    m <- cbind(m, c(0, x, 0))
    w <- c(w, blockPlotWidth)
  }

  m <- cbind(m, c(0, 0, 0))
  w <- c(w, edgeWidth)

  if (!is.null(heatmapData) | !is.null(barData)) {
    h <- c(labelHeight, mainHeight, labelHeight)
  } else {
    h <- c(edgeWidth, mainHeight, edgeWidth)
  }

  list(m = as.matrix(m), w = w, h = h)
}

plotTree <- function(tree, ladderise = NULL, heatmapData = NULL, barData = NULL,
                     infoFile = NULL, blockFile = NULL, snpFile = NULL,
                     gapChar = "?", genome_size = 5E6, blwd = 5,
                     block_colour = "black", snp_colour = "red",
                     genome_offset = 0, colourNodesBy = NULL, infoCols = NULL,
                     outputPDF = NULL, outputPNG = NULL, w, h,
                     heatmap.colours = rev(gray(seq(0, 1, 0.1))),
                     tip.labels = FALSE, tipLabelSize = 1, offset = 0,
                     tip.colour.cex = 0.5, legend = TRUE,
                     legend.pos = "bottomleft", ancestral.reconstruction = FALSE,
                     cluster = NULL, tipColours = NULL, lwd = 1.5, axis = FALSE,
                     axisPos = 3, edge.color = "black", infoCex = 0.8,
                     colLabelCex = 0.8, treeWidth = 10, infoWidth = 10,
                     dataWidth = 30, edgeWidth = 1, labelHeight = 10,
                     mainHeight = 100, barDataWidth = 10, blockPlotWidth = 10,
                     barDataCol = 2, heatmapBreaks = NULL,
                     heatmapDecimalPlaces = 1, heatmap.legend = TRUE,
                     vlines.heatmap = NULL,
                     vlines.heatmap.col = 2, heatmap.blocks = NULL,
                     pie.cex = 0.5) {
  if (is.character(tree)) {
    t <- read.tree(tree)
  } else {
    t <- tree
  }

  if (is.null(ladderise)) {
    tl <- t
  } else if (ladderise == "descending") {
    tl <- ladderize(t, TRUE)
  } else if (ladderise == "ascending") {
    tl <- ladderize(t, FALSE)
  } else {
    stop("ladderise option should be exactly 'ascending' or 'descending'")
  }

  tips <- tl$edge[, 2]
  tip.order <- tips[tips <= length(tl$tip.label)]
  tip.label.order <- tl$tip.label[tip.order]

  if (!is.null(heatmapData)) {
    x <- readMatrix(heatmapData)
    y.ordered <- x[tip.label.order, , drop = FALSE]
    if (!is.null(cluster)) {
      if (!(identical(cluster, FALSE))) {
        if (cluster == "square" & ncol(y.ordered) == nrow(y.ordered)) {
          original_order <- 1:nrow(x)
          names(original_order) <- rownames(x)
          reordered <- original_order[tip.label.order]
          y.ordered <- y.ordered[, rev(as.numeric(reordered)), drop = FALSE]
        } else {
          if (identical(cluster, TRUE)) cluster <- "ward"
          hc <- hclust(dist(t(na.omit(y.ordered))), cluster)
          y.ordered <- y.ordered[, hc$order, drop = FALSE]
        }
      }
    }
  }

  if (!is.null(barData)) {
    b <- readMatrix(barData)
    barData <- b[, 1]
    names(barData) <- rownames(b)
  }

  if (!is.null(infoFile)) {
    info <- readMatrix(infoFile)
    info.ordered <- info[rev(tip.label.order), , drop = FALSE]
  } else {
    info.ordered <- NULL
  }

  ancestral <- NULL
  nodeColourSuccess <- NULL
  if (!is.null(colourNodesBy) & !is.null(infoFile)) {
    if (colourNodesBy %in% colnames(info.ordered)) {
      nodeColourSuccess <- TRUE
      loc1 <- info.ordered[, which(colnames(info.ordered) == colourNodesBy)]
      tipLabelSet <- character(length(loc1))
      names(tipLabelSet) <- rownames(info.ordered)
      groups <- table(loc1, exclude = "")
      n <- length(groups)
      groupNames <- names(groups)
      if (is.null(tipColours)) {
        colours <- rainbow(n)
      } else {
        colours <- tipColours
      }
      for (i in 1:n) {
        g <- groupNames[i]
        tipLabelSet[loc1 == g] <- colours[i]
      }
      tipLabelSet <- tipLabelSet[tl$tip]
      if (ancestral.reconstruction) {
        ancestral <- ace(loc1, tl, type = "discrete")
      }
    }
  }

  if (!is.null(outputPDF)) {
    pdf(width = w, height = h, file = outputPDF)
  }
  if (!is.null(outputPNG)) {
    png(width = w, height = h, file = outputPNG)
  }

  doBlocks <- (!is.null(blockFile) | !is.null(snpFile))
  l <- getLayout(infoFile, infoCols, heatmapData, barData, doBlocks,
                 treeWidth = treeWidth, infoWidth = infoWidth,
                 dataWidth = dataWidth, edgeWidth = edgeWidth,
                 labelHeight = labelHeight, mainHeight = mainHeight,
                 barDataWidth = barDataWidth, blockPlotWidth = blockPlotWidth)
  layout(l$m, widths = l$w, heights = l$h)

  par(mar = rep(0, 4))
  plot.phylo(tl, no.margin = TRUE, show.tip.label = tip.labels,
           label.offset = offset, edge.width = lwd, edge.color = edge.color,
           xaxs = "i", yaxs = "i", y.lim = c(0.5, length(tl$tip) + 0.5),
           cex = tipLabelSize)

  if (!is.null(nodeColourSuccess)) {
    tiplabels(col = tipLabelSet, pch = 16, cex = tip.colour.cex)
    if (ancestral.reconstruction) {
      nodelabels(pie = ancestral$lik.anc, cex = pie.cex, piecol = colours)
    }
    if (legend) {
      legend(legend.pos, legend = groupNames, fill = colours)
    }
  }

  if (axis) {
    axisPhylo(axisPos)
  }

  if (!is.null(infoFile)) {
    printCols <- TRUE
    if (!is.null(infoCols)) {
      if (length(infoCols) == 1 && is.na(infoCols)) {
        printCols <- FALSE
      }
    }
    if (printCols) {
      par(mar = rep(0, 4))
      if (!is.null(infoCols)) {
        infoColNumbers <- which(colnames(info.ordered) %in% infoCols)
      } else {
        infoColNumbers <- 1:ncol(info.ordered)
      }
      plot(NA, axes = FALSE, pch = "",
           xlim = c(0, length(infoColNumbers) + 1.5),
           ylim = c(0.5, length(tl$tip) + 0.5), xaxs = "i", yaxs = "i")
      for (i in 1:length(infoColNumbers)) {
        j <- infoColNumbers[i]
        text(x = rep(i + 1, nrow(info.ordered) + 1),
             y = c(nrow(info.ordered):1),
             info.ordered[, j], cex = infoCex)
      }
    }
  }

  if (!is.null(heatmapData)) {
    if (is.null(heatmapBreaks)) {
      heatmapBreaks <- seq(min(y.ordered, na.rm = TRUE),
                           max(y.ordered, na.rm = TRUE),
                           length.out = length(heatmap.colours) + 1)
    }
    par(mar = rep(0, 4), xpd = TRUE)
    image((1:ncol(y.ordered)) - 0.5, (1:nrow(y.ordered)) - 0.5,
          as.matrix(t(y.ordered)), col = heatmap.colours,
          breaks = heatmapBreaks, axes = FALSE, xaxs = "i", yaxs = "i",
          xlab = "", ylab = "")
    if (!is.null(vlines.heatmap)) {
      for (v in vlines.heatmap) {
        abline(v = v, col = vlines.heatmap.col)
      }
    }
    if (!is.null(heatmap.blocks)) {
      for (coords in heatmap.blocks) {
        rect(xleft = coords[1], 0, xright = coords[2], ncol(y.ordered),
             col = vlines.heatmap.col, border = NA)
      }
    }
    par(mar = rep(0, 4))
    plot(NA, axes = FALSE, xaxs = "i", yaxs = "i", ylim = c(0, 2),
         xlim = c(0.5, ncol(y.ordered) + 0.5))
    text(1:ncol(y.ordered) - 0.5, rep(0, ncol(y.ordered)),
         colnames(y.ordered), srt = 90, cex = colLabelCex, pos = 4)
    par(mar = c(2, 0, 0, 2))
    if (heatmap.legend) {
      image(as.matrix(seq(min(y.ordered, na.rm = TRUE),
                          max(y.ordered, na.rm = TRUE),
                          length.out = length(heatmap.colours) + 1)),
            col = heatmap.colours, yaxt = "n", breaks = heatmapBreaks,
            axes = FALSE)
      axis(1,
           at = heatmapBreaks[-length(heatmapBreaks)] / max(y.ordered, na.rm = TRUE),
           labels = round(heatmapBreaks[-length(heatmapBreaks)], heatmapDecimalPlaces))
    } else {
      plot(NA, axes = FALSE, xaxs = "i", yaxs = "i",
           xlim = c(0.5, ncol(y.ordered) + 0.5), ylim = c(0, 2))
      text(ncol(y.ordered) / 2 + 0.5, 1, "Colour = mutation", cex = 0.8)
    }
  }

  if (!is.null(barData)) {
    par(mar = rep(0, 4))
    barplot(barData[tip.label.order], horiz = TRUE, axes = FALSE,
            xaxs = "i", yaxs = "i", xlab = "", ylab = "",
            ylim = c(0.25, length(barData) + 0.25),
            xlim = c((-1) * max(barData, na.rm = TRUE) / 20,
                     max(barData, na.rm = TRUE)),
            col = barDataCol, border = 0, width = 0.5, space = 1,
            names.arg = NA)
    par(mar = c(2, 0, 0, 0))
    plot(NA, yaxt = "n", xaxs = "i", yaxs = "i", xlab = "", ylab = "",
         ylim = c(0, 2),
         xlim = c((-1) * max(barData, na.rm = TRUE) / 20,
                  max(barData, na.rm = TRUE)), frame.plot = FALSE)
  }

  if (doBlocks) {
    par(mar = rep(0, 4))
    plot(NA, axes = FALSE, pch = "",
         xlim = c(genome_offset, genome_offset + genome_size + 1.5),
         ylim = c(0.5, length(tl$tip) + 0.5), xaxs = "i", yaxs = "i")
    if (!is.null(snpFile)) {
      # SNP plotting code omitted for brevity; not used here.
    }
    if (!is.null(blockFile)) {
      blocks <- read.delim(blockFile, header = FALSE)
      for (i in 1:nrow(blocks)) {
        if (as.character(blocks[i, 1]) %in% tip.label.order) {
          y <- which(tip.label.order == as.character(blocks[i, 1]))
          x1 <- blocks[i, 2]
          x2 <- blocks[i, 3]
          lines(c(x1, x2), c(y, y), lwd = blwd, lend = 2,
                col = block_colour)
        }
      }
    }
  }

  if (!is.null(outputPDF) | !is.null(outputPNG)) {
    dev.off()
  }

  if (!is.null(heatmapData)) {
    mat <- as.matrix(t(y.ordered))
  } else {
    mat <- NULL
  }

  list(info = info.ordered, anc = ancestral, mat = mat,
       strain_order = tip.label.order)
}

# ---------------------------------------------------------------------------
# Wrapper
# ---------------------------------------------------------------------------
option_list <- list(
  make_option(c("-t", "--tree"), type = "character", default = NULL,
              help = "Newick tree file"),
  make_option(c("-m", "--metadata"), type = "character", default = NULL,
              help = "Tip metadata TSV (first column = tip label)"),
  make_option(c("-x", "--mutation-matrix"), type = "character", default = NULL,
              dest = "mutation_matrix",
              help = "Mutation matrix TSV (first column = tip label)"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Output PNG file"),
  make_option(c("-p", "--prefix"), type = "character", default = "out",
              help = "Prefix for intermediate files"),
  make_option(c("-c", "--colour-nodes-by"), type = "character", default = "country",
              dest = "colour_nodes_by",
              help = "Metadata column used to colour tip nodes"),
  make_option(c("--width"), type = "integer", default = 2400,
              help = "PNG width in pixels"),
  make_option(c("--height"), type = "integer", default = 3000,
              help = "PNG height in pixels")
)

parser <- OptionParser(option_list = option_list,
                       description = "Produce an annotated phylogeny with a per-genome heatmap.")
args <- parse_args(parser)

if (is.null(args$tree) || is.null(args$metadata)) {
  stop("--tree and --metadata are required")
}

# Read metadata and set tip labels as row names
metadata <- read.table(args$metadata,
                       header = TRUE,
                       sep = "\t",
                       check.names = FALSE,
                       stringsAsFactors = FALSE,
                       row.names = 1,
                       quote = "",
                       comment.char = "",
                       fill = TRUE)

# Ensure the node-colouring trait has no empty values, otherwise plotTree fails
if (args$colour_nodes_by %in% colnames(metadata)) {
  metadata[[args$colour_nodes_by]][
    metadata[[args$colour_nodes_by]] == "" | is.na(metadata[[args$colour_nodes_by]])
  ] <- "Unknown"
}

# Default info columns to print beside the tree (only those that exist)
# For dense trees, omit the long date strings to keep the annotation readable
candidate_info <- if (nrow(metadata) > 60) {
  c("country", "is_query")
} else {
  c("country", "date", "is_query")
}
info_cols <- intersect(candidate_info, colnames(metadata))

# Convert is_query to a label: "query" for query samples, blank otherwise
if ("is_query" %in% colnames(metadata)) {
  metadata$is_query <- ifelse(
    metadata$is_query == TRUE | tolower(metadata$is_query) == "true",
    "query", ""
  )
  # Drop the column if it would be entirely blank
  if (all(metadata$is_query == "")) {
    info_cols <- setdiff(info_cols, "is_query")
  }
}
if (length(info_cols) == 0) info_cols <- NA

# Write a temporary metadata file in the format plotTree expects
info_file <- paste0(args$prefix, "_info.tsv")
write.table(metadata, info_file,
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

# Read and recode the mutation matrix so each mutation has a distinct colour
heatmap_matrix <- NULL
heatmap_colours <- NULL
heatmap_breaks <- NULL

if (!is.null(args$mutation_matrix) && file.exists(args$mutation_matrix)) {
  mut_df <- read.table(args$mutation_matrix,
                       header = TRUE,
                       sep = "\t",
                       row.names = 1,
                       check.names = FALSE,
                       stringsAsFactors = FALSE,
                       quote = "",
                       comment.char = "")
  mut_df <- as.data.frame(lapply(mut_df, function(x) as.numeric(as.character(x))),
                          row.names = rownames(mut_df))
  mut_mat <- as.matrix(mut_df)
  ncols <- ncol(mut_mat)
  mut_numeric <- mut_mat
  for (j in seq_len(ncols)) {
    mut_numeric[, j] <- ifelse(mut_mat[, j] == 1, j, 0)
  }
  colnames(mut_numeric) <- colnames(mut_mat)

  base_col <- "grey95"
  if (ncols <= 9) {
    mut_cols <- brewer.pal(max(3, ncols), "Set1")[1:ncols]
  } else if (ncols <= 12) {
    mut_cols <- brewer.pal(ncols, "Paired")
  } else {
    mut_cols <- rainbow(ncols, start = 0, end = 0.85)
  }
  heatmap_matrix <- mut_numeric
  heatmap_colours <- c(base_col, mut_cols)
  heatmap_breaks <- seq(-0.5, ncols + 0.5, length.out = ncols + 2)
}

# Colour nodes only if the requested column exists
 colour_by <- args$colour_nodes_by
if (!colour_by %in% colnames(metadata)) {
  message("Column '", colour_by, "' not found in metadata; node colouring disabled")
  colour_by <- NULL
}

output_png <- args$output
if (is.null(output_png)) {
  output_png <- paste0(args$prefix, "_tree_heatmap.png")
}

# Build the figure
plotTree(
  tree = args$tree,
  ladderise = "descending",
  heatmapData = heatmap_matrix,
  heatmap.colours = heatmap_colours,
  heatmapBreaks = heatmap_breaks,
  heatmap.legend = FALSE,
  infoFile = info_file,
  infoCols = info_cols,
  colourNodesBy = colour_by,
  outputPNG = output_png,
  w = args$width,
  h = args$height,
  tip.labels = FALSE,
  tip.colour.cex = 1.5,
  lwd = 2,
  infoCex = 1.5,
  colLabelCex = 1.2,
  legend = TRUE,
  legend.pos = "bottomleft",
  treeWidth = 20,
  infoWidth = 15,
  dataWidth = 12,
  mainHeight = 150,
  labelHeight = 12
)

message("Tree + heatmap figure written to: ", output_png)
