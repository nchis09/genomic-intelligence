#!/usr/bin/env Rscript
#
# plot_geo_map.R
#
# Static geographic visualisation of the samples represented in a phylogenetic
# tree. Given a tip-metadata TSV (first column = tip label) it plots the
# country-level distribution of samples on a world map using ggplot2 + maps.
#
# Optional: provide a lat/long TSV with columns location, latitude, longitude
# to override the automatic centroid lookup.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(maps)
  library(optparse)
  library(RColorBrewer)
})

option_list <- list(
  make_option(c("-m", "--metadata"), type = "character", default = NULL,
              help = "Tip metadata TSV (first column = tip label)"),
  make_option(c("-l", "--lat-longs"), type = "character", default = NULL,
              dest = "lat_longs",
              help = "Optional lat/long TSV: location\\tlatitude\\tlongitude"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Output PNG file"),
  make_option(c("-p", "--prefix"), type = "character", default = "out",
              help = "Output prefix if --output not given"),
  make_option(c("-c", "--colour-by"), type = "character", default = "country",
              dest = "colour_by",
              help = "Metadata column used to colour points"),
  make_option(c("--width"), type = "integer", default = 2000,
              help = "PNG width in pixels"),
  make_option(c("--height"), type = "integer", default = 1400,
              help = "PNG height in pixels")
)

parser <- OptionParser(option_list = option_list,
                       description = "Static geographic map of tree tips.")
args <- parse_args(parser)

if (is.null(args$metadata)) {
  stop("--metadata is required")
}

output_png <- args$output
if (is.null(output_png)) {
  output_png <- paste0(args$prefix, "_geo_map.png")
}

metadata <- read.table(args$metadata,
                       header = TRUE,
                       sep = "\t",
                       check.names = FALSE,
                       stringsAsFactors = FALSE,
                       row.names = 1,
                       quote = "",
                       comment.char = "",
                       fill = TRUE)

if (!args$colour_by %in% colnames(metadata)) {
  stop("Colour column '", args$colour_by, "' not found in metadata")
}

# Keep rows that have a non-empty geographic value
metadata <- metadata[metadata[[args$colour_by]] != "" & !is.na(metadata[[args$colour_by]]), , drop = FALSE]

if (nrow(metadata) == 0) {
  stop("No geographic information found in metadata")
}

# Summarise sample counts per location
loc_counts <- metadata %>%
  group_by(across(all_of(args$colour_by))) %>%
  summarise(n = n(), .groups = "drop") %>%
  rename(location = !!args$colour_by)

# ---------------------------------------------------------------------------
# Build / load coordinates
# ---------------------------------------------------------------------------
if (!is.null(args$lat_longs) && file.exists(args$lat_longs)) {
  coords <- read.table(args$lat_longs,
                       header = TRUE,
                       sep = "\t",
                       check.names = FALSE,
                       stringsAsFactors = FALSE)
  required <- c("location", "latitude", "longitude")
  missing <- setdiff(required, colnames(coords))
  if (length(missing) > 0) {
    stop("lat/long file must contain columns: location, latitude, longitude")
  }
  coords <- coords %>%
    select(location, lat = latitude, lon = longitude) %>%
    mutate(across(c(lat, lon), as.numeric))
} else {
  # Use maps::world.cities for country centroids (mean lat/lon per country)
  data("world.cities", package = "maps", envir = environment())
  cities <- world.cities
  coords <- as.data.frame(cities) %>%
    group_by(country = country.etc) %>%
    summarise(lat = mean(lat, na.rm = TRUE),
              lon = mean(long, na.rm = TRUE),
              .groups = "drop") %>%
    rename(location = country)

  unmatched <- setdiff(loc_counts$location, coords$location)
  if (length(unmatched) > 0) {
    warning("Could not map the following locations to centroids: ",
            paste(unmatched, collapse = ", "))
  }
}

loc_points <- loc_counts %>%
  left_join(coords, by = "location") %>%
  filter(!is.na(lat), !is.na(lon))

if (nrow(loc_points) == 0) {
  stop("No coordinates available for the locations in the metadata")
}

# ---------------------------------------------------------------------------
# Zoom the map to the region containing the samples
# ---------------------------------------------------------------------------
lon_range <- range(loc_points$lon, na.rm = TRUE)
lat_range <- range(loc_points$lat, na.rm = TRUE)
x_pad <- max(5, 0.25 * diff(lon_range))
y_pad <- max(5, 0.25 * diff(lat_range))
xlim <- c(max(-180, lon_range[1] - x_pad), min(180, lon_range[2] + x_pad))
ylim <- c(max(-90, lat_range[1] - y_pad), min(90, lat_range[2] + y_pad))

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
n_groups <- length(unique(loc_points$location))
col_palette <- if (n_groups <= 8) {
  "Set1"
} else if (n_groups <= 12) {
  "Set3"
} else {
  "Paired"
}
custom_cols <- colorRampPalette(brewer.pal(8, col_palette))(n_groups)

world_map <- map_data("world")

png(output_png, width = args$width, height = args$height, res = 200)
p <- ggplot() +
  geom_polygon(data = world_map,
               aes(x = long, y = lat, group = group),
               fill = "grey90", color = "white", linewidth = 0.2) +
  geom_point(data = loc_points,
             aes(x = lon, y = lat, color = location, size = n),
             alpha = 0.85) +
  scale_color_manual(values = custom_cols, name = "Country") +
  scale_size_continuous(range = c(3, 16), name = "Samples") +
  coord_fixed(ratio = 1.3, xlim = xlim, ylim = ylim, expand = FALSE) +
  labs(title = paste(args$prefix, "geographic distribution"),
       x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10)
  )
print(p)
dev.off()

message("Geographic map written to: ", output_png)
