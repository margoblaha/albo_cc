# =============================================================================
# Script 3b: Map Mean Future-Period Predictions
# =============================================================================
# Description:
#   For each future scenario prediction file (Script 2b), aggregates weekly
#   spatial predictions into per-cell, per-week mean estimates, converts
#   them to multi-band rasters, and produces faceted maps showing 10 selected
#   weeks across the season.
#
#   The loop iterates over all future prediction RDS files found in data_dir,
#   producing one aggregated data frame, one raster object, and one figure
#   per scenario.
#
# Inputs:
#   - Predictions_<model>_<period>_<date>.RDS : Future spatial predictions
#     (Script 2b)
#
# Outputs (per scenario):
#   - Mean_Predictions_<model>_<period>_<date>.RDS        : Aggregated data frame
#   - Mean_Predictions_Rasters_<model>_<period>_<date>.RDS: Multi-band SpatRaster
#   - Mean_Predictions_<model>_<period>_<date>.png        : Faceted map figure
#
# Dependencies: terra, sf, tidyverse, stars, geodata, ggplot2, dplyr,
#               lubridate, ISOweek, eurostat, ggspatial, scales
# =============================================================================

# ---- Libraries ----

library(terra)
library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stars)
library(geodata)
library(lubridate)
library(ISOweek)
library(eurostat)
library(ggspatial)
library(scales)

# ---- Paths ----
# data_dir : folder containing future prediction RDS files and for saving outputs
# fig_dir  : destination for figure outputs
# geo_dir  : cache directory for geodata::gadm() downloads

data_dir <- "data/"      # <-- update to your local data directory
fig_dir  <- "figures/"   # <-- update to your local figures directory
geo_dir  <- "geodata/"   # <-- update to your preferred geodata cache path

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(geo_dir, showWarnings = FALSE, recursive = TRUE)

# ---- List future prediction files ----

fpreds <- list.files(data_dir,
                     pattern = "^Predictions_.*\\.RDS$",
                     full.names = TRUE)

if (length(fpreds) == 0) {
  stop("No future prediction files found in: ", data_dir,
       "\nExpected pattern: Predictions_<model>_<period>_<date>.RDS")
}

# ---- Build shared map elements ----

# ISO weeks to display and their letter codes used as facet labels
week_codes <- data.frame(
  band    = LETTERS[1:10],
  isoweek = c(10, 15, 18, 20, 25, 30, 35, 40, 43, 45)
)

ref_year <- 2015  # reference year for computing calendar date ranges per week

week_df <- week_codes %>%
  mutate(
    isoweek_str = sprintf("%d-W%02d", ref_year, isoweek),
    start_date  = ISOweek2date(paste0(isoweek_str, "-1")),  # Monday
    end_date    = ISOweek2date(paste0(isoweek_str, "-7")),  # Sunday
    label       = paste0("Week ", isoweek, " (",
                          format(start_date, "%d %b"), " – ",
                          format(end_date,   "%d %b"), ")")
  )

labeller_vector <- setNames(week_df$label, week_df$band)

# NUTS-2 sub-national boundaries
nuts2_sub <- eurostat::get_eurostat_geospatial(
  output_class = "sf",
  resolution   = "01",
  nuts_level   = 2,
  year         = 2024
) %>%
  filter(CNTR_CODE %in% c("IT", "FR", "CH", "LI", "DE", "AT", "SI", "HR"))

# Country-level borders
countries_all <- geodata::gadm(
  country = c("Italy", "France", "Switzerland", "Liechtenstein",
               "Germany", "Austria", "Slovenia", "Croatia"),
  level = 0, path = geo_dir
) %>% st_as_sf()

box_sfc <- st_as_sfc(
  st_bbox(c(xmin = 5.07, xmax = 18, ymin = 43, ymax = 48.82),
          crs = st_crs(countries_all))
)

# Clip borders to bounding box to avoid edge artefacts
countries_lines <- st_intersection(st_boundary(countries_all), box_sfc)
nuts2_lines     <- st_intersection(st_boundary(nuts2_sub),     box_sfc)

# ---- Process each future scenario ----

for (file_path in fpreds) {

  message("Processing: ", basename(file_path))

  fpreds_data <- readRDS(file_path)

  file_info <- strsplit(basename(file_path), "_")[[1]]
  model     <- file_info[2]
  years     <- file_info[3]

  # Aggregate across years: per cell × week mean and uncertainty
  mean_fpreds <- fpreds_data %>%
    group_by(tileID, lon, lat, week) %>%
    summarise(
      mean_pred   = mean(spatial_pred,   na.rm = TRUE),
      median_pred = median(spatial_pred, na.rm = TRUE),
      pred_se     = sd(spatial_pred,     na.rm = TRUE) / sqrt(n()),
      pred_lci    = mean_pred - (1.96 * pred_se),
      pred_uci    = mean_pred + (1.96 * pred_se),
      .groups     = "drop"
    )

  saveRDS(mean_fpreds,
          file.path(data_dir,
                    paste0("Mean_Predictions_", model, "_", years, "_",
                           Sys.Date(), ".RDS")))

  # Convert to multi-band raster
  mean_wide  <- mean_fpreds %>%
    ungroup() %>%
    select(lon, lat, week, mean_pred) %>%
    pivot_wider(names_from = "week", values_from = "mean_pred")

  spat_pred              <- terra::rast(mean_wide, type = "xyz")
  spat_pred[spat_pred < 0] <- 0  # clip negative predictions to zero

  saveRDS(spat_pred,
          file.path(data_dir,
                    paste0("Mean_Predictions_Rasters_", model, "_", years, "_",
                           Sys.Date(), ".RDS")))

  # Extract selected bands and convert to stars
  stars_pred        <- spat_pred[[week_codes$isoweek]]
  names(stars_pred) <- LETTERS[seq_len(nlyr(stars_pred))]
  stars_pred        <- st_as_stars(stars_pred)

  # Faceted map
  albo <- ggplot() +
    geom_stars(data = stars_pred, inherit.aes = FALSE) +
    geom_sf(data = countries_lines, color = "black", linewidth = 0.3) +
    geom_sf(data = nuts2_lines,     color = "black", linewidth = 0.18, alpha = 0.4) +
    scale_fill_viridis_b(
      option    = "rocket",
      direction = -1,
      limits    = c(0, 500),
      breaks    = c(0, 10, 25, 50, 75, 100, 150, 200, 300, 400),
      oob       = scales::squish,
      na.value  = "transparent"
    ) +
    ggspatial::annotation_scale(
      data       = data.frame(band = "J", location = "br"),
      aes(location = location),
      width_hint = 0.25
    ) +
    scale_y_continuous(limits = c(43, 49)) +
    scale_x_continuous(limits = c(5, 18), breaks = seq(5, 18, by = 3)) +
    facet_wrap(~band, ncol = 5, labeller = labeller(band = labeller_vector)) +
    labs(x = "Longitude", y = "Latitude", fill = "Mean number of eggs") +
    theme_classic() +
    theme(
      legend.position   = "bottom",
      text              = element_text(size = 12),
      legend.text       = element_text(size = 10),
      legend.key.width  = unit(1.5, "cm"),
      legend.key.height = unit(0.35, "cm")
    )

  ggsave(
    file.path(fig_dir, paste0("Mean_Predictions_", model, "_", years, "_",
                               Sys.Date(), ".png")),
    albo, width = 3200, height = 1300, units = "px", dpi = 300
  )

  rm(fpreds_data, mean_fpreds, mean_wide, spat_pred, stars_pred, albo)
  gc()

  message("Saved figure for: ", model, " | ", years)
}
