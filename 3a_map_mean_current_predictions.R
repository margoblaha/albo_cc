# =============================================================================
# Script 3a: Map Mean Current-Period Predictions (2005–2024)
# =============================================================================
# Description:
#   Aggregates weekly spatial predictions from the current climate period
#   (2005–2024) into mean (and uncertainty) estimates per grid cell and ISO
#   week. The aggregated means are converted to a multi-band raster and
#   plotted as a faceted map showing 10 selected weeks across the season.
#
#   Country and NUTS-2 region boundaries are overlaid for context.
#
# Inputs:
#   - predictions_current_scenario_<date>.RDS : Spatial predictions (Script 1)
#
# Outputs:
#   - mean_predictions_current_scenario_<date>.RDS    : Aggregated predictions
#     per grid cell × week (data frame)
#   - Mean_Predictions_Current_2005-2024_<date>.RDS   : Multi-band SpatRaster
#     of mean predictions
#   - Mean_Predictions_current_2005-2024_<date>.png   : Faceted map figure
#
# Dependencies: sf, dplyr, ggplot2, tidyr, stars, terra, stringr, lubridate,
#               ISOweek, eurostat, geodata, ggspatial, scales
# =============================================================================

# ---- Libraries ----

library(sf)
library(dplyr)
library(ggplot2)
library(tidyr)
library(stars)
library(terra)
library(stringr)
library(lubridate)
library(ISOweek)
library(eurostat)
library(geodata)
library(ggspatial)
library(scales)

# ---- Paths ----
# data_dir : folder containing input prediction RDS files and for saving outputs
# fig_dir  : destination for figure outputs
# geo_dir  : cache directory for geodata::gadm() downloads

setwd("Dati/albo_cc/")      # <-- update to your local data directory
geo_dir  <- "/home/geodata/"   # <-- update to your preferred geodata cache path

#dir.create(fig_dir,  showWarnings = FALSE, recursive = TRUE)
#dir.create(geo_dir,  showWarnings = FALSE, recursive = TRUE)

# ---- 1. Load and aggregate spatial predictions ----

my_spPred_all_df <- readRDS("predictions_current_scenario_2026-07-14.RDS")

# Compute per-cell, per-week summary statistics across years
mean_preds <- my_spPred_all_df %>%
  group_by(ID, lon, lat, week) %>%
  summarise(
    mean_pred   = mean(spatial_pred, na.rm = TRUE),
    median_pred = median(spatial_pred, na.rm = TRUE),
    pred_se     = sd(spatial_pred, na.rm = TRUE) / sqrt(n()),
    .groups     = "drop"
  ) %>%
  mutate(
    pred_lci = mean_pred - (1.96 * pred_se),
    pred_uci = mean_pred + (1.96 * pred_se)
  )

saveRDS(mean_preds,
        file.path(paste0("mean_predictions_current_scenario_", Sys.Date(), ".RDS")))

mean_preds <- readRDS("~/Dati/albo_cc/mean_predictions_current_scenario_2025-08-19.RDS")

# ---- 2.1 Convert to multi-band raster ----

mean_wide <- mean_preds %>%
  mutate(week = paste0("w_", str_pad(week, 2, pad = "0"))) %>%
  select(lon, lat, week, mean_pred) %>%
  pivot_wider(names_from = "week", values_from = "mean_pred")

mean_rast              <- terra::rast(mean_wide, type = "xyz")
mean_rast[mean_rast < 0] <- 0  # clip negative predictions to zero

saveRDS(mean_rast,
        file.path(paste0("Mean_Predictions_Current_2005-2024_", Sys.Date(), ".RDS")))

# ---- 2.2 Create LCI and UCI rasters ----

lci_wide <- mean_preds %>%
  mutate(week = paste0("w_", str_pad(week, 2, pad = "0"))) %>%
  select(lon, lat, week, pred_lci) %>%
  pivot_wider(names_from = "week", values_from = "pred_lci")

lci_rast              <- terra::rast(lci_wide, type = "xyz")
lci_rast[lci_rast < 0] <- 0  # clip negative predictions to zero

saveRDS(lci_rast,
        file.path(paste0("LCI_Predictions_Current_2005-2024_", Sys.Date(), ".RDS")))

uci_wide <- mean_preds %>%
  mutate(week = paste0("w_", str_pad(week, 2, pad = "0"))) %>%
  select(lon, lat, week, pred_uci) %>%
  pivot_wider(names_from = "week", values_from = "pred_uci")

uci_rast              <- terra::rast(uci_wide, type = "xyz")
uci_rast[uci_rast < 0] <- 0  # clip negative predictions to zero

saveRDS(uci_rast,
        file.path(paste0("UCI_Predictions_Current_2005-2024_", Sys.Date(), ".RDS")))

# ---- 3. Build faceted maps for Mean, LCI, and UCI ----

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
  year         = 2024,
  cache = FALSE
)

# Country-level borders
countries_all <- geodata::gadm(
  country = c("Italy", "France", "Switzerland", "Liechtenstein",
              "Germany", "Austria", "Slovenia", "Croatia",
              "Slovakia", "Hungary", "Bosnia and Herzegovina", "Czech Republic"),
  level = 0, path = "~/geodata/"
) %>% st_as_sf()

box_sfc <- st_as_sfc(
  st_bbox(c(xmin = 5.07, xmax = 17.8, ymin = 43, ymax = 48.82),
          crs = st_crs(countries_all))
)

# Clip borders to bounding box to avoid edge artefacts
countries_lines <- st_intersection(st_boundary(countries_all), box_sfc)
nuts2_lines     <- st_intersection(st_boundary(nuts2_sub),     box_sfc)

# ---- Loop over the three prediction rasters ----

raster_list <- list(
  Mean = list(rast = mean_rast, label = "Mean number of eggs",       out_prefix = "Mean_Predictions_current_2005-2024"),
  LCI  = list(rast = lci_rast,  label = "Lower CI (number of eggs)", out_prefix = "LCI_Predictions_current_2005-2024"),
  UCI  = list(rast = uci_rast,  label = "Upper CI (number of eggs)", out_prefix = "UCI_Predictions_current_2005-2024")
)

for (nm in names(raster_list)) {
  
  entry <- raster_list[[nm]]
  
  # Extract selected bands and convert to stars
  stars_pred        <- entry$rast[[week_codes$isoweek]]
  names(stars_pred) <- LETTERS[seq_len(nlyr(stars_pred))]
  stars_pred        <- st_as_stars(stars_pred)
  
  p <- ggplot() +
    geom_stars(data = stars_pred, inherit.aes = FALSE) +
    geom_sf(data = countries_lines, color = "black", linewidth = 0.25, alpha = 0.7) +
    geom_sf(data = nuts2_lines,     color = "black", linewidth = 0.18, alpha = 0.3) +
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
    labs(x = "Longitude", y = "Latitude", fill = entry$label) +
    theme_classic() +
    theme(
      legend.position   = "bottom",
      text              = element_text(size = 12),
      legend.text       = element_text(size = 10),
      legend.key.width  = unit(1.5, "cm"),
      legend.key.height = unit(0.35, "cm")
    )
  
  ggsave(
    file.path(paste0("Figures/", entry$out_prefix, "_", Sys.Date(), ".png")),
    p, width = 3200, height = 1300, units = "px", dpi = 300
  )
}

# old version
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

# Extract selected bands and convert to stars
stars_pred              <- mean_rast[[week_codes$isoweek]]
names(stars_pred)       <- LETTERS[seq_len(nlyr(stars_pred))]
stars_pred              <- st_as_stars(stars_pred)

# NUTS-2 sub-national boundaries
nuts2_sub <- eurostat::get_eurostat_geospatial(
  output_class = "sf",
  resolution   = "01",
  nuts_level   = 2,
  year         = 2024,
  cache = FALSE
)

# nuts2_sub <- gisco_get_nuts(
#   year = "2024",
#   nuts_level = 2,
#   resolution = "01"
# )

# %>%
#   filter(CNTR_CODE %in% c("IT", "FR", "CH", "LI", "DE", "AT", "SI", "HR",
#                            "SK", "HU", "BA", "CZ"))

# Country-level borders
countries_all <- geodata::gadm(
  country = c("Italy", "France", "Switzerland", "Liechtenstein",
               "Germany", "Austria", "Slovenia", "Croatia",
               "Slovakia", "Hungary", "Bosnia and Herzegovina", "Czech Republic"),
  level = 0, path = "~/geodata/"
) %>% st_as_sf()

box_sfc <- st_as_sfc(
  st_bbox(c(xmin = 5.07, xmax = 17.8, ymin = 43, ymax = 48.82),
          crs = st_crs(countries_all))
)

# Clip borders to bounding box to avoid edge artefacts
countries_lines <- st_intersection(st_boundary(countries_all), box_sfc)
nuts2_lines     <- st_intersection(st_boundary(nuts2_sub),     box_sfc)

# Faceted map
albo <- ggplot() +
  geom_stars(data = stars_pred, inherit.aes = FALSE) +
  geom_sf(data = countries_lines, color = "black", linewidth = 0.25, alpha = 0.7) +
  geom_sf(data = nuts2_lines,     color = "black", linewidth = 0.18, alpha = 0.3) +
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
    legend.position  = "bottom",
    text             = element_text(size = 12),
    legend.text      = element_text(size = 10),
    legend.key.width  = unit(1.5, "cm"),
    legend.key.height = unit(0.35, "cm")
  )

ggsave(
  file.path(paste0("Figures/Mean_Predictions_current_2005-2024_", Sys.Date(), ".png")),
  albo, width = 3200, height = 1300, units = "px", dpi = 300
)
