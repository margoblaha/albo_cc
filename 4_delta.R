# =============================================================================
# Script 4: Delta Maps (Future − Current Change in Egg Counts)
# =============================================================================
# Description:
#   Computes the difference (delta) between future and current mean spatial
#   predictions for each scenario (model × period combination) and produces
#   faceted maps showing projected changes in weekly mosquito egg counts for
#   10 selected weeks across the season.
#
#   The current baseline raster is loaded once; each future raster is cropped
#   to the common spatial extent before differencing to ensure grid alignment.
#
# Inputs:
#   - Mean_Predictions_Current_2005-2024_<date>.RDS       : Current-period
#     mean raster (Script 3a)
#   - Mean_Predictions_Rasters_<model>_<period>_<date>.RDS: Future-period
#     mean rasters (Script 3b)
#
# Outputs (per scenario):
#   - Delta_Predictions_<model>_<period>_<date>.png : Faceted delta map figure
#
# Dependencies: terra, sf, ggplot2, geodata, stars, ggspatial, dplyr, tidyr,
#               ISOweek, eurostat, scales
# =============================================================================

# ---- Libraries ----

library(terra)
library(sf)
library(ggplot2)
library(geodata)
library(stars)
library(ggspatial)
library(dplyr)
library(tidyr)
library(ISOweek)
library(eurostat)
library(scales)

# ---- Paths ----
# data_dir : folder containing mean prediction raster RDS files
# fig_dir  : destination for figure outputs
# geo_dir  : cache directory for geodata::gadm() downloads

data_dir <- "data/"      # <-- update to your local data directory
fig_dir  <- "figures/"   # <-- update to your local figures directory
geo_dir  <- "geodata/"   # <-- update to your preferred geodata cache path

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(geo_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Load current-period baseline raster ----

spat.pred.c <- terra::unwrap(readRDS(
  tail(sort(list.files(data_dir,
                       pattern = "^Mean_Predictions_Current_2005-2024_.*\\.RDS$",
                       full.names = TRUE)), 1)
))

# ---- List future mean prediction raster files ----

spat.pred.f.files <- list.files(data_dir,
                                 pattern = "^Mean_Predictions_Rasters_.*\\.RDS$",
                                 full.names = TRUE)

# Filter to expected model names
spat.pred.f.files <- spat.pred.f.files[
  grepl("MOHC|ICHEC|IPSL", spat.pred.f.files)
]

if (length(spat.pred.f.files) == 0) {
  stop("No future raster files found in: ", data_dir,
       "\nExpected pattern: Mean_Predictions_Rasters_<MOHC|ICHEC|IPSL>_*.RDS")
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
  filter(CNTR_CODE %in% c("IT", "FR", "CH", "LI", "DE", "AT", "SI", "HR",
                           "SK", "HU", "BA", "CZ"))

# Country-level borders
countries_all <- geodata::gadm(
  country = c("Italy", "France", "Switzerland", "Liechtenstein",
               "Germany", "Austria", "Slovenia", "Croatia",
               "Slovakia", "Hungary", "Bosnia and Herzegovina", "Czech Republic"),
  level = 0, path = geo_dir
) %>% st_as_sf()

box_sfc <- st_as_sfc(
  st_bbox(c(xmin = 5.07, xmax = 17.8, ymin = 43, ymax = 48.82),
          crs = st_crs(countries_all))
)

# Clip borders to bounding box to avoid edge artefacts
countries_lines <- st_intersection(st_boundary(countries_all), box_sfc)
nuts2_lines     <- st_intersection(st_boundary(nuts2_sub),     box_sfc)

# ---- Compute and map delta for each future scenario ----

for (f in spat.pred.f.files) {

  message("Processing: ", basename(f))

  file_info     <- strsplit(basename(f), "_")[[1]]
  model         <- file_info[4]
  years         <- file_info[5]

  spat.pred.fut <- readRDS(f)

  # Crop both rasters to their common spatial extent
  common.extent  <- terra::intersect(ext(spat.pred.c), ext(spat.pred.fut))
  current.cropped <- terra::crop(spat.pred.c,   common.extent)
  future.cropped  <- terra::crop(spat.pred.fut, common.extent)

  delta <- future.cropped - current.cropped

  # Extract selected bands and convert to stars
  stars.pred        <- delta[[week_codes$isoweek]]
  names(stars.pred) <- LETTERS[seq_len(nlyr(stars.pred))]
  stars.pred        <- st_as_stars(stars.pred)

  # Faceted delta map with diverging colour scale
  delta.plot <- ggplot() +
    geom_stars(data = stars.pred, inherit.aes = FALSE) +
    geom_sf(data = countries_lines, color = "black", linewidth = 0.3) +
    geom_sf(data = nuts2_lines,     color = "black", linewidth = 0.18, alpha = 0.4) +
    ggspatial::annotation_scale(
      data       = data.frame(band = "J", location = "br"),
      aes(location = location),
      width_hint = 0.25
    ) +
    scale_fill_stepsn(
      colours = c(
        "#ABD8DE",  # (−25, −10]
        "#C6EDE3",  # (−10,  −5]
        "#FFFFFF",  # ( −5,   5]
        "#FEEDDE",  # (  5,  10]
        "#FDD0A2",  # ( 10,  25]
        "#FDAE6B",  # ( 25,  50]
        "#FD8D3C",  # ( 50,  75]
        "#F16913",  # ( 75, 100]
        "#D94801",  # (100, 200]
        "#8C2D04"   # (200, 300]
      ),
      values = scales::rescale(c(-25, -10, -5, 5, 10, 25, 50, 75, 100, 200, 300)),
      breaks = c(-25, -10, -5, 5, 10, 25, 50, 75, 100, 200, 300),
      limits = c(-25, 300),
      na.value = "transparent",
      name = "\u0394 mean number of eggs"
    ) +
    labs(x = "Longitude", y = "Latitude") +
    scale_y_continuous(limits = c(43, 49)) +
    scale_x_continuous(limits = c(5, 18), breaks = seq(5, 18, by = 3)) +
    facet_wrap(~band, ncol = 5, labeller = labeller(band = labeller_vector)) +
    theme_classic() +
    theme(
      legend.position   = "bottom",
      text              = element_text(size = 12),
      legend.text       = element_text(size = 10),
      legend.key.width  = unit(1.5, "cm"),
      legend.key.height = unit(0.35, "cm")
    )

  ggsave(
    file.path(fig_dir, paste0("Delta_Predictions_", model, "_", years, "_",
                               Sys.Date(), ".png")),
    delta.plot, width = 3200, height = 1300, units = "px", dpi = 300
  )

  rm(spat.pred.fut, current.cropped, future.cropped, delta, stars.pred, delta.plot)
  gc()

  message("Saved delta map for: ", model, " | ", years)
}
