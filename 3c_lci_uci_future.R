# =============================================================================
# Script 3c: Map LCI and UCI Future-Period Predictions
# =============================================================================
# Description:
#   For each future scenario's aggregated prediction file (Script 3b output,
#   which already contains pred_lci/pred_uci per cell x week), builds LCI and
#   UCI multi-band rasters and produces faceted maps showing 10 selected
#   weeks across the season — mirroring the mean-prediction maps from 3b.
#
# Inputs:
#   - Mean_Predictions_<model>_<period>_<date>.RDS : Aggregated data frame
#     with mean_pred, pred_se, pred_lci, pred_uci (Script 3b)
#
# Outputs (per scenario):
#   - LCI_Predictions_Rasters_<model>_<period>_<date>.RDS
#   - UCI_Predictions_Rasters_<model>_<period>_<date>.RDS
#   - LCI_Predictions_<model>_<period>_<date>.png
#   - UCI_Predictions_<model>_<period>_<date>.png
#
# Dependencies: terra, sf, tidyverse, stars, geodata, ggplot2, dplyr,
#               lubridate, ISOweek, eurostat, ggspatial, scales, stringr
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
library(stringr)

# ---- Paths ----
# data_dir : folder containing Mean_Predictions_*.RDS files (Script 3b output) and for saving outputs
# fig_dir  : destination for figure outputs

data_dir <- "~/Dati/albo_cc/"      # <-- update to your local data directory
fig_dir  <- "Figures/"   # <-- update to your local figures directory

# ---- List aggregated scenario prediction files (Script 3b output) ----
# Excludes the "..._Rasters_..." files also written by 3b

mpreds <- list.files(data_dir,
                     pattern = "^Mean_Predictions_.*\\.RDS$",
                     full.names = TRUE)
mpreds <- mpreds[grepl("2025-08-19.*\\.RDS$", basename(mpreds))]

mpreds <- mpreds[!grepl("Rasters", basename(mpreds))]

mpreds <- mpreds[grep("ICHEC|IPSL|MOHC", basename(mpreds))]

if (length(mpreds) == 0) {
  stop("No aggregated scenario prediction files found in: ", data_dir,
       "\nExpected pattern: Mean_Predictions_<model>_<period>_<date>.RDS")
}

# ---- Build shared map elements (identical to Script 3b) ----

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
  level = 0, path = "~/geodata/"
) %>% st_as_sf()

box_sfc <- st_as_sfc(
  st_bbox(c(xmin = 5.07, xmax = 18, ymin = 43, ymax = 48.82),
          crs = st_crs(countries_all))
)

# Clip borders to bounding box to avoid edge artefacts
countries_lines <- st_intersection(st_boundary(countries_all), box_sfc)
nuts2_lines     <- st_intersection(st_boundary(nuts2_sub),     box_sfc)

# ---- Which CI bound to process, and how to label/name each ----

ci_specs <- list(
  LCI = list(col = "pred_lci", label = "Lower 95% CI (number of eggs)", prefix = "LCI_Predictions"),
  UCI = list(col = "pred_uci", label = "Upper 95% CI (number of eggs)", prefix = "UCI_Predictions")
)

# ---- Process each future scenario ----

for (file_path in mpreds) {
  
  message("Processing: ", basename(file_path))
  
  fname <- basename(file_path)
  years <- str_extract(fname, "\\d{4}-\\d{4}")
  model <- str_extract(fname, "(?<=Mean_Predictions_)[A-Za-z0-9\\-]+(?=_\\d{4}-\\d{4})")
  
  mean_fpreds <- readRDS(file_path)
  
  for (ci_name in names(ci_specs)) {
    
    spec <- ci_specs[[ci_name]]
    
    # Convert to multi-band raster
    ci_wide <- mean_fpreds %>%
      ungroup() %>%
      select(lon, lat, week, all_of(spec$col)) %>%
      pivot_wider(names_from = "week", values_from = spec$col)
    
    ci_rast              <- terra::rast(ci_wide, type = "xyz")
    ci_rast[ci_rast < 0] <- 0  # clip negative predictions to zero
    
    saveRDS(ci_rast,
            file.path(data_dir,
                      paste0(spec$prefix, "_Rasters_", model, "_", years, "_",
                             Sys.Date(), ".RDS")))
    
    # Extract selected bands and convert to stars
    stars_pred        <- ci_rast[[week_codes$isoweek]]
    names(stars_pred) <- LETTERS[seq_len(nlyr(stars_pred))]
    stars_pred        <- st_as_stars(stars_pred)
    
    # Faceted map
    p <- ggplot() +
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
      labs(x = "Longitude", y = "Latitude", fill = spec$label) +
      theme_classic() +
      theme(
        legend.position   = "bottom",
        text              = element_text(size = 12),
        legend.text       = element_text(size = 10),
        legend.key.width  = unit(1.5, "cm"),
        legend.key.height = unit(0.35, "cm")
      )
    
    ggsave(
      file.path(fig_dir, paste0(spec$prefix, "_", model, "_", years, "_",
                                Sys.Date(), ".png")),
      p, width = 3200, height = 1300, units = "px", dpi = 300
    )
    
    message("  Saved ", ci_name, " for: ", model, " | ", years)
  }
  
  rm(mean_fpreds)
  gc()
}
