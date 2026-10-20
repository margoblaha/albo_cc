# =============================================================================
# Script 2b: Convert Adjusted Climate Rasters to Daily Covariate Files
# =============================================================================
# Description:
#   Reads the adjusted daily climate rasters produced by Script 2a (one pair
#   of temperature and precipitation .tif files per model × period) and
#   converts them to long-format daily data frames, one file per scenario.
#
#   Each output file contains one row per grid cell per day with columns:
#   lon, lat, label (layer name), date, temperature, precipitation.
#
#   These daily files feed into Script 2c for weekly aggregation, photoperiod
#   computation, and lag generation.
#
# Inputs (from data_dir):
#   - <Model>_<index>_temp_<period>.tif : Adjusted daily temperature rasters
#   - <Model>_<index>_prec_<period>.tif : Adjusted daily precipitation rasters
#     (both produced by Script 2a)
#
# Outputs (written to out_dir):
#   - daily_covariates_<period>_<model>.RDS : Long-format daily data frame
#     per model × period combination
#
# Dependencies: terra, dplyr, tidyr, stringr, lubridate, purrr
# =============================================================================

# ---- Libraries ----

library(terra)
library(dplyr)
library(tidyr)
library(stringr)
library(lubridate)
library(purrr)

# ---- Paths ----
# data_dir : folder containing adjusted .tif rasters from Script 2a
# out_dir  : destination for daily covariate RDS output files

data_dir <- "data/"   # <-- update to your local data directory
out_dir  <- "data/"   # <-- update to your desired output directory

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Discover and parse available raster files ----

all_files <- list.files(data_dir, pattern = "\\.tif$",
                         full.names = FALSE, recursive = FALSE) %>%
  keep(~ !str_detect(.x, "\\.aux\\.json$|2005-2024"))  # exclude baseline rasters

if (length(all_files) == 0) {
  stop("No adjusted raster .tif files found in: ", data_dir,
       "\nRun Script 2a first to generate change-factor-adjusted rasters.")
}

# Parse model, scenario index, variable, and period from filename
# Expected filename format: <Model>_<index>_<var>_<period>.tif
info <- tibble(file = all_files) %>%
  mutate(
    parts  = str_split(file, "_"),
    model  = map_chr(parts, 1),
    var    = map_chr(parts, 3),
    period = str_extract(file, "(?<=_)[0-9]{4}-[0-9]{4}(?=\\.tif)")
  ) %>%
  select(-parts)

models  <- unique(info$model)
periods <- unique(info$period)

# ---- Helper: convert raster to long-format data frame ----

raster_to_long <- function(r, var_prefix) {
  # r          : SpatRaster with one layer per day and terra time attributes
  # var_prefix : column prefix for the raster variable (e.g., "tg_" or "rr_")
  #
  # Returns a data frame with columns: lon, lat, label, date, <var_prefix>

  dates  <- as.Date(terra::time(r))
  labels <- names(r)

  date_df <- tibble(label = labels, date = dates)
  values  <- as.data.frame(r, xy = TRUE)

  values %>%
    pivot_longer(
      cols      = starts_with(var_prefix),
      names_to  = "label",
      values_to = var_prefix
    ) %>%
    left_join(date_df, by = "label") %>%
    rename(lon = x, lat = y)
}

# ---- Convert rasters to daily data frames ----

for (mod in models) {
  for (per in periods) {
    message("Processing model: ", mod, " | period: ", per)

    subset_info <- filter(info, model == mod, period == per)

    fn_temp <- filter(subset_info, var == "temp")$file
    fn_prec <- filter(subset_info, var == "prec")$file

    if (length(fn_temp) == 0 || length(fn_prec) == 0) {
      warning("Missing temp or prec file for ", mod, " | ", per, " — skipping.")
      next
    }

    r_temp <- rast(file.path(data_dir, fn_temp))
    r_prec <- rast(file.path(data_dir, fn_prec))

    long_temp <- raster_to_long(r_temp, "tg_")
    long_prec <- raster_to_long(r_prec, "rr_") %>% select(rr_)

    covs <- cbind(long_temp, long_prec) %>%
      rename(temperature   = tg_,
             precipitation = rr_)

    out_name <- file.path(out_dir,
                           sprintf("daily_covariates_%s_%s.RDS", per, mod))
    saveRDS(covs, out_name)

    rm(r_temp, r_prec, long_temp, long_prec, covs)
    gc()

    message("Saved: ", basename(out_name))
  }
}
