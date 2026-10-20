# =============================================================================
# Script 2a: Climate Change Factor Application (Delta Method)
# =============================================================================
# Description:
#   Applies monthly change factors (delta method) derived from EURO-CORDEX
#   regional climate model (RCM) ensembles to E-OBS observed daily
#   temperature and precipitation grids. Adjusted climate fields are saved
#   as GeoTIFF rasters for use in downstream spatial covariate generation.
#
#   Three RCMs are processed:
#     1. MOHC-HadGEM2-ES   – represents average expected changes
#     2. IPSL-IPSL-CM5A-MR – warmer winters, wetter summers
#     3. ICHEC-EC-EARTH     – much warmer and drier summers
#
#   Two future periods are produced per model: 2036–2055 and 2066–2085.
#   Temperature is adjusted by addition; precipitation by multiplication.
#
# Inputs:
#   - E-OBS NetCDF ensemble files (precipitation "rr", temperature "tg")
#     matching pattern "ens" in data_dir
#   - Monthly change-factor NetCDF files in cf_dir/<var>_<period>/
#     (filenames contain EUR-11_<ModelName>)
#
# Outputs (written to data_dir):
#   - <Model>_<index>_prec_<period>.tif  : Adjusted precipitation rasters
#   - <Model>_<index>_temp_<period>.tif  : Adjusted temperature rasters
#   - prec_2005-2024.tif                 : Baseline E-OBS precipitation
#   - temp_2005-2024.tif                 : Baseline E-OBS temperature
#
# Dependencies: terra, dplyr, tidyr
# =============================================================================

# ---- Libraries ----

library(terra)
library(dplyr)
library(tidyr)

# ---- Paths ----
# data_dir : folder containing E-OBS ensemble NetCDF files
# cf_dir   : root folder for change-factor NetCDF files (subfolders per variable/period)
# out_dir  : destination for output GeoTIFF rasters

data_dir <- "data/"          # <-- update to your local data directory
cf_dir   <- "data/change_factors/"  # <-- update to change-factor root directory
out_dir  <- "data/"          # <-- update to your desired output directory

# ---- Load E-OBS baseline grids ----

file_list <- list.files(data_dir, pattern = "ens", full.names = TRUE)

covs <- lapply(file_list, terra::rast)

# Spatial bounding box covering the study domain
bbox <- ext(2, 20, 42, 50)

prec <- terra::rast(covs[[1]])
temp <- terra::rast(covs[[2]])

# Subset to the study period starting from 2005
time_vector    <- terra::time(temp)
from_2005_idx  <- which(time_vector >= as.Date("2005-01-01"))

temp <- temp[[from_2005_idx]]
prec <- prec[[from_2005_idx]]

# Crop to study domain
temp <- crop(temp, bbox)
prec <- crop(prec, bbox)

# Save baseline rasters for reference
writeRaster(prec, file.path(out_dir, "prec_2005-2024.tif"), overwrite = TRUE)
writeRaster(temp, file.path(out_dir, "temp_2005-2024.tif"), overwrite = TRUE)

# Extract month index for each layer (used to select the correct monthly change factor)
time_vector <- terra::time(temp)
months      <- as.integer(format(time_vector, "%m"))

# ---- Define models and future periods ----

# Model name as it appears in change-factor filenames (EUR-11_<ModelName>_...)
model_centers <- list(
  "MOHC-HadGEM2-ES"   = 1,  # average expected changes
  "IPSL-IPSL-CM5A-MR" = 2,  # warmer winters, wetter summers
  "ICHEC-EC-EARTH"    = 3   # much warmer and drier summers
)

future_periods <- c("2036-2055", "2066-2085")

# ---- Apply change factors ----

for (model in names(model_centers)) {

  for (period in future_periods) {

    prec_files <- list.files(file.path(cf_dir, paste0("pr_",  period)),
                             pattern = "\\.nc$", full.names = TRUE)
    temp_files <- list.files(file.path(cf_dir, paste0("tas_", period)),
                             pattern = "\\.nc$", full.names = TRUE)

    # --- Precipitation (multiplicative adjustment) ---
    for (prec_file in prec_files) {
      model_name <- gsub(".*EUR-11_([^_]+)_.*", "\\1", prec_file)
      if (model_name != model) next

      prec_cf           <- rast(prec_file)
      prec_cf_resampled <- resample(prec_cf, prec[[1]])

      adjusted_prec <- vector("list", nlyr(prec))
      for (i in seq_len(nlyr(prec))) {
        month_index        <- as.integer(months[i])
        adjusted_prec[[i]] <- prec[[i]] * prec_cf_resampled[[month_index]]
      }

      adjusted_prec <- do.call(c, adjusted_prec)
      writeRaster(
        adjusted_prec,
        file.path(out_dir, paste0(model, "_", model_centers[[model]],
                                  "_prec_", period, ".tif")),
        overwrite = TRUE
      )
      rm(adjusted_prec)
      gc()
    }

    # --- Temperature (additive adjustment) ---
    for (temp_file in temp_files) {
      model_name <- gsub(".*EUR-11_([^_]+)_.*", "\\1", temp_file)
      if (model_name != model) next

      temp_cf           <- rast(temp_file)
      temp_cf_resampled <- resample(temp_cf, temp[[1]])

      adjusted_temp <- vector("list", nlyr(temp))
      for (i in seq_len(nlyr(temp))) {
        month_index        <- as.integer(months[i])
        adjusted_temp[[i]] <- temp[[i]] + temp_cf_resampled[[month_index]]
      }

      adjusted_temp <- do.call(c, adjusted_temp)
      writeRaster(
        adjusted_temp,
        file.path(out_dir, paste0(model, "_", model_centers[[model]],
                                  "_temp_", period, ".tif")),
        overwrite = TRUE
      )
      rm(adjusted_temp)
      gc()
    }

    message("Completed: ", model, " | ", period)
  }
}
