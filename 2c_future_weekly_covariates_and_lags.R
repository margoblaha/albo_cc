# =============================================================================
# Script 2c: Weekly Aggregation, Photoperiod, and Lagged Covariates
#            for Future Climate Scenarios
# =============================================================================
# Description:
#   Processes the daily future covariate files (Script 2b) into weekly
#   aggregates with lagged climate predictors matching the format expected
#   by the trained model (Script 1). For each model × period scenario:
#
#   1. Aggregates daily temperature and precipitation to weekly medians/sums.
#   2. Shifts internal model years to the target future period by extracting
#      the period from the filename (e.g. 2036-2055 → year offset = +31 from
#      the 2005-2024 E-OBS baseline used in Script 2a).
#   3. Computes weekly median photoperiod (day length) per grid cell using
#      the geosphere package.
#   4. Computes rolling-mean lagged predictors (lags 2–10 weeks) for
#      temperature, precipitation, and photoperiod.
#   5. Applies a minimum temperature threshold (values below threshold set
#      to 0) to exclude biologically inactive periods.
#   6. Back-fills leading NA values introduced by the rolling-lag windows.
#
# Inputs (from data_dir):
#   - daily_covariates_<period>_<model>.RDS : Daily covariate files (Script 2b)
#
# Outputs (written to out_dir):
#   - covariates_<model>_<period>_with-lags_<date>.RDS : Weekly covariates
#     with all lagged predictors, ready for spatial prediction (Script 2d)
#
# Dependencies: data.table, zoo, stringr, geosphere, lubridate
# =============================================================================

# ---- Libraries ----

library(data.table)
library(zoo)
library(stringr)
library(geosphere)
library(lubridate)

# ---- Paths ----
# data_dir : folder containing daily_covariates_*.RDS files from Script 2b
# out_dir  : destination for weekly covariate output files

data_dir <- "data/"   # <-- update to your local data directory
out_dir  <- "data/"   # <-- update to your desired output directory

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Parameters ----

# Spatial domain (must match the domain used in Scripts 1 and 2a)
lat_min <- 43; lat_max <- 49
lon_min <- 5;  lon_max <- 18

# Temperature threshold: weekly median temperatures below this value (°C)
# are set to 0 to exclude biologically inactive periods for Ae. albopictus.
temp_threshold <- 15

# Lag window sizes (weeks) for rolling-mean predictors
lag_windows <- 2:10

# ---- Discover daily covariate files ----

covs_files <- list.files(data_dir,
                          pattern = "^daily_covariates_.*\\.RDS$",
                          full.names = TRUE)

if (length(covs_files) == 0) {
  stop("No daily covariate files found in: ", data_dir,
       "\nRun Script 2b first.")
}

# ---- Helper: extract year offset from filename ----
# The adjusted rasters in Script 2a retain E-OBS dates (2005-2024).
# The year offset maps these to the correct future period:
#   period 2036-2055 → offset = 2036 - 2005 = 31
#   period 2066-2085 → offset = 2066 - 2005 = 61
# If the E-OBS baseline in Script 2a differs, adjust eobs_start_year below.

eobs_start_year <- 2005  # first year of the E-OBS subset used in Script 2a

get_year_offset <- function(filename, eobs_start = eobs_start_year) {
  period      <- str_extract(basename(filename), "[0-9]{4}-[0-9]{4}")
  future_start <- as.integer(str_split(period, "-")[[1]][1])
  future_start - eobs_start
}

# ---- Pre-compute photoperiod across all files ----
# Photoperiod depends only on latitude and date, so it is computed once
# across all unique (ID, date) combinations and then joined to each file.

message("Pre-computing photoperiod across all scenario files ...")

all_covs <- rbindlist(lapply(covs_files, function(f) {
  dt <- as.data.table(readRDS(f))
  dt <- dt[lat >= lat_min & lat <= lat_max & lon >= lon_min & lon <= lon_max]
  year_offset <- get_year_offset(f)
  dt[, `:=`(
    date       = as.Date(date) + lubridate::years(year_offset),
    year       = lubridate::isoyear(as.Date(date) + lubridate::years(year_offset)),
    week       = lubridate::isoweek(as.Date(date) + lubridate::years(year_offset)),
    ID         = paste0("ID_", sprintf("%.4f", lat), "_", sprintf("%.4f", lon))
  )]
  dt[, temporalID := sprintf("%d_%02d", year, week)]
  dt[, .(ID, lat, date, year, week, temporalID)]
}))

lat_ids      <- unique(all_covs[, .(ID, lat)])
temporal_info <- unique(all_covs[, .(date, year, week, temporalID)])

phOut <- rbindlist(lapply(seq_len(nrow(lat_ids)), function(i) {
  if (i %% 100 == 0) message("  Photoperiod: ", i, " / ", nrow(lat_ids))
  tmp <- copy(temporal_info)
  tmp[, photop := geosphere::daylength(lat = lat_ids$lat[i], doy = date)]
  tmp[, ID := lat_ids$ID[i]]
  tmp <- tmp[, .(medianPhotoweek = median(photop, na.rm = TRUE)),
             by = .(year, week, ID)]
  tmp[, temporalID := sprintf("%d_%02d", year, week)]
  tmp[, join.id := paste0(ID, "_", week)]
  tmp[, .(join.id, medianPhotoweek)]
}))

message("Photoperiod computation complete.")

# ---- Process each scenario file ----

for (f in covs_files) {
  message("Processing: ", basename(f))

  # Parse scenario metadata from filename
  fname   <- basename(f)
  period  <- str_extract(fname, "[0-9]{4}-[0-9]{4}")
  model   <- str_remove(fname, "^daily_covariates_[0-9]{4}-[0-9]{4}_")
  model   <- tools::file_path_sans_ext(model)
  year_offset <- get_year_offset(f)

  # Load and filter to study domain
  covariates <- as.data.table(readRDS(f))
  covariates <- covariates[lat >= lat_min & lat <= lat_max &
                             lon >= lon_min & lon <= lon_max]

  # Shift dates from E-OBS baseline to target future period
  covariates[, `:=`(
    date = as.Date(date) + lubridate::years(year_offset),
    year = lubridate::isoyear(as.Date(date) + lubridate::years(year_offset)),
    week = lubridate::isoweek(as.Date(date) + lubridate::years(year_offset))
  )]

  # Weekly aggregation: median temperature, sum precipitation
  covariates <- covariates[, .(
    date     = date[1],
    med_TAVG = median(temperature,   na.rm = TRUE),
    sum_PREC = sum(precipitation,    na.rm = TRUE)
  ), by = .(lon, lat, week, year)]

  # Create stable grid cell ID and join keys
  covariates[, `:=`(
    ID         = paste0("ID_", sprintf("%.4f", lat), "_", sprintf("%.4f", lon)),
    temporalID = sprintf("%d_%02d", year, week),
    join.id    = paste0("ID_", sprintf("%.4f", lat), "_", sprintf("%.4f", lon),
                         "_", week)
  )]

  # Join precomputed weekly photoperiod
  covariates <- merge(covariates, phOut, by = "join.id", all.x = TRUE)

  # ---- Compute rolling-mean lagged predictors ----
  setorder(covariates, ID, year, week)

  for (col in c("med_TAVG", "sum_PREC", "medianPhotoweek")) {
    for (k in lag_windows) {
      new_col <- paste0(col, ".lag", k)
      covariates[, (new_col) := zoo::rollmean(get(col), k,
                                               align = "right",
                                               fill = NA),
                 by = ID]
    }
  }

  # ---- Apply minimum temperature threshold ----
  # Set lagged temperature covariates to 0 when below the biological threshold.
  temp_lag_cols <- grep("^med_TAVG\\.lag", names(covariates), value = TRUE)
  for (col in temp_lag_cols) {
    covariates[, (col) := fifelse(get(col) < temp_threshold, 0, get(col))]
  }

  # ---- Back-fill leading NAs from rolling windows ----
  # The first (k-1) rows per cell will be NA due to the rolling window.
  # These are filled with the next valid value (backward fill) to avoid
  # dropping the earliest weeks of each time series.
  lagged_cols <- grep("\\.lag[0-9]", names(covariates), value = TRUE)
  for (col in lagged_cols) {
    covariates[, (col) := zoo::na.locf(get(col),
                                        na.rm    = FALSE,
                                        fromLast = TRUE),
               by = ID]
  }

  # ---- Save output ----
  out_name <- file.path(out_dir,
                         sprintf("covariates_%s_%s_with-lags_%s.RDS",
                                 model, period, Sys.Date()))
  saveRDS(covariates, out_name)

  rm(covariates)
  gc()

  message("Saved: ", basename(out_name))
}
