# =============================================================================
# Script 2b: Spatial Predictions on Future Climate Scenarios
# =============================================================================
# Description:
#   Applies the trained stacked ensemble model (Script 1) to future-period
#   spatial covariate datasets (produced from change-factor adjusted climate
#   grids) to generate weekly mosquito egg-count predictions for each
#   future scenario and model combination.
#
#   Predictions are made year-by-year to manage memory, then combined and
#   saved per scenario file.
#
# Inputs:
#   - stacked_Model_train_<date>.RDS  : Trained stacked learner (Script 1)
#   - covariates_<period>_<date>_<model>_with-lags.RDS : Future spatial
#     covariates (one file per model/period combination)
#
# Outputs (written to out_dir):
#   - Predictions_<model>_<period>_<date>.RDS : Weekly spatial predictions
#     for each scenario (response + SE)
#
# Dependencies: mlr3, dplyr, stringr
# =============================================================================

# ---- Libraries ----

library(mlr3)
library(dplyr)
library(stringr)

# ---- Paths ----
# data_dir   : folder containing future covariate RDS files
# model_dir  : folder containing the trained model RDS file
# out_dir    : destination for prediction output files

data_dir  <- "data/"   # <-- update to your local data directory
model_dir <- "data/"   # <-- update if model file is in a different location
out_dir   <- "data/"   # <-- update to your desired output directory

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Load trained model ----

best.learner <- readRDS(
  tail(sort(list.files(model_dir, pattern = "^stacked_Model_train_.*\\.RDS$",
                       full.names = TRUE)), 1)
)

# ---- List future covariate files ----

future_covs <- list.files(
  path        = data_dir,
  pattern     = "covariates_.*_with-lags\\.RDS$",
  full.names  = TRUE
)

if (length(future_covs) == 0) {
  stop("No future covariate files found in: ", data_dir,
       "\nExpected pattern: covariates_*_with-lags.RDS")
}

# ---- Generate predictions for each scenario ----

for (file in future_covs) {

  message("Loading: ", basename(file))

  spatial_cov <- readRDS(file) %>%
    select(tileID, join.id, lon, lat, temporalID, date, year, week,
           med_TAVG.lag2:med_TAVG.lag10,
           sum_PREC.lag2:sum_PREC.lag10,
           medianPhotoweek.lag2:medianPhotoweek.lag10) %>%
    na.omit()

  # Predict year by year to limit memory usage
  yearly_preds <- spatial_cov %>%
    split(.$year) %>%
    lapply(function(df_year) {
      message("  Processing year: ", unique(df_year$year))

      test.sp <- df_year %>%
        ungroup() %>%
        select(-lon, -lat, -tileID, -week, -temporalID, -date, -year, -join.id)

      spPred <- predict(best.learner, newdata = test.sp)

      df_year %>%
        select(lon, lat, tileID, week, year, temporalID) %>%
        bind_cols(
          spatial_pred    = round(spPred$response),
          spatial_pred_se = round(spPred$se)
        )
    })

  all_preds <- bind_rows(yearly_preds)

  # Extract model name and period from filename for output naming
  file_name <- basename(file)
  years     <- str_extract(file_name, "\\d{4}-\\d{4}")
  model     <- str_extract(file_name, "(?<=_20[0-9]{2}-[0-9]{4}_)[A-Za-z0-9\\-]+(?=_)")

  outname <- file.path(out_dir,
                        paste0("Predictions_", model, "_", years, "_", Sys.Date(), ".RDS"))
  saveRDS(all_preds, outname)

  rm(all_preds, spatial_cov, yearly_preds)
  gc()

  message("Saved: ", outname)
}
