# =============================================================================
# Script 1: Model Training and Current-Period Spatial Predictions
# =============================================================================
# Description:
#   Trains a stacked ensemble model (XGBoost + Random Forest + Cubist → LM)
#   using mlr3 on observed mosquito egg-count data with lagged climate
#   covariates (temperature, precipitation, photoperiod). The script performs
#   spatial and temporal train/test splitting, hyperparameter tuning, internal
#   cross-validation, and spatial predictions over the current period
#   (2005–2024) using gridded E-OBS covariates.
#
# Inputs:
#   - biomatrix_<date>.RDS       : Observation data with lagged climate covariates
#   - eobs_covariates_<date>.RDS : Gridded spatial covariates for current period
#
# Outputs:
#   - test-ids_<date>.RDS                    : IDs used for spatial test set
#   - proportions_<date>.RDS                 : Train/test set proportions
#   - TunedHyperParam_Model_<date>.RDS       : Tuning results
#   - stacked_Model_train_<date>.RDS         : Trained stacked learner
#   - stacked_Model_train_internalCV_<date>.RDS : Cross-validation results
#   - mlr3_stack_Pred_<date>.RDS             : Predictions (train + test sets)
#   - mlr3_stack_Obs_<date>.RDS              : Observations (train + test sets)
#   - predictions_current_scenario_<date>.RDS: Spatial predictions 2005–2024
#
# Dependencies: See library calls below.
# R version: >= 4.2.0
# Note: mlr3extralearners must be installed from GitHub:
#   remotes::install_github("mlr-org/mlr3extralearners@*release")
# =============================================================================

# ---- Libraries ----

library(mlr3)
library(mlr3pipelines)
library(mlr3tuning)
library(mlr3viz)
library(mlr3extralearners)
library(mlr3learners)
library(mlr3fselect)
library(paradox)
library(data.table)
library(dplyr)
library(lubridate)
library(stringr)
library(skimr)
library(Cubist)
library(ggplot2)
library(future)

# ---- Paths ----
# Set data_dir to the root folder containing all input/output data files.
# Set fig_dir to the folder where figures will be saved.

data_dir <- "data/"          # <-- update to your local data directory
fig_dir  <- "figures/"       # <-- update to your local figures directory

dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir,  showWarnings = FALSE, recursive = TRUE)

# ---- 1. Load observations ----

bio.matrix <- readRDS(
  tail(sort(list.files(data_dir, pattern = "^biomatrix_.*\\.RDS$",
                       full.names = TRUE)), 1)
)

bio.matrix <- bio.matrix %>%
  select(ID, Country, Region, value, lon, lat, temporalID, date, year, week,
         med_TAVG.lag2:med_TAVG.lag10,
         sum_PREC.lag2:sum_PREC.lag10,
         medianPhotoweek.lag2:medianPhotoweek.lag10)

# Exploratory plot: mean egg counts by week and region
bio.matrix %>%
  group_by(week, Region) %>%
  summarize(mean_value = mean(value, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = as.integer(week), y = mean_value)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~Region) +
  labs(x = "Week", y = "Mean egg count") +
  theme_light()

# ---- 2. Define predictors and prepare modelling dataset ----

lag_vars <- c(
  paste0("med_TAVG.lag",        2:10),
  paste0("sum_PREC.lag",        2:10),
  paste0("medianPhotoweek.lag", 2:10)
)

vars <- unique(make.names(unlist(
  sapply(lag_vars, function(i) names(bio.matrix)[grep(i, names(bio.matrix))])
)))

formula <- as.formula(paste("value ~", paste(vars, collapse = "+")))
r.sel   <- stats::complete.cases(bio.matrix[, all.vars(formula)])

df.s <- bio.matrix[which(r.sel),
                   c(all.vars(formula), "ID", "year", "date", "Country", "Region")]

# ---- 3. Spatial and temporal train/test split ----

# Sample one site per region (with >= 3 sites) for spatial hold-out test set
test.ids <- df.s %>%
  group_by(Region) %>%
  filter(n_distinct(ID) >= 3) %>%
  slice_sample(n = 1) %>%
  ungroup() %>%
  pull(ID)

saveRDS(test.ids, file.path(data_dir, paste0("test-ids_", Sys.Date(), ".RDS")))

# Years used as temporal hold-out test set
temporal_test_years <- c(2012, 2014, 2022)

# Spatial test set (held-out sites, excluding temporal test years)
test.df       <- df.s %>% filter(ID %in% test.ids,
                                  !(year %in% temporal_test_years))
date.test     <- test.df[, c("ID", "Country", "Region", "year", "date")]
test.df       <- test.df %>% select(-ID, -year, -date, -Country, -Region)

# Temporal test set (held-out years)
test.yrs      <- df.s %>% filter(year %in% temporal_test_years)
date.test.yrs <- test.yrs[, c("ID", "Country", "Region", "year", "date")]
test.yrs      <- test.yrs %>% select(-ID, -year, -date, -Country, -Region)

# Training set (excluding held-out sites and years)
train.df      <- df.s %>% filter(!(ID %in% test.ids),
                                   !(year %in% temporal_test_years))
date.train    <- train.df[, c("ID", "Country", "Region", "year", "date")]
id            <- train.df$ID
train.df      <- train.df %>% select(-ID, -year, -date, -Country, -Region)

# Integrity check
if (nrow(test.df) + nrow(train.df) + nrow(test.yrs) != nrow(df.s)) {
  stop("Row count mismatch between subsets and full dataset — check split logic.")
}

proportions <- data.frame(
  Dataset    = c("Training set", "Spatial test", "Temporal test"),
  Proportion = c(round(nrow(train.df) / nrow(df.s), 3),
                 round(nrow(test.df)  / nrow(df.s), 3),
                 round(nrow(test.yrs) / nrow(df.s), 3))
)
print(proportions)
saveRDS(proportions, file.path(data_dir, paste0("proportions_", Sys.Date(), ".RDS")))

# ---- 4. Define the stacked ensemble learner ----

tv <- all.vars(formula)[1]  # response variable name
fts <- names(train.df[, 2:ncol(train.df)])  # predictor names

# Base learners run in cross-validation, predictions fed to meta-learner (LM)
stacked_graph <- gunion(list(
  po("learner_cv", lrn("regr.xgboost", predict_type = "response", nthread = 5)),
  po("learner_cv", lrn("regr.ranger",  predict_type = "response",
                        num.threads = 5, importance = "permutation")),
  po("learner_cv", lrn("regr.cubist",  predict_type = "response"))
)) %>>%
  po("featureunion") %>>%
  lrn("regr.lm", id = "lm_ensemble", predict_type = "se")

# stacked_graph$plot(html = FALSE)  # uncomment to visualise pipeline graph

stacked_graph$keep_results <- TRUE
stacked_learner <- as_learner(stacked_graph)

# ---- 5. Hyperparameter search space ----

stacked.search_space <- paradox::ps(
  regr.ranger.num.trees        = p_int(lower = 200, upper = 500),
  regr.ranger.mtry             = p_int(lower = round((length(fts) - 1) * 0.5),
                                        upper = (length(fts) - 3)),
  regr.ranger.min.node.size    = p_int(lower = 5, upper = 10),
  regr.xgboost.nrounds         = p_int(10, 30),
  regr.xgboost.max_depth       = p_int(4, 10),
  regr.xgboost.eta             = p_dbl(0.2, 0.4),
  regr.xgboost.subsample       = p_dbl(0.9, 1),
  regr.xgboost.min_child_weight = p_int(1, 4),
  regr.xgboost.colsample_bytree = p_dbl(0.5, 0.6),
  regr.cubist.committees       = p_int(lower = 1, upper = 10)
)

# ---- 6. Tuning ----

task <- mlr3::as_task_regr(train.df, target = tv, id = "eggs_regr")
task$cbind(data.table("ID" = id))
task$set_col_roles("ID", roles = "group")

tune.eml <- TuningInstanceBatchSingleCrit$new(
  task        = task,
  learner     = stacked_learner,
  resampling  = rsmp("cv", folds = 10),
  measure     = msr("regr.rmse"),
  search_space = stacked.search_space,
  terminator  = trm("evals", n_evals = 5)
)

set.seed(123)
tuner <- tnr("random_search")
future::plan("multisession", workers = 5)
tuner$optimize(tune.eml)

saveRDS(tune.eml,
        file.path(data_dir, paste0("TunedHyperParam_Model_", Sys.Date(), ".RDS")))

# ---- 7. Train final model with best hyperparameters ----

best.learner <- stacked_learner
best.learner$param_set$values <- tune.eml$result$learner_param_vals[[1]]

# Re-create task (ensures clean state after tuning)
task <- as_task_regr(train.df, target = "value", id = "eggs_regr")
task$cbind(data.table("ID" = id))
task$set_col_roles("ID", roles = "group")

best.learner$train(task)

summary(best.learner$model$lm_ensemble$model)
AIC(best.learner$model$lm_ensemble$model)

saveRDS(best.learner,
        file.path(data_dir, paste0("stacked_Model_train_", Sys.Date(), ".RDS")))

# ---- 8. Internal cross-validation ----

myresmpl     <- rsmp("cv", folds = 10L)
cv_predictions <- mlr3::resample(task     = task,
                                  learner  = best.learner,
                                  resampling = myresmpl)

cv_predictions$aggregate(msr("regr.rmse"))

saveRDS(cv_predictions,
        file.path(data_dir, paste0("stacked_Model_train_internalCV_", Sys.Date(), ".RDS")))

# ---- 9. Predict on training and test sets ----

pred.train    <- best.learner$predict(task)
pred.test     <- predict(best.learner, newdata = test.df)
pred.test.yrs <- predict(best.learner, newdata = test.yrs)

fitted.stack <- data.frame(
  ID        = date.train$ID,
  Region    = date.train$Region,
  Country   = date.train$Country,
  year      = date.train$year,
  date      = date.train$date,
  value     = round(pred.train$response, 3),
  se        = round(pred.train$se, 3),
  model     = "stacked",
  partition = "Train"
)

test.stack <- data.frame(
  ID        = date.test$ID,
  Region    = date.test$Region,
  Country   = date.test$Country,
  year      = date.test$year,
  date      = date.test$date,
  value     = round(pred.test$response, 3),
  se        = round(pred.test$se, 3),
  model     = "stacked",
  partition = "Test_spatial"
)

test.stack.yrs <- data.frame(
  ID        = date.test.yrs$ID,
  Region    = date.test.yrs$Region,
  Country   = date.test.yrs$Country,
  year      = date.test.yrs$year,
  date      = date.test.yrs$date,
  value     = round(pred.test.yrs$response, 3),
  se        = round(pred.test.yrs$se, 3),
  model     = "stacked",
  partition = "Test_temporal"
)

myPred <- rbind(fitted.stack, test.stack, test.stack.yrs)
myPred$week       <- lubridate::isoweek(myPred$date)
myPred$year       <- lubridate::isoyear(myPred$date)
myPred$temporalID <- paste0(myPred$year, "_", str_pad(myPred$week, 2, pad = "0"))

saveRDS(myPred,
        file.path(data_dir, paste0("mlr3_stack_Pred_", Sys.Date(), ".RDS")))

# Export observations for comparison
myObs <- bind_rows(
  cbind(date.train,    value = train.df$value, partition = "Train",         model = "Observations"),
  cbind(date.test,     value = test.df$value,  partition = "Test_spatial",  model = "Observations"),
  cbind(date.test.yrs, value = test.yrs$value, partition = "Test_temporal", model = "Observations")
)
myObs$week       <- lubridate::isoweek(myObs$date)
myObs$year       <- lubridate::isoyear(myObs$date)
myObs$temporalID <- paste0(myObs$year, "_", str_pad(myObs$week, 2, pad = "0"))

saveRDS(myObs,
        file.path(data_dir, paste0("mlr3_stack_Obs_", Sys.Date(), ".RDS")))

# ---- 10. Spatial predictions over the current period (2005–2024) ----

spatial_cov_full <- readRDS(
  tail(sort(list.files(data_dir,
                       pattern = "^eobs_covariates_.*\\.RDS$",
                       full.names = TRUE)), 1)
)

my_spPred_all <- list()

for (yr in sort(unique(spatial_cov_full$year))) {
  message("Processing year: ", yr)

  spatial_cov <- spatial_cov_full %>% filter(year == yr)

  test.sp <- spatial_cov %>%
    ungroup() %>%
    select(-lon, -lat, -ID, -week, -year, -date, -temporalID, -join.id)

  gc()

  prediction <- predict(best.learner, newdata = test.sp)

  my_spPred_all[[as.character(yr)]] <- spatial_cov %>%
    select(lon, lat, ID, week, year, date, temporalID, join.id) %>%
    bind_cols(
      spatial_pred    = round(prediction$response),
      spatial_pred_se = round(prediction$se)
    )
}

my_spPred_all_df <- bind_rows(my_spPred_all)

saveRDS(my_spPred_all_df,
        file.path(data_dir, paste0("predictions_current_scenario_", Sys.Date(), ".RDS")))
