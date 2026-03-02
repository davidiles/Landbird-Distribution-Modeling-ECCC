# ============================================================
# 11b_crossvalidation_brt.R
#
# Purpose:
#   Spatial block cross-validation for Poisson boosted regression trees (BRT)
#   to compare with INLA CV outputs from 11a.
#
# Outputs (mirrors 11a):
#   data_clean/model_output/cv_brt/predictions/<sp>.rds
#   data_clean/model_output/cv_brt/summaries/<sp>.rds
#   data_clean/model_output/cv_brt/block_summaries/<sp>.rds
#   data_clean/model_output/cv_brt/meta/cv_plan.rds (reused from 11a if you point to same)
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(sf)
  
  library(tidymodels)
  library(xgboost)   # engine used by parsnip
})

# ------------------------------------------------------------
# Config (match 11a wherever possible)
# ------------------------------------------------------------

in_file <- "data_clean/birds/data_ready_for_analysis.rds"

out_dir <- "data_clean/model_output/xval_BRT"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "predictions"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "summaries"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "meta"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "block_summaries"), recursive = TRUE, showWarnings = FALSE)

# Source helpers (reuse your existing ones)
source("R/functions/cv_utils.R")       # score_cv_predictions(), read_cv_plan_or_create(), etc
source("R/functions/spatial_utils.R")  # make_spatial_blocks_grid(), make_block_polygons_grid()
source("R/functions/inla_model_utils.R")

# Species selection (same as 11a)
min_detections_obba3 <- 100
min_squares_obba3 <- 50

# CV design (same as 11a)
block_size_km <- 30
n_folds <- 10
n_repeats <- 1
cv_seed <- 123

# Covariates: same list used in 11a for INLA fixed effects
base_covars <- c(
  "prec",
  "tmax",
  "on_river",
  "on_road",
  "urban_3",
  "lc_1",
  "lc_4",
  "lc_5",
  "lc_8",
  "lc_9",
  "lc_10",
  "lc_11",
  "lc_12",
  "lc_14",
  "lc_17",
  "insect_broadleaf",
  "insect_needleleaf"
)

# BRT hyperparameters (start conservative; adjust later)
brt_spec <- boost_tree(
  trees = 2000,
  tree_depth = 6,
  min_n = 10,
  loss_reduction = 0,
  sample_size = 0.8,
  learn_rate = 0.05,
  mtry = NULL,
  stop_iter = Inf
) %>%
  set_engine("xgboost", objective = "count:poisson") %>%
  set_mode("regression")

# ------------------------------------------------------------
# Tuning spec + function (run once per species per repeat)
# ------------------------------------------------------------

brt_tune_spec <- boost_tree(
  trees = tune(),
  tree_depth = tune(),
  min_n = tune(),
  loss_reduction = tune(),
  sample_size = tune(),
  mtry = tune(),
  learn_rate = tune(),
  stop_iter = Inf
) %>%
  set_engine("xgboost", objective = "count:poisson") %>%
  set_mode("regression")

tune_brt_once <- function(dat_train, rec, seed = 123, v_inner = 5, grid_size = 25) {
  # dat_train: training data for outer fold 1 (must include block_id + count)
  # rec: recipe built from dat_train (your make_brt_recipe output)
  
  stopifnot("count" %in% names(dat_train))
  stopifnot("block_id" %in% names(dat_train))
  
  # Prep recipe on training data to learn dummy columns etc.
  rec_prep <- prep(rec, training = dat_train, verbose = FALSE)
  
  # Get number of post-processed predictors (for mtry range)
  x_train <- juice(rec_prep) %>% select(-count)
  p <- ncol(x_train)
  if (p < 1) stop("After preprocessing there are no predictors (p < 1).")
  
  # Inner CV (spatial grouping)
  set.seed(seed)
  inner_folds <- group_vfold_cv(dat_train, group = block_id, v = v_inner)
  
  wf_tune <- workflow() %>%
    add_recipe(rec) %>%
    add_model(brt_tune_spec)
  
  # Define parameter ranges (sane-ish defaults; adjust later if needed)
  param_set <- wf_tune %>%
    extract_parameter_set_dials() %>%
    update(
      trees = trees(range = c(300L, 3000L)),
      tree_depth = tree_depth(range = c(2L, 8L)),
      min_n = min_n(range = c(5L, 50L)),
      learn_rate = learn_rate(range = c(-3, -1)),  
      sample_size = sample_prop(range = c(0.6, 1.0)),
      loss_reduction = loss_reduction(range = c(-5, 1)),
      mtry = mtry(range = c(1L, min(p, 250L)))
    )
  
  # Latin hypercube grid in that parameter space
  set.seed(seed + 1)
  grid <- grid_latin_hypercube(param_set, size = grid_size)
  
  # Tune
  set.seed(seed + 2)
  res <- tune_grid(
    wf_tune,
    resamples = inner_folds,
    grid = grid,
    metrics = yardstick::metric_set(
      yardstick::poisson_log_loss,
      yardstick::rmse,
      yardstick::mae
    ),
    control = control_grid(verbose = FALSE, parallel_over = "resamples")
  )
  
  best <- select_best(res, metric = "poisson_log_loss")
  
  list(
    best = best,
    p = p,
    tuning_metrics = collect_metrics(res),
    res = res
  )
}

# ------------------------------------------------------------
# Load data
# ------------------------------------------------------------

stopifnot(file.exists(in_file))
dat <- readRDS(in_file)

all_surveys <- dat$all_surveys %>%
  mutate(BCR_numeric = as.numeric(as.factor(BCR)))

counts      <- dat$counts
study_boundary <- dat$study_boundary %>% st_as_sf()
species_to_model <- dat$species_to_model

stopifnot(nrow(all_surveys) == nrow(counts))
stopifnot(inherits(all_surveys, "sf"))

species_run <- species_to_model %>%
  filter(total_detections_OBBA3 >= min_detections_obba3,
         total_squares_OBBA3 >= min_squares_obba3) %>%
  na.omit()

message("Species queued for BRT CV: ", nrow(species_run))

# ------------------------------------------------------------
# Create/reuse same CV plan as 11a (block IDs must match)
# ------------------------------------------------------------

cv_plan_path <- file.path("data_clean/model_output/xval_INLA", "meta", paste0("cv_plan_blocks_", block_size_km, "km.rds"))

plan_obj <- read_cv_plan_or_create(
  cv_plan_path = cv_plan_path,
  all_surveys = all_surveys,
  block_size_km = block_size_km,
  n_folds = n_folds,
  n_repeats = n_repeats,
  cv_seed = cv_seed,
  out_dir = out_dir
)

block_ids   <- plan_obj$block_ids
block_polys <- plan_obj$block_polys
cv_plan     <- plan_obj$cv_plan

message("Using CV plan: ", cv_plan_path)

# ------------------------------------------------------------
# Helper: build a recipe that mirrors INLA’s “available predictors”
# ------------------------------------------------------------

make_brt_recipe <- function(df_train, base_covars) {
  stopifnot("count" %in% names(df_train))
  
  desired_predictors <- c(
    base_covars,
    "days_since_june15",
    "Hours_Since_Sunrise",
    "BCR_numeric",
    "ARU",
    "Atlas3"
  )
  
  predictors_present <- intersect(desired_predictors, names(df_train))
  
  missing_preds <- setdiff(desired_predictors, predictors_present)
  if (length(missing_preds) > 0) {
    message(
      "BRT recipe: predictors missing from df_train and will be ignored: ",
      paste(missing_preds, collapse = ", ")
    )
  }
  
  recipe(count ~ ., data = df_train) %>%
    # mark ONLY these columns as predictors; everything else gets "ignore"
    update_role(all_of(predictors_present), new_role = "predictor") %>%
    update_role(-all_of(c("count", predictors_present)), new_role = "ignore") %>%
    step_rm(has_role("ignore")) %>%
    
    # preprocessing
    step_unknown(all_nominal_predictors()) %>%
    step_novel(all_nominal_predictors()) %>%
    step_impute_median(all_numeric_predictors()) %>%
    step_dummy(all_nominal_predictors(), one_hot = TRUE) %>%
    step_zv(all_predictors())
}

# ------------------------------------------------------------
# Main CV loop
# ------------------------------------------------------------

species_to_check <- c("Bobolink",
                      "Blue Jay",
                      "Canada Jay",
                      "Olive-sided Flycatcher",
                      "Winter Wren",
                      "Lesser Yellowlegs",
                      "Blackpoll Warbler",
                      "Connecticut Warbler",
                      "Palm Warbler",
                      "Lincoln's Sparrow",
                      "Fox Sparrow",
                      "Common Nighthawk",
                      "Long-eared Owl",
                      "American Tree Sparrow",
                      "LeConte's Sparrow",
                      "Nelson's Sparrow",
                      "Boreal Chickadee",
                      "Rusty Blackbird",
                      "Yellow-bellied Flycatcher",
                      "Greater Yellowlegs",
                      "Hudsonian Godwit",
                      "Canada Warbler",
                      "Eastern Wood-Peewee",
                      "Grasshopper Sparrow",
                      "Solitary Sandpiper",
                      "White-throated Sparrow",
                      "Bay-breasted Warbler")


species_run <- species_run %>%
  subset(english_name %in% species_to_check)

for (i in seq_len(nrow(species_run))) {
  
  sp_english <- species_run$english_name[i]
  sp_code <- as.character(species_run$species_id[i])
  sp_file <- sp_filename(sp_english)
  
  message("\n====================\n", i, "/", nrow(species_run), ": ", sp_english, " (", sp_code, ")\n====================")
  
  if (!(sp_code %in% names(counts))) {
    message("Skipping (species_id not found in counts columns): ", sp_code)
    next
  }
  
  pred_path  <- file.path(out_dir, "predictions", paste0(sp_file, ".rds"))
  summ_path  <- file.path(out_dir, "summaries",   paste0(sp_file, ".rds"))
  block_path <- file.path(out_dir, "block_summaries", paste0(sp_file, ".rds"))
  
  tune_dir <- file.path(out_dir, "meta", "tuned_params")
  dir.create(tune_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Build species modeling data (same as 11a)
  sp_dat <- all_surveys %>%
    mutate(
      count = counts[[sp_code]],
      days_since_june15 = DayOfYear - 166,
      BCR_factor = as.numeric(factor(BCR)),
      Atlas3 = ifelse(Atlas == "OBBA2", 0, 1),
      block_id = block_ids
    )
  
  covars_present <- intersect(base_covars, names(sp_dat))
  
  for (r in seq_len(n_repeats)) {
    
    plan_r <- cv_plan %>% filter(rep == r)
    tune_path <- file.path(tune_dir, paste0(sp_file, "_rep", r, "_tuned_params.rds"))
    
    for (f in seq_len(n_folds)) {
      
      set.seed(999 + 100*r + f)
      
      test_blocks <- plan_r %>% filter(fold == f) %>% pull(block_id)
      
      is_test <- sp_dat$block_id %in% test_blocks
      dat_train_sf <- sp_dat[!is_test, ]
      dat_test_sf  <- sp_dat[ is_test, ]
      
      expected_pairs <- sf::st_drop_geometry(dat_test_sf) %>%
        dplyr::distinct(block_id, Atlas)
      
      if (fold_already_completed(block_path, sp_code, block_size_km, r, f, expected_pairs)) {
        message("BRT CV: repeat ", r, "/", n_repeats, " fold ", f, "/", n_folds,
                " already completed for ", sp_english, " — skipping.")
        next
      }
      
      if (nrow(dat_test_sf) == 0 || nrow(dat_train_sf) == 0) {
        message("Repeat ", r, " fold ", f, ": empty train/test split; skipping.")
        next
      }
      
      message("BRT CV: repeat ", r, "/", n_repeats, " fold ", f, "/", n_folds,
              " | train n=", nrow(dat_train_sf), " test n=", nrow(dat_test_sf))
      
      # Drop geometry for modeling
      dat_train <- sf::st_drop_geometry(dat_train_sf)
      dat_test  <- sf::st_drop_geometry(dat_test_sf)
      
      # Build workflow
      rec <- make_brt_recipe(dat_train, covars_present)
      
      # --- Tune once per species per repeat using outer fold 1 training data
      tune_dir <- file.path(out_dir, "meta", "tuned_params")
      dir.create(tune_dir, recursive = TRUE, showWarnings = FALSE)
      tune_path <- file.path(tune_dir, paste0(sp_file, "_rep", r, "_tuned_params.rds"))
      
      if (!file.exists(tune_path)) {
        if (f == 1) {
          message("Tuning BRT hyperparameters (once) for ", sp_english, " rep=", r, " using outer fold 1 training data...")
          
          # Important: dat_train must include block_id for group_vfold_cv
          # and count must be present for training
          tune_out <- tune_brt_once(
            dat_train = dat_train,
            rec = rec,
            seed = 1000 + r,
            v_inner = 5,
            grid_size = 15
          )
          
          saveRDS(tune_out, tune_path)
          message("Saved tuned params: ", tune_path)
        } else {
          # If fold 1 was skipped (already completed) but tuning wasn't created,
          # we force tuning now using the current fold's training data.
          message("Tuned params missing but fold != 1; tuning now using current training data for ", sp_english, " rep=", r)
          tune_out <- tune_brt_once(
            dat_train = dat_train,
            rec = rec,
            seed = 1000 + r,
            v_inner = 5,
            grid_size = 25
          )
          saveRDS(tune_out, tune_path)
        }
      }
      
      tune_out <- readRDS(tune_path)
      best_params <- tune_out$best
      
      # --- Finalize workflow with tuned params
      wf_final <- workflow() %>%
        add_recipe(rec) %>%
        add_model(brt_tune_spec) %>%   # model with tune() params
        finalize_workflow(best_params)
      
      # Fit final workflow on this fold's training data
      fit_wf <- try(workflows::fit(wf_final, dat_train), silent = TRUE)
      
      if (inherits(fit_wf, "try-error")) {
        message("BRT fit error (rep=", r, ", fold=", f, "):")
        message(attr(fit_wf, "condition")$message)
        fit_ok <- FALSE
      } else {
        fit_ok <- TRUE
      }
      
      if (!fit_ok) {
        message("BRT fit failed (rep=", r, ", fold=", f, "). Recording NA summary + continuing.")
        
        summ_row <- tibble(
          sp_english = sp_english,
          sp_code = sp_code,
          block_size_km = block_size_km,
          rep = r,
          fold = f,
          Atlas = NA_character_,
          n = nrow(dat_test),
          rmse = NA_real_,
          mae  = NA_real_,
          spearman = NA_real_,
          pearson  = NA_real_,
          poisson_dev = NA_real_,
          mean_obs = NA_real_,
          mean_pred_mean = NA_real_,
          mean_pred_q025 = NA_real_,
          mean_pred_q50  = NA_real_,
          mean_pred_q975 = NA_real_,
          mean_pred_yrep_mean = NA_real_,
          mean_pred_yrep_q025 = NA_real_,
          mean_pred_yrep_q50  = NA_real_,
          mean_pred_yrep_q975 = NA_real_,
          coverage_95 = NA_real_,
          mean_log_score = NA_real_,
          brier_ge1 = NA_real_
        )
        
        append_rds_rows(
          summ_path,
          summ_row,
          key_cols = c("sp_code","rep","fold","block_size_km","Atlas")
        )
        next
      }
      
      # # Optional: extract fitted model
      # fit_parsnip <- workflows::extract_fit_parsnip(fit_wf)
      # 
      # # extract the raw xgboost engine (xgb.Booster)
      # xgb_fit <- workflows::extract_fit_engine(fit_wf)
      # xgb_fit
      # 
      # # Variable importance
      # imp <- xgboost::xgb.importance(model = xgb_fit)
      # head(imp, 20)
      # 
      # # Partial dependence plots
      # pred_fun <- function(object, newdata) {
      #   predict(object, new_data = newdata) %>%
      #     dplyr::pull(.pred)
      # }
      # 
      # pdp::partial(
      #   object = fit_wf,
      #   pred.var = "tmax",
      #   train = dat_train,
      #   pred.fun = pred_fun,
      #   type = "regression",
      #   plot = TRUE
      # )
      
      # Predict mean (mu)
      mu_pred <- predict(fit_wf, new_data = dat_test) %>% pull(.pred)
      mu_pred <- pmax(as.numeric(mu_pred), 1e-12)
      
      # “Probabilistic” summaries using Poisson observation noise only
      yrep_q025 <- qpois(0.025, lambda = mu_pred)
      yrep_q50  <- qpois(0.50,  lambda = mu_pred)
      yrep_q975 <- qpois(0.975, lambda = mu_pred)
      p_ge1 <- 1 - exp(-mu_pred)
      
      pred_rows <- dat_test_sf %>%
        mutate(
          sp_english = sp_english,
          sp_code = sp_code,
          block_size_km = block_size_km,
          rep = r,
          fold = f,
          
          survey_id = if ("survey_id" %in% names(dat_test_sf)) dat_test_sf$survey_id else NA_character_,
          Atlas = Atlas,
          count_obs = count,
          mu_pred = mu_pred,
          yrep_q025 = yrep_q025,
          yrep_q50  = yrep_q50,
          yrep_q975 = yrep_q975,
          p_ge1 = p_ge1
        ) %>%
        select(
          sp_english, sp_code, block_size_km, rep, fold,
          block_id, survey_id, Atlas,
          count_obs, mu_pred,
          yrep_q025, yrep_q50, yrep_q975, p_ge1,
          geometry
        )
      
      append_rds_rows(
        pred_path,
        pred_rows,
        key_cols = c("sp_code","rep","fold","block_id","survey_id","Atlas")
      )
      
      # Fold-level scoring per Atlas (same structure as 11a)
      test_df <- sf::st_drop_geometry(pred_rows)
      
      atlases_in_test <- sort(unique(test_df$Atlas))
      for (atl in atlases_in_test) {
        test_df_atl <- test_df %>% filter(Atlas == atl)
        if (nrow(test_df_atl) == 0) next
        
        fold_score <- score_cv_predictions(test_df_atl, y_col = "count_obs", mu_col = "mu_pred")
        
        cov95 <- interval_coverage(
          y = test_df_atl$count_obs,
          lo = test_df_atl$yrep_q025,
          hi = test_df_atl$yrep_q975
        )
        
        mean_log_score <- mean(dpois(test_df_atl$count_obs, lambda = test_df_atl$mu_pred, log = TRUE), na.rm = TRUE)
        
        brier_ge1 <- brier_score(
          y01 = as.numeric(test_df_atl$count_obs >= 1),
          p = test_df_atl$p_ge1
        )
        
        fold_mean_obs <- mean(test_df_atl$count_obs, na.rm = TRUE)
        fold_mean_pred <- mean(test_df_atl$mu_pred, na.rm = TRUE)
        
        summ_row <- fold_score %>%
          mutate(
            sp_english = sp_english,
            sp_code = sp_code,
            block_size_km = block_size_km,
            rep = r,
            fold = f,
            Atlas = atl,
            mean_obs = fold_mean_obs,
            mean_pred_mean = fold_mean_pred,
            mean_pred_q025 = NA_real_,
            mean_pred_q50  = NA_real_,
            mean_pred_q975 = NA_real_,
            mean_pred_yrep_mean = fold_mean_pred,
            mean_pred_yrep_q025 = NA_real_,
            mean_pred_yrep_q50  = NA_real_,
            mean_pred_yrep_q975 = NA_real_,
            coverage_95 = cov95,
            mean_log_score = mean_log_score,
            brier_ge1 = brier_ge1
          ) %>%
          select(
            sp_english, sp_code,
            block_size_km, rep, fold, Atlas,
            n, rmse, mae, spearman, pearson, poisson_dev,
            mean_obs, mean_pred_mean, mean_pred_q025, mean_pred_q50, mean_pred_q975,
            mean_pred_yrep_mean, mean_pred_yrep_q025, mean_pred_yrep_q50, mean_pred_yrep_q975,
            coverage_95, mean_log_score, brier_ge1
          )
        
        append_rds_rows(
          summ_path,
          summ_row,
          key_cols = c("sp_code","rep","fold","block_size_km","Atlas")
        )
      }
      
      # Block-level summaries (analogous to 11a)
      test_df0 <- sf::st_drop_geometry(dat_test_sf) %>%
        transmute(block_id = block_id,
                  Atlas = as.character(Atlas),
                  count_obs = count)
      
      atlases_in_test <- sort(unique(test_df0$Atlas))
      block_summ_list <- list()
      
      for (atl in atlases_in_test) {
        idx_atl <- which(test_df0$Atlas == atl)
        if (!length(idx_atl)) next
        
        blocks_here <- unique(test_df0$block_id[idx_atl])
        
        for (bid in blocks_here) {
          idx <- which(test_df0$block_id == bid & test_df0$Atlas == atl)
          
          obs_all_zero <- as.integer(all(test_df0$count_obs[idx] == 0, na.rm = TRUE))
          obs_any_detect <- 1L - obs_all_zero
          
          mu_sum <- sum(mu_pred[idx], na.rm = TRUE)
          n_in_block <- length(idx)
          
          # Under independent Poisson per survey:
          p_all_zero <- exp(-mu_sum)
          p_any_detect <- 1 - p_all_zero
          
          mean_obs_block <- mean(test_df0$count_obs[idx], na.rm = TRUE)
          mean_pred_block <- mean(mu_pred[idx], na.rm = TRUE)
          
          # Distribution of mean count in block via sum ~ Poisson(mu_sum)
          mean_yrep_q025 <- qpois(0.025, lambda = mu_sum) / n_in_block
          mean_yrep_q50  <- qpois(0.50,  lambda = mu_sum) / n_in_block
          mean_yrep_q975 <- qpois(0.975, lambda = mu_sum) / n_in_block
          
          covered_95 <- as.integer(mean_obs_block >= mean_yrep_q025 & mean_obs_block <= mean_yrep_q975)
          
          block_summ_list[[length(block_summ_list) + 1]] <- tibble(
            sp_english = sp_english,
            sp_code = sp_code,
            block_size_km = block_size_km,
            rep = r,
            fold = f,
            Atlas = atl,
            block_id = bid,
            n = n_in_block,
            mean_obs = mean_obs_block,
            mean_pred_mean = mean_pred_block,
            mean_pred_q025 = NA_real_,
            mean_pred_q50  = NA_real_,
            mean_pred_q975 = NA_real_,
            mean_pred_yrep_mean = mean_pred_block,
            mean_pred_yrep_q025 = mean_yrep_q025,
            mean_pred_yrep_q50  = mean_yrep_q50,
            mean_pred_yrep_q975 = mean_yrep_q975,
            covered_95 = covered_95,
            p_all_zero = p_all_zero,
            p_any_detect = p_any_detect,
            obs_any_detect = obs_any_detect,
            rmse = sqrt(mean((mu_pred[idx] - test_df0$count_obs[idx])^2, na.rm = TRUE)),
            mae  = mean(abs(mu_pred[idx] - test_df0$count_obs[idx]), na.rm = TRUE),
            bias = mean(mu_pred[idx] - test_df0$count_obs[idx], na.rm = TRUE)
          )
        }
      }
      
      if (length(block_summ_list) > 0) {
        block_summ <- bind_rows(block_summ_list) %>%
          left_join(block_polys %>% st_drop_geometry(), by = "block_id") %>%
          left_join(block_polys %>% select(block_id, geometry), by = "block_id") %>%
          st_as_sf()
        
        append_rds_rows(
          block_path,
          block_summ,
          key_cols = c("sp_code","rep","fold","block_size_km","Atlas","block_id")
        )
      }
      
    } # folds
  } # repeats
  
  message("Done BRT CV: ", sp_english)
}

message("\n11b_crossvalidation_brt.R complete.")