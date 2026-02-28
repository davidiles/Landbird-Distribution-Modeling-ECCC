# ============================================================
# 09b_fit_models_and_predict_separately.R
#
# Purpose:
#   Fit INLA/inlabru models for selected species by FITTING EACH
#   ATLAS PERIOD SEPARATELY (OBBA2 and OBBA3), predicting onto each
#   atlas grid, and constructing a change surface by differencing
#   posterior predictions pixel-by-pixel.
#
#   This script exists mainly as a demonstration that the joint model
#   (09a) produces similar mean predictions, but with improved precision.
#
# Inputs:
#   data_clean/birds/data_ready_for_analysis.rds  (from script 07)
#   R/functions/inla_model_utils.R                (helper functions)
#
# Outputs (incrementally per species):
#   data_clean/model_output/models_separate/<sp>__OBBA2.rds
#   data_clean/model_output/models_separate/<sp>__OBBA3.rds
#   data_clean/model_output/predictions_separate/<sp>.rds
#
# Outputs (incrementally updated caches):
#   data_clean/model_output/summaries_separate/model_summaries.rds
#   data_clean/model_output/summaries_separate/change_summaries.rds
#   data_clean/model_output/summaries_separate/HSS_DOY_summaries.rds
#
# Assumptions:
#   - grid_OBBA2 and grid_OBBA3 have identical row ordering so that
#     row i corresponds to the same pixel/location across periods.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(sf)
  
  library(INLA)
  library(inlabru)
  library(fmesher)
})

# Helper functions (must provide: sp_filename, save_atomic, load_or_empty_list,
# make_cov_df, predict_inla, summarize_posterior, fit_inla_single_atlas)
source("R/functions/inla_model_utils.R")

# ------------------------------------------------------------
# Config
# ------------------------------------------------------------

rerun_models <- TRUE

in_file <- "data_clean/birds/data_ready_for_analysis.rds"

out_dir <- "data_clean/model_output"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "models_separate"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "predictions_separate"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "summaries_separate"), recursive = TRUE, showWarnings = FALSE)

# Species eligibility filters (OBBA3)
min_detections_obba3 <- 100
min_squares_obba3 <- 50

# Posterior draw count for predictions
n_samples_predict <- 1000

# If model takes more than this long to fit, assume it stalled and retry
timeout_min <- 15

# Spatial field priors (same for both separate fits)
prior_range_abund <- c(100, 0.50)   # 50% chance range < 100 km
prior_sigma_abund <- c(0.1, 0.05)   # 5% chance SD > 0.1

# Reference prediction settings (match 09a convention)
ref_HSS <- 0
ref_days_since_june15 <- -7

# Covariates to *potentially* include (subset to those present per species)
base_covars <- c(
  "on_river",
  "on_road",
  "urban_3",
  "lc_1","lc_4","lc_5",
  "lc_8S","lc_8N",
  "lc_9S","lc_9N",
  "lc_10S","lc_10N",
  "lc_11","lc_12","lc_14","lc_17"
)

# ------------------------------------------------------------
# Load data
# ------------------------------------------------------------

stopifnot(file.exists(in_file))
dat <- readRDS(in_file)

all_surveys <- dat$all_surveys
counts      <- dat$counts
grid_OBBA2  <- dat$grid_OBBA2
grid_OBBA3  <- dat$grid_OBBA3
study_boundary <- dat$study_boundary %>% st_as_sf()
species_to_model <- dat$species_to_model

stopifnot(nrow(all_surveys) == nrow(counts))

# For separate-fit differencing we require aligned grids
stopifnot(nrow(grid_OBBA2) == nrow(grid_OBBA3))

# Add atlas labels (useful for downstream plots/debugging)
grid_OBBA2 <- grid_OBBA2 %>% mutate(Atlas = "OBBA2")
grid_OBBA3 <- grid_OBBA3 %>% mutate(Atlas = "OBBA3")

# ------------------------------------------------------------
# Summary caches
# ------------------------------------------------------------

model_summaries_path  <- file.path(out_dir, "summaries_separate", "model_summaries.rds")
change_summaries_path <- file.path(out_dir, "summaries_separate", "change_summaries.rds")
hss_doy_path          <- file.path(out_dir, "summaries_separate", "HSS_DOY_summaries.rds")

model_summaries  <- load_or_empty_list(model_summaries_path)
change_summaries <- load_or_empty_list(change_summaries_path)
hss_doy_cache    <- load_or_empty_list(hss_doy_path)

# ------------------------------------------------------------
# Species list
# ------------------------------------------------------------

species_run <- species_to_model %>%
  filter(total_detections_OBBA3 >= min_detections_obba3,
         total_squares_OBBA3 >= min_squares_obba3) %>%
  tidyr::drop_na()

message("Species meeting thresholds: ", nrow(species_run))

# Optional: restrict to a test set (demo-only)
# species_to_check <- c("Canada Warbler")
# species_run <- species_run %>% filter(english_name %in% species_to_check)

message("Species queued: ", nrow(species_run))

# ------------------------------------------------------------
# Main loop
# ------------------------------------------------------------

for (i in seq_len(nrow(species_run))) {
  
  t0 <- Sys.time()
  
  sp_english <- species_run$english_name[i]
  sp_code    <- as.character(species_run$species_id[i])
  sp_file    <- sp_filename(sp_english)
  
  message("\n====================\n",
          i, "/", nrow(species_run), ": ", sp_english, " (", sp_code, ")\n",
          "====================")
  
  model_path2 <- file.path(out_dir, "models_separate", paste0(sp_file, "__OBBA2.rds"))
  model_path3 <- file.path(out_dir, "models_separate", paste0(sp_file, "__OBBA3.rds"))
  pred_path   <- file.path(out_dir, "predictions_separate", paste0(sp_file, ".rds"))
  
  # --- Build species analysis data
  if (!(sp_code %in% names(counts))) {
    message("Skipping (species_id not found in counts columns): ", sp_code)
    next
  }
  
  sp_dat_all <- all_surveys %>%
    mutate(
      count = counts[[sp_code]],
      days_since_june15 = DayOfYear - 166,
      BCR_factor = as.numeric(factor(BCR))
    )
  
  sp_dat2 <- sp_dat_all %>% filter(Atlas == "OBBA2")
  sp_dat3 <- sp_dat_all %>% filter(Atlas == "OBBA3")
  
  # --- Covariates spec (per species, based on survey table)
  covars_present <- intersect(base_covars, names(sp_dat_all))
  cov_df_sp <- make_cov_df(covars_present)
  
  # ----------------------------------------------------------
  # Fit/load models (OBBA2 and OBBA3)
  # ----------------------------------------------------------
  
  mod2 <- NULL
  if (file.exists(model_path2) && !rerun_models) {
    mod2 <- readRDS(model_path2)
    message("Loaded existing OBBA2 model: ", model_path2)
  } else {
    mod2 <- try(
      fit_inla_single_atlas(
        sp_dat = sp_dat2,
        study_boundary = study_boundary,
        covariates = cov_df_sp,
        timeout_min = timeout_min,
        prior_range_abund = prior_range_abund,
        prior_sigma_abund = prior_sigma_abund,
        family = "poisson"
      ),
      silent = TRUE
    )
    if (inherits(mod2, "try-error") || is.null(mod2)) {
      message("OBBA2 model failed for ", sp_english, "; continuing.")
      next
    }
    save_atomic(mod2, model_path2)
    message("Saved OBBA2 model: ", model_path2)
  }
  
  mod3 <- NULL
  if (file.exists(model_path3) && !rerun_models) {
    mod3 <- readRDS(model_path3)
    message("Loaded existing OBBA3 model: ", model_path3)
  } else {
    mod3 <- try(
      fit_inla_single_atlas(
        sp_dat = sp_dat3,
        study_boundary = study_boundary,
        covariates = cov_df_sp,
        timeout_min = timeout_min,
        prior_range_abund = prior_range_abund,
        prior_sigma_abund = prior_sigma_abund,
        family = "poisson"
      ),
      silent = TRUE
    )
    if (inherits(mod3, "try-error") || is.null(mod3)) {
      message("OBBA3 model failed for ", sp_english, "; continuing.")
      next
    }
    save_atomic(mod3, model_path3)
    message("Saved OBBA3 model: ", model_path3)
  }
  
  # --- Save minimal model summary (fast) for BOTH atlases
  model_summaries[[sp_english]] <- list(
    sp_english = sp_english,
    sp_code = sp_code,
    OBBA2 = list(summary_fixed = mod2$summary.fixed, summary_hyperpar = mod2$summary.hyperpar),
    OBBA3 = list(summary_fixed = mod3$summary.fixed, summary_hyperpar = mod3$summary.hyperpar)
  )
  save_atomic(model_summaries, model_summaries_path)
  
  # ----------------------------------------------------------
  # Predictions (or load)
  # ----------------------------------------------------------
  
  if (file.exists(pred_path) && !rerun_models) {
    preds_out <- readRDS(pred_path)
    message("Loaded existing predictions: ", pred_path)
  } else {
    
    message("Generating predictions (separate fits for each atlas)...")
    
    pred_grid2 <- grid_OBBA2 %>%
      mutate(Hours_Since_Sunrise = ref_HSS,
             days_since_june15   = ref_days_since_june15)
    
    pred_grid3 <- grid_OBBA3 %>%
      mutate(Hours_Since_Sunrise = ref_HSS,
             days_since_june15   = ref_days_since_june15)
    
    # NOTE: this assumes your single-atlas model uses the same component
    # naming as the joint model for HSS/DOY (HSS + DOY_global) AND that
    # predict_inla() uses the provided pred_formula appropriately.
    #
    # If you have make_pred_formula() in utils for covariates, you can
    # use it here too; otherwise this mirrors your original intent by
    # predicting the full linear predictor.
    pred_formula <- make_pred_formula_single(cov_df_sp)
    
    pred2 <- predict_inla(
      mod = mod2,
      grid = pred_grid2,
      pred_formula = pred_formula,
      n.samples = n_samples_predict,
      seed = 123
    )
    
    pred3 <- predict_inla(
      mod = mod3,
      grid = pred_grid3,
      pred_formula = pred_formula,
      n.samples = n_samples_predict,
      seed = 123
    )
    
    mu2 <- exp(pred2$eta)  # n_pixels x n_draws
    mu3 <- exp(pred3$eta)
    
    stopifnot(nrow(mu2) == nrow(mu3), ncol(mu2) == ncol(mu3))
    abs_change <- mu3 - mu2
    
    # Posterior summaries for mapping (no geometry assumptions inside)
    OBBA2_summary <- summarize_posterior(mu2, CI_probs = c(0.05, 0.95), prefix = "OBBA2")
    OBBA3_summary <- summarize_posterior(mu3, CI_probs = c(0.05, 0.95), prefix = "OBBA3")
    abs_change_summary <- summarize_posterior(abs_change, CI_probs = c(0.05, 0.95), prefix = "abs_change")
    
    # Overall change summaries (optionally mask water)
    mu2_use <- mu2
    mu3_use <- mu3
    if ("on_water" %in% names(grid_OBBA2) && any(grid_OBBA2$on_water, na.rm = TRUE)) {
      mu2_use[grid_OBBA2$on_water, ] <- 0
    }
    if ("on_water" %in% names(grid_OBBA3) && any(grid_OBBA3$on_water, na.rm = TRUE)) {
      mu3_use[grid_OBBA3$on_water, ] <- 0
    }
    
    sum2 <- colSums(mu2_use)
    sum3 <- colSums(mu3_use)
    pct_change <- ifelse(sum2 > 0, (sum3 - sum2) / sum2 * 100, NA_real_)
    
    overall_summary <- tibble(
      mean_change = mean(pct_change, na.rm = TRUE),
      median_change = median(pct_change, na.rm = TRUE),
      lower_change = unname(quantile(pct_change, 0.05, na.rm = TRUE)),
      upper_change = unname(quantile(pct_change, 0.95, na.rm = TRUE)),
      prob_decline = mean(pct_change <= 0, na.rm = TRUE),
      prob_decline_30 = mean(pct_change <= -30, na.rm = TRUE),
      n_draws_used = sum(!is.na(pct_change))
    )
    
    # Update cached change summaries (used by later scripts if desired)
    change_summaries[[sp_english]] <- list(
      sp_english = sp_english,
      sp_code = sp_code,
      overall_summary = overall_summary
    )
    save_atomic(change_summaries, change_summaries_path)
    
    # --------------------------------------------------------
    # HSS / DOY curves (for both atlases)
    # --------------------------------------------------------
    
    make_hss_doy_grid <- function(sp_dat, atlas_label) {
      HSS_pred <- tibble(
        Hours_Since_Sunrise = seq(min(sp_dat$Hours_Since_Sunrise, na.rm = TRUE),
                                  max(sp_dat$Hours_Since_Sunrise, na.rm = TRUE),
                                  length.out = 100),
        days_since_june15 = ref_days_since_june15,
        effect = "HSS",
        atlas = atlas_label
      )
      
      DOY_pred <- tibble(
        days_since_june15 = seq(min(sp_dat$days_since_june15, na.rm = TRUE),
                                max(sp_dat$days_since_june15, na.rm = TRUE),
                                by = 1),
        Hours_Since_Sunrise = ref_HSS,
        effect = "DOY",
        atlas = atlas_label
      )
      
      bind_rows(HSS_pred, DOY_pred)
    }
    
    curve_grid2 <- make_hss_doy_grid(sp_dat2, "OBBA2")
    curve_grid3 <- make_hss_doy_grid(sp_dat3, "OBBA3")
    
    curve_pred2 <- predict_inla(
      mod = mod2,
      grid = curve_grid2,
      pred_formula = as.formula("~ data.frame(eta = HSS + DOY_global)"),
      n.samples = n_samples_predict,
      seed = 123
    )
    
    curve_pred3 <- predict_inla(
      mod = mod3,
      grid = curve_grid3,
      pred_formula = as.formula("~ data.frame(eta = HSS + DOY_global)"),
      n.samples = n_samples_predict,
      seed = 123
    )
    
    q2 <- t(apply(exp(curve_pred2$eta), 1, quantile, probs = c(0.025, 0.5, 0.975)))
    q3 <- t(apply(exp(curve_pred3$eta), 1, quantile, probs = c(0.025, 0.5, 0.975)))
    
    curve_summ2 <- bind_cols(curve_grid2, as_tibble(q2, .name_repair = "minimal"))
    curve_summ3 <- bind_cols(curve_grid3, as_tibble(q3, .name_repair = "minimal"))
    
    names(curve_summ2)[(ncol(curve_summ2)-2):ncol(curve_summ2)] <- c("q025","q50","q975")
    names(curve_summ3)[(ncol(curve_summ3)-2):ncol(curve_summ3)] <- c("q025","q50","q975")
    
    hss_doy_cache[[sp_english]] <- list(
      sp_english = sp_english,
      sp_code = sp_code,
      ref_HSS = ref_HSS,
      ref_days_since_june15 = ref_days_since_june15,
      curves = bind_rows(curve_summ2, curve_summ3)
    )
    save_atomic(hss_doy_cache, hss_doy_path)
    
    # --------------------------------------------------------
    # Save per-species outputs (match 09a style)
    # --------------------------------------------------------
    
    preds_out <- list(
      sp_english = sp_english,
      sp_code = sp_code,
      n_draws = ncol(mu2),
      
      pred_grid_meta = list(
        ref_HSS = ref_HSS,
        ref_days_since_june15 = ref_days_since_june15
      ),
      
      # posterior draws (grids)
      OBBA2 = mu2,
      OBBA3 = mu3,
      abs_change = abs_change,
      
      # posterior summaries (NOTE: geometry is not bound here; downstream can bind to grid)
      OBBA2_summary = OBBA2_summary,
      OBBA3_summary = OBBA3_summary,
      abs_change_summary = abs_change_summary,
      
      # Ontario-wide summary (sums across pixels per draw)
      overall_summary = overall_summary,
      
      # provenance
      change_method = "separate_fit_difference_pixelwise"
    )
    
    save_atomic(preds_out, pred_path)
    message("Saved predictions: ", pred_path)
  }
  
  t1 <- Sys.time()
  elapsed_min <- as.numeric(difftime(t1, t0, units = "mins"))
  message("Done: ", sp_english, " - ", sprintf("%.2f", elapsed_min), " min")
}

message("\n09b_fit_models_and_predict_separately.R complete.")