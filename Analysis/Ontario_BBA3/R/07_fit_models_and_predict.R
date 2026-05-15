# ============================================================
# 07_fit_models_and_predict.R
#
# Purpose
#   Fit joint OBBA2/OBBA3 INLA/inlabru models for selected species and
#   generate prediction products on the full 1-km prediction grids.
#
# Main outputs
#   - predictions_<model_name>/<species>_1km.rds
#   - summaries_<model_name>/model_summaries.rds
#   - data_used_<model_name>/<species>_1km.rds
# ============================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(sf)
  library(ggplot2)
  library(INLA)
  library(inlabru)
  library(fmesher)
  library(here)
  library(mgcv)
})

# ============================================================
# 1. Paths, utilities, and global configuration
# ============================================================

source(here::here("R", "00_config_paths.R"))
source(file.path(paths$functions, "inla_model_utils_testing.R"))

model_name <- "model_CLeffort"

rerun_models <- FALSE
rerun_predictions <- FALSE

n_prediction_draws <- 500
prediction_seed <- 123

min_detections <- 250
min_squares <- 50

# Candidate fixed-effect covariates. Species-specific filtering below removes
# covariates that are absent from the data or have no variation after filtering.
base_covars <- c(
  "ForestNeedleleaf", "ForestBroadleaf", "ForestMixed", "Wetland", "Cropland",
  "Urban", "On_Road",
  "Grassland_BCR7_8", "Grassland_BCR12_13",
  "Shrubland_BCR13", "Shrubland_BCR7_8_12",
  "Lake_Lg", "Lake_Sm", "GreatLakes", "HudsonBayCoast",
  "River_Lg_BCR7", "River_Lg_BCR13_12_8",
  "River_Sm_BCR7", "River_Sm_BCR13_12_8"
)

# Priors passed to fit_inla_multi_atlas().
priors_list <- list(
  prior_range_abund = c(200, 0.1),      # 10% chance range is less than 200 km
  prior_sigma_abund = c(0.5, 0.1),      # 10% chance SD is larger than 0.5
  
  prior_range_change = c(200, 0.1),     # 10% chance range is less than 200 km
  prior_sigma_change = c(0.1, 0.1),     # 10% chance SD is larger than 0.1
  
  prior_HSS_range = c(2, 0.1),          # 10% chance range is less than 2 hours
  prior_HSS_sigma = c(1, 0.1),          # 10% chance SD is larger than 1
  
  prior_DOY_range_global = c(3, 0.1),  # 10% chance range is less than 3 days
  prior_DOY_sigma_global = c(1, 0.1),   # 10% chance SD is larger than 1
  
  kappa_pcprec_diff = c(log(2), 0.1)    # 10% chance square iid random effect SD is larger than log(2)
)

# INLA approximation settings used inside fit_inla_multi_atlas().
int_strategy <- "ccd"
strategy <- "laplace"

# ============================================================
# 2. Input/output locations
# ============================================================

in_file <- file.path(paths$data_clean, "birds", "data_ready_for_analysis.rds")

if (!file.exists(in_file)) {
  stop(
    "Cannot find input at: ", in_file,
    "\nHave you run 06_filter_and_finalize_surveys.R?"
  )
}

out_dir <- paths$model_output

pred_dir <- file.path(out_dir, paste0("predictions_", model_name))
summary_dir <- file.path(out_dir, paste0("summaries_", model_name))
data_used_dir <- file.path(out_dir, paste0("data_used_", model_name))

purrr::walk(
  c(out_dir, pred_dir, summary_dir, data_used_dir),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
)

model_summaries_path <- file.path(summary_dir, "model_summaries.rds")
model_summaries <- load_or_empty_list(model_summaries_path)

# ============================================================
# 3. Load finalized data
# ============================================================

dat <- readRDS(in_file)

all_surveys <- dat$all_surveys
counts <- dat$counts
grid_OBBA2 <- dat$grid_OBBA2
grid_OBBA3 <- dat$grid_OBBA3
study_boundary <- dat$study_boundary %>% sf::st_as_sf()
species_to_model <- dat$species_to_model
safe_dates_breeding <- dat$safe_dates_breeding

# ============================================================
# 4. Select species to model
# ============================================================

species_run <- species_to_model %>%
  filter(
    (detections_safe_OBBA2 >= min_detections |
       detections_safe_OBBA3 >= min_detections) &
      (n_squares_safe_OBBA2 >= min_squares |
         n_squares_safe_OBBA3 >= min_squares)
  ) %>%
  filter(
    !is.na(detections_safe_OBBA2),
    !is.na(detections_safe_OBBA3),
    !is.na(n_squares_safe_OBBA2),
    !is.na(n_squares_safe_OBBA3)
  ) %>%
  mutate(
    delta_dets_safe = log(detections_safe_OBBA3 / detections_safe_OBBA2),
    delta_squares_safe = log(n_squares_safe_OBBA3 / n_squares_safe_OBBA2)
  ) %>%
  arrange(abs(delta_squares_safe)) %>%
  # 
  # # Temporary/manual species subset for development.
  # # Remove this block to run all eligible species.
  filter(
    english_name %in% c(
      # "American Goldfinch",
      # "Bank Swallow",
      # "Belted Kingfisher",
      # "Bobolink",
      # "Palm Warbler",
      # "Belted Kingfisher",
      # "Winter Wren",
      # "Spotted Sandpiper",
      # "Wilson's Warbler",
      # "Yellow-rumped Warbler",
      # "Wilson's Snipe",
      # "Solitary Sandpiper",
      # "Savannah Sparrow",
      # "Philadelphia Vireo",
      # "Hermit Thrush",
      # "Connecticut Warbler",
      # "Black-and-white Warbler",
      # "Bobolink",
      # "Palm Warbler",
      # "Sandhill Crane",
      # "Bank Swallow",
      # "Eastern Meadowlark",
      # "Northern Cardinal",
      "Alder Flycatcher",
      "American Crow",
      "Black-and-white Warbler"
      # "Olive-sided Flycatcher"
      # "Dark-eyed Junco",
      # "Osprey",
      # "Rock Pigeon (Feral Pigeon)",
      # "Double-crested Cormorant",
      # "Swainson's Thrush",
      # "American Crow",
      # "American Goldfinch",
      # "Bald Eagle",
      # "Boreal Chickadee",
      # "Common Yellowthroat",
      # "Common Nighthawk"
    )
  )

print(species_run)
message("Species queued: ", nrow(species_run))

# ============================================================
# 5. Main species loop
# ============================================================

for (i in seq_len(nrow(species_run))) {
  
  # ----------------------------------------------------------
  # 5.1 Species identifiers and output paths
  # ----------------------------------------------------------
  
  sp_name <- species_run$english_name[i]
  sp_code <- as.character(species_run$species_id[i])
  sp_file <- sp_filename(sp_name)
  
  message(
    "\n====================\n",
    i, "/", nrow(species_run), ": ", sp_name,
    " (species_id = ", sp_code, ")\n",
    "===================="
  )
  
  pred_path <- file.path(pred_dir, paste0(sp_file, "_1km.rds"))
  dat_path <- file.path(data_used_dir, paste0(sp_file, "_1km.rds"))
  
  if (!(sp_code %in% names(counts))) {
    message("Skipping; species_id not found in counts columns: ", sp_code)
    next
  }
  
  if (file.exists(pred_path) && !rerun_predictions) {
    message("Predictions already exist for ", sp_name, "; skipping.")
    next
  }
  
  sp_dat <- all_surveys %>%
    mutate(count = counts[[sp_code]])
  
  # ----------------------------------------------------------
  # 5.2 Apply species-specific safe-date filtering
  # ----------------------------------------------------------
  # For each species, keep only surveys within the species-specific
  # safe-date window for each biological region.
  #
  # If a region has no species-specific safe dates, use the intersection of
  # safe dates across available regions as a fallback. If no valid fallback
  # exists, skip the species.
  
  sp_safe_dates <- safe_dates_breeding %>%
    filter(sp_english == sp_name)
  
  if (nrow(sp_safe_dates) == 0) {
    message("Skipping ", sp_name, "; no safe-date records found.")
    next
  }
  
  fallback_start <- max(sp_safe_dates$start_doy, na.rm = TRUE)
  fallback_end <- min(sp_safe_dates$end_doy, na.rm = TRUE)
  
  if (!is.finite(fallback_start) ||
      !is.finite(fallback_end) ||
      fallback_start > fallback_end) {
    message("Skipping ", sp_name, "; no valid shared safe-date window.")
    next
  }
  
  all_regions <- sp_dat %>%
    distinct(Biol_Region)
  
  safe_regions <- sp_safe_dates %>%
    distinct(Biol_Region)
  
  missing_regions <- all_regions %>%
    anti_join(safe_regions, by = "Biol_Region")
  
  if (nrow(missing_regions) > 0) {
    sp_safe_dates <- bind_rows(
      sp_safe_dates,
      missing_regions %>%
        mutate(
          sp_english = sp_name,
          start_doy = fallback_start,
          end_doy = fallback_end,
          midpoint = floor((fallback_start + fallback_end) / 2)
        ) %>%
        select(sp_english, Biol_Region, start_doy, end_doy, midpoint)
    )
  }
  
  # Prediction reference date: midpoint of the shared safe-date window.
  pred_doy <- floor((fallback_start + fallback_end) / 2)
  
  sp_dat <- sp_dat %>%
    left_join(sp_safe_dates, by = "Biol_Region") %>%
    filter(
      !is.na(start_doy),
      !is.na(end_doy),
      DayOfYear >= start_doy,
      DayOfYear <= end_doy
    ) %>%
    mutate(
      days_midpoint = DayOfYear - pred_doy,
      duration_rescaled = Survey_Duration_Minutes - 5
    )
  
  if (nrow(sp_dat) == 0) {
    message("Skipping ", sp_name, "; no surveys remain after safe-date filtering.")
    next
  }
  
  # ----------------------------------------------------------
  # 5.3 Remove extreme counts and choose likelihood family
  # ----------------------------------------------------------
  # Very large counts are removed before model fitting. The negative-binomial
  # switch is currently disabled, but the diagnostic is retained for review.
  
  large_counts <- sp_dat %>%
    filter(count > 50)
  
  sp_dat <- sp_dat %>%
    filter(count <= 50)
  
  n_det <- sum(sp_dat$count > 0)
  y_pos <- sp_dat$count[sp_dat$count > 0]
  
  if (length(y_pos) == 0) {
    message("Skipping ", sp_name, "; no positive counts remain after filtering.")
    next
  }
  
  n_top <- max(1, ceiling(0.01 * n_det))
  
  prop_total_top1pct_nonzero <- sum(
    sort(y_pos, decreasing = TRUE)[seq_len(n_top)]
  ) / sum(y_pos)
  
  error_family <- "poisson"
  
  # Optional future switch:
  if (prop_total_top1pct_nonzero >= 0.10 || sum(y_pos > 15) > 10) {
    error_family <- "nbinomial"
  }
  
  # ----------------------------------------------------------
  # 5.4 Estimate reference Hours Since Sunrise
  # ----------------------------------------------------------
  # Fit a simple GAM to estimate the survey time with highest expected count.
  # This value is used as the standardized HSS value for prediction.
  
  hss_dat <- sp_dat %>%
    filter(
      Survey_Type %in% c("Point_Count", "ARU"),
      is.finite(Hours_Since_Sunrise)
    )
  
  if (nrow(hss_dat) < 50 || dplyr::n_distinct(hss_dat$Hours_Since_Sunrise) < 5) {
    optimal_HSS <- median(sp_dat$Hours_Since_Sunrise, na.rm = TRUE)
    message("Using median HSS for ", sp_name, ": ", round(optimal_HSS, 2))
  } else {
    hss_gam <- mgcv::gam(
      count ~ s(Hours_Since_Sunrise),
      data = hss_dat,
      family = "poisson"
    )
    
    hss_min <- quantile(hss_dat$Hours_Since_Sunrise, 0.05, na.rm = TRUE)
    hss_max <- quantile(hss_dat$Hours_Since_Sunrise, 0.95, na.rm = TRUE)
    
    hss_pred <- data.frame(
      Hours_Since_Sunrise = seq(hss_min, hss_max, length.out = 100)
    )
    
    hss_pred$pred <- predict(hss_gam, newdata = hss_pred, type = "response")
    
    optimal_HSS <- hss_pred$Hours_Since_Sunrise[which.max(hss_pred$pred)]
    
    plot(
      pred ~ Hours_Since_Sunrise,
      data = hss_pred,
      type = "l",
      main = paste0(
        sp_name,
        " - effect of Hours Since Sunrise\n\nOptimal survey timing = ",
        round(optimal_HSS, 2),
        " hours since sunrise"
      ),
      lwd = 2
    )
    abline(v = optimal_HSS, lty = 2, col = "blue")
  }
  
  # ----------------------------------------------------------
  # 5.5 Save compact data record for review
  # ----------------------------------------------------------
  
  dat_for_review <- sp_dat %>%
    select(
      Date_Time,
      Survey_Type,
      count,
      Hours_Since_Sunrise,
      Survey_Duration_Minutes,
      Distance_Traveled_m,
      DayOfYear,
      days_midpoint,
      Atlas,
      Atlas3_c,
      square_id,
      square_atlas,
      BCR,
      Biol_Region
    )
  
  saveRDS(
    list(
      sp_english = sp_name,
      sp_code = sp_code,
      sp_safe_dates = sp_safe_dates,
      pred_doy = pred_doy,
      large_counts = large_counts,
      prop_total_top1pct_nonzero = prop_total_top1pct_nonzero,
      error_family = error_family,
      sp_dat = dat_for_review,
      optimal_HSS = optimal_HSS
    ),
    dat_path
  )
  
  # ----------------------------------------------------------
  # 5.6 Build species-specific covariate table
  # ----------------------------------------------------------
  # Retain only candidate covariates that:
  #   1. exist in this species dataset; and
  #   2. have more than one unique finite value after filtering.
  
  covars_present <- intersect(base_covars, names(sp_dat))
  
  if (length(covars_present) > 0) {
    cov_n_unique <- sp_dat %>%
      sf::st_drop_geometry() %>%
      select(all_of(covars_present)) %>%
      summarise(
        across(
          everything(),
          ~ dplyr::n_distinct(.x[is.finite(.x)], na.rm = TRUE)
        )
      ) %>%
      tidyr::pivot_longer(
        cols = everything(),
        names_to = "covariate",
        values_to = "n_unique"
      )
    
    covars_present <- cov_n_unique %>%
      filter(n_unique > 1) %>%
      pull(covariate)
  }
  
  cov_df_sp <- make_cov_df(
    covars_present,
    mean = 0,
    sd_linear = 0.5
  )
  
  # ----------------------------------------------------------
  # 5.7 Fit or load model
  # ----------------------------------------------------------
  # Compatible with the revised fit_inla_multi_atlas(), which creates meshes
  # internally and uses int_strategy / strategy via control.inla.
  
  start_model <- Sys.time()
  mod <- NULL
  
  mod <- try(
    fit_inla_multi_atlas(
      sp_dat = sp_dat,
      study_boundary = study_boundary,
      covariates = cov_df_sp,
      
      prior_range_abund = priors_list$prior_range_abund,
      prior_sigma_abund = priors_list$prior_sigma_abund,
      
      prior_range_change = priors_list$prior_range_change,
      prior_sigma_change = priors_list$prior_sigma_change,
      
      prior_HSS_range = priors_list$prior_HSS_range,
      prior_HSS_sigma = priors_list$prior_HSS_sigma,
      
      prior_DOY_range_global = priors_list$prior_DOY_range_global,
      prior_DOY_sigma_global = priors_list$prior_DOY_sigma_global,
      
      kappa_pcprec_diff = priors_list$kappa_pcprec_diff,
      
      int_strategy = int_strategy,
      strategy = strategy,
      family = error_family
    ),
    silent = TRUE
  )
  
  if (inherits(mod, "try-error") || is.null(mod)) {
    message("Model failed for ", sp_name, "; skipping this species.")
    print(mod)
    next
  }
  
  end_model <- Sys.time()
  fit_minutes <- round(as.numeric(end_model - start_model, units = "mins"), 1)
  
  model_summaries[[sp_name]] <- list(
    sp_name = sp_name,
    sp_code = sp_code,
    error_family = error_family,
    priors = priors_list,
    int_strategy = int_strategy,
    strategy = strategy,
    n_surveys = nrow(sp_dat),
    n_detections = n_det,
    n_covariates = nrow(cov_df_sp),
    covariates = cov_df_sp,
    pred_doy = pred_doy,
    optimal_HSS = optimal_HSS,
    fit_minutes = fit_minutes,
    summary_fixed = mod$summary.fixed,
    summary_hyperpar = mod$summary.hyperpar
  )
  
  save_atomic(model_summaries, model_summaries_path)
  
  print(summary(mod))
  
  message(
    "\n====================\n",
    i, "/", nrow(species_run), ": ", sp_name,
    " (species_id = ", sp_code, "); ", fit_minutes, " min to fit model\n",
    "===================="
  )
  
  # ----------------------------------------------------------
  # 5.8 Generate full-grid predictions
  # ----------------------------------------------------------
  # Predictions are standardized to:
  #   - optimal_HSS for Hours_Since_Sunrise
  #   - days_midpoint = 0, i.e. the shared safe-date midpoint
  
  start_prediction <- Sys.time()
  message("Generating predictions for full 1-km grid")
  
  pred_formula <- make_pred_formula_multiatlas(cov_df_sp)
  
  pred_grid <- make_pred_grid(grid_OBBA2, grid_OBBA3) %>%
    mutate(
      Hours_Since_Sunrise = optimal_HSS,
      days_midpoint = 0
    )
  
  preds <- predict_all_pixels(
    mod = mod,
    pred_grid = pred_grid,
    pred_formula = pred_formula,
    n.samples = n_prediction_draws,
    seed = prediction_seed
  )
  
  # ----------------------------------------------------------
  # 5.9 Apply terrestrial open-water correction
  # ----------------------------------------------------------
  # This scales predicted abundance by the non-open-water fraction of each
  # pixel. It assumes the model estimates terrestrial abundance/density.
  
  preds$mu2_Corrected_for_Water <- preds$mu2 * (1 - grid_OBBA2$open_water)
  preds$mu3_Corrected_for_Water <- preds$mu3 * (1 - grid_OBBA3$open_water)
  
  # ----------------------------------------------------------
  # 5.10 Summarize pixel predictions and hex draws
  # ----------------------------------------------------------
  
  pred_summary <- summarize_predictions(preds$mu2, preds$mu3)
  
  pred_summary_Corrected_for_Water <- summarize_predictions(
    preds$mu2_Corrected_for_Water,
    preds$mu3_Corrected_for_Water
  )
  
  g2 <- pred_grid %>% filter(Atlas == "OBBA2")
  g3 <- pred_grid %>% filter(Atlas == "OBBA3")
  
  preds_OBBA2_summary <- bind_cols(
    g2 %>% st_drop_geometry() %>% select(pixel_id, hex_id),
    pred_summary$OBBA2
  )
  
  preds_OBBA3_summary <- bind_cols(
    g3 %>% st_drop_geometry() %>% select(pixel_id, hex_id),
    pred_summary$OBBA3
  )
  
  preds_abs_change_summary <- bind_cols(
    g2 %>% st_drop_geometry() %>% select(pixel_id, hex_id),
    pred_summary$abs_change
  )
  
  hex_draws <- make_hex_draws(
    g2 = g2,
    mu2 = preds$mu2,
    mu3 = preds$mu3
  )
  
  preds_OBBA2_summary_Corrected_for_Water <- bind_cols(
    g2 %>% st_drop_geometry() %>% select(pixel_id, hex_id),
    pred_summary_Corrected_for_Water$OBBA2
  )
  
  preds_OBBA3_summary_Corrected_for_Water <- bind_cols(
    g3 %>% st_drop_geometry() %>% select(pixel_id, hex_id),
    pred_summary_Corrected_for_Water$OBBA3
  )
  
  preds_abs_change_summary_Corrected_for_Water <- bind_cols(
    g2 %>% st_drop_geometry() %>% select(pixel_id, hex_id),
    pred_summary_Corrected_for_Water$abs_change
  )
  
  hex_draws_Corrected_for_Water <- make_hex_draws(
    g2 = g2,
    mu2 = preds$mu2_Corrected_for_Water,
    mu3 = preds$mu3_Corrected_for_Water
  )
  
  end_prediction <- Sys.time()
  pred_minutes <- round(as.numeric(end_prediction - start_prediction, units = "mins"), 1)
  
  # ----------------------------------------------------------
  # 5.11 Summarize observed survey coverage by atlas square
  # ----------------------------------------------------------
  
  sp_square_summary <- sp_dat %>%
    as.data.frame() %>%
    group_by(Atlas, square_id) %>%
    summarise(
      n_surveys = n(),
      total_count = sum(count),
      n_detections = sum(count > 0),
      BCR = names(which.max(table(BCR))),
      .groups = "drop"
    )
  
  # ----------------------------------------------------------
  # 5.12 Save prediction products
  # ----------------------------------------------------------
  
  save_atomic(
    list(
      sp_name = sp_name,
      sp_code = sp_code,
      sp_safe_dates = sp_safe_dates,
      sp_square_summary = sp_square_summary,
      
      error_family = error_family,
      priors = priors_list,
      int_strategy = int_strategy,
      strategy = strategy,
      
      fit_minutes = fit_minutes,
      pred_minutes = pred_minutes,
      
      summary_fixed = mod$summary.fixed,
      summary_hyperpar = mod$summary.hyperpar,
      
      pred_doy = pred_doy,
      optimal_HSS = optimal_HSS,
      prediction_seed = prediction_seed,
      n_prediction_draws = n_prediction_draws,
      
      # Predictions for terrestrial habitats
      OBBA2 = preds_OBBA2_summary,
      OBBA3 = preds_OBBA3_summary,
      abs_change = preds_abs_change_summary,
      hex_draws = hex_draws,
      
      # Predictions that correct for open water proportions in each pixel
      OBBA2_Corrected_for_Water = preds_OBBA2_summary_Corrected_for_Water,
      OBBA3_Corrected_for_Water = preds_OBBA3_summary_Corrected_for_Water,
      abs_change_Corrected_for_Water = preds_abs_change_summary_Corrected_for_Water,
      hex_draws_Corrected_for_Water = hex_draws_Corrected_for_Water
    ),
    pred_path
  )
  
  message(
    "\n====================\n",
    i, "/", nrow(species_run), ": ", sp_name,
    " (species_id = ", sp_code, "); ", pred_minutes,
    " min to generate predictions\n",
    "===================="
  )
  
  rm(preds, pred_grid, g2, g3, hex_draws, hex_draws_Corrected_for_Water)
  gc(verbose = FALSE)
}

message("\n07_fit_models_and_predict.R complete.")