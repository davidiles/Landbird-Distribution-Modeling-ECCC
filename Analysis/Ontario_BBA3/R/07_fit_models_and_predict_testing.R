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

# ------------------------------------------------------------
# Paths and utilities
# ------------------------------------------------------------

source(here::here("R", "00_config_paths.R"))
source(file.path(paths$functions, "inla_model_utils_testing.R"))

# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------

model_name <- "m_multiscale"
rerun_models <- FALSE
rerun_predictions <- FALSE

n_prediction_draws <- 500
prediction_seed <- 123

min_detections <- 250
min_squares <- 50

# Fixed-effect covariates considered for each species. Covariates absent from
# the species data, or with no variation, are dropped inside the species loop.
base_covars <- c(
  "ForestNeedleleaf", "ForestBroadleaf", "ForestMixed", "Wetland", "Cropland",
  "Urban", "On_Road",
  "Grassland_BCR7_8", "Grassland_BCR12_13",
  "Shrubland_BCR13", "Shrubland_BCR7_8_12",
  "Lake_Lg", "Lake_Sm","GreatLakes", "HudsonBayCoast",
  "River_Lg_BCR7", "River_Lg_BCR13_12_8",
  "River_Sm_BCR7", "River_Sm_BCR13_12_8"
)

# Default priors for the model components. These are passed unchanged to
# fit_inla_multi_atlas().
priors_list <- list(
  prior_range_abund_coarse  = c(200, 0.1),   # 10% chance range < 200 km
  prior_sigma_abund_coarse  = c(3,   0.1),   # 10% chance SD > 3
  
  prior_range_abund_fine  = c(200, 0.9),     # 90% chance range < 200 km
  prior_sigma_abund_finer  = c(0.5,   0.1),  # 10% chance SD > 0.5
  
  prior_range_change = c(200, 0.1),   # 10% chance range < 200 km
  prior_sigma_change = c(0.1, 0.1),   # 10% chance SD > 0.1
  
  prior_HSS_range = c(3, 0.9),
  prior_HSS_sigma = c(3, 0.1),
  prior_DOY_range_global = c(7, 0.9),
  prior_DOY_sigma_global = c(3, 0.1),
  kappa_pcprec_diff = c(log(2), 0.1)
)

# INLA settings. These are relatively fast defaults; use fuller Laplace settings
# if convergence or approximation quality becomes a concern.
int_strategy <- "ccd" # eb
strategy <- "laplace" # simplified.laplace

# ------------------------------------------------------------
# Input and output locations
# ------------------------------------------------------------

in_file <- file.path(paths$data_clean, "birds", "data_ready_for_analysis.rds")
if (!file.exists(in_file)) {
  stop(
    "Cannot find input at: ", in_file,
    "\nHave you run 07_filter_and_finalize_surveys.R?"
  )
}

out_dir <- paths$model_output

model_dir <- file.path(out_dir, paste0("models_", model_name))
pred_dir <- file.path(out_dir, paste0("predictions_", model_name))
summary_dir <- file.path(out_dir, paste0("summaries_", model_name))
data_used_dir <- file.path(out_dir, paste0("data_used_", model_name))

purrr::walk(
  c(out_dir, model_dir, pred_dir, summary_dir, data_used_dir),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
)

model_summaries_path <- file.path(summary_dir, "model_summaries.rds")
model_summaries <- load_or_empty_list(model_summaries_path)

# ------------------------------------------------------------
# Load finalized survey and prediction-grid data
# ------------------------------------------------------------

dat <- readRDS(in_file)

all_surveys <- dat$all_surveys
counts <- dat$counts
grid_OBBA2 <- dat$grid_OBBA2
grid_OBBA3 <- dat$grid_OBBA3
study_boundary <- dat$study_boundary %>% st_as_sf()
species_to_model <- dat$species_to_model
safe_dates_breeding <- dat$safe_dates_breeding

# ------------------------------------------------------------
# Select species to model
# ------------------------------------------------------------

set.seed(123)

species_run <- species_to_model %>%
  filter(
    (detections_safe_OBBA2 >= min_detections |
       detections_safe_OBBA3 >= min_detections) &
      (n_squares_safe_OBBA2 >= min_squares |
         n_squares_safe_OBBA3 >= min_squares)
  ) %>%
  na.omit() %>%
  mutate(
    delta_dets_safe = log(detections_safe_OBBA3 / detections_safe_OBBA2),
    delta_squares_safe = log(n_squares_safe_OBBA3 / n_squares_safe_OBBA2)
  ) %>%
  arrange(abs(delta_squares_safe))

# Temporary/manual species subset for development and review runs.
# Remove this block to run all species passing the filters above.
species_run <- species_run %>%
  subset(english_name %in% c(
    "Belted Kingfisher",
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
    # "Eastern Meadlowlark",
    # "Northern Cardinal",
    # "Olive-sided Flycatcher",
    # "Dark-eyed Junco",
    # "Osprey",
    # "Rock Pigeon (Feral Pigeon)",
    # "Double-crested Cormorant",
    # "Swainson's Thrush",
    # "American Crow",
    "American Goldfinch"
    # "Bald Eagle",
    # "Boreal Chickadee",
    # "Common Yellowthroat",
    # "Common Nighthawk"
  ))

print(species_run)
message("Species queued: ", nrow(species_run))

# ============================================================
# Main species loop
# ============================================================

for (i in seq_len(nrow(species_run))) {
  
  sp_name <- species_run$english_name[i]
  sp_code <- as.character(species_run$species_id[i])
  sp_file <- sp_filename(sp_name)
  
  message(
    "\n====================\n",
    i, "/", nrow(species_run), ": ", sp_name,
    " (species_id = ", sp_code, ")\n",
    "===================="
  )
  
  model_path <- file.path(model_dir, paste0(sp_file, "_1km.rds"))
  pred_path <- file.path(pred_dir, paste0(sp_file, "_1km.rds"))
  dat_path <- file.path(data_used_dir, paste0(sp_file, "_1km.rds"))
  
  if (!(sp_code %in% names(counts))) {
    message("Skipping (species_id not found in counts columns): ", sp_code)
    next
  }
  
  sp_dat <- all_surveys %>%
    mutate(count = counts[[sp_code]])
  
  # ----------------------------------------------------------
  # Skip if prediction products already exist
  # ----------------------------------------------------------
  
  if (file.exists(pred_path) && !rerun_predictions) {
    message("Predictions already exist for ", sp_name, "; skipping")
    next
  }
  
  # ----------------------------------------------------------
  # Apply species-specific safe-date filtering
  # ----------------------------------------------------------
  # Candidate for utility function:
  #   prepare_species_safe_dates(sp_name, sp_dat, safe_dates_breeding)
  # It would return sp_safe_dates and pred_doy.
  
  sp_safe_dates <- safe_dates_breeding %>%
    filter(sp_english == sp_name)
  
  if (nrow(sp_safe_dates) == 0) {
    pred_doy <- NA_real_
    warning(
      paste0("Species '", sp_name, "' has no BCR safe dates at all."),
      call. = FALSE
    )
  } else {
    fallback_start <- max(sp_safe_dates$start_doy, na.rm = TRUE)
    fallback_end <- min(sp_safe_dates$end_doy, na.rm = TRUE)
    
    if (fallback_start <= fallback_end) {
      all_Regions <- sp_dat %>% distinct(Biol_Region)
      safe_Regions <- sp_safe_dates %>% distinct(Biol_Region)
      missing_Regions <- all_Regions %>% anti_join(safe_Regions, by = "Biol_Region")
      
      if (nrow(missing_Regions) > 0) {
        sp_safe_dates <- bind_rows(
          sp_safe_dates,
          missing_Regions %>%
            mutate(
              sp_english = sp_name,
              start_doy = fallback_start,
              end_doy = fallback_end,
              midpoint = floor((fallback_start + fallback_end) / 2)
            ) %>%
            select(sp_english, Biol_Region, start_doy, end_doy, midpoint)
        )
      }
      
      # Predict to a shared safe-date midpoint across all regions.
      pred_doy <- floor((fallback_start + fallback_end) / 2)
    } else {
      pred_doy <- NA_real_
    }
  }
  
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
  
  # ----------------------------------------------------------
  # Estimate species-specific optimal survey timing
  # ----------------------------------------------------------
  # Candidate for utility function:
  #   estimate_optimal_hss(sp_dat, sp_name, plot = TRUE)
  
  hss_gam <- mgcv::gam(
    count ~ s(Hours_Since_Sunrise),
    data = sp_dat %>% subset(Survey_Type %in% c("Point_Count", "ARU")),
    family = "poisson"
  )
  
  hss_pred <- data.frame(
    Hours_Since_Sunrise = seq(
      min(sp_dat$Hours_Since_Sunrise) + 0.5,
      max(sp_dat$Hours_Since_Sunrise) - 3,
      length.out = 100
    ) %>%
      round(2)
  )
  hss_pred$pred <- predict(hss_gam, newdata = hss_pred, type = "response")
  optimal_HSS <- subset(hss_pred, pred == max(hss_pred$pred))$Hours_Since_Sunrise
  
  plot(
    pred ~ Hours_Since_Sunrise,
    data = hss_pred,
    type = "l",
    main = paste0(
      sp_name,
      " - effect of Hours Since Sunrise\n\nOptimal survey timing = ",
      optimal_HSS,
      " hours since sunrise"
    ),
    lwd = 2
  )
  abline(v = optimal_HSS, lty = 2, col = "blue")
  
  # ----------------------------------------------------------
  # Count screening and likelihood-family choice
  # ----------------------------------------------------------
  # Counts > 50 are removed before modeling. The negative-binomial switch is
  # currently disabled, so all fitted models use Poisson unless changed below.
  
  large_counts <- sp_dat %>% filter(count > 50)
  sp_dat <- sp_dat %>% filter(count <= 50)
  
  n_det <- sum(sp_dat$count > 0)
  y_pos <- sp_dat$count[sp_dat$count > 0]
  
  if (length(y_pos) == 0) {
    message("Skipping ", sp_name, " because there are no positive counts after filtering.")
    next
  }
  
  n_top <- max(1, ceiling(0.01 * n_det))
  prop_total_top1pct_nonzero <- sum(sort(y_pos, decreasing = TRUE)[seq_len(n_top)]) / sum(y_pos)
  
  error_family <- "poisson"
  # if (prop_total_top1pct_nonzero >= 0.10 || sum(y_pos > 15) > 10) {
  #   error_family <- "nbinomial"
  # }
  
  # Save a compact record of the data used for this species model.
  dat_for_review <- sp_dat %>%
    select(Date_Time, Survey_Type, count, Hours_Since_Sunrise, Atlas)
  
  saveRDS(
    list(
      sp_english = sp_name,
      sp_code = sp_code,
      sp_safe_dates = sp_safe_dates,
      sp_dat = dat_for_review,
      optimal_HSS = optimal_HSS
    ),
    dat_path
  )
  
  # ----------------------------------------------------------
  # Species-specific covariate table
  # ----------------------------------------------------------
  
  covars_present <- intersect(base_covars, names(sp_dat))
  
  cov_sd <- sp_dat %>%
    as.data.frame() %>%
    select(covars_present) %>%
    apply(2, function(x) length(unique(x)))
  
  covars_present <- names(cov_sd)[which(cov_sd > 0)]
  
  # Fixed-effect priors: mean 0, SD 0.5 on the log scale.
  cov_df_sp <- make_cov_df(covars_present, mean = 0, sd_linear = 0.5)
  
  # ----------------------------------------------------------
  # Fit or load model
  # ----------------------------------------------------------
  
  start_model <- Sys.time()
  mod <- NULL
  
  if (file.exists(model_path) && !rerun_models) {
    mod <- readRDS(model_path)
    message("Loaded existing model: ", model_path)
  } else {
    mod <- try(
      fit_inla_multi_atlas(
        sp_dat = sp_dat,
        study_boundary = study_boundary,
        covariates = cov_df_sp,
        mesh_abund = mesh_abund,
        mesh_chg = mesh_chg,
        
        prior_range_abund_coarse = priors_list$prior_range_abund_coarse,
        prior_sigma_abund_coarse = priors_list$prior_sigma_abund_coarse,
        prior_range_abund_fine = priors_list$prior_range_abund_fine,
        prior_sigma_abund_fine = priors_list$prior_sigma_abund_fine,
        
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
      message("Model failed for ", sp_name, "; continuing.")
      next
    }
    
    # Models are not saved because fitted inlabru objects are very large.
    # save_atomic(mod, model_path)
  }
  
  end_model <- Sys.time()
  fit_minutes <- round(as.numeric(end_model - start_model, units = "mins"))
  
  model_summaries[[sp_name]] <- list(
    sp_name = sp_name,
    sp_code = sp_code,
    error_family = error_family,
    priors = priors_list,
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
  # Predict on the full 1-km grid
  # ----------------------------------------------------------
  # Candidate for utility function:
  #   build_prediction_products(...)
  # This block creates all saved prediction-summary objects and should be moved
  # only after testing object equality with the current script outputs.
  
  if (file.exists(pred_path) && !rerun_predictions) {
    message("Loaded existing predictions: ", pred_path)
  } else {
    
    start_prediction <- Sys.time()
    message("Generating predictions for full 1-km grid")
    
    pred_formula <- make_pred_formula_multiatlas(cov_df_sp)
    pred_grid <- make_pred_grid(grid_OBBA2, grid_OBBA3) %>%
      mutate(Hours_Since_Sunrise = optimal_HSS)
    
    preds <- predict_all_pixels(
      mod = mod,
      pred_grid = pred_grid,
      pred_formula = pred_formula,
      n.samples = n_prediction_draws,
      seed = prediction_seed
    )
    
    # Terrestrial correction: multiply predictions by non-open-water area.
    # This assumes survey counts represent terrestrial density only. It can make
    # observed counts and pixel predictions diverge in cells with many small lakes.
    preds$mu2_Corrected_for_Water <- preds$mu2 * (1 - grid_OBBA2$open_water)
    preds$mu3_Corrected_for_Water <- preds$mu3 * (1 - grid_OBBA3$open_water)
    
    pred_summary <- summarize_predictions(preds$mu2, preds$mu3)
    pred_summary_Corrected_for_Water <- summarize_predictions(
      preds$mu2_Corrected_for_Water,
      preds$mu3_Corrected_for_Water
    )
    
    g2 <- pred_grid %>% filter(Atlas == "OBBA2")
    g3 <- pred_grid %>% filter(Atlas == "OBBA3")
    
    # Pixel-level summaries, uncorrected for open water (should be more similar
    # to observed counts, since surveys are only collected on land)
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
    
    # Pixel-level summaries after correcting for open water.
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
    pred_minutes <- round(as.numeric(end_prediction - start_prediction, units = "mins"))
    
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
    
    save_atomic(
      list(
        sp_name = sp_name,
        sp_code = sp_code,
        sp_safe_dates = sp_safe_dates,
        sp_square_summary = sp_square_summary,
        error_family = error_family,
        priors = priors_list,
        fit_minutes = fit_minutes,
        summary_fixed = mod$summary.fixed,
        summary_hyperpar = mod$summary.hyperpar,
        pred_minutes = pred_minutes,
        pred_doy = pred_doy,
        prediction_seed = prediction_seed,
        n_prediction_draws = n_prediction_draws,
        
        OBBA2 = preds_OBBA2_summary,
        OBBA3 = preds_OBBA3_summary,
        abs_change = preds_abs_change_summary,
        hex_draws = hex_draws,
        
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
    
    rm(preds, pred_grid, g2, g3, hex_draws)
    gc(verbose = FALSE)
  }
}

message("\n07_fit_models_and_predict.R complete.")
