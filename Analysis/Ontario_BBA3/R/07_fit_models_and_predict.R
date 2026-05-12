# ============================================================
# 09_fit_models_and_predict_PC_ARU.R
#
# Purpose:
#   Fit INLA/inlabru models for selected species and generate
#   predictions on the full OBBA2/OBBA3 1-km grids.
#
# Inputs:
#   data_clean/birds/data_ready_for_analysis.rds
#
# Outputs:
#   data_clean/model_output/predictions_<model_name>/<sp>_1km.rds
#   data_clean/model_output/summaries_<model_name>/model_summaries.rds
#   data_clean/model_output/summaries_<model_name>/change_summaries.rds
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
# Centralized paths
# ------------------------------------------------------------

source(here::here("R", "00_config_paths.R"))
source(file.path(paths$functions, "inla_model_utils2.R"))

# ------------------------------------------------------------
# Config
# ------------------------------------------------------------

model_name <- "m1"
rerun_models <- FALSE
rerun_predictions <- FALSE

n_prediction_draws <- 500
prediction_seed <- 123

in_file <- file.path(paths$data_clean, "birds", "data_ready_for_analysis.rds")
if (!file.exists(in_file)) {
  stop("Cannot find input at: ", in_file,
       "\nHave you run 07_filter_and_finalize_surveys.R?")
}

out_dir <- paths$model_output

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, paste0("models_", model_name)), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, paste0("predictions_", model_name)), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, paste0("summaries_", model_name)), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, paste0("data_used_", model_name)), recursive = TRUE, showWarnings = FALSE)

model_summaries_path  <- file.path(out_dir, paste0("summaries_", model_name), "model_summaries.rds")
change_summaries_path <- file.path(out_dir, paste0("summaries_", model_name), "change_summaries.rds")

model_summaries  <- load_or_empty_list(model_summaries_path)
change_summaries <- load_or_empty_list(change_summaries_path)

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

make_pred_grid <- function(grid_obba2, grid_obba3) {
  bind_rows(
    grid_obba2 %>% mutate(Atlas3 = 0L),
    grid_obba3 %>% mutate(Atlas3 = 1L)
  ) %>%
    mutate(
      Hours_Since_Sunrise = 0,
      days_midpoint = 0,
      Atlas3_c = Atlas3 - 0.5,
      BCR_idx = as.integer(BCR)
    )
}

predict_all_pixels <- function(mod, pred_grid, pred_formula,
                               n.samples, seed,
                               on_water_col = "on_water") {
  
  preds <- predict_inla(
    mod = mod,
    grid = pred_grid,
    pred_formula = pred_formula,
    n.samples = n.samples,
    seed = seed
  )
  
  # # Force deterministic zero abundance on open water
  # if (on_water_col %in% names(pred_grid)) {
  #   water_idx <- which(!is.na(pred_grid[[on_water_col]]) & pred_grid[[on_water_col]])
  #   if (length(water_idx) > 0) {
  #     preds$eta[water_idx, ] <- -Inf
  #   }
  # }
  
  idx2 <- which(pred_grid$Atlas == "OBBA2")
  idx3 <- which(pred_grid$Atlas == "OBBA3")
  
  eta <- preds$eta
  mu2 <- exp(eta[idx2, , drop = FALSE])
  mu3 <- exp(eta[idx3, , drop = FALSE])
  
  if (nrow(mu2) != nrow(mu3)) {
    stop("Prediction grid does not contain matched OBBA2/OBBA3 rows.")
  }
  
  list(
    eta = eta,
    mu2 = mu2,
    mu3 = mu3
  )
}

summarize_predictions <- function(mu2, mu3) {
  abs_change <- mu3 - mu2
  list(
    OBBA2 = summarize_posterior(mu2, CI_probs = c(0.05, 0.95), prefix = "OBBA2"),
    OBBA3 = summarize_posterior(mu3, CI_probs = c(0.05, 0.95), prefix = "OBBA3"),
    abs_change = summarize_posterior(abs_change, CI_probs = c(0.05, 0.95), prefix = "abs_change")
  )
}

make_hex_draws <- function(g2, mu2, mu3) {
  
  if (!("hex_id" %in% names(g2))) {
    stop("Cannot create hex_draws because `hex_id` is not present in prediction grid.")
  }
  
  hex_ids <- g2$hex_id
  u_hex <- unique(hex_ids)
  
  bind_rows(
    lapply(u_hex, function(hx) {
      
      idx_hex <- which(hex_ids == hx)
      
      mu2_hex <- colMeans(mu2[idx_hex, , drop = FALSE])
      mu3_hex <- colMeans(mu3[idx_hex, , drop = FALSE])
      
      tibble(
        hex_id = hx,
        n_pixels = length(idx_hex),
        mu_OBBA2 = list(mu2_hex),
        mu_OBBA3 = list(mu3_hex),
        abs_change = list(mu3_hex - mu2_hex)
      )
    })
  )
}

# ------------------------------------------------------------
# Load data
# ------------------------------------------------------------

dat <- readRDS(in_file)

all_surveys       <- dat$all_surveys
counts            <- dat$counts
grid_OBBA2        <- dat$grid_OBBA2
grid_OBBA3        <- dat$grid_OBBA3
study_boundary    <- dat$study_boundary %>% st_as_sf()
species_to_model  <- dat$species_to_model
hex_grid          <- dat$hex_grid

# Safe dates expanded to BCR should already exist from script 07
if ("safe_dates_bcr" %in% names(dat)) {
  safe_dates_bcr <- dat$safe_dates_bcr
} else {
  safe_dates_bcr <- dat$safe_dates %>%
    mutate(
      BCR = case_when(
        ecoregion == "Mixedwood Plains" ~ "13",
        ecoregion == "Hudson Plains"    ~ "7",
        ecoregion == "Boreal Shield"    ~ "8, 12",
        TRUE ~ NA_character_
      )
    ) %>%
    separate_rows(BCR, sep = ",") %>%
    mutate(BCR = as.integer(trimws(BCR))) %>%
    select(sp_english, BCR, start_doy, end_doy) %>%
    arrange(sp_english, BCR) %>%
    mutate(midpoint = (start_doy + end_doy) / 2)
}

# ------------------------------------------------------------
# Main loop
# ------------------------------------------------------------

base_covars <- c(
  "ForestNeedleleaf", "ForestBroadleaf", "ForestMixed", "Wetland", "Cropland",
  "Urban", "On_Road",
  "Grassland_BCR7_8", "Grassland_BCR12_13",
  "Shrubland_BCR13", "Shrubland_BCR7_8_12",
  "Lake_Lg", "Lake_Sm",
  "GreatLakes", "HudsonBayCoast",
  "River_Lg_BCR7", "River_Lg_BCR13_12_8",
  "River_Sm_BCR7", "River_Sm_BCR13_12_8"
)

min_detections <- 250
min_squares <- 50

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

species_run <- species_run %>%
  subset(english_name %in%
           c(
             "Belted Kingfisher",
             "Winter Wren",
             "Spotted Sandpiper",           # BAM makes different predictions in HBL
             "Wilson's Warbler",            # eBird, BAM, and Atlas make different predictions in HBL
             "Yellow-rumped Warbler",       # eBird and BAM make very different predictions in HBL
             "Wilson's Snipe",              # Atlas and eBird make different predictions in HBL
             "Solitary Sandpiper",          # Good example of eBird and Atlas showing different hotspot locations; BAM and eBird also make very different predictions
             "Savannah Sparrow",            # Check for correspondence with eBird in revised Atlas model
             "Philadelphia Vireo",          # BAM predicts some extreme hotspots in Northern Ontario
             "Hermit Thrush",               # Atlas shows HBL as hotspot, BAM shows almost none, eBird intermediate
             "Dark-eyed Junco",             # Atlas shows HBL as hotspot
             "Connecticut Warbler",
             "Black-and-white Warbler",
             "Bobolink",
             "Palm Warbler",
             "Sandhill Crane",
             "Bank Swallow",
             "Eastern Meadlowlark",
             "Northern Cardinal",
             "Olive-sided Flycatcher",
             "Dark-eyed Junco",
             "Osprey",
             "Rock Pigeon (Feral Pigeon)",
             "Double-crested Cormorant",
             "Swainson's Thrush",
             "American Crow",
             "American Goldfinch",
             "Bald Eagle",
             "Boreal Chickadee",
             "Common Yellowthroat",
             "Common Nighthawk"
           ))

print(species_run)
message("Species queued: ", nrow(species_run))

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
  
  model_path <- file.path(out_dir, paste0("models_", model_name), paste0(sp_file, "_1km.rds"))
  pred_path  <- file.path(out_dir, paste0("predictions_", model_name), paste0(sp_file, "_1km.rds"))
  dat_path   <- file.path(out_dir, paste0("data_used_", model_name), paste0(sp_file, "_1km.rds"))
  
  if (!(sp_code %in% names(counts))) {
    message("Skipping (species_id not found in counts columns): ", sp_code)
    next
  }
  
  sp_dat <- all_surveys %>%
    mutate(count = counts[[sp_code]])
  
  # ------------------------------------------------------------
  # Filter species data to species- and BCR-specific safe dates
  # ------------------------------------------------------------
  
  sp_safe_dates <- safe_dates_bcr %>%
    filter(sp_english == sp_name)
  
  if (nrow(sp_safe_dates) == 0) {
    pred_doy <- NA_real_
    warning(paste0("Species '", sp_name, "' has no BCR safe dates at all."), call. = FALSE)
  } else {
    fallback_start <- max(sp_safe_dates$start_doy, na.rm = TRUE)
    fallback_end   <- min(sp_safe_dates$end_doy,   na.rm = TRUE)
    
    if (fallback_start <= fallback_end) {
      all_bcrs <- sp_dat %>% distinct(BCR)
      safe_bcrs <- sp_safe_dates %>% distinct(BCR)
      missing_bcrs <- all_bcrs %>% anti_join(safe_bcrs, by = "BCR")
      
      if (nrow(missing_bcrs) > 0) {
        sp_safe_dates <- bind_rows(
          sp_safe_dates,
          missing_bcrs %>%
            mutate(
              sp_english = sp_name,
              start_doy = fallback_start,
              end_doy   = fallback_end,
              midpoint  = floor((fallback_start + fallback_end) / 2)
            ) %>%
            select(sp_english, BCR, start_doy, end_doy, midpoint)
        )
      }
      
      pred_doy <- floor((fallback_start + fallback_end) / 2)
    } else {
      pred_doy <- NA_real_
    }
  }
  
  # 
  sp_dat <- sp_dat %>%
    left_join(sp_safe_dates, by = "BCR") %>%
    filter(
      !is.na(start_doy),
      !is.na(end_doy),
      DayOfYear >= start_doy,
      DayOfYear <= end_doy
    ) %>%
    mutate(
      
      # days from safe date midpoint
      days_midpoint = DayOfYear - pred_doy,
      
      # difference in duration from standard 5-min survey
      duration_rescaled = Survey_Duration_Minutes - 5)
  
  # ------------------------------------------------------------
  # Determine approximately optimal time of day for detection
  # ------------------------------------------------------------
  
  hss_gam <- mgcv::gam(count ~ s(Hours_Since_Sunrise), data = sp_dat %>% subset(Survey_Type %in% c("Point_Count","ARU")), family = "poisson")
  
  hss_pred <- data.frame(Hours_Since_Sunrise = seq(min(sp_dat$Hours_Since_Sunrise)+0.5,max(sp_dat$Hours_Since_Sunrise)-0.5, length.out = 100) %>% round(2))
  hss_pred$pred <- predict(hss_gam, newdata = hss_pred, type = "response")
  optimal_HSS <- subset(hss_pred, pred == max(hss_pred$pred))$Hours_Since_Sunrise
  
  plot(pred~Hours_Since_Sunrise, data = hss_pred, type = "l", main = paste0(sp_name," - effect of Hours Since Sunrise\n\nOptimal survey timing = ",optimal_HSS," hours since sunrise"), lwd = 2)
  abline(v = optimal_HSS, lty = 2, col = "blue")
  
  # ------------------------------------------------------------
  # Error family triage
  # ------------------------------------------------------------
  
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
  
  # Use poisson as default.  Only consider nbinomial if posterior predictive checks suggest overdispersion
  error_family <- "poisson"
  # if (prop_total_top1pct_nonzero >= 0.10 || sum(y_pos > 15) > 10) {
  #   error_family <- "nbinomial"
  # }
  
  dat_for_review <- sp_dat %>%
    select(Date_Time, Survey_Type, count, Hours_Since_Sunrise, Atlas)
  
  # Save data that was used for this species model
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
  
  # ------------------------------------------------------------
  # Fit model
  # ------------------------------------------------------------
  
  if (file.exists(pred_path) && !rerun_predictions) {
    message("Predictions already exist for ", sp_name, "; skipping")
    next
  }
  
  covars_present <- intersect(base_covars, names(sp_dat))
  
  cov_sd <- sp_dat %>%
    as.data.frame() %>%
    select(covars_present) %>%
    apply(2, function(x) length(unique(x)))
  
  covars_present <- names(cov_sd)[which(cov_sd > 0)]
  
  # Priors for fixed effects (mean 0 and SD = 0.5 on log scale)
  cov_df_sp <- make_cov_df(covars_present, mean = 0, sd_linear = 0.5)
  
  # Priors for random effects; these are good defaults
  priors_list <- list(
    prior_range_abund = c(150,0.1),   #c(250, 0.5),
    prior_sigma_abund = c(0.5, 0.1),
    prior_range_change = c(100, 0.1), # c(100,0.1)
    prior_sigma_change = c(0.1, 0.1),
    prior_HSS_range = c(3, 0.9),
    prior_HSS_sigma = c(3, 0.1),
    prior_DOY_range_global = c(7, 0.9),
    prior_DOY_sigma_global = c(3, 0.1),
    kappa_pcprec_diff = c(log(2), 0.1)
  )
  
  # ------------------------------------------------------------
  # Prepare meshes
  # ------------------------------------------------------------
  
  hull <- fmesher::fm_extensions(
    study_boundary,
    convex  = c(50, 200),
    concave = c(10, 200)
  )
  
  # Finer mesh for shared abundance surface
  mesh_abund <- fmesher::fm_mesh_2d_inla(
    loc      = sf::st_as_sfc(sp_dat),
    boundary = hull,
    max.edge = c(40, 100),
    cutoff   = 40,
    crs      = sf::st_crs(sp_dat)
  )
  
  # Coarser mesh for change surface
  mesh_chg <- fmesher::fm_mesh_2d_inla(
    loc      = sf::st_as_sfc(sp_dat),
    boundary = hull,
    max.edge = c(40, 100),
    cutoff   = 40,
    crs      = sf::st_crs(sp_dat)
  )
  
  
  start_model <- Sys.time()
  
  # These allow the model to fit more quickly.  Use alternate settings if struggling with convergence
  int_strategy <- "eb"
  strategy <- "simplified.laplace" # alternate: laplace
  
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
      message("Model failed for ", sp_name, "; continuing.")
      next
    }
    
    # Do not save model, as file sizes are prohibitive
    # save_atomic(mod, model_path)
  }
  
  end_model <- Sys.time()
  fit_minutes <- round(as.numeric(end_model - start_model, units = "mins"))
  
  # Save summary of model effects
  model_summaries[[sp_name]] <- list(
    sp_name = sp_name,
    sp_code = sp_code,
    error_family = error_family,
    priors = priors_list,
    summary_fixed = mod$summary.fixed,
    summary_hyperpar = mod$summary.hyperpar
  )
  
  save_atomic(model_summaries, model_summaries_path)
  
  message(
    "\n====================\n",
    i, "/", nrow(species_run), ": ", sp_name,
    " (species_id = ", sp_code, "); ", fit_minutes, " min to fit model\n",
    "===================="
  )
  
  print(summary(mod))
  
  # ------------------------------------------------------------
  # Predictions onto full 1-km pixel grid for Atlas2 and Atlas3
  # ------------------------------------------------------------
  
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
    
    # ---- Multiply each prediction by (1-proportion open water) in pixel
    #      NOTE: THIS IMPLICITY ASSUMES THAT SURVEYS ARE NOT REPRESENTATIVE OF 
    #      OPEN-WATER DENSITIES; THESE ARE TERRESTRIAL PREDICTIONS ONLY
    #      NOTE2: THIS WILL ALSO MAKE OBSERVED COUNTS LOOK LESS LIKE PREDICTIONS
    #             IN AREAS WITH LOTS OF SMALL SCATTERED LAKES BECAUSE OBSERVED 
    #             COUNTS ARE ONLY TERRESTRIAL, BUT PREDICTIONS INCLUDE THE UNSAMPLED WATERBODIES
    
    preds$mu2 <- preds$mu2 * (1-grid_OBBA2$open_water)
    preds$mu3 <- preds$mu3 * (1-grid_OBBA3$open_water)
    
    # ---- Summarize predictions
    
    pred_summary <- summarize_predictions(preds$mu2,preds$mu3)
    
    g2 <- pred_grid %>%
      filter(Atlas == "OBBA2")
    
    g3 <- pred_grid %>%
      filter(Atlas == "OBBA3")
    
    preds_OBBA2_summary <- bind_cols(
      g2 %>%
        st_drop_geometry() %>%
        select(pixel_id, hex_id),
      pred_summary$OBBA2
    )
    
    preds_OBBA3_summary <- bind_cols(
      g3 %>%
        st_drop_geometry() %>%
        select(pixel_id, hex_id),
      pred_summary$OBBA3
    )
    
    preds_abs_change_summary <- bind_cols(
      g2 %>%
        st_drop_geometry() %>%
        select(pixel_id, hex_id),
      pred_summary$abs_change
    )
    
    # Store individual posterior draws only at the larger hexagon level
    # (i.e., sum pixel-level predictions within each hexagon within each posterior draw)
    hex_draws <- make_hex_draws(
      g2 = g2,
      mu2 = preds$mu2,
      mu3 = preds$mu3
    )
    
    end_prediction <- Sys.time()
    pred_minutes <- round(as.numeric(end_prediction - start_prediction, units = "mins"))
    
    sp_square_summary <- sp_dat %>%
      as.data.frame() %>%
      group_by(Atlas, square_id) %>%
      summarise(
        n_surveys    = n(),
        total_count  = sum(count),
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
        hex_draws = hex_draws
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

message("\n09_fit_models_and_predict_PC_ARU.R complete.")