# ============================================================
# 09_fit_models_and_predict_jointly.R
#
# Purpose:
#   Fit INLA/inlabru models for selected species and generate
#   predictions on the OBBA2/OBBA3 grids.
#
# Inputs:
#   data_clean/birds/data_ready_for_analysis.rds  (from script 07)
#
# Outputs (written incrementally per species):
#   data_clean/model_output/models/<sp>.rds
#   data_clean/model_output/predictions/<sp>.rds
#   data_clean/model_output/summaries/model_summaries.rds
#   data_clean/model_output/summaries/change_summaries.rds
#   data_clean/model_output/summaries/HSS_DOY_summaries.rds
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
})

# ------------------------------------------------------------
# Centralized paths
# ------------------------------------------------------------
source(here::here("R", "00_config_paths.R"))

# helper-functions
inla_utils_path <- file.path(paths$functions, "inla_model_utils.R")
source(inla_utils_path)

# ------------------------------------------------------------
# Config
# ------------------------------------------------------------

rerun_models <- FALSE
rerun_predictions <- FALSE

in_file <- file.path(paths$data_clean, "birds", "data_ready_for_analysis.rds")
if (!file.exists(in_file)) {
  stop("Cannot find input at: ", in_file,
       "\nHave you run 07_filter_and_finalize_surveys.R?")
}

out_dir <- paths$model_output
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "models_joint"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "predictions_joint"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "summaries_joint"), recursive = TRUE, showWarnings = FALSE)

# Summary caches
model_summaries_path  <- file.path(out_dir, "summaries_joint", "model_summaries.rds")
change_summaries_path <- file.path(out_dir, "summaries_joint", "change_summaries.rds")
hss_doy_path          <- file.path(out_dir, "summaries_joint", "HSS_DOY_summaries.rds")

model_summaries   <- load_or_empty_list(model_summaries_path)
change_summaries  <- load_or_empty_list(change_summaries_path)
hss_doy_summaries <- load_or_empty_list(hss_doy_path)

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

# Add Atlas labels to grids (helps later)
grid_OBBA2 <- grid_OBBA2 %>% mutate(Atlas = "OBBA2")
grid_OBBA3 <- grid_OBBA3 %>% mutate(Atlas = "OBBA3")

# ------------------------------------------------------------
# Create a hexagon grid across the study area, for assessing statistical significance
# of trends within discrete regions
# ------------------------------------------------------------

# Build hex grid
hex_grid <- make_hex_grid(study_boundary, width_km = 25) # Each hex ~ 540 km^2

# ------------------------------------------------------------
# Main loop (~ 45 min per species)
# ------------------------------------------------------------

# Covariates to *potentially* include; per-species covariate variance screening can happen later
base_covars <- c(
  "on_river",
  "on_road",
  "urban_3",
  "lc_1S","lc_1N",
  "lc_4S","lc_4N",
  "lc_5S","lc_5N",
  "lc_8S","lc_8N",
  "lc_9S","lc_9N",
  "lc_10S","lc_10N",
  "lc_11","lc_12","lc_14","lc_17",
  "insect_broadleaf","insect_needleleaf"
)

# If model fails to update after 15 min, assume it stalled and retry
timeout_min <- 15

# Select species to model
min_detections <- 150     # at least 150 detections in one atlas
min_squares <- 50         # at least 50 squares in one atlas

species_run <- species_to_model %>%
  filter((total_detections_OBBA3 >= min_detections |
           total_detections_OBBA2 >= min_detections )&
           (total_squares_OBBA3 >= min_squares |
           total_squares_OBBA2 >= min_squares))%>%
  na.omit()

species_to_check <- c(
  "Palm Warbler",
  "Bobolink",
  # "Blue Jay",
  "Canada Jay",
  "Olive-sided Flycatcher",
  # "Winter Wren",
  "Lesser Yellowlegs",
  "Blackpoll Warbler",
  "Connecticut Warbler",
  # "Lincoln's Sparrow",
  # "Fox Sparrow",
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
  # "White-throated Sparrow",
  "Bay-breasted Warbler"
)

species_run <- species_run %>%
  subset(english_name %in% species_to_check)

message("Species queued: ", nrow(species_run))
species_run

for (i in seq_len(nrow(species_run))) {
  
  start <- Sys.time()
  
  sp_english <- species_run$english_name[i]
  sp_code <- as.character(species_run$species_id[i])
  sp_file <- sp_filename(sp_english)
  
  message("\n====================\n", i, "/", nrow(species_run), ": ", sp_english, " (species_id = ", sp_code, ")\n====================")
  
  model_path <- file.path(out_dir, "models_joint", paste0(sp_file, ".rds"))
  pred_path  <- file.path(out_dir, "predictions_joint", paste0(sp_file, ".rds"))
  
  # --- Build species data
  if (!(sp_code %in% names(counts))) {
    message("Skipping (species_id not found in counts columns): ", sp_code)
    next
  }
  
  # Append counts
  sp_dat <- all_surveys %>%mutate(count = counts[[sp_code]])
  
  # --- Covariates spec
  covars_present <- intersect(base_covars, names(sp_dat))
  cov_df_sp <- make_cov_df(covars_present)
  
  # Priors controlling spatial autocorrelation fields
  prior_range_abund  <- c(500, 0.90) # 90% chance spatial autocorrelation range of shared field is less than 500 km
  prior_sigma_abund  <- c(3, 0.1)    # 10% chance SD is larger than 3
  prior_range_change <- c(500, 0.5)  # 50% chance spatial autocorrelation range of change field is less than 500 km
  prior_sigma_change <- c(0.1, 0.1)  # 10% chance SD is larger than 0.1
  
  # square identifier (persistent across atlases)
  sp_dat$square_idx <- as.numeric(as.factor(sp_dat$square_id))

  # --- Fit (or load existing)
  mod <- NULL
  if (file.exists(model_path) & rerun_models == FALSE) {
    mod <- readRDS(model_path)
    message("Loaded existing model: ", model_path)
  } else {
    mod <- try(
      
      fit_inla_multi_atlas(
        sp_dat = sp_dat,
        study_boundary = study_boundary,
        covariates = cov_df_sp,
        timeout_min = timeout_min,
        prior_range_abund = prior_range_abund,
        prior_sigma_abund = prior_sigma_abund,
        prior_range_change = prior_range_change,
        prior_sigma_change = prior_sigma_change,
        family = "poisson"),
      
      silent = TRUE
    )
    
    if (inherits(mod, "try-error") || is.null(mod)) {
      message("Model failed for ", sp_english, "; continuing.")
      next
    }
    
    save_atomic(mod, model_path)
    message("Saved model: ", model_path)
  }
  
  # --- Save minimal model summary (fast)
  model_summaries[[sp_english]] <- list(
    sp_english = sp_english,
    sp_code = sp_code,
    summary_fixed = mod$summary.fixed,
    summary_hyperpar = mod$summary.hyperpar
  )
  save_atomic(model_summaries, model_summaries_path)
  
  print(summary(mod))
  
  # --- Predictions (or load existing)
  
  preds_obj <- NULL
  if (file.exists(pred_path)  & rerun_predictions == FALSE) {
    preds_obj <- readRDS(pred_path)
    message("Loaded existing predictions: ", pred_path)
    
  } else {
    
    message("Generating predictions")
    
    pred_grid <- bind_rows(
      grid_OBBA2 %>% mutate(Atlas3 = 0),
      grid_OBBA3 %>% mutate(Atlas3 = 1)
    ) %>%
      mutate(
        # reference values for predictions
        Hours_Since_Sunrise = 0,
        days_rescaled = -7,
        Atlas3_c = Atlas3 - 0.5,
        BCR_factor = as.numeric(factor(BCR)),
        BCR_idx = as.integer(BCR)
      )
    
    pred_formula <- make_pred_formula_multiatlas(cov_df_sp)
    
    preds <- predict_inla(
      mod = mod,
      grid = pred_grid,
      pred_formula = pred_formula,
      n.samples = 1000,
      seed = 123
    )
    
    # preds$eta is n_grid x n_draw
    eta <- preds$eta
    
    idx2 <- which(pred_grid$Atlas == "OBBA2")
    idx3 <- which(pred_grid$Atlas == "OBBA3")
    
    mu2 <- exp(eta[idx2, , drop = FALSE])
    mu3 <- exp(eta[idx3, , drop = FALSE])
    abs_change <- mu3 - mu2
    
    # --- Summarize predictions for each pixel
    preds_OBBA2_summary <- summarize_posterior(mu2, CI_probs = c(0.05, 0.95), prefix = "OBBA2")
    preds_OBBA3_summary <- summarize_posterior(mu3, CI_probs = c(0.05, 0.95), prefix = "OBBA3")
    preds_abs_change_summary <- summarize_posterior(abs_change, CI_probs = c(0.05, 0.95), prefix = "abs_change")
    
    # -------------------------------------------------
    # Assess support for meaningful population change within larger polygons,
    # (e.g., for drawing "significance" boundaries on maps, or for summarizing results into BCRs)
    # -------------------------------------------------
    
    # Associate each pixel with a particular hexagon
    idx <- build_pixel_polygon_index(
      grid_sf      = grid_OBBA2,         # the full pixel grid
      polygons_sf  = hex_grid,           # e.g., BCRs/ecoregions/hexes
      poly_id_col  = "hex_id",
      join         = "within"
    )
    
    # Assess evidence of change within each polygon
    eta_draws_per_hex <- aggregate_polygon_draws_from_eta(
      eta         = preds$eta,
      pred_grid   = pred_grid,
      pix_poly_id = idx$pix_poly_id,
      poly_ids    = idx$poly_ids,
      poly_id_col = "hex_id",
      atlas_col   = "Atlas",
      atlas_levels = c("OBBA2", "OBBA3")
    )
    
    # Save posterior summaries
    save_atomic(list(
      sp_english = sp_english,
      sp_code = sp_code,
      OBBA2 = preds_OBBA2_summary,
      OBBA3 = preds_OBBA3_summary,
      abs_change = preds_abs_change_summary,
      
      # For plotting boundaries of important/significant population change
      hexagon_sf = hex_grid,
      eta_draws_per_hex = eta_draws_per_hex,
      
      meta = list(n_draws = ncol(eta))
    ), pred_path)
  }
  
  end <- Sys.time()
  message("Done: ", sp_english," - ",round(end-start,2)," min")
}

message("\n09_fit_models_and_predict_jointly.R complete.")
