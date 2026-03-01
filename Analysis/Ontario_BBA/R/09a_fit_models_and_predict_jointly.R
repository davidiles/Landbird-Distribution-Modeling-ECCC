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
})

# helper-functions
source("R/functions/0_inla_model_utils.R")  

# ------------------------------------------------------------
# Config
# ------------------------------------------------------------
rerun_models = TRUE

in_file <- "data_clean/birds/data_ready_for_analysis.rds"

out_dir <- "data_clean/model_output"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "models_joint"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "predictions_joint"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "summaries_joint"), recursive = TRUE, showWarnings = FALSE)

# Which species to run
min_detections_obba3 <- 100
min_squares_obba3 <- 50
n_samples_predict <- 1000

# If model takes more than 15 min to fit, assume it stalled and retry
timeout_min <- 15

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

# Summary caches
model_summaries_path  <- file.path(out_dir, "summaries_joint", "model_summaries.rds")
change_summaries_path <- file.path(out_dir, "summaries_joint", "change_summaries.rds")
hss_doy_path          <- file.path(out_dir, "summaries_joint", "HSS_DOY_summaries.rds")

model_summaries  <- load_or_empty_list(model_summaries_path)
change_summaries <- load_or_empty_list(change_summaries_path)
hss_doy_summaries <- load_or_empty_list(hss_doy_path)

# ------------------------------------------------------------
# Create a hexagon grid across the study area, for assessing statistical significance
# of trends
# ------------------------------------------------------------

make_hex_grid <- function(study_boundary, width_km = 25, seed = 123) {
  set.seed(123)
  hex <- sf::st_make_grid(
    study_boundary,
    cellsize = width_km ,
    square = FALSE
  )
  
  hex_sf <- sf::st_sf(
    hex_id = seq_along(hex),
    geometry = hex
  )
  
  # Keep only hexes that overlap the study boundary
  hex_sf <- hex_sf[sf::st_intersects(hex_sf, study_boundary, sparse = FALSE), , drop = FALSE]
  
  hex_sf
}

# Build hex grid
hex_grid <- make_hex_grid(study_boundary, width_km = 50) # Each hex ~ 2000 km^2
ggplot(hex_grid)+geom_sf()

# Assign each pixel in grids to a hexagon
# 3) Assign each pixel to a hex
#    Using st_within; pixels outside hexes become NA (should be rare if grids cover boundary)
hex_grid <- hex_grid %>% st_transform(st_crs(grid_OBBA2))
pix_hex <- sf::st_join(
  sf::st_sf(pixel_id = seq_len(nrow(grid_OBBA2)), geometry = sf::st_geometry(grid_OBBA2)),
  hex_grid[, c("hex_id")],
  join = sf::st_within,
  left = TRUE
)
hex_id_by_pixel <- pix_hex$hex_id

# Which pixels are valid
valid_pix <- which(!is.na(hex_id_by_pixel))
hex_ids <- sort(unique(hex_id_by_pixel[valid_pix]))

# ------------------------------------------------------------
# Main loop (~ 45 min per species)
# ------------------------------------------------------------

# Select species list
species_run <- species_to_model %>%
  filter(total_detections_OBBA3 >= min_detections_obba3 &
           total_squares_OBBA3 >= min_squares_obba3) %>%
  na.omit()

message("Species queued: ", nrow(species_run))

# species_run <- sample_n(species_run, nrow(species_run))

species_to_check <- c(
  #"Palm Warbler",
  #"Bobolink",
  #"Blue Jay",
  "Canada Jay",
  "Olive-sided Flycatcher",
  #"Winter Wren",
  #"Lesser Yellowlegs",
  #"Blackpoll Warbler",
  "Connecticut Warbler",
  "Palm Warbler",
  #"Lincoln's Sparrow",
  #"Fox Sparrow",
  "Common Nighthawk",
  #"Long-eared Owl",
  #"American Tree Sparrow",
  #"LeConte's Sparrow",
  #"Nelson's Sparrow",
  #"Boreal Chickadee",
  #"Rusty Blackbird",
  #"Yellow-bellied Flycatcher",
  "Greater Yellowlegs",
  #"Hudsonian Godwit",
  "Canada Warbler",
  #"Eastern Wood-Peewee",
  #"Grasshopper Sparrow",
  #"Solitary Sandpiper",
  #"White-throated Sparrow",
  "Bay-breasted Warbler"
)


species_run <- species_run %>%
  subset(english_name %in% species_to_check)

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

message("Species queued: ", nrow(species_run))

for (i in seq_len(nrow(species_run))) {
  
  start <- Sys.time()
  
  sp_english <- species_run$english_name[i]
  sp_code <- as.character(species_run$species_id[i])
  sp_file <- sp_filename(sp_english)
  
  message("\n====================\n", i, "/", nrow(species_run), ": ", sp_english, " (", sp_code, ")\n====================")
  
  model_path <- file.path(out_dir, "models_joint", paste0(sp_file, ".rds"))
  pred_path  <- file.path(out_dir, "predictions_joint", paste0(sp_file, ".rds"))
  
  # --- Build species data
  if (!(sp_code %in% names(counts))) {
    message("Skipping (species_id not found in counts columns): ", sp_code)
    next
  }
  
  sp_dat <- all_surveys %>%
    mutate(
      count = counts[[sp_code]],
      days_since_june15 = DayOfYear - 166,
      BCR_factor = as.numeric(factor(BCR)),
      BCR_idx = as.integer(BCR),
      Atlas3 = ifelse(Atlas == "OBBA2", 0, 1),
      Atlas3_c = Atlas3 - 0.5
    )
  
  # --- Covariates spec
  covars_present <- intersect(base_covars, names(sp_dat))
  cov_df_sp <- make_cov_df(covars_present)
  
  # Priors controlling spatial autocorrelation fields
  prior_range_abund  <- c(250, 0.50) # 50% chance spatial autocorrelation range of shared field is less than 250 km
  prior_sigma_abund  <- c(3, 0.1)    # 10% chance SD is larger than 3
  prior_range_change <- c(250, 0.1)  # 10% chance spatial autocorrelation range of change field is less than 250 km
  prior_sigma_change <- c(0.1, 0.1)  # 10% chance SD is larger than 0.1
  
  # Modify square_atlas so that it uses the same identifier for both atlases
  sp_dat$square_idx <- as.numeric(as.factor(sp_dat$square_id))
  
  # --- Fit (or load existing)
  mod <- NULL
  if (file.exists(model_path) & rerun_models == FALSE) {
    mod <- readRDS(model_path)
    message("Loaded existing model: ", model_path)
  } else {
    mod <- try(
      
      fit_inla_multi_atlas3(
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
  if (file.exists(pred_path)  & rerun_models == FALSE) {
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
        days_since_june15 = -7,
        Atlas3_c = Atlas3 - 0.5,
        BCR_factor = as.numeric(factor(BCR)),
        BCR_idx = as.integer(BCR)
      )
    
    pred_formula <- make_pred_formula_multiatlas(cov_df_sp)
    
    preds <- predict_inla(
      mod = mod,
      grid = pred_grid,
      pred_formula = pred_formula,
      n.samples = n_samples_predict,
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
    
    # Calculate proportional change in each polygon
    eps <- 1e-9  # prevents division by 0
    polygon_prop_by_draw <- lapply(hex_ids, function(h) {
      rows <- which(hex_id_by_pixel == h)
      Mu2 <- colSums(mu2[rows, , drop = FALSE])
      Mu3 <- colSums(mu3[rows, , drop = FALSE])
      (Mu3 - Mu2) / pmax(Mu2, eps)
    })
    
    # Compute probabilities + summaries per hex
    delta_pct <- 0.10 # Threshold for "significant effect"
    polygon_change_summary <- bind_rows(lapply(seq_along(hex_ids), function(k) {
      v <- polygon_prop_by_draw[[k]]
      tibble(
        hex_id = hex_ids[k],
        n_pixels = sum(hex_id_by_pixel == hex_ids[k], na.rm = TRUE),
        p_inc = mean(v >  delta_pct),
        p_dec = mean(v < -delta_pct),
        q50 = median(v),
        q05 = unname(quantile(v, 0.05)),
        q95 = unname(quantile(v, 0.95))
      )
    }))
    
    # Join summaries back onto hex geometry
    polygon_change_sf <- dplyr::left_join(hex_grid, polygon_change_summary, by = "hex_id")
    
    # # --- Overall + BCR change summaries (mask water if available)
    # # If on_water exists, mask
    # mu2_use <- preds_obj$OBBA2
    # mu3_use <- preds_obj$OBBA3
    # if ("on_water" %in% names(grid_OBBA2) && any(grid_OBBA2$on_water, na.rm = TRUE)) {
    #   mu2_use[grid_OBBA2$on_water, ] <- 0
    # }
    # if ("on_water" %in% names(grid_OBBA3) && any(grid_OBBA3$on_water, na.rm = TRUE)) {
    #   mu3_use[grid_OBBA3$on_water, ] <- 0
    # }
    # 
    # # Overall percent change across all pixels (using sums across pixels per draw)
    # sum2 <- colSums(mu2_use)
    # sum3 <- colSums(mu3_use)
    # pct_change <- (sum3 - sum2) / sum2 * 100
    # 
    # overall_summary <- tibble(
    #   mean_change = mean(pct_change),
    #   median_change = median(pct_change),
    #   lower_change = unname(quantile(pct_change, 0.05)),
    #   upper_change = unname(quantile(pct_change, 0.95)),
    #   prob_decline = mean(pct_change <= 0),
    #   prob_decline_30 = mean(pct_change <= -30)
    # )
    # 
    # # --- HSS / DOY curves (lightweight)
    # ref_doy <- -7
    # ref_hss <- 0
    # 
    # HSS_pred <- tibble(
    #   Hours_Since_Sunrise = seq(min(sp_dat$Hours_Since_Sunrise, na.rm = TRUE),
    #                             max(sp_dat$Hours_Since_Sunrise, na.rm = TRUE),
    #                             length.out = 100),
    #   days_since_june15 = ref_doy,
    #   effect = "HSS"
    # )
    # 
    # DOY_pred <- tibble(
    #   days_since_june15 = seq(min(sp_dat$days_since_june15, na.rm = TRUE),
    #                           max(sp_dat$days_since_june15, na.rm = TRUE),
    #                           by = 1),
    #   Hours_Since_Sunrise = ref_hss,
    #   effect = "DOY"
    # )
    # 
    # pred_df <- bind_rows(HSS_pred, DOY_pred)
    # 
    # # predict only these components
    # pred_all <- predict_inla(
    #   mod = mod,
    #   grid = pred_df,
    #   pred_formula = as.formula("~ data.frame(eta = HSS + DOY_global)"),
    #   n.samples = n_samples_predict,
    #   seed = 123
    # )
    # 
    # # summarize on response scale
    # q <- t(apply(exp(pred_all$eta), 1, quantile, probs = c(0.025, 0.5, 0.975)))
    # hss_doy_summaries <- bind_cols(pred_df, as_tibble(q, .name_repair = "minimal"))
    # names(hss_doy_summaries)[(ncol(hss_doy_summaries)-2):ncol(hss_doy_summaries)] <- c("q025","q50","q975")
    # 
    # Save posterior summaries
    save_atomic(list(
      sp_english = sp_english,
      sp_code = sp_code,
      OBBA2 = preds_OBBA2_summary,
      OBBA3 = preds_OBBA3_summary,
      abs_change = preds_abs_change_summary,
      
      # For plotting boundaries of important/significant population change
      polygon_change_sf = polygon_change_sf,
      
      # summary of overall change across Ontario
      # overall_summary = overall_summary,
      
      # summary of hss and doy effects
      #hss_doy_summaries = hss_doy_summaries,
      
      meta = list(n_draws = ncol(eta))
    ), pred_path)
  }
  
  end <- Sys.time()
  message("Done: ", sp_english," - ",round(end-start,2)," min")
}

message("\n08_fit_models_and_predict_jointly.R complete.")
