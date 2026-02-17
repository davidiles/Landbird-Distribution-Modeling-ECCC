# ============================================================
# 08_fit_models_and_predict_separately.R
#
# Purpose:
#   Fit INLA/inlabru models for selected species, FITTING EACH
#   ATLAS SEPARATELY, predicting onto the corresponding atlas
#   grid (with time-specific covariates), and constructing a
#   change map by differencing posterior predictions pixel-by-pixel.
#
# Inputs:
#   data_clean/birds/data_ready_for_analysis.rds  (from script 07)
#
# Outputs (written incrementally per species):
#   data_clean/model_output/models_separate/<sp>__OBBA2.rds
#   data_clean/model_output/models_separate/<sp>__OBBA3.rds
#   data_clean/model_output/predictions_separate/<sp>.rds
#   data_clean/model_output/summaries_separate/model_summaries.rds
#   data_clean/model_output/summaries_separate/change_summaries.rds
#   data_clean/model_output/summaries_separate/HSS_DOY_summaries.rds
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

# ------------------------------------------------------------
# Config
# ------------------------------------------------------------
rerun_models = TRUE

in_file <- "data_clean/birds/data_ready_for_analysis.rds"

out_dir <- "data_clean/model_output"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "models_separate"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "predictions_separate"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "summaries_separate"), recursive = TRUE, showWarnings = FALSE)

# helper-function source
source("R/functions/inla_model_utils.R")

# Which species to run
min_detections_obba3 <- 100
min_squares_obba3 <- 50
n_samples_predict <- 1000

# If model takes more than 15 min to fit, assume it stalled
timeout_min <- 15

# Derived-covariate logic
south_bcr <- c(12, 13)
north_bcr <- c(7, 8)

# Covariates to *potentially* include
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
# Helper functions (script-local; move to utils if reused)
# ------------------------------------------------------------

sp_filename <- function(sp_english) {
  sp_english %>%
    str_to_lower() %>%
    str_replace_all("[^a-z0-9]+", "_") %>%
    str_replace_all("^_|_$", "")
}

load_or_empty_list <- function(path) {
  if (file.exists(path)) readRDS(path) else list()
}

# Save RDS atomically
save_atomic <- function(obj, path) {
  tmp <- paste0(path, ".tmp")
  saveRDS(obj, tmp)
  file.rename(tmp, path)
}

# Create “derived” covariates consistently across surveys + grids
add_derived_covariates <- function(surveys, grid2, grid3,
                                   south_bcr = c(12, 13),
                                   north_bcr = c(7, 8)) {
  
  make_bits <- function(df) {
    if (!("water_river" %in% names(df))) df$water_river <- NA_real_
    if (!("BCR" %in% names(df))) stop("BCR missing; required for lc_* split.")
    
    df %>%
      mutate(
        on_river = as.integer(!is.na(water_river) & water_river > 0.1),
        on_road  = as.integer(!is.na(road) & road > 0.1),
        lc_8S  = ifelse(BCR %in% south_bcr, lc_8, 0),
        lc_8N  = ifelse(BCR %in% north_bcr, lc_8, 0),
        lc_9S  = ifelse(BCR %in% south_bcr, lc_9, 0),
        lc_9N  = ifelse(BCR %in% north_bcr, lc_9, 0),
        lc_10S = ifelse(BCR %in% south_bcr, lc_10, 0),
        lc_10N = ifelse(BCR %in% north_bcr, lc_10, 0)
      )
  }
  
  surveys <- make_bits(surveys)
  grid2   <- make_bits(grid2)
  grid3   <- make_bits(grid3)
  
  list(surveys = surveys, grid2 = grid2, grid3 = grid3)
}

# Build a covariate specification dataframe (simple linear priors)
make_cov_df <- function(covars) {
  tibble(
    covariate = covars,
    beta = 1,
    sd_linear = 3,
    model = "linear",
    mean = 0,
    prec = 1 / (sd_linear^2)
  )
}

# Prediction formula builder for the SINGLE-ATLAS model
# (no Atlas3 terms, no spde_change)
make_pred_formula_single <- function(cov_df = NULL, include_kappa = FALSE, include_aru = FALSE) {
  cov_terms <- character(0)
  
  if (!is.null(cov_df) && nrow(cov_df) > 0) {
    cov_terms <- cov_df %>%
      dplyr::mutate(term = paste0("Beta", beta, "_", covariate, "*I(", covariate, "^", beta, ")")) %>%
      dplyr::pull(term)
  }
  
  base_terms <- c(
    "Intercept",
    "spde_abund",
    "DOY_global",
    "HSS"
  )
  
  if (include_aru) {
    base_terms <- c(base_terms, "ARU * effect_ARU")
  }
  
  eta_expr <- paste(c(base_terms, cov_terms), collapse = " + ")
  
  if (include_kappa) {
    eta_expr <- paste0(eta_expr, " + kappa")
  }
  
  as.formula(paste0("~ data.frame(eta = ", eta_expr, ")"))
}

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
stopifnot(nrow(grid_OBBA2) == nrow(grid_OBBA3))  # user guarantee of alignment

# Add Atlas labels to grids (optional convenience)
grid_OBBA2 <- grid_OBBA2 %>% mutate(Atlas = "OBBA2")
grid_OBBA3 <- grid_OBBA3 %>% mutate(Atlas = "OBBA3")

# Derived covariates
tmp <- add_derived_covariates(all_surveys, grid_OBBA2, grid_OBBA3, south_bcr, north_bcr)
all_surveys <- tmp$surveys
grid_OBBA2  <- tmp$grid2
grid_OBBA3  <- tmp$grid3

# Summary caches
model_summaries_path  <- file.path(out_dir, "summaries_separate", "model_summaries.rds")
change_summaries_path <- file.path(out_dir, "summaries_separate", "change_summaries.rds")
hss_doy_path          <- file.path(out_dir, "summaries_separate", "HSS_DOY_summaries.rds")

model_summaries   <- load_or_empty_list(model_summaries_path)
change_summaries  <- load_or_empty_list(change_summaries_path)
hss_doy_summaries <- load_or_empty_list(hss_doy_path)

# ------------------------------------------------------------
# Main loop
# ------------------------------------------------------------

# Select species list
species_run <- species_to_model %>%
  filter(total_detections_OBBA3 >= min_detections_obba3 &
           total_squares_OBBA3 >= min_squares_obba3) %>%
  na.omit()

message("Species queued: ", nrow(species_run))

# species_run <- sample_n(species_run, nrow(species_run))

species_to_check <- c(
  #"Bobolink",
  #"Blue Jay",
  "Canada Jay",
  "Olive-sided Flycatcher",
  "Winter Wren",
  "Lesser Yellowlegs",
  "Blackpoll Warbler",
  #"Connecticut Warbler",
  "Palm Warbler",
  "Lincoln's Sparrow",
  "Fox Sparrow",
  "Common Nighthawk",
  # "Long-eared Owl",
  "American Tree Sparrow",
  # "LeConte's Sparrow",
  # "Nelson's Sparrow",
  "Boreal Chickadee",
  "Rusty Blackbird",
  "Yellow-bellied Flycatcher",
  # "Greater Yellowlegs",
  # "Hudsonian Godwit",
  # "Canada Warbler",
  "Eastern Wood-Peewee",
  "Grasshopper Sparrow",
  "Solitary Sandpiper"
  #"White-throated Sparrow"
  # "Bay-breasted Warbler"
)

species_to_check <- "Canada Warbler"

species_run <- species_run %>%
  subset(english_name %in% species_to_check)

for (i in seq_len(nrow(species_run))) {
  
  start <- Sys.time()
  
  sp_english <- species_run$english_name[i]
  sp_code <- as.character(species_run$species_id[i])
  sp_file <- sp_filename(sp_english)
  
  message("\n====================\n",
          i, "/", nrow(species_run), ": ", sp_english, " (", sp_code, ")\n====================")
  
  # per-atlas model paths
  model_path2 <- file.path(out_dir, "models_separate", paste0(sp_file, "__OBBA2.rds"))
  model_path3 <- file.path(out_dir, "models_separate", paste0(sp_file, "__OBBA3.rds"))
  
  # combined predictions path (stores both surfaces + change)
  pred_path  <- file.path(out_dir, "predictions_separate", paste0(sp_file, ".rds"))
  
  # --- Build species data
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
  
  # Separate the data from each atlas, so they can be analyzed independently
  sp_dat2 <- sp_dat_all %>% filter(Atlas == "OBBA2")
  sp_dat3 <- sp_dat_all %>% filter(Atlas == "OBBA3")
  
  # --- Covariates spec (based on what exists in the survey table)
  covars_present <- intersect(base_covars, names(sp_dat_all))
  cov_df_sp <- make_cov_df(covars_present)
  
  # Priors controlling spatial autocorrelation field (same for both fits)
  # Priors controlling spatial autocorrelation fields
  prior_range_abund  <- c(100, 0.50) # 50% chance spatial autocorrelation is smaller than 100 km
  prior_sigma_abund  <- c(0.1, 0.05)

  # --- Fit/load OBBA2
  mod2 <- NULL
  if (file.exists(model_path2)  & rerun_models == FALSE) {
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
  
  # --- Fit/load OBBA3
  mod3 <- NULL
  if (file.exists(model_path3)  & rerun_models == FALSE) {
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
  
  # --- Predictions (or load existing)
  preds_obj <- NULL
  if (file.exists(pred_path) & rerun_models == FALSE) {
    preds_obj <- readRDS(pred_path)
    message("Loaded existing predictions: ", pred_path)
    
  } else {
    
    message("Generating predictions (separate fits for each atlas)")
    
    # Reference values for effect surfaces (same as joint script convention)
    ref_hss <- 0
    ref_doy <- -7
    
    pred_grid2 <- grid_OBBA2 %>%
      mutate(
        Hours_Since_Sunrise = ref_hss,
        days_since_june15 = ref_doy
      )
    
    pred_grid3 <- grid_OBBA3 %>%
      mutate(
        Hours_Since_Sunrise = ref_hss,
        days_since_june15 = ref_doy
      )
    
    # Single-atlas prediction formula
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
    
    # Pixel-by-pixel change (user guarantee: same row = same location)
    stopifnot(nrow(mu2) == nrow(mu3), ncol(mu2) == ncol(mu3))
    abs_change <- mu3 - mu2
    
    preds_obj <- list(
      sp_english = sp_english,
      sp_code = sp_code,
      
      # Keep both grids because their covariates differ by atlas
      #grid_OBBA2 = pred_grid2,
      #grid_OBBA3 = pred_grid3,
      
      n_draws = ncol(mu2),
      pred_grid_meta = list(ref_HSS = ref_hss, ref_days_since_june15 = ref_doy),
      
      OBBA2 = mu2,
      OBBA3 = mu3,
      abs_change = abs_change,
      
      change_method = "separate_fit_difference_pixelwise"
    )
    
    # save_atomic(preds_obj, pred_path)
    # message("Saved predictions: ", pred_path)
  }
  
  # --- Posterior summaries for mapping
  preds_OBBA2_summary <- summarize_posterior(preds_obj$OBBA2, prefix = "OBBA2")
  preds_OBBA3_summary <- summarize_posterior(preds_obj$OBBA3, prefix = "OBBA3")
  preds_abs_change_summary <- summarize_posterior(preds_obj$abs_change, prefix = "abs_change")
  
  # Attach summaries back onto grids
  preds_OBBA2_summary <- bind_cols(preds_obj$grid_OBBA2, as_tibble(preds_OBBA2_summary))
  preds_OBBA3_summary <- bind_cols(preds_obj$grid_OBBA3, as_tibble(preds_OBBA3_summary))
  
  # For change, use either grid geometry (they are identical). Use OBBA3 by default.
  preds_abs_change_summary <- bind_cols(preds_obj$grid_OBBA3, as_tibble(preds_abs_change_summary))
  
  # --- Overall percent change across all pixels (using sums across pixels per draw)
  # mimic joint script approach of summing pixels per draw
  mu2_use <- preds_obj$OBBA2
  mu3_use <- preds_obj$OBBA3
  
  if ("on_water" %in% names(preds_obj$grid_OBBA2) && any(preds_obj$grid_OBBA2$on_water, na.rm = TRUE)) {
    mu2_use[preds_obj$grid_OBBA2$on_water, ] <- 0
  }
  if ("on_water" %in% names(preds_obj$grid_OBBA3) && any(preds_obj$grid_OBBA3$on_water, na.rm = TRUE)) {
    mu3_use[preds_obj$grid_OBBA3$on_water, ] <- 0
  }
  
  sum2 <- colSums(mu2_use)
  sum3 <- colSums(mu3_use)
  pct_change <- (sum3 - sum2) / sum2 * 100
  
  overall_summary <- tibble(
    mean_change = mean(pct_change),
    median_change = median(pct_change),
    lower_change = unname(quantile(pct_change, 0.05)),
    upper_change = unname(quantile(pct_change, 0.95)),
    prob_decline = mean(pct_change <= 0),
    prob_decline_30 = mean(pct_change <= -30)
  )
  
  # --- HSS / DOY curves (lightweight) for both atlases
  # (optional but useful: you get two curves that can differ across periods)
  ref_doy <- preds_obj$pred_grid_meta$ref_days_since_june15
  ref_hss <- preds_obj$pred_grid_meta$ref_HSS
  
  make_hss_doy_preds <- function(sp_dat, which_atlas_label = c("OBBA2","OBBA3")) {
    which_atlas_label <- match.arg(which_atlas_label)
    
    HSS_pred <- tibble(
      Hours_Since_Sunrise = seq(min(sp_dat$Hours_Since_Sunrise, na.rm = TRUE),
                                max(sp_dat$Hours_Since_Sunrise, na.rm = TRUE),
                                length.out = 100),
      days_since_june15 = ref_doy,
      effect = "HSS",
      atlas = which_atlas_label
    )
    
    DOY_pred <- tibble(
      days_since_june15 = seq(min(sp_dat$days_since_june15, na.rm = TRUE),
                              max(sp_dat$days_since_june15, na.rm = TRUE),
                              by = 1),
      Hours_Since_Sunrise = ref_hss,
      effect = "DOY",
      atlas = which_atlas_label
    )
    
    bind_rows(HSS_pred, DOY_pred)
  }
  
  pred_df2 <- make_hss_doy_preds(sp_dat2, "OBBA2")
  pred_df3 <- make_hss_doy_preds(sp_dat3, "OBBA3")
  
  pred_all2 <- predict_inla(
    mod = mod2,
    grid = pred_df2,
    pred_formula = as.formula("~ data.frame(eta = HSS + DOY_global)"),
    n.samples = n_samples_predict,
    seed = 123
  )
  
  pred_all3 <- predict_inla(
    mod = mod3,
    grid = pred_df3,
    pred_formula = as.formula("~ data.frame(eta = HSS + DOY_global)"),
    n.samples = n_samples_predict,
    seed = 123
  )
  
  q2 <- t(apply(exp(pred_all2$eta), 1, quantile, probs = c(0.025, 0.5, 0.975)))
  q3 <- t(apply(exp(pred_all3$eta), 1, quantile, probs = c(0.025, 0.5, 0.975)))
  
  hss_doy_df2 <- bind_cols(pred_df2, as_tibble(q2, .name_repair = "minimal"))
  hss_doy_df3 <- bind_cols(pred_df3, as_tibble(q3, .name_repair = "minimal"))
  
  names(hss_doy_df2)[(ncol(hss_doy_df2)-2):ncol(hss_doy_df2)] <- c("q025","q50","q975")
  names(hss_doy_df3)[(ncol(hss_doy_df3)-2):ncol(hss_doy_df3)] <- c("q025","q50","q975")
  
  hss_doy_summaries[[sp_english]] <- list(
    sp_english = sp_english,
    sp_code = sp_code,
    curves = bind_rows(hss_doy_df2, hss_doy_df3)
  )
  save_atomic(hss_doy_summaries, hss_doy_path)
  
  # --- Save per-species summary bundle 
  save_atomic(
    list(
      sp_english = sp_english,
      sp_code = sp_code,
      
      OBBA2 = preds_OBBA2_summary,
      OBBA3 = preds_OBBA3_summary,
      abs_change = preds_abs_change_summary,
      
      # summary of overall change across Ontario
      overall_summary = overall_summary,
      
      # summary of hss and doy effects
      hss_doy_summaries = hss_doy_summaries,
      
      meta = list(
        n_draws = preds_obj$n_draws,
        ref_HSS = preds_obj$pred_grid_meta$ref_HSS,
        ref_days_since_june15 = preds_obj$pred_grid_meta$ref_days_since_june15
      )
    ),pred_path)
  
  end <- Sys.time()
  message("Done: ", sp_english, " - ", round(as.numeric(difftime(end, start, units = "mins")), 2), " min")
  
}

message("\n08_fit_models_and_predict_separately.R complete.")
