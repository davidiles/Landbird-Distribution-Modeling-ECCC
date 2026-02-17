# ============================================================
# 08_fit_models_and_predict.R
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
dir.create(file.path(out_dir, "models_joint"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "predictions_joint"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "summaries_joint"), recursive = TRUE, showWarnings = FALSE)

# helper-function source (you should create this)
source("R/functions/inla_model_utils.R")  # fit_inla_multi_atlas(), predict_inla(), summarize_posterior(), etc.

# Which species to run
min_detections_obba3 <- 100
min_squares_obba3 <- 50
n_samples_predict <- 1000

# If model takes more than 15 min to fit, assume it stalled and retry
timeout_min <- 15

# Derived-covariate logic
south_bcr <- c(12, 13)
north_bcr <- c(7, 8)

# Covariates to *potentially* include; per-species covariate variance screening can happen later
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

# Save RDS
save_atomic <- function(obj, path) {
  tmp <- paste0(path, ".tmp")
  saveRDS(obj, tmp)
  file.rename(tmp, path)
}

# Create “derived” covariates consistently across surveys + grids
add_derived_covariates <- function(surveys, grid2, grid3,
                                   south_bcr = c(12, 13), north_bcr = c(7, 8)) {
  
  make_bits <- function(df) {
    if (!("water_river" %in% names(df))) df$water_river <- NA_real_
    if (!("BCR" %in% names(df))) stop("BCR missing; required for lc_* split.")
    
    df %>%
      mutate(
        on_river = as.integer(!is.na(water_river) & water_river > 0.1),
        on_road = as.integer(!is.na(road) & road > 0.1),
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

# Prediction formula builder: returns an inlabru generate() formula that outputs one column "eta"
make_pred_formula <- function(cov_df = NULL, include_kappa = FALSE, include_aru = FALSE) {
  cov_terms <- character(0)
  
  if (!is.null(cov_df) && nrow(cov_df) > 0) {
    cov_terms <- cov_df %>%
      dplyr::mutate(term = paste0("Beta", beta, "_", covariate, "*I(", covariate, "^", beta, ")")) %>%
      dplyr::pull(term)
  }
  
  base_terms <- c(
    "Intercept",
    "spde_abund",
    "Atlas3 * effect_Atlas3",
    "Atlas3 * spde_change",
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

# Add Atlas labels to grids (helps later)
grid_OBBA2 <- grid_OBBA2 %>% mutate(Atlas = "OBBA2")
grid_OBBA3 <- grid_OBBA3 %>% mutate(Atlas = "OBBA3")

# Derived covariates
tmp <- add_derived_covariates(all_surveys, grid_OBBA2, grid_OBBA3, south_bcr, north_bcr)
all_surveys <- tmp$surveys
grid_OBBA2  <- tmp$grid2
grid_OBBA3  <- tmp$grid3

# Summary caches
model_summaries_path  <- file.path(out_dir, "summaries_joint", "model_summaries.rds")
change_summaries_path <- file.path(out_dir, "summaries_joint", "change_summaries.rds")
hss_doy_path          <- file.path(out_dir, "summaries_joint", "HSS_DOY_summaries.rds")

model_summaries  <- load_or_empty_list(model_summaries_path)
change_summaries <- load_or_empty_list(change_summaries_path)
hss_doy_summaries <- load_or_empty_list(hss_doy_path)

# ------------------------------------------------------------
# Main loop (~ 30 min per species)
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
      Atlas3 = ifelse(Atlas == "OBBA2", 0, 1)
    )
  
  # --- Covariates spec
  covars_present <- intersect(base_covars, names(sp_dat))
  cov_df_sp <- make_cov_df(covars_present)
  
  # Note: previous defaults were stable as:
  # prior_range_abund  <- c(500, 0.90)
  # prior_sigma_abund  <- c(3,   0.05)
  # prior_range_change <- c(500, 0.10)
  # prior_sigma_change <- c(0.1, 0.05)
  
  # Priors controlling spatial autocorrelation fields
  prior_range_abund  <- c(100, 0.50) # 50% chance spatial autocorrelation is smaller than 100 km
  prior_sigma_abund  <- c(0.1, 0.50)
  prior_range_change <- c(250, 0.10) # 10 chance spatial autocorrelation is less than 250 km
  prior_sigma_change <- c(0.1, 0.05)
  
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
        family = "poisson"
      ),
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
        days_since_june15 = -7
      )
    
    pred_formula <- make_pred_formula(cov_df_sp)
    
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
    
    preds_obj <- list(
      sp_english = sp_english,
      sp_code = sp_code,
      n_draws = ncol(eta),
      pred_grid_meta = list(
        ref_HSS = 0,
        ref_days_since_june15 = -7
      ),
      OBBA2 = mu2,
      OBBA3 = mu3,
      abs_change = abs_change
    )
    
    
    # --- Summarize predictions for output tables
    preds_OBBA2_summary <- summarize_posterior(preds_obj$OBBA2, CI_probs = c(0.05, 0.95), prefix = "OBBA2")
    preds_OBBA3_summary <- summarize_posterior(preds_obj$OBBA3, CI_probs = c(0.05, 0.95), prefix = "OBBA3")
    preds_abs_change_summary <- summarize_posterior(preds_obj$abs_change, CI_probs = c(0.05, 0.95), prefix = "abs_change")
    
    # --- Overall + BCR change summaries (mask water if available)
    # If on_water exists, mask
    mu2_use <- preds_obj$OBBA2
    mu3_use <- preds_obj$OBBA3
    if ("on_water" %in% names(grid_OBBA2) && any(grid_OBBA2$on_water, na.rm = TRUE)) {
      mu2_use[grid_OBBA2$on_water, ] <- 0
    }
    if ("on_water" %in% names(grid_OBBA3) && any(grid_OBBA3$on_water, na.rm = TRUE)) {
      mu3_use[grid_OBBA3$on_water, ] <- 0
    }
    
    # Overall percent change across all pixels (using sums across pixels per draw)
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
    
    # --- HSS / DOY curves (lightweight)
    ref_doy <- -7
    ref_hss <- 0
    
    HSS_pred <- tibble(
      Hours_Since_Sunrise = seq(min(sp_dat$Hours_Since_Sunrise, na.rm = TRUE),
                                max(sp_dat$Hours_Since_Sunrise, na.rm = TRUE),
                                length.out = 100),
      days_since_june15 = ref_doy,
      effect = "HSS"
    )
    
    DOY_pred <- tibble(
      days_since_june15 = seq(min(sp_dat$days_since_june15, na.rm = TRUE),
                              max(sp_dat$days_since_june15, na.rm = TRUE),
                              by = 1),
      Hours_Since_Sunrise = ref_hss,
      effect = "DOY"
    )
    
    pred_df <- bind_rows(HSS_pred, DOY_pred)
    
    # predict only these components
    pred_all <- predict_inla(
      mod = mod,
      grid = pred_df,
      pred_formula = as.formula("~ data.frame(eta = HSS + DOY_global)"),
      n.samples = n_samples_predict,
      seed = 123
    )
    
    # summarize on response scale
    q <- t(apply(exp(pred_all$eta), 1, quantile, probs = c(0.025, 0.5, 0.975)))
    hss_doy_summaries <- bind_cols(pred_df, as_tibble(q, .name_repair = "minimal"))
    names(hss_doy_summaries)[(ncol(hss_doy_summaries)-2):ncol(hss_doy_summaries)] <- c("q025","q50","q975")
    
    # Save posterior summaries
    save_atomic(list(
      sp_english = sp_english,
      sp_code = sp_code,
      OBBA2 = preds_OBBA2_summary,
      OBBA3 = preds_OBBA3_summary,
      abs_change = preds_abs_change_summary,
      
      # summary of overall change across Ontario
      overall_summary = overall_summary,
      
      # summary of hss and doy effects
      hss_doy_summaries = hss_doy_summaries,
      
      meta = list(n_draws = preds_obj$n_draws,
                  ref_HSS = preds_obj$pred_grid_meta$ref_HSS,
                  ref_days_since_june15 = preds_obj$pred_grid_meta$ref_days_since_june15)
    ), pred_path)
  }
  
  end <- Sys.time()
  message("Done: ", sp_english," - ",round(end-start,2)," min")
}

message("\n08_fit_models_and_predict_jointly.R complete.")
