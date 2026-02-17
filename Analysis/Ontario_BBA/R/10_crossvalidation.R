# ============================================================
# 10_cross_validate_spatiotemporal_models.R
#
# Purpose:
#   Spatial block cross-validation for INLA/inlabru species models.
#   Refit models on training folds and predict to withheld blocks.
#
# Inputs:
#   data_clean/birds/data_ready_for_analysis.rds (from script 07)
#
# Outputs:
#   data_clean/model_output/cv/predictions/<sp>.rds  (row-level predictions)
#   data_clean/model_output/cv/summaries/<sp>.rds    (fold-level metrics)
#   data_clean/model_output/cv/meta/cv_plan.rds      (block->fold plan)
#   data_clean/model_output/cv/block_summaries/<sp>.rds (per-block metrics with square polygons)
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

in_file <- "data_clean/birds/data_ready_for_analysis.rds"

out_dir <- "data_clean/model_output/cv"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "predictions"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "summaries"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "meta"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "block_summaries"), recursive = TRUE, showWarnings = FALSE)

# Source helpers
source("R/functions/inla_model_utils.R")   # fit_inla_testing(), predict_inla(), summarize_posterior() 
source("R/functions/spatial_utils.R")      # make_spatial_blocks_grid(), make_block_cv_plan(), make_block_polygons_grid() 
source("R/functions/cv_utils.R")           # score_cv_predictions()

# Species selection (same spirit as script 08)
min_detections_obba3 <- 100
min_squares_obba3 <- 50

# CV design
block_size_km <- 30
n_folds <- 10
n_repeats <- 1
cv_seed <- 123

# Runtime controls
timeout_min <- 3
n_samples_predict_cv <- 1000  
retry_fit <- 1                # set to 1 to retry fitting if it fails; 0 to only try fitting a single time

# Note: previous defaults were stable as:
# prior_range_abund  <- c(500, 0.90)
# prior_sigma_abund  <- c(3,   0.05)
# prior_range_change <- c(500, 0.10)
# prior_sigma_change <- c(0.1, 0.05)

# Priors controlling spatial autocorrelation fields
prior_range_abund  <- c(100, 0.50) # 50% chance spatial autocorrelation is smaller than 100 km
prior_sigma_abund  <- c(3,   0.05) # 5% chance SD is larger than 3
prior_range_change <- c(250, 0.10) # 10 chance spatial autocorrelation is less than 250 km
prior_sigma_change <- c(0.1, 0.05) # 5% chance SD is larger than 0.1

# Derived-covariate logic (match script 08)
south_bcr <- c(12, 13)
north_bcr <- c(7, 8)

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

# Create “derived” covariates consistently across surveys + (optionally) grids
add_derived_covariates <- function(surveys,
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
  
  make_bits(surveys)
}

make_cov_df <- function(covars) {
  tibble(
    covariate = covars,
    beta = 1,
    sd_linear = 1,
    model = "linear",
    mean = 0,
    prec = 1 / (sd_linear^2)
  )
}

# Prediction formula builder (reuse exact logic from script 08)
make_pred_formula <- function(cov_df = NULL, include_kappa = FALSE, include_aru = FALSE) {
  cov_terms <- character(0)
  
  if (!is.null(cov_df) && nrow(cov_df) > 0) {
    cov_terms <- cov_df %>%
      mutate(term = paste0("Beta", beta, "_", covariate, "*I(", covariate, "^", beta, ")")) %>%
      pull(term)
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
study_boundary <- dat$study_boundary %>% st_as_sf()
species_to_model <- dat$species_to_model

stopifnot(nrow(all_surveys) == nrow(counts))

# Derived covariates (surveys only for CV)
all_surveys <- add_derived_covariates(all_surveys, south_bcr, north_bcr)

# Species list (match script 08 selection)
species_run <- species_to_model %>%
  filter(total_detections_OBBA3 >= min_detections_obba3 &
           total_squares_OBBA3 >= min_squares_obba3) %>%
  na.omit()

message("Species queued for CV: ", nrow(species_run))

# ------------------------------------------------------------
# Create CV blocks + plan (saved once)
# ------------------------------------------------------------

# Ensure we have geometry
stopifnot(inherits(all_surveys, "sf"))

# Create or load a persistent CV plan (do not overwrite once created)
cv_plan_path <- file.path(out_dir, "meta", paste0("cv_plan_blocks_", block_size_km, "km.rds"))

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
  
  start <- Sys.time()
  
  sp_english <- species_run$english_name[i]
  sp_code <- as.character(species_run$species_id[i])
  sp_file <- sp_filename(sp_english)
  
  message("\n====================\n", i, "/", nrow(species_run), ": ", sp_english, " (", sp_code, ")\n====================")
  
  if (!(sp_code %in% names(counts))) {
    message("Skipping (species_id not found in counts columns): ", sp_code)
    next
  }
  
  # Output paths (per species)
  pred_path  <- file.path(out_dir, "predictions", paste0(sp_file, ".rds"))
  summ_path  <- file.path(out_dir, "summaries",   paste0(sp_file, ".rds"))
  block_path <- file.path(out_dir, "block_summaries", paste0(sp_file, ".rds"))
  
  # Build species modeling data (exactly like script 08)
  sp_dat <- all_surveys %>%
    mutate(
      count = counts[[sp_code]],
      days_since_june15 = DayOfYear - 166,
      BCR_factor = as.numeric(factor(BCR)),
      Atlas3 = ifelse(Atlas == "OBBA2", 0, 1),
      block_id = block_ids
    )
  
  # Covariates spec
  covars_present <- intersect(base_covars, names(sp_dat))
  cov_df_sp <- make_cov_df(covars_present)
  
  # Prediction formula; include_kappa = TRUE generates new square-level random effects in predictions (turn off to remove this variance component from predictions)
  pred_formula <- make_pred_formula(cov_df_sp, include_kappa = TRUE, include_aru = TRUE)
  
  # Loop repeats
  for (r in seq_len(n_repeats)) {
    
    # Each repeat "reshuffles" the folds in which each block appears
    plan_r <- cv_plan %>% filter(rep == r)
    
    # Loop folds
    for (f in seq_len(n_folds)) {
      
      # determine withheld blocks for this fold
      test_blocks <- plan_r %>% filter(fold == f) %>% pull(block_id)
      
      is_test <- sp_dat$block_id %in% test_blocks
      dat_train <- sp_dat[!is_test, ]
      dat_test  <- sp_dat[ is_test, ]
      
      # Skip if this species/rep/fold has already been completed (i.e., all its withheld block_id x Atlas combinations are present)
      expected_pairs <- sf::st_drop_geometry(dat_test) %>%
        dplyr::distinct(block_id, Atlas)
      
      if (fold_already_completed(block_path, sp_code, block_size_km, r, f, expected_pairs)) {
        message("CV: repeat ", r, "/", n_repeats, " fold ", f, "/", n_folds,
                " already completed for ",
                sp_english, " — skipping.")
        next
      }
      
      # Guard against degenerate splits
      if (nrow(dat_test) == 0 || nrow(dat_train) == 0) {
        message("Repeat ", r, " fold ", f, ": empty train/test split; skipping.")
        next
      }
      
      message("CV: repeat ", r, "/", n_repeats, " fold ", f, "/", n_folds,
              " | train n=", nrow(dat_train), " test n=", nrow(dat_test))
      
      start <- Sys.time()
      # Attempt to fit the model
      fit_res <- fit_with_fallback(
        dat_train = dat_train,
        study_boundary = study_boundary,
        cov_df_sp = cov_df_sp,
        timeout_min = timeout_min,
        prior_range_abund = prior_range_abund,
        prior_sigma_abund = prior_sigma_abund,
        prior_range_change = prior_range_change,
        prior_sigma_change = prior_sigma_change,
        retry_fit = retry_fit,
        nb_pc_target_prob = 0.5,
        nb_pc_threshold_theta = 5,
        family_order = "poisson",
        max_tries_per_family = 1,
        verbose = TRUE
      )
      end <- Sys.time()
      mod <- fit_res$mod
      family_used <- fit_res$family_used
      
      if (is.null(mod)) {
        message("Fit failed after NB(2) + Poisson(2) attempts (rep=", r, ", fold=", f, "). Recording NA summary + continuing.")
        
        summ_row <- tibble(
          sp_english = sp_english,
          sp_code = sp_code,
          block_size_km = block_size_km,
          rep = r,
          fold = f,
          fit_ok = FALSE,
          family_used = family_used,   # will be NA here
          n_train = nrow(dat_train),
          n_test  = nrow(dat_test),
          n_test_blocks = length(unique(test_blocks)),
          n_draws = NA_integer_,
          rmse = NA_real_,
          mae  = NA_real_,
          spearman = NA_real_,
          pearson  = NA_real_,
          poisson_dev = NA_real_,
          mean_obs = NA_real_,
          mean_pred_mean = NA_real_,
          mean_pred_q025 = NA_real_,
          mean_pred_q50  = NA_real_,
          mean_pred_q975 = NA_real_
        )
        
        append_rds_rows(
          summ_path,
          summ_row,
          key_cols = c("sp_code","rep","fold","block_size_km")
        )
        next
      }
      
      # Predict to withheld surveys
      pred_draws <- predict_inla(
        mod = mod,
        grid = dat_test,
        pred_formula = pred_formula,
        n.samples = n_samples_predict_cv,
        seed = 1000 + 100*r + f
      )
      
      # posterior draws: eta is n_test x n_draws
      eta <- pred_draws$eta
      mu_draws <- exp(eta)  
      mu_mean  <- rowMeans(mu_draws)     # posterior mean per survey
      
      # ------------------------------------------------------------
      # Posterior predictive draws for counts (Poisson)
      # These represent *new* observations, conditional on posterior draws of mu.
      # ------------------------------------------------------------
      
      yrep_draws <- sim_yrep_draws(mu_draws, dat_test, mod)
      
      # Per-survey predictive summaries from posterior predictive draws
      yrep_q025 <- apply(yrep_draws, 1, stats::quantile, probs = 0.025, na.rm = TRUE)
      yrep_q50  <- apply(yrep_draws, 1, stats::quantile, probs = 0.50,  na.rm = TRUE)
      yrep_q975 <- apply(yrep_draws, 1, stats::quantile, probs = 0.975, na.rm = TRUE)
      
      # Event probabilities (example: at least one detection)
      p_ge1 <- rowMeans(yrep_draws >= 1, na.rm = TRUE)
      
      # ------------------------------------------------------------
      # Save survey-level predictions, retaining geometry (sf)
      # ------------------------------------------------------------
      
      pred_rows <- dat_test %>%
        mutate(
          sp_english = sp_english,
          sp_code = sp_code,
          block_size_km = block_size_km,
          rep = r,
          fold = f,
          block_id = block_id,
          
          # optional id fields if present
          survey_id = if ("survey_id" %in% names(dat_test)) dat_test$survey_id else NA_character_,
          Atlas = Atlas,
          count_obs = count,
          mu_pred = mu_mean,
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
      
      # Save predictions incrementally
      append_rds_rows(
        pred_path,
        pred_rows,
        key_cols = c("sp_code","rep","fold","block_id","survey_id","Atlas")
      )
      
      # ------------------------------------------------------------
      # Fold-level scoring (separately for each atlas)
      # ------------------------------------------------------------
      
      test_df <- sf::st_drop_geometry(pred_rows)
      
      atlases_in_test <- sort(unique(test_df$Atlas))
      
      for (atl in atlases_in_test) {
        
        test_df_atl <- test_df %>% dplyr::filter(Atlas == atl)
        if (nrow(test_df_atl) == 0) next
        
        # Point metrics based on posterior mean of mu
        fold_score <- score_cv_predictions(test_df_atl, y_col = "count_obs", mu_col = "mu_pred")
        
        # Probabilistic metrics
        idx_atl <- which(test_df$Atlas == atl)
        
        # Proportion of observations within 95% prediction interval within the fold
        cov95 <- interval_coverage(
          y = test_df$count_obs[idx_atl],
          lo = test_df$yrep_q025[idx_atl],
          hi = test_df$yrep_q975[idx_atl]
        )
        
        mean_log_score <- poisson_log_score_mc(
          y = test_df$count_obs[idx_atl],
          mu_draws = mu_draws[idx_atl, , drop = FALSE]
        )
        
        brier_ge1 <- brier_score(
          y01 = as.numeric(test_df$count_obs[idx_atl] >= 1),
          p = test_df$p_ge1[idx_atl]
        )
        
        # Fold mean observed count
        fold_mean_obs <- mean(test_df$count_obs[idx_atl], na.rm = TRUE)
        
        # Fold mean predicted mu (process mean), per draw
        fold_mu_draw_means <- colMeans(mu_draws[idx_atl, , drop = FALSE], na.rm = TRUE)
        
        fold_mean_pred_mean <- mean(fold_mu_draw_means, na.rm = TRUE)
        fold_mean_pred_q025 <- unname(stats::quantile(fold_mu_draw_means, 0.025, na.rm = TRUE))
        fold_mean_pred_q50  <- unname(stats::quantile(fold_mu_draw_means, 0.50,  na.rm = TRUE))
        fold_mean_pred_q975 <- unname(stats::quantile(fold_mu_draw_means, 0.975, na.rm = TRUE))
        
        # Fold mean predicted counts (posterior predictive), per draw
        fold_yrep_draw_means <- colMeans(yrep_draws[idx_atl, , drop = FALSE], na.rm = TRUE)
        
        fold_mean_pred_yrep_mean <- mean(fold_yrep_draw_means, na.rm = TRUE)
        fold_mean_pred_yrep_q025 <- unname(stats::quantile(fold_yrep_draw_means, 0.025, na.rm = TRUE))
        fold_mean_pred_yrep_q50  <- unname(stats::quantile(fold_yrep_draw_means, 0.50,  na.rm = TRUE))
        fold_mean_pred_yrep_q975 <- unname(stats::quantile(fold_yrep_draw_means, 0.975, na.rm = TRUE))
        
        summ_row <- fold_score %>%
          mutate(
            sp_english = sp_english,
            sp_code = sp_code,
            block_size_km = block_size_km,
            rep = r,
            fold = f,
            Atlas = atl,
            mean_obs = fold_mean_obs,
            mean_pred_mean = fold_mean_pred_mean,
            mean_pred_q025 = fold_mean_pred_q025,
            mean_pred_q50  = fold_mean_pred_q50,
            mean_pred_q975 = fold_mean_pred_q975,
            mean_pred_yrep_mean = fold_mean_pred_yrep_mean,
            mean_pred_yrep_q025 = fold_mean_pred_yrep_q025,
            mean_pred_yrep_q50  = fold_mean_pred_yrep_q50,
            mean_pred_yrep_q975 = fold_mean_pred_yrep_q975,
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
      
      # ------------------------------------------------------------
      # Per-block summaries (with polygons), separately for each atlas
      # ------------------------------------------------------------
      
      # Per-block summaries for this fold (use actual square polygon geometry per block)
      # We compute mean predicted count per block WITH 95% credible intervals by:
      #   - taking mu_draws = exp(eta) for each survey in the test set
      #   - for each posterior draw: average mu across surveys in the block
      #   - summarize those draw-level block means
      
      test_df <- sf::st_drop_geometry(dat_test) %>%
        transmute(block_id = block_id,
                  Atlas = as.character(Atlas),
                  count_obs = count)
      
      atlases_in_test <- sort(unique(test_df$Atlas))
      
      for (atl in atlases_in_test) {
        
        test_df_atl <- test_df %>% filter(Atlas == atl)
        if (nrow(test_df_atl) == 0) next
        
        block_ids_in_test <- unique(test_df_atl$block_id)
        
        block_summ_list <- lapply(block_ids_in_test, function(bid) {
          idx <- which(test_df$block_id == bid & test_df$Atlas == atl)
          
          # Observed presence/absence at the block level
          obs_all_zero  <- as.integer(all(test_df$count_obs[idx] == 0, na.rm = TRUE))
          obs_any_detect <- 1L - obs_all_zero
          
          # Posterior predictive probability that ALL surveys in the block are zero
          # yrep_draws is n_test x n_draws
          # For each draw (column), are all surveys in this block zero?
          p_all_zero <- mean(colSums(yrep_draws[idx, , drop = FALSE] > 0, na.rm = TRUE) == 0)
          
          # Probability of at least one detection in the block
          p_any_detect <- 1 - p_all_zero
          
          # observed mean in this block (for this atlas)
          mean_obs_block <- mean(test_df$count_obs[idx], na.rm = TRUE)
          
          # predicted mean per draw in this block (for this atlas)
          block_mu_draw_means <- colMeans(mu_draws[idx, , drop = FALSE], na.rm = TRUE)
          
          # predicted count mean per draw in this block (posterior predictive)
          block_yrep_draw_means <- colMeans(yrep_draws[idx, , drop = FALSE], na.rm = TRUE)
          
          tibble(
            sp_english = sp_english,
            sp_code = sp_code,
            block_size_km = block_size_km,
            rep = r,
            fold = f,
            Atlas = atl,
            block_id = bid,
            
            n = length(idx),
            
            # per-block mean observed and predicted with 95% CrI
            mean_obs = mean_obs_block,
            mean_pred_mean = mean(block_mu_draw_means, na.rm = TRUE),
            mean_pred_q025 = unname(stats::quantile(block_mu_draw_means, 0.025, na.rm = TRUE)),
            mean_pred_q50  = unname(stats::quantile(block_mu_draw_means, 0.50,  na.rm = TRUE)),
            
            # Posterior predictive (count) distribution for the block mean
            mean_pred_yrep_mean = mean(block_yrep_draw_means, na.rm = TRUE),
            mean_pred_yrep_q025 = unname(stats::quantile(block_yrep_draw_means, 0.025, na.rm = TRUE)),
            mean_pred_yrep_q50  = unname(stats::quantile(block_yrep_draw_means, 0.50,  na.rm = TRUE)),
            mean_pred_yrep_q975 = unname(stats::quantile(block_yrep_draw_means, 0.975, na.rm = TRUE)),
            covered_95 = as.integer(mean_obs_block >= unname(stats::quantile(block_yrep_draw_means, 0.025, na.rm = TRUE)) &
                                      mean_obs_block <= unname(stats::quantile(block_yrep_draw_means, 0.975, na.rm = TRUE))),
            
            # For AUC assessments
            p_all_zero   = p_all_zero,
            p_any_detect = p_any_detect,
            obs_any_detect = obs_any_detect,
            
            # Optional (useful for “where is it worst?”)
            rmse = sqrt(mean((mu_mean[idx] - test_df$count_obs[idx])^2, na.rm = TRUE)),
            mae  = mean(abs(mu_mean[idx] - test_df$count_obs[idx]), na.rm = TRUE),
            bias = mean(mu_mean[idx] - test_df$count_obs[idx], na.rm = TRUE)
          )
        }) # close lapply
        
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
      
    } # close folds loop
    
  } # close repeats loop
  
  end <- Sys.time()
  message("Done CV: ", sp_english, " - ", round(end - start, 2), " min")
}

message("\n10_cross_validate_spatiotemporal_models.R complete.")
