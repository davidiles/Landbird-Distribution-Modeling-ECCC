# ============================================================
# 09_fit_models_and_predict_PC_ARU.R
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

model_name <- "PC_ARU"
rerun_models <- FALSE
rerun_predictions <- FALSE

in_file <- file.path(paths$data_clean, "birds", "data_ready_for_analysis.rds")
if (!file.exists(in_file)) {
  stop("Cannot find input at: ", in_file,
       "\nHave you run 07_filter_and_finalize_surveys.R?")
}

out_dir <- paths$model_output
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, paste0("models_",model_name)), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, paste0("predictions_",model_name)), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, paste0("summaries_",model_name)), recursive = TRUE, showWarnings = FALSE)

# Summary caches
model_summaries_path  <- file.path(out_dir, paste0("summaries_",model_name), "model_summaries.rds")
change_summaries_path <- file.path(out_dir, paste0("summaries_",model_name), "change_summaries.rds")
hss_doy_path          <- file.path(out_dir, paste0("summaries_",model_name), "HSS_DOY_summaries.rds")

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
safe_dates <- dat$safe_dates
species_to_model <- dat$species_to_model

stopifnot(nrow(all_surveys) == nrow(counts))

# Add Atlas labels to grids (helps later)
grid_OBBA2 <- grid_OBBA2 %>% mutate(Atlas = "OBBA2")
grid_OBBA3 <- grid_OBBA3 %>% mutate(Atlas = "OBBA3")

# Hexagon grid for summarizing statistical support for population change
hex_grid <- dat$hex_grid

# Associate each pixel with a particular hexagon
idx <- build_pixel_polygon_index(
  grid_sf      = grid_OBBA2,         # the full pixel grid
  polygons_sf  = hex_grid,           # e.g., BCRs/ecoregions/hexes
  poly_id_col  = "hex_id",
  join         = "within"
)

# ------------------------------------------------------------
# Assign safe dates to BCRs, rather than ecoregions
# ------------------------------------------------------------

safe_dates_bcr <- safe_dates %>%
  mutate(
    BCR = case_when(
      ecoregion == "Mixedwood Plains" ~ "13",
      ecoregion == "Hudson Plains"    ~ "7",
      ecoregion == "Boreal Shield"    ~ "8,12",
      TRUE ~ NA_character_
    )
  ) %>%
  separate_rows(BCR, sep = ",") %>%
  mutate(BCR = as.integer(BCR)) %>%
  select(sp_english, BCR, start_doy, end_doy) %>%
  arrange(sp_english, BCR) %>%
  mutate(midpoint = (start_doy + end_doy)/2)

# ------------------------------------------------------------
# Main loop
# ------------------------------------------------------------

# Covariates to *potentially* include; per-species covariate variance screening can happen later
base_covars <- c(
  "on_river",
  "on_road",
  "urban_3",
  "lc_1",
  "lc_4",
  "lc_5",
  "lc_8S","lc_8N",   # Woody savannas
  "lc_9S","lc_9N",   # Savannas
  "lc_10S","lc_10N", # Grasslands
  "lc_11","lc_12",
  "lc_14","lc_17",
  "insect_broadleaf","insect_needleleaf"
)

species_to_check <- c(
  "Palm Warbler",
  "Bobolink",
  "Blue Jay",
  "Canada Jay",
  "Olive-sided Flycatcher",
  "Winter Wren",
  "Lesser Yellowlegs",
  "Blackpoll Warbler",
  "Connecticut Warbler",
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
  "Bay-breasted Warbler",
  "Marsh Wren",
  "Eastern Towhee",
  "Northern Cardinal",
  "Carolina Wren",
  "Black-and-white Warbler",
  "Common Yellowthroat",
  "Gray Catbird",
  "Mourning Warbler",
  "Red-eyed Vireo",
  "Veery",
  "Tennessee Warbler",
  "Chestnut-sided Warbler",
  "Brown Thrasher",
  "American Redstart",
  "Hairy Woodpecker",
  "Belted Kingfisher",
  "Yellow-bellied Sapsucker",
  "Ruffed Grouse",
  "Red-winged Blackbird",
  "White-breasted Nuthatch",
  "Black-billed Cuckoo",
  "Northern Harrier",
  "European Starling",
  "House Finch",
  "Savannah Sparrow",
  "Alder Flycatcher",
  "Evening Grosbeak",
  "Northern Parula",
  "Pine Siskin",
  "Common Raven",
  "Sandhill Crane",
  "Wild Turkey",
  "Bald Eagle",
  "Blue-winged Teal",
  "Western Meadowlark",
  "Bank Swallow",
  "Tree Swallow",
  "Black-backed Woodpecker",
  "Black-capped Chickadee",
  "Black-throated Blue Warbler",
  "Black-throated Green Warbler",
  "Blackburnian Warbler",
  "Blue-headed Vireo",
  "Brown Creeper",
  "Chipping Sparrow",
  "Field Sparrow",
  "Golden-crowned Kinglet",
  "Golden-winged Warbler",
  "Least Flycatcher",
  "Nashville Warbler",
  "Northern Waterthrush",
  "Ovenbird",
  "Purple Finch",
  "Ruby-crowned Kinglet",
  "Swainson's Thrush",
  "Vesper Sparrow",
  "Warbling Vireo",
  "White-winged Crossbill",
  "Wilson's Warbler",
  "Orchard Oriole",
  "Turkey Vulture",
  "Sedge Wren",
  # "Herring Gull",
  "Black Tern",
  "Wood Thrush",
  "Purple Martin",
  "Brown-headed Cowbird",
  "Cliff Swallow",
  "Northern Rough-winged Swallow",
  # "Ring-billed Gull",
  "Chimney Swift",
  "Killdeer",
  "Great Blue Heron",
  "American Black Duck",
  "Yellow-throated Vireo",
  "House Wren",
  "American Robin"
) %>% unique()

# Data requirements
min_detections <- 150     # at least 150 detections in one atlas
min_squares <- 50         # at least 50 squares in one atlas

set.seed(123)
species_run <- species_to_model %>%
  filter((detections_safe_OBBA2 >= min_detections |
            detections_safe_OBBA3 >= min_detections )&
           (n_squares_safe_OBBA2 >= min_squares |
              n_squares_safe_OBBA3 >= min_squares))%>%
  na.omit() %>%
  mutate(delta_dets_safe = log(detections_safe_OBBA3/detections_safe_OBBA2),
         delta_squares_safe = log(n_squares_safe_OBBA3/n_squares_safe_OBBA2)) %>%
  
  # Arrange by species with largest changes
  arrange(desc(abs(delta_squares_safe)))

#species_run <- species_run %>% subset(english_name %in% species_to_check)

message("Species queued: ", nrow(species_run))
species_run

for (i in seq_len(nrow(species_run))) {
  
  # To choose a particular species:
  # i = which(species_run$english_name == "Purple Martin")
  
  sp_name <- species_run$english_name[i]
  sp_code <- as.character(species_run$species_id[i])
  sp_file <- sp_filename(sp_name)

  int_strategy = "auto" # can be faster with "eb", but more likley to fail
  strategy = "laplace"  # can be faster with "simplified.laplace" but more approximate
  
  message("\n====================\n", i, "/", nrow(species_run), ": ", sp_name, " (species_id = ", sp_code, ")\n====================")
  
  model_path <- file.path(out_dir,paste0("models_",model_name), paste0(sp_file, ".rds"))
  pred_path  <- file.path(out_dir, paste0("predictions_",model_name), paste0(sp_file, ".rds"))
  
  if (file.exists(pred_path) & rerun_predictions == FALSE){
    message("Predictions already exist for ", sp_name, "; skipping")
    next
  }
  
  # --- Build species data
  if (!(sp_code %in% names(counts))) {
    message("Skipping (species_id not found in counts columns): ", sp_code)
    next
  }
  
  # Append counts
  sp_dat <- all_surveys %>%
    dplyr::mutate(count = counts[[sp_code]])
  

  # ------------------------------------------------------------
  # Filter data to only surveys conducted within the safe dates for this species
  # ------------------------------------------------------------
  
  # Safe dates for this species
  sp_safe_dates <- safe_dates_bcr %>%
    dplyr::filter(sp_english == sp_name)
  
  if (nrow(sp_safe_dates) == 0) {
    pred_doy <- NA_real_
    
    warning(
      paste0("Species '", sp_name, "' has no BCR safe dates at all."),
      call. = FALSE
    )
    
  } else {
    
    # Narrowest overlap across known BCR safe-date windows
    fallback_start <- max(sp_safe_dates$start_doy, na.rm = TRUE)
    fallback_end   <- min(sp_safe_dates$end_doy,   na.rm = TRUE)
    
    if (fallback_start <= fallback_end) {
      
      all_bcrs <- sp_dat %>%
        dplyr::distinct(BCR)
      
      safe_bcrs <- sp_safe_dates %>%
        dplyr::distinct(BCR)
      
      missing_bcrs <- all_bcrs %>%
        dplyr::anti_join(safe_bcrs, by = "BCR")
      
      if (nrow(missing_bcrs) > 0) {
        sp_safe_dates <- dplyr::bind_rows(
          sp_safe_dates,
          missing_bcrs %>%
            dplyr::mutate(
              sp_english = sp_name,
              start_doy = fallback_start,
              end_doy   = fallback_end,
              midpoint  = floor((fallback_start + fallback_end) / 2)
            ) %>%
            dplyr::select(sp_english, BCR, start_doy, end_doy, midpoint)
        )
        
        message(
          paste0(
            "Species '", sp_name,
            "' was missing safe dates for BCR(s): ",
            paste(missing_bcrs$BCR, collapse = ", "),
            ". Assigned fallback dates ", fallback_start, " to ", fallback_end, "."
          )
        )
      }
      
      pred_doy <- floor((fallback_start + fallback_end) / 2)
      
      message(
        paste0(
          "Species '", sp_name,
          "': province-wide prediction DOY = ", pred_doy,
          " (overlap window ", fallback_start, " to ", fallback_end, ")."
        )
      )
      
    } else {
      pred_doy <- NA_real_
      
      message(
        paste0(
          "Species '", sp_name,
          "' has no single DOY that falls within all BCR safe-date windows. ",
          "Latest start_doy = ", fallback_start,
          ", earliest end_doy = ", fallback_end, "."
        )
      )
    }
  }
  
  # Filter to within-safe-date surveys
  sp_dat <- sp_dat %>%
    dplyr::left_join(sp_safe_dates, by = "BCR") %>%
    dplyr::filter(
      !is.na(start_doy),
      !is.na(end_doy),
      DayOfYear >= start_doy,
      DayOfYear <= end_doy
    ) %>%
    dplyr::mutate(
      days_midpoint = DayOfYear - pred_doy
    )
  
  # ------------------------------------------------------------
  # Use raw counts as a basis for deciding whether to fit poisson or nbinomial
  # ------------------------------------------------------------
  
  # ---- Remove extreme counts (assume these arise from a different observation process such as colony detection)
  large_counts <- sp_dat %>% filter(count > 50)
  sp_dat <- sp_dat %>% filter(count <= 50)
  
  # Number of detections
  n_det <- sum(sp_dat$count>0)
  y_pos <- sp_dat$count[sp_dat$count>0]
  max_count <- max(y_pos)
  
  # Proportion of total count in top 1% of counts
  n_top <- max(1, ceiling(0.01 * n_det))
  prop_total_top1pct_nonzero <- sum(sort(y_pos, decreasing = TRUE)[seq_len(n_top)]) / sum(y_pos)
  
  error_family <- "poisson"
  
  if (
    
    # More than 10% of total in top 1% of counts
    prop_total_top1pct_nonzero >= 0.10 |
    
    # More than 20 counts greater than 15
      sum(y_pos > 15) > 10){
    error_family = "nbinomial"
  }
  
  message("\n====================\n", i, "/", nrow(species_run), ": ", sp_name, " (species_id = ", sp_code, "): fitting model with ",error_family," error model\n====================")
  
  # ------------------------------------------------------------
  # Fit model
  # ------------------------------------------------------------
  
  # --- Covariates spec
  covars_present <- intersect(base_covars, names(sp_dat))
  cov_df_sp <- make_cov_df(covars_present)
  
  # --- Specify priors
  priors_list <- list(
    
    # SPDE priors
    prior_range_abund = c(250, 0.5),           # 10% chance spatial autocorrelation range of shared field is less than 200 km
    prior_sigma_abund = c(3, 0.1),             # 10% chance SD is larger than 3
    prior_range_change = c(250, 0.5),          # 10% chance spatial autocorrelation range of change field is less than 200 km
    prior_sigma_change = c(0.1, 0.1),          # 10% chance SD is larger than 0.1
    
    # Priors for HSS effect
    prior_HSS_range = c(5, 0.9),
    prior_HSS_sigma = c(3, 0.1),
    
    # Priors for DOY effects,
    prior_DOY_range_global = c(7, 0.9),
    prior_DOY_sigma_global = c(3, 0.1),        # global curve captures large deviations
    
    # atlas square iid effect
    kappa_pcprec_diff   = c(1, 0.1)           # 10% chance SD of iid random effect is larger than 1
    
  )
  
  # --- Fit (or load existing)
  start_model <- Sys.time()
  
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
        
        # SPDE priors
        prior_range_abund = priors_list$prior_range_abund,          
        prior_sigma_abund = priors_list$prior_sigma_abund,            
        prior_range_change = priors_list$prior_range_change,          
        prior_sigma_change = priors_list$prior_sigma_change,          
        
        # Priors for HSS effect
        prior_HSS_range = priors_list$prior_HSS_range,
        prior_HSS_sigma = priors_list$prior_HSS_sigma,
        
        # Priors for DOY effects,
        prior_DOY_range_global = priors_list$prior_DOY_range_global,
        prior_DOY_sigma_global = priors_list$prior_DOY_sigma_global,        
        
        #prior_DOY_range_BCR = priors_list$prior_DOY_range_BCR,
        #prior_DOY_sigma_BCR = priors_list$prior_DOY_sigma_BCR,         
        
        # atlas square iid effect
        kappa_pcprec_diff   = priors_list$kappa_pcprec_diff,
        
        int_strategy = int_strategy,
        strategy = strategy,
        
        family = error_family),
      
      silent = TRUE
    )
    
    if (inherits(mod, "try-error") || is.null(mod)) {
      message("Model failed for ", sp_name, "; continuing.")
      next
    }
    
    
    #save_atomic(mod, model_path)
    #message("Saved model: ", model_path)
  }
  end_model <- Sys.time()
  
  fit_minutes <- round(as.numeric(end_model - start_model, units = "mins"))
  
  message(
    "\n====================\n",
    i, "/", nrow(species_run), ": ", sp_name,
    " (species_id = ", sp_code, "): ",
    fit_minutes, " minutes to fit model\n",
    "===================="
  )
  
  # --- Save minimal model summary (fast)
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
  
  # --- Predictions (or load existing)
  preds_obj <- NULL
  if (file.exists(pred_path)  & rerun_predictions == FALSE) {
    preds_obj <- readRDS(pred_path)
    message("Loaded existing predictions: ", pred_path)
    
  } else {
    
    start_prediction <- Sys.time()
    message("Generating predictions")
    
    pred_grid <- bind_rows(
      grid_OBBA2 %>% mutate(Atlas3 = 0),
      grid_OBBA3 %>% mutate(Atlas3 = 1)
    ) %>%
      mutate(
        
        # reference values for predictions
        Hours_Since_Sunrise = 0,
        days_midpoint = 0,
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
    
    # Set "water" pixels to 0 expected abundance
    preds$eta[which(pred_grid$on_water),] <- -Inf  # Note this is log-scale predicted abundance, so log(0) = -Inf
    
    idx2 <- which(pred_grid$Atlas == "OBBA2")
    idx3 <- which(pred_grid$Atlas == "OBBA3")
    
    # preds$eta is n_grid x n_draw
    eta <- preds$eta
    mu2 <- exp(eta[idx2, , drop = FALSE])
    mu3 <- exp(eta[idx3, , drop = FALSE])
    
    # Mask out pixels that are considered "open water"
    
    # Calculate change
    abs_change <- mu3 - mu2
    
    # --- Summarize predictions for each pixel
    preds_OBBA2_summary <- summarize_posterior(mu2, CI_probs = c(0.05, 0.95), prefix = "OBBA2")
    preds_OBBA3_summary <- summarize_posterior(mu3, CI_probs = c(0.05, 0.95), prefix = "OBBA3")
    preds_abs_change_summary <- summarize_posterior(abs_change, CI_probs = c(0.05, 0.95), prefix = "abs_change")
    
    # -------------------------------------------------
    # Assess support for meaningful population change within larger polygons,
    # (e.g., for drawing "significance" boundaries on maps, or for summarizing results into BCRs)
    # -------------------------------------------------
    
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
    
    end_prediction <- Sys.time()
    
    pred_minutes <- round(as.numeric(end_prediction - start_prediction, units = "mins"))
    
    message(
      "\n====================\n",
      i, "/", nrow(species_run), ": ", sp_name,
      " (species_id = ", sp_code, "): ",
      pred_minutes, " minutes to generate predictions\n",
      "===================="
    )
    
    sp_square_summary <- sp_dat %>%
      as.data.frame() %>%
      group_by(Atlas,square_id) %>%
      summarize(n_surveys = n(),
                total_count = sum(count),
                n_detections = sum(count>0))
    
    sp_square_summary <- sp_dat %>%
      as.data.frame() %>%
      dplyr::group_by(Atlas, square_id) %>%
      dplyr::summarise(
        n_surveys    = dplyr::n(),
        total_count  = sum(count),
        n_detections = sum(count > 0),
        
        # BCR with the most surveys in this square
        BCR = names(which.max(table(BCR))),
        
        .groups = "drop"
      )
    
    # Save model summaries and predictions
    save_atomic(list(
      
      # Species data
      sp_name = sp_name,
      sp_code = sp_code,
      sp_safe_dates = sp_safe_dates,
      sp_square_summary = sp_square_summary,
      
      # Model priors and error structure
      error_family = error_family,
      priors = priors_list,
      
      # Fitted model summary
      fit_minutes = fit_minutes,
      summary_fixed = mod$summary.fixed,
      summary_hyperpar = mod$summary.hyperpar,
      
      # Predictions
      pred_minutes = pred_minutes,
      pred_doy = pred_doy,
      OBBA2 = preds_OBBA2_summary, # Per-pixel summaries
      OBBA3 = preds_OBBA3_summary, # Per-pixel summaries
      abs_change = preds_abs_change_summary, # Per-pixel summaries
      eta_draws_per_hex = eta_draws_per_hex # Posterior draws within each 25-km hexagon across Ontario
      
    ), pred_path)
  }
  
}

message("\n09_fit_models_and_predict_jointly.R complete.")
