# ============================================================
# 10b_generate_crossvalidation_figures with cross-validation flags.R
# ============================================================

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(cowplot)
  library(pROC)
})

# ---- utils
source("R/functions/figure_utils.R")

# ------------------------------------------------------------
# Cross-validation overlay settings
# ------------------------------------------------------------
cv_dir_INLA <- "data_clean/model_output/xval_INLA"
cv_dir_BRT <- "data_clean/model_output/xval_BRT"

include_cv_flags_default <- TRUE

aggregate_cv_blocks <- function(block_summ) {
  req <- c("sp_english","sp_code","block_size_km","block_id","Atlas",
           "rep","fold","n","mean_obs","bias","rmse","mae","geometry",
           "mean_pred_yrep_q50","mean_pred_yrep_q025","mean_pred_yrep_q975","covered_95",
           "p_any_detect","obs_any_detect")
  stopifnot(all(req %in% names(block_summ)))
  
  block_summ %>%
    dplyr::group_by(sp_english, sp_code, block_size_km, block_id, Atlas) %>%
    dplyr::summarize(
      n_eval = dplyr::n_distinct(interaction(rep, fold, drop = TRUE)),
      n_surveys_total = mean(n, na.rm = TRUE),
      
      mean_obs = mean(mean_obs, na.rm = TRUE),
      
      mean_pred_yrep_q025 = mean(mean_pred_yrep_q025, na.rm = TRUE),
      mean_pred_yrep_q50  = mean(mean_pred_yrep_q50,  na.rm = TRUE),
      mean_pred_yrep_q975 = mean(mean_pred_yrep_q975, na.rm = TRUE),
      mean_pred_yrep_mean = mean(mean_pred_mean, na.rm = TRUE),
      
      covered_95 = mean(covered_95, na.rm = TRUE),
      
      rmse = mean(rmse, na.rm = TRUE),
      mae  = mean(mae,  na.rm = TRUE),
      bias = mean(bias, na.rm = TRUE),
      
      # presence/absence channel
      p_any_detect   = mean(p_any_detect, na.rm = TRUE),         # model probability of ≥1 detection
      obs_any_detect = as.integer(mean(obs_any_detect, na.rm = TRUE) >= 0.5),  # should be constant; this is defensive
      
      geometry = dplyr::first(geometry),
      .groups = "drop"
    ) %>%
    sf::st_as_sf()
}

# ------------------------------------------------------------
# Config
# ------------------------------------------------------------

in_data <- "data_clean/birds/data_ready_for_analysis.rds"
pred_dir <- "data_clean/model_output/predictions"
fig_dir  <- "data_clean/model_output/cv/figures/"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

in_bcr <- "../../Data/Spatial/BCR/BCR_Terrestrial_master.shp"
in_water <- "data_clean/spatial/water_filtered.shp"
in_atlas_squares <- "../../Data/Spatial/National/AtlasSquares/NationalSquares_FINAL.shp"

# ------------------------------------------------------------
# Load base data
# ------------------------------------------------------------
stopifnot(file.exists(in_data))
dat <- readRDS(in_data)

all_surveys <- dat$all_surveys
counts      <- dat$counts

grid2 <- dat$grid_OBBA2
grid3 <- dat$grid_OBBA3
study_boundary <- dat$study_boundary %>% sf::st_as_sf()

stopifnot(inherits(grid2, "sf"), inherits(grid3, "sf"))

crs_use <- st_crs(all_surveys)
study_boundary <- st_transform(study_boundary, crs_use)
grid3 <- st_transform(grid3, crs_use)
grid2 <- st_transform(grid2, crs_use)

bcr_sf <- st_read(in_bcr, quiet = TRUE) %>%
  st_make_valid() %>%
  st_transform(crs_use) %>%
  dplyr::filter(PROVINCE_S %in% c("ONTARIO", "ON", "Ontario") | is.na(PROVINCE_S)) %>%
  dplyr::select(BCR, BCRNAME, PROVINCE_S) %>%
  group_by(BCR, BCRNAME) %>%
  summarise(geometry = st_union(geometry), .groups = "drop")

# ------------------------------------------------------------
# Loop species
# ------------------------------------------------------------

fitted_species <- c("Olive-sided Flycatcher")

for (sp_english in fitted_species) {
  
  sp_file <- sp_english %>%
    str_to_lower() %>%
    str_replace_all("[^a-z0-9]+", "_") %>%
    str_replace_all("^_|_$", "")
  
  cv_path_survey_summary_INLA <- file.path(cv_dir_INLA, "predictions", paste0(sp_file, ".rds"))     # Predictions at individual survey level
  cv_path_block_summary_INLA <- file.path(cv_dir_INLA, "block_summaries", paste0(sp_file, ".rds"))  # Predictions summarized at block level
  if (!file.exists(cv_path_block_summary_INLA)) next
  
  cv_path_survey_summary_BRT <- file.path(cv_dir_BRT, "predictions", paste0(sp_file, ".rds"))     # Predictions at individual survey level
  cv_path_block_summary_BRT <- file.path(cv_dir_BRT, "block_summaries", paste0(sp_file, ".rds"))  # Predictions summarized at block level
  if (!file.exists(cv_path_block_summary_BRT)) next
  
  # ------------------------------------------------------------
  # Load block-level summaries and average across reps
  # ------------------------------------------------------------
  
  cv_block_summ_INLA <- readRDS(cv_path_block_summary_INLA)
  cv_block_avg_INLA  <- aggregate_cv_blocks(cv_block_summ_INLA)
  
  cv_block_summ_BRT <- readRDS(cv_path_block_summary_BRT)
  cv_block_avg_BRT  <- aggregate_cv_blocks(cv_block_summ_BRT)
  
  shared_blocks <- intersect(cv_block_avg_INLA$block_id,cv_block_avg_BRT$block_id)
  cv_block_avg_INLA <- cv_block_avg_INLA %>% subset(block_id %in% shared_blocks) %>% arrange(block_id)
  cv_block_avg_BRT <- cv_block_avg_BRT %>% subset(block_id %in% shared_blocks) %>% arrange(block_id)
  
  mean(cv_block_avg_INLA$mae)
  mean(cv_block_avg_BRT$mae)
  
  # ------------------------------------------------------------
  # Compare measures of cross-validation accuracy
  # ------------------------------------------------------------
  mean(cv_block_avg_INLA$mae)
  mean(cv_block_avg_BRT$mae)
  
  cor(cv_block_avg_INLA$mean_obs,cv_block_avg_INLA$mean_pred_yrep_q50)
  cor(cv_block_avg_INLA$mean_obs,cv_block_avg_BRT$mean_pred_yrep_q50)
  
}

message("09_make_species_figures.R complete: ", fig_dir)
