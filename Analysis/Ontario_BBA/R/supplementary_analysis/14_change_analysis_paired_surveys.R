# ============================================================
# 
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
model_name <- "joint_paired"
rerun_models <- TRUE
rerun_predictions <- TRUE

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
# RESTRICT COMPARISON TO SHARED FOOTPRINT
# ------------------------------------------------------------

all_surveys$row_id <- 1:nrow(all_surveys)

buffer_km <- 0.5

# Create buffer around atlas 2 points
a2_pts <- all_surveys %>% filter(Atlas == "OBBA2")
a3_pts <- all_surveys %>% filter(Atlas == "OBBA3")

# For each A2 point: which A3 points are within buffer?
a2_to_a3 <- st_is_within_distance(a2_pts, a3_pts, dist = buffer_km)

# Keep only A2 points that have at least one nearby A3
a2_keep <- lengths(a2_to_a3) > 0
a2_matched <- a2_pts[a2_keep, ]

# Now keep A3 points that are within buffer of at least one kept A2 point
a3_to_a2keep <- st_is_within_distance(a3_pts, a2_matched, dist = buffer_km)
a3_keep <- lengths(a3_to_a2keep) > 0

a3_matched <- a3_pts[a3_keep, ]

# Combined "shared footprint" dataset (matched A2 + nearby A3)
paired_surveys <- bind_rows(a2_matched, a3_matched) %>%
  mutate(shared_pairwise = TRUE) %>%
  arrange(row_id)

paired_counts <- counts[paired_surveys$row_id,] 

# ------------------------------------------------------------
# Main loop
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
  #"Palm Warbler",
  #"Bobolink",
  # "Blue Jay",
  #"Canada Jay",
  "Olive-sided Flycatcher"
  # "Winter Wren",
  #"Lesser Yellowlegs",
  #"Blackpoll Warbler",
  #"Connecticut Warbler",
  # "Lincoln's Sparrow",
  # "Fox Sparrow",
  #"Common Nighthawk"
  # "Long-eared Owl",
  # "American Tree Sparrow",
  # "LeConte's Sparrow",
  # "Nelson's Sparrow",
  # "Boreal Chickadee",
  # "Rusty Blackbird",
  # "Yellow-bellied Flycatcher",
  # "Greater Yellowlegs",
  # "Hudsonian Godwit",
  # "Canada Warbler",
  # "Eastern Wood-Peewee",
  # "Grasshopper Sparrow",
  # "Solitary Sandpiper",
  # # "White-throated Sparrow",
  # "Bay-breasted Warbler"
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
  
  model_path <- file.path(out_dir, paste0("models_",model_name), paste0(sp_file, ".rds"))
  pred_path  <- file.path(out_dir, paste0("predictions_",model_name), paste0(sp_file, ".rds"))
  
  # --- Build species data
  if (!(sp_code %in% names(counts))) {
    message("Skipping (species_id not found in counts columns): ", sp_code)
    next
  }
  
  # Append counts
  sp_dat_paired <- paired_surveys %>% mutate(count = paired_counts[[sp_code]])
  
  # ------------------------------------------------------
  # Mapping
  # ------------------------------------------------------
  
  ggplot() +
    # Plot all survey locations
    geom_sf(data = sp_dat_paired,
            col = "gray50",
            size = 2,
            shape = 18,
            stroke = 2.2) +
    
    # Overlay detection locations
    geom_sf(
      data = sp_dat_paired %>% filter(count > 0),
      shape = 4,
      stroke = 2.2,
      size = 2,
      aes(color = factor(ARU))
    ) +
    
    # Overlay BCR boundaries
    geom_sf(data = dat$bcr_sf, fill = NA, linewidth = 0.6) +
    
    scale_color_manual(values = c("black","blue"), name = "ARU")+
    coord_sf(expand = FALSE) +
    facet_grid(. ~ Atlas,
               labeller = labeller(
                 Atlas = c(
                   OBBA2 = "Atlas 2",
                   OBBA3 = "Atlas 3"
                 )
               )) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      strip.text = element_text(size = 12, face = "bold")
    )+
    ggtitle(sp_english)
  
  # ------------------------------------------------------
  # Consider evidence in "shared footprint" (near_A2) vs not,
  # Separately for each BCR
  # ------------------------------------------------------
  
  evidence_paired <- sp_dat_paired %>%
    st_drop_geometry() %>%
    group_by(Atlas, ARU, BCR) %>%
    summarise(
      n_surveys = n(),
      n_squares = n_distinct(square_id),
      frac_ARU  = mean(ARU == 1, na.rm = TRUE),
      n_svy_detected = sum(count>0,na.rm = TRUE),
      mean_svy_count = mean(count, na.rm = TRUE),
      prop_svy_detected = mean(count > 0, na.rm = TRUE),
      mean_duration_min = mean(Survey_Duration_Minutes, na.rm = TRUE),
      median_DOY = median(DayOfYear, na.rm = TRUE),
      q10_DOY = quantile(DayOfYear, 0.10, na.rm = TRUE),
      q90_DOY = quantile(DayOfYear, 0.90, na.rm = TRUE),
      median_HSS = median(Hours_Since_Sunrise, na.rm = TRUE),
      q10_HSS = quantile(Hours_Since_Sunrise, 0.10, na.rm = TRUE),
      q90_HSS = quantile(Hours_Since_Sunrise, 0.90, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(BCR,Atlas,ARU)
  
  evidence_paired %>% as.data.frame()
  
  ggplot(evidence_paired,aes(x = BCR,y = n_surveys,
                             fill = Atlas, 
                             linetype = factor(ARU),
                             col = factor(ARU)))+
    geom_bar(position = position_dodge(width = 1), stat = "identity", linewidth = 1)+
    scale_color_manual(values=c("transparent","black"),name = "ARU")+
    scale_linetype_manual(values=c(1,5),name = "ARU")+
    scale_fill_manual(values=c("black","dodgerblue"),name = "ARU")+
    theme_bw()+
    facet_grid(.~BCR, scales = "free_x")+
    scale_x_continuous(breaks = seq(1,100))+
    ggtitle("Day of Year")
  
  ggplot(evidence_paired,aes(x = BCR,y = median_DOY, ymin = q10_DOY, ymax = q90_DOY,col = Atlas, 
                                shape = factor(ARU),
                                linetype = factor(ARU)))+
    geom_point(position = position_dodge(width = 1), size = 3)+
    geom_errorbar(position = position_dodge(width = 1), linewidth = 1, width = 2)+
    scale_color_manual(values=c("black","dodgerblue"))+
    scale_shape_manual(values = c(19,17),name = "ARU")+
    scale_linetype_manual(values=c(1,5),name = "ARU")+
    theme_bw()+
    facet_grid(.~BCR, scales = "free_x")+
    scale_x_continuous(breaks = seq(1,100))+
    ggtitle("Day of Year")
  
  ggplot(evidence_paired,aes(x = BCR,y = median_HSS, ymin = q10_HSS, ymax = q90_HSS,col = Atlas, 
                                shape = factor(ARU),
         linetype = factor(ARU)))+
    geom_point(position = position_dodge(width = 1), size = 3)+
    geom_errorbar(position = position_dodge(width = 1), linewidth = 1, width = 2)+
    scale_color_manual(values=c("black","dodgerblue"))+
    scale_shape_manual(values = c(19,17),name = "ARU")+
    scale_linetype_manual(values=c(1,5),name = "ARU")+
    theme_bw()+
    facet_grid(.~BCR, scales = "free_x")+
    scale_x_continuous(breaks = seq(1,100))+
    ggtitle("Hours Since Sunrise")
  
  ggplot(evidence_paired,aes(x = BCR,y = mean_svy_count,col = Atlas, 
                                shape = factor(ARU),
                                linetype = factor(ARU)))+
    geom_point(position = position_dodge(width = 1), size = 3)+
    scale_color_manual(values=c("black","dodgerblue"))+
    scale_shape_manual(values = c(19,17),name = "ARU")+
    scale_linetype_manual(values=c(1,5),name = "ARU")+
    theme_bw()+
    facet_grid(.~BCR, scales = "free_x")+
    scale_x_continuous(breaks = seq(1,100))+
    ggtitle(paste0(sp_english,"\n\nMean count per 5-min survey"))
  
  
  
  # ------------------------------------------------------
  # Fit Bayesian model
  # ------------------------------------------------------
  
  # --- Covariates spec
  covars_present <- intersect(base_covars, names(sp_dat_paired))
  cov_df_sp <- make_cov_df(covars_present)
  
  # Priors controlling spatial autocorrelation fields
  prior_range_abund  <- c(250, 0.90) # 90% chance spatial autocorrelation range of shared field is less than 500 km
  prior_sigma_abund  <- c(3, 0.1)    # 10% chance SD is larger than 3
  prior_range_change <- c(250, 0.1)  # 10% chance spatial autocorrelation range of change field is less than 500 km
  prior_sigma_change <- c(0.1, 0.1)  # 10% chance SD is larger than 0.1
  
  # square identifier (persistent across atlases)
  sp_dat_paired$square_idx <- as.numeric(as.factor(sp_dat_paired$square_id))
  length(unique(sp_dat_paired$square_atlas))
  length(unique(sp_dat_paired$square_idx))
  
  # --- Fit (or load existing)
  mod <- NULL
  if (file.exists(model_path) & rerun_models == FALSE) {
    mod <- readRDS(model_path)
    message("Loaded existing model: ", model_path)
  } else {
    mod <- try(
      
      fit_inla_multi_atlas(
        sp_dat = sp_dat_paired,
        study_boundary = study_boundary,
        covariates = cov_df_sp,
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

message("\nfit_models_and_predict_jointly_paired.R complete.")
