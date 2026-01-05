# ===============================================================================
# BAYESIAN ANALYSIS / SPECIES DISTRIBUTION MODELS FOR ONTARIO BREEDING BIRD ATLAS
# ===============================================================================

library(tidyverse)
library(sf)
library(fmesher)
library(INLA)
library(inlabru)
library(readxl)
library(patchwork)
library(DHARMa)
library(pROC)

rm(list=ls())

setwd("C:/Users/IlesD/OneDrive - EC-EC/Iles/Projects/Landbirds/Landbird-Distribution-Modeling-ECCC/")

# Load prepared data
dat <- readRDS(file = "Analysis/Ontario_BBA/Data_Cleaned/data_ready_for_analysis.rds")
all_surveys <- dat$all_surveys
counts <- dat$counts
grid_OBBA2 <- dat$grid_OBBA2 %>% mutate(Atlas = "OBBA2") %>% na.omit()
grid_OBBA3 <- dat$grid_OBBA3 %>% mutate(Atlas = "OBBA3") %>% na.omit()
species_square_summary <- dat$species_square_summary
species_to_model <- dat$species_to_model
study_boundary <- dat$study_boundary
water = dat$water

# ---------------------------------------------
# ---- Load several shapefiles needed for plotting
# ---------------------------------------------

# Load shapefile containing atlas squares (10km x 10km)
atlas_squares <- st_read("Data/Spatial/National/AtlasSquares/NationalSquares_FINAL.shp")  %>%
  subset(prov == "ON") %>%
  st_transform(st_crs(all_surveys))

# BCR boundaries
allBCR <- st_read("Data/Spatial/BCR/BCR_Terrestrial_master.shp") %>%
  st_transform(st_crs(study_boundary)) %>%
  st_make_valid() %>%
  st_buffer(0) %>%
  group_by(BCR,BCRNAME,PROVINCE_S) %>%
  summarise() %>%
  st_transform(st_crs(all_surveys))

BCR <- st_read("Data/Spatial/BCR/BCR_Terrestrial_master.shp") %>%
  subset(PROVINCE_S == "ONTARIO") %>%
  st_transform(st_crs(study_boundary)) %>%
  st_make_valid() %>%
  st_buffer(0) %>%
  group_by(BCR,BCRNAME) %>%
  summarise() %>%
  st_transform(st_crs(all_surveys))

# ---------------------------------------------
# ---- Re-classify pixels into either "near river" or "not near river"
# ---------------------------------------------

grid_OBBA2$near_river <- as.numeric(grid_OBBA2$water_river > 0)
grid_OBBA3$near_river <- as.numeric(grid_OBBA3$water_river > 0)
all_surveys$near_river <- as.numeric(all_surveys$water_river>0)

# ---------------------------------------------
# ---- Re-classify woody savanna (lc_8), savannah (lc_9), and grassland (lc_10) into southern and northern components 
# ---------------------------------------------

grid_OBBA2$lc_8S <- ifelse(grid_OBBA2$BCR %in% c(13,12),grid_OBBA2$lc_8,0)
grid_OBBA2$lc_8N <- ifelse(grid_OBBA2$BCR %in% c(7,8),grid_OBBA2$lc_8,0)
grid_OBBA3$lc_8S <- ifelse(grid_OBBA3$BCR %in% c(13,12),grid_OBBA3$lc_8,0)
grid_OBBA3$lc_8N <- ifelse(grid_OBBA3$BCR %in% c(7,8),grid_OBBA3$lc_8,0)
all_surveys$lc_8S <- ifelse(all_surveys$BCR %in% c(13,12),all_surveys$lc_8,0)
all_surveys$lc_8N <- ifelse(all_surveys$BCR %in% c(7,8),all_surveys$lc_8,0)

grid_OBBA2$lc_9S <- ifelse(grid_OBBA2$BCR %in% c(13,12),grid_OBBA2$lc_9,0)
grid_OBBA2$lc_9N <- ifelse(grid_OBBA2$BCR %in% c(7,8),grid_OBBA2$lc_9,0)
grid_OBBA3$lc_9S <- ifelse(grid_OBBA3$BCR %in% c(13,12),grid_OBBA3$lc_9,0)
grid_OBBA3$lc_9N <- ifelse(grid_OBBA3$BCR %in% c(7,8),grid_OBBA3$lc_9,0)
all_surveys$lc_9S <- ifelse(all_surveys$BCR %in% c(13,12),all_surveys$lc_9,0)
all_surveys$lc_9N <- ifelse(all_surveys$BCR %in% c(7,8),all_surveys$lc_9,0)

grid_OBBA2$lc_10S <- ifelse(grid_OBBA2$BCR %in% c(13,12),grid_OBBA2$lc_10,0)
grid_OBBA2$lc_10N <- ifelse(grid_OBBA2$BCR %in% c(7,8),grid_OBBA2$lc_10,0)
grid_OBBA3$lc_10S <- ifelse(grid_OBBA3$BCR %in% c(13,12),grid_OBBA3$lc_10,0)
grid_OBBA3$lc_10N <- ifelse(grid_OBBA3$BCR %in% c(7,8),grid_OBBA3$lc_10,0)
all_surveys$lc_10S <- ifelse(all_surveys$BCR %in% c(13,12),all_surveys$lc_10,0)
all_surveys$lc_10N <- ifelse(all_surveys$BCR %in% c(7,8),all_surveys$lc_10,0)

# ---------------------------------------------
# ---- Variance in covariates within each BCR
# ---------------------------------------------

cov_dat <- bind_rows(grid_OBBA2, grid_OBBA3) %>%
  as.data.frame() %>%
  dplyr::select(-geometry,-urban_1,-urban_2)

covar_names <- cov_dat %>%
  select(where(is.numeric), -BCR) %>%
  names()

cov_summary <- cov_dat %>%
  group_by(BCR) %>%
  summarize(
    # SD
    across(
      all_of(covar_names),
      sd,
      na.rm = TRUE,
      .names = "sd__{.col}"
    ),
    
    # Median
    across(
      all_of(covar_names),
      median,
      na.rm = TRUE,
      .names = "med__{.col}"
    ),
    
    # Proportion > 0.1 from median
    across(
      all_of(covar_names),
      ~ mean(abs(.x - median(.x, na.rm = TRUE)) > 0.1, na.rm = TRUE),
      .names = "gt_0.1__{.col}"
    ),
    
    .groups = "drop"
  )

cov_summary_long <- cov_summary %>%
  pivot_longer(
    -BCR,
    names_to = c("stat", "covariate"),
    names_sep = "__",
    values_to = "value"
  )

# ---------------------------------------------
# ---- Covariates to include in model
# ---------------------------------------------

covars <- c(#"insect_broadleaf",
  # "insect_needleleaf",
  # "road",
  "near_river",
  #"prec",
  # "tmax",
  #"urban_2",
  "urban_3",
  "lc_1",  # Evergreen Needleaf Forest
  "lc_4",  # Deciduous Broadleaf Forest
  "lc_5",  # Mixed Forest
  
  
  "lc_8S",  # Woody Savanna; Tree cover 30-60% (canopy >2m
  "lc_8N",  # Woody Savanna; Tree cover 30-60% (canopy >2m
  
  "lc_9S",  # Savanna; Tree cover 10-30% (canopy >2m)
  "lc_9N",  # Savanna; Tree cover 10-30% (canopy >2m)
  
  "lc_10S", # Grassland
  "lc_10N", # Grassland
  
  "lc_11", # Permanent Wetland
  "lc_12", # Cropland
  "lc_14", # Cropland / Natural Vegetation Mosaic
  "lc_17"  # Water
)
#covars <- covars[1:2]

# Priors for covariate effects

cov_df <- data.frame(covariate = covars, 
                     beta = 1,
                     sd_linear = 1,
                     model = "linear",
                     mean = 0) %>%
  mutate(prec = ifelse(covariate %in% c("prec","tmax"),1/1^2,1/sd_linear^2))

#cov_df <- cov_df %>% bind_rows(cov_df %>% mutate(beta = 2, prec = 1 / (sd_linear / 2)^2))

# ---------------------------------------------
# ---- Loop through species and conduct analysis
# ---------------------------------------------

species_to_model <- dat$species_to_model %>%
  subset(english_name %in% c("Bobolink",
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
                             "White-throated Sparrow"
  )) %>%
  subset(total_detections_OBBA3 > 100)


for (i in c(17,18,19,20)){
  
  model_summaries <- list()
  if (file.exists("Analysis/Ontario_BBA/Output/Tables_Summaries/model_summaries.rds")) model_summaries <- readRDS("Analysis/Ontario_BBA/Output/Tables_Summaries/model_summaries.rds")
  
  # ---------------------------------------------
  # ---- Prepare data for this species
  # ---------------------------------------------
  
  sp_english <- species_to_model$english_name[i]
  sp_filename <- sp_english %>%
    str_to_lower() %>%                       # lower case
    str_replace_all("[^a-z0-9]+", "_") %>%   # replace non-alphanumeric with _
    str_replace_all("^_|_$", "")             # trim leading/trailing underscores
  
  sp_code <- subset(species_to_model, english_name == sp_english)$species_id[1] %>% as.character()
  
  print(sp_english)
  
  # Prepare species survey data
  sp_dat <- all_surveys %>% 
    mutate(count = counts[[sp_code]],
           days_since_june15 = DayOfYear - 166,
           BCR_factor = as.numeric(factor(BCR))) 
  
  # Optional: skip species if already run
  # if (length(model_summaries )>0 & sp_english %in% names(model_summaries )) next
  
  # ---------------------------------------------
  # ---- Define area for modeling
  # ---------------------------------------------
  
  # Identify BCRs in which species has been observed
  BCRs_to_model <- sp_dat %>%
    as.data.frame() %>%
    group_by(BCR,square_atlas) %>%
    summarize(det = sum(count>0)) %>%
    group_by(BCR) %>%
    summarize(n_sq = sum(det>0)) %>%
    mutate(BCR_factor = as.numeric(factor(BCR))) %>%
    subset(n_sq > 10)
  
  model_boundary_sp <- BCR %>% 
    st_union() %>%
    st_set_precision(1) %>%
    st_make_valid() %>%
    st_collection_extract("POLYGON") %>%
    st_cast("MULTIPOLYGON") %>%
    st_as_sf()
  
  # Restrict covariate list to those with useful variance within the BCRs in which species was observed
  valid_covariates <- cov_summary_long %>%
    filter(stat == "gt_0.1" & value >= 0.2) %>%
    semi_join(BCRs_to_model, by = "BCR") %>%
    pull(covariate)
  
  valid_covariates <- c("near_river",valid_covariates) %>% unique()
  
  cov_df_sp <- cov_df %>%
    filter(covariate %in% valid_covariates)
  
  dim(cov_df_sp)
  cov_df_sp
  
  # ---------------------------------------------
  # ---- Fit models
  # ---------------------------------------------
  
  sp_dat <- subset(sp_dat, BCR %in% BCRs_to_model$BCR)
  
  # n_sq <- min(c(species_to_model$total_squares_OBBA2[i],species_to_model$total_squares_OBBA3[i]))
  # 
  # if (n_sq < 100) {
  #   prior_range_abund  <- c(200, 0.1)
  #   prior_range_change <- c(200, 0.1)
  # }
  
  # Fit model to both atlas periods
  source("Analysis/Ontario_BBA/Code/inla_model_functions.R")
  mod <- fit_inla_testing(sp_dat,  
                          model_boundary_sp, 
                          covariates = cov_df_sp,
                          prior_range_abund = c(500,0.99),  # 90% chance range is smaller than 250 km 
                          prior_sigma_abund = c(5,0.1),     # 10% chance SD is larger than 5  c(5,0.1)
                          prior_range_change = c(500,0.1),  # 10% chance range is smaller than 500 km c(500,0.1)
                          prior_sigma_change = c(0.1,0.01)  # 10% chance SD is larger than 0.1 c(0.1,0.1)
  )
  
  print(sp_english)
  print(summary(mod)) # 14 min to fit
  
  # ---------------------------------------------
  # ---- Save model summary
  # ---------------------------------------------
  
  model_summaries[[sp_english]]$summary.fixed <- mod$summary.fixed
  model_summaries[[sp_english]]$summary.hyperpar <- mod$summary.hyperpar
  saveRDS(model_summaries, "Analysis/Ontario_BBA/Output/Tables_Summaries/model_summaries.rds")
  
  # ---------------------------------------------
  # ---- Predictions onto grid across study area
  # ---------------------------------------------
  
  # Collapse covariate formulas into one string
  covariates <- cov_df_sp %>%
    mutate(formula = paste0("Beta", beta, "_", covariate, "*", covariate, "^", beta))
  
  covariate_terms <- paste(covariates$formula, collapse = " + ")
  
  # Prediction expression
  pred_expr <- paste(
    "Intercept",
    "spde_abund",
    "Atlas3 * effect_Atlas3",
    "Atlas3 * spde_change",
    "DOY_global",
    "DOY_BCR",
    covariate_terms,
    sep = " + "
  )
  
  # Assemble into single prediction formula
  pred_formula <- as.formula(
    paste0("~ data.frame(pred = ", pred_expr,")")
  )
  
  
  grid_OBBA2_model <- grid_OBBA2 
  grid_OBBA3_model <- grid_OBBA3
  
  pred_grid <- bind_rows(grid_OBBA2_model,grid_OBBA3_model) %>%
    mutate(Hours_Since_Sunrise = 0,
           days_since_june15 = -7,
           Atlas3 = ifelse(Atlas == "OBBA2",0,1)) %>%
    left_join(BCRs_to_model)
  
  # Predict relative abundance for every pixel
  start_pred <- Sys.time()
  preds <- predict_inla(mod = mod,
                        grid = pred_grid,
                        pred_formula = pred_formula)
  end_pred <- Sys.time()
  end_pred - start_pred # ~ 20 min to generate predictions
  
  # Reformat predictions
  OBBA2_indices <- which(pred_grid$Atlas == "OBBA2")
  OBBA3_indices <- which(pred_grid$Atlas == "OBBA3")
  
  preds$OBBA2 <- exp(preds$pred[OBBA2_indices,])
  preds$OBBA3 <- exp(preds$pred[OBBA3_indices,])
  preds$pred <- NULL
  
  ## Save full posterior of prediction if desired
  #if (sp_english %in% example_species){
  #  saveRDS(preds,file = paste0("model_output/species_predictions_full/",sp_english,".rds"))
  #}
  
  preds$abs_change <- preds$OBBA3 - preds$OBBA2
  
  # ---------------------------------------------
  # ---- Summarize predictions for every pixel
  # ---------------------------------------------
  
  # ---- OBBA2
  preds_OBBA2_summary <- summarize_posterior(preds$OBBA2,CI_probs = c(0.05,0.95),prefix = "OBBA2")
  
  # ---- OBBA3
  preds_OBBA3_summary <- summarize_posterior(preds$OBBA3,CI_probs = c(0.05,0.95),prefix = "OBBA3")
  
  # ---- Change between atlases (OBBA3 - OBBA2)
  preds_abs_change_summary <- summarize_posterior(preds$abs_change, prefix = "abs_change")
  
  saveRDS(list(sp_english = sp_english,
               sp_code = sp_code,
               preds_OBBA2_summary = preds_OBBA2_summary,
               preds_OBBA3_summary = preds_OBBA3_summary,
               preds_abs_change_summary = preds_abs_change_summary),
          file = paste0("Analysis/Ontario_BBA/Output/Species_Predictions/",sp_filename,".rds"))
  
  
  # ---------------------------------------------
  # ---- Create maps
  # ---------------------------------------------
  
  map_res <- max(st_distance(grid_OBBA2[1:2,]))
  
  # ---- Relative abundance during each atlas period
  
  upper_bound <- quantile(c(preds_OBBA2_summary$OBBA2_q50,preds_OBBA2_summary$OBBA3_q50),0.99,na.rm = TRUE) %>% as.numeric()
  lower_bound <- 0.01
  if (lower_bound >= upper_bound) lower_bound <- upper_bound/100
  if (lower_bound >= upper_bound/100) lower_bound <- upper_bound/100
  
  # Atlas 3
  map_relabund(species_name = sp_english,
               species_filename = sp_filename,
               obs_dat = sp_dat %>% subset(Atlas == "OBBA3"),
               grid = grid_OBBA3,
               preds_summarized = preds_OBBA3_summary,
               atlas_squares = atlas_squares,
               study_boundary = study_boundary,
               BCR_shapefile = allBCR,
               water_shapefile = water,
               map_dir = "Analysis/Ontario_BBA/Output/Prediction_Maps/Relative_Abundance/",
               train_dat_filter = "TRUE",
               prefix = "OBBA3",
               plot_obs_data = TRUE,
               title = "Atlas 3",
               subtitle = "Relative Abundance",
               upper_bound = upper_bound,
               lower_bound = lower_bound,
               res = map_res*1.01)
  
  # Atlas 2
  map_relabund(species_name = sp_english,
               species_filename = sp_filename,
               obs_dat = sp_dat %>% subset(Atlas == "OBBA2"),
               grid = grid_OBBA2,
               preds_summarized = preds_OBBA2_summary,
               atlas_squares = atlas_squares,
               study_boundary = study_boundary,
               BCR_shapefile = allBCR,
               water_shapefile = water,
               map_dir = "Analysis/Ontario_BBA/Output/Prediction_Maps/Relative_Abundance/",
               train_dat_filter = "TRUE",
               prefix = "OBBA2",
               plot_obs_data = TRUE,
               title = "Atlas 2",
               subtitle = "Relative Abundance",
               upper_bound = upper_bound,
               lower_bound = lower_bound,
               res = map_res*1.01)
  
  # Change between atlas periods
  map_change(species_name = sp_english,
             species_filename = sp_filename,
             grid = grid_OBBA3,
             preds_summarized = preds_abs_change_summary,
             study_boundary = study_boundary,
             BCR_shapefile = allBCR,
             water_shapefile = water,
             map_dir = "Analysis/Ontario_BBA/Output/Prediction_Maps/Relative_Abundance/",
             res = map_res*1.01,
             upper_bound = upper_bound,
             lower_bound = -upper_bound,
             change_type = "Absolute")
  
  # ---------------------------------------------
  # ---- Mask out water pixels prior to summarizing overall change
  # ---------------------------------------------
  
  preds$OBBA2[grid_OBBA2$on_water,] <- 0
  preds$OBBA3[grid_OBBA3$on_water,] <- 0
  
  # ---------------------------------------------
  # ---- Summarize percent change in each BCR
  # ---------------------------------------------
  
  # Initialize a list to store results
  BCR_vec <- unique(sp_dat$BCR)
  bcr_change_list <- vector("list", length(BCR_vec))
  for (j in seq_len(length(BCR_vec))) {
    
    bcr_id <- BCR_vec[j]
    
    # -----------------------------------------
    # Using all pixels in BCR
    # -----------------------------------------
    
    # Identify pixels in this BCR
    pixels_in_bcr <- which(grid_OBBA2$BCR == bcr_id)
    
    # Sum posterior draws across pixels
    sum_OBBA2 <- colSums(preds$OBBA2[pixels_in_bcr, , drop = FALSE])
    sum_OBBA3 <- colSums(preds$OBBA3[pixels_in_bcr, , drop = FALSE])
    
    # Percent change:
    percent_change <- (sum_OBBA3 - sum_OBBA2) / sum_OBBA2 * 100
    
    # -----------------------------------------
    # Probability change estimate is lower in bbs area
    # -----------------------------------------
    
    # Store results in a tibble
    bcr_change_list[[j]] <- tibble(
      BCR = bcr_id,
      
      # Overall across the BCR
      mean_change = mean(percent_change),
      median_change = median(percent_change),
      lower_change = quantile(percent_change, 0.05),
      upper_change = quantile(percent_change, 0.95),
      prob_decline = mean(percent_change<=0),
      prob_decline_30 = mean(percent_change<=-30))
    
  }
  
  # ---------------------------------------------
  # ---- Summarize overall percent change across all pixels
  # ---------------------------------------------
  
  # Sum posterior draws across all pixels
  sum_all_OBBA2 <- colSums(preds$OBBA2)
  sum_all_OBBA3 <- colSums(preds$OBBA3)
  
  # Percent change across the full study area
  overall_percent_change <- (sum_all_OBBA3 - sum_all_OBBA2) / sum_all_OBBA2 * 100
  
  overall_summary <- tibble(
    mean_change = mean(overall_percent_change),
    median_change = median(overall_percent_change),
    lower_change = quantile(overall_percent_change, 0.05),
    upper_change = quantile(overall_percent_change, 0.95),
    prob_decline = mean(overall_percent_change <= 0),
    prob_decline_30 = mean(overall_percent_change <= -30)
  )
  
  # ---------------------------------------------
  # ---- Save change results
  # ---------------------------------------------
  
  # Initialize a list for this species
  species_change_summary <- list()
  species_change_summary[[sp_english]]$sp_english <- sp_english
  species_change_summary[[sp_english]]$sp_code <- sp_code
  
  # --- Overall change ---
  species_change_summary[[sp_english]]$overall <- overall_summary
  
  # --- BCR-level summaries ---
  species_change_summary[[sp_english]]$BCR <- bcr_change_list
  
  # Save BCR-level summaries
  change_summaries <- list()
  if (file.exists("Analysis/Ontario_BBA/Output/Tables_Summaries/change_summaries.rds")) change_summaries <- readRDS("Analysis/Ontario_BBA/Output/Tables_Summaries/change_summaries.rds")
  
  # Append species summary using species English name as the list name
  change_summaries[[sp_english]] <- species_change_summary[[sp_english]]
  saveRDS(change_summaries,"Analysis/Ontario_BBA/Output/Tables_Summaries/change_summaries.rds")
  
  # # ---------------------------------------------
  # # ---- Examine HSS and DOY effects
  # # ---------------------------------------------
  # 
  # # Reference values
  # ref_doy <- 0
  # ref_hss <- 0
  # ref_bcr <- 1
  # 
  # # HSS prediction block
  # HSS_pred <- tibble(
  #   Hours_Since_Sunrise = seq(
  #     min(sp_dat$Hours_Since_Sunrise, na.rm = TRUE),
  #     max(sp_dat$Hours_Since_Sunrise, na.rm = TRUE),
  #     length.out = 100
  #   ),
  #   days_since_june15 = ref_doy,
  #   BCR_factor = ref_bcr,
  #   effect = "HSS"
  # )
  # 
  # # DOY-by-BCR prediction block
  # DOY_seq <- seq(
  #   min(sp_dat$days_since_june15, na.rm = TRUE),
  #   max(sp_dat$days_since_june15, na.rm = TRUE),
  #   by = 1
  # )
  # 
  # DOY_pred <- tidyr::expand_grid(
  #   days_since_june15 = DOY_seq,
  #   BCR_factor = sort(unique(sp_dat$BCR_factor))
  # ) %>%
  #   mutate(
  #     Hours_Since_Sunrise = ref_hss,
  #     effect = "DOY"
  #   )
  # 
  # # Combine into a single dataframe
  # pred_df <- bind_rows(HSS_pred, DOY_pred) %>%
  #   mutate(
  #     BCR_factor = if_else(
  #       effect == "HSS",
  #       ref_bcr,
  #       BCR_factor
  #     )
  #   )
  # 
  # pred_all <- predict_inla(
  #   mod,
  #   pred_df,
  #   as.formula(~ data.frame(pred = HSS + DOY_global + DOY_BCR))
  # )
  # 
  # pred_summary <- apply(pred_all$pred, 1, function(x)
  #   quantile(exp(x), c(0.025, 0.5, 0.975))) %>%
  #   t() %>%
  #   as.data.frame()
  # 
  # pred_df <- cbind(pred_df,pred_summary)
  # 
  # # Simple method for calculating peak DOY
  # peak_DOY <- pred_df %>%
  #   filter(effect == "DOY" & days_since_june15 >= -15 & days_since_june15 <= 25) %>%
  #   group_by(BCR_factor) %>%
  #   slice_max(`50%`, n = 1, with_ties = FALSE) %>%
  #   ungroup() %>%
  #   transmute(
  #     BCR_factor,
  #     days_since_june15_peak = days_since_june15
  #   )
  # 
  # # Plots
  # HSS_plot <- ggplot(subset(pred_df, effect == "HSS"), aes(x = Hours_Since_Sunrise, y = `50%`, ymin = `2.5%`, ymax = `97.5%`))+
  #   geom_ribbon(fill = "gray90") +
  #   geom_line() +
  #   theme_bw()+
  #   xlab("Hours Since Sunrise")+
  #   ylab("Relative abundance")+
  #   ggtitle(sp_english)
  # 
  # print(HSS_plot)
  # 
  # DOY_plot <- ggplot(subset(pred_df, effect == "DOY" & days_since_june15 >= -15 & days_since_june15 <= 25), aes(x = days_since_june15, y = `50%`, ymin = `2.5%`, ymax = `97.5%`))+
  #   geom_ribbon(fill = "gray90") +
  #   geom_line() +
  #   geom_vline(data = peak_DOY, aes(xintercept = days_since_june15_peak), col = "dodgerblue", linewidth = 2, alpha = 0.5)+
  #   facet_wrap(BCR_factor~., scales = "free_y")+
  #   theme_bw()+
  #   xlab("Days since June 15")+
  #   ylab("Relative abundance")+
  #   ggtitle(sp_english)
  # 
  # print(DOY_plot)
  # 
  # # Save time-of-day and day-of-year results
  # HSS_DOY_summaries <- list()
  # if (file.exists("Analysis/Ontario_BBA/Output/Tables_Summaries/HSS_DOY_summaries.rds")) HSS_DOY_summaries <- readRDS("Analysis/Ontario_BBA/Output/Tables_Summaries/HSS_DOY_summaries.rds")
  # 
  # # Append species summary using species English name as the list name
  # HSS_DOY_summaries[[sp_english]] <- pred_df
  # saveRDS(HSS_DOY_summaries,"Analysis/Ontario_BBA/Output/Tables_Summaries/HSS_DOY_summaries.rds")
  # 
  
}









# # ---------------------------------------------
# # ---- Fitted habitat covariate effects
# # ---------------------------------------------
# 
# make_stacked_df <- function(data, covariates, length_out = 50, ref_values = NULL) {
#   if (is.null(ref_values)) {
#     ref_values <- apply(data[, covariates, drop = FALSE], 2, mean, na.rm = TRUE)
#   }
#   
#   grids <- lapply(covariates, function(cv) {
#     grid <- as.data.frame(as.list(ref_values))[rep(1, length_out), ]
#     grid[[cv]] <- seq(quantile(data[[cv]],0.01, na.rm = TRUE),
#                       quantile(data[[cv]],0.99, na.rm = TRUE),
#                       length.out = length_out)
#     grid$covariate <- cv
#     grid
#   })
#   
#   df <- dplyr::bind_rows(grids)
#   rownames(df) <- NULL
#   return(df)
# }
# 
# # --- List of covariates to plot
# cov_to_plot <- unique(cov_df$covariate)
# 
# # --- prepare dataframe for which predictions must be generated
# # Optional: define reference values; if not, the mean of each covariate is used
# ref_vals <- sapply(cov_to_plot, function(cv) mean(sp_dat[[cv]], na.rm = TRUE))
# 
# stacked_df <- make_stacked_df(data = sp_dat,
#                               covariates = cov_to_plot,
#                               length_out = 50,
#                               ref_values = ref_vals)
# 
# # --- Build prediction formula for inlabru::generate
# cov_df <- cov_df %>%
#   mutate(formula_term = paste0("Beta", beta, "_", covariate, "*", covariate, "^", beta))
# covar_pred_string <- paste(cov_df$formula_term, collapse = " + ")
# covar_pred_formula <- as.formula(paste("~", covar_pred_string))
# 
# # --- Generate posterior draws
# covar_preds <- inlabru::generate(mod,
#                                  stacked_df,
#                                  formula = covar_pred_formula,
#                                  n.samples = 1000,
#                                  seed = 123)
# 
# # --- Combine stacked_df and posterior draws
# draws_df <- cbind(stacked_df, covar_preds)
# 
# # --- Initialize an empty list to store per-covariate data
# draws_list <- list()
# 
# # --- Loop over each covariate
# for (cv in cov_to_plot) {
#   
#   # --- Get numeric covariate values for the prediction grid
#   cov_values <- stacked_df[[cv]]
#   
#   # --- Flatten the posterior draws matrix for this covariate
#   # covar_preds is a matrix: rows = nrow(stacked_df), columns = draws
#   # we need to repeat cov_values for each draw
#   n_draws <- ncol(covar_preds)
#   n_values <- length(cov_values)
#   
#   df <- data.frame(
#     covar_to_plot = rep(cov_values, times = n_draws),
#     draw = rep(paste0("draw", 1:n_draws), each = n_values),
#     pred = as.vector(covar_preds)  # flatten matrix
#   )
#   
#   df$covar_to_plot <- round(df$covar_to_plot,5)
#   
#   # Remove covariate values exactly at zero
#   df <- subset(df, covar_to_plot != 0)
#   
#   # --- Add covariate name and definition/label
#   df$covariate <- cv
#   df <- df %>%
#     left_join(cov_labels, by = c("covariate" = "covariate"))
#   
#   # --- Apply exponential transformation
#   df$pred <- exp(df$pred)
#   
#   # --- Append to list
#   draws_list[[cv]] <- df
# }
# 
# # --- Combine all covariates into a single data frame
# draws_long <- dplyr::bind_rows(draws_list)
# 
# # --- Summarize for plotting
# draws_summary <- draws_long %>%
#   group_by(covar_to_plot, covariate, definition) %>%
#   summarise(
#     q05 = quantile(pred, 0.05),
#     q50 = quantile(pred, 0.5),
#     q95 = quantile(pred, 0.95),
#     .groups = "drop"
#   )
# 
# # --- Plot each covariate
# plots <- list()
# for (cv in unique(draws_summary$covariate)) {
#   cov_data <- filter(draws_summary, covariate == cv)
#   plots[[cv]] <- ggplot(cov_data, aes(x = covar_to_plot, y = q50)) +
#     geom_ribbon(aes(ymin = q05, ymax = q95), alpha = 0.3) +
#     geom_line() +
#     theme_bw() +
#     labs(
#       x = cov_data$definition[1],
#       y = "Predicted abundance",
#       title = cov_data$definition[1]
#     )
# }
# 
# # --- Combine all plots in a grid
# covar_plot <- patchwork::wrap_plots(plots, ncol = 3)

# # ---------------------------------------------
# # ---- Goodness-of-fit evaluation
# # ---------------------------------------------
# 
# nb_size <- mod$summary.hyperpar["size for the nbinomial observations (1/overdispersion)","mean"]
# kappa_var <- 1/mod$summary.hyperpar["Precision for kappa","mean"]
# 
# # -----------
# # ---- Step 1: define several types of predictions to generate for each observation
# # -----------
# 
# # Predictions using fixed effects and SPDEs (but omitting square-level random effects)
# 
# # This extracts posterior draws for every observation, omitting random effects of atlas_square
# pred_expr <- paste(
#   "Intercept",
#   "ARU * effect_ARU",
#   "spde_abund",
#   "Atlas3 * effect_Atlas3",
#   "Atlas3 * spde_change",
#   "HSS",
#   covariate_terms,
#   sep = " + "
# )
# 
# # This extracts posterior draws for every observation, simulating new random effects of atlas_square
# pred_expr_kappa <- paste(
#   "Intercept",
#   "ARU * effect_ARU",
#   "spde_abund",
#   "Atlas3 * effect_Atlas3",
#   "Atlas3 * spde_change",
#   "HSS",
#   "kappa",
#   covariate_terms,
#   sep = " + "
# )
# 
# pred_formula <- as.formula(
#   paste0("~ data.frame(preds = exp(", pred_expr,"),",
#          "preds_kappa = exp(",pred_expr_kappa,"))")
# )
# 
# # -----------
# # ---- Step 2: Generate posterior predictive draws for every point count observation
# # -----------
# 
# point_count_df <- sp_dat %>% 
#   subset(ARU == 0) %>%
#   
#   # Ensures new values of kappa (square_atlas) are simulated (if values outside the training data are used, INLA will generate new samples for them)
#   mutate(square_atlas = square_atlas + max(sp_dat$square_atlas))
# 
# ppc_preds <- predict_inla(mod = mod,
#                           grid = point_count_df,
#                           pred_formula = pred_formula)
# 
# # Simulate counts based on negative binomial error
# nb_size <- mod$summary.hyperpar["size for the nbinomial observations (1/overdispersion)","mean"]
# 
# nb_sim <- function(mat,theta){
#   matrix(rnbinom(length(mat), mu = as.vector(mat), size = theta),
#          nrow = nrow(mat),
#          ncol = ncol(mat))
# }
# 
# ppc_preds_nb <- nb_sim(ppc_preds$preds,nb_size)
# ppc_preds_kappa_nb <- nb_sim(ppc_preds$preds_kappa,nb_size)
# 
# # -----------
# # ---- Step 3: Compare observed vs predicted mean counts per atlas
# # -----------
# 
# # Compute observed mean per Atlas
# obs_means <- point_count_df %>%
#   as.data.frame() %>%
#   group_by(Atlas) %>%
#   summarize(obs_mean = mean(count), .groups = "drop")
# 
# # Combine Atlas column with posterior predictive draws
# ppc_mat <- as.data.frame(ppc_preds$preds_kappa)
# ppc_df <- cbind(as.data.frame(point_count_df)["Atlas"], ppc_mat)
# 
# # Split posterior draws by Atlas
# ppc_list <- split(ppc_df, ppc_df$Atlas)
# 
# # Function to compute posterior mean quantiles
# posterior_mean_summary <- function(df) {
#   mat <- as.matrix(df[ , -1]) # drop Atlas column
#   draw_means <- apply(mat, 2, mean) # mean per posterior draw
#   quantile(draw_means, c(0.05, 0.5, 0.95))
# }
# 
# # Apply to each Atlas
# ppc_summary <- map_dfr(names(ppc_list), function(atlas) {
#   quants <- posterior_mean_summary(ppc_list[[atlas]])
#   tibble(
#     Atlas = atlas,
#     predicted_mean_q05 = quants[1],
#     predicted_mean_q50 = quants[2],
#     predicted_mean_q95 = quants[3]
#   )
# })
# 
# # Join observed means
# obs_vs_pred_counts <- obs_means %>%
#   left_join(ppc_summary, by = "Atlas")
# 
# # -----------
# # ---- Step 4: DHARMa statistics
# # -----------
# 
# yobs <- point_count_df$count # observed counts
# pred_median <- apply(ppc_preds_kappa_nb,1,median)
# 
# # Full model posterior predictive simulations
# simObj_full <- DHARMa::createDHARMa(simulatedResponse = ppc_preds_kappa_nb, 
#                                     observedResponse = yobs,
#                                     fittedPredictedResponse = pred_median,
#                                     integerResponse = TRUE)
# 
# zero_inf <- testZeroInflation(simObj_full)          # Should be close to 1
# uniformity <- testUniformity(simObj_full)           # Should be close to 0
# 
# # overall Bayesian p-value, using Freeman-Tukey test statistic
# mu_hat <- rowMeans(ppc_preds_kappa_nb)         # length nobs
# T_obs_FT <- sum((sqrt(yobs) - sqrt(mu_hat))^2)
# T_sims_FT <- apply(ppc_preds_kappa_nb, 2, function(y_sim) sum((sqrt(y_sim) - sqrt(mu_hat))^2))
# pvalue_FT <- mean(T_sims_FT >= T_obs_FT)
# 
# # -----------
# # ---- Step 5: Correlation between predicted and observed counts (omitting square-level random effects)
# # -----------
# 
# ppc_mat <- as.data.frame(exp(log(ppc_preds$preds) + 0.5*kappa_var))
# 
# point_count_df$pred_median <- apply(ppc_mat,1,median) # Posterior medians for every observation
# 
# # Raw correlation per observation (pointwise)
# pointwise_cor <- point_count_df %>% 
#   as.data.frame() %>%
#   group_by(Atlas) %>%
#   summarize(cor_pointwise = cor(count,pred_median))
# 
# # --- Calculate square_level predictions (one calculation for each draw from posterior)
# square_preds <- point_count_df %>%
#   as.data.frame() %>%
#   select(square_id, Atlas) %>%
#   bind_cols(as.data.frame(ppc_preds$preds)) %>%
#   pivot_longer(-c(square_id, Atlas),
#                names_to="draw", values_to="pred") %>%
#   group_by(square_id, Atlas, draw) %>%
#   summarize(pred_square = mean(pred), .groups="drop")
# 
# # summarize across draws
# square_summary <- square_preds %>%
#   group_by(square_id, Atlas) %>%
#   summarize(pred_median = median(pred_square))
# 
# # compute observed square means
# square_obs <- point_count_df %>%
#   as.data.frame() %>%
#   group_by(square_id, Atlas) %>%
#   summarize(mean_obs = mean(count))
# 
# # compute correlation at square level (using median)
# squarelevel_cor <- square_summary %>%
#   left_join(square_obs, by=c("square_id","Atlas")) %>%
#   group_by(Atlas) %>%
#   summarize(cor_square = cor(mean_obs, pred_median))
# 
# # combine into single dataframe
# cor_obs_pred <- full_join(pointwise_cor,squarelevel_cor)
# 
# # -----------
# # ---- Step 6: Ability to predict presence/absence (note this omits square-level random effects)
# # -----------
# 
# point_count_df$presence <- as.numeric(point_count_df$count>0)
# 
# #  square-level random effects
# pobs <- rowMeans(nb_sim(exp(log(ppc_preds$preds) + 0.5*kappa_var),nb_size) > 0)
# roc <- pROC::roc(point_count_df$presence,pobs)
# auc <- pROC::auc(roc)
# 
# # Tjur's R-squared
# tjur_r2 <- function(y, p) {
#   if (!all(y %in% c(0, 1))) stop("y must be binary (0/1).")
#   mean(p[y == 1]) - mean(p[y == 0])
# }
# 
# # Apply to your data
# tjur <- tjur_r2(point_count_df$presence, pobs)
# tjur
# 
# # ---- Optional: show predicted probability distributions for presences vs absences
# # # Put data into a dataframe for plotting
# # plot_df <- data.frame(
# #   presence = factor(point_count_df$presence, labels = c("Absence", "Presence")),
# #   pobs = pobs
# # )
# # 
# # # Density plot
# # ggplot(plot_df, aes(x = pobs, fill = presence)) +
# #   geom_density(alpha = 0.4) +
# #   labs(
# #     x = "Predicted probability of presence",
# #     y = "Density",
# #     fill = "Observed"
# #   ) +
# #   theme_minimal()
# 
# # ---- Package goodness of fit statistics and save
# 
# gof_list <- list()
# if (file.exists("model_output/gof_linear.rds")) gof_list <- readRDS("model_output/gof_linear.rds")
# gof_list[[sp_english]] <- list(obs_vs_pred_counts = obs_vs_pred_counts,
#                                cor_obs_pred = cor_obs_pred,
#                                
#                                # DHARMA statistics
#                                zero_inf = zero_inf,
#                                uniformity = uniformity,
#                                pvalue_FT = pvalue_FT,
#                                
#                                auc = auc,
#                                tjur = tjur,
#                                
#                                pred_df = as.data.frame(point_count_df) %>% dplyr::select(count,pred_median))
# 
# saveRDS(gof_list,"model_output/gof_linear.rds")
