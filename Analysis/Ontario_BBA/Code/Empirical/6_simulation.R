# ===============================================================================
# BAYESIAN ANALYSIS / SPECIES DISTRIBUTION MODELS FOR ONTARIO BREEDING BIRD ATLAS
# ===============================================================================

library(tidyverse)
library(sf)
library(fmesher)
library(INLA)
library(MASS)  # for mvrnorm
library(terra)
library(ebirdst)
library(viridis)
library(patchwork)
library(dplyr)
library(purrr)

# Load birdDistribution package
devtools::load_all("C:/Users/IlesD/OneDrive - EC-EC/Iles/Projects/Landbirds/birdDistribution")

rm(list=ls())

# ===============================================================================
# Function to simulate spatially autocorrelated process with desired range, mean, and sd
# ===============================================================================

sim_field_fn <- function(grid,
                         surveys,
                         range_km,
                         mean_field,
                         sigma_field){
  
  grid_coords <- st_coordinates(grid)
  survey_coords <- st_coordinates(surveys)
  
  # ---- Build mesh with max edge size ~ range/5
  mesh <- inla.mesh.2d(
    loc = grid_coords[sample(1:nrow(grid_coords), 5000), ], # subsample for efficiency
    max.edge = c(range_km/5, range_km),
    cutoff = range_km/10
  )
  
  # ---- Simulate surface
  kappa <- sqrt(8) / range_km
  tau <- 1 / (sqrt(4 * pi) * kappa * sigma_field)
  
  spde <- inla.spde2.matern(
    mesh,
    alpha = 2
  )
  
  Q <- inla.spde.precision(spde, theta = c(log(tau), log(kappa)))
  x <- inla.qsample(n = 1, Q = Q)  # one realization
  
  # ---- Projection onto grid
  A_grid <- inla.spde.make.A(mesh, loc = grid_coords) # mesh projector
  grid_vals <- as.numeric(A_grid %*% x)
  
  grid_mean <- mean(grid_vals)
  grid_sd <- sd(grid_vals)
  
  # ---- Projection onto surveys
  A_surveys <- inla.spde.make.A(mesh, loc = survey_coords) # mesh projector
  survey_vals <- as.numeric(A_surveys %*% x)
  
  # ---- Ensure fields have correct mean/sd
  grid_vals <- (grid_vals/grid_sd - grid_mean) * sigma_field + mean_field
  survey_vals <- (survey_vals/grid_sd - grid_mean) * sigma_field + mean_field
  
  field_vals <- list(grid_vals = grid_vals,
                     survey_vals = survey_vals)
  return(field_vals)
}

# ===============================================================================
# Function to simulate Hours Since Sunrise effect
# ===============================================================================

simulate_activity_curve <- function(hours = seq(-2, 10, length.out = 200),
                                    peak = 0,
                                    spread = 3) {
  # Gaussian curve
  y <- exp(-((hours - peak)^2) / (2 * spread^2))
  # rescale to max = 1
  y <- y / max(y)
  y
}

# ===============================================================================
# ----  Load prepared data
# ===============================================================================

dat <- readRDS(file = "data_clean/data_ready_for_analysis.rds")
all_surveys <- dat$all_surveys

counts <- dat$counts
grid_OBBA2 <- dat$grid_OBBA2
grid_OBBA3 <- dat$grid_OBBA3
species_square_summary <- dat$species_square_summary
species_to_model <- dat$species_to_model
sp_list <- dat$sp_list
study_boundary <- dat$study_boundary
water = dat$water

# ---------------------------------------------
# ---- Load several shapefiles needed for plotting
# ---------------------------------------------

# Load shapefile containing atlas squares (10km x 10km)
atlas_squares <- st_read("C:/Users/IlesD/OneDrive - EC-EC/Iles/Projects/Landbirds/Landbird-Distribution-Modeling-ECCC/Data/Spatial/National/AtlasSquares/NationalSquares_FINAL.shp")  %>%
  subset(prov == "ON") %>%
  st_transform(st_crs(all_surveys))

# BCR boundaries
allBCR <- st_read("data/spatial/BCR/BCR_Terrestrial_master.shp") %>%
  st_transform(st_crs(study_boundary)) %>%
  st_make_valid() %>%
  st_buffer(0) %>%
  group_by(BCR,BCRNAME,PROVINCE_S) %>%
  summarise() %>%
  st_transform(st_crs(all_surveys))

BCR <- st_read("data/spatial/BCR/BCR_Terrestrial_master.shp") %>%
  subset(PROVINCE_S == "ONTARIO" & BCR %in% c(13,12,8)) %>%
  st_transform(st_crs(study_boundary)) %>%
  st_make_valid() %>%
  st_buffer(0) %>%
  group_by(BCR,BCRNAME) %>%
  summarise() %>%
  st_transform(st_crs(all_surveys))

# ---------------------------------------------
# ---- Covariates to include in model
# ---------------------------------------------

covars <- c("insect_broadleaf",
            "insect_needleleaf",
            "road",
            "water_river",
            "prec",
            "tmax",
            "urban_2",
            "urban_3",
            "lc_1",  # Evergreen Needleaf Forest
            "lc_4",  # Deciduous Broadleaf Forest
            "lc_5",  # Mixed Forest
            "lc_8",  # Woody Savanna; Tree cover 30-60% (canopy >2m
            "lc_9",  # Savanna; Tree cover 10-30% (canopy >2m)
            "lc_10", # Grassland
            "lc_11", # Permanent Wetland
            "lc_12", # Cropland
            "lc_14",  # Cropland / Natural Vegetation Mosaic'
            "lc_17"  # Water
)

# ---------------------------------------------
# ---- Ensure all_surveys uses covariates from grid_sf
# ---------------------------------------------

grid_sf <- grid_OBBA3
all_surveys <- all_surveys %>% dplyr::select(-all_of(covars))

# find index of nearest grid point for each survey point
nearest_idx <- st_nearest_feature(all_surveys, grid_sf)

# pull covariate values from grid_sf
covar_vals <- grid_sf %>%
  st_drop_geometry() %>%   # drop geometry, keep attributes
  slice(nearest_idx) %>%
  select(all_of(covars))

all_surveys <- all_surveys %>%
  bind_cols(covar_vals)

# ---------------------------------------------
# ---- Specify priors for covariate effects
# ---------------------------------------------

sd_linear <- 0.5
cov_df <- data.frame(covariate = covars, 
                     beta = 1,
                     model = "linear",
                     mean = 0, prec = 1 / sd_linear^2) 
sp_list <- dat$sp_list

sim_results <- list()
for (sp_english in sp_list){
  
  if (sp_english == "Yellow-bellied Sapsucker") next
  
  print(sp_english)
  if (file.exists("model_output/sim_results_linear.rds")) sim_results <- readRDS("model_output/sim_results_linear.rds")
  if (sp_english %in% names(sim_results)) next
  sp_code <- subset(species_to_model, english_name == sp_english)$species_id[1] %>% as.character()
  if (is.na(sp_code)) next
  
  # Parameters controlling the spatial field for population change
  range_km = runif(1,100,500)
  mean_field = runif(1,log(0.5),log(2))
  sigma_field = 0.5
  nb_size = 1
  
  # ===============================================================================
  # Simulate initial abundance using eBird surface
  # ===============================================================================
  
  # ---- Use species relative abundance surface from eBird
  
  # Load ebird species relative abundance surface (3 km resolution)
  sp_code_ebird <- subset(ebirdst::ebirdst_runs, common_name == sp_english)$species_code
  
  if ("yebsap" %in% sp_code_ebird) sp_code_ebird == "yebsap"
  
  # Check if species is available
  check <- get_species(sp_code_ebird)
  if (length(check)==0){
    print(paste0(sp_english, " not available in eBird"))
    next
  }
  
  if (length(check)>0 & is.na(check)){
    print(paste0(species_name, " not available in eBird"))
    next
  }
  
  # ebirdst_download_status(sp_code_ebird,
  #                         download_abundance = TRUE,
  #                         download_ranges = FALSE,
  #                         pattern = "abundance_median_3km_")
  
  path <- get_species_path(sp_code_ebird)
  
  ebird_raster <- load_raster(sp_code_ebird, product = "abundance", period = "seasonal",resolution = "3km")
  
  if ("breeding" %in% names(ebird_raster)){
    ebird_raster <- ebird_raster[["breeding"]]
  } else{
    ebird_raster <- ebird_raster[["resident"]]
  }
  
  # ---- Sample the raster at each location in grid_sf and all_surveys
  # Make sure CRS are consistent
  if (crs(grid_sf) != crs(ebird_raster)) grid_sf <- st_transform(grid_OBBA3, crs(ebird_raster))
  if (crs(all_surveys) != crs(ebird_raster)) all_surveys <- st_transform(all_surveys, crs(ebird_raster))
  
  sp_dat <- st_transform(all_surveys, crs(ebird_raster))
  
  # Extract values for grid_sf
  grid_vals <- terra::extract(ebird_raster, vect(grid_sf))
  grid_sf$ebird <- grid_vals[,2]   # first column is ID, second is value
  
  # Extract values for sp_dat
  survey_vals <- terra::extract(ebird_raster, vect(sp_dat))
  sp_dat$ebird <- survey_vals[,2]
  
  # Convert both objects back to original crs
  grid_sf <- st_transform(grid_sf, st_crs(all_surveys))
  sp_dat <- st_transform(sp_dat, st_crs(all_surveys))
  
  grid_sf$ebird[is.na(grid_sf$ebird)] <- 0
  sp_dat$ebird[is.na(sp_dat$ebird)] <- 0
  
  # Rescale relative abundances surfaces so that counts will be realistic
  target_mean_count <- mean(counts[[sp_code]],na.rm = TRUE)
  ebird_mean_count <- mean(sp_dat$ebird,na.rm = TRUE)
  
  sp_dat$ebird  <- sp_dat$ebird/ ebird_mean_count * target_mean_count
  grid_sf$ebird <- grid_sf$ebird/ ebird_mean_count * target_mean_count
  
  # ===============================================================================
  # Simulate change surface across the province, with random intercept and range
  # ===============================================================================
  
  change_sim <- sim_field_fn(grid = grid_sf %>% st_transform(st_crs(grid_OBBA3)),
                             surveys = sp_dat %>% st_transform(st_crs(grid_OBBA3)),
                             range_km = range_km,
                             mean_field = mean_field,
                             sigma_field = sigma_field)
  
  # ===============================================================================
  # Simulate data for this species
  # ===============================================================================
  
  # Survey counts
  sp_dat$log_mu <- log(sp_dat$ebird) + simulate_activity_curve(sp_dat$Hours_Since_Sunrise, peak = 0, spread = 6) + ifelse(sp_dat$Atlas == "OBBA3",1,0) * change_sim$survey_vals
  sp_dat$expected_count <- exp(sp_dat$log_mu)
  sp_dat$count <- rnbinom(nrow(sp_dat),mu = sp_dat$expected_count, size = nb_size)
  
  plot_sim_surveys <- ggplot(sp_dat)+
    geom_sf(data = study_boundary, col = "gray80", fill = "transparent")+
    geom_sf(aes(col = count))+
    scale_color_gradientn(colours = c("transparent","black"))+
    facet_grid(.~Atlas)+
    theme_bw()
  
  # Expected counts (at sunrise) in each 1 km pixel
  grid_sf$log_mu_OBBA2 <- log(grid_sf$ebird)
  grid_sf$log_mu_OBBA3 <- log(grid_sf$ebird) + change_sim$grid_vals
  
  grid_sf$expected_count_OBBA2 <- exp(grid_sf$log_mu_OBBA2)
  grid_sf$expected_count_OBBA3 <- exp(grid_sf$log_mu_OBBA3)
  
  lim <- range(c(grid_sf$expected_count_OBBA2, grid_sf$expected_count_OBBA3))
  plot_sim_OBBA2 <- ggplot(grid_sf)+
    geom_sf(data = study_boundary, col = "gray80", fill = "transparent")+
    geom_sf(aes(col = expected_count_OBBA2))+
    scale_color_gradientn(colours = c("transparent","black"), limits = lim, name = "Expected\ncount")+
    theme_bw()+
    ggtitle("OBBA2")
  
  plot_sim_OBBA3 <- ggplot(grid_sf)+
    geom_sf(data = study_boundary, col = "gray80", fill = "transparent")+
    geom_sf(aes(col = expected_count_OBBA3))+
    scale_color_gradientn(colours = c("transparent","black"), limits = lim, name = "Expected\ncount")+
    theme_bw()+
    ggtitle("OBBA3")
  
  plot_grid_combined <- plot_sim_OBBA2 + plot_sim_OBBA3
  # print(plot_grid_combined)
  
  # ===============================================================================
  # Analyze simulated data
  # ===============================================================================
  
  # Fit model to both atlas periods
  sp_dat <- st_transform(sp_dat, st_crs(grid_OBBA3))
  mod <- fit_inla(sp_dat,  
                  study_boundary, 
                  covariates = cov_df,
                  prior_range_abund = c(250,0.5),  # 50% chance range is smaller than 250 km
                  prior_sigma_abund = c(1,0.1),    # 10% chance SD is larger than 1
                  prior_range_change = c(250,0.5), # 50% chance range is smaller than 250 km
                  prior_sigma_change = c(1,0.1)    # 10% chance SD is larger than 1
  ) #~ 20 min
  
  # ---------------------------------------------
  # ---- Predictions onto grid
  # ---------------------------------------------
  
  # Combine the grids for both atlases into a single object and generate predictions for both
  # simultaneously (required to preserve dependencies in joint posterior for change mapping)
  grid_sf <- st_transform(grid_sf, st_crs(grid_OBBA3))
  grid_sf_Atlas2 <- grid_sf %>% mutate(Atlas = "OBBA2")
  grid_sf_Atlas3 <- grid_sf %>% mutate(Atlas = "OBBA3")
  
  pred_grid <- bind_rows(grid_sf_Atlas2,grid_sf_Atlas3) %>% 
    mutate(Hours_Since_Sunrise = 0,
           Atlas3 = ifelse(Atlas == "OBBA2",0,1))
  
  # ---------------------------------------------
  # ---- Predictions onto grid across study area
  # ---------------------------------------------
  
  # Collapse covariate formulas into one string
  covariates <- cov_df %>%
    mutate(formula = paste0("Beta", beta, "_", covariate, "*", covariate, "^", beta))
  
  covariate_terms <- paste(covariates$formula, collapse = " + ")
  
  # Prediction expression
  pred_expr <- paste(
    "Intercept",
    "spde_abund",
    "Atlas3 * effect_Atlas3",
    "Atlas3 * spde_change",
    covariate_terms,
    sep = " + "
  )
  
  # Assemble into single prediction formula
  pred_formula <- as.formula(
    paste0("~ data.frame(pred = ", pred_expr,")")
  )
  
  
  # Predict relative abundance and change for every pixel
  preds <- predict_inla(mod = mod,
                        grid = pred_grid,
                        pred_formula = pred_formula)
  
  # Reformat predictions
  OBBA2_indices <- which(pred_grid$Atlas == "OBBA2")
  OBBA3_indices <- which(pred_grid$Atlas == "OBBA3")
  
  preds$OBBA2 <- exp(preds$pred[OBBA2_indices,])
  preds$OBBA3 <- exp(preds$pred[OBBA3_indices,])
  preds$pred <- NULL
  preds$abs_change <- preds$OBBA3 - preds$OBBA2
  
  
  # ---------------------------------------------
  # ---- Summarize percent change in each BCR
  # ---------------------------------------------
  
  # Initialize a list to store results
  bcr_summary_list <- vector("list", nrow(BCR))
  
  for (i in seq_len(nrow(BCR))) {
    
    bcr_id <- BCR$BCR[i]
    
    # -----------------------------------------
    # Using all pixels in BCR
    # -----------------------------------------
    
    # Identify pixels in this BCR
    pixels_in_bcr <- st_intersects(grid_OBBA2, BCR[i, ], sparse = FALSE)[,1]
    
    # Truth (simulated)
    sum_OBBA2_true <- sum(grid_sf$expected_count_OBBA2[pixels_in_bcr])
    sum_OBBA3_true <- sum(grid_sf$expected_count_OBBA3[pixels_in_bcr])
    percent_change_true <- (sum_OBBA3_true - sum_OBBA2_true) / sum_OBBA2_true * 100
    
    # Model-based estimates
    sum_OBBA2_est <- colSums(preds$OBBA2[pixels_in_bcr, , drop = FALSE])
    sum_OBBA3_est <- colSums(preds$OBBA3[pixels_in_bcr, , drop = FALSE])
    percent_change_est <- (sum_OBBA3_est - sum_OBBA2_est) / sum_OBBA2_est * 100
    
    # Store results in a tibble
    bcr_summary_list[[i]] <- tibble(
      
      BCR = bcr_id,
      
      # Atlas 2
      sum_OBBA2_true = sum_OBBA2_true,
      sum_OBBA2_est_q50 = quantile(sum_OBBA2_est,0.5),
      sum_OBBA2_est_q05 = quantile(sum_OBBA2_est,0.05),
      sum_OBBA2_est_q95 = quantile(sum_OBBA2_est,0.95),
      
      # Atlas 3
      sum_OBBA3_true = sum_OBBA3_true,
      sum_OBBA3_est_q50 = quantile(sum_OBBA3_est,0.5),
      sum_OBBA3_est_q05 = quantile(sum_OBBA3_est,0.05),
      sum_OBBA3_est_q95 = quantile(sum_OBBA3_est,0.95),
      
      # Percent change
      percent_change_true =  percent_change_true,
      percent_change_est_q50 = quantile(percent_change_est,0.5),
      percent_change_est_q05 = quantile(percent_change_est,0.05),
      percent_change_est_q95 = quantile(percent_change_est,0.95)
      )
    
  }
  
  # ---------------------------------------------
  # ---- Summarize overall percent change across all pixels
  # ---------------------------------------------
  
  # Truth (simulated)
  sum_overall_OBBA2_true <- sum(grid_sf$expected_count_OBBA2)
  sum_overall_OBBA3_true <- sum(grid_sf$expected_count_OBBA3)
  overall_percent_change_true <- (sum_overall_OBBA3_true - sum_overall_OBBA2_true) / sum_overall_OBBA2_true * 100
  
  # Model-based estimates
  sum_overall_OBBA2_est <- colSums(preds$OBBA2)
  sum_overall_OBBA3_est <- colSums(preds$OBBA3)
  overall_percent_change_est <- (sum_overall_OBBA3_est - sum_overall_OBBA2_est) / sum_overall_OBBA2_est * 100
  
  overall_summary <- tibble(
    
    # Atlas 2
    sum_overall_OBBA2_true = sum_overall_OBBA2_true,
    sum_overall_OBBA2_est_q50 = quantile(sum_overall_OBBA2_est,0.5),
    sum_overall_OBBA2_est_q05 = quantile(sum_overall_OBBA2_est,0.05),
    sum_overall_OBBA2_est_q95 = quantile(sum_overall_OBBA2_est,0.95),
    
    # Atlas 3
    sum_overall_OBBA3_true = sum_overall_OBBA3_true,
    sum_overall_OBBA3_est_q50 = quantile(sum_overall_OBBA3_est,0.5),
    sum_overall_OBBA3_est_q05 = quantile(sum_overall_OBBA3_est,0.05),
    sum_overall_OBBA3_est_q95 = quantile(sum_overall_OBBA3_est,0.95),
    
    # Percent change
    overall_percent_change_true =  overall_percent_change_true,
    overall_percent_change_est_q50 = quantile(overall_percent_change_est,0.5),
    overall_percent_change_est_q05 = quantile(overall_percent_change_est,0.05),
    overall_percent_change_est_q95 = quantile(overall_percent_change_est,0.95)
  )
  
  # ---------------------------------------------
  # ---- Save change results
  # ---------------------------------------------
  
  # Initialize a list for this species
  species_results <- list(range_km = range_km,
                          mean_field = mean_field,
                          sigma_field = sigma_field,
                          nb_size = nb_size,
                          bcr_summary_list = bcr_summary_list,
                          overall_summary = overall_summary)

  # Save BCR-level summaries
  if (file.exists("model_output/sim_results_linear.rds")) sim_results <- readRDS("model_output/sim_results_linear.rds")
  
  # Append species summary using species English name as the list name
  sim_results[[sp_english]] <- species_results
  saveRDS(sim_results,"model_output/sim_results_linear.rds")
  
}


# ---------------------------------------
# Combine results and plot
# ---------------------------------------

library(ggplot2)
library(dplyr)
library(scales)  # for minor_breaks

sim_results <- readRDS("model_output/sim_results_linear.rds")
combined_df <- map_dfr(
  names(sim_results),
  ~ sim_results[[.x]]$overall_summary %>%
    as.data.frame() %>%
    mutate(sp_english = .x)
) %>% relocate(sp_english)

#combined_df <- subset(combined_df, overall_percent_change_est_q95 < 500)

# Signed log transform function
signed_log <- function(x) {
  sign(x) * log1p(abs(x)/100)   # divide by 100 if x is in percent
}

# Inverse transform for axis labels
inv_signed_log <- function(x) {
  sign(x) * (exp(abs(x)) - 1) * 100
}

# Define pretty breaks for the log scale
pretty_breaks_signed_log <- function(limits) {
  # Compute breaks on the transformed scale
  raw <- c(-90, -50, -20, 0, 20, 50, 100, 200, 400)
  signed_log(raw)
}

combined_df <- combined_df %>%
  mutate(
    overall_percent_change_true_log = signed_log(overall_percent_change_true),
    overall_percent_change_est_q50_log = signed_log(overall_percent_change_est_q50),
    overall_percent_change_est_q05_log = signed_log(overall_percent_change_est_q05),
    overall_percent_change_est_q95_log = signed_log(overall_percent_change_est_q95),
    
    cov = overall_percent_change_est_q95_log > overall_percent_change_true_log & overall_percent_change_est_q05_log < overall_percent_change_true_log
  )

# Define reasonable axis limits on the signed-log scale
xlim_log <- c(signed_log(-90), signed_log(400))  # example: -90% to +400%
ylim_log <- c(signed_log(-90), signed_log(400))

ggplot(combined_df) +
  geom_abline(slope = 1, intercept = 0, col = "gray80") +
  geom_point(aes(x = overall_percent_change_true_log, y = overall_percent_change_est_q50_log, col = cov)) +
  geom_errorbar(aes(x = overall_percent_change_true_log,
                    ymin = overall_percent_change_est_q05_log,
                    ymax = overall_percent_change_est_q95_log,col = cov), width = 0) +
  scale_x_continuous(
    breaks = pretty_breaks_signed_log(xlim_log),
    labels = function(x) paste0(round(inv_signed_log(x)), "%")
  ) +
  scale_y_continuous(
    breaks = pretty_breaks_signed_log(ylim_log),
    labels = function(x) paste0(round(inv_signed_log(x)), "%")
  ) +
  coord_cartesian(xlim = xlim_log, ylim = ylim_log) +  # clip the axes but keep error bars
  xlab("True (simulated) percent population change") +
  ylab("Estimated percent population change") +
  ggtitle("True (simulated) vs model-estimated percent population change\nfrom Atlas 2 to Atlas 3") +
  theme_bw() +
  theme(
    panel.grid.minor = element_line(color = "gray90", linetype = "dotted")
  )+
  scale_color_manual(values = c("red","black"), guide = "none")

mean(combined_df$overall_percent_change_est_q50 > combined_df$overall_percent_change_true)
mean(combined_df$overall_percent_change_est_q95 > combined_df$overall_percent_change_true & combined_df$overall_percent_change_est_q05 < combined_df$overall_percent_change_true)
mean(combined_df$overall_percent_change_est_q50_log - combined_df$overall_percent_change_true_log)

# lim <- range(combined_df[,c("overall_percent_change_true","overall_percent_change_est_q05","overall_percent_change_est_q95")])
# ggplot(combined_df)+
#   geom_abline(slope = 1, intercept = 0, col = "gray80")+
#   geom_point(aes(x = overall_percent_change_true, y = overall_percent_change_est_q50))+
#   geom_errorbar(aes(x = overall_percent_change_true, ymin = overall_percent_change_est_q05, ymax = overall_percent_change_est_q95), width = 0)+
#   ylab("Estimated percent population change")+
#   xlab("True (simulated) percent population change")+
#   ggtitle("True (simulated) vs model-estimated percent population change\nfrom Atlas 2 to Atlas 3")+
#   coord_cartesian(ylim=lim,xlim=lim)+
#   theme_bw()