# ============================================================
# 13_evaluate_sampling_artifacts.R
#
# Purpose:
#
# Inputs:
#   data_clean/birds/data_ready_for_analysis.rds  (from script 07)
#
# Outputs (written incrementally per species):
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
figure_utils_path <- file.path(paths$functions, "figure_utils.R")

source(inla_utils_path)
source(figure_utils_path)


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

sim_dir <- file.path(out_dir, paste0("simulations_",model_name))
dir.create(sim_dir, recursive = TRUE, showWarnings = FALSE)

fig_dir <- file.path(sim_dir,"simulation_figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

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
# PREPARE SIMULATION GRID WHERE NOTHING ON LANDSCAPE HAS CHANGED
# Ensure surveys match this
# ------------------------------------------------------------

grid2_sim <- grid_OBBA2 %>% mutate(Atlas = "OBBA2")
grid3_sim <- grid_OBBA2 %>% mutate(Atlas = "OBBA3")

# Pick covariate columns that exist on grid2
covars_present <- intersect(base_covars, names(grid2_sim))






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
  sample_n(nrow(.)) %>%
  filter(english_name %in% c("Palm Warbler","Common Nighthawk","Connecticut Warbler",
                             "Greater Yellowlegs","Savannah Sparrow","Olive-sided Flycatcher",
                             "Evening Grosbeak","Bobolink","Tree Swallow"))

message("Species queued: ", nrow(species_run))
species_run

max_runs <- 1

for (i in 1:nrow(species_run)) {
  
  sp_name <- species_run$english_name[i]
  sp_code <- as.character(species_run$species_id[i])
  sp_file <- sp_filename(sp_name)
  
  int_strategy = "auto" # can be faster with "eb", but more likley to fail
  strategy = "laplace"  # can be faster with "simplified.laplace" but more approximate
  
  message("\n====================\n", i, "/", nrow(species_run), ": ", sp_name, " (species_id = ", sp_code, ")\n====================")
  
  if (file.exists(file.path(sim_dir,"simulation_summaries", paste0(sp_file,"_",max_runs,".rds")))) next
  
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
  
  mod <- fit_inla_multi_atlas(
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
        
        family = error_family)
    
  end_model <- Sys.time()
  
  fit_minutes <- round(as.numeric(end_model - start_model, units = "mins"))
  
  message(
    "\n====================\n",
    i, "/", nrow(species_run), ": ", sp_name,
    " (species_id = ", sp_code, "): ",
    fit_minutes, " minutes to fit model\n",
    "===================="
  )
  
  
  
  # ------------------------------------------------------------
  # Simulate new data (where no population has occurred) and fit model to simulated data
  # ------------------------------------------------------------
  
  crs_dat <- st_crs(sp_dat)
  sp_sim <- sp_dat
  
  # Find nearest grid2_sim index for each survey point
  # (Requires same CRS; if not, transform one to match the other.)
  if (st_crs(sp_sim) != st_crs(grid2_sim)) {
    sp_sim <- st_transform(sp_sim, st_crs(grid2_sim))
  }

  nn_idx <- st_nearest_feature(sp_sim, grid2_sim)

  # Build a covariate table aligned to surveys
  cov_tbl <- st_drop_geometry(grid2_sim[nn_idx, ]) %>%
    select(all_of(covars_present))

  # Overwrite / attach covariates on surveys
  sp_sim <- sp_sim %>%
    # drop any existing versions of those covariates to avoid duplicates
    select(-any_of(covars_present)) %>%
    bind_cols(cov_tbl)

  # Simulate new data (assuming no change between atlases)
  sp_sim <- sp_sim %>%
    mutate(
      BCR_factor = as.numeric(factor(BCR)),
      BCR_idx = as.integer(BCR),

      Atlas3 = 0,  # Set to 0 for simulation so that all data are simulated for atlas 2
      Atlas3_c = Atlas3 - 0.5,

      # New levels for squares (so that new square-level random effects are also simulated)
      square_atlas = square_atlas + max(square_atlas) # This ensures new random square effects are generated
    ) %>%
    st_transform(crs_dat)
  
  sim_formula <- make_pred_formula_multiatlas(cov_df_sp, include_kappa = FALSE, include_aru = TRUE)
  simdat <- generate(mod, newdata = sp_sim, formula = sim_formula, n.samples = 10)
  
  # ------------------------------------------------
  # Simulation runs
  # ------------------------------------------------
  
  run_number <- 1 # Only doing 1 run for now
  
  sim_path <- file.path(sim_dir,"simulation_summaries", paste0(sp_file,"_",run_number,".rds"))
  
  if (file.exists(sim_path)){
    message("Simulation already exists for ", sp_name, "; skipping")
    next
  }
  
  # Simulate new square-level random effects
  kappa_sd <- sqrt(1/mod$summary.hyperpar["Precision for kappa_diff","mean"])
  kappa_sim <- rnorm(max(sp_sim$square_atlas),0,kappa_sd)
  sp_sim$kappa_new <- kappa_sim[sp_sim$square_atlas]
  
  # Simulate either poisson or nbinomial error in counts
  if (error_family == "poisson"){
    new_counts = rpois(nrow(sp_dat),exp(simdat[[1]]$eta + sp_sim$kappa_new))
  }
  
  if (error_family == "nbinomial"){
    size_PC <- summary(mod)$inla$hyperpar[1,"mean"]
    size_ARU <- summary(mod)$inla$hyperpar[2,"mean"]
    size_vec <- ifelse(sp_dat$ARU==0,size_PC,size_ARU)
    new_counts = rnbinom(nrow(sp_dat),mu = exp(simdat[[1]]$eta + sp_sim$kappa_new), size = size_vec)
  }
  
  # Prepare for refitting
  sp_dat_refit <- sp_sim %>% 
    mutate(count = new_counts,
           Atlas3 = ifelse(Atlas == "OBBA2",0,1), 
           Atlas3_c = Atlas3 - 0.5,
           square_atlas = sp_dat$square_atlas)
  
  # Plot new simulated data
  sim_plot <- ggplot(sp_dat_refit) + 
    geom_sf(col = "gray50")+
    geom_sf(data = subset(sp_dat_refit,count>0), col = "black")+
    facet_grid(Atlas~.)+
    theme_bw()+
    ggtitle(sp_name)
  
  print(sim_plot)

  mod_refit <- fit_inla_multi_atlas(
      sp_dat = sp_dat_refit,
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
      
      # atlas square iid effect
      kappa_pcprec_diff   = priors_list$kappa_pcprec_diff,
      
      int_strategy = int_strategy,
      strategy = strategy,
      
      family = error_family)
  
  summary(mod)
  summary(mod_refit)
  
  # ------------------------------------------------------------
  # Generate predictions for every pixel based on mod_refit
  # ------------------------------------------------------------
  
    message("Generating simulated data")
    
    # Create an 
    pred_grid <- bind_rows(
      grid2_sim %>% mutate(Atlas3 = 0),
      grid3_sim %>% mutate(Atlas3 = 1)
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
      mod = mod_refit,
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
    
    message(
      "\n====================\n",
      i, "/", nrow(species_run), ": ", sp_name,
      " (species_id = ", sp_code, "): Simulation complete",
      "===================="
    )
    
    # Calculate zlim (rather than storing OBBA2 and OBBA3 relative abundance predictions)
    zmax2 <- as.numeric(stats::quantile(preds_OBBA2_summary$OBBA2_q50, 0.99, na.rm = TRUE))
    zmax3 <- as.numeric(stats::quantile(preds_OBBA3_summary$OBBA3_q50, 0.99, na.rm = TRUE))
    zmax <- max(zmax2, zmax3, na.rm = TRUE)
    
    # Save model summaries and predictions
    save_atomic(list(

      # Species data
      sp_name = sp_name,
      sp_code = sp_code,
      run_number = run_number,
      
      mod_summary_empirical = list(hyperpar = summary(mod)$inla$hyperpar,
                              fixed = summary(mod)$inla$fixed),
      mod_summary_simulation = list(hyperpar = summary(mod_refit)$inla$hyperpar,
                                   fixed = summary(mod_refit)$inla$fixed),
      
      zmax = zmax,
      abs_change = preds_abs_change_summary,   # Per-pixel summaries
      eta_draws_per_hex = eta_draws_per_hex    # Posterior draws within each 25-km hexagon across Ontario

    ), sim_path)
  
}

message("\nSimulations complete.")


# ------------------------------------------------------------
# Process / summarize simulation results
# ------------------------------------------------------------

# Plot export settings
dpi <- 1000
width_in <- 10
height_in <- 8
ggsave_type <- "cairo"

# Raster resolution in map units (meters in EPSG:3978)
plot_res <- 1001

sim_files <- list.files(file.path(sim_dir,"simulation_summaries"), pattern = "\\.rds$", full.names = TRUE)
stopifnot(length(sim_files) > 0)

study_boundary <- st_transform(study_boundary, st_crs(grid_OBBA2))
message("Simulation files found: ", length(sim_files))


for (i in seq_len(length(sim_files))) {
  
  preds <- readRDS(sim_files[i])
  
  # ----------------------------------------------------------
  # Basic identifiers
  # ----------------------------------------------------------
  sp_english <- preds$sp_name
  sp_file <- sp_english |>
    stringr::str_to_lower() |>
    stringr::str_replace_all("[^a-z0-9]+", "_") |>
    stringr::str_replace_all("^_|_$", "")
  
  sp_code <- dat$species_to_model |>
    dplyr::filter(english_name == sp_english) |>
    dplyr::pull(species_id)
  
  run_number <- preds$run_number
  
  fig_path <- file.path(fig_dir, paste0(sp_file, "_",run_number,".png"))
  
  # Skip if already run
  if (file.exists(fig_path)) next
  
  message("Mapping: ", sp_english," simulation run ",run_number)
  
  # ----------------------------------------------------------
  # Shared relative-abundance/change plotting scale
  # ----------------------------------------------------------
  
  zmax <- preds$zmax
  rel_bounds <- list(lower = 0, upper = zmax)
  
  # ----------------------------------------------------------
  # Meaningful change polygons
  # ----------------------------------------------------------
  
  signif_change_polys <- build_meaningful_change_polys(
    eta_draws_per_hex = preds$eta_draws_per_hex,
    hexagon_sf = hex_grid,
    study_boundary = study_boundary,
    min_area_km2 = 10000,  # Area of contiguous polygons with "significant" change to be plotted
    smoothing_bandwidth_m = 75000,
    param = "abs_change",
    threshold = zmax / 20,
    prob_level = 0.90,
    direction = c("two_sided", "increase", "decrease"),
    ci_probs = c(0.05, 0.95),
    include_summary = TRUE,
    drop_holes = TRUE
  )
  
  chg <- make_abs_change_maps(
    species_name = sp_english,
    grid_sf = grid_OBBA2,
    abs_change_summary = preds$abs_change,
    study_boundary = study_boundary,
    bcr_sf = NULL,
    water_sf = NULL,
    signif_poly = signif_change_polys,
    res = plot_res,
    max_abs = zmax,
    title = "Absolute change",
    subtitle = "OBBA2 to OBBA3"
  )
  

  ggsave(
    filename =fig_path,
    plot = chg$chg_plot, width = width_in, height = height_in, units = "in",
    dpi = dpi, type = ggsave_type, limitsize = FALSE
  )
  
  
}

message("10_make_species_figures.R complete: ", fig_dir)