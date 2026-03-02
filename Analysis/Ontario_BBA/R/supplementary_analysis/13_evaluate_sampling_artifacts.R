# ============================================================
# 13_evaluate_sampling_artifacts.R
#
# Purpose:
#   Use already-fitted JOINT (multi-atlas) models as data-generating truth,
#   but simulate a "no change" dataset at the REAL survey locations/timing/design,
#   with covariates fixed to OBBA2 (grid2) values.
#
# Implements:
#   1) Overwrite grid3 covariates with grid2 covariates
#   2) Re-attach covariates to all_surveys via nearest grid2 point
#   3) Detect fitted models already on disk
#   4) Loop through species, load fitted model
#   5) Simulate new observations at all survey events
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

source("R/functions/inla_model_utils.R")

# ------------------------------------------------------------
# Config
# ------------------------------------------------------------
in_file <- "data_clean/birds/data_ready_for_analysis.rds"

out_dir <- "data_clean/model_output"
models_dir <- file.path(out_dir, "models_joint2")
dir.create(file.path(out_dir, "sim_no_change"), recursive = TRUE, showWarnings = FALSE)

# number of posterior draws used for simulation
n_draws_sim <- 1  # set to 1 for one simulated dataset per species; increase for repeats
seed_base <- 123

# Covariates that may be included
base_covars <- c(
  "on_river_N","on_river_S","on_road","urban_3",
  "lc_1S","lc_1N","lc_4S","lc_4N","lc_5S","lc_5N",
  "lc_8S","lc_8N","lc_9S","lc_9N","lc_10S","lc_10N",
  "lc_11","lc_12","lc_14","lc_17",
  "insect_broadleaf","insect_needleleaf"
)

# ------------------------------------------------------------
# Load data
# ------------------------------------------------------------
stopifnot(file.exists(in_file))
dat <- readRDS(in_file)

all_surveys <- dat$all_surveys
counts      <- dat$counts
grid2       <- dat$grid_OBBA2
grid3       <- dat$grid_OBBA3
species_to_model <- dat$species_to_model
study_boundary <- dat$study_boundary %>% st_as_sf()

# Ensure sf
all_surveys <- st_as_sf(all_surveys) %>% st_transform(st_crs(study_boundary))
grid2 <- st_as_sf(grid2) %>% st_transform(st_crs(study_boundary))
grid3 <- st_as_sf(grid3) %>% st_transform(st_crs(study_boundary))

# ------------------------------------------------------------
# 1) "Convert grid3 to grid2"
#    - simplest interpretation: use grid2 for both atlases in simulation
#    - If you need a grid3 object downstream, just duplicate grid2
# ------------------------------------------------------------
grid3_sim <- grid2

# (Optional) keep Atlas label if useful later
grid2 <- grid2 %>% mutate(Atlas = "OBBA2")
grid3_sim <- grid3_sim %>% mutate(Atlas = "OBBA3")

# ------------------------------------------------------------
# 2) Re-attach covariates to all_surveys from nearest grid2 point
#    - This overwrites atlas-specific covariate columns in all_surveys
# ------------------------------------------------------------

# Pick covariate columns that exist on grid2
covars_present <- intersect(base_covars, names(grid2))

# Find nearest grid2 index for each survey point
# (Requires same CRS; if not, transform one to match the other.)
if (st_crs(all_surveys) != st_crs(grid2)) {
  all_surveys <- st_transform(all_surveys, st_crs(grid2))
}

nn_idx <- st_nearest_feature(all_surveys, grid2)

# Build a covariate table aligned to surveys
cov_tbl <- st_drop_geometry(grid2[nn_idx, ]) %>%
  select(all_of(covars_present))

# Overwrite / attach covariates on surveys
surveys_sim_design <- all_surveys %>%
  # drop any existing versions of those covariates to avoid duplicates
  select(-any_of(covars_present)) %>%
  bind_cols(cov_tbl)

# Keep derived covariates
surveys_sim_design <- surveys_sim_design %>%
  mutate(
    days_since_june15 = DayOfYear - 166,
    BCR_factor = as.numeric(factor(BCR)),
    BCR_idx = as.integer(BCR),
    Atlas3 = ifelse(Atlas == "OBBA2", 0, 1),
    Atlas3_c = Atlas3 - 0.5,
    square_idx = as.numeric(as.factor(square_id))
    
  )

# ------------------------------------------------------------
# 3) Identify which species already have fitted models
# ------------------------------------------------------------
model_files <- list.files(models_dir, pattern = "\\.rds$", full.names = TRUE)
if (length(model_files) == 0) stop("No fitted models found in: ", models_dir)

# Map model file back to species english names using your sp_filename() convention
model_index <- tibble(
  model_path = model_files,
  sp_file = tools::file_path_sans_ext(basename(model_files))
)

# Join to your species table to get species_id, english_name
# (Assumes sp_filename(english_name) matches what 09a wrote to disk.)
sp_lookup <- species_to_model %>%
  mutate(sp_file = sp_filename(english_name)) %>%
  select(species_id, english_name, sp_file)

model_index <- model_index %>%
  left_join(sp_lookup, by = "sp_file") %>%
  filter(!is.na(species_id))

model_index
message("Models found (matched to species_to_model): ", nrow(model_index))

# ------------------------------------------------------------
# 4-5) Loop: load fitted model, simulate new observations at survey locations
# ------------------------------------------------------------

for (k in seq_len(nrow(model_index))) {
  
  sp_english <- model_index$english_name[k]
  sp_code    <- as.character(model_index$species_id[k])
  model_path <- model_index$model_path[k]
  
  message("\n--- Simulating: ", sp_english, " (", sp_code, ")")
  
  # Load fitted joint model (from empirical data)
  mod <- readRDS(model_path)
  
  # Build species-specific covariate spec
  covars_sp <- intersect(base_covars, names(surveys_sim_design))
  cov_df_sp <- make_cov_df(covars_sp)
  pred_formula <- make_pred_formula3(cov_df_sp, include_kappa = TRUE, include_aru = TRUE)
  
  # Build newdata for simulation:
  #   keep real design (timing, effort, ARU indicator, etc.)
  #   BUT force Atlas3=0 to simulate "no change truth"
  sim_newdata <- surveys_sim_design %>%
    mutate(Atlas3 = 0,  # <-- critical no-change truth
           Atlas3_c = Atlas3 - 0.5)
  
  # Get posterior draws of eta at survey events
  # preds$eta is n_obs x n_draw
  preds <- predict_inla(
    mod = mod,
    grid = sim_newdata,
    pred_formula = pred_formula,
    n.samples = n_draws_sim,
    seed = seed_base + k
  )
  
  eta <- preds$eta
  mu  <- exp(eta)
  
  # Simulate counts (Poisson here; match your family)
  # If n_draws_sim == 1, eta is (n_obs x 1), so this yields one dataset.
  y_sim <- matrix(rpois(length(mu), lambda = as.vector(mu)),
                  nrow = nrow(mu), ncol = ncol(mu))
  
  # Choose which simulated draw to store as "the" dataset (if multiple)
  draw_to_keep <- 1L
  count_sim <- as.integer(y_sim[, draw_to_keep])
  
  # Create per-species simulated dataset:
  #   - keep Atlas label (design)
  #   - restore Atlas3 design indicator for later refit
  sp_dat_sim <- surveys_sim_design %>%
    mutate(
      count = count_sim,
      Atlas3 = if_else(Atlas == "OBBA2", 0, 1)
    )
  
  ggplot()+
    geom_sf(data = sp_dat_sim, col = "gray90")+
    geom_sf(data = sp_dat_sim %>% subset(count>0), col = "black", size = 0.1)+
    theme_bw()+
    facet_wrap(.~Atlas)
  
  # ------------------------------------------------------------
  # 6: fit joint atlas model to simulated data
  # ------------------------------------------------------------
  
  # --- Covariates spec (per species)
  covars_present <- intersect(base_covars, names(sp_dat_sim))
  cov_df_sp <- make_cov_df(covars_present)
  
  timeout_min <- 15
  
  # Covariates that may be included (subset to those present for each species)
  base_covars <- c(
    "on_river_N",
    "on_river_S",
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
  
  # Spatial field priors
  prior_range_abund  <- c(250, 0.50) # 50% chance spatial autocorrelation is less than 250 km
  prior_sigma_abund  <- c(3, 0.1)
  prior_range_change <- c(250, 0.50) # 10% chance spatial autocorrelation is less than 250 km
  prior_sigma_change <- c(0.1, 0.1)  # 10% chance SD is larger than 0.1
  
  # Fit model
  mod_refit <- NULL
  mod_refit <- try(
    fit_inla_multi_atlas3(
      sp_dat = sp_dat_sim,
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
  
  # ------------------------------------------------------------
  # Generate predictions across full study area
  # ------------------------------------------------------------
  
  
  
  
  message("Saved simulated dataset: ", sim_path)
}

message("\nSimulation complete.")