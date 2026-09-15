# ============================================================
# 07_fit_models_and_predict.R
#
# Fit joint OBBA2/OBBA3 INLA/inlabru models for selected species and generate
# prediction products on the full 1-km grids. The `survey_set` switch selects a
# mode from mode_settings; everything else is shared.
#
# Modes differ only in:
#   model_name    -> output folder names
#   survey_types  -> which survey types are kept per species
#   fit_function  -> fit_PC_ARU()    (point counts + ARUs; site + square iid)
#                    fit_PC_ARU_CL() (adds BBA + linear-transect checklists)
#   include_site / include_square -> the two structural iid switches (independent)
#
# Error family is chosen per species from the pooled count distribution -- see
# choose_stream_family() and the loop block in 6.3 -- and ONE family is shared
# across all of a species' surveys. Negative binomial is used only where the
# non-zero counts have a real tail; otherwise Poisson, stable on near-binary data.
#
# Prediction scale: the prediction formula omits the iid terms (kappa_diff,
# site_iid; both fitted constr = TRUE), so raw exp(eta) is the GEOMETRIC mean over
# squares/sites. An optional lognormal variance correction (vc_* below, off by
# default) rescales to the arithmetic mean. Predictions are standardized to
# optimal_TOD and days_midpoint = 0 -- a standardized index.
#
# Outputs (model_output/):
#   predictions_<model_name>/<species>_1km.rds
#   summaries_<model_name>/model_summaries.rds
#   data_used_<model_name>/<species>_1km.rds
# ============================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(sf)
  library(INLA)
  library(inlabru)
  library(fmesher)
  library(here)
  library(ggplot2)
})

source(here::here("R", "00_config_paths.R"))
source(file.path(paths$functions, "inla_model_utils_SS.R"))
source(file.path(paths$functions, "survey_processing_utils.R"))

# ============================================================
# Choose which model to run (normally the only lines you edit)
# ============================================================

survey_set         <- "PC_ARU_CL_SS"

# ============================================================
# Additional settings that could be modified
# ============================================================

rerun_predictions  <- FALSE
n_prediction_draws <- 500
prediction_seed    <- 0

# Modelability thresholds (point-count/ARU detections and squares, per atlas).
min_detections <- 100
min_squares    <- 20

# Special (targeted single-protocol) surveys. Used only by the PC_ARU_CL_SS mode
# below: each type gets its OWN intercept and NO effort term, and its intercept
# is fit only where the species is detected in >= its square threshold (types
# below threshold have their rows dropped). See fit_PC_ARU_CL() for the mechanics.
special_survey_types <- c(
  "Eastern Screech-Owl Survey", "Long-eared Owl Survey", "Marshbird Survey",
  "Nightjar Survey Protocol",   "Central Ontario Owls",  "Northern Hawk Owl Survey"
)
# One threshold (atlas squares with a detection) applied to every special type,
# or a named vector to override per type, e.g.
#   c(default = 10, "Marshbird Survey" = 20, "Nightjar Survey Protocol" = 15)
min_special_squares <- 10

# Optional: restrict to a subset for testing (NULL = run all).
species_test <- c("Eastern Whip-poor-will","Virginia Rail")

# Each mode points at a base fit function and carries the two structural switches
# explicitly, so any site/square combination is reachable without a wrapper.
# include_site / include_square drive EVERY output folder via model_name.
# Both fit functions take a single `family` arg: one shared error family spans all
# of a species' surveys. The family itself is chosen per species in the loop (6.3).
mode_settings <- list(
  
  # ---- Point counts + ARUs (1, 3, 5 min) ----
  PC_ARU = list(
    model_name       = "PC_ARU",
    fit_function     = fit_PC_ARU,
    survey_types     = c("Point_Count", "ARU"),
    survey_durations = c(1, 3, 5),
    include_site     = TRUE,
    include_square   = TRUE
  ),
  
  # ---- Point counts + ARUs + checklists (BBA stationary + linear transects) ----
  # OBBA3-only. PC/ARU/BBA and linear transects share ONE observation model with a
  # single shared error family across all surveys; linear transects keep their own
  # effort terms and site-year iid (snap_m_lt) but not a separate family.
  PC_ARU_CL = list(
    model_name       = "PC_ARU_CL",
    fit_function     = fit_PC_ARU_CL,
    survey_types     = c("Point_Count", "ARU",
                         "Breeding Bird Atlas", "Linear transect"),
    survey_durations = c(1, 3, 5),
    include_site     = TRUE,
    include_square   = TRUE
  ),

  # ---- PC + ARU + checklists + special single-protocol surveys ----
  # As PC_ARU_CL, plus the targeted owl/nightjar/marshbird protocols. Each special
  # type keeps its own intercept (a protocol detectability offset) but NO effort
  # term. A type's intercept is fit only where the species is detected in >= its
  # square threshold; types below threshold have their rows dropped inside the fit
  # function. A species that qualifies for NO special survey is skipped entirely
  # (see the loop), since its model would just duplicate PC_ARU_CL.
  #
  # special_share_tod / special_share_doy control, PER TYPE, whether a special
  # survey shares the diurnal TOD / seasonal DOY smooths. FALSE keeps every
  # special type out (timing/season absorbed by its intercept) -- the safe default
  # given the deep-night owl and dawn+dusk marshbird protocols. To let only the
  # near-dawn hawk-owl survey share TOD (it sits inside the PC/ARU diurnal
  # window), pass e.g.
  #   special_share_tod = c("Northern Hawk Owl Survey" = TRUE)
  # Shares fit_PC_ARU_CL(): the special machinery is inert unless
  # special_survey_types is non-empty.
  PC_ARU_CL_SS = list(
    model_name       = "PC_ARU_CL_SS",
    fit_function     = fit_PC_ARU_CL,
    survey_types     = c("Point_Count", "ARU",
                         "Breeding Bird Atlas", "Linear transect",
                         special_survey_types),
    survey_durations = c(1, 3, 5),
    include_site     = TRUE,
    include_square   = TRUE,
    special_survey_types = special_survey_types,
    min_special_squares  = min_special_squares,
    special_share_tod    = FALSE,
    special_share_doy    = FALSE
  )
)


stopifnot(survey_set %in% names(mode_settings))
mode       <- mode_settings[[survey_set]]
model_name <- mode$model_name
fit_model  <- mode$fit_function

message("Mode: ", survey_set, "  ->  model_name = ", model_name)

# ============================================================
# 1. Global configuration
# ============================================================

# Per-stream negative-binomial gate (see choose_stream_family()). A stream is fit
# as nbinomial only when it has >= nb_min_positive positive counts AND >=
# nb_min_tail_n of them >= nb_tail_value -- i.e. a genuine non-zero tail with which
# to estimate overdispersion. Otherwise Poisson (the NB size is unidentified on
# near-binary counts and destabilises the fit).

nb_min_positive <- 50
nb_tail_value   <- 5
nb_min_tail_n   <- 10

# Decide one stream's error family from its count distribution.
choose_stream_family <- function(count,
                                 min_positive = 50,
                                 tail_value   = 10,
                                 min_tail_n   = 10) {
  y_pos <- count[is.finite(count) & count > 0]
  if (length(y_pos) >= min_positive && sum(y_pos >= tail_value) >= min_tail_n) {
    "nbinomial"
  } else {
    "poisson"
  }
}

# Candidate fixed-effect covariates. Species-specific filtering below drops any
# absent from the data or with no variation after filtering.
base_covars <- c(
  
  "ForestNeedleleaf", "ForestBroadleaf", "ForestMixed", "Cropland",
  "Urban", "On_Road",
  "Lake_Lg", "Lake_Sm", "GreatLakes", "HudsonBayCoast", "River_Lg", "River_Sm",
  
  # Modeled with separate north/south effects.
  "Grassland_South", "Grassland_North",
  "Shrubland_South", "Shrubland_North",
  "Wetland_South",   "Wetland_North"
)

priors_list <- list(
  # Large-scale abundance & change fields (SPDE, PC-Matern). 
  # Range prior favours long (300-1000 km) ranges
  prior_range_abund  = c(300, 0.1),   # P(range < 300 km) = 0.1
  prior_sigma_abund  = c(1, 0.1),     # P(sigma > 1)    = 0.1
  prior_range_change = c(300, 0.1),   # P(range < 300 km) = 0.1
  prior_sigma_change = c(0.3, 0.1),   # P(sigma > 0.3)    = 0.1
  
  # iid sinks: atlas-square (~10 km) and site (repeated visits).
  kappa_pcprec_diff  = c(1, 0.1),     # P(sigma_square > 1) = 0.1
  prior_site_pcprec  = c(1, 0.1),     # P(sigma_site  > 1) = 0.1
  
  # Time-of-day smooth: FIXED range + estimated sigma.
  fixed_TOD_range = 4.5,
  prior_TOD_sigma = c(1, 0.1),        # P(sigma > 1) = 0.1 
  
  # Day-of-year smooth: FIXED range + estimated sigma.
  fixed_DOY_range = 30,
  prior_DOY_sigma = c(1, 0.1),        # P(sigma > 1) = 0.1
  
  # Checklist effort slopes on log(effort / median effort). A single log-linear
  # coefficient per protocol (BBA duration, LT distance), NOT a 1-D SPDE. Prior
  # mean > 0 favours more detections with more effort (diminishing returns).
  BBA_log_duration_prior_mean = 0.30,
  BBA_log_duration_prior_sd   = 0.20,
  LT_log_distance_prior_mean  = 0.30,
  LT_log_distance_prior_sd    = 0.20
)

# INLA approximation settings used inside the fit function.
int_strategy <- "ccd"
strategy     <- "simplified.laplace"

# ============================================================
# 2. Input/output locations
# ============================================================

in_file     <- file.path(paths$data_clean, "birds", "data_ready_for_analysis.rds")
species_dir <- file.path(paths$data_clean, "birds", "species_data")   # written by 06

if (!file.exists(in_file)) {
  stop("Cannot find input at: ", in_file,
       "\nHave you run 06_filter_and_finalize_surveys.R?")
}

out_dir       <- paths$model_output
pred_dir      <- file.path(out_dir, paste0("predictions_", model_name))
summary_dir   <- file.path(out_dir, paste0("summaries_", model_name))
data_used_dir <- file.path(out_dir, paste0("data_used_", model_name))

purrr::walk(c(out_dir, pred_dir, summary_dir, data_used_dir),
            dir.create, recursive = TRUE, showWarnings = FALSE)

model_summaries_path <- file.path(summary_dir, "model_summaries.rds")
model_summaries      <- load_or_empty_list(model_summaries_path)

# ============================================================
# 3. Load finalized data
# ============================================================

dat <- readRDS(in_file)

all_surveys                 <- dat$all_surveys
grid_OBBA2                  <- dat$grid_OBBA2
grid_OBBA3                  <- dat$grid_OBBA3
study_boundary              <- dat$study_boundary %>% sf::st_as_sf()
species_detection_summaries <- dat$species_detection_summaries
checklist_candidates        <- dat$checklist_candidates

if (!dir.exists(species_dir)) {
  stop("Per-species data directory not found: ", species_dir,
       "\nRerun 06_filter_and_finalize_surveys.R.")
}

# ============================================================
# 4. List of species to consider
# ============================================================

species_all <- species_detection_summaries %>%
  dplyr::select(sp_english, species_id) %>%
  distinct()

# Apply the test filter to both sets when supplied (NULL means "all").
if (!is.null(species_test)) {
  species_all <- species_all %>% filter(sp_english %in% species_test)
}

message("Species in list: ", nrow(species_all))

# ============================================================
# 5. Build spatial meshes once
# ============================================================

bndry <- study_boundary %>%
  st_make_valid() %>%
  st_cast("POLYGON", warn = FALSE) %>%
  mutate(area_km2 = as.numeric(st_area(geometry))) %>%
  filter(area_km2 >= 25) %>%
  summarise(geometry = st_union(geometry)) %>%
  st_simplify(dTolerance = 10, preserveTopology = TRUE)

hull <- fmesher::fm_extensions(bndry, convex = c(50, 600), concave = c(50, 600))

mesh_abund <- fmesher::fm_mesh_2d_inla(
  loc = sf::st_as_sfc(all_surveys), boundary = hull,
  max.edge = c(30, 150), cutoff = 30, crs = sf::st_crs(all_surveys)
)

mesh_chg <- fmesher::fm_mesh_2d_inla(
  loc = sf::st_as_sfc(all_surveys), boundary = hull,
  max.edge = c(30, 150), cutoff = 30, crs = sf::st_crs(all_surveys)
)

# ---- 1-D time-of-day mesh ----
# Built ONCE so every species shares one basis and the fixed range means the same
# thing for all of them. Knot spacing must be <= fixed_TOD_range / 5 or the basis
# (not the range) sets smoothness. Domain is cyclic (Hours_After_Reference is
# measured from 3 h before sunrise).
TOD_knot_spacing <- priors_list$fixed_TOD_range / 6
mesh_TOD <- make_TOD_mesh(span = c(0, 24), knot_spacing = TOD_knot_spacing,
                          pad_range = priors_list$fixed_TOD_range)

message("TOD mesh: ", length(mesh_TOD$loc), " knots at ",
        signif(TOD_knot_spacing, 3), " h spacing; range fixed at ",
        priors_list$fixed_TOD_range, " h")

# Survey coverage by hour and protocol (informs the reference-hour logic in 6.6).
TOD_coverage <- all_surveys %>%
  sf::st_drop_geometry() %>%
  dplyr::mutate(hour = floor(Hours_After_Reference)) %>%
  dplyr::count(hour, Survey_Type) %>%
  tidyr::pivot_wider(names_from = Survey_Type, values_from = n, values_fill = 0)
print(TOD_coverage, n = 25)

# ============================================================
# 6. Main species loop
# ============================================================

for (i in seq_len(nrow(species_all))) {
  
  # ---- 6.1 Identifiers, paths, modelability ----
  sp_name <- species_all$sp_english[i]
  sp_code <- as.character(species_all$species_id[i])
  sp_file <- sp_filename(sp_name)
  
  message("\n====================\n",
          i, "/", nrow(species_all), ": ", sp_name,
          " (species_id = ", sp_code, ")\n====================")
  
  pred_path <- file.path(pred_dir, paste0(sp_file, "_1km.rds"))
  
  if (file.exists(pred_path) && !rerun_predictions) {
    message("Predictions already exist for ", sp_name, "; skipping model fit.")
    next
  }
  
  dat_path <- file.path(data_used_dir, paste0(sp_file, "_1km.rds"))
  sp_path  <- sp_data_path(species_dir, sp_name)
  
  if (!file.exists(sp_path)) {
    message("Skipping; no per-species data file from script 06.")
    next
  }
  
  # ---- 6.2 Load species data (safe-date filtering + pred_doy from 06) ----
  sp <- load_sp_dat(sp_path, all_surveys,
                    survey_types = mode$survey_types,
                    survey_durations = mode$survey_durations)
  
  sp_dat        <- sp$sp_dat
  sp_safe_dates <- sp$sp_safe_dates_unmodified
  pred_doy      <- sp$pred_doy
  
  if (nrow(sp_dat) == 0) {
    message("Skipping ", sp_name, "; no surveys of the required types remain.")
    next
  }
  
  n_det <- sum(sp_dat$count > 0)
  if (n_det == 0) {
    message("Skipping ", sp_name, "; no positive counts remain after filtering.")
    next
  }
  
  # ---- 6.2b Restrict to the special surveys actually analysed ----
  # load_sp_dat() returns every type in mode$survey_types, but fit_PC_ARU_CL()
  # keeps a special type's intercept only where the species is detected in >= its
  # square threshold and DROPS all rows of the types below it (fit section 3a).
  # Mirror that here on the SAME resolver so sp_dat holds only the rows the model
  # will see BEFORE it feeds the family choice, sp_det_summary, is_modelable, the
  # saved survey_counts and the fit -- otherwise the dat_path record (used by the
  # downstream "analysed data" plots) overstates the data. Also folds in the old
  # "fit only if some special type qualifies" gate.
  special_keep   <- character(0)
  special_report <- list()
  if (length(mode$special_survey_types) > 0) {
    if (!"square_id" %in% names(sp_dat)) {
      stop("square_id required to threshold special surveys for ", sp_name, ".")
    }
    special_present <- intersect(mode$special_survey_types, unique(sp_dat$Survey_Type))
    for (type in special_present) {
      det  <- sp_dat$Survey_Type == type & sp_dat$count > 0
      n_sq <- dplyr::n_distinct(sp_dat$square_id[det])
      thr  <- resolve_special_threshold(type, min_special_squares = mode$min_special_squares)
      keep <- n_sq >= thr
      if (keep) special_keep <- c(special_keep, type)
      special_report[[type]] <- c(n_det_squares = n_sq, threshold = thr,
                                  kept = as.numeric(keep))
    }
    special_drop <- setdiff(special_present, special_keep)

    if (length(special_keep) == 0) {
      message("Skipping ", sp_name, "; no special survey meets its square ",
              "threshold (mode ", survey_set, ").")
      next
    }
    if (length(special_drop) > 0) {
      message("  dropping below-threshold special rows: ",
              paste(special_drop, collapse = ", "))
      sp_dat <- sp_dat %>% filter(!Survey_Type %in% special_drop)
      n_det  <- sum(sp_dat$count > 0)   # recompute on the analysed row set
    }
  }
  
  
  # ---- 6.3 Error family from the count distribution ----
  # ONE error family is shared across all of a species' surveys, so it is chosen
  # from the pooled count distribution (every retained row, LT included). NB
  # overdispersion is estimable only where the non-zero counts have a real tail
  # (choose_stream_family); near-binary or sparse data fall back to Poisson.
  shared_family <- choose_stream_family(
    sp_dat$count,
    min_positive = nb_min_positive,
    tail_value   = nb_tail_value,
    min_tail_n   = nb_min_tail_n
  )
  
  family_args  <- list(family = shared_family)
  stat_arg     <- "family"
  error_family <- shared_family
  
  y_pos <- sp_dat$count[sp_dat$count > 0]
  count_dist_summary <- tibble::tibble(
    n         = nrow(sp_dat),
    n_pos     = length(y_pos),
    det_rate  = round(mean(sp_dat$count > 0), 4),
    max_count = if (length(y_pos)) max(y_pos) else 0L,
    n_ge2     = sum(y_pos >= 2),
    n_tail    = sum(y_pos >= nb_tail_value),   # counts >= nb_tail_value
    family    = shared_family
  )
  
  # Determine if species clears sample size threshold for modeling
  sp_det_summary <- sp_dat %>%
    as.data.frame() %>%
    subset(count>0) %>%
    
    group_by(Atlas) %>%
    summarize(n_squares    = length(unique(square_id)),
              n_detections = sum(count>0),
              n_detections_PC = sum(Survey_Type == "Point_Count"),
              n_detections_ARU = sum(Survey_Type == "ARU"),
              n_detections_BBA = sum(Survey_Type == "Breeding Bird Atlas"),
              n_detections_LT = sum(Survey_Type == "Linear transect"))
  
  # Special-survey detection squares (pooled across atlases), to audit / tune the
  # per-type min_special_squares thresholds. Matches the count fit_PC_ARU_CL()
  # thresholds on: unique square_id with count > 0 for each special type. Empty
  # for modes that don't load the special surveys.
  sp_special_summary <- sp_dat %>%
    st_drop_geometry() %>%
    filter(Survey_Type %in% special_survey_types, count > 0) %>%
    group_by(Survey_Type) %>%
    summarize(n_det_squares = n_distinct(square_id),
              n_detections  = n(),
              .groups       = "drop")
  
  # Special-survey keep/drop is resolved in 6.2b; below-threshold rows (and any
  # species that qualifies for no special type) are already handled there.
  
  if (max(sp_det_summary$n_squares) > min_squares & max(sp_det_summary$n_detections) > min_detections){
    is_modelable <- TRUE
  } else{
    is_modelable <- FALSE
  }
  
  # Site-level REs: 500-m site-year grid.
  sp_dat <- add_site_ids(sp_dat, snap_m = 500, tolerance_m = 0)
  
  # ---- 6.4 Save data record for every species ----
  save_atomic(
    list(
      sp_english    = sp_name,
      sp_code       = sp_code,
      is_modelable  = is_modelable,
      sp_safe_dates = sp_safe_dates,
      pred_doy      = pred_doy,
      survey_counts = sp_dat %>% st_drop_geometry() %>% select(survey_id, count, site_id),
      sp_det_summary = sp_det_summary,
      sp_special_summary = sp_special_summary,
      special_report = special_report,     # per-type n_det_squares / threshold / kept
      error_family  = error_family
    ),
    dat_path
  )
  
  if (!is_modelable) {
    message("Honeycomb data saved for ", sp_name, "; not enough data to fit a model.")
    next
  }
  
  
  # OPTIONAL: THIN TO A SINGLE OBSERVATION PER SITE
  # sp_dat <- thin_one_per_site(sp_dat, group_cols = "site_id")
  
  # ---- 6.5 Species-specific covariate table ----
  # Keep candidate covariates present in this species' data with > 1 unique finite
  # value after filtering. sd_linear kept tight (0.5): loose habitat betas on
  # sparse 0/1 detections drive checkerboard artefacts in the abundance maps.
  covars_present <- intersect(base_covars, names(sp_dat))
  if (length(covars_present) > 0) {
    covars_present <- sp_dat %>%
      sf::st_drop_geometry() %>%
      select(all_of(covars_present)) %>%
      summarise(across(everything(),
                       ~ dplyr::n_distinct(.x[is.finite(.x)], na.rm = TRUE))) %>%
      tidyr::pivot_longer(everything(), names_to = "covariate", values_to = "n_unique") %>%
      filter(n_unique > 1) %>%
      pull(covariate)
  }
  cov_df_sp <- make_cov_df(covars_present, mean = 0, sd_linear = 0.5)
  
  # ---- 6.6 Candidate reference hours + empirical starting value ----
  obs_dat <- sp_dat %>%
    st_drop_geometry() %>%
    filter(Survey_Type %in% c("Point_Count", "ARU")) %>%
    mutate(hour = round(Hours_After_Reference)) %>%
    group_by(hour, site_id) %>%
    summarize(mean_count = mean(count > 0), n = n(), .groups = "drop") %>%
    group_by(hour) %>%
    summarize(mean_count = mean(mean_count), n = n(), .groups = "drop") %>%
    filter(n > 1000)
  
  hours_supported <- obs_dat$hour
  
  if (nrow(obs_dat) == 0) {
    message("  too few PC/ARU surveys to set a reference hour; defaulting to 3 (sunrise)")
    optimal_TOD_empirical <- 3
    hours_supported       <- 3
  } else {
    optimal_TOD_empirical <- obs_dat$hour[which.max(obs_dat$mean_count)]
    plot(mean_count ~ hour, data = obs_dat, type = "p", pch = 19, lwd = 2,
         main = paste0(sp_name, " - detection rate by time of day\n\nEmpirical peak = ",
                       round(optimal_TOD_empirical, 2), " hours since reference"))
    abline(v = optimal_TOD_empirical, col = "dodgerblue", lwd = 2)
    abline(v = 3, col = "gray80", lwd = 2, lty = 2)
  }
  
  # ---- 1-D day-of-year mesh (linear domain, per species) ----
  DOY_range        <- range(sp_dat$days_midpoint)
  DOY_knot_spacing <- priors_list$fixed_DOY_range / 6
  mesh_DOY <- make_DOY_mesh(span = c(DOY_range[1] - 7, DOY_range[2] + 7),
                            knot_spacing = DOY_knot_spacing,
                            pad_range = priors_list$fixed_DOY_range)
  message("DOY mesh: ", length(mesh_DOY$loc), " knots at ",
          signif(DOY_knot_spacing, 3), " d spacing; range fixed at ",
          priors_list$fixed_DOY_range, " d")
  
  # ---- 6.7 Fit model ----
  start_model <- Sys.time()
  
  # Shared arguments; the family arg(s) come from the per-stream decision (6.3)
  # and are appended via family_args in the do.call below.
  fit_args <- list(
    sp_dat     = sp_dat,
    mesh_abund = mesh_abund,
    mesh_chg   = mesh_chg,
    mesh_TOD   = mesh_TOD,
    covariates = cov_df_sp,
    
    prior_range_abund  = priors_list$prior_range_abund,
    prior_sigma_abund  = priors_list$prior_sigma_abund,
    prior_range_change = priors_list$prior_range_change,
    prior_sigma_change = priors_list$prior_sigma_change,
    
    fixed_TOD_range = priors_list$fixed_TOD_range,
    prior_TOD_sigma = priors_list$prior_TOD_sigma,
    
    mesh_DOY        = mesh_DOY,
    fixed_DOY_range = priors_list$fixed_DOY_range,
    prior_DOY_sigma = priors_list$prior_DOY_sigma,
    
    kappa_pcprec_diff = priors_list$kappa_pcprec_diff,
    prior_site_pcprec = priors_list$prior_site_pcprec,
    
    int_strategy = int_strategy,
    strategy     = strategy,
    
    include_site   = mode$include_site,
    include_square = mode$include_square
  )
  
  # Mode-specific extras: checklist effort priors (fit_PC_ARU_CL only) and the
  # special-survey settings (PC_ARU_CL_SS only). Appended only when the selected
  # fit function DECLARES the argument, so PC/ARU-only modes don't error on unused
  # ones; NULL entries (fields a mode omits) are dropped so absent modes keep the
  # fit-function defaults.
  mode_extra_args <- list(
    BBA_log_duration_prior_mean = priors_list$BBA_log_duration_prior_mean,
    BBA_log_duration_prior_sd   = priors_list$BBA_log_duration_prior_sd,
    LT_log_distance_prior_mean  = priors_list$LT_log_distance_prior_mean,
    LT_log_distance_prior_sd    = priors_list$LT_log_distance_prior_sd,
    
    special_survey_types = mode$special_survey_types,
    min_special_squares  = mode$min_special_squares,
    special_share_tod    = mode$special_share_tod,
    special_share_doy    = mode$special_share_doy
  )
  
  mode_extra_args <- mode_extra_args[!vapply(mode_extra_args, is.null, logical(1))]
  fit_args <- c(
    fit_args,
    mode_extra_args[names(mode_extra_args) %in% names(formals(fit_model))]
  )
  
  mod <- try(do.call(fit_model, c(fit_args, family_args)), silent = TRUE)
  
  # Fallback: a chosen nbinomial stationary family that fails to fit -> retry once
  # with Poisson (rare, since nbinomial is only chosen when the tail supports it).
  if ((inherits(mod, "try-error") || is.null(mod)) &&
      identical(family_args[[stat_arg]], "nbinomial")) {
    message("  stationary nbinomial fit failed for ", sp_name, "; retrying with Poisson.")
    print(mod)
    family_args[[stat_arg]] <- "poisson"
    mod <- try(do.call(fit_model, c(fit_args, family_args)), silent = TRUE)
  }
  
  if (inherits(mod, "try-error") || is.null(mod)) {
    message("Model failed for ", sp_name, "; skipping this species.")
    print(mod)
    next
  }
  
  # Stationary family actually used (after any fallback).
  stationary_family_used <- family_args[[stat_arg]]
  family_fallback        <- !identical(stationary_family_used, shared_family)
  if (family_fallback) {
    message("  -> fit succeeded with Poisson stationary error (nbinomial fallback).")
  }
  
  end_model   <- Sys.time()
  fit_minutes <- round(as.numeric(end_model - start_model, units = "mins"), 1)
  
  print(summary(mod))
  message("\n====================\n", i, "/", nrow(species_all), ": ", sp_name,
          " (species_id = ", sp_code, "); ", fit_minutes, " min to fit model\n====================")
  
  
  # ---- 6.7b Reference time of day, from the fitted TOD curve ----
  # Read off TOD_global (estimated alongside the protocol intercepts), restricted
  # to hours_supported so the reference is never placed where the curve
  # extrapolates. Falls back to the empirical peak if unavailable.
  optimal_TOD_fitted <- optimal_TOD_from_fit(mod, hours_allowed = hours_supported)
  
  grid  <- seq(min(hours_supported), max(hours_supported), by = 0.05)
  curve <- TOD_curve_from_fit(mod, grid_hours = grid)
  plot(pred ~ Hours_After_Reference, data = curve, type = "l", lwd = 2,
       xlab = "Hours after reference", ylab = "TOD_global (log scale)")
  rect(hours_supported - 0.5, par("usr")[3], hours_supported + 0.5, par("usr")[4],
       col = adjustcolor("dodgerblue", 0.3), border = NA)
  abline(h = 0, col = "grey70")
  
  # ---- 6.7c Fitted day-of-year curve, over a survey-density histogram ----
  # The histogram shows where days_midpoint actually carries data. If the dip in
  # the DOY curve lines up with a density trough at day 0, it is a centring/
  # constraint artifact of pred_doy, not phenology.
  DOY_BIN_DAYS <- 3   # histogram bin width (days)
  
  doy_range <- range(sp_dat$days_midpoint, na.rm = TRUE)
  doy_grid  <- seq(doy_range[1], doy_range[2], by = 0.5)
  doy_curve <- DOY_curve_from_fit(mod, grid_days = doy_grid)
  
  if (!is.null(doy_curve)) {
    breaks <- seq(floor(doy_range[1]), ceiling(doy_range[2]) + DOY_BIN_DAYS,
                  by = DOY_BIN_DAYS)
    h_srv  <- hist(sp_dat$days_midpoint, breaks = breaks, plot = FALSE)
    
    op <- par(mar = c(5, 4, 4, 4) + 0.1)   # room for a right-hand axis
    
    # Histogram first (left axis). Inflate ylim so the bars fill only the lower
    # ~third of the panel and sit UNDER the curve rather than over it.
    plot(h_srv, col = "grey88", border = "grey70",
         xlim = doy_range, ylim = c(0, max(h_srv$counts) * 3),
         main = paste0(sp_name, " - day-of-year effect"),
         xlab = "Days from safe-date midpoint", ylab = "n surveys")
    
    # Curve on top, on its own (right) axis.
    par(new = TRUE)
    plot(pred ~ days_midpoint, data = doy_curve, type = "l", lwd = 2,
         xlim = doy_range, axes = FALSE, xlab = "", ylab = "")
    axis(4)
    mtext("DOY_global (log scale)", side = 4, line = 2.5)
    
    abline(h = 0, col = "grey70")                       # curve scale (right axis)
    abline(v = 0, col = "dodgerblue", lwd = 2, lty = 2) # day 0 = pred_doy
    
    par(op)
  }
  
  if (is.na(optimal_TOD_fitted)) {
    optimal_TOD     <- optimal_TOD_empirical
    optimal_TOD_src <- "empirical (fitted curve unavailable)"
  } else {
    optimal_TOD     <- optimal_TOD_fitted
    optimal_TOD_src <- "fitted TOD curve"
  }
  
  message("  reference hour = ", round(optimal_TOD, 2),
          " (", optimal_TOD_src, "; empirical peak was ",
          round(optimal_TOD_empirical, 2), ")")
  
  model_summaries[[sp_name]] <- list(
    sp_name          = sp_name,
    sp_code          = sp_code,
    error_family           = error_family,             # chosen error family
    stationary_family_used = stationary_family_used,   # after any fallback
    family_fallback        = family_fallback,
    priors           = priors_list,
    int_strategy     = int_strategy,
    strategy         = strategy,
    n_surveys        = nrow(sp_dat),
    n_detections     = n_det,
    n_covariates     = nrow(cov_df_sp),
    covariates       = cov_df_sp,
    pred_doy         = pred_doy,
    fitted_TOD_curve = curve,      # no uncertainty
    fitted_DOY_curve = doy_curve,  # no uncertainty
    
    optimal_TOD           = optimal_TOD,
    optimal_TOD_empirical = optimal_TOD_empirical,
    optimal_TOD_fitted    = optimal_TOD_fitted,
    optimal_TOD_source    = optimal_TOD_src,
    hours_supported       = hours_supported,
    tod_fixed_range       = mod$tod_fixed_range,
    tod_knot_gap          = mod$tod_knot_gap,
    fit_minutes           = fit_minutes,
    summary_fixed         = mod$summary.fixed,
    summary_hyperpar      = mod$summary.hyperpar
  )
  save_atomic(model_summaries, model_summaries_path)
  
  # ---- 6.8 Full-grid predictions ----
  # Standardized to optimal_TOD, days_midpoint = 0, with kappa_diff and site_iid
  # omitted from pred_formula -> exp(eta) is a geometric mean (6.9 optionally
  # rescales to the arithmetic mean).
  start_prediction <- Sys.time()
  message("Generating predictions for full 1-km grid")
  
  pred_formula <- make_pred_formula_multiatlas(cov_df_sp)
  
  pred_grid <- make_pred_grid(grid_OBBA2, grid_OBBA3) %>%
    mutate(Hours_After_Reference = optimal_TOD, days_midpoint = 0)
  
  preds <- predict_all_pixels(
    mod = mod, pred_grid = pred_grid, pred_formula = pred_formula,
    n.samples = n_prediction_draws, seed = prediction_seed
  )
  
  # ---- 6.9 Optional lognormal variance correction ----
  # Note this is a carryover from an earlier version; preserved in case it is of use in the future
  # Script 09 will incorporate variance corrections for mapping
  vc_terms        <- character(0)   # Disabled by default
  vc_propagate    <- FALSE          # Propogate random effect error on a draw-by-draw basis
  vc_cap_quantile <- 0.99           # Remove extreme values
  
  vc <- make_lognormal_correction(
    mod = mod, terms = vc_terms, propagate = vc_propagate,
    n_draws = n_prediction_draws, cap_quantile = vc_cap_quantile
  )
  message("  ", describe_lognormal_correction(vc))
  
  preds$mu2 <- apply_lognormal_correction(preds$mu2, vc$factor)
  preds$mu3 <- apply_lognormal_correction(preds$mu3, vc$factor)
  
  # Stored without the per-draw factor vector; factor_summary keeps its quantiles.
  vc_record <- vc[setdiff(names(vc), "factor")]
  model_summaries[[sp_name]]$variance_correction <- vc_record
  save_atomic(model_summaries, model_summaries_path)
  
  # ---- 6.10 Terrestrial open-water correction ----
  preds$mu2_Corrected_for_Water <- preds$mu2 * (1 - grid_OBBA2$open_water)
  preds$mu3_Corrected_for_Water <- preds$mu3 * (1 - grid_OBBA3$open_water)
  
  # ---- 6.11 Summarize pixel predictions and hex draws ----
  pred_summary <- summarize_predictions(preds$mu2, preds$mu3)
  pred_summary_Corrected_for_Water <- summarize_predictions(
    preds$mu2_Corrected_for_Water, preds$mu3_Corrected_for_Water
  )
  
  g2 <- pred_grid %>% filter(Atlas == "OBBA2")
  g3 <- pred_grid %>% filter(Atlas == "OBBA3")
  
  preds_OBBA2_summary <- bind_cols(
    g2 %>% st_drop_geometry() %>% select(pixel_id, hex_id), pred_summary$OBBA2
  )
  preds_OBBA3_summary <- bind_cols(
    g3 %>% st_drop_geometry() %>% select(pixel_id, hex_id), pred_summary$OBBA3
  )
  preds_abs_change_summary <- bind_cols(
    g2 %>% st_drop_geometry() %>% select(pixel_id, hex_id), pred_summary$abs_change
  )
  
  hex_draws <- make_hex_draws(g2 = g2, mu2 = preds$mu2, mu3 = preds$mu3)
  
  preds_OBBA2_summary_Corrected_for_Water <- bind_cols(
    g2 %>% st_drop_geometry() %>% select(pixel_id, hex_id),
    pred_summary_Corrected_for_Water$OBBA2
  )
  preds_OBBA3_summary_Corrected_for_Water <- bind_cols(
    g3 %>% st_drop_geometry() %>% select(pixel_id, hex_id),
    pred_summary_Corrected_for_Water$OBBA3
  )
  preds_abs_change_summary_Corrected_for_Water <- bind_cols(
    g2 %>% st_drop_geometry() %>% select(pixel_id, hex_id),
    pred_summary_Corrected_for_Water$abs_change
  )
  
  hex_draws_Corrected_for_Water <- make_hex_draws(
    g2 = g2, mu2 = preds$mu2_Corrected_for_Water, mu3 = preds$mu3_Corrected_for_Water
  )
  
  end_prediction <- Sys.time()
  pred_minutes   <- round(as.numeric(end_prediction - start_prediction, units = "mins"), 1)
  
  # ---- 6.12 Observed survey coverage by atlas square ----
  sp_square_summary <- sp_dat %>%
    st_drop_geometry() %>%
    group_by(Atlas, square_id) %>%
    summarise(
      n_surveys    = n(),
      total_count  = sum(count),
      n_detections = sum(count > 0),
      BCR          = names(which.max(table(BCR))),
      .groups      = "drop"
    )
  
  # ---- 6.13 Save prediction products ----
  save_atomic(
    list(
      sp_name           = sp_name,
      sp_code           = sp_code,
      sp_safe_dates     = sp_safe_dates,
      sp_square_summary = sp_square_summary,
      
      error_family           = error_family,             # chosen error family
      stationary_family_used = stationary_family_used,   # after any fallback
      family_fallback        = family_fallback,
      priors       = priors_list,
      int_strategy = int_strategy,
      strategy     = strategy,
      
      fit_minutes  = fit_minutes,
      pred_minutes = pred_minutes,
      
      summary_fixed    = mod$summary.fixed,
      summary_hyperpar = mod$summary.hyperpar,
      
      pred_doy           = pred_doy,
      optimal_TOD        = optimal_TOD,
      prediction_seed    = prediction_seed,
      n_prediction_draws = n_prediction_draws,
      
      # Scale metadata (factor is 1 when vc_terms is empty, i.e. geometric scale).
      variance_correction = vc_record,
      
      # Predictions for terrestrial habitats.
      OBBA2      = preds_OBBA2_summary,
      OBBA3      = preds_OBBA3_summary,
      abs_change = preds_abs_change_summary,
      hex_draws  = hex_draws,
      
      # As above, additionally scaled by the non-open-water pixel fraction.
      OBBA2_Corrected_for_Water      = preds_OBBA2_summary_Corrected_for_Water,
      OBBA3_Corrected_for_Water      = preds_OBBA3_summary_Corrected_for_Water,
      abs_change_Corrected_for_Water = preds_abs_change_summary_Corrected_for_Water,
      hex_draws_Corrected_for_Water  = hex_draws_Corrected_for_Water
    ),
    pred_path
  )
  
  message("\n====================\n", i, "/", nrow(species_all), ": ", sp_name,
          " (species_id = ", sp_code, "); ", pred_minutes,
          " min to generate predictions\n====================")
  
  rm(preds, pred_grid, g2, g3, hex_draws, hex_draws_Corrected_for_Water)
  gc(verbose = FALSE)
}

message("\n07_fit_models_and_predict.R complete  (mode: ", survey_set, ").")