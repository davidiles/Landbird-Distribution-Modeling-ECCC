# ============================================================
# 07_fit_models_and_predict.R   (PC+ARU  or  PC+ARU+checklists)
#
# Purpose
#   Fit joint OBBA2/OBBA3 INLA/inlabru models for selected species and generate
#   prediction products on the full 1-km prediction grids. Choose the analysis
#   with the `survey_set` switch; everything else is shared across modes.
#   Each survey set has a *_nosite twin that fits the identical model with the
#   site-level iid random effect forced off, for A/B comparison of estimates,
#   interval widths, and fit/predict time.
#
# What differs between the two modes (and ONLY these things):
#   1. model_name    -> output folder names.
#   2. survey_types  -> which survey types are kept from the per-species data:
#        "PC_ARU"     point counts + ARUs (checklists excluded)
#        "PC_ARU_CL"  everything (NULL = no survey-type filter)
#   3. fit_function:
#        "PC_ARU"     fit_PC_ARU()     (point counts + ARUs; site + square iid)
#        "PC_ARU_CL"  fit_PC_ARU_CL()  (adds checklists + effort corrections)
#
# Per-species data
#   06 writes one lightweight record per species (survey_id + count for every
#   survey inside the safe-date window, across ALL survey types). load_sp_dat()
#   rejoins all_surveys, filters to this mode's survey types, and rebuilds
#   days_midpoint / duration_rescaled / ARU_1min / ARU_3min.
#
# Prediction scale
#   The prediction formula omits the iid terms (kappa_diff, site_iid), both
#   fitted with constr = TRUE, so raw exp(eta) is the GEOMETRIC mean over
#   squares/sites (a typical/median site). An OPTIONAL lognormal variance
#   correction rescales that to the arithmetic mean; it is disabled by default
#   (vc_terms = character(0) -> factor 1, a no-op) and controlled by the vc_*
#   config below. Predictions are otherwise standardized to optimal_TOD and
#   days_midpoint = 0 -- a standardized index, not a prediction of a typical outing.
#
# Main outputs (model_output/)
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
source(file.path(paths$functions, "inla_model_utils.R"))

# ============================================================
# CHOOSE WHICH MODEL TO RUN  (normally the only line you edit)
# ============================================================
#   "PC_ARU"            = point counts + ARUs
#   "PC_ARU_CL"         = point counts + ARUs + checklists (effort corrections)
#   "PC_ARU_nosite"     = as PC_ARU,    site-level iid random effect OFF
#   "PC_ARU_CL_nosite"  = as PC_ARU_CL, site-level iid random effect OFF

survey_set <- "PC_ARU_CL_nosite"

mode_settings <- list(
  PC_ARU = list(
    model_name       = "PC_ARU",
    fit_function     = fit_PC_ARU,
    survey_types     = c("Point_Count", "ARU"),
    survey_durations = c(5)
  ),
  PC_ARU_CL = list(
    model_name       = "PC_ARU_CL",
    fit_function     = fit_PC_ARU_CL,
    survey_types     = NULL,   # NULL = no survey-type filter
    survey_durations = NULL
  ),

  # ---- site-level iid FORCED OFF (identical models otherwise) ----
  # The *_nosite fit functions are thin wrappers that call their parent with
  # include_site = FALSE. Distinct model_name -> outputs land in their own
  # predictions_/summaries_/data_used_ folders, so the two arms never collide.
  PC_ARU_nosite = list(
    model_name       = "PC_ARU_nosite",
    fit_function     = fit_PC_ARU_nosite,
    survey_types     = c("Point_Count", "ARU"),
    survey_durations = c(5)
  ),
  PC_ARU_CL_nosite = list(
    model_name       = "PC_ARU_CL_nosite",
    fit_function     = fit_PC_ARU_CL_nosite,
    survey_types     = NULL,
    survey_durations = NULL
  )
)

stopifnot(survey_set %in% names(mode_settings))
mode       <- mode_settings[[survey_set]]
model_name <- mode$model_name
fit_model  <- mode$fit_function

message("Mode: ", survey_set, "  ->  model_name = ", model_name)

# ============================================================
# 1. Global configuration shared by both modes
# ============================================================

rerun_predictions  <- FALSE
n_prediction_draws <- 500
prediction_seed    <- 0

# ---- Optional lognormal variance correction (see "Prediction scale" above) --
# vc_terms       iid components whose variance is reinstated on the mean scale.
#                Terms absent from a species' model are skipped automatically.
#                character(0) disables the correction (factor = 1, a no-op);
#                "kappa_diff" applies only the square-level correction.
# vc_propagate   TRUE draws one factor per posterior sample so variance-component
#                uncertainty widens the mapped intervals; FALSE uses the plug-in
#                factor from the median variances. FALSE is the recommended default.
# vc_cap_quantile
#                exp(sigma^2/2) is convex with a long right tail for sparse
#                species; sampled variances are capped at this quantile. NA disables.
vc_terms        <- character(0)
vc_propagate    <- FALSE
vc_cap_quantile <- 0.99

# Modelability thresholds (point-count/ARU detections and squares, per atlas).
min_PC_detections <- 20
min_squares       <- 20

# Negative-binomial switch thresholds (see choose_error_family()).
nb_switch_count <- 50
nb_switch_min_n <- 50

# Candidate fixed-effect covariates. Species-specific filtering below removes
# covariates absent from the data or with no variation after filtering.
base_covars <- c(
  "ForestNeedleleaf", "ForestBroadleaf", "ForestMixed", "Wetland", "Cropland",
  "Urban", "On_Road",
  "Grassland_South", "Grassland_North",
  "Shrubland_South", "Shrubland_North",
  "Lake_Lg", "Lake_Sm",
  "GreatLakes", "HudsonBayCoast",
  "River_Lg_South", "River_Lg_North",
  "River_Sm_South", "River_Sm_North"
)

priors_list <- list(
  # Large-scale spatial fields (SPDEs). The PC prior favours long ranges; 300 km
  # strengthens that without hitting the constr=TRUE degeneracy ceiling (~domain/2).
  prior_range_abund = c(300, 0.1),   # P(range < 300 km) = 0.1
  prior_sigma_abund = c(1.5, 0.1),   # P(sigma > 1.5)   = 0.1  (prior median ~0.45, log)

  prior_range_change = c(300, 0.1),  # P(range < 300 km) = 0.1
  prior_sigma_change = c(0.3, 0.1),  # P(sigma > 0.3)    = 0.1  (tighter; shrinks change)

  kappa_pcprec_diff = c(0.5, 0.1),   # P(sigma_square > 0.5) = 0.1 (iid soaks up fine-scale)

  # NOTE: two further partitioning priors live OUTSIDE this list, by design:
  #   * site (100-m) iid   -> prior_site_pcprec, default c(0.5, 0.1) in the fit_*
  #     functions (== recommended value, so left at default and NOT passed below).
  #   * habitat covariate SD -> sd_linear in make_cov_df() (kept tight at 0.5)

  # TOD range is FIXED, not estimated (see fit_PC_ARU).
  fixed_TOD_range = 9,               # hours; held fixed, no prior
  prior_TOD_sigma = c(2, 0.1),       # P(sigma > 2) = 0.1 (prior median ~0.6)

  # Day-of-year smooth: FIXED range (days) + estimated sigma.
  fixed_DOY_range = 45,              # days; held fixed, no prior
  prior_DOY_sigma = c(2, 0.1)        # P(sigma > 2) = 0.1 (prior median ~0.6)
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
# 4. Select species to model
# ============================================================
#   species_all : EVERY species (honeycomb survey data is saved for all).
#   species_sel : the subset with enough data to fit a model (fitted + predicted).

# Optional: restrict to a handful of species for testing (NULL = run all).
species_test <- NULL

species_all <- species_detection_summaries %>%
  dplyr::select(sp_english, species_id) %>%
  distinct()

species_sel <- species_detection_summaries %>%
  subset(Survey_Type == "Point_Count" & n_det >= min_PC_detections & n_sq >= min_squares) %>%
  group_by(sp_english, species_id) %>%
  summarize(n_PCdet_BothAtlas = sum(n_det),
            n_PCsq_BothAtlas  = sum(n_sq), .groups = "drop")

# Apply the test filter to BOTH sets when supplied. (NULL truly means "all"; the
# previous `%in% NULL` form silently selected zero species.)
if (!is.null(species_test)) {
  species_all <- species_all %>% filter(sp_english %in% species_test)
  species_sel <- species_sel %>% filter(sp_english %in% species_test)
}

message("Species total (honeycomb data saved for all): ", nrow(species_all))
message("Species to model: ", nrow(species_sel))

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

hull <- fmesher::fm_extensions(bndry, convex = c(50, 200), concave = c(50, 200))

mesh_abund <- fmesher::fm_mesh_2d_inla(
  loc = sf::st_as_sfc(all_surveys), boundary = hull,
  max.edge = c(40, 80), cutoff = 40, crs = sf::st_crs(all_surveys)
)

mesh_chg <- fmesher::fm_mesh_2d_inla(
  loc = sf::st_as_sfc(all_surveys), boundary = hull,
  max.edge = c(60, 100), cutoff = 60, crs = sf::st_crs(all_surveys)
)

# ---- 1-D mesh for the time-of-day smooth ----
# Built ONCE so every species shares one basis and the fixed range means the same
# thing for all of them. Knot spacing must be <= fixed_TOD_range / 5 or the basis,
# not the range, sets the smoothness. The domain is cyclic (Hours_After_Reference
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

  # ---- 6.1 Species identifiers and output paths ----
  sp_name <- species_all$sp_english[i]
  sp_code <- as.character(species_all$species_id[i])
  sp_file <- sp_filename(sp_name)

  message("\n====================\n",
          i, "/", nrow(species_all), ": ", sp_name,
          " (species_id = ", sp_code, ")\n====================")

  pred_path <- file.path(pred_dir, paste0(sp_file, "_1km.rds"))
  dat_path  <- file.path(data_used_dir, paste0(sp_file, "_1km.rds"))
  sp_path   <- sp_data_path(species_dir, sp_name)

  if (!file.exists(sp_path)) {
    message("Skipping; no per-species data file from script 06.")
    next
  }

  # Enough data to fit a model? Honeycomb data is saved either way.
  is_modelable <- sp_code %in% species_sel$species_id

  # ---- 6.2 Load species data (safe-date filtering + pred_doy come from 06) ----
  sp <- load_sp_dat(sp_path, all_surveys,
                    survey_types = mode$survey_types,
                    survey_durations = mode$survey_durations)

  sp_dat        <- sp$sp_dat
  sp_safe_dates <- sp$sp_safe_dates
  pred_doy      <- sp$pred_doy

  sp_dat <- sp_dat %>% subset(!(Survey_Type == "ARU" & Survey_Duration_Minutes != 5))
  
  if (nrow(sp_dat) == 0) {
    message("Skipping ", sp_name, "; no surveys of the required types remain.")
    next
  }

  n_det <- sum(sp_dat$count > 0)
  if (n_det == 0) {
    message("Skipping ", sp_name, "; no positive counts remain after filtering.")
    next
  }

  # ---- 6.3 Likelihood family and honeycomb survey data ----
  # Only survey_id + count are stored; geometry/attributes live once in
  # all_surveys and are rejoined by survey_id downstream (08/09).
  error_family <- choose_error_family(
    count = subset(sp_dat, Survey_Type %in% c("ARU","Point_Count"))$count,
    nb_switch_count = nb_switch_count,
    nb_switch_min_n = nb_switch_min_n
  )

  save_atomic(
    list(
      sp_english    = sp_name,
      sp_code       = sp_code,
      is_modelable  = is_modelable,
      sp_safe_dates = sp_safe_dates,
      pred_doy      = pred_doy,
      survey_counts = sp_dat %>% st_drop_geometry() %>% select(survey_id, count),
      error_family  = error_family
    ),
    dat_path
  )

  # ---- 6.4 Decide whether to fit a model ----
  if (!is_modelable) {
    message("Honeycomb data saved for ", sp_name, "; not enough data to fit a model.")
    next
  }
  if (file.exists(pred_path) && !rerun_predictions) {
    message("Predictions already exist for ", sp_name, "; skipping model fit.")
    next
  }

  # ---- 6.5 Species-specific covariate table ----
  # Retain only candidate covariates that exist in this species dataset AND have
  # more than one unique finite value after filtering.
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

  # sd_linear kept tight (0.5) on purpose: loose habitat betas on sparse (0/1)
  # detections drive checkerboard artefacts in the abundance maps.
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
    filter(n > 500)

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

  # ---- Construct the 1-D DOY mesh (linear domain, per species) ----
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

  mod <- try(
    fit_model(
      sp_dat     = sp_dat,
      mesh_abund = mesh_abund,
      mesh_chg   = mesh_chg,
      mesh_TOD   = mesh_TOD,
      covariates = cov_df_sp,

      prior_range_abund = priors_list$prior_range_abund,
      prior_sigma_abund = priors_list$prior_sigma_abund,
      prior_range_change = priors_list$prior_range_change,
      prior_sigma_change = priors_list$prior_sigma_change,

      fixed_TOD_range = priors_list$fixed_TOD_range,
      prior_TOD_sigma = priors_list$prior_TOD_sigma,

      mesh_DOY        = mesh_DOY,
      fixed_DOY_range = priors_list$fixed_DOY_range,
      prior_DOY_sigma = priors_list$prior_DOY_sigma,

      kappa_pcprec_diff = priors_list$kappa_pcprec_diff,

      int_strategy = int_strategy,
      strategy     = strategy,
      family       = error_family
    ),
    silent = TRUE
  )

  if (inherits(mod, "try-error") || is.null(mod)) {
    message("Model failed for ", sp_name, "; skipping this species.")
    print(mod)
    next
  }

  end_model   <- Sys.time()
  fit_minutes <- round(as.numeric(end_model - start_model, units = "mins"), 1)

  # ---- 6.7b Reference time of day, from the fitted TOD curve ----
  # Read off TOD_global (estimated alongside the protocol intercepts, so not
  # confounded), restricted to hours_supported so the reference is never placed
  # where the curve extrapolates. Falls back to the empirical peak if unavailable.
  optimal_TOD_fitted <- optimal_TOD_from_fit(mod, hours_allowed = hours_supported)

  grid  <- seq(min(hours_supported), max(hours_supported), by = 0.05)
  curve <- TOD_curve_from_fit(mod, grid_hours = grid)
  plot(pred ~ Hours_After_Reference, data = curve, type = "l", lwd = 2,
       xlab = "Hours after reference", ylab = "TOD_global (log scale)")
  rect(hours_supported - 0.5, par("usr")[3], hours_supported + 0.5, par("usr")[4],
       col = adjustcolor("dodgerblue", 0.3), border = NA)
  abline(h = 0, col = "grey70")

  # ---- 6.7c Fitted day-of-year curve ----
  doy_range <- range(sp_dat$days_midpoint, na.rm = TRUE)
  doy_grid  <- seq(doy_range[1], doy_range[2], by = 0.5)
  doy_curve <- DOY_curve_from_fit(mod, grid_days = doy_grid)
  if (!is.null(doy_curve)) {
    plot(pred ~ days_midpoint, data = doy_curve, type = "l", lwd = 2,
         xlab = "Days from safe-date midpoint", ylab = "DOY_global (log scale)",
         main = paste0(sp_name, " - day-of-year effect"))
    abline(h = 0, col = "grey70")
    abline(v = 0, col = "dodgerblue", lwd = 2, lty = 2)
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
    error_family     = error_family,
    priors           = priors_list,
    int_strategy     = int_strategy,
    strategy         = strategy,
    n_surveys        = nrow(sp_dat),
    n_detections     = n_det,
    n_covariates     = nrow(cov_df_sp),
    covariates       = cov_df_sp,
    pred_doy         = pred_doy,
    fitted_TOD_curve = curve,      # does not show uncertainty
    fitted_DOY_curve = doy_curve,  # does not show uncertainty

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

  print(summary(mod))
  message("\n====================\n", i, "/", nrow(species_all), ": ", sp_name,
          " (species_id = ", sp_code, "); ", fit_minutes, " min to fit model\n====================")

  # ---- 6.8 Full-grid predictions ----
  # Standardized to optimal_TOD, days_midpoint = 0, kappa_diff = 0, site_iid = 0
  # (the last two omitted from pred_formula -> exp(eta) is a geometric mean;
  # 6.9 optionally rescales to the arithmetic mean).
  start_prediction <- Sys.time()
  message("Generating predictions for full 1-km grid")

  pred_formula <- make_pred_formula_multiatlas(cov_df_sp)

  pred_grid <- make_pred_grid(grid_OBBA2, grid_OBBA3) %>%
    mutate(Hours_After_Reference = optimal_TOD, days_midpoint = 0)

  preds <- predict_all_pixels(
    mod = mod, pred_grid = pred_grid, pred_formula = pred_formula,
    n.samples = n_prediction_draws, seed = prediction_seed
  )

  # ---- 6.9 Optional lognormal variance correction (draw matrices) ----
  # Rescales exp(eta) from the geometric to the arithmetic mean over the omitted
  # iid levels. Disabled by default: vc_terms = character(0) -> factor 1 (no-op).
  # Applied HERE so pixel summaries, abs_change and hex_draws all come from the
  # same (corrected or uncorrected) draws. NB overdispersion is intentionally not
  # corrected (E[y] = mu already).
  vc <- make_lognormal_correction(
    mod = mod, terms = vc_terms, propagate = vc_propagate,
    n_draws = n_prediction_draws, cap_quantile = vc_cap_quantile
  )
  message("  ", describe_lognormal_correction(vc))

  preds$mu2 <- apply_lognormal_correction(preds$mu2, vc$factor)
  preds$mu3 <- apply_lognormal_correction(preds$mu3, vc$factor)

  # Stored without the large per-draw factor vector; factor_summary keeps its
  # median and 5th/95th percentiles.
  vc_record <- vc[setdiff(names(vc), "factor")]
  model_summaries[[sp_name]]$variance_correction <- vc_record
  save_atomic(model_summaries, model_summaries_path)

  # ---- 6.10 Terrestrial open-water correction ----
  # Scales predicted abundance by the non-open-water fraction of each pixel.
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

      error_family = error_family,
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

      # Scale metadata (present when the variance correction was applied; its
      # factor is 1 when vc_terms is empty, i.e. surfaces are on the geometric scale).
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
