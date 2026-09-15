# ============================================================
# inla_model_utils.R
#
# Helper functions for the INLA/inlabru multi-atlas modelling workflow.
#
# Contents, in the order the pipeline uses them:
#
#   1. Small utilities        ensure_numeric, load_or_empty_list, save_atomic,
#                             sp_filename, make_cov_df, calibrate_pc_lambda
#   2. Per-species datasets   sp_data_path, make_sp_dat_record, load_sp_dat,
#                             select_modelable_species, choose_error_family
#   3. Model fitting          make_TOD_mesh, optimal_TOD_from_fit,
#                             fit_PC_ARU, fit_PC_ARU_CL
#   4. Prediction             make_pred_formula_multiatlas, make_pred_grid,
#                             predict_inla, predict_all_pixels
#   5. Variance correction    get_var, has_hyper, sample_hyper_vars,
#                             make_lognormal_correction,
#                             apply_lognormal_correction,
#                             describe_lognormal_correction
#   6. Posterior summaries    summarize_posterior, summarize_predictions,
#                             make_hex_draws
#   7. Shared-footprint       make_shared_footprint_dataset,
#      (paired) analysis      fit_inla_shared_footprint_change,
#                             summarize_shared_footprint_change
#
# Sourced by 06 (section 2 only), 07, 07b and 08.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
})

# ============================================================
# 1. Small utilities
# ============================================================

# Coerce to numeric (thin wrapper)
#
# Used to make intent explicit when preparing inputs for INLA/inlabru,
# where factors/characters can otherwise silently propagate and break
# model fitting.
#
# x Vector to coerce.
# Numeric vector.
ensure_numeric <- function(x) as.numeric(x)

# Load an RDS object or return an empty list if missing
#
# Convenience helper for incremental pipelines that cache intermediate
# results. If `path` doesn't exist, returns `list()` rather than erroring.
#
# path RDS file path.
# Loaded object, or empty list.
load_or_empty_list <- function(path) {
  if (file.exists(path)) readRDS(path) else list()
}

# Save an object atomically (local copy)
#
# Same idea as in `cv_utils.R`: write to a temp file then rename, preventing
# corrupted `.rds` files. This local definition keeps `inla_model_utils.R`
# self-contained.
#
# obj Object to save.
# path Output path.
# Invisibly, the path.
save_atomic <- function(obj, path) {
  tmp <- paste0(path, ".tmp")
  saveRDS(obj, tmp)
  file.rename(tmp, path)
}

# Convert a species common name into a safe filename slug
#
# Replaces spaces/punctuation with underscores and normalizes case so
# file outputs are consistent and portable across operating systems.
#
# sp_english Species common name.
# Filename-safe character string.
sp_filename <- function(sp_english) {
  sp_english %>%
    str_to_lower() %>%
    str_replace_all("[^a-z0-9]+", "_") %>%
    str_replace_all("^_|_$", "")
}

# Build a covariate prior specification dataframe (simple linear priors)
make_cov_df <- function(covars, mean = 0, sd_linear = 1) {
  tibble(
    covariate = covars,
    beta = 1,
    sd_linear = sd_linear,
    model = "linear",
    mean = 0,
    prec = 1 / (sd_linear^2)
  )
}

# Resolve the per-type special-survey square threshold from min_special_squares.
#
# min_special_squares is either a single number (applied to every type) or a
# named vector: named types use their value; unlisted types fall back to a
# "default" entry if present, else the first entry. Used by BOTH the main loop
# (to decide whether to fit a species at all in a special-survey mode) and
# fit_PC_ARU_CL() (to decide which intercepts to keep), so the two cannot drift.
resolve_special_threshold <- function(type, min_special_squares) {
  v <- min_special_squares
  if (is.null(names(v)) || all(names(v) == "")) return(as.numeric(v[[1]]))
  if (type %in% names(v))      return(as.numeric(v[[type]]))
  if ("default" %in% names(v)) return(as.numeric(v[["default"]]))
  as.numeric(v[[1]])
}

#' Calibrate lambda for a PC prior via Monte Carlo tail probability
#'
#' Some PC priors are parameterized via a rate/scale (lambda) chosen so that
#' `P(theta > threshold) = target_prob`. This helper finds lambda by root finding
#' using a Monte Carlo approximation to the tail probability.
#'
#' @param target_prob Desired tail probability.
#' @param threshold_theta Threshold defining the tail event.
#' @param n_mc Number of Monte Carlo draws used to approximate the probability.
#' @param lower,upper Search bounds for lambda.
#' @return Calibrated lambda value.
calibrate_pc_lambda <- function(target_prob,
                                threshold_theta,
                                n_mc = 20000,
                                lower = 1e-6,
                                upper = 20,
                                tol = 1e-3) {
  stopifnot(is.numeric(target_prob), length(target_prob) == 1, target_prob > 0, target_prob < 1)
  stopifnot(is.numeric(threshold_theta), length(threshold_theta) == 1, threshold_theta > 0)
  
  #' Monte Carlo estimate of PC prior tail probability
  #'
  #' Internal helper used by `calibrate_pc_lambda()`; estimates
  #' `P(theta > threshold)` for a given lambda by simulation.
  #'
  #' @param lambda Candidate lambda.
  #' @param n Number of draws.
  #' @param threshold Threshold for tail event.
  #' @return Estimated probability.
  pc_tail_prob <- function(lambda, n = n_mc, threshold = threshold_theta) {
    samp <- INLA::inla.pc.rgamma(n = n, lambda = lambda)
    mean(samp > threshold)
  }
  
  #' Root function for PC lambda calibration
  #'
  #' Returns `pc_tail_prob(lambda) - target_prob` so `uniroot()` can solve
  #' for the lambda where the tail probability matches the target.
  #'
  #' @param lambda Candidate lambda.
  #' @return Signed difference from target probability.
  f_for_root <- function(lambda) pc_tail_prob(lambda) - target_prob
  
  sol <- try(uniroot(f_for_root, lower = lower, upper = upper, tol = tol), silent = TRUE)
  if (inherits(sol, "try-error")) {
    stop("Could not find lambda in the specified interval. Try increasing `upper` or `n_mc`.")
  }
  sol$root
}


# ============================================================
# 2. Per-species analysis datasets
# ============================================================

# ============================================================
# Per-species analysis datasets
#
# Script 06 writes one lightweight record per species; scripts 07 and 07b read
# them back with load_sp_dat() and rehydrate a full sf `sp_dat`.
#
# WHAT IS STORED
#   Only `survey_id` + `count`, plus species metadata. Geometry, covariates and
#   every other survey attribute live once in `all_surveys` and are rejoined by
#   `survey_id` on load. Storing full sf objects per species would duplicate the
#   geometry and covariate columns 400+ times (tens of GB); the natural key
#   costs a few MB per species and cannot drift out of sync with all_surveys.
#
# WHAT IS FILTERED (in 06, once)
#   - Biological regions with no safe-date window for the species are dropped.
#   - Surveys outside the species x region safe-date window are dropped.
#   Survey types are NOT filtered: every record holds point counts, ARUs and
#   checklists, i.e. the PC_ARU_CL superset. Downstream scripts narrow this:
#     07  PC_ARU     -> survey_types = c("Point_Count", "ARU")
#     07  PC_ARU_CL  -> survey_types = NULL   (keep everything)
#     07b paired     -> survey_types = "Point_Count"
# ============================================================

# Path to one species' record. Kept in one place so 06/07/07b cannot disagree.
sp_data_path <- function(dir, sp_english) {
  file.path(dir, paste0(sp_filename(sp_english), "_sp_dat.rds"))
}


# Build the per-species record written by script 06.
#
# surveys_f      Filtered survey sf; must carry survey_id, Biol_Region, DayOfYear.
# count_vec      Species count vector, row-aligned with surveys_f.
# sp_safe_dates  MODELLING windows: one row per biological region, every window
#                filled. Regions the species breeds in keep their own safe-date
#                window; regions with no safe dates borrow the narrowest
#                cross-regional window (assigned in 06). Used HERE to decide which
#                surveys enter the model -- including the (mostly zero) surveys in
#                no-safe-date regions that bend the SPDE toward zero.
#                Must have Biol_Region, start_doy, end_doy (all non-NA).
# sp_safe_dates_unmodified
#                ORIGINAL safe-date definitions: same rows, but start_doy/end_doy
#                are NA in regions where the species has no safe dates. Stored
#                verbatim and NOT used for filtering here; script 07 ships it into
#                each prediction file and script 08 uses it to zero predictions in
#                the no-safe-date regions. Keeping both objects on the record is
#                what lets "fit everywhere, predict only where safe" stay coherent.
#
# Returns a list, or NULL if no surveys survive the safe-date filter.
make_sp_dat_record <- function(surveys_f,
                               count_vec,
                               sp_english,
                               species_id,
                               sp_safe_dates,
                               sp_safe_dates_unmodified) {
  
  stopifnot("survey_id" %in% names(surveys_f))
  stopifnot(nrow(surveys_f) == length(count_vec))
  
  if (nrow(sp_safe_dates) == 0) return(NULL)
  
  # Prediction reference day: midpoint of the window shared by all regions the
  # species occupies. Computed once here so 07 and 07b cannot diverge. na.rm is
  # defensive; sp_safe_dates windows are filled (non-NA) by the time we get here.
  pred_doy <- (min(sp_safe_dates$end_doy,   na.rm = TRUE) +
               max(sp_safe_dates$start_doy, na.rm = TRUE)) / 2
  
  windows <- sp_safe_dates %>%
    dplyr::select(Biol_Region, start_doy, end_doy)
  
  sp_dat <- surveys_f %>%
    sf::st_drop_geometry() %>%
    dplyr::mutate(count = count_vec) %>%
    dplyr::filter(Biol_Region %in% windows$Biol_Region) %>%
    dplyr::left_join(windows, by = "Biol_Region") %>%
    dplyr::filter(
      !is.na(start_doy),
      !is.na(end_doy),
      DayOfYear >= start_doy,
      DayOfYear <= end_doy
    )
  
  if (nrow(sp_dat) == 0) return(NULL)
  
  list(
    sp_english               = sp_english,
    species_id               = as.character(species_id),
    sp_safe_dates            = sp_safe_dates,             # filled: what 06 filtered surveys on
    sp_safe_dates_unmodified = sp_safe_dates_unmodified,  # original: what 08 zeros predictions with
    pred_doy                 = pred_doy,
    survey_counts            = dplyr::select(sp_dat, survey_id, count),
    n_surveys                = nrow(sp_dat),
    n_detections             = sum(sp_dat$count > 0),
    date_created             = Sys.time()
  )
}


# Read a species record and rehydrate the full modelling dataset.
#
# path          Record written by script 06 (see sp_data_path()).
# all_surveys   The `all_surveys` sf from data_ready_for_analysis.rds.
# survey_types  Character vector of Survey_Type values to keep, or NULL for all.
#
# Returns the record with an `sp_dat` sf attached, ready for the fit functions.
load_sp_dat <- function(path, all_surveys, survey_types = NULL,survey_durations = NULL) {
  
  rec <- readRDS(path)
  
  sp_dat <- all_surveys %>%
    dplyr::inner_join(rec$survey_counts, by = "survey_id")
  
  if (nrow(sp_dat) != nrow(rec$survey_counts)) {
    stop(
      "survey_id join is not 1:1 for ", rec$sp_english, " (",
      nrow(rec$survey_counts), " stored rows -> ", nrow(sp_dat),
      " joined rows). Was all_surveys rebuilt without rerunning 06?"
    )
  }
  
  if (!is.null(survey_types)) {
    sp_dat <- sp_dat %>% dplyr::filter(Survey_Type %in% survey_types)
  }
  
  # Derived covariates depend on pred_doy, so they are rebuilt on load rather
  # than stored. days_midpoint and the ARU duration indicators are cheap and
  # keep the record free of anything the fit functions might later redefine.
  #
  # ARU_1min / ARU_3min / ARU_5min are MUTUALLY EXCLUSIVE: durations are
  # restricted to {1, 3, 5} upstream (06 + the survey_durations filter below),
  # so every ARU row loads exactly one of them and point counts load none.
  # Paired with the matching effect_ARU_* components, each effect is a DIRECT
  # contrast of that ARU protocol against a 5-minute point count (rather than a
  # deviation from the 5-minute ARU), so the three ARU effects are directly
  # comparable in the model output.
  sp_dat <- sp_dat %>%
    dplyr::mutate(
      days_midpoint = DayOfYear - rec$pred_doy,
      ARU_1min      = as.numeric(ARU == 1 & Survey_Duration_Minutes == 1),
      ARU_3min      = as.numeric(ARU == 1 & Survey_Duration_Minutes == 3),
      ARU_5min      = as.numeric(ARU == 1 & Survey_Duration_Minutes == 5)
    )
  
  if (!is.null(survey_durations)){
    # Duration filter applies ONLY to structured surveys (point counts + ARUs).
    # Checklists ("Breeding Bird Atlas", "Linear transect") have variable
    # duration by design and must never be dropped by the ARU/PC duration filter.
    sp_dat <- sp_dat %>%
      filter(!(Survey_Type %in% c("Point_Count", "ARU")) |
               Survey_Duration_Minutes %in% survey_durations)
  }
  
  rec$sp_dat <- sp_dat
  rec
}


# Species that have enough data to support a model fit.
#
# Applied identically by 07 and 07b so the two analyses run on the same species
# list. The detection thresholds are always evaluated on point counts and ARUs,
# regardless of which survey types a given mode ends up fitting.
select_modelable_species <- function(species_detection_summaries,
                                     min_detections = 50,
                                     min_squares    = 20,
                                     survey_types   = c("Point_Count", "ARU")) {
  
  species_all <- species_detection_summaries %>%
    dplyr::group_by(sp_english, species_id) %>%
    dplyr::summarise(n_det = sum(n_det), .groups = "drop")
  
  model_ids <- species_detection_summaries %>%
    dplyr::filter(
      Survey_Type %in% survey_types,
      n_det >= min_detections,
      n_sq  >  min_squares
    ) %>%
    dplyr::pull(species_id) %>%
    as.character() %>%
    unique()
  
  list(species_all = species_all, model_ids = model_ids)
}


# Choose poisson vs nbinomial from the observed counts.
#
# No counts are trimmed. Negative binomial is used when the top 1% of nonzero
# counts hold >= `top1pct_share` of all individuals, or when more than
# `nb_switch_min_n` surveys exceed `nb_switch_count`.
choose_error_family <- function(count,
                                nb_switch_count = 20,
                                nb_switch_min_n = 20,
                                top1pct_share   = 0.20) {
  
  y_pos <- count[count > 0]
  if (length(y_pos) == 0) return(NA_character_)
  
  n_top <- max(1, ceiling(0.01 * length(y_pos)))
  prop_top1pct <- sum(sort(y_pos, decreasing = TRUE)[seq_len(n_top)]) / sum(y_pos)
  
  if (prop_top1pct >= top1pct_share ||
      sum(y_pos > nb_switch_count) > nb_switch_min_n) {
    "nbinomial"
  } else {
    "poisson"
  }
}

# ============================================================
# 3. Model fitting
# ============================================================

# ==============================================================================
# make_TOD_mesh()  --  ONE global 1-D mesh for the time-of-day smooth
#
# Built once and reused for every species. Previously each fit derived its own
# mesh from range(sp_dat$Hours_After_Reference), so a species with only
# point counts got a mesh spanning ~1-10 h while another spanned 0-24 h. Because
# the SPDE range is an ABSOLUTE length, its meaning relative to mesh resolution
# then differed between species and the smooths were not comparable.
#
# Two properties matter:
#
#   1. KNOT SPACING <= range/5. The finite-element basis cannot represent
#      variation finer than its knot spacing, so a coarse mesh silently
#      overrides the range: the old `length.out = 11` over 0-24 h gave 2.6 h
#      spacing, which made any range prior below ~3 h inoperative.
#
#   2. CYCLIC DOMAIN. Hours_After_Reference is measured from 3 h before
#      sunrise, so hour 24 is ~1 h before the next day's hour 0. On a free
#      mesh those are maximally distant and both ends carry inflated boundary
#      variance -- exactly where nocturnal detections and the ARU evening peak
#      (16-23 h) live. A cyclic mesh joins the seam.
#
# If the installed fmesher/INLA does not support boundary = "cyclic", falls
# back to a free mesh padded by one full range at each end, which at least
# keeps the boundary effect away from the data.
#
# span         Domain of Hours_After_Reference. Default c(0, 24).
# knot_spacing Mesh resolution in hours. Keep <= fixed_TOD_range / 5.
# pad_range    Padding used only by the non-cyclic fallback.
# ==============================================================================

make_TOD_mesh <- function(span         = c(0, 24),
                          knot_spacing = 0.25,
                          pad_range    = 1.5,
                          cyclic       = TRUE) {
  
  stopifnot(length(span) == 2, diff(span) > 0, knot_spacing > 0)
  
  if (cyclic) {
    # Domain is [span[1], span[2]); the final knot is one step short of the
    # upper bound, otherwise the wrapped endpoint is duplicated.
    loc <- seq(span[1], span[2] - knot_spacing, by = knot_spacing)
    
    mesh <- try(
      INLA::inla.mesh.1d(
        loc      = loc,
        interval = span,
        degree   = 2,
        boundary = "cyclic"
      ),
      silent = TRUE
    )
    
    if (!inherits(mesh, "try-error")) {
      return(mesh)
    }
    
    warning(
      "Cyclic 1-D mesh unsupported by this INLA version; ",
      "falling back to a padded free mesh. Dusk and pre-dawn will not ",
      "share information across the 0/24 h seam.",
      call. = FALSE
    )
  }
  
  INLA::inla.mesh.1d(
    loc = seq(span[1] - pad_range,
              span[2] + pad_range,
              by = knot_spacing),
    degree   = 2,
    boundary = "free"
  )
}


# ==============================================================================
# optimal_TOD_from_fit()  --  reference hour from the FITTED TOD curve
#
# The empirical alternative -- the hour with the highest raw detection rate --
# is confounded with protocol, because survey types are not evenly distributed
# across the day (hours 16-23 are effectively ARU-only, 10-16 are checklist
# -dominated). A raw peak in a single-protocol stretch of the day may be a
# protocol effect wearing a time-of-day label.
#
# TOD_global is estimated with the protocol intercepts already in the model, so
# reading the reference hour off the fitted curve removes that confounding.
#
# `hours_allowed` restricts the search to hours with adequate survey support, so
# the reference is never placed where the curve is extrapolating. Returns NA on
# failure; the caller should fall back to the empirical estimate.
#
# mod           Fitted inlabru model.
# hours_allowed Numeric vector of admissible hours (see coverage table in 07).
# by            Grid resolution for the search, in hours.
# ==============================================================================

TOD_curve_from_fit <- function(mod, grid_hours, comp_name = "TOD_global") {
  w <- mod$summary.random[[comp_name]]$mean
  if (is.null(w)) return(NULL)
  mapper <- tryCatch(mod$bru_info$model$effects[[comp_name]]$main$mapper,
                     error = function(e) NULL)
  if (is.null(mapper)) return(NULL)
  vals <- try(as.vector(inlabru::ibm_eval(mapper, input = grid_hours, state = w)),
              silent = TRUE)
  if (inherits(vals, "try-error") || length(vals) != length(grid_hours)) return(NULL)
  data.frame(Hours_After_Reference = grid_hours, pred = vals)
}

optimal_TOD_from_fit <- function(mod, hours_allowed, by = 0.25,
                                 tol = 0.5, comp_name = "TOD_global") {
  
  hours_allowed <- sort(unique(hours_allowed[is.finite(hours_allowed)]))
  if (!length(hours_allowed)) return(NA_real_)
  
  grid_hours <- seq(min(hours_allowed), max(hours_allowed), by = by)
  
  # Restore the gap filter: never search hours with no survey support.
  keep <- vapply(grid_hours,
                 function(h) any(abs(h - hours_allowed) <= tol),
                 logical(1))
  grid_hours <- grid_hours[keep]
  if (!length(grid_hours)) return(NA_real_)
  
  curve <- TOD_curve_from_fit(mod, grid_hours, comp_name)
  if (is.null(curve)) return(NA_real_)
  
  curve$Hours_After_Reference[which.max(curve$pred)]
}



# ==============================================================================
# make_DOY_mesh()  --  ONE global 1-D mesh for the day-of-year smooth
#
# The day-of-year analogue of make_TOD_mesh(). Built once and reused for every
# species so the fixed range means the same thing everywhere, exactly as for
# TOD. The one structural difference is that the DOY domain is NOT cyclic:
# days_midpoint is a linear axis centred on the shared safe-date midpoint (0),
# so a free 1-D mesh is used (no 0/24-style seam to join). The domain is padded
# by one full range at each end to keep boundary variance away from the data.
#
# span         Domain of days_midpoint (days from the safe-date midpoint).
#              Default c(-60, 60) comfortably covers any breeding safe-date
#              window; widen it if a species' window is unusually long.
# knot_spacing Mesh resolution in days. Keep <= fixed_DOY_range / 5.
# pad_range    Padding added at each end (defaults to one full range).
# ==============================================================================

make_DOY_mesh <- function(span         = c(-60, 60),
                          knot_spacing = 5,
                          pad_range    = 30) {
  
  stopifnot(length(span) == 2, diff(span) > 0, knot_spacing > 0)
  
  INLA::inla.mesh.1d(
    loc = seq(span[1] - pad_range,
              span[2] + pad_range,
              by = knot_spacing),
    degree   = 2,
    boundary = "free"
  )
}


# ==============================================================================
# DOY_curve_from_fit()  --  evaluate the fitted DOY smooth on a grid of days
#
# The day-of-year analogue of TOD_curve_from_fit(). Reads the posterior-mean
# weights of the DOY_global component and maps them to a fine grid of
# days_midpoint values via the same inlabru mapper, returning a data.frame with
# columns days_midpoint / pred (the log-scale seasonal effect). Returns NULL if
# the component is absent (e.g. a model fitted before DOY was added).
#
# mod        Fitted inlabru model.
# grid_days  Numeric vector of days_midpoint values to evaluate.
# ==============================================================================

DOY_curve_from_fit <- function(mod, grid_days, comp_name = "DOY_global") {
  w <- mod$summary.random[[comp_name]]$mean
  if (is.null(w)) return(NULL)
  mapper <- tryCatch(mod$bru_info$model$effects[[comp_name]]$main$mapper,
                     error = function(e) NULL)
  if (is.null(mapper)) return(NULL)
  vals <- try(as.vector(inlabru::ibm_eval(mapper, input = grid_days, state = w)),
              silent = TRUE)
  if (inherits(vals, "try-error") || length(vals) != length(grid_days)) return(NULL)
  data.frame(days_midpoint = grid_days, pred = vals)
}



# ==============================================================================
# fit_PC_ARU()  --  point counts + ARUs, with a site-level iid random effect
#
# Repeated visits to one station (especially short ARU recordings) are
# pseudoreplicates. The site term quarantines that within-site correlation. It
# is built AFTER species-specific filtering, and is dropped entirely when too
# few multi-visit sites remain to identify sigma_site.
#
# Checklists are not modelled here -- use fit_PC_ARU_CL() for those.
#
# Requires sp_dat columns: count, Atlas3, Atlas3_c, Hours_After_Reference,
# days_midpoint, ARU, Survey_Type, square_atlas, site_id, geometry.
# ==============================================================================

fit_PC_ARU <- function(
    sp_dat,
    mesh_abund,
    mesh_chg,
    covariates,
    mesh_TOD = NULL,
    timeout_min = 40,
    
    # Spatial SPDE priors
    prior_range_abund  = c(300, 0.1),   # P(range < 300 km) = 0.1  -> favours ranges > 300 km
    prior_sigma_abund  = c(0.5, 0.1),   # P(sigma > 0.5) = 0.1     -> prior median ~0.15 (log); caps field amplitude
    prior_range_change = c(300, 0.1),   # P(range < 300 km) = 0.1  -> favours ranges > 300 km
    prior_sigma_change = c(0.3, 0.1),   # P(sigma > 0.3) = 0.1     -> prior median ~0.09; tighter than abundance (change is a noisy difference)
    
    # 1-D SPDE: survey timing. FIXED range + estimated sigma, as for the
    # effort smooths. Range and marginal variance are only jointly identified
    # through sigma^2 * kappa^(2*nu), so on a bounded 1-D domain the likelihood
    # is a ridge and the range collapses for sparse species. Fixing it also
    # keeps the smoothness of the detectability correction COMPARABLE ACROSS
    # SPECIES, which an estimated range would not.
    
    fixed_TOD_range = 1.5,                # hours; NOT estimated
    prior_TOD_sigma = c(2, 0.1),          # P(sigma > 2) = 0.1 -> prior median ~0.6
    
    # 1-D SPDE: day of year (seasonal phenology / detectability). Modelled
    # exactly like TOD -- shared mesh (mesh_DOY), FIXED range (default 30 days),
    # estimated sigma.
    mesh_DOY = NULL,
    fixed_DOY_range = 30,                 # days; NOT estimated
    prior_DOY_sigma = c(2, 0.1),          # P(sigma > 2) = 0.1 -> prior median ~0.6

    # ARU protocol effects (ARU_1min / ARU_3min / ARU_5min). Each is a contrast
    # of that ARU protocol against a 5-minute point count on the log link; all
    # three share ONE prior: mean 0, SD = aru_prior_sd. The default log(5)/2
    # (SD ~= 0.8) puts ~95% of prior mass within a factor of ~5 of the PC
    # reference -- wide enough to admit the ARU-vs-observer method offset plus
    # the shortening from 5 to 1 minute. That shortening is modest, not the naive
    # log(1/5): most birds at a location are detected in the first minute, so a
    # 1-minute window recovers the bulk of a 5-minute count. mean 0 stays
    # agnostic about sign.
    aru_prior_sd = log(5) / 2,
    #
    # Atlas-square iid prior
    kappa_pcprec_diff = c(0.5, 0.1),    # P(sigma_square > 0.5) = 0.1 -> fine-scale (10-km square) sink

    # include_square: master switch for the atlas-square iid random effect
    #   (kappa_diff). TRUE keeps it in the model; FALSE forces it OFF for
    #   every species. Symmetric to include_site below, so the site- and
    #   square-level effects can be toggled independently.
    include_square    = TRUE,

    # Site-level iid prior and inclusion rules.
    # min_site_surveys: minimum RETAINED surveys for a site to get a level.
    # min_site_levels:  minimum number of such sites, else the term is dropped
    #   (sigma_site cannot be identified from a handful of sites).
    prior_site_pcprec = c(1, 0.5),    # P(sigma_site > 1) = 0.5  (default; script 07 overrides with c(1, 0.1))
    
    # include_site: master switch for the site-level iid random effect.
    #   TRUE  -> keep the data-driven gate below (the term is used only
    #            where it is identifiable: n_site >= min_site_levels).
    #   FALSE -> force the term OFF for every species, regardless of how
    #            many multi-visit sites exist. This is the switch that the
    #            fit_*_nosite() wrappers set; everything else is unchanged.
    include_site      = TRUE,

    min_site_surveys  = 2,
    min_site_levels   = 10,
    
    # Likelihood
    family = c("poisson", "nbinomial"),
    
    # Negative-binomial PC prior settings
    nb_pc_target_prob     = 0.5,
    nb_pc_threshold_theta = 5,
    
    # inlabru / INLA options
    inla_mode    = "experimental",
    int_strategy = "ccd",
    strategy     = "simplified.laplace",
    bru_verbose  = 4,
    waic         = FALSE,
    cpo          = FALSE,
    retry        = 0
) {
  
  # ------------------------------------------------------------
  # 1. Validate inputs
  # ------------------------------------------------------------
  
  family <- match.arg(family)
  
  if (!inherits(sp_dat, "sf")) {
    stop("sp_dat must be an sf object with a geometry column.")
  }
  
  if (!inherits(study_boundary, "sf")) {
    stop("study_boundary must be an sf object.")
  }
  
  required_cols <- c(
    "count",
    "Atlas3",
    "Atlas3_c",
    "Hours_After_Reference",
    "days_midpoint",
    "ARU",
    "ARU_1min",
    "ARU_3min",
    "ARU_5min",
    "Survey_Type",
    "square_atlas",
    "site_id",
    "geometry"
  )
  
  missing_cols <- setdiff(required_cols, names(sp_dat))
  
  if (length(missing_cols) > 0) {
    stop(
      "sp_dat is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  if (!is.data.frame(covariates)) {
    stop("covariates must be a data.frame.")
  }
  
  needed_cov_cols <- c("covariate", "beta", "model", "mean", "prec")
  cov_missing <- setdiff(needed_cov_cols, names(covariates))
  
  if (length(cov_missing) > 0) {
    stop(
      "covariates is missing required columns: ",
      paste(cov_missing, collapse = ", ")
    )
  }
  
  # ------------------------------------------------------------
  # 2. Set INLA options and align CRS
  # ------------------------------------------------------------
  
  INLA::inla.setOption(inla.timeout = 60 * timeout_min)
  
  if (sf::st_crs(sp_dat) != sf::st_crs(study_boundary)) {
    study_boundary <- sf::st_transform(study_boundary, sf::st_crs(sp_dat))
  }
  
  # ------------------------------------------------------------
  # 3. Ensure model variables have expected types
  # ------------------------------------------------------------
  
  sp_dat$Hours_After_Reference <- ensure_numeric(sp_dat$Hours_After_Reference)
  sp_dat$days_midpoint       <- ensure_numeric(sp_dat$days_midpoint)
  sp_dat$Atlas3_c            <- ensure_numeric(sp_dat$Atlas3_c)
  sp_dat$ARU                 <- ensure_numeric(sp_dat$ARU)
  sp_dat$ARU_1min            <- ensure_numeric(sp_dat$ARU_1min)
  sp_dat$ARU_3min            <- ensure_numeric(sp_dat$ARU_3min)
  sp_dat$ARU_5min            <- ensure_numeric(sp_dat$ARU_5min)

  # The three duration indicators must PARTITION the ARU rows: every ARU row
  # carries exactly one (duration in {1, 3, 5}) and every non-ARU row carries
  # none. This is what makes each effect_ARU_* a clean contrast against the
  # 5-minute point-count reference; a row with all-zero indicators would fold
  # silently back into that reference.
  stopifnot(all(sp_dat$ARU_1min + sp_dat$ARU_3min + sp_dat$ARU_5min == sp_dat$ARU))
  
  # ----------------------------------------------------------------------------
  # 3b. Site index for the repeated-measures random effect
  # ----------------------------------------------------------------------------
  # Computed HERE, not upstream: sp_dat has already been through species-specific
  # safe-date and survey-type filtering, so a site that had repeats in the raw
  # data may be a singleton for THIS species. Only sites with >= min_site_surveys
  # RETAINED surveys get a level.
  #
  # Singletons are excluded because they carry no within-site correlation to
  # absorb; an iid level informed by a single observation merely duplicates the
  # nbinomial overdispersion parameter. Division of labour:
  #     nbinomial  -> observation-level overdispersion, everywhere
  #     site_iid   -> shared correlation, only where repeats exist
  #
  # Singleton rows are pointed at level 1 and carry site_w = 0, so
  # `site_w * site_iid` gives them a zero Jacobian row: they contribute nothing
  # to any site node. This avoids depending on inlabru's NA-index handling,
  # which is mapper-dependent and has varied across versions.
  # ----------------------------------------------------------------------------
  
  site_tab    <- table(sp_dat$site_id)
  multi_sites <- names(site_tab)[site_tab >= min_site_surveys]
  is_multi    <- sp_dat$site_id %in% multi_sites
  
  sp_dat$site_idx_safe <- 1L
  sp_dat$site_idx_safe[is_multi] <-
    as.integer(factor(sp_dat$site_id[is_multi]))
  sp_dat$site_w <- as.numeric(is_multi)
  
  n_site       <- length(multi_sites)
  enough_sites <- n_site >= min_site_levels
  
  # The term is used only when the caller ALLOWED it (include_site) AND it is
  # identifiable (enough multi-visit sites). fit_*_nosite() sets the former
  # FALSE, which forces use_site_re FALSE no matter how many sites there are.
  use_site_re  <- isTRUE(include_site) && enough_sites
  
  site_re_note <- if (!isTRUE(include_site)) {
    "  -- TERM DISABLED (include_site = FALSE)"
  } else if (!enough_sites) {
    "  -- TERM DROPPED (too few sites)"
  } else {
    ""
  }
  
  message(
    "  site RE: ", n_site, " sites with >= ", min_site_surveys,
    " surveys, covering ", sum(is_multi), "/", nrow(sp_dat), " surveys",
    site_re_note
  )
  # ----------------------------------------------------------------------------
  # 3c. Square-level iid switch (kappa_diff)
  # ----------------------------------------------------------------------------
  # No data-driven gate: atlas squares are always plentiful, so the master
  # switch is the only control. n_square_levels is recorded on the fit object.
  use_square_re   <- isTRUE(include_square)
  n_square_levels <- length(unique(sp_dat$square_atlas))
  message(
    "  square RE: ", n_square_levels, " atlas-square levels",
    if (!use_square_re) "  -- TERM DISABLED (include_square = FALSE)" else ""
  )
  # ----------------------------------------------------------------------------
  
  # ------------------------------------------------------------
  # 6. Build spatial SPDE models
  # ------------------------------------------------------------
  
  matern_mean <- INLA::inla.spde2.pcmatern(
    mesh        = mesh_abund,
    prior.range = prior_range_abund,
    prior.sigma = prior_sigma_abund,
    constr      = TRUE
  )
  
  matern_diff <- INLA::inla.spde2.pcmatern(
    mesh        = mesh_chg,
    prior.range = prior_range_change,
    prior.sigma = prior_sigma_change,
    constr      = TRUE
  )
  
  # ------------------------------------------------------------
  # 7. Build 1-D SPDE smoothers for timing and effort corrections
  # ------------------------------------------------------------
  
  # --- Time-of-day smooth: SHARED mesh, FIXED range ---------------------------
  # The mesh is built once (make_TOD_mesh) and passed in, so every species uses
  # the same basis and the same absolute range means the same thing everywhere.
  # Deriving it from range(sp_dat$Hours_After_Reference) made the resolution
  # species-dependent.
  #
  # prior.range = c(range0, NA) -> range HELD FIXED at range0 (no prior).
  # Only sigma is estimated. Same idiom as the effort SPDEs.
  if (is.null(mesh_TOD)) {
    mesh_TOD <- make_TOD_mesh(knot_spacing = fixed_TOD_range / 6,
                              pad_range    = fixed_TOD_range)
  }
  
  # Guard the requirement that makes the fixed range meaningful: a basis coarser
  # than range/5 cannot represent the correlation length being asked for, and
  # the mesh, not the range, then sets the smoothness.
  TOD_knot_gap <- max(diff(sort(mesh_TOD$loc)))
  if (TOD_knot_gap > fixed_TOD_range / 5) {
    warning(
      "TOD mesh knot spacing (", signif(TOD_knot_gap, 3),
      " h) is coarser than fixed_TOD_range/5 (",
      signif(fixed_TOD_range / 5, 3),
      " h). The mesh, not the range, is setting the smoothness.",
      call. = FALSE
    )
  }
  
  TOD_spde <- INLA::inla.spde2.pcmatern(
    mesh        = mesh_TOD,
    alpha       = 2,
    prior.range = c(fixed_TOD_range, NA),
    prior.sigma = prior_TOD_sigma,
    constr      = TRUE
  )
  
  # --- Day-of-year smooth: SHARED mesh, FIXED range --------------------------
  # Modelled exactly like the time-of-day smooth above: a 1-D SPDE on a shared
  # mesh (make_DOY_mesh, passed in via mesh_DOY) with the RANGE HELD FIXED at
  # fixed_DOY_range and only sigma estimated. Unlike TOD the domain is NOT
  # cyclic -- days_midpoint is a linear axis centred on the safe-date midpoint.
  #
  # prior.range = c(range0, NA) -> range HELD FIXED at range0 (no prior).
  if (is.null(mesh_DOY)) {
    mesh_DOY <- make_DOY_mesh(knot_spacing = fixed_DOY_range / 6,
                              pad_range    = fixed_DOY_range)
  }
  
  # Same knot-spacing guard as for TOD: a basis coarser than range/5 cannot
  # represent the correlation length being asked for.
  DOY_knot_gap <- max(diff(sort(mesh_DOY$loc)))
  if (DOY_knot_gap > fixed_DOY_range / 5) {
    warning(
      "DOY mesh knot spacing (", signif(DOY_knot_gap, 3),
      " d) is coarser than fixed_DOY_range/5 (",
      signif(fixed_DOY_range / 5, 3),
      " d). The mesh, not the range, is setting the smoothness.",
      call. = FALSE
    )
  }
  
  DOY_spde <- INLA::inla.spde2.pcmatern(
    mesh        = mesh_DOY,
    alpha       = 2,
    prior.range = c(fixed_DOY_range, NA),
    prior.sigma = prior_DOY_sigma,
    constr      = TRUE
  )
  
  # ------------------------------------------------------------
  # 8. Define iid priors
  # ------------------------------------------------------------
  
  pc_prec_diff <- list(
    prior = "pcprec",
    param = kappa_pcprec_diff
  )
  
  pc_prec_site <- list(
    prior = "pcprec",
    param = prior_site_pcprec
  )
  
  # ------------------------------------------------------------
  # 9. Build covariate components and formula terms
  # ------------------------------------------------------------
  
  covariates <- covariates %>%
    dplyr::mutate(
      component = paste0(
        "Beta", beta, "_", covariate,
        '(1, model="', model,
        '", mean.linear=', mean,
        ", prec.linear=", prec,
        ")"
      ),
      term = paste0(
        "Beta", beta, "_", covariate,
        " * I(", covariate, "^", beta, ")"
      )
    )
  
  covar_components_str <- if (nrow(covariates) > 0) {
    paste(covariates$component, collapse = " + ")
  } else {
    ""
  }
  
  covar_terms_str <- if (nrow(covariates) > 0) {
    paste(covariates$term, collapse = " + ")
  } else {
    ""
  }
  
  # ------------------------------------------------------------
  # 10. Define inlabru components
  # ------------------------------------------------------------
  
  # Shared ARU protocol-effect prior (ARU_1min / ARU_3min / ARU_5min): mean 0,
  # SD = aru_prior_sd on the log link (see argument). The precision is injected
  # as a literal below so the fitted formula carries no environment dependency.
  aru_prec <- 1 / (aru_prior_sd^2)
  
  components_str <- paste0(
    "~ Intercept(1) + ",
    
    "effect_Atlas3(
       1,
       model = 'linear',
       mean.linear = 0,
       prec.linear = 1 / ((log(1.5) / 2)^2)
     ) + ",
    
    "effect_ARU_1min(
       1,
       model = 'linear',
       mean.linear = 0,
       prec.linear = ", aru_prec, "
     ) + ",
    
    "effect_ARU_3min(
       1,
       model = 'linear',
       mean.linear = 0,
       prec.linear = ", aru_prec, "
     ) + ",
    
    "effect_ARU_5min(
       1,
       model = 'linear',
       mean.linear = 0,
       prec.linear = ", aru_prec, "
     ) + ",
    
    "spde_mean(main = geometry, model = matern_mean) + ",
    "spde_diff(main = geometry, model = matern_diff) + ",
    
    "TOD_global(main = Hours_After_Reference, model = TOD_spde) + ",
    "DOY_global(main = days_midpoint, model = DOY_spde)",

    # Atlas-square iid, only when include_square = TRUE.
    if (use_square_re) {
      "+ kappa_diff(
         square_atlas,
         model = 'iid',
         constr = TRUE,
         hyper = list(prec = pc_prec_diff)
       )"
    } else {
      ""
    },
    
    # Site-level iid, only when there are enough multi-visit sites.
    # constr = TRUE (as for kappa_diff) makes the site effects sum to zero, so
    # dropping the term at prediction predicts the median site
    if (use_site_re) {
      "+ site_iid(
         site_idx_safe,
         model = 'iid',
         n = n_site,
         constr = TRUE,
         hyper = list(prec = pc_prec_site)
       )"
    } else {
      ""
    },
    
    if (nchar(covar_components_str) > 0) {
      paste0(" + ", covar_components_str)
    } else {
      ""
    }
  )
  
  model_components <- stats::as.formula(components_str)
  
  # ------------------------------------------------------------
  # 11. Define linear predictor
  # ------------------------------------------------------------
  
  # DOY_global +
  
  formula_str <- paste0(
    "count ~
      Intercept +
      TOD_global +
      DOY_global +
      

      ARU_1min * effect_ARU_1min +
      ARU_3min * effect_ARU_3min +
      ARU_5min * effect_ARU_5min +

      spde_mean +

      Atlas3_c * spde_diff +
      Atlas3_c * effect_Atlas3",

    # kappa_diff only when include_square = TRUE.
    if (use_square_re) " + kappa_diff" else "",

    # site_w is 0 at singleton sites, so those rows contribute nothing.
    # Same idiom as `ARU_5min * effect_ARU_5min` above.
    if (use_site_re) " + site_w * site_iid" else "",
    
    if (nchar(covar_terms_str) > 0) {
      paste0(" + ", covar_terms_str)
    } else {
      ""
    }
  )
  
  model_formula <- stats::as.formula(formula_str)
  
  # ------------------------------------------------------------
  # 12. Build likelihood
  # ------------------------------------------------------------
  
  if (family == "poisson") {
    
    lik <- inlabru::like(
      family  = "poisson",
      formula = model_formula,
      data    = sp_dat
    )
    
  } else {
    
    lambda <- calibrate_pc_lambda(
      target_prob     = nb_pc_target_prob,
      threshold_theta = nb_pc_threshold_theta
    )
    
    lik <- inlabru::like(
      family  = "nbinomial",
      formula = model_formula,
      data    = sp_dat,
      control.family = list(
        hyper = list(
          theta = list(
            prior = "pc.gamma",
            param = c(lambda)
          )
        )
      )
    )
  }
  
  # ------------------------------------------------------------
  # 13. Set inlabru / INLA options
  # ------------------------------------------------------------
  
  bru_opts <- list(
    inla.mode = inla_mode,
    
    control.inla = list(
      int.strategy = int_strategy,
      strategy     = strategy
    ),
    
    control.compute = list(
      waic = waic,
      cpo  = cpo
    ),
    
    bru_verbose = bru_verbose
  )
  
  # ------------------------------------------------------------
  # 14. Fit model, with optional retries
  # ------------------------------------------------------------
  
  attempt <- 0
  last_err <- NULL
  fit <- NULL
  
  while (attempt <= retry) {
    
    attempt <- attempt + 1
    
    fit <- try(
      inlabru::bru(
        components = model_components,
        lik,
        options = bru_opts
      ),
      silent = TRUE
    )
    
    if (!inherits(fit, "try-error")) {
      break
    }
    
    last_err <- fit
    fit <- NULL
  }
  
  if (is.null(fit)) {
    stop(
      "bru() failed after ", attempt, " attempt(s).\nLast error:\n",
      as.character(last_err)
    )
  }
  
  # Record what the site and square terms actually did, for model_summaries.
  fit$use_site_re     <- use_site_re
  fit$use_square_re   <- use_square_re
  fit$n_site_levels   <- if (use_site_re) n_site else 0L
  fit$n_site_surveys  <- if (use_site_re) sum(is_multi) else 0L
  fit$n_square_levels <- if (use_square_re) n_square_levels else 0L
  fit$tod_fixed_range <- fixed_TOD_range
  fit$tod_knot_gap    <- TOD_knot_gap
  fit$doy_fixed_range <- fixed_DOY_range
  fit$doy_knot_gap    <- DOY_knot_gap
  fit$aru_prior_sd    <- aru_prior_sd   # shared SD for the ARU_1min/3min/5min contrasts
  
  fit
}

# fit_PC_ARU_nosite()  --  fit_PC_ARU() with the site-level iid FORCED OFF
#
# Thin wrapper around fit_PC_ARU(): identical data handling, priors, meshes
# and defaults, except the 100-m site random effect is disabled for every
# species (include_site = FALSE) rather than gated on n_site. Use it as the
# comparison arm when checking (a) whether the site term shifts estimates or
# only interval widths, and (b) how much fit/predict time the term costs.
#
# All other arguments pass straight through via `...`, so this wrapper never
# needs updating when fit_PC_ARU()'s signature changes. Passing include_site
# here is an error -- the whole point of the wrapper is that it is fixed.
fit_PC_ARU_nosite <- function(sp_dat, mesh_abund, mesh_chg, covariates, ...) {
  if ("include_site" %in% names(list(...))) {
    stop("fit_PC_ARU_nosite() always sets include_site = FALSE; do not pass it.")
  }
  fit_PC_ARU(
    sp_dat       = sp_dat,
    mesh_abund   = mesh_abund,
    mesh_chg     = mesh_chg,
    covariates   = covariates,
    include_site = FALSE,
    ...
  )
}

fit_PC_ARU_CL <- function(
    sp_dat,
    mesh_abund,
    mesh_chg,
    covariates,
    mesh_TOD = NULL,
    timeout_min = 120,
    
    # Spatial SPDE priors
    prior_range_abund  = c(300, 0.1),   # P(range < 300 km) = 0.1  -> favours ranges > 300 km
    prior_sigma_abund  = c(0.5, 0.1),   # P(sigma > 0.5) = 0.1     -> prior median ~0.15; caps field amplitude
    prior_range_change = c(300, 0.1),   # P(range < 300 km) = 0.1  -> favours ranges > 300 km
    prior_sigma_change = c(0.3, 0.1),   # P(sigma > 0.3) = 0.1     -> prior median ~0.09; tighter than abundance (change is a noisy difference)
    
    # 1-D SPDE: survey timing. FIXED range + estimated sigma (see fit_PC_ARU).
    # Fixing the range keeps the smoothness of the detectability correction
    # comparable across species and stops the range collapsing for sparse ones.
    fixed_TOD_range = 1.5,                # hours; NOT estimated
    prior_TOD_sigma = c(2, 0.1),          # P(sigma > 2) = 0.1 -> prior median ~0.6
    
    # 1-D SPDE: day of year (seasonal phenology / detectability). Modelled
    # exactly like TOD -- shared mesh (mesh_DOY), FIXED range (default 30 days),
    # estimated sigma.
    mesh_DOY = NULL,
    fixed_DOY_range = 30,                 # days; NOT estimated
    prior_DOY_sigma = c(2, 0.1),          # P(sigma > 2) = 0.1 -> prior median ~0.6

    # ARU protocol effects (ARU_1min / ARU_3min / ARU_5min). Each is a contrast
    # of that ARU protocol against a 5-minute point count on the log link; all
    # three share ONE prior: mean 0, SD = aru_prior_sd. The default log(5)/2
    # (SD ~= 0.8) puts ~95% of prior mass within a factor of ~5 of the PC
    # reference -- wide enough for the ARU-vs-observer method offset plus the
    # shortening from 5 to 1 minute. That shortening is modest, not the naive
    # log(1/5): most birds at a location are detected in the first minute, so a
    # 1-minute window recovers the bulk of a 5-minute count.
    aru_prior_sd = log(5) / 2,

    # Checklist effort effects: log-linear corrections on centred log-effort.
    #   BBA:  beta_BBA * log(Survey_Duration_Minutes / median_BBA_duration)
    #   LT:   beta_LT  * log(Distance_Traveled_m     / median_LT_distance)
    # Each beta is a SINGLE fixed effect with a Gaussian prior (mean, sd on the
    # log link) -- not a 1-D SPDE. Centring the covariate on the protocol's
    # median effort keeps effect_SC / effect_LT as protocol contrasts evaluated
    # at typical checklist effort. Prior mean > 0 favours more detections with
    # more effort (with diminishing returns, since it is linear in log-effort).
    BBA_log_duration_prior_mean = 0.30,
    BBA_log_duration_prior_sd   = 0.20,
    LT_log_distance_prior_mean  = 0.30,
    LT_log_distance_prior_sd    = 0.20,

    # --- Special surveys (targeted single-protocol surveys) --------------------
    # Survey types given their OWN intercept (a protocol detectability offset vs
    # the 5-min point-count reference) and NO effort covariate, because each is
    # standardized within its own protocol (e.g. the marshbird quiet+playback
    # sequence). They still load onto the shared abundance/change surface, habitat
    # covariates and the square/site iid, so they inform distribution/abundance
    # like any other survey. Leave as character(0) for the plain checklist model.
    special_survey_types  = character(0),
    # A special type's intercept is fit only where it is identifiable: the species
    # must be detected in >= this many atlas squares (unique square_id, count > 0)
    # by that survey type. Types below threshold have ALL their rows dropped -- no
    # intercept, and no contamination of the reference protocol. Either one number
    # for every type, or a named vector to override per type; types not named fall
    # back to a "default" entry if present, else the first entry, e.g.
    #   c(default = 10, "Marshbird Survey" = 20)
    min_special_squares   = 10,
    # SD (log link) of each special-survey intercept's Gaussian prior; matches the
    # factor-of-5 detectability prior used for the checklist intercepts.
    special_intercept_prior_sd = log(5) / 2,
    # Whether special-survey rows share the diurnal TOD / seasonal DOY smooths.
    # Either one logical for every special type, or a NAMED logical vector keyed
    # by survey type (e.g. c("Northern Hawk Owl Survey" = TRUE)); types not named
    # fall back to a "default" entry if present, else FALSE. FALSE (the default)
    # treats timing/season as fixed by protocol and absorbed by the intercept, and
    # keeps nocturnal owl/nightjar rows from bending the PC/ARU-based smooths where
    # those carry no data. Per-type control lets, say, the near-dawn hawk-owl
    # survey opt into TOD without pulling the deep-night owl surveys in with it.
    special_share_tod     = FALSE,
    special_share_doy     = FALSE,

    # Atlas-square iid prior
    kappa_pcprec_diff = c(0.25, 0.1),    # P(sigma_square > 0.25) = 0.1; fine-scale (10-km square) sink
    
    # include_square: master switch for the atlas-square iid random effect
    #   (kappa_diff). TRUE keeps it; FALSE forces it OFF. Symmetric to
    #   include_site, so the two effects toggle independently.
    include_square    = TRUE,
    
    # Site-level iid prior and inclusion rules. ONE site term, SHARED by every
    # survey type (see 3b) -- there is no separate linear-transect site effect.
    # min_site_surveys: minimum RETAINED surveys (any survey type) for a site to
    #   get a level.
    # min_site_levels:  minimum number of such sites, else the term is dropped
    #   (sigma_site cannot be identified from a handful of sites).
    prior_site_pcprec = c(0.25, 0.1),    # P(sigma_site > 0.25) = 0.1; site sink (kept tight; low counts barely inform it)
    
    # include_site: master switch for the site-level iid random effect.
    #   TRUE  -> keep the data-driven gate below (the term is used only
    #            where it is identifiable: n_site >= min_site_levels).
    #   FALSE -> force the term OFF for every species, regardless of how
    #            many multi-visit sites exist.
    include_site      = TRUE,
    
    min_site_surveys  = 2,
    min_site_levels   = 10,
    
    # Error family (NO fallback -- a failure is raised). ONE family, shared by
    # every survey type (PC + ARU + BBA + LT), and when nbinomial ONE size -- a
    # single observation model, far easier to identify than a per-stream split.
    # Only the effort terms differ by stream; the family and the site-level
    # random effect are shared (3b).
    family = c("nbinomial", "poisson"),
    
    # Negative-binomial PC prior settings (shared across all surveys).
    nb_pc_target_prob     = 0.5,
    nb_pc_threshold_theta = 5,
    
    # inlabru / INLA options
    inla_mode    = "experimental",
    int_strategy = "ccd",
    strategy     = "simplified.laplace",
    bru_verbose  = 4,
    waic         = FALSE,
    cpo          = FALSE,
    retry        = 0
) {
  
  # ------------------------------------------------------------
  # 1. Validate inputs
  # ------------------------------------------------------------
  
  family <- match.arg(family)
  
  if (!inherits(sp_dat, "sf")) {
    stop("sp_dat must be an sf object with a geometry column.")
  }
  
  if (!inherits(study_boundary, "sf")) {
    stop("study_boundary must be an sf object.")
  }
  
  required_cols <- c(
    "count",
    "Atlas3",
    "Atlas3_c",
    "Hours_After_Reference",
    "days_midpoint",
    "ARU",
    "ARU_1min",
    "ARU_3min",
    "ARU_5min",
    "SC",
    "LT",
    "Survey_Type",
    "Survey_Duration_Minutes",
    "Distance_Traveled_m",
    "square_atlas",
    "site_id",
    "geometry"
  )
  
  missing_cols <- setdiff(required_cols, names(sp_dat))
  
  if (length(missing_cols) > 0) {
    stop(
      "sp_dat is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  if (!is.data.frame(covariates)) {
    stop("covariates must be a data.frame.")
  }
  
  needed_cov_cols <- c("covariate", "beta", "model", "mean", "prec")
  cov_missing <- setdiff(needed_cov_cols, names(covariates))
  
  if (length(cov_missing) > 0) {
    stop(
      "covariates is missing required columns: ",
      paste(cov_missing, collapse = ", ")
    )
  }
  
  stopifnot(
    is.numeric(fixed_TOD_range), length(fixed_TOD_range) == 1,
    fixed_TOD_range > 0,
    is.numeric(BBA_log_duration_prior_mean), length(BBA_log_duration_prior_mean) == 1,
    is.numeric(BBA_log_duration_prior_sd),   length(BBA_log_duration_prior_sd)   == 1,
    BBA_log_duration_prior_sd > 0,
    is.numeric(LT_log_distance_prior_mean),  length(LT_log_distance_prior_mean)  == 1,
    is.numeric(LT_log_distance_prior_sd),    length(LT_log_distance_prior_sd)    == 1,
    LT_log_distance_prior_sd > 0
  )
  
  # ------------------------------------------------------------
  # 2. Set INLA options and align CRS
  # ------------------------------------------------------------
  
  INLA::inla.setOption(inla.timeout = 60 * timeout_min)
  
  if (sf::st_crs(sp_dat) != sf::st_crs(study_boundary)) {
    study_boundary <- sf::st_transform(study_boundary, sf::st_crs(sp_dat))
  }
  
  # ------------------------------------------------------------
  # 3. Ensure model variables have expected types
  # ------------------------------------------------------------
  
  sp_dat$Hours_After_Reference   <- ensure_numeric(sp_dat$Hours_After_Reference)
  sp_dat$days_midpoint           <- ensure_numeric(sp_dat$days_midpoint)
  sp_dat$Atlas3_c                <- ensure_numeric(sp_dat$Atlas3_c)
  sp_dat$ARU                     <- ensure_numeric(sp_dat$ARU)
  sp_dat$ARU_1min                <- ensure_numeric(sp_dat$ARU_1min)
  sp_dat$ARU_3min                <- ensure_numeric(sp_dat$ARU_3min)
  sp_dat$ARU_5min                <- ensure_numeric(sp_dat$ARU_5min)
  sp_dat$SC                      <- ensure_numeric(sp_dat$SC)
  sp_dat$LT                      <- ensure_numeric(sp_dat$LT)
  sp_dat$Survey_Duration_Minutes <- ensure_numeric(sp_dat$Survey_Duration_Minutes)
  sp_dat$Distance_Traveled_m     <- ensure_numeric(sp_dat$Distance_Traveled_m)

  # Special-survey rows are, by definition, neither ARU nor BBA (SC) nor LT, so
  # force those protocol indicators to 0 on them. This keeps the ARU partition
  # check below robust (a stray NA/1 would break it) and guarantees a special row
  # never loads onto the BBA/LT intercepts -- its own intercept is added in 3a.
  if (length(special_survey_types) > 0) {
    is_special_raw <- sp_dat$Survey_Type %in% special_survey_types
    if (any(is_special_raw)) {
      for (col in c("ARU", "ARU_1min", "ARU_3min", "ARU_5min", "SC", "LT")) {
        sp_dat[[col]][is_special_raw] <- 0
      }
    }
  }

  # The three duration indicators must PARTITION the ARU rows: every ARU row
  # carries exactly one (duration in {1, 3, 5}) and every non-ARU row carries
  # none. This is what makes each effect_ARU_* a clean contrast against the
  # 5-minute point-count reference; a row with all-zero indicators would fold
  # silently back into that reference.
  stopifnot(all(sp_dat$ARU_1min + sp_dat$ARU_3min + sp_dat$ARU_5min == sp_dat$ARU))
  
  # ----------------------------------------------------------------------------
  # 3a. Resolve special surveys, then build survey-type masks
  # ----------------------------------------------------------------------------
  # Special (targeted single-protocol) surveys each get their own intercept and
  # no effort term (sections 10/11). A type's intercept is fit only where the
  # species is detected in >= its square threshold; types below threshold have
  # ALL their rows dropped here, so they neither get an intercept nor contaminate
  # the reference protocol. Row-dropping happens BEFORE the masks, the site index
  # and the effort medians below, all of which must see the final row set.
  if (length(special_survey_types) > 0 && !"square_id" %in% names(sp_dat)) {
    stop("fit_PC_ARU_CL(): square_id is required to threshold special surveys.")
  }

  special_types_req <- intersect(special_survey_types, unique(sp_dat$Survey_Type))

  special_keep   <- character(0)
  special_report <- list()
  for (type in special_types_req) {
    det  <- sp_dat$Survey_Type == type & sp_dat$count > 0
    n_sq <- dplyr::n_distinct(sp_dat$square_id[det])
    thr  <- resolve_special_threshold(type, min_special_squares)
    keep <- n_sq >= thr
    if (keep) special_keep <- c(special_keep, type)
    special_report[[type]] <- c(n_det_squares = n_sq, threshold = thr,
                                kept = as.numeric(keep))
    message(sprintf(
      "  special survey '%s': detected in %d square(s) (threshold %g) -> %s",
      type, n_sq, thr, if (keep) "KEEP intercept" else "DROP rows"))
  }

  special_drop <- setdiff(special_types_req, special_keep)
  if (length(special_drop) > 0) {
    sp_dat <- sp_dat[!(sp_dat$Survey_Type %in% special_drop), , drop = FALSE]
  }

  # PC + ARU (structured) and BBA (stationary checklist) share one predictor; LT
  # and each KEPT special type add their own intercept. All share the site-level
  # random effect built in 3b. Masks are computed on the FINAL row set.
  is_pc_aru     <- sp_dat$Survey_Type %in% c("Point_Count", "ARU")
  is_bba        <- sp_dat$Survey_Type == "Breeding Bird Atlas"
  is_lt         <- sp_dat$Survey_Type == "Linear transect"
  is_special    <- sp_dat$Survey_Type %in% special_keep
  is_stationary <- is_pc_aru | is_bba
  
  n_pc_aru <- sum(is_pc_aru)
  n_bba    <- sum(is_bba)
  n_lt     <- sum(is_lt)
  
  if (n_pc_aru == 0) {
    stop("fit_PC_ARU_CL(): no Point_Count/ARU rows to fit.")
  }
  if (any(!(is_pc_aru | is_bba | is_lt | is_special))) {
    bad <- sort(unique(sp_dat$Survey_Type[!(is_pc_aru | is_bba | is_lt | is_special)]))
    stop("fit_PC_ARU_CL(): rows with unrecognised Survey_Type (",
         paste(bad, collapse = ", "),
         "). Expected Point_Count, ARU, Breeding Bird Atlas, Linear transect, ",
         "or a declared special_survey_types value.")
  }
  
  # ----------------------------------------------------------------------------
  # 3b. Site index -- ONE term, SHARED by every survey type
  # ----------------------------------------------------------------------------
  # Keyed on site_id (the site-year grid built in script 07) across ALL rows, so
  # linear transects load onto the SAME site nodes as PC/ARU/BBA: repeat visits
  # to a place are repeat visits whatever the protocol, and one shared variance
  # is far better identified than two. The error family is shared across streams
  # too; only the effort terms in each predictor differ.
  #
  # Computed HERE, not upstream: sp_dat has already been through species-specific
  # safe-date and survey-type filtering, so a site with repeats in the raw data
  # may be a singleton for THIS species. Only sites with >= min_site_surveys
  # RETAINED surveys -- counted across all survey types -- get a level.
  #
  # Singletons are excluded because they carry no within-site correlation to
  # absorb; an iid level informed by a single observation merely duplicates the
  # nbinomial overdispersion parameter. Division of labour:
  #     nbinomial  -> observation-level overdispersion, everywhere
  #     site_iid   -> shared correlation, only where repeats exist
  #
  # Singleton rows are pointed at level 1 and carry site_w = 0, so
  # `site_w * site_iid` gives them a zero Jacobian row: they contribute nothing
  # to any site node.
  # ----------------------------------------------------------------------------
  
  site_tab    <- table(sp_dat$site_id)
  multi_sites <- names(site_tab)[site_tab >= min_site_surveys]
  is_multi    <- !is.na(sp_dat$site_id) & (sp_dat$site_id %in% multi_sites)
  
  sp_dat$site_idx_safe <- 1L
  sp_dat$site_idx_safe[is_multi] <-
    as.integer(factor(sp_dat$site_id[is_multi]))
  sp_dat$site_w <- as.numeric(is_multi)
  
  n_site       <- length(multi_sites)
  enough_sites <- n_site >= min_site_levels
  
  # The term is used only when the caller ALLOWED it (include_site) AND it is
  # identifiable (enough multi-visit sites).
  use_site_re <- isTRUE(include_site) && enough_sites
  
  site_re_note <- if (!isTRUE(include_site)) {
    "  -- TERM DISABLED (include_site = FALSE)"
  } else if (!enough_sites) {
    "  -- TERM DROPPED (too few sites)"
  } else {
    ""
  }
  
  message(
    "  site RE (shared by all survey types): ", n_site, " sites with >= ",
    min_site_surveys, " surveys, covering ", sum(is_multi), "/", nrow(sp_dat),
    " surveys (", sum(is_multi & is_stationary), " stationary, ",
    sum(is_multi & is_lt), " LT)", site_re_note
  )
  
  # ----------------------------------------------------------------------------
  # Square-level iid switch (kappa_diff)
  # ----------------------------------------------------------------------------
  # No data-driven gate: atlas squares are always plentiful, so the master
  # switch is the only control. n_square_levels is recorded on the fit object.
  use_square_re   <- isTRUE(include_square)
  n_square_levels <- length(unique(sp_dat$square_atlas))
  message(
    "  square RE: ", n_square_levels, " atlas-square levels",
    if (!use_square_re) "  -- TERM DISABLED (include_square = FALSE)" else ""
  )
  
  # ----------------------------------------------------------------------------
  # 3c. Two observation models
  # ----------------------------------------------------------------------------
  # The stationary model (PC + ARU + BBA) and the LT model each get their OWN
  # likelihood -- and, when nbinomial, their OWN size -- so their very different
  # overdispersion is never forced onto a shared parameter. Both share the SAME
  # latent predictor, so abundance/change, the site effect and every other latent
  # effect are informed by all rows.
  message(
    "  sample sizes: ", n_pc_aru + n_bba, " stationary (",
    n_pc_aru, " PC/ARU + ", n_bba, " BBA)",
    if (n_lt > 0) paste0(" + ", n_lt, " LT") else "",
    ""
  )
  
  # ------------------------------------------------------------
  # 4. Create centred log-effort covariates for the effort corrections
  # ------------------------------------------------------------
  # Each checklist effort correction is a single log-linear term:
  #   BBA: beta_BBA * log(Survey_Duration_Minutes / bba_ref_duration)
  #   LT:  beta_LT  * log(Distance_Traveled_m     / lt_ref_distance)
  # The reference is the protocol's MEDIAN effort, so the covariate is 0 at
  # typical effort and effect_SC / effect_LT stay protocol contrasts there.
  # Non-applicable rows are set to EXACTLY 0, so `<var> * effect_<var>` drops out
  # of every other protocol without needing a survey-indicator multiplier.
  
  sp_dat$BBA <- as.numeric(sp_dat$Survey_Type == "Breeding Bird Atlas")
  
  # BBA reference duration (median over BBA rows). Log-effort needs finite,
  # strictly-positive durations on BBA rows.
  if (n_bba > 0) {
    bba_ref_duration <- stats::median(
      sp_dat$Survey_Duration_Minutes[is_bba],
      na.rm = TRUE
    )
    if (!is.finite(bba_ref_duration)) {
      stop("Could not determine a finite BBA reference duration.")
    }
    if (any(!is.finite(sp_dat$Survey_Duration_Minutes[is_bba]))) {
      stop("All BBA Survey_Duration_Minutes values must be finite.")
    }
    if (any(sp_dat$Survey_Duration_Minutes[is_bba] <= 0)) {
      stop("BBA Survey_Duration_Minutes must be > 0 for log-effort correction.")
    }
  } else {
    bba_ref_duration <- NA_real_
  }
  
  # LT reference distance (median over LT rows), same positivity requirement.
  if (n_lt > 0) {
    lt_ref_distance <- stats::median(
      sp_dat$Distance_Traveled_m[is_lt],
      na.rm = TRUE
    )
    if (!is.finite(lt_ref_distance)) {
      stop("Could not determine a finite linear-transect reference distance.")
    }
    if (any(!is.finite(sp_dat$Distance_Traveled_m[is_lt]))) {
      stop("All linear-transect Distance_Traveled_m values must be finite.")
    }
    if (any(sp_dat$Distance_Traveled_m[is_lt] <= 0)) {
      stop("Linear-transect Distance_Traveled_m must be > 0 for log-effort correction.")
    }
  } else {
    lt_ref_distance <- NA_real_
  }
  
  sp_dat$BBA_log_duration <- 0
  if (n_bba > 0) {
    sp_dat$BBA_log_duration[is_bba] <-
      log(sp_dat$Survey_Duration_Minutes[is_bba] / bba_ref_duration)
  }
  
  sp_dat$LT_log_distance <- 0
  if (n_lt > 0) {
    sp_dat$LT_log_distance[is_lt] <-
      log(sp_dat$Distance_Traveled_m[is_lt] / lt_ref_distance)
  }

  # ------------------------------------------------------------
  # 4b. Special-survey indicators and TOD/DOY participation weights
  # ------------------------------------------------------------
  # One 0/1 indicator per surviving special type (the column name is reused as
  # the effect name, so a fit carries readable effect_SS_<slug> terms). tod_w /
  # doy_w switch the diurnal TOD / seasonal DOY smooths on or off PER TYPE: a
  # special type participates only if special_share_tod / special_share_doy opts
  # it in (see the resolver below); every non-special row always participates.
  # This is per type, not global, so a near-dawn protocol can share the TOD smooth
  # while deep-night ones stay out. When a special type is excluded, its rows'
  # timing/season covariate is pinned to an in-domain reference so the 1-D
  # projector never has to evaluate an out-of-mesh (e.g. nocturnal) value.
  special_slug <- stats::setNames(paste0("SS_", sp_filename(special_keep)),
                                  special_keep)
  for (type in special_keep) {
    sp_dat[[special_slug[[type]]]] <- as.numeric(sp_dat$Survey_Type == type)
  }

  # Resolve a scalar-or-named logical flag for one special type. Mirrors
  # resolve_special_threshold(): named types use their value; unlisted types use a
  # "default" entry if present, else FALSE.
  special_flag <- function(type, flags) {
    if (is.null(names(flags)) || all(names(flags) == "")) {
      return(isTRUE(unname(flags)[1]))
    }
    if (type %in% names(flags))      return(isTRUE(flags[[type]]))
    if ("default" %in% names(flags)) return(isTRUE(flags[["default"]]))
    FALSE
  }

  share_tod_type <- vapply(special_keep, special_flag, logical(1),
                           flags = special_share_tod)
  share_doy_type <- vapply(special_keep, special_flag, logical(1),
                           flags = special_share_doy)
  names(share_tod_type) <- special_keep
  names(share_doy_type) <- special_keep

  tod_off_types <- special_keep[!share_tod_type]
  doy_off_types <- special_keep[!share_doy_type]

  # A row is gated out of a smooth iff its (special) type opted out; all PC/ARU/
  # BBA/LT rows keep weight 1.
  sp_dat$tod_w <- as.numeric(!(sp_dat$Survey_Type %in% tod_off_types))
  sp_dat$doy_w <- as.numeric(!(sp_dat$Survey_Type %in% doy_off_types))
  if (length(tod_off_types) > 0) {
    sp_dat$Hours_After_Reference[sp_dat$Survey_Type %in% tod_off_types] <-
      stats::median(sp_dat$Hours_After_Reference[is_pc_aru], na.rm = TRUE)
  }
  if (length(doy_off_types) > 0) {
    sp_dat$days_midpoint[sp_dat$Survey_Type %in% doy_off_types] <- 0
  }
  if (length(special_keep) > 0) {
    message("  special surveys share smooths: TOD = {",
            paste(special_keep[share_tod_type], collapse = ", "), "}, DOY = {",
            paste(special_keep[share_doy_type], collapse = ", "), "}")
  }

  # ------------------------------------------------------------
  # 6. Build spatial SPDE models
  # ------------------------------------------------------------
  
  matern_mean <- INLA::inla.spde2.pcmatern(
    mesh        = mesh_abund,
    prior.range = prior_range_abund,
    prior.sigma = prior_sigma_abund,
    constr      = TRUE
  )
  
  matern_diff <- INLA::inla.spde2.pcmatern(
    mesh        = mesh_chg,
    prior.range = prior_range_change,
    prior.sigma = prior_sigma_change,
    constr      = TRUE
  )
  
  # ------------------------------------------------------------
  # 7. Build 1-D SPDE smoothers for timing and effort corrections
  # ------------------------------------------------------------
  
  # --- Time-of-day smooth: SHARED mesh, FIXED range ---------------------------
  # The mesh is built once (make_TOD_mesh) and passed in, so every species uses
  # the same basis and the same absolute range means the same thing everywhere.
  #
  # prior.range = c(range0, NA) -> range HELD FIXED at range0 (no prior).
  # Only sigma is estimated. Same idiom as the effort SPDEs.
  if (is.null(mesh_TOD)) {
    mesh_TOD <- make_TOD_mesh(knot_spacing = fixed_TOD_range / 6,
                              pad_range    = fixed_TOD_range)
  }
  
  # Guard the requirement that makes the fixed range meaningful: a basis coarser
  # than range/5 cannot represent the correlation length being asked for, and
  # the mesh, not the range, then sets the smoothness.
  TOD_knot_gap <- max(diff(sort(mesh_TOD$loc)))
  if (TOD_knot_gap > fixed_TOD_range / 5) {
    warning(
      "TOD mesh knot spacing (", signif(TOD_knot_gap, 3),
      " h) is coarser than fixed_TOD_range/5 (",
      signif(fixed_TOD_range / 5, 3),
      " h). The mesh, not the range, is setting the smoothness.",
      call. = FALSE
    )
  }
  
  TOD_spde <- INLA::inla.spde2.pcmatern(
    mesh        = mesh_TOD,
    alpha       = 2,
    prior.range = c(fixed_TOD_range, NA),
    prior.sigma = prior_TOD_sigma,
    constr      = TRUE
  )
  
  # --- Day-of-year smooth: SHARED mesh, FIXED range --------------------------
  # Modelled exactly like the time-of-day smooth above, but the domain is NOT
  # cyclic -- days_midpoint is a linear axis centred on the safe-date midpoint.
  if (is.null(mesh_DOY)) {
    mesh_DOY <- make_DOY_mesh(knot_spacing = fixed_DOY_range / 6,
                              pad_range    = fixed_DOY_range)
  }
  
  # Same knot-spacing guard as for TOD.
  DOY_knot_gap <- max(diff(sort(mesh_DOY$loc)))
  if (DOY_knot_gap > fixed_DOY_range / 5) {
    warning(
      "DOY mesh knot spacing (", signif(DOY_knot_gap, 3),
      " d) is coarser than fixed_DOY_range/5 (",
      signif(fixed_DOY_range / 5, 3),
      " d). The mesh, not the range, is setting the smoothness.",
      call. = FALSE
    )
  }
  
  DOY_spde <- INLA::inla.spde2.pcmatern(
    mesh        = mesh_DOY,
    alpha       = 2,
    prior.range = c(fixed_DOY_range, NA),
    prior.sigma = prior_DOY_sigma,
    constr      = TRUE
  )
  
  # --- Checklist effort corrections: log-linear (no 1-D SPDE) -----------------
  # The effort effects are single fixed effects on the centred log-effort
  # covariates built in section 4; their Gaussian priors are injected as literals
  # into the components below. Nothing mesh-based is needed here.
  message(
    "  effort log-linear: BBA ref duration = ",
    if (n_bba > 0) round(bba_ref_duration, 1) else NA, " min",
    ", LT ref distance = ",
    if (n_lt > 0) round(lt_ref_distance, 1) else NA, " m"
  )
  
  # ------------------------------------------------------------
  # 8. Define iid priors
  # ------------------------------------------------------------
  
  pc_prec_diff <- list(
    prior = "pcprec",
    param = kappa_pcprec_diff
  )
  
  pc_prec_site <- list(
    prior = "pcprec",
    param = prior_site_pcprec
  )
  
  # ------------------------------------------------------------
  # 9. Build covariate components and formula terms
  # ------------------------------------------------------------
  
  covariates <- covariates %>%
    dplyr::mutate(
      component = paste0(
        "Beta", beta, "_", covariate,
        '(1, model="', model,
        '", mean.linear=', mean,
        ", prec.linear=", prec,
        ")"
      ),
      term = paste0(
        "Beta", beta, "_", covariate,
        " * I(", covariate, "^", beta, ")"
      )
    )
  
  covar_components_str <- if (nrow(covariates) > 0) {
    paste(covariates$component, collapse = " + ")
  } else {
    ""
  }
  
  covar_terms_str <- if (nrow(covariates) > 0) {
    paste(covariates$term, collapse = " + ")
  } else {
    ""
  }
  
  # ------------------------------------------------------------
  # 10. Define inlabru components
  # ------------------------------------------------------------
  
  # Shared ARU protocol-effect prior (ARU_1min / ARU_3min / ARU_5min): mean 0,
  # SD = aru_prior_sd on the log link (see argument). The precision is injected
  # as a literal below so the fitted formula carries no environment dependency.
  aru_prec <- 1 / (aru_prior_sd^2)
  
  components_str <- paste0(
    "~ Intercept(1) + ",
    
    "effect_Atlas3(
       1,
       model = 'linear',
       mean.linear = 0,
       prec.linear = 1 / ((log(1.5) / 2)^2)
     ) + ",
    
    "effect_ARU_1min(
       1,
       model = 'linear',
       mean.linear = 0,
       prec.linear = ", aru_prec, "
     ) + ",
    
    "effect_ARU_3min(
       1,
       model = 'linear',
       mean.linear = 0,
       prec.linear = ", aru_prec, "
     ) + ",
    
    "effect_ARU_5min(
       1,
       model = 'linear',
       mean.linear = 0,
       prec.linear = ", aru_prec, "
     ) + ",
    
    "spde_mean(main = geometry, model = matern_mean) + ",
    "spde_diff(main = geometry, model = matern_diff) + ",
    
    "TOD_global(main = Hours_After_Reference, model = TOD_spde) + ",
    "DOY_global(main = days_midpoint, model = DOY_spde)",
    
    # BBA (stationary checklist) intercept + log-duration effort coefficient,
    # defined ONLY when BBA rows exist -- i.e. iff those terms enter a predictor
    # below -- so no component is ever left referenced by zero likelihoods. The
    # effort effect is a single linear coefficient (multiplied by BBA_log_duration
    # in the formula), NOT a 1-D SPDE.
    if (n_bba > 0) {
      paste0(
        " + effect_SC(
               1,
               model = 'linear',
               mean.linear = 0,
               prec.linear = 1 / ((log(5) / 2)^2)
             ) +
            effect_BBA_log_duration(
               1,
               model = 'linear',
               mean.linear = ", BBA_log_duration_prior_mean, ",
               prec.linear = ", 1 / (BBA_log_duration_prior_sd^2), "
             )"
      )
    } else {
      ""
    },
    
    # LT (linear transect) intercept + log-distance effort coefficient, defined
    # ONLY when LT rows exist (iff the LT likelihood is built below).
    if (n_lt > 0) {
      paste0(
        " + effect_LT(
               1,
               model = 'linear',
               mean.linear = 0,
               prec.linear = 1 / ((log(5) / 2)^2)
             ) +
            effect_LT_log_distance(
               1,
               model = 'linear',
               mean.linear = ", LT_log_distance_prior_mean, ",
               prec.linear = ", 1 / (LT_log_distance_prior_sd^2), "
             )"
      )
    } else {
      ""
    },

    # Special-survey intercepts: one protocol offset per surviving special type,
    # NO effort term. Built only for types that cleared their square threshold, so
    # no component is ever left referenced by zero likelihood rows.
    if (length(special_keep) > 0) {
      paste0(" + ", paste(sprintf(
        "effect_%s(1, model = 'linear', mean.linear = 0, prec.linear = %s)",
        special_slug[special_keep],
        1 / (special_intercept_prior_sd^2)
      ), collapse = " + "))
    } else {
      ""
    },
    
    # Atlas-square iid, only when include_square = TRUE.
    if (use_square_re) {
      "+ kappa_diff(
         square_atlas,
         model = 'iid',
         constr = TRUE,
         hyper = list(prec = pc_prec_diff)
       )"
    } else {
      ""
    },
    
    # Site-level iid: ONE component, entered in the shared predictor below, so
    # every survey type draws on the same site nodes and the same variance.
    # constr = TRUE (as for kappa_diff) makes the site effects sum to zero, so
    # dropping the term at prediction predicts the AVERAGE site exactly rather
    # than leaving the intercept and the site mean weakly confounded.
    if (use_site_re) {
      "+ site_iid(
         site_idx_safe,
         model = 'iid',
         n = n_site,
         constr = TRUE,
         hyper = list(prec = pc_prec_site)
       )"
    } else {
      ""
    },
    
    if (nchar(covar_components_str) > 0) {
      paste0(" + ", covar_components_str)
    } else {
      ""
    }
  )
  
  model_components <- stats::as.formula(components_str)
  
  # ------------------------------------------------------------
  # 11. Define linear predictor
  # ------------------------------------------------------------
  
  # Terms shared by every survey type: the latent abundance/change surface,
  # shared detectability (TOD, DOY), the square and SITE random effects, and the
  # covariates. The site term lives here, which is what makes LT and stationary
  # rows share one site effect.
  shared_terms <- paste0(
    "count ~
      Intercept +

      spde_mean +

      Atlas3_c * spde_diff +
      Atlas3_c * effect_Atlas3",
    # Special rows are gated out of the diurnal/seasonal smooths (via tod_w/doy_w)
    # unless the caller opted them in; with no special surveys these are bare
    # TOD_global + DOY_global, identical to the plain checklist model.
    if (length(special_keep) > 0)
         " + tod_w * TOD_global + doy_w * DOY_global"
    else " + TOD_global + DOY_global",
    if (use_square_re) " + kappa_diff" else "",
    if (use_site_re)   " + site_w * site_iid" else "",
    if (nchar(covar_terms_str) > 0) paste0(" + ", covar_terms_str) else ""
  )
  
  # ONE predictor over all rows. Each protocol/effort term is gated by its own
  # indicator (ARU_*, BBA/SC, LT), so it is exactly zero on rows of the other
  # protocols; every survey type therefore shares this single mean structure.
  # Terms for an absent protocol are dropped entirely so no zero-weight projector
  # is built for a component with no active rows.
  bba_terms <- if (n_bba > 0) {
    " +
      SC  * effect_SC +
      BBA_log_duration * effect_BBA_log_duration"
  } else {
    ""
  }
  
  lt_terms <- if (n_lt > 0) {
    " +
      LT  * effect_LT +
      LT_log_distance * effect_LT_log_distance"
  } else {
    ""
  }

  # One intercept term per surviving special type, gated by its own indicator so
  # it is exactly zero on every other protocol's rows (mirrors SC/LT). No effort
  # term accompanies it.
  special_terms <- if (length(special_keep) > 0) {
    paste0(" +\n      ",
           paste(sprintf("%s * effect_%s",
                         special_slug[special_keep],
                         special_slug[special_keep]),
                 collapse = " +\n      "))
  } else {
    ""
  }

  model_formula <- stats::as.formula(paste0(
    shared_terms,
    " +
      ARU_1min * effect_ARU_1min +
      ARU_3min * effect_ARU_3min +
      ARU_5min * effect_ARU_5min",
    bba_terms,
    lt_terms,
    special_terms
  ))
  
  # ------------------------------------------------------------
  # 12. Build likelihood
  # ------------------------------------------------------------
  # ONE observation model over all rows (PC + ARU + BBA + LT): a single shared
  # family and, when nbinomial, a single size. Mirrors fit_PC_ARU().
  
  if (family == "poisson") {
    
    lik <- inlabru::like(
      family  = "poisson",
      formula = model_formula,
      data    = sp_dat
    )
    
  } else {
    
    lambda <- calibrate_pc_lambda(
      target_prob     = nb_pc_target_prob,
      threshold_theta = nb_pc_threshold_theta
    )
    
    lik <- inlabru::like(
      family         = "nbinomial",
      formula        = model_formula,
      data           = sp_dat,
      control.family = list(
        hyper = list(
          theta = list(
            prior = "pc.gamma",
            param = c(lambda)
          )
        )
      )
    )
  }
  
  # ------------------------------------------------------------
  # 13. Set inlabru / INLA options
  # ------------------------------------------------------------
  
  bru_opts <- list(
    inla.mode = inla_mode,
    
    control.inla = list(
      int.strategy = int_strategy,
      strategy     = strategy
    ),
    
    control.compute = list(
      waic = waic,
      cpo  = cpo
    ),
    
    bru_verbose = bru_verbose
  )
  
  # ------------------------------------------------------------
  # 14. Fit model, with optional retries
  # ------------------------------------------------------------
  
  attempt <- 0
  last_err <- NULL
  fit <- NULL
  
  while (attempt <= retry) {
    
    attempt <- attempt + 1
    
    fit <- try(
      inlabru::bru(
        components = model_components,
        lik,
        options = bru_opts
      ),
      silent = TRUE
    )
    
    if (!inherits(fit, "try-error")) {
      break
    }
    
    last_err <- fit
    fit <- NULL
  }
  
  if (is.null(fit)) {
    stop(
      "bru() failed after ", attempt, " attempt(s).\nLast error:\n",
      as.character(last_err)
    )
  }
  
  # Record what the site/square terms and fixed effort ranges did, for summaries.
  # site_re_shared flags that ONE site effect spans all survey types (there is no
  # separate LT site term), so downstream code reads a single sigma_site.
  fit$use_site_re        <- use_site_re
  fit$site_re_shared     <- TRUE
  fit$use_square_re      <- use_square_re
  fit$n_site_levels      <- if (use_site_re) n_site else 0L
  fit$n_site_surveys     <- if (use_site_re) sum(is_multi) else 0L
  fit$n_site_surveys_stationary <- if (use_site_re) sum(is_multi & is_stationary) else 0L
  fit$n_site_surveys_lt  <- if (use_site_re) sum(is_multi & is_lt) else 0L
  fit$n_square_levels    <- if (use_square_re) n_square_levels else 0L
  fit$tod_fixed_range    <- fixed_TOD_range
  fit$tod_knot_gap       <- TOD_knot_gap
  fit$doy_fixed_range    <- fixed_DOY_range
  fit$doy_knot_gap       <- DOY_knot_gap
  fit$bba_ref_duration   <- bba_ref_duration
  fit$lt_ref_distance    <- lt_ref_distance
  fit$bba_log_duration_prior <- c(mean = BBA_log_duration_prior_mean,
                                  sd   = BBA_log_duration_prior_sd)
  fit$lt_log_distance_prior  <- c(mean = LT_log_distance_prior_mean,
                                  sd   = LT_log_distance_prior_sd)
  fit$aru_prior_sd       <- aru_prior_sd   # shared SD for the ARU_1min/3min/5min contrasts
  
  # One shared family/size across all survey types; record it and the per-stream
  # observation counts (data composition, not a family split).
  fit$family            <- family
  fit$n_stationary_obs  <- n_pc_aru + n_bba
  fit$n_pc_aru_obs      <- n_pc_aru
  fit$n_bba_obs         <- n_bba
  fit$n_lt_obs          <- n_lt

  # Special-survey structure: which targeted protocols kept an intercept, the
  # per-type detection-square audit (n_det_squares / threshold / kept), whether
  # they shared the TOD/DOY smooths, and how many rows they contributed.
  fit$special_survey_types       <- special_keep
  fit$special_survey_report      <- special_report
  fit$special_share_tod          <- share_tod_type
  fit$special_share_doy          <- share_doy_type
  fit$special_intercept_prior_sd <- special_intercept_prior_sd
  fit$n_special_obs              <- sum(is_special)
  
  fit
}

# fit_PC_ARU_CL_nosite()  --  fit_PC_ARU_CL() with the site-level iid FORCED OFF
#
# Thin wrapper around fit_PC_ARU_CL(): identical checklist/effort structure,
# priors, meshes and defaults, except the 100-m site random effect is disabled
# for every species (include_site = FALSE) rather than gated on n_site. This
# is the arm to compare against fit_PC_ARU_CL() when deciding whether the site
# term earns its (substantial) compute cost for the checklist model.
#
# All other arguments pass through via `...`; passing include_site is an error.
fit_PC_ARU_CL_nosite <- function(sp_dat, mesh_abund, mesh_chg, covariates, ...) {
  if ("include_site" %in% names(list(...))) {
    stop("fit_PC_ARU_CL_nosite() always sets include_site = FALSE; do not pass it.")
  }
  fit_PC_ARU_CL(
    sp_dat       = sp_dat,
    mesh_abund   = mesh_abund,
    mesh_chg     = mesh_chg,
    covariates   = covariates,
    include_site = FALSE,
    ...
  )
}

# ============================================================
# 4. Prediction
# ============================================================

# Build a prediction expression/formula for multi-atlas outputs
#
# Produces the formula/expression used in `predict_inla()` to compute
# posterior draws for quantities of interest (e.g., mean abundance per atlas,
# change, relative change), optionally including random effects components.
#
# cov_df Optional covariate reference values (for marginal predictions).
# include_kappa Include atlas-square random effects in predictions.
# include_aru Include ARU effect terms (if present in the model).
# A formula/expression suitable for inlabru prediction.
make_pred_formula_multiatlas <- function(cov_df = NULL, include_kappa = FALSE, include_aru = FALSE) {
  cov_terms <- character(0)
  
  if (!is.null(cov_df) && nrow(cov_df) > 0) {
    cov_terms <- cov_df %>%
      dplyr::mutate(term = paste0("Beta", beta, "_", covariate, "*I(", covariate, "^", beta, ")")) %>%
      dplyr::pull(term)
  }
  
  base_terms <- c(
    "Intercept",
    "spde_mean",
    "Atlas3_c * spde_diff",
    "Atlas3_c * effect_Atlas3",
    "DOY_global",
    "TOD_global"
  )
  
  if (include_aru) {
    base_terms <- c(base_terms, "ARU * effect_ARU")
  }
  
  eta_expr <- paste(c(base_terms, cov_terms), collapse = " + ")
  
  if (include_kappa) {
    eta_expr <- paste0(eta_expr, " + kappa_diff")
  }
  
  as.formula(paste0("~ data.frame(eta = ", eta_expr, ")"))
}

make_pred_grid <- function(grid_obba2, grid_obba3) {
  bind_rows(
    grid_obba2 %>% mutate(Atlas3 = 0L),
    grid_obba3 %>% mutate(Atlas3 = 1L)
  ) %>%
    mutate(
      days_midpoint = 0,
      Atlas3_c = Atlas3 - 0.5,
      BCR_idx = as.integer(factor(BCR))
    )
}

# Generate posterior draws for predictions from a fitted inlabru model
#
# Wrapper around inlabru prediction tools that:
# - evaluates an expression/formula at new locations,
# - returns posterior draws for the linear predictor and/or transformed mean,
# - and packages results in a consistent format used by downstream scripts.
#
# mod Fitted `bru` model.
# newdata Data frame/sf of prediction locations and covariates.
# formula Prediction formula/expression.
# n_samples Number of posterior draws to generate.
# A list containing draws and metadata suitable for summarization and mapping.
predict_inla <- function(mod,
                         grid,
                         pred_formula,
                         n.samples = 1000,
                         seed = 123) {
  stopifnot(!is.null(mod))
  stopifnot(is.data.frame(grid) || inherits(grid, "sf"))
  stopifnot(inherits(pred_formula, "formula"))
  
  pred <- inlabru::generate(
    mod,
    grid,
    formula = pred_formula,
    n.samples = n.samples,
    seed = seed
  )
  
  # pred is a list (length = n.samples) of 1-row data.frames
  vars <- names(pred[[1]])
  out <- lapply(vars, function(v) sapply(pred, function(x) x[[v]]))
  names(out) <- vars
  out
}

predict_all_pixels <- function(mod, pred_grid, pred_formula,
                               n.samples, seed,
                               on_water_col = "on_water") {
  
  preds <- predict_inla(
    mod = mod,
    grid = pred_grid,
    pred_formula = pred_formula,
    n.samples = n.samples,
    seed = seed
  )
  
  idx2 <- which(pred_grid$Atlas == "OBBA2")
  idx3 <- which(pred_grid$Atlas == "OBBA3")
  
  eta <- preds$eta
  mu2 <- exp(eta[idx2, , drop = FALSE])
  mu3 <- exp(eta[idx3, , drop = FALSE])
  
  if (nrow(mu2) != nrow(mu3)) {
    stop("Prediction grid does not contain matched OBBA2/OBBA3 rows.")
  }
  
  list(
    eta = eta,
    mu2 = mu2,
    mu3 = mu3
  )
}


# ============================================================
# 5. Lognormal variance correction
# ============================================================
#
# Why this exists
# ---------------
# make_pred_formula_multiatlas() deliberately omits the iid terms (kappa_diff,
# and site_iid where it was fitted), so exp(eta) is the conditional mean given
# kappa = 0 and site = 0. Both terms are fitted with constr = TRUE, and that
# sum-to-zero constraint acts on the LOG scale, so exp(eta) is the GEOMETRIC
# mean of the level-specific means: an "average square at an average site".
#
# The ARITHMETIC mean over the population of squares/sites is
#
#     E[exp(eta + z)] = exp(eta) * exp(sigma^2 / 2),      z ~ N(0, sigma^2)
#
# so the surfaces are rescaled by exp(0.5 * sum(sigma^2)) over the omitted iid
# terms. The mapped value then answers "what would a survey at a randomly
# chosen site return on average?" rather than "what would a survey at a typical
# (median) site return?". For species with large square- or site-level variance
# these differ substantially, and only the former aggregates coherently to
# population totals or compares sensibly against an absolute absence floor.
#
# What is deliberately NOT corrected
# ----------------------------------
#   - negative-binomial overdispersion: E[y] = mu already, so observation-level
#     overdispersion does not shift the mean and needs no correction.
#   - spde_mean, spde_diff, TOD_global, and the fixed effects: these ARE in the
#     prediction formula, so their posterior draws are already inside exp(eta).
#
# Known limitation
# ----------------
# site_iid is only active where site_w = 1, i.e. at sites with enough repeat
# visits to support a level -- in practice ARU stations. Its variance is
# therefore estimated from a non-random subdesign and then applied to every
# pixel, including regions surveyed only by single-visit point counts (whose
# site-level heterogeneity is absorbed by the nbinomial overdispersion instead,
# where it does not shift the mean). Whether that over- or under-corrects is an
# empirical question; set terms = "kappa_diff" to apply the square-level
# correction only.

# Safely extract a variance component (1 / precision) from a summary.hyperpar
# table.
#
# INLA names iid hyperparameters "Precision for <term>". A term that was
# dropped for a given species -- site_iid is dropped when n_site is below
# min_site_levels -- simply has no row in the table. Returning 0 for an absent
# term (rather than raising a subscript error) lets callers use one code path
# for species with and without a site random effect.
#
# hyper     summary.hyperpar table: mod$summary.hyperpar, or the
#           summary_hyperpar element saved by script 07.
# name      Component name as written in the model components, e.g. "site_iid".
# quantile  Column to read. "0.5quant" gives the posterior median precision;
#           because medians are invariant to monotone transformation,
#           1 / median(precision) is exactly the median variance.
# Numeric scalar >= 0.
get_var <- function(hyper, name, quantile = "0.5quant") {
  
  if (is.null(hyper) || is.null(rownames(hyper))) return(0)
  
  i <- match(paste0("Precision for ", name), rownames(hyper))
  if (is.na(i)) return(0)
  
  if (!quantile %in% colnames(hyper)) {
    stop("Column '", quantile, "' not found in summary.hyperpar; available: ",
         paste(colnames(hyper), collapse = ", "))
  }
  
  v <- 1 / hyper[i, quantile]
  if (!is.finite(v) || v < 0) 0 else v
}

# Is a hyperparameter present in a summary.hyperpar table?
#
# Distinguishes "term absent from this species' model" from "term present but
# with a variance estimated near zero", which get_var() alone cannot.
has_hyper <- function(hyper, name) {
  if (is.null(hyper) || is.null(rownames(hyper))) return(FALSE)
  !is.na(match(paste0("Precision for ", name), rownames(hyper)))
}

# Posterior draws of one or more variance components.
#
# Draws come from INLA::inla.hyperpar.sample(), which samples the joint
# posterior of the hyperparameters. Those draws are NOT paired with the
# latent-field draws produced by inlabru::generate(); treating the two as
# independent is a mild approximation that affects the joint dependence
# structure but not the marginal calibration of either piece. Pairing them
# properly would require driving both from a single inla.posterior.sample()
# call, which is more fragile across INLA versions than it is worth here.
#
# mod        Fitted bru/inla model.
# var_names  Character vector of component names.
# n_draws    Number of draws.
# Numeric matrix, n_draws x length(var_names), columns named after var_names.
# Components that cannot be located return a column of zeros.
sample_hyper_vars <- function(mod, var_names, n_draws) {
  
  out <- matrix(0, nrow = n_draws, ncol = length(var_names),
                dimnames = list(NULL, var_names))
  
  if (length(var_names) == 0) return(out)
  
  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop("INLA is required to sample hyperparameters.")
  }
  
  hyp <- INLA::inla.hyperpar.sample(n = n_draws, result = mod)
  
  if (is.null(dim(hyp)) || nrow(hyp) != n_draws) {
    stop("inla.hyperpar.sample() returned an unexpected shape.")
  }
  
  for (nm in var_names) {
    target <- paste0("Precision for ", nm)
    
    j <- match(target, colnames(hyp))
    if (is.na(j)) {
      # Some INLA versions decorate hyperparameter names; fall back to a fixed
      # (non-regex) partial match before treating the term as absent.
      j_alt <- which(grepl(target, colnames(hyp), fixed = TRUE))
      j <- if (length(j_alt) > 0) j_alt[1] else NA_integer_
    }
    if (is.na(j)) next
    
    v <- 1 / hyp[, j]
    v[!is.finite(v) | v < 0] <- 0
    out[, nm] <- v
  }
  
  out
}

# Build the multiplicative correction that rescales exp(eta) from the
# geometric mean over the omitted iid levels to the arithmetic mean.
#
# mod           Fitted bru/inla model.
# terms         Component names whose variance is being reinstated. Terms not
#               present in this species' model are skipped silently.
# propagate     If TRUE, return one correction factor per posterior draw so
#               that uncertainty in the variance components propagates into the
#               mapped intervals. If FALSE, return the single plug-in factor
#               exp(0.5 * sum(median variances)).
# n_draws       Number of posterior draws; required when propagate = TRUE, and
#               must equal the number of draws used for the predictions.
# cap_quantile  exp(sigma^2 / 2) is convex in sigma^2 and sigma^2 has a long
#               right tail for sparse species, so a handful of draws can
#               dominate the upper credible limit. Sampled variances are capped
#               at this quantile of their own posterior. Set to NA to disable.
#
# List with:
#   factor          numeric, length 1 (plug-in) or n_draws (propagated)
#   terms, present  which components were requested / actually found
#   var_median      median variance per requested term
#   factor_plugin   the plug-in factor, always returned for reference
#   factor_summary  median and 5th/95th percentiles of `factor`
#   propagated      whether draw-level propagation was actually used
#   cap_quantile    the cap applied, or NA
make_lognormal_correction <- function(mod,
                                      terms        = c("kappa_diff", "site_iid"),
                                      propagate    = TRUE,
                                      n_draws      = NULL,
                                      cap_quantile = 0.99) {
  
  stopifnot(!is.null(mod), is.character(terms))
  
  hyper      <- mod$summary.hyperpar
  present    <- vapply(terms, function(tm) has_hyper(hyper, tm), logical(1))
  var_median <- vapply(terms, function(tm) get_var(hyper, tm), numeric(1))
  
  names(present)    <- terms
  names(var_median) <- terms
  
  factor_plugin <- exp(0.5 * sum(var_median))
  
  make_result <- function(f_vec, propagated, cap) {
    list(
      factor         = f_vec,
      terms          = terms,
      present        = present,
      var_median     = var_median,
      factor_plugin  = factor_plugin,
      factor_summary = stats::quantile(f_vec, c(0.05, 0.5, 0.95), names = TRUE),
      propagated     = propagated,
      cap_quantile   = cap
    )
  }
  
  # Nothing to correct, or propagation not requested.
  if (length(terms) == 0 || !any(present) || !isTRUE(propagate)) {
    return(make_result(factor_plugin, propagated = FALSE, cap = NA_real_))
  }
  
  if (is.null(n_draws) || !is.finite(n_draws) || n_draws < 1) {
    stop("n_draws must be a positive integer when propagate = TRUE.")
  }
  
  var_draws <- try(
    sample_hyper_vars(mod, terms[present], n_draws),
    silent = TRUE
  )
  
  if (inherits(var_draws, "try-error")) {
    warning("Could not sample hyperparameters (",
            conditionMessage(attr(var_draws, "condition")),
            "); falling back to the plug-in variance correction.",
            call. = FALSE)
    return(make_result(factor_plugin, propagated = FALSE, cap = NA_real_))
  }
  
  cap_used <- NA_real_
  if (!is.null(cap_quantile) && length(cap_quantile) == 1 &&
      is.finite(cap_quantile) && cap_quantile > 0 && cap_quantile < 1) {
    for (j in seq_len(ncol(var_draws))) {
      cap <- stats::quantile(var_draws[, j], cap_quantile, na.rm = TRUE)
      var_draws[, j] <- pmin(var_draws[, j], cap)
    }
    cap_used <- cap_quantile
  }
  
  f <- exp(0.5 * rowSums(var_draws))
  
  if (!all(is.finite(f)) || any(f <= 0)) {
    warning("Non-finite variance correction draws; falling back to the ",
            "plug-in variance correction.", call. = FALSE)
    return(make_result(factor_plugin, propagated = FALSE, cap = NA_real_))
  }
  
  make_result(f, propagated = TRUE, cap = cap_used)
}

# Apply a variance correction to a matrix of posterior draws.
#
# mat         Numeric matrix, rows = prediction locations, cols = draws.
# correction  Length 1 (constant) or ncol(mat) (one factor per draw).
# log_scale   If TRUE, `mat` holds linear-predictor draws and log(correction)
#             is ADDED rather than the factor being multiplied.
apply_lognormal_correction <- function(mat, correction, log_scale = FALSE) {
  
  stopifnot(is.matrix(mat), is.numeric(correction))
  
  if (length(correction) != 1 && length(correction) != ncol(mat)) {
    stop("correction must have length 1 or ncol(mat) (= number of posterior ",
         "draws); got ", length(correction), " for ", ncol(mat), " draws.")
  }
  
  if (log_scale) {
    adj <- log(correction)
    if (length(adj) == 1) return(mat + adj)
    return(sweep(mat, 2, adj, "+"))
  }
  
  if (length(correction) == 1) return(mat * correction)
  sweep(mat, 2, correction, "*")
}

# One-line description of a correction object, for log messages and for
# reporting the scale of a saved prediction file.
describe_lognormal_correction <- function(vc, digits = 3) {
  
  if (is.null(vc)) return("no variance correction recorded")
  
  used <- vc$terms[vc$present]
  
  if (length(used) == 0) {
    return("variance correction = 1.000 (no iid terms in this model)")
  }
  
  vars_txt <- paste0(
    used, " sigma^2 = ", round(vc$var_median[used], digits),
    collapse = "; "
  )
  
  fs <- vc$factor_summary
  
  paste0(
    "variance correction x", round(fs[["50%"]], digits),
    " [", round(fs[["5%"]], digits), ", ", round(fs[["95%"]], digits), "]",
    if (isTRUE(vc$propagated)) " (propagated" else " (plug-in",
    if (isTRUE(vc$propagated) && is.finite(vc$cap_quantile)) {
      paste0(", capped at q", vc$cap_quantile, ")")
    } else ")",
    "; ", vars_txt
  )
}

# ============================================================
# 6. Posterior summaries
# ============================================================

# Summarize posterior draws row-wise (median + credible intervals)
#
# Given a matrix of posterior draws (rows = prediction locations, columns = draws),
# compute common summaries used in mapping and reporting: mean, median, sd,
# and quantiles (e.g., 2.5%, 50%, 97.5%).
#
# mat Numeric matrix of posterior draws.
# probs Quantiles to compute.
# A data.frame of summaries with one row per prediction location.
summarize_posterior <- function(mat,
                                CI_probs = c(0.05, 0.95),
                                prefix = "var") {
  stopifnot(is.matrix(mat))
  stopifnot(length(CI_probs) == 2, CI_probs[1] < CI_probs[2])
  
  # avoid dependency surprises: use matrixStats if installed, otherwise base
  if (requireNamespace("matrixStats", quietly = TRUE)) {
    mean_vals   <- matrixStats::rowMeans2(mat, na.rm = TRUE)
    median_vals <- matrixStats::rowMedians(mat, na.rm = TRUE)
    sd_vals     <- matrixStats::rowSds(mat, na.rm = TRUE)
    lower_vals  <- matrixStats::rowQuantiles(mat, probs = CI_probs[1], na.rm = TRUE)
    upper_vals  <- matrixStats::rowQuantiles(mat, probs = CI_probs[2], na.rm = TRUE)
  } else {
    mean_vals   <- rowMeans(mat, na.rm = TRUE)
    median_vals <- apply(mat, 1, median, na.rm = TRUE)
    sd_vals     <- apply(mat, 1, sd, na.rm = TRUE)
    lower_vals  <- apply(mat, 1, quantile, probs = CI_probs[1], na.rm = TRUE)
    upper_vals  <- apply(mat, 1, quantile, probs = CI_probs[2], na.rm = TRUE)
  }
  
  cv_vals <- sd_vals / median_vals
  
  data.frame(
    setNames(list(mean_vals),   paste0(prefix, "_mean")),
    setNames(list(median_vals), paste0(prefix, "_q50")),
    setNames(list(sd_vals),     paste0(prefix, "_sd")),
    setNames(list(cv_vals),     paste0(prefix, "_cv_median")),
    setNames(list(lower_vals),  paste0(prefix, "_lower")),
    setNames(list(upper_vals),  paste0(prefix, "_upper"))
  )
}

summarize_predictions <- function(mu2, mu3) {
  abs_change <- mu3 - mu2
  list(
    OBBA2 = summarize_posterior(mu2, CI_probs = c(0.05, 0.95), prefix = "OBBA2"),
    OBBA3 = summarize_posterior(mu3, CI_probs = c(0.05, 0.95), prefix = "OBBA3"),
    abs_change = summarize_posterior(abs_change, CI_probs = c(0.05, 0.95), prefix = "abs_change")
  )
}

# 
make_hex_draws <- function(g2, mu2, mu3) {
  
  if (!("hex_id" %in% names(g2))) {
    stop("Cannot create hex_draws because `hex_id` is not present in prediction grid.")
  }
  
  hex_ids <- g2$hex_id
  u_hex <- unique(hex_ids)
  
  bind_rows(
    lapply(u_hex, function(hx) {
      
      idx_hex <- which(hex_ids == hx)
      
      mu2_hex <- colMeans(mu2[idx_hex, , drop = FALSE])
      mu3_hex <- colMeans(mu3[idx_hex, , drop = FALSE])
      
      tibble(
        hex_id = hx,
        n_pixels = length(idx_hex),
        mu_OBBA2 = list(mu2_hex),
        mu_OBBA3 = list(mu3_hex),
        abs_change = list(mu3_hex - mu2_hex)
      )
    })
  )
}


# ============================================================
# 7. Shared-footprint (paired) analysis
# ============================================================

# Function to only pull surveys from a shared footprint
make_shared_footprint_dataset <- function(
    dat,
    atlas_col = "Atlas",
    atlas_levels = c("OBBA2", "OBBA3"),
    buffer_km = 0.5,
    row_id_col = "shared_footprint_row_id",
    flag_col = "shared_footprint"
) {
  
  # ---- Input checks ----
  stopifnot(inherits(dat, "sf"))
  stopifnot(atlas_col %in% names(dat))
  stopifnot(length(atlas_levels) == 2)
  
  # ---- Add stable row ID for rejoining counts or other external objects ----
  dat <- dat |>
    dplyr::mutate(
      "{row_id_col}" := dplyr::row_number()
    )
  
  # ---- Split surveys by atlas period ----
  a1 <- dat |>
    dplyr::filter(.data[[atlas_col]] == atlas_levels[1])
  
  a2 <- dat |>
    dplyr::filter(.data[[atlas_col]] == atlas_levels[2])
  
  if (nrow(a1) == 0 || nrow(a2) == 0) {
    out <- dat[0, ]
    out[[flag_col]] <- logical(0)
    return(out)
  }
  
  # ---- Find cross-period neighbours within buffer distance ----
  # CRS units are assumed to be kilometres, so buffer_km is passed directly.
  a1_to_a2 <- sf::st_is_within_distance(
    a1,
    a2,
    dist = buffer_km
  )
  
  # ---- Symmetric shared-footprint logic ----
  # Keep:
  #   1. atlas-1 surveys that are near at least one atlas-2 survey
  #   2. atlas-2 surveys that are near at least one retained atlas-1 survey
  #
  # This avoids the previous second-pass asymmetry.
  
  a1_keep <- lengths(a1_to_a2) > 0
  
  a2_keep_idx <- sort(unique(unlist(a1_to_a2)))
  
  a1_matched <- a1[a1_keep, ]
  a2_matched <- a2[a2_keep_idx, ]
  
  # ---- Combine retained surveys ----
  shared_dat <- dplyr::bind_rows(a1_matched, a2_matched) |>
    dplyr::mutate(
      "{flag_col}" := TRUE
    ) |>
    dplyr::arrange(.data[[row_id_col]])
  
  shared_dat
}

# Fits only point counts and ARUs
fit_inla_shared_footprint_change <- function(
    sp_dat_shared,
    region_col = "BCR",
    family = c("poisson", "nbinomial"),
    
    mesh_TOD = NULL,
    fixed_TOD_range = 1.5,                # hours; NOT estimated
    prior_TOD_sigma = c(3, 0.1),

    mesh_DOY = NULL,
    fixed_DOY_range = 30,                 # days; NOT estimated
    prior_DOY_sigma = c(2, 0.1),          # P(sigma > 2) = 0.1 -> prior median ~0.6
    
    square_pcprec = c(log(2), 0.1),
    
    timeout_min = 5,
    int_strategy = "eb",
    strategy = "simplified.laplace",
    inla_mode = "experimental",
    bru_verbose = 4,
    waic = TRUE,
    cpo = FALSE
) {
  
  family <- match.arg(family)
  
  # ------------------------------------------------------------
  # 1. Validate inputs
  # ------------------------------------------------------------
  
  if (!inherits(sp_dat_shared, "sf")) {
    stop("sp_dat_shared must be an sf object.")
  }
  
  required_cols <- c(
    "count",
    "Atlas",
    "Hours_After_Reference",
    "days_midpoint",
    "Survey_Type",
    "square_atlas",
    region_col
  )
  
  missing_cols <- setdiff(required_cols, names(sp_dat_shared))
  
  if (length(missing_cols) > 0) {
    stop("Missing columns: ", paste(missing_cols, collapse = ", "))
  }
  
  INLA::inla.setOption(inla.timeout = 60 * timeout_min)
  
  # ------------------------------------------------------------
  # 2. Prepare data
  # ------------------------------------------------------------
  
  sp_dat_shared <- sp_dat_shared |>
    dplyr::filter(Survey_Type %in% c("Point_Count", "ARU")) |>
    dplyr::mutate(
      count = as.numeric(count),
      Atlas3 = as.numeric(Atlas == "OBBA3"),
      ARU = as.numeric(Survey_Type == "ARU"),
      Hours_After_Reference = as.numeric(Hours_After_Reference),
      days_midpoint = as.numeric(days_midpoint),
      region_factor = factor(.data[[region_col]])
    )
  
  if (nrow(sp_dat_shared) == 0) {
    stop("No point-count or ARU surveys remain after filtering.")
  }
  
  # ------------------------------------------------------------
  # 3. Numeric iid index for atlas-square random effect
  # ------------------------------------------------------------
  
  square_levels <- sort(unique(sp_dat_shared$square_atlas))
  
  square_lookup <- tibble::tibble(
    square_atlas = square_levels,
    square_id_iid = seq_along(square_levels)
  )
  
  sp_dat_shared <- sp_dat_shared |>
    dplyr::left_join(square_lookup, by = "square_atlas")
  
  n_square <- nrow(square_lookup)
  
  # ------------------------------------------------------------
  # 4. Cell-means fixed effects for regions and change
  # ------------------------------------------------------------
  # region_* = region-specific OBBA2 log-scale intercept
  # change_* = region-specific OBBA3 - OBBA2 log-scale change
  
  # Built by hand (rather than via stats::model.matrix(~0 + region_factor))
  # because model.matrix() always routes factor terms through R's
  # `contrasts<-` machinery -- even for a no-intercept, cell-means formula --
  # and that machinery errors ("contrasts can be applied only to factors
  # with 2 or more levels") whenever the shared footprint happens to contain
  # only a single region. This explicit 0/1 encoding has no such
  # restriction and is numerically identical to model.matrix()'s output
  # whenever there are 2+ levels.
  region_levels <- levels(sp_dat_shared$region_factor)
  region_chr    <- as.character(sp_dat_shared$region_factor)
  
  region_mm <- matrix(
    0,
    nrow = nrow(sp_dat_shared),
    ncol = length(region_levels),
    dimnames = list(NULL, paste0("region_", region_levels))
  )
  
  for (j in seq_along(region_levels)) {
    region_mm[, j] <- as.numeric(region_chr == region_levels[j])
  }
  
  change_mm <- region_mm * sp_dat_shared$Atlas3
  colnames(change_mm) <- paste0("change_", region_levels)
  
  sp_dat_shared <- dplyr::bind_cols(
    sp_dat_shared,
    as.data.frame(region_mm),
    as.data.frame(change_mm)
  )
  
  region_terms <- paste(colnames(region_mm), collapse = " + ")
  change_terms <- paste(colnames(change_mm), collapse = " + ")
  
  # ------------------------------------------------------------
  # 5. Build 1-D SPDE smoothers
  # ------------------------------------------------------------
  
  # --- Time-of-day smooth: SHARED mesh, FIXED range ---------------------------
  # The mesh is built once (make_TOD_mesh) and passed in, so every species uses
  # the same basis and the same absolute range means the same thing everywhere.
  # Deriving it from range(sp_dat$Hours_After_Reference) made the resolution
  # species-dependent.
  #
  # prior.range = c(range0, NA) -> range HELD FIXED at range0 (no prior).
  # Only sigma is estimated. Same idiom as the effort SPDEs.
  if (is.null(mesh_TOD)) {
    mesh_TOD <- make_TOD_mesh(knot_spacing = fixed_TOD_range / 6,
                              pad_range    = fixed_TOD_range)
  }
  
  # Guard the requirement that makes the fixed range meaningful: a basis coarser
  # than range/5 cannot represent the correlation length being asked for, and
  # the mesh, not the range, then sets the smoothness.
  TOD_knot_gap <- max(diff(sort(mesh_TOD$loc)))
  if (TOD_knot_gap > fixed_TOD_range / 5) {
    warning(
      "TOD mesh knot spacing (", signif(TOD_knot_gap, 3),
      " h) is coarser than fixed_TOD_range/5 (",
      signif(fixed_TOD_range / 5, 3),
      " h). The mesh, not the range, is setting the smoothness.",
      call. = FALSE
    )
  }
  
  TOD_spde <- INLA::inla.spde2.pcmatern(
    mesh        = mesh_TOD,
    alpha       = 2,
    prior.range = c(fixed_TOD_range, NA),
    prior.sigma = prior_TOD_sigma,
    constr      = TRUE
  )
  
  # --- Day-of-year smooth: FIXED range, sigma estimated ----------------------
  # Added so the paired analysis carries the SAME seasonal detectability
  # adjustment as the full model in script 07 (fit_PC_ARU). Modelled exactly like
  # the TOD smooth above: a 1-D SPDE with the RANGE HELD FIXED at fixed_DOY_range
  # and only sigma estimated. The domain is NOT cyclic -- days_midpoint is a
  # linear axis centred on the species' safe-date midpoint (0), so a free mesh is
  # used. DOY_global is constr = TRUE (sum-to-zero over the data), so it absorbs
  # systematic OBBA2-vs-OBBA3 timing differences within the shared footprint
  # without competing with the change_* fixed effects.
  #
  # When mesh_DOY is not supplied it is built to span this species' observed
  # days_midpoint range (+/- 7 d), matching how script 07 sizes its per-species
  # DOY mesh. Pass mesh_DOY explicitly to share one basis across species.
  if (is.null(mesh_DOY)) {
    DOY_obs_range <- range(sp_dat_shared$days_midpoint, na.rm = TRUE)
    if (!all(is.finite(DOY_obs_range))) DOY_obs_range <- c(-60, 60)
    mesh_DOY <- make_DOY_mesh(
      span         = c(DOY_obs_range[1] - 7, DOY_obs_range[2] + 7),
      knot_spacing = fixed_DOY_range / 6,
      pad_range    = fixed_DOY_range
    )
  }

  # Same knot-spacing guard as for TOD: a basis coarser than range/5 cannot
  # represent the correlation length being asked for.
  DOY_knot_gap <- max(diff(sort(mesh_DOY$loc)))
  if (DOY_knot_gap > fixed_DOY_range / 5) {
    warning(
      "DOY mesh knot spacing (", signif(DOY_knot_gap, 3),
      " d) is coarser than fixed_DOY_range/5 (",
      signif(fixed_DOY_range / 5, 3),
      " d). The mesh, not the range, is setting the smoothness.",
      call. = FALSE
    )
  }

  DOY_spde <- INLA::inla.spde2.pcmatern(
    mesh        = mesh_DOY,
    alpha       = 2,
    prior.range = c(fixed_DOY_range, NA),
    prior.sigma = prior_DOY_sigma,
    constr      = TRUE
  )

  # ------------------------------------------------------------
  # 6. Priors
  # ------------------------------------------------------------
  
  pc_prec_square <- list(
    prior = "pcprec",
    param = square_pcprec
  )
  
  # ------------------------------------------------------------
  # 7. Fixed-effect components
  # ------------------------------------------------------------
  
  region_coef_names <- paste0("b_region_", region_levels)
  change_coef_names <- paste0("b_change_", region_levels)
  
  region_components <- paste0(
    region_coef_names,
    "(1, model = 'linear', mean.linear = 0.2, prec.linear = 1 / (3^2))"
  )
  
  change_components <- paste0(
    change_coef_names,
    "(1, model = 'linear', mean.linear = 0, prec.linear = 1 / (log(2)^2))"
  )
  
  fixed_components <- paste(
    c(region_components, change_components),
    collapse = " + "
  )
  
  region_terms <- paste(
    paste0(colnames(region_mm), " * ", region_coef_names),
    collapse = " + "
  )
  
  change_terms <- paste(
    paste0(colnames(change_mm), " * ", change_coef_names),
    collapse = " + "
  )
  
  # ------------------------------------------------------------
  # 8. Define inlabru components
  # ------------------------------------------------------------
  # NOTE: no global Intercept(1), because region_* terms are cell means.
  
  components <- stats::as.formula(paste0(
    "~ ",
    fixed_components,
    
    " + effect_ARU(
        1,
        model = 'linear',
        mean.linear = 0,
        prec.linear = 1 / ((log(1.5) / 2)^2)
      )",
    
    " + square_effect(
        main = square_id_iid,
        model = 'iid',
        n = n_square,
        constr = TRUE,
        hyper = list(prec = pc_prec_square)
      )",
    
    " + TOD_global(
        main = Hours_After_Reference,
        model = TOD_spde
      )",
    
    " + DOY_global(
        main = days_midpoint,
        model = DOY_spde
      )"
  ))
  
  # ------------------------------------------------------------
  # 9. Define likelihood formula
  # ------------------------------------------------------------
  # NOTE: no global intercept term here either.
  
  formula <- stats::as.formula(paste0(
    "count ~ ",
    region_terms,
    " + ",
    change_terms,
    " + ARU * effect_ARU",
    " + TOD_global",
    " + DOY_global",
    " + square_effect"
  ))
  
  # ------------------------------------------------------------
  # 10. Fit model
  # ------------------------------------------------------------
  
  lik <- inlabru::like(
    family = family,
    formula = formula,
    data = sp_dat_shared
  )
  
  fit <- inlabru::bru(
    components = components,
    lik,
    options = list(
      inla.mode = inla_mode,
      control.inla = list(
        int.strategy = int_strategy,
        strategy = strategy
      ),
      control.compute = list(
        waic = waic,
        cpo = cpo
      ),
      bru_verbose = bru_verbose
    )
  )
  
  fit$region_levels   <- region_levels
  fit$square_lookup   <- square_lookup
  fit$tod_fixed_range <- fixed_TOD_range
  fit$tod_knot_gap    <- TOD_knot_gap
  fit$doy_fixed_range <- fixed_DOY_range
  fit$doy_knot_gap    <- DOY_knot_gap
  
  fit
}

summarize_shared_footprint_change <- function(mod_shared, region_levels = NULL) {
  
  # Region ordering. Levels are taken from the fitted model, which stores them
  # in fit_inla_shared_footprint_change(), rather than being hardcoded here.
  # The previous hardcoded set c("12","13","74","76_77") silently converted any
  # unlisted region code to NA -- so splitting BCR 76_77 into 76 and 77 would
  # have dropped the region labels from every paired change summary without
  # error. Falling back to the codes actually present keeps this correct for
  # any future BCR definition.
  if (is.null(region_levels)) {
    region_levels <- mod_shared$region_levels
  }
  
  probs <- c(
    q025 = 0.025,
    q05  = 0.05,
    q10  = 0.10,
    q50  = 0.50,
    q90  = 0.90,
    q95  = 0.95,
    q975 = 0.975
  )
  
  change_names <- grep(
    "^b_change_",
    rownames(mod_shared$summary.fixed),
    value = TRUE
  )
  
  change_summary <- purrr::map_dfr(change_names, function(term) {
    
    marg <- mod_shared$marginals.fixed[[term]]
    
    q_log <- INLA::inla.qmarginal(
      p = unname(probs),
      marginal = marg
    )
    
    names(q_log) <- names(probs)
    
    q_pct <- 100 * (exp(q_log) - 1)
    
    log_change_mean <- mod_shared$summary.fixed[term, "mean"]
    
    tibble::tibble(
      BCR  = sub("^b_change_", "", term),
      term = term,
      
      log_change_mean = log_change_mean,
      log_change_sd   = mod_shared$summary.fixed[term, "sd"],
      
      pct_change_mean = 100 * (exp(log_change_mean) - 1),
      pct_change_q50  = unname(q_pct["q50"]),
      
      pct_change_q025 = unname(q_pct["q025"]),
      pct_change_q05  = unname(q_pct["q05"]),
      pct_change_q10  = unname(q_pct["q10"]),
      pct_change_q90  = unname(q_pct["q90"]),
      pct_change_q95  = unname(q_pct["q95"]),
      pct_change_q975 = unname(q_pct["q975"])
    )
  })
  
  # Any code present in the fit but absent from region_levels is appended rather
  # than turned into NA, so a stale level set can never silently discard a
  # region (this also covers region_levels being NULL on legacy fits).
  region_levels <- union(region_levels, sort(unique(change_summary$BCR)))
  
  change_summary |>
    dplyr::mutate(
      BCR = factor(BCR, levels = region_levels)
    ) |>
    dplyr::arrange(BCR) |>
    dplyr::mutate(
      BCR = as.character(BCR)
    )
}
