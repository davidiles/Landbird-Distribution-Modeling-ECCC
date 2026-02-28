# ============================================================
# R/functions/inla_model_utils.R
#
# Purpose:
#   Helper functions for INLA/inlabru species model workflow:
#     - fit_inla_testing(): fit a spatio-temporal abundance model
#     - predict_inla(): posterior draws for arbitrary expressions
#     - summarize_posterior(): row-wise posterior summaries
#     - calibrate_pc_lambda(): PC prior calibration helper (optional)
#
# Notes:
#   - Designed to be called from 08_fit_models_and_predict.R
#   - No reliance on global variables (species name, sp_code, etc.)
#   - Uses inlabru::bru() + INLA SPDEs
# ============================================================

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(stringr)
  
  library(INLA)
  library(inlabru)
  library(fmesher)
})

# ------------------------------------------------------------
# Small utilities
# ------------------------------------------------------------

`%||%` <- function(x, y) if (is.null(x)) y else x

ensure_numeric <- function(x) as.numeric(x)

# ------------------------------------------------------------
# PC prior lambda calibration for inla.pc.rgamma
# (Useful if you want a PC prior on NB "size"/theta)
# ------------------------------------------------------------

calibrate_pc_lambda <- function(target_prob,
                                threshold_theta,
                                n_mc = 20000,
                                lower = 1e-6,
                                upper = 20,
                                tol = 1e-3) {
  stopifnot(is.numeric(target_prob), length(target_prob) == 1, target_prob > 0, target_prob < 1)
  stopifnot(is.numeric(threshold_theta), length(threshold_theta) == 1, threshold_theta > 0)
  
  pc_tail_prob <- function(lambda, n = n_mc, threshold = threshold_theta) {
    samp <- INLA::inla.pc.rgamma(n = n, lambda = lambda)
    mean(samp > threshold)
  }
  
  f_for_root <- function(lambda) pc_tail_prob(lambda) - target_prob
  
  sol <- try(uniroot(f_for_root, lower = lower, upper = upper, tol = tol), silent = TRUE)
  if (inherits(sol, "try-error")) {
    stop("Could not find lambda in the specified interval. Try increasing `upper` or `n_mc`.")
  }
  sol$root
}

# ------------------------------------------------------------
# Posterior summaries (row-wise) for n_row x n_draw matrices
# ------------------------------------------------------------

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

# ------------------------------------------------------------
# Prediction draws via inlabru::generate()
# Returns a named list of matrices, one per output variable.
# Each matrix is nrow(grid) x n.samples
# ------------------------------------------------------------

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

# ------------------------------------------------------------
# Core fitter
# ------------------------------------------------------------

fit_inla_multi_atlas <- function(sp_dat,
                                 study_boundary,
                                 covariates,
                                 timeout_min = 15,
                                 # SPDE priors
                                 prior_range_abund  = c(150, 0.1),
                                 prior_sigma_abund  = c(0.5, 0.1),
                                 prior_range_change = c(500, 0.1),
                                 prior_sigma_change = c(0.1, 0.1),
                                 # Mesh controls
                                 mesh_max_edge = c(100, 200),
                                 mesh_cutoff = 50,
                                 mesh_convex = c(50, 150),
                                 mesh_concave = c(50, 150),
                                 
                                 # 1D smoothers
                                 HSS_prior_range  = c(2, 0.1),
                                 HSS_prior_sigma  = c(1, 0.5),
                                 
                                 DOY_prior_range  = c(30, 0.1),
                                 DOY_prior_sigma  = c(1, 0.5),
                                 
                                 # atlas square iid effect
                                 kappa_pcprec = c(1, 0.1),
                                 # likelihood
                                 family = c("poisson", "nbinomial"),
                                 # NB PC prior config (only used if family="nbinomial")
                                 nb_pc_target_prob = 0.5,
                                 nb_pc_threshold_theta = 5,
                                 # inlabru/INLA options
                                 inla_mode = "experimental",
                                 int_strategy = "eb",
                                 strategy = "simplified.laplace",
                                 bru_verbose = 4,
                                 waic = FALSE,
                                 cpo = FALSE,
                                 retry = 0) {
  
  family <- match.arg(family)
  
  stopifnot(is.data.frame(sp_dat) || inherits(sp_dat, "sf"))
  stopifnot(inherits(study_boundary, "sf"))
  
  required_cols <- c("count", "Atlas3", "Hours_Since_Sunrise", "days_since_june15", "square_atlas", "geometry")
  missing <- setdiff(required_cols, names(sp_dat))
  if (length(missing) > 0) stop("sp_dat is missing required columns: ", paste(missing, collapse = ", "))
  
  # covariates is a data.frame with at least: covariate, beta, model, mean, prec
  if (!is.data.frame(covariates)) stop("covariates must be a data.frame")
  needed_cov_cols <- c("covariate", "beta", "model", "mean", "prec")
  cov_missing <- setdiff(needed_cov_cols, names(covariates))
  if (length(cov_missing) > 0) stop("covariates is missing columns: ", paste(cov_missing, collapse = ", "))
  
  # Timeout (seconds)
  INLA::inla.setOption(inla.timeout = 60 * timeout_min)
  
  # Make sure geometry exists and CRS alignments are consistent
  if (!inherits(sp_dat, "sf")) {
    stop("sp_dat must be an sf object with a geometry column.")
  }
  if (sf::st_crs(sp_dat) != sf::st_crs(study_boundary)) {
    study_boundary <- sf::st_transform(study_boundary, sf::st_crs(sp_dat))
  }
  
  # ------------------------------------------------------------
  # Mesh + SPDEs
  # ------------------------------------------------------------
  
  hull <- fmesher::fm_extensions(
    study_boundary,
    convex  = mesh_convex,
    concave = mesh_concave
  )
  
  mesh_spatial <- fmesher::fm_mesh_2d_inla(
    loc = sf::st_as_sfc(sp_dat),
    boundary = hull,
    max.edge = mesh_max_edge,
    cutoff = mesh_cutoff,
    crs = sf::st_crs(sp_dat)
  )
  
  matern_abund <- INLA::inla.spde2.pcmatern(
    mesh_spatial,
    prior.range = prior_range_abund,
    prior.sigma = prior_sigma_abund,
    constr = TRUE
  )
  
  matern_change <- INLA::inla.spde2.pcmatern(
    mesh_spatial,
    prior.range = prior_range_change,
    prior.sigma = prior_sigma_change,
    constr = TRUE
  )
  
  # ------------------------------------------------------------
  # 1D effects: Hours since sunrise + Day-of-year
  # ------------------------------------------------------------
  
  sp_dat$Hours_Since_Sunrise <- ensure_numeric(sp_dat$Hours_Since_Sunrise)
  sp_dat$days_since_june15   <- ensure_numeric(sp_dat$days_since_june15)
  sp_dat$Atlas3              <- ensure_numeric(sp_dat$Atlas3)
  
  # Hours-since-sunrise mesh
  HSS_range <- range(sp_dat$Hours_Since_Sunrise, na.rm = TRUE)
  HSS_meshpoints <- seq(HSS_range[1] - 1, HSS_range[2] + 1, length.out = 41)
  HSS_mesh1D <- INLA::inla.mesh.1d(HSS_meshpoints, boundary = "free")
  HSS_spde <- INLA::inla.spde2.pcmatern(
    HSS_mesh1D,
    prior.range = HSS_prior_range,
    prior.sigma = HSS_prior_sigma,
    constr = TRUE
  )
  
  # DOY mesh (days since June 15)
  DOY_range <- range(sp_dat$days_since_june15, na.rm = TRUE)
  DOY_meshpoints <- seq(DOY_range[1] - 5, DOY_range[2] + 5, length.out = 41)
  DOY_mesh1D <- INLA::inla.mesh.1d(DOY_meshpoints, boundary = "free")
  DOY_spde_global <- INLA::inla.spde2.pcmatern(
    DOY_mesh1D,
    prior.range = DOY_prior_range,
    prior.sigma = DOY_prior_sigma,
    constr = TRUE
  )
  
  # ------------------------------------------------------------
  # IID effect for atlas squares
  # ------------------------------------------------------------
  
  pc_prec <- list(prior = "pcprec", param = kappa_pcprec)
  
  # ------------------------------------------------------------
  # Covariate components (linear by default)
  # covariates df must have: covariate, beta, model, mean, prec
  # Creates components: Beta<beta>_<cov>(1, model="linear", ...)
  # And formula terms: Beta... * I(cov^beta)
  # ------------------------------------------------------------
  
  covariates <- covariates %>%
    mutate(
      components = paste0(
        "Beta", beta, "_", covariate,
        '(1, model="', model,
        '", mean.linear=', mean,
        ", prec.linear=", prec,
        ")"
      ),
      formula = paste0("Beta", beta, "_", covariate, "*I(", covariate, "^", beta, ")")
    )
  
  covar_components_str <- if (nrow(covariates) > 0) paste(covariates$components, collapse = " + ") else ""
  covar_terms_str <- if (nrow(covariates) > 0) paste(covariates$formula, collapse = " + ") else ""
  
  # ------------------------------------------------------------
  # Model components (bru)
  # ------------------------------------------------------------
  
  components_str <- paste0(
    '~ Intercept(1) +
       effect_Atlas3(1, model="linear", mean.linear=0, prec.linear=1) +
       effect_ARU(1, model="linear", mean.linear=0, prec.linear=100) +
       spde_abund(main = geometry, model = matern_abund) +
       spde_change(main = geometry, model = matern_change) +
       HSS(main = Hours_Since_Sunrise, model = HSS_spde) +
       DOY_global(main = days_since_june15, model = DOY_spde_global) +
       kappa(square_atlas, model="iid", constr=TRUE, hyper=list(prec=pc_prec))',
    if (nchar(covar_components_str) > 0) paste0(" + ", covar_components_str) else ""
  )
  
  model_components <- as.formula(components_str)
  
  # ------------------------------------------------------------
  # Likelihood formula (shared by ARU / point counts)
  # ------------------------------------------------------------
  
  # You can keep ARU * effect_ARU term even if ARU is all zeros; it just becomes irrelevant.
  base_formula_str <- paste0(
    "count ~
      Intercept +
      ARU * effect_ARU +
      HSS +
      DOY_global +
      kappa +
      spde_abund +
      Atlas3 * spde_change +
      Atlas3 * effect_Atlas3",
    if (nchar(covar_terms_str) > 0) paste0(" + ", covar_terms_str) else ""
  )
  model_formula <- as.formula(base_formula_str)
  
  # ------------------------------------------------------------
  # Likelihood setup (two likelihoods: ARU==0 and ARU==1)
  # ------------------------------------------------------------
  
  make_like <- function(data_subset) {
    if (family == "poisson") {
      inlabru::like(
        family = "poisson",
        formula = model_formula,
        data = data_subset
      )
    } else {
      # negative binomial with PC prior on theta (size)
      lambda <- calibrate_pc_lambda(
        target_prob = nb_pc_target_prob,
        threshold_theta = nb_pc_threshold_theta
      )
      
      inlabru::like(
        family = "nbinomial",
        formula = model_formula,
        data = data_subset,
        control.family = list(
          hyper = list(theta = list(
            prior = "pc.gamma",
            param = c(lambda)
          ))
        )
      )
    }
  }
  
  dat_pc  <- sp_dat %>% dplyr::filter(ARU == 0)
  dat_aru <- sp_dat %>% dplyr::filter(ARU == 1)
  
  # If one of the subsets is empty, bru() can still run with a single likelihood
  likes <- list()
  if (nrow(dat_pc)  > 0) likes <- c(likes, list(make_like(dat_pc)))
  if (nrow(dat_aru) > 0) likes <- c(likes, list(make_like(dat_aru)))
  
  if (length(likes) == 0) stop("No data rows for ARU==0 or ARU==1 after filtering; cannot fit model.")
  
  # ------------------------------------------------------------
  # Fit with optional retries
  # ------------------------------------------------------------
  
  bru_opts <- list(
    inla.mode = inla_mode,
    control.compute = list(waic = waic, cpo = cpo),
    control.inla = list(
      int.strategy = int_strategy,
      strategy = strategy
    ),
    bru_verbose = bru_verbose
  )
  
  attempt <- 0
  last_err <- NULL
  fit <- NULL
  
  while (attempt <= retry) {
    attempt <- attempt + 1
    
    fit <- try(
      do.call(
        inlabru::bru,
        c(
          list(components = model_components),
          likes,
          list(options = bru_opts)
        )
      ),
      silent = TRUE
    )
    
    if (!inherits(fit, "try-error")) break
    last_err <- fit
    fit <- NULL
  }
  
  if (is.null(fit)) {
    stop("bru() failed after ", attempt, " attempt(s).\nLast error:\n", as.character(last_err))
  }
  
  fit
}

fit_inla_single_atlas <- function(sp_dat,
                                  study_boundary,
                                  covariates,
                                  timeout_min = 15,
                                  
                                  # SPDE priors
                                  prior_range_abund  = c(150, 0.1),
                                  prior_sigma_abund  = c(0.5, 0.1),
                                  
                                  # Mesh controls
                                  mesh_max_edge = c(100, 200),
                                  mesh_cutoff = 50,
                                  mesh_convex = c(50, 150),
                                  mesh_concave = c(50, 150),
                                  
                                  # 1D smoothers
                                  HSS_prior_range  = c(2, 0.1),
                                  HSS_prior_sigma  = c(1, 0.5),
                                  
                                  DOY_prior_range  = c(30, 0.1),
                                  DOY_prior_sigma  = c(1, 0.5),
                                  
                                  # atlas square iid effect
                                  kappa_pcprec = c(1, 0.1),
                                  
                                  # likelihood
                                  family = c("poisson", "nbinomial"),
                                  
                                  # NB PC prior config (only used if family="nbinomial")
                                  nb_pc_target_prob = 0.5,
                                  nb_pc_threshold_theta = 5,
                                  # inlabru/INLA options
                                  inla_mode = "experimental",
                                  int_strategy = "eb",
                                  strategy = "simplified.laplace",
                                  bru_verbose = 4,
                                  waic = FALSE,
                                  cpo = FALSE,
                                  retry = 0) {
  
  family <- match.arg(family)
  
  stopifnot(is.data.frame(sp_dat) || inherits(sp_dat, "sf"))
  stopifnot(inherits(study_boundary, "sf"))
  
  required_cols <- c("count", "Hours_Since_Sunrise", "days_since_june15", "square_atlas", "geometry")
  missing <- setdiff(required_cols, names(sp_dat))
  if (length(missing) > 0) stop("sp_dat is missing required columns: ", paste(missing, collapse = ", "))
  
  # covariates is a data.frame with at least: covariate, beta, model, mean, prec
  if (!is.data.frame(covariates)) stop("covariates must be a data.frame")
  needed_cov_cols <- c("covariate", "beta", "model", "mean", "prec")
  cov_missing <- setdiff(needed_cov_cols, names(covariates))
  if (length(cov_missing) > 0) stop("covariates is missing columns: ", paste(cov_missing, collapse = ", "))
  
  # Timeout (seconds)
  INLA::inla.setOption(inla.timeout = 60 * timeout_min)
  
  # Make sure geometry exists and CRS alignments are consistent
  if (!inherits(sp_dat, "sf")) {
    stop("sp_dat must be an sf object with a geometry column.")
  }
  if (sf::st_crs(sp_dat) != sf::st_crs(study_boundary)) {
    study_boundary <- sf::st_transform(study_boundary, sf::st_crs(sp_dat))
  }
  
  # ------------------------------------------------------------
  # Mesh + SPDEs
  # ------------------------------------------------------------
  
  hull <- fmesher::fm_extensions(
    study_boundary,
    convex  = mesh_convex,
    concave = mesh_concave
  )
  
  mesh_spatial <- fmesher::fm_mesh_2d_inla(
    loc = sf::st_as_sfc(sp_dat),
    boundary = hull,
    max.edge = mesh_max_edge,
    cutoff = mesh_cutoff,
    crs = sf::st_crs(sp_dat)
  )
  
  matern_abund <- INLA::inla.spde2.pcmatern(
    mesh_spatial,
    prior.range = prior_range_abund,
    prior.sigma = prior_sigma_abund,
    constr = TRUE
  )
  
  # ------------------------------------------------------------
  # 1D effects: Hours since sunrise + Day-of-year
  # ------------------------------------------------------------
  
  sp_dat$Hours_Since_Sunrise <- ensure_numeric(sp_dat$Hours_Since_Sunrise)
  sp_dat$days_since_june15   <- ensure_numeric(sp_dat$days_since_june15)
  
  # Hours-since-sunrise mesh
  HSS_range <- range(sp_dat$Hours_Since_Sunrise, na.rm = TRUE)
  HSS_meshpoints <- seq(HSS_range[1] - 1, HSS_range[2] + 1, length.out = 41)
  HSS_mesh1D <- INLA::inla.mesh.1d(HSS_meshpoints, boundary = "free")
  HSS_spde <- INLA::inla.spde2.pcmatern(
    HSS_mesh1D,
    prior.range = HSS_prior_range,
    prior.sigma = HSS_prior_sigma,
    constr = TRUE
  )
  
  # DOY mesh (days since June 15)
  DOY_range <- range(sp_dat$days_since_june15, na.rm = TRUE)
  DOY_meshpoints <- seq(DOY_range[1] - 5, DOY_range[2] + 5, length.out = 41)
  DOY_mesh1D <- INLA::inla.mesh.1d(DOY_meshpoints, boundary = "free")
  DOY_spde_global <- INLA::inla.spde2.pcmatern(
    DOY_mesh1D,
    prior.range = DOY_prior_range,
    prior.sigma = DOY_prior_sigma,
    constr = TRUE
  )
  
  # ------------------------------------------------------------
  # IID effect for atlas squares
  # ------------------------------------------------------------
  
  pc_prec <- list(prior = "pcprec", param = kappa_pcprec)
  
  # ------------------------------------------------------------
  # Covariate components (linear by default)
  # covariates df must have: covariate, beta, model, mean, prec
  # Creates components: Beta<beta>_<cov>(1, model="linear", ...)
  # And formula terms: Beta... * I(cov^beta)
  # ------------------------------------------------------------
  
  covariates <- covariates %>%
    mutate(
      components = paste0(
        "Beta", beta, "_", covariate,
        '(1, model="', model,
        '", mean.linear=', mean,
        ", prec.linear=", prec,
        ")"
      ),
      formula = paste0("Beta", beta, "_", covariate, "*I(", covariate, "^", beta, ")")
    )
  
  covar_components_str <- if (nrow(covariates) > 0) paste(covariates$components, collapse = " + ") else ""
  covar_terms_str <- if (nrow(covariates) > 0) paste(covariates$formula, collapse = " + ") else ""
  
  # ------------------------------------------------------------
  # Model components (bru)
  # ------------------------------------------------------------
  
  components_str <- paste0(
    '~ Intercept(1) +
       effect_ARU(1, model="linear", mean.linear=0, prec.linear=100) +
       spde_abund(main = geometry, model = matern_abund) +
       HSS(main = Hours_Since_Sunrise, model = HSS_spde) +
       DOY_global(main = days_since_june15, model = DOY_spde_global) +
       kappa(square_atlas, model="iid", constr=TRUE, hyper=list(prec=pc_prec))',
    if (nchar(covar_components_str) > 0) paste0(" + ", covar_components_str) else ""
  )
  
  model_components <- as.formula(components_str)
  
  # ------------------------------------------------------------
  # Likelihood formula (shared by ARU / point counts)
  # ------------------------------------------------------------
  
  # You can keep ARU * effect_ARU term even if ARU is all zeros; it just becomes irrelevant.
  base_formula_str <- paste0(
    "count ~
      Intercept +
      ARU * effect_ARU +
      HSS +
      DOY_global +
      kappa +
      spde_abund",
    if (nchar(covar_terms_str) > 0) paste0(" + ", covar_terms_str) else ""
  )
  model_formula <- as.formula(base_formula_str)
  
  # ------------------------------------------------------------
  # Likelihood setup (two likelihoods: ARU==0 and ARU==1)
  # ------------------------------------------------------------
  
  make_like <- function(data_subset) {
    if (family == "poisson") {
      inlabru::like(
        family = "poisson",
        formula = model_formula,
        data = data_subset
      )
    } else {
      # negative binomial with PC prior on theta (size)
      lambda <- calibrate_pc_lambda(
        target_prob = nb_pc_target_prob,
        threshold_theta = nb_pc_threshold_theta
      )
      
      inlabru::like(
        family = "nbinomial",
        formula = model_formula,
        data = data_subset,
        control.family = list(
          hyper = list(theta = list(
            prior = "pc.gamma",
            param = c(lambda)
          ))
        )
      )
    }
  }
  
  dat_pc  <- sp_dat %>% dplyr::filter(ARU == 0)
  dat_aru <- sp_dat %>% dplyr::filter(ARU == 1)
  
  # If one of the subsets is empty, bru() can still run with a single likelihood
  likes <- list()
  if (nrow(dat_pc)  > 0) likes <- c(likes, list(make_like(dat_pc)))
  if (nrow(dat_aru) > 0) likes <- c(likes, list(make_like(dat_aru)))
  
  if (length(likes) == 0) stop("No data rows for ARU==0 or ARU==1 after filtering; cannot fit model.")
  
  # ------------------------------------------------------------
  # Fit with optional retries
  # ------------------------------------------------------------
  
  bru_opts <- list(
    inla.mode = inla_mode,
    control.compute = list(waic = waic, cpo = cpo),
    control.inla = list(
      int.strategy = int_strategy,
      strategy = strategy
    ),
    bru_verbose = bru_verbose
  )
  
  attempt <- 0
  last_err <- NULL
  fit <- NULL
  
  while (attempt <= retry) {
    attempt <- attempt + 1
    
    fit <- try(
      do.call(
        inlabru::bru,
        c(
          list(components = model_components),
          likes,
          list(options = bru_opts)
        )
      ),
      silent = TRUE
    )
    
    if (!inherits(fit, "try-error")) break
    last_err <- fit
    fit <- NULL
  }
  
  if (is.null(fit)) {
    stop("bru() failed after ", attempt, " attempt(s).\nLast error:\n", as.character(last_err))
  }
  
  fit
}

# Attempt to fit the model twice using family = nbinomial, if that fails default to poisson twice (and then stop if it fails)
fit_with_fallback_single_atlas <- function(dat_train,
                                           study_boundary,
                                           cov_df_sp,
                                           timeout_min,
                                           prior_range_abund,
                                           prior_sigma_abund,
                                           retry_fit = TRUE,
                                           nb_pc_target_prob = 0.5,
                                           nb_pc_threshold_theta = 5,
                                           family_order = c("nbinomial", "poisson"),
                                           max_tries_per_family = 2,
                                           verbose = TRUE) {
  
  last_error <- NULL
  
  for (fam in family_order) {
    for (attempt in seq_len(max_tries_per_family)) {
      if (verbose) message("  Fitting family=", fam, " attempt ", attempt, "/", max_tries_per_family)
      
      mod_try <- try(
        fit_inla_testing(
          sp_dat = dat_train,
          study_boundary = study_boundary,
          covariates = cov_df_sp,
          timeout_min = timeout_min,
          family = fam,
          
          # NB-only args (safe to pass even if fam = poisson if your function ignores them;
          # if not, wrap them in if(fam=="nbinomial") logic inside fit_inla_testing)
          nb_pc_target_prob = nb_pc_target_prob,
          nb_pc_threshold_theta = nb_pc_threshold_theta,
          
          prior_range_abund = prior_range_abund,
          prior_sigma_abund = prior_sigma_abund,
          retry = retry_fit
        ),
        silent = TRUE
      )
      
      if (!inherits(mod_try, "try-error") && !is.null(mod_try)) {
        return(list(mod = mod_try, family_used = fam, attempts = attempt))
      }
      
      # store last error for debugging/logging
      last_error <- mod_try
    }
  }
  
  # complete failure
  list(mod = NULL, family_used = NA_character_, attempts = NA_integer_, last_error = last_error)
}


# Extract the error family for a fitted model object
get_bru_family <- function(mod) {
  fam <- tryCatch(
    summary(mod)$bru_info$lhoods[[1]]$family,
    error = function(e) NA_character_
  )
  fam <- tolower(fam)
  if (is.na(fam) || !nzchar(fam)) {
    stop("Could not extract family from model via summary(mod)$bru_info$lhoods[[1]]$family")
  }
  fam
}


# Extract size parameters
get_nb_sizes_from_summary <- function(mod) {
  
  hp <- summary(mod)$inla$hyperpar
  
  idx <- grep("^size for the nbinomial observations", rownames(hp))
  if (length(idx) == 0) stop("No nbinomial size hyperparameter found in summary(mod)$hyperpar")
  
  # use posterior median by default (robust)
  sizes <- hp[idx, "0.5quant"]
  as.numeric(sizes)
}

# Simulate new observations (e.g., for posterior predictive checks)
sim_yrep_draws <- function(mu_draws, dat_test, mod) {
  
  fam <- get_bru_family(mod)
  
  # Poisson observations
  if (fam == "poisson") {
    
    return(matrix(
      stats::rpois(length(mu_draws), lambda = as.vector(mu_draws)),
      nrow = nrow(mu_draws),
      ncol = ncol(mu_draws)
    ))
  }
  
  # Negative binomial observations
  if (startsWith(fam, "nbinomial")) {
    nb_sizes <- get_nb_sizes_from_summary(mod)
    
    # ASSUMPTION: size[1] = ARU==0, size[2] = ARU==1
    if (length(nb_sizes) == 1) {
      size_vec <- rep(nb_sizes, nrow(dat_test))
    } else {
      stopifnot(all(dat_test$ARU %in% c(0,1)), !anyNA(dat_test$ARU))
      size_vec <- nb_sizes[dat_test$ARU + 1]
    }
    
    return(
      matrix(
        stats::rnbinom(n = length(mu_draws),
                       mu = as.vector(mu_draws),
                       size = rep(size_vec, times = ncol(mu_draws))),
        nrow = nrow(mu_draws),
        ncol = ncol(mu_draws)
      )
      
    )
  }
  
  stop("Unsupported family: ", fam)
}


sp_filename <- function(sp_english) {
  sp_english %>%
    str_to_lower() %>%
    str_replace_all("[^a-z0-9]+", "_") %>%
    str_replace_all("^_|_$", "")
}

load_or_empty_list <- function(path) {
  if (file.exists(path)) readRDS(path) else list()
}

# Save RDS
save_atomic <- function(obj, path) {
  tmp <- paste0(path, ".tmp")
  saveRDS(obj, tmp)
  file.rename(tmp, path)
}

# Create “derived” covariates consistently across surveys + grids
add_derived_covariates <- function(surveys, grid2, grid3,
                                   south_bcr = c(12, 13), north_bcr = c(7, 8)) {
  
  make_bits <- function(df) {
    if (!("water_river" %in% names(df))) df$water_river <- NA_real_
    if (!("BCR" %in% names(df))) stop("BCR missing; required for lc_* split.")
    
    df %>%
      mutate(
        on_river = as.integer(!is.na(water_river) & water_river > 0),
        on_road = as.integer(!is.na(road) & road > 0),
        lc_8S  = ifelse(BCR %in% south_bcr, lc_8, 0),
        lc_8N  = ifelse(BCR %in% north_bcr, lc_8, 0),
        lc_9S  = ifelse(BCR %in% south_bcr, lc_9, 0),
        lc_9N  = ifelse(BCR %in% north_bcr, lc_9, 0),
        lc_10S = ifelse(BCR %in% south_bcr, lc_10, 0),
        lc_10N = ifelse(BCR %in% north_bcr, lc_10, 0)
      )
  }
  
  surveys <- make_bits(surveys)
  grid2   <- make_bits(grid2)
  grid3   <- make_bits(grid3)
  
  list(surveys = surveys, grid2 = grid2, grid3 = grid3)
}

# Build a covariate specification dataframe (simple linear priors)
make_cov_df <- function(covars) {
  tibble(
    covariate = covars,
    beta = 1,
    sd_linear = 3,
    model = "linear",
    mean = 0,
    prec = 1 / (sd_linear^2)
  )
}

# Prediction formula builder: returns an inlabru generate() formula that outputs one column "eta"
make_pred_formula <- function(cov_df = NULL, include_kappa = FALSE, include_aru = FALSE) {
  cov_terms <- character(0)
  
  if (!is.null(cov_df) && nrow(cov_df) > 0) {
    cov_terms <- cov_df %>%
      dplyr::mutate(term = paste0("Beta", beta, "_", covariate, "*I(", covariate, "^", beta, ")")) %>%
      dplyr::pull(term)
  }
  
  base_terms <- c(
    "Intercept",
    "spde_abund",
    "Atlas3 * effect_Atlas3",
    "Atlas3 * spde_change",
    "DOY_global",
    "HSS"
  )
  
  if (include_aru) {
    base_terms <- c(base_terms, "ARU * effect_ARU")
  }
  
  eta_expr <- paste(c(base_terms, cov_terms), collapse = " + ")
  
  if (include_kappa) {
    eta_expr <- paste0(eta_expr, " + kappa")
  }
  
  as.formula(paste0("~ data.frame(eta = ", eta_expr, ")"))
}

# ------------------------------------------------------------------------------
# Updated functions that fit a centered atlas model
# ------------------------------------------------------------------------------

fit_inla_multi_atlas2 <- function(sp_dat,
                                  study_boundary,
                                  covariates,
                                  timeout_min = 15,
                                  
                                  # SPDE priors
                                  prior_range_abund  = c(150, 0.1),
                                  prior_sigma_abund  = c(0.5, 0.1),
                                  prior_range_change = c(500, 0.1),
                                  prior_sigma_change = c(0.1, 0.1),
                                  
                                  # Mesh controls
                                  mesh_max_edge = c(100, 200),
                                  mesh_cutoff = 50,
                                  mesh_convex = c(50, 150),
                                  mesh_concave = c(50, 150),
                                  
                                  # 1D smoothers
                                  HSS_prior_range  = c(2, 0.1),
                                  HSS_prior_sigma  = c(1, 0.5),
                                  
                                  DOY_prior_range  = c(30, 0.1),
                                  DOY_prior_sigma  = c(1, 0.5),
                                  
                                  # atlas square iid effect
                                  kappa_pcprec_shared = c(3, 0.1),
                                  kappa_pcprec_diff   = c(0.1, 0.1),
                                  
                                  # likelihood
                                  family = c("poisson", "nbinomial"),
                                  
                                  # NB PC prior config (only used if family="nbinomial")
                                  nb_pc_target_prob = 0.5,
                                  nb_pc_threshold_theta = 5,
                                  
                                  # inlabru/INLA options
                                  inla_mode = "experimental",
                                  int_strategy = "eb",
                                  strategy = "simplified.laplace",
                                  bru_verbose = 4,
                                  waic = FALSE,
                                  cpo = FALSE,
                                  retry = 0) {
  
  family <- match.arg(family)
  
  stopifnot(is.data.frame(sp_dat) || inherits(sp_dat, "sf"))
  stopifnot(inherits(study_boundary, "sf"))
  
  required_cols <- c("count", "Atlas3", "Hours_Since_Sunrise", "days_since_june15",
                     "square_atlas", "geometry")
  missing <- setdiff(required_cols, names(sp_dat))
  if (length(missing) > 0) stop("sp_dat is missing required columns: ", paste(missing, collapse = ", "))
  
  # covariates is a data.frame with at least: covariate, beta, model, mean, prec
  if (!is.data.frame(covariates)) stop("covariates must be a data.frame")
  needed_cov_cols <- c("covariate", "beta", "model", "mean", "prec")
  cov_missing <- setdiff(needed_cov_cols, names(covariates))
  if (length(cov_missing) > 0) stop("covariates is missing columns: ", paste(cov_missing, collapse = ", "))
  
  # Timeout (seconds)
  INLA::inla.setOption(inla.timeout = 60 * timeout_min)
  
  # Make sure geometry exists and CRS alignments are consistent
  if (!inherits(sp_dat, "sf")) {
    stop("sp_dat must be an sf object with a geometry column.")
  }
  if (sf::st_crs(sp_dat) != sf::st_crs(study_boundary)) {
    study_boundary <- sf::st_transform(study_boundary, sf::st_crs(sp_dat))
  }
  
  # ------------------------------------------------------------
  # Mesh + SPDEs
  # ------------------------------------------------------------
  
  hull <- fmesher::fm_extensions(
    study_boundary,
    convex  = mesh_convex,
    concave = mesh_concave
  )
  
  mesh_spatial <- fmesher::fm_mesh_2d_inla(
    loc = sf::st_as_sfc(sp_dat),
    boundary = hull,
    max.edge = mesh_max_edge,
    cutoff = mesh_cutoff,
    crs = sf::st_crs(sp_dat)
  )
  
  # Mean/shared field (persistent abundance surface)
  matern_mean <- INLA::inla.spde2.pcmatern(
    mesh_spatial,
    prior.range = prior_range_abund,
    prior.sigma = prior_sigma_abund,
    constr = TRUE
  )
  
  # Diff field (Atlas3 - Atlas2 on log scale, under centered coding)
  matern_diff <- INLA::inla.spde2.pcmatern(
    mesh_spatial,
    prior.range = prior_range_change,
    prior.sigma = prior_sigma_change,
    constr = TRUE
  )
  
  # ------------------------------------------------------------
  # 1D effects: Hours since sunrise + Day-of-year
  # ------------------------------------------------------------
  
  sp_dat$Hours_Since_Sunrise <- ensure_numeric(sp_dat$Hours_Since_Sunrise)
  sp_dat$days_since_june15   <- ensure_numeric(sp_dat$days_since_june15)
  sp_dat$Atlas3_c              <- ensure_numeric(sp_dat$Atlas3_c)
  
  # Hours-since-sunrise mesh
  HSS_range <- range(sp_dat$Hours_Since_Sunrise, na.rm = TRUE)
  HSS_meshpoints <- seq(HSS_range[1] - 1, HSS_range[2] + 1, length.out = 41)
  HSS_mesh1D <- INLA::inla.mesh.1d(HSS_meshpoints, boundary = "free")
  HSS_spde <- INLA::inla.spde2.pcmatern(
    HSS_mesh1D,
    prior.range = HSS_prior_range,
    prior.sigma = HSS_prior_sigma,
    constr = TRUE
  )
  
  # DOY mesh (days since June 15)
  DOY_range <- range(sp_dat$days_since_june15, na.rm = TRUE)
  DOY_meshpoints <- seq(DOY_range[1] - 5, DOY_range[2] + 5, length.out = 41)
  DOY_mesh1D <- INLA::inla.mesh.1d(DOY_meshpoints, boundary = "free")
  DOY_spde_global <- INLA::inla.spde2.pcmatern(
    DOY_mesh1D,
    prior.range = DOY_prior_range,
    prior.sigma = DOY_prior_sigma,
    constr = TRUE
  )
  
  # ------------------------------------------------------------
  # IID effects for atlas squares
  # ------------------------------------------------------------
  
  pc_prec_shared <- list(prior = "pcprec", param = kappa_pcprec_shared)
  pc_prec_diff <- list(prior = "pcprec", param = kappa_pcprec_diff)
  
  # ------------------------------------------------------------
  # Covariate components (linear by default)
  # ------------------------------------------------------------
  
  covariates <- covariates %>%
    dplyr::mutate(
      components = paste0(
        "Beta", beta, "_", covariate,
        '(1, model="', model,
        '", mean.linear=', mean,
        ", prec.linear=", prec,
        ")"
      ),
      formula = paste0("Beta", beta, "_", covariate, "*I(", covariate, "^", beta, ")")
    )
  
  covar_components_str <- if (nrow(covariates) > 0) paste(covariates$components, collapse = " + ") else ""
  covar_terms_str <- if (nrow(covariates) > 0) paste(covariates$formula, collapse = " + ") else ""
  
  # ------------------------------------------------------------
  # Model components (bru)
  # ------------------------------------------------------------
  
  components_str <- paste0(
    '~ Intercept(1) +
    
       # NOTE: effect_Atlas3 now multiplies Atlas3_c in the likelihood,
       # making this parameter directly interpretable as global (Atlas3 - Atlas2).
       effect_Atlas3(1, model="linear", mean.linear=0, prec.linear=1) +
       effect_ARU(1, model="linear", mean.linear=0, prec.linear=100) +
       spde_mean(main = geometry, model = matern_mean) +
       spde_diff(main = geometry, model = matern_diff) +
       HSS(main = Hours_Since_Sunrise, model = HSS_spde) +
       DOY_global(main = days_since_june15, model = DOY_spde_global) +
       kappa_diff(square_atlas, model="iid", constr=TRUE, hyper=list(prec=pc_prec_diff)) +
       kappa_shared(square_idx, model="iid", constr=TRUE, hyper=list(prec=pc_prec_shared))',
    if (nchar(covar_components_str) > 0) paste0(" + ", covar_components_str) else ""
  )
  
  model_components <- as.formula(components_str)
  
  # ------------------------------------------------------------
  # Likelihood formula (shared by ARU / point counts)
  # ------------------------------------------------------------
  
  base_formula_str <- paste0(
    "count ~
      Intercept +
      ARU * effect_ARU +
      HSS +
      DOY_global +
      kappa_diff +
      kappa_shared +
      spde_mean +
      Atlas3_c * spde_diff +
      Atlas3_c * effect_Atlas3",
    if (nchar(covar_terms_str) > 0) paste0(" + ", covar_terms_str) else ""
  )
  model_formula <- as.formula(base_formula_str)
  
  # ------------------------------------------------------------
  # Likelihood setup (two likelihoods: ARU==0 and ARU==1)
  # ------------------------------------------------------------
  
  make_like <- function(data_subset) {
    if (family == "poisson") {
      inlabru::like(
        family = "poisson",
        formula = model_formula,
        data = data_subset
      )
    } else {
      lambda <- calibrate_pc_lambda(
        target_prob = nb_pc_target_prob,
        threshold_theta = nb_pc_threshold_theta
      )
      
      inlabru::like(
        family = "nbinomial",
        formula = model_formula,
        data = data_subset,
        control.family = list(
          hyper = list(theta = list(
            prior = "pc.gamma",
            param = c(lambda)
          ))
        )
      )
    }
  }
  
  dat_pc  <- sp_dat %>% dplyr::filter(ARU == 0)
  dat_aru <- sp_dat %>% dplyr::filter(ARU == 1)
  
  likes <- list()
  if (nrow(dat_pc)  > 0) likes <- c(likes, list(make_like(dat_pc)))
  if (nrow(dat_aru) > 0) likes <- c(likes, list(make_like(dat_aru)))
  
  if (length(likes) == 0) stop("No data rows for ARU==0 or ARU==1 after filtering; cannot fit model.")
  
  # ------------------------------------------------------------
  # Fit with optional retries
  # ------------------------------------------------------------
  
  bru_opts <- list(
    inla.mode = inla_mode,
    control.compute = list(waic = waic, cpo = cpo),
    control.inla = list(
      int.strategy = int_strategy,
      strategy = strategy
    ),
    bru_verbose = bru_verbose
  )
  
  attempt <- 0
  last_err <- NULL
  fit <- NULL
  
  while (attempt <= retry) {
    attempt <- attempt + 1
    
    fit <- try(
      do.call(
        inlabru::bru,
        c(
          list(components = model_components),
          likes,
          list(options = bru_opts)
        )
      ),
      silent = TRUE
    )
    
    if (!inherits(fit, "try-error")) break
    last_err <- fit
    fit <- NULL
  }
  
  if (is.null(fit)) {
    stop("bru() failed after ", attempt, " attempt(s).\nLast error:\n", as.character(last_err))
  }
  
  fit
}

# Prediction formula builder: returns an inlabru generate() formula that outputs one column "eta"
make_pred_formula2 <- function(cov_df = NULL, include_kappa = FALSE, include_aru = FALSE) {
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
    "HSS"
  )
  
  if (include_aru) {
    base_terms <- c(base_terms, "ARU * effect_ARU")
  }
  
  eta_expr <- paste(c(base_terms, cov_terms), collapse = " + ")
  
  if (include_kappa) {
    eta_expr <- paste0(eta_expr, " + kappa")
  }
  
  as.formula(paste0("~ data.frame(eta = ", eta_expr, ")"))
}

# BCR-specific day and hour smooths
fit_inla_multi_atlas3 <- function(sp_dat,
                                  study_boundary,
                                  covariates,
                                  timeout_min = 15,
                                  
                                  # SPDE priors
                                  prior_range_abund  = c(150, 0.1),
                                  prior_sigma_abund  = c(0.5, 0.1),
                                  prior_range_change = c(500, 0.1),
                                  prior_sigma_change = c(0.1, 0.1),
                                  
                                  # Mesh controls
                                  mesh_max_edge = c(100, 200),
                                  mesh_cutoff = 50,
                                  mesh_convex = c(50, 150),
                                  mesh_concave = c(50, 150),
                                  
                                  # 1D smoothers (global + deviations share these priors)
                                  HSS_prior_range  = c(2, 0.1),
                                  HSS_prior_sigma  = c(1, 0.5),
                                  
                                  DOY_prior_range  = c(30, 0.1),
                                  DOY_prior_sigma  = c(1, 0.5),
                                  
                                  # atlas square iid effect
                                  kappa_pcprec_shared = c(3, 0.1),
                                  kappa_pcprec_diff   = c(0.1, 0.1),
                                  
                                  # likelihood
                                  family = c("poisson", "nbinomial"),
                                  
                                  # NB PC prior config (only used if family="nbinomial")
                                  nb_pc_target_prob = 0.5,
                                  nb_pc_threshold_theta = 5,
                                  
                                  # inlabru/INLA options
                                  inla_mode = "experimental",
                                  int_strategy = "eb",
                                  strategy = "simplified.laplace",
                                  bru_verbose = 4,
                                  waic = FALSE,
                                  cpo = FALSE,
                                  retry = 0) {
  
  family <- match.arg(family)
  
  stopifnot(is.data.frame(sp_dat) || inherits(sp_dat, "sf"))
  stopifnot(inherits(study_boundary, "sf"))
  
  # ---- Required columns
  required_cols <- c(
    "count",
    "Atlas3",
    "Atlas3_c",
    "Hours_Since_Sunrise",
    "days_since_june15",
    "BCR_idx",          # <-- NEW (integer 1..nBCR)
    "square_atlas",
    "square_idx",
    "geometry"
  )
  missing <- setdiff(required_cols, names(sp_dat))
  if (length(missing) > 0) stop("sp_dat is missing required columns: ", paste(missing, collapse = ", "))
  
  # covariates is a data.frame with at least: covariate, beta, model, mean, prec
  if (!is.data.frame(covariates)) stop("covariates must be a data.frame")
  needed_cov_cols <- c("covariate", "beta", "model", "mean", "prec")
  cov_missing <- setdiff(needed_cov_cols, names(covariates))
  if (length(cov_missing) > 0) stop("covariates is missing columns: ", paste(cov_missing, collapse = ", "))
  
  # ---- Timeout (seconds)
  INLA::inla.setOption(inla.timeout = 60 * timeout_min)
  
  # ---- Make sure geometry exists and CRS alignments are consistent
  if (!inherits(sp_dat, "sf")) {
    stop("sp_dat must be an sf object with a geometry column.")
  }
  if (sf::st_crs(sp_dat) != sf::st_crs(study_boundary)) {
    study_boundary <- sf::st_transform(study_boundary, sf::st_crs(sp_dat))
  }
  
  # ------------------------------------------------------------
  # Mesh + SPDEs
  # ------------------------------------------------------------
  
  hull <- fmesher::fm_extensions(
    study_boundary,
    convex  = mesh_convex,
    concave = mesh_concave
  )
  
  mesh_spatial <- fmesher::fm_mesh_2d_inla(
    loc = sf::st_as_sfc(sp_dat),
    boundary = hull,
    max.edge = mesh_max_edge,
    cutoff = mesh_cutoff,
    crs = sf::st_crs(sp_dat)
  )
  
  # Mean/shared field (persistent abundance surface)
  matern_mean <- INLA::inla.spde2.pcmatern(
    mesh_spatial,
    prior.range = prior_range_abund,
    prior.sigma = prior_sigma_abund,
    constr = TRUE
  )
  
  # Diff field (Atlas3 - Atlas2 on log scale, under centered coding)
  matern_diff <- INLA::inla.spde2.pcmatern(
    mesh_spatial,
    prior.range = prior_range_change,
    prior.sigma = prior_sigma_change,
    constr = TRUE
  )
  
  # ------------------------------------------------------------
  # 1D effects: Hours since sunrise + Day-of-year
  #   Now: GLOBAL smooth + BCR-SPECIFIC DEVIATIONS (replicates)
  # ------------------------------------------------------------
  
  sp_dat$Hours_Since_Sunrise <- ensure_numeric(sp_dat$Hours_Since_Sunrise)
  sp_dat$days_since_june15   <- ensure_numeric(sp_dat$days_since_june15)
  sp_dat$Atlas3_c            <- ensure_numeric(sp_dat$Atlas3_c)
  sp_dat$BCR_idx             <- ensure_numeric(sp_dat$BCR_idx)
  
  # Hours-since-sunrise mesh
  HSS_range <- range(sp_dat$Hours_Since_Sunrise, na.rm = TRUE)
  HSS_meshpoints <- seq(HSS_range[1] - 1, HSS_range[2] + 1, length.out = 41)
  HSS_mesh1D <- INLA::inla.mesh.1d(HSS_meshpoints, boundary = "free")
  HSS_spde <- INLA::inla.spde2.pcmatern(
    HSS_mesh1D,
    prior.range = HSS_prior_range,
    prior.sigma = HSS_prior_sigma,
    constr = TRUE
  )
  
  # DOY mesh (days since June 15)
  DOY_range <- range(sp_dat$days_since_june15, na.rm = TRUE)
  DOY_meshpoints <- seq(DOY_range[1] - 5, DOY_range[2] + 5, length.out = 41)
  DOY_mesh1D <- INLA::inla.mesh.1d(DOY_meshpoints, boundary = "free")
  DOY_spde <- INLA::inla.spde2.pcmatern(
    DOY_mesh1D,
    prior.range = DOY_prior_range,
    prior.sigma = DOY_prior_sigma,
    constr = TRUE
  )
  
  # ------------------------------------------------------------
  # IID effects for atlas squares
  # ------------------------------------------------------------
  
  pc_prec_shared <- list(prior = "pcprec", param = kappa_pcprec_shared)
  pc_prec_diff   <- list(prior = "pcprec", param = kappa_pcprec_diff)
  
  # ------------------------------------------------------------
  # Covariate components (linear by default)
  # ------------------------------------------------------------
  
  covariates <- covariates %>%
    dplyr::mutate(
      components = paste0(
        "Beta", beta, "_", covariate,
        '(1, model="', model,
        '", mean.linear=', mean,
        ", prec.linear=", prec,
        ")"
      ),
      formula = paste0("Beta", beta, "_", covariate, "*I(", covariate, "^", beta, ")")
    )
  
  covar_components_str <- if (nrow(covariates) > 0) paste(covariates$components, collapse = " + ") else ""
  covar_terms_str      <- if (nrow(covariates) > 0) paste(covariates$formula, collapse = " + ") else ""
  
  # ------------------------------------------------------------
  # Model components (bru)
  #   Key change: HSS_global + HSS_dev(BCR) and DOY_global + DOY_dev(BCR)
  # ------------------------------------------------------------
  # HSS_dev(main = Hours_Since_Sunrise, model = HSS_spde,group = BCR_idx, control.group = list(model = "iid")) +
  components_str <- paste0(
    '~ Intercept(1) +
       
       # Atlas effects
       effect_Atlas3(1, model="linear", mean.linear=0, prec.linear=1) +
       effect_ARU(1, model="linear", mean.linear=0, prec.linear=100) +
       
       # Spatial fields
       spde_mean(main = geometry, model = matern_mean) +
       spde_diff(main = geometry, model = matern_diff) +
       
       # Detectability smoothers: global mean + BCR-specific deviations
       HSS_global(main = Hours_Since_Sunrise, model = HSS_spde) +
       
       DOY_global(main = days_since_june15, model = DOY_spde) +
       DOY_dev(main = days_since_june15, model = DOY_spde,
               group = BCR_idx, control.group = list(model = "iid")) +
       
       # Atlas square iid effects
       kappa_diff(square_atlas, model="iid", constr=TRUE, hyper=list(prec=pc_prec_diff)) +
       kappa_shared(square_idx, model="iid", constr=TRUE, hyper=list(prec=pc_prec_shared))',
    if (nchar(covar_components_str) > 0) paste0(" + ", covar_components_str) else ""
  )
  
  model_components <- as.formula(components_str)
  
  # ------------------------------------------------------------
  # Likelihood formula (shared by ARU / point counts)
  # ------------------------------------------------------------
  # + HSS_dev +
  base_formula_str <- paste0(
    "count ~
      Intercept +
      ARU * effect_ARU +
      
      # Global + deviation parameterization
      HSS_global +
      DOY_global + DOY_dev +
      
      kappa_diff +
      kappa_shared +
      spde_mean +
      Atlas3_c * spde_diff +
      Atlas3_c * effect_Atlas3",
    if (nchar(covar_terms_str) > 0) paste0(" + ", covar_terms_str) else ""
  )
  model_formula <- as.formula(base_formula_str)
  
  # ------------------------------------------------------------
  # Likelihood setup (two likelihoods: ARU==0 and ARU==1)
  # ------------------------------------------------------------
  
  make_like <- function(data_subset) {
    if (family == "poisson") {
      inlabru::like(
        family  = "poisson",
        formula = model_formula,
        data    = data_subset
      )
    } else {
      lambda <- calibrate_pc_lambda(
        target_prob      = nb_pc_target_prob,
        threshold_theta  = nb_pc_threshold_theta
      )
      
      inlabru::like(
        family  = "nbinomial",
        formula = model_formula,
        data    = data_subset,
        control.family = list(
          hyper = list(theta = list(
            prior = "pc.gamma",
            param = c(lambda)
          ))
        )
      )
    }
  }
  
  dat_pc  <- sp_dat %>% dplyr::filter(ARU == 0)
  dat_aru <- sp_dat %>% dplyr::filter(ARU == 1)
  
  likes <- list()
  if (nrow(dat_pc)  > 0) likes <- c(likes, list(make_like(dat_pc)))
  if (nrow(dat_aru) > 0) likes <- c(likes, list(make_like(dat_aru)))
  if (length(likes) == 0) stop("No data rows for ARU==0 or ARU==1 after filtering; cannot fit model.")
  
  # ------------------------------------------------------------
  # Fit with optional retries
  # ------------------------------------------------------------
  
  bru_opts <- list(
    inla.mode = inla_mode,
    control.compute = list(waic = waic, cpo = cpo),
    control.inla = list(
      int.strategy = int_strategy,
      strategy = strategy
    ),
    bru_verbose = bru_verbose
  )
  
  attempt <- 0
  last_err <- NULL
  fit <- NULL
  
  while (attempt <= retry) {
    attempt <- attempt + 1
    
    fit <- try(
      do.call(
        inlabru::bru,
        c(
          list(components = model_components),
          likes,
          list(options = bru_opts)
        )
      ),
      silent = TRUE
    )
    
    if (!inherits(fit, "try-error")) break
    last_err <- fit
    fit <- NULL
  }
  
  if (is.null(fit)) {
    stop("bru() failed after ", attempt, " attempt(s).\nLast error:\n", as.character(last_err))
  }
  
  fit
}

# Prediction formula builder: returns an inlabru generate() formula that outputs one column "eta"
make_pred_formula3 <- function(cov_df = NULL, include_kappa = FALSE, include_aru = FALSE) {
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
    "DOY_dev",
    "HSS_global"
  )
  
  if (include_aru) {
    base_terms <- c(base_terms, "ARU * effect_ARU")
  }
  
  eta_expr <- paste(c(base_terms, cov_terms), collapse = " + ")
  
  if (include_kappa) {
    eta_expr <- paste0(eta_expr, " + kappa_diff + kappa_shared")
  }
  
  as.formula(paste0("~ data.frame(eta = ", eta_expr, ")"))
}

