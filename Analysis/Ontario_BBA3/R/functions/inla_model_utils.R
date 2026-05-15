# ============================================================
# inla_model_utils.R
#
# Purpose
#   Helper functions for the INLA/inlabru multi-atlas modeling workflow.
#   Functions cover small file/name utilities, model fitting, prediction,
#   posterior summaries, and hex-level posterior draw aggregation.
#
# Notes
#   - Active functions are kept above.
#   - Fully commented-out experimental/downstream functions are preserved
#     in an archive block at the bottom of the file.
# ============================================================

# ------------------------------------------------------------
# Small utilities
# ------------------------------------------------------------

`%||%` <- function(x, y) if (is.null(x)) y else x

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


# Save RDS

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

# Fit the joint (multi-atlas) INLA/inlabru model used in the OBBA workflow
#
# This is the main model-fitting routine for the joint Atlas 2 vs Atlas 3
# analysis.
fit_inla_multi_atlas <- function(
    sp_dat,
    study_boundary,
    covariates,
    timeout_min = 15,
    
    # Spatial SPDE priors
    prior_range_abund  = c(200, 0.1),
    prior_sigma_abund  = c(3,   0.1),
    prior_range_change = c(200, 0.1),
    prior_sigma_change = c(0.1, 0.1),
    
    # 1-D SPDE priors
    prior_HSS_range = c(5, 0.9),
    prior_HSS_sigma = c(3, 0.1),
    
    prior_DOY_range_global = c(7, 0.9),
    prior_DOY_sigma_global = c(3, 0.1),
    
    # Atlas-square iid prior
    kappa_pcprec_diff = c(1, 0.1),
    
    # Likelihood
    family = c("poisson", "nbinomial"),
    
    # Negative-binomial PC prior settings
    nb_pc_target_prob     = 0.5,
    nb_pc_threshold_theta = 5,
    
    # inlabru / INLA options
    inla_mode    = "experimental",
    int_strategy = "eb",
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
    "Hours_Since_Sunrise",
    "days_midpoint",
    "ARU",
    "SC",
    "LT",
    "square_atlas",
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
  
  sp_dat$Hours_Since_Sunrise <- ensure_numeric(sp_dat$Hours_Since_Sunrise)
  sp_dat$days_midpoint       <- ensure_numeric(sp_dat$days_midpoint)
  sp_dat$Atlas3_c            <- ensure_numeric(sp_dat$Atlas3_c)
  sp_dat$ARU                 <- ensure_numeric(sp_dat$ARU)
  sp_dat$SC                  <- ensure_numeric(sp_dat$SC)
  sp_dat$LT                  <- ensure_numeric(sp_dat$LT)
  
  # ------------------------------------------------------------
  # 4. Build spatial meshes internally
  # ------------------------------------------------------------
  
  hull <- fmesher::fm_extensions(
    study_boundary,
    convex  = c(50, 200),
    concave = c(10, 200)
  )
  
  mesh_abund <- fmesher::fm_mesh_2d_inla(
    loc      = sf::st_as_sfc(sp_dat),
    boundary = hull,
    max.edge = c(40, 100),
    cutoff   = 40,
    crs      = sf::st_crs(sp_dat)
  )
  
  mesh_chg <- fmesher::fm_mesh_2d_inla(
    loc      = sf::st_as_sfc(sp_dat),
    boundary = hull,
    max.edge = c(40, 100),
    cutoff   = 40,
    crs      = sf::st_crs(sp_dat)
  )
  
  # ------------------------------------------------------------
  # 5. Build spatial SPDE models
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
  # 6. Build 1-D SPDE smoothers for detectability corrections
  # ------------------------------------------------------------
  
  HSS_range <- range(sp_dat$Hours_Since_Sunrise, na.rm = TRUE)
  
  HSS_mesh1D <- INLA::inla.mesh.1d(
    loc = seq(HSS_range[1] - 1, HSS_range[2] + 1, length.out = 11),
    boundary = "free"
  )
  
  HSS_spde <- INLA::inla.spde2.pcmatern(
    mesh        = HSS_mesh1D,
    prior.range = prior_HSS_range,
    prior.sigma = prior_HSS_sigma,
    constr      = TRUE
  )
  
  DOY_range <- range(sp_dat$days_midpoint, na.rm = TRUE)
  
  DOY_mesh1D <- INLA::inla.mesh.1d(
    loc = seq(DOY_range[1] - 5, DOY_range[2] + 5, length.out = 15),
    boundary = "free"
  )
  
  DOY_spde_global <- INLA::inla.spde2.pcmatern(
    mesh        = DOY_mesh1D,
    prior.range = prior_DOY_range_global,
    prior.sigma = prior_DOY_sigma_global,
    constr      = TRUE
  )
  
  # ------------------------------------------------------------
  # 7. Define iid atlas-square prior
  # ------------------------------------------------------------
  
  pc_prec_diff <- list(
    prior = "pcprec",
    param = kappa_pcprec_diff
  )
  
  # ------------------------------------------------------------
  # 8. Build covariate components and formula terms
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
  # 9. Define inlabru components
  # ------------------------------------------------------------
  
  components_str <- paste0(
    "~ Intercept(1) + ",
    
    "effect_Atlas3(
       1,
       model = 'linear',
       mean.linear = 0,
       prec.linear = 1 / ((log(1.5) / 2)^2)
     ) + ",
    
    "effect_ARU(
       1,
       model = 'linear',
       mean.linear = 0,
       prec.linear = 1 / ((log(1.25) / 2)^2)
     ) + ",
    
    "effect_SC(
       1,
       model = 'linear',
       mean.linear = 0,
       prec.linear = 1 / ((log(1.5) / 2)^2)
     ) + ",
    
    "effect_LT(
       1,
       model = 'linear',
       mean.linear = 0,
       prec.linear = 1 / ((log(1.5) / 2)^2)
     ) + ",
    
    "spde_mean(main = geometry, model = matern_mean) + ",
    "spde_diff(main = geometry, model = matern_diff) + ",
    
    "HSS_global(main = Hours_Since_Sunrise, model = HSS_spde) + ",
    "DOY_global(main = days_midpoint, model = DOY_spde_global) + ",
    
    "kappa_diff(
       square_atlas,
       model = 'iid',
       constr = TRUE,
       hyper = list(prec = pc_prec_diff)
     )",
    
    if (nchar(covar_components_str) > 0) {
      paste0(" + ", covar_components_str)
    } else {
      ""
    }
  )
  
  model_components <- stats::as.formula(components_str)
  
  # ------------------------------------------------------------
  # 10. Define linear predictor
  # ------------------------------------------------------------
  
  formula_str <- paste0(
    "count ~
      Intercept +
      HSS_global +
      DOY_global +

      ARU * effect_ARU +
      SC  * effect_SC +
      LT  * effect_LT +

      kappa_diff +
      spde_mean +

      Atlas3_c * spde_diff +
      Atlas3_c * effect_Atlas3",
    
    if (nchar(covar_terms_str) > 0) {
      paste0(" + ", covar_terms_str)
    } else {
      ""
    }
  )
  
  model_formula <- stats::as.formula(formula_str)
  
  # ------------------------------------------------------------
  # 11. Build likelihood
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
  # 12. Set inlabru / INLA options
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
  # 13. Fit model, with optional retries
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
  
  fit
}

summarize_predictions <- function(mu2, mu3) {
  abs_change <- mu3 - mu2
  list(
    OBBA2 = summarize_posterior(mu2, CI_probs = c(0.05, 0.95), prefix = "OBBA2"),
    OBBA3 = summarize_posterior(mu3, CI_probs = c(0.05, 0.95), prefix = "OBBA3"),
    abs_change = summarize_posterior(abs_change, CI_probs = c(0.05, 0.95), prefix = "abs_change")
  )
}


# Prediction formula builder: returns an inlabru generate() formula that outputs one column "eta"

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
    # "DOY_dev",
    "HSS_global"
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
      Hours_Since_Sunrise = 0,
      days_midpoint = 0,
      Atlas3_c = Atlas3 - 0.5,
      BCR_idx = as.integer(BCR)
    )
}

# ------------------------------------------------------------
# Prediction draws via inlabru::generate()
# Returns a named list of matrices, one per output variable.
# Each matrix is nrow(grid) x n.samples
# ------------------------------------------------------------

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




# 

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

# ------------------------------------------------------------
# PC prior lambda calibration for inla.pc.rgamma
# (Useful if you want a PC prior on NB "size"/theta)
# ------------------------------------------------------------

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
# Archived commented-out functions
#
# These functions were already fully commented out in the source file.
# They are intentionally preserved here for possible downstream reuse.
# ============================================================

# ------------------------------------------------------------------------------
# Updated functions that fit a centered atlas model
# ------------------------------------------------------------------------------




#' #' Aggregate posterior draws of pixel-level eta into polygon-level draws
#' #'
#' #' Purpose:
#' #'   Convert pixel-level linear predictor draws (eta) into polygon-level
#' #'   draws of total expected abundance for each atlas period (Mu2, Mu3),
#' #'   plus derived draws (absolute and proportional change).
#' #'
#' #' Returns only the draw matrices (no summaries, no geometry), so this
#' #' can be called inside a species loop and passed to downstream
#' #' summarizers for different hypotheses.
#' #'
#' #' @param eta matrix [n_grid x n_draw] linear predictor draws for the full pred_grid
#' #' @param pred_grid data.frame/tibble with n_grid rows and an atlas indicator column
#' #' @param pix_poly_id vector length n_pixels assigning each pixel to a polygon ID (NA allowed)
#' #' @param poly_ids vector of polygon IDs to include (typically unique(pix_poly_id[!is.na()]))
#' #' @param poly_id_col name used for polygon IDs (only used for dimnames + metadata)
#' #' @param atlas_col name of atlas indicator column in pred_grid
#' #' @param atlas_levels length-2 vector, e.g. c("OBBA2","OBBA3") (order defines Mu2/Mu3)
#' #' @param eps small numeric to prevent division by 0 when computing proportional change
#' #' @return list with:
#' #'   - Mu2, Mu3, abs_change, rel_change matrices [n_poly x n_draw]
#' #'   - meta: list with poly_ids, n_pixels, atlas_levels, atlas_col, eps
#' aggregate_polygon_draws_from_eta <- function(eta,
#'                                              pred_grid,
#'                                              pix_poly_id,
#'                                              poly_ids,
#'                                              poly_id_col = "poly_id",
#'                                              atlas_col = "Atlas",
#'                                              atlas_levels = c("OBBA2", "OBBA3"),
#'                                              eps = 1e-9) {
#'   stopifnot(is.matrix(eta))
#'   stopifnot(nrow(eta) == nrow(pred_grid))
#'   stopifnot(length(atlas_levels) == 2)
#'   
#'   # Indices for atlas periods in the combined prediction grid
#'   idx2 <- which(pred_grid[[atlas_col]] == atlas_levels[1])
#'   idx3 <- which(pred_grid[[atlas_col]] == atlas_levels[2])
#'   
#'   # Convert eta -> mu for each atlas period (pixel x draw)
#'   mu2 <- exp(eta[idx2, , drop = FALSE])
#'   mu3 <- exp(eta[idx3, , drop = FALSE])
#'   
#'   # Sanity: pix_poly_id should correspond to "one atlas worth" of pixels
#'   n_pix <- nrow(mu2)
#'   if (length(pix_poly_id) != n_pix) {
#'     stop("pix_poly_id length (", length(pix_poly_id),
#'          ") must equal number of pixel rows per atlas (", n_pix, ").")
#'   }
#'   
#'   # Precompute pixel row indices per polygon
#'   rows_by_poly <- lapply(poly_ids, function(pid) which(pix_poly_id == pid))
#'   n_poly <- length(poly_ids)
#'   n_draw <- ncol(mu2)
#'   
#'   # Allocate polygon-level draws
#'   Mu2 <- matrix(NA_real_, nrow = n_poly, ncol = n_draw,
#'                 dimnames = list(as.character(poly_ids), NULL))
#'   Mu3 <- matrix(NA_real_, nrow = n_poly, ncol = n_draw,
#'                 dimnames = list(as.character(poly_ids), NULL))
#'   n_pixels <- integer(n_poly)
#'   
#'   # Aggregate: sum pixel-level mu within polygon for each draw
#'   for (i in seq_len(n_poly)) {
#'     rows <- rows_by_poly[[i]]
#'     n_pixels[i] <- length(rows)
#'     if (length(rows) == 0) next
#'     Mu2[i, ] <- colMeans(mu2[rows, , drop = FALSE])
#'     Mu3[i, ] <- colMeans(mu3[rows, , drop = FALSE])
#'   }
#'   
#'   abs_change <- Mu3 - Mu2
#'   rel_change <- abs_change / pmax(Mu2, eps)
#'   
#'   list(
#'     Mu2 = Mu2,
#'     Mu3 = Mu3,
#'     abs_change = abs_change,
#'     rel_change = rel_change,
#'     meta = list(
#'       poly_ids = poly_ids,
#'       poly_id_col = poly_id_col,
#'       n_pixels = n_pixels,
#'       atlas_col = atlas_col,
#'       atlas_levels = atlas_levels,
#'       eps = eps
#'     )
#'   )
#' }
#' 
#' # ============================================================
#' # 2) Summarize evidence for hypotheses using polygon-level draw matrices
#' # ============================================================
#' 
#' #' Summarize posterior evidence for threshold-based hypotheses from polygon-level draws
#' #'
#' #' Purpose:
#' #'   Given polygon-level posterior draws (e.g., from aggregate_polygon_draws_from_eta()),
#' #'   compute posterior probabilities and (optional) interval summaries for hypotheses
#' #'   such as:
#' #'     - abs_change > T
#' #'     - abs_change < -T
#' #'     - rel_change > t
#' #'     - rel_change < -t
#' #'     - Mu3 / Mu2 > r   (if you add ratio draws)
#' #'
#' #' @param draws list produced by aggregate_polygon_draws_from_eta()
#' #' @param param which draw matrix to test: one of c("abs_change","rel_change","Mu2","Mu3")
#' #' @param threshold numeric threshold T (same units as the chosen param)
#' #' @param prob_level posterior probability cutoff for "support", e.g. 0.95 or 0.975
#' #' @param direction one of c("increase","decrease","two_sided")
#' #'   - "increase": tests Pr(param > +threshold)
#' #'   - "decrease": tests Pr(param < -threshold)
#' #'   - "two_sided": tests both sides and reports p_inc, p_dec and a classification
#' #' @param ci_probs numeric length-2 vector for posterior interval summaries (optional)
#' #' @param include_summary logical; if TRUE include mean/median/quantiles of the param draws
#' #' @return tibble with one row per polygon id, including posterior probabilities and flags
#' summarize_polygon_hypothesis <- function(draws,
#'                                          param = c("abs_change", "rel_change", "Mu2", "Mu3"),
#'                                          threshold,
#'                                          prob_level = 0.95,
#'                                          direction = c("two_sided", "increase", "decrease"),
#'                                          ci_probs = c(0.05, 0.95),
#'                                          include_summary = TRUE) {
#'   param <- match.arg(param)
#'   direction <- match.arg(direction)
#'   
#'   if (!is.list(draws) || is.null(draws[[param]])) {
#'     stop("draws must be a list containing a matrix named '", param, "'.")
#'   }
#'   X <- draws[[param]]
#'   stopifnot(is.matrix(X))
#'   stopifnot(is.numeric(threshold) && length(threshold) == 1)
#'   stopifnot(is.numeric(prob_level) && prob_level > 0 && prob_level < 1)
#'   
#'   poly_ids <- rownames(X)
#'   if (is.null(poly_ids)) {
#'     # fall back if dimnames weren't set
#'     poly_ids <- draws$meta$poly_ids %||% seq_len(nrow(X))
#'   }
#'   
#'   # Posterior probabilities for exceeding thresholds
#'   p_inc <- apply(X, 1, function(v) mean(v >  threshold))
#'   p_dec <- apply(X, 1, function(v) mean(v < -threshold))
#'   
#'   # Summaries of the draw distribution (optional)
#'   out <- tibble::tibble(
#'     poly_id = poly_ids
#'   )
#'   
#'   if (include_summary) {
#'     out <- out %>%
#'       dplyr::mutate(
#'         mean = apply(X, 1, mean),
#'         q50  = apply(X, 1, stats::median),
#'         q_lo = apply(X, 1, function(v) unname(stats::quantile(v, ci_probs[1]))),
#'         q_hi = apply(X, 1, function(v) unname(stats::quantile(v, ci_probs[2])))
#'       )
#'   }
#'   
#'   if (direction == "increase") {
#'     out <- out %>%
#'       dplyr::mutate(
#'         p = p_inc,
#'         support = p >= prob_level
#'       )
#'   } else if (direction == "decrease") {
#'     out <- out %>%
#'       dplyr::mutate(
#'         p = p_dec,
#'         support = p >= prob_level
#'       )
#'   } else {
#'     # two-sided: report both + classify
#'     out <- out %>%
#'       dplyr::mutate(
#'         p_inc = p_inc,
#'         p_dec = p_dec,
#'         support_increase = p_inc >= prob_level,
#'         support_decrease = p_dec >= prob_level,
#'         classification = dplyr::case_when(
#'           support_increase ~ "Increase",
#'           support_decrease ~ "Decrease",
#'           TRUE ~ "None"
#'         )
#'       )
#'   }
#'   
#'   # If metadata includes n_pixels, attach it for downstream filtering
#'   if (!is.null(draws$meta$n_pixels)) {
#'     out$n_pixels <- draws$meta$n_pixels
#'   }
#'   
#'   # Make param/threshold explicit (useful if you bind results for multiple tests)
#'   out$param <- param
#'   out$threshold <- threshold
#'   out$prob_level <- prob_level
#'   
#'   out
#' }
#' 
#' # helper: provide %||% without importing rlang
#' `%||%` <- function(x, y) if (!is.null(x)) x else y
#' 
#' 
#' 
#' 
#' #' Fit a single-atlas INLA/inlabru model (used for stress tests / comparisons)
#' #'
#' #' Fits the same general model structure as the joint model, but for one
#' #' atlas period at a time. This is useful for:
#' #' - sanity checks against the joint model,
#' #' - comparisons of uncertainty,
#' #' - and debugging model components without the change field.
#' #'
#' #' @param sp_dat Survey data (`sf`) for a single atlas.
#' #' @param study_boundary Boundary polygon for mesh.
#' #' @param covariates Fixed-effect covariate table and priors.
#' #' @param timeout_min INLA timeout.
#' #' @param prior_range_abund,prior_sigma_abund SPDE priors.
#' #' @param mesh_* Mesh construction controls.
#' #' @param family Likelihood family.
#' #' @param nb_* NB PC prior settings (if applicable).
#' #' @param inla_mode,int_strategy,strategy,bru_verbose,waic,cpo,retry INLA controls.
#' #' @return A list containing the fitted model and components needed for prediction.
#' fit_inla_single_atlas <- function(sp_dat,
#'                                   study_boundary,
#'                                   covariates,
#'                                   timeout_min = 15,
#'                                   
#'                                   # SPDE priors
#'                                   prior_range_abund  = c(150, 0.1),
#'                                   prior_sigma_abund  = c(0.5, 0.1),
#'                                   
#'                                   # Mesh controls
#'                                   mesh_max_edge = c(100, 200),
#'                                   mesh_cutoff = 50,
#'                                   mesh_convex = c(50, 150),
#'                                   mesh_concave = c(50, 150),
#'                                   
#'                                   # 1D smoothers
#'                                   HSS_prior_range  = c(2, 0.1),
#'                                   HSS_prior_sigma  = c(1, 0.5),
#'                                   
#'                                   DOY_prior_range  = c(30, 0.1),
#'                                   DOY_prior_sigma  = c(1, 0.5),
#'                                   
#'                                   # atlas square iid effect
#'                                   kappa_pcprec = c(1, 0.1),
#'                                   
#'                                   # likelihood
#'                                   family = c("poisson", "nbinomial"),
#'                                   
#'                                   # NB PC prior config (only used if family="nbinomial")
#'                                   nb_pc_target_prob = 0.5,
#'                                   nb_pc_threshold_theta = 5,
#'                                   # inlabru/INLA options
#'                                   inla_mode = "experimental",
#'                                   int_strategy = "eb",
#'                                   strategy = "simplified.laplace",
#'                                   bru_verbose = 4,
#'                                   waic = FALSE,
#'                                   cpo = FALSE,
#'                                   retry = 0) {
#'   
#'   family <- match.arg(family)
#'   
#'   stopifnot(is.data.frame(sp_dat) || inherits(sp_dat, "sf"))
#'   stopifnot(inherits(study_boundary, "sf"))
#'   
#'   required_cols <- c("count", "Hours_Since_Sunrise", "days_rescaled", "square_atlas", "geometry")
#'   missing <- setdiff(required_cols, names(sp_dat))
#'   if (length(missing) > 0) stop("sp_dat is missing required columns: ", paste(missing, collapse = ", "))
#'   
#'   # covariates is a data.frame with at least: covariate, beta, model, mean, prec
#'   if (!is.data.frame(covariates)) stop("covariates must be a data.frame")
#'   needed_cov_cols <- c("covariate", "beta", "model", "mean", "prec")
#'   cov_missing <- setdiff(needed_cov_cols, names(covariates))
#'   if (length(cov_missing) > 0) stop("covariates is missing columns: ", paste(cov_missing, collapse = ", "))
#'   
#'   # Timeout (seconds)
#'   INLA::inla.setOption(inla.timeout = 60 * timeout_min)
#'   
#'   # Make sure geometry exists and CRS alignments are consistent
#'   if (!inherits(sp_dat, "sf")) {
#'     stop("sp_dat must be an sf object with a geometry column.")
#'   }
#'   if (sf::st_crs(sp_dat) != sf::st_crs(study_boundary)) {
#'     study_boundary <- sf::st_transform(study_boundary, sf::st_crs(sp_dat))
#'   }
#'   
#'   # ------------------------------------------------------------
#'   # Mesh + SPDEs
#'   # ------------------------------------------------------------
#'   
#'   hull <- fmesher::fm_extensions(
#'     study_boundary,
#'     convex  = mesh_convex,
#'     concave = mesh_concave
#'   )
#'   
#'   mesh_spatial <- fmesher::fm_mesh_2d_inla(
#'     loc = sf::st_as_sfc(sp_dat),
#'     boundary = hull,
#'     max.edge = mesh_max_edge,
#'     cutoff = mesh_cutoff,
#'     crs = sf::st_crs(sp_dat)
#'   )
#'   
#'   matern_abund <- INLA::inla.spde2.pcmatern(
#'     mesh_spatial,
#'     prior.range = prior_range_abund,
#'     prior.sigma = prior_sigma_abund,
#'     constr = TRUE
#'   )
#'   
#'   # ------------------------------------------------------------
#'   # 1D effects: Hours since sunrise + Day-of-year
#'   # ------------------------------------------------------------
#'   
#'   sp_dat$Hours_Since_Sunrise <- ensure_numeric(sp_dat$Hours_Since_Sunrise)
#'   sp_dat$days_rescaled   <- ensure_numeric(sp_dat$days_rescaled)
#'   
#'   # Hours-since-sunrise mesh
#'   HSS_range <- range(sp_dat$Hours_Since_Sunrise, na.rm = TRUE)
#'   HSS_meshpoints <- seq(HSS_range[1] - 1, HSS_range[2] + 1, length.out = 41)
#'   HSS_mesh1D <- INLA::inla.mesh.1d(HSS_meshpoints, boundary = "free")
#'   HSS_spde <- INLA::inla.spde2.pcmatern(
#'     HSS_mesh1D,
#'     prior.range = HSS_prior_range,
#'     prior.sigma = HSS_prior_sigma,
#'     constr = TRUE
#'   )
#'   
#'   # DOY mesh (days since June 15)
#'   DOY_range <- range(sp_dat$days_rescaled, na.rm = TRUE)
#'   DOY_meshpoints <- seq(DOY_range[1] - 5, DOY_range[2] + 5, length.out = 41)
#'   DOY_mesh1D <- INLA::inla.mesh.1d(DOY_meshpoints, boundary = "free")
#'   DOY_spde_global <- INLA::inla.spde2.pcmatern(
#'     DOY_mesh1D,
#'     prior.range = DOY_prior_range,
#'     prior.sigma = DOY_prior_sigma,
#'     constr = TRUE
#'   )
#'   
#'   # ------------------------------------------------------------
#'   # IID effect for atlas squares
#'   # ------------------------------------------------------------
#'   
#'   pc_prec <- list(prior = "pcprec", param = kappa_pcprec)
#'   
#'   # ------------------------------------------------------------
#'   # Covariate components (linear by default)
#'   # covariates df must have: covariate, beta, model, mean, prec
#'   # Creates components: Beta<beta>_<cov>(1, model="linear", ...)
#'   # And formula terms: Beta... * I(cov^beta)
#'   # ------------------------------------------------------------
#'   
#'   covariates <- covariates %>%
#'     mutate(
#'       components = paste0(
#'         "Beta", beta, "_", covariate,
#'         '(1, model="', model,
#'         '", mean.linear=', mean,
#'         ", prec.linear=", prec,
#'         ")"
#'       ),
#'       formula = paste0("Beta", beta, "_", covariate, "*I(", covariate, "^", beta, ")")
#'     )
#'   
#'   covar_components_str <- if (nrow(covariates) > 0) paste(covariates$components, collapse = " + ") else ""
#'   covar_terms_str <- if (nrow(covariates) > 0) paste(covariates$formula, collapse = " + ") else ""
#'   
#'   # ------------------------------------------------------------
#'   # Model components (bru)
#'   # ------------------------------------------------------------
#'   
#'   components_str <- paste0(
#'     '~ Intercept(1) +
#'        effect_ARU(1, model="linear", mean.linear=0, prec.linear=100) +
#'        spde_abund(main = geometry, model = matern_abund) +
#'        HSS(main = Hours_Since_Sunrise, model = HSS_spde) +
#'        DOY_global(main = days_rescaled, model = DOY_spde_global) +
#'        kappa(square_atlas, model="iid", constr=TRUE, hyper=list(prec=pc_prec))',
#'     if (nchar(covar_components_str) > 0) paste0(" + ", covar_components_str) else ""
#'   )
#'   
#'   model_components <- as.formula(components_str)
#'   
#'   # ------------------------------------------------------------
#'   # Likelihood formula (shared by ARU / point counts)
#'   # ------------------------------------------------------------
#'   
#'   # You can keep ARU * effect_ARU term even if ARU is all zeros; it just becomes irrelevant.
#'   base_formula_str <- paste0(
#'     "count ~
#'       Intercept +
#'       ARU * effect_ARU +
#'       HSS +
#'       DOY_global +
#'       kappa +
#'       spde_abund",
#'     if (nchar(covar_terms_str) > 0) paste0(" + ", covar_terms_str) else ""
#'   )
#'   model_formula <- as.formula(base_formula_str)
#'   
#'   # ------------------------------------------------------------
#'   # Likelihood setup (two likelihoods: ARU==0 and ARU==1)
#'   # ------------------------------------------------------------
#'   
#'   #' Build the inlabru likelihood object for a subset of data
#'   #'
#'   #' Helper that constructs the `like()` call used by `inlabru::bru()`,
#'   #' including coordinate mapping and the chosen likelihood family.
#'   #' @param data_subset Data frame/sf containing response and coordinates.
#'   #' @return An `inlabru` likelihood specification.
#'   make_like <- function(data_subset) {
#'     if (family == "poisson") {
#'       inlabru::like(
#'         family = "poisson",
#'         formula = model_formula,
#'         data = data_subset
#'       )
#'     } else {
#'       # negative binomial with PC prior on theta (size)
#'       lambda <- calibrate_pc_lambda(
#'         target_prob = nb_pc_target_prob,
#'         threshold_theta = nb_pc_threshold_theta
#'       )
#'       
#'       inlabru::like(
#'         family = "nbinomial",
#'         formula = model_formula,
#'         data = data_subset,
#'         control.family = list(
#'           hyper = list(theta = list(
#'             prior = "pc.gamma",
#'             param = c(lambda)
#'           ))
#'         )
#'       )
#'     }
#'   }
#'   
#'   dat_pc  <- sp_dat %>% dplyr::filter(ARU == 0)
#'   dat_aru <- sp_dat %>% dplyr::filter(ARU == 1)
#'   
#'   # If one of the subsets is empty, bru() can still run with a single likelihood
#'   likes <- list()
#'   if (nrow(dat_pc)  > 0) likes <- c(likes, list(make_like(dat_pc)))
#'   if (nrow(dat_aru) > 0) likes <- c(likes, list(make_like(dat_aru)))
#'   
#'   if (length(likes) == 0) stop("No data rows for ARU==0 or ARU==1 after filtering; cannot fit model.")
#'   
#'   # ------------------------------------------------------------
#'   # Fit with optional retries
#'   # ------------------------------------------------------------
#'   
#'   bru_opts <- list(
#'     inla.mode = inla_mode,
#'     control.compute = list(waic = waic, cpo = cpo),
#'     control.inla = list(
#'       int.strategy = int_strategy,
#'       strategy = strategy
#'     ),
#'     bru_verbose = bru_verbose
#'   )
#'   
#'   attempt <- 0
#'   last_err <- NULL
#'   fit <- NULL
#'   
#'   while (attempt <= retry) {
#'     attempt <- attempt + 1
#'     
#'     fit <- try(
#'       do.call(
#'         inlabru::bru,
#'         c(
#'           list(components = model_components),
#'           likes,
#'           list(options = bru_opts)
#'         )
#'       ),
#'       silent = TRUE
#'     )
#'     
#'     if (!inherits(fit, "try-error")) break
#'     last_err <- fit
#'     fit <- NULL
#'   }
#'   
#'   if (is.null(fit)) {
#'     stop("bru() failed after ", attempt, " attempt(s).\nLast error:\n", as.character(last_err))
#'   }
#'   
#'   fit
#' }
#' 

#' 
#' # ------------------------------------------------------------
#' # Posterior summaries (row-wise) for n_row x n_draw matrices
#' # ------------------------------------------------------------
#' 


#' # Extract the error family for a fitted model object
#' 
#' #' Extract the likelihood family used in a bru model
#' #'
#' #' Convenience helper to inspect a fitted inlabru model and return a
#' #' standardized family string ("poisson" / "nbinomial") used by other helpers.
#' #'
#' #' @param mod Fitted bru model.
#' #' @return Character family name (or NA on failure).
#' get_bru_family <- function(mod) {
#'   fam <- tryCatch(
#'     summary(mod)$bru_info$lhoods[[1]]$family,
#'     error = function(e) NA_character_
#'   )
#'   fam <- tolower(fam)
#'   if (is.na(fam) || !nzchar(fam)) {
#'     stop("Could not extract family from model via summary(mod)$bru_info$lhoods[[1]]$family")
#'   }
#'   fam
#' }
#' 
#' 
#' # Extract size parameters
#' 
#' #' Extract Negative Binomial size/overdispersion parameters from model summary
#' #'
#' #' When fitting NB models, INLA stores parameterization details in the summary.
#' #' This helper pulls those values so posterior predictive simulations can use
#' #' the correct NB parameters.
#' #'
#' #' @param mod Fitted bru model.
#' #' @return A list/data.frame of NB parameter values (size/theta, etc.).
#' get_nb_sizes_from_summary <- function(mod) {
#'   
#'   hp <- summary(mod)$inla$hyperpar
#'   
#'   idx <- grep("^size for the nbinomial observations", rownames(hp))
#'   if (length(idx) == 0) stop("No nbinomial size hyperparameter found in summary(mod)$hyperpar")
#'   
#'   # use posterior median by default (robust)
#'   sizes <- hp[idx, "0.5quant"]
#'   as.numeric(sizes)
#' }
#' 
#' # Simulate new observations (e.g., for posterior predictive checks)
#' 
#' #' Simulate posterior predictive replicate counts (y-rep)
#' #'
#' #' Generates posterior predictive draws of observed counts for held-out data,
#' #' given posterior draws of the mean (`mu_draws`) and the fitted likelihood.
#' #' Used in CV scoring and posterior predictive checks.
#' #'
#' #' @param mu_draws Matrix of posterior draws for the mean.
#' #' @param dat_test Held-out data frame containing observed counts and metadata.
#' #' @param mod Fitted bru model (for family/parameters).
#' #' @return A matrix of y-rep draws aligned to rows of `dat_test`. 
#' sim_yrep_draws <- function(mu_draws, dat_test, mod) {
#'   
#'   fam <- get_bru_family(mod)
#'   
#'   # Poisson observations
#'   if (fam == "poisson") {
#'     
#'     return(matrix(
#'       stats::rpois(length(mu_draws), lambda = as.vector(mu_draws)),
#'       nrow = nrow(mu_draws),
#'       ncol = ncol(mu_draws)
#'     ))
#'   }
#'   
#'   # Negative binomial observations
#'   if (startsWith(fam, "nbinomial")) {
#'     nb_sizes <- get_nb_sizes_from_summary(mod)
#'     
#'     # ASSUMPTION: size[1] = ARU==0, size[2] = ARU==1
#'     if (length(nb_sizes) == 1) {
#'       size_vec <- rep(nb_sizes, nrow(dat_test))
#'     } else {
#'       stopifnot(all(dat_test$ARU %in% c(0,1)), !anyNA(dat_test$ARU))
#'       size_vec <- nb_sizes[dat_test$ARU + 1]
#'     }
#'     
#'     return(
#'       matrix(
#'         stats::rnbinom(n = length(mu_draws),
#'                        mu = as.vector(mu_draws),
#'                        size = rep(size_vec, times = ncol(mu_draws))),
#'         nrow = nrow(mu_draws),
#'         ncol = ncol(mu_draws)
#'       )
#'       
#'     )
#'   }
#'   
#'   stop("Unsupported family: ", fam)
#' }
#' 
#' 




#' 
#' # Prediction formula builder: returns an inlabru generate() formula that outputs one column "eta"
#' make_pred_formula <- function(cov_df = NULL, include_kappa = FALSE, include_aru = FALSE) {
#'   cov_terms <- character(0)
#'   
#'   if (!is.null(cov_df) && nrow(cov_df) > 0) {
#'     cov_terms <- cov_df %>%
#'       dplyr::mutate(term = paste0("Beta", beta, "_", covariate, "*I(", covariate, "^", beta, ")")) %>%
#'       dplyr::pull(term)
#'   }
#'   
#'   base_terms <- c(
#'     "Intercept",
#'     "spde_abund",
#'     "Atlas3 * effect_Atlas3",
#'     "Atlas3 * spde_change",
#'     "DOY_global",
#'     "HSS"
#'   )
#'   
#'   if (include_aru) {
#'     base_terms <- c(base_terms, "ARU * effect_ARU")
#'   }
#'   
#'   eta_expr <- paste(c(base_terms, cov_terms), collapse = " + ")
#'   
#'   if (include_kappa) {
#'     eta_expr <- paste0(eta_expr, " + kappa")
#'   }
#'   
#'   as.formula(paste0("~ data.frame(eta = ", eta_expr, ")"))
#' }
