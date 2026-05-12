# ============================================================
# R/functions/cv_utils.R
#
# Purpose:
#   Helpers for cross-validation runs:
#     - scoring (RMSE/MAE/correlations + Poisson deviance)
#     - probabilistic metrics from posterior predictive summaries
#     - safe incremental saving / append
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
})

# Atomic save RDS (same pattern as script 08)

#' Save an R object to disk atomically (write-then-rename)
#'
#' Writes to a temporary file in the same directory and then renames into
#' place. This prevents partially-written `.rds` files if the process is
#' interrupted (common in long model runs).
#'
#' @param obj R object to save.
#' @param path Output path for the `.rds`.
#' @return Invisibly, the output path.
save_atomic <- function(obj, path) {
  tmp <- paste0(path, ".tmp")
  saveRDS(obj, tmp)
  file.rename(tmp, path)
}

# Append rows to an on-disk RDS data.frame without blowing up memory.
# If key_cols provided, will drop duplicates based on those keys.

#' Append rows to an on-disk RDS table (with optional de-duplication)
#'
#' Loads an existing data.frame from `path` (if present), binds `new_rows`,
#' and optionally drops duplicates based on `key_cols`. Saves back using
#' atomic write to avoid corruption.
#'
#' @param path `.rds` file path.
#' @param new_rows data.frame/tibble of rows to append.
#' @param key_cols Optional character vector of columns defining uniqueness.
#' @return The combined table (invisibly), and writes it to `path`.
append_rds_rows <- function(path, new_rows, key_cols = NULL) {
  stopifnot(is.data.frame(new_rows))
  
  if (file.exists(path)) {
    old <- readRDS(path)
    out <- dplyr::bind_rows(old, new_rows)
  } else {
    out <- new_rows
  }
  
  if (!is.null(key_cols)) {
    missing <- setdiff(key_cols, names(out))
    if (length(missing) > 0) stop("append_rds_rows(): missing key_cols: ", paste(missing, collapse = ", "))
    out <- out %>% distinct(across(all_of(key_cols)), .keep_all = TRUE)
  }
  
  save_atomic(out, path)
  invisible(TRUE)
}

#' Root mean squared error
#'
#' Standard regression error metric (lower is better).
#' `rmse = sqrt(mean((y - yhat)^2))`, ignoring NAs.
#'
#' @param y Observed values.
#' @param yhat Predicted values.
#' @return Numeric scalar RMSE.
rmse <- function(y, yhat) sqrt(mean((y - yhat)^2, na.rm = TRUE))

#' Mean absolute error
#'
#' `mae = mean(|y - yhat|)`, ignoring NAs.
#'
#' @param y Observed values.
#' @param yhat Predicted values.
#' @return Numeric scalar MAE.
mae  <- function(y, yhat) mean(abs(y - yhat), na.rm = TRUE)

# Poisson deviance contribution (good baseline even if you fit NB sometimes)
# D = 2 * sum( y*log(y/mu) - (y-mu) ), with convention y=0 -> 0 for y*log(y/mu)

#' Poisson deviance for count predictions
#'
#' Computes the (scaled) deviance contribution comparing observed counts `y`
#' to predicted mean `mu` under a Poisson likelihood. Useful for comparing
#' probabilistic count forecasts.
#'
#' @param y Observed counts (non-negative).
#' @param mu Predicted mean counts (positive).
#' @return Numeric scalar deviance (lower is better).
poisson_deviance <- function(y, mu) {
  y <- as.numeric(y)
  mu <- pmax(as.numeric(mu), 1e-12)
  
  term <- ifelse(y == 0, 0, y * log(y / mu))
  2 * sum(term - (y - mu), na.rm = TRUE)
}



# Coverage indicator for an interval [lo, hi]

#' Prediction interval coverage
#'
#' Fraction of observations that fall within the interval [lo, hi].
#'
#' @param y Observed values.
#' @param lo Lower bound of interval.
#' @param hi Upper bound of interval.
#' @return Numeric scalar in [0,1].
interval_coverage <- function(y, lo, hi) {
  y <- as.numeric(y); lo <- as.numeric(lo); hi <- as.numeric(hi)
  mean(!is.na(y) & !is.na(lo) & !is.na(hi) & (y >= lo) & (y <= hi))
}

# Brier score for binary events (lower is better)

#' Brier score for binary outcomes
#'
#' Mean squared error between observed binary outcomes and predicted
#' probabilities. Lower is better.
#'
#' @param y01 Binary outcomes (0/1).
#' @param p Predicted probabilities in [0,1].
#' @return Numeric scalar Brier score.
brier_score <- function(y01, p) {
  y01 <- as.numeric(y01)
  p <- pmin(pmax(as.numeric(p), 0), 1)
  mean((y01 - p)^2, na.rm = TRUE)
}

# Monte Carlo log score for Poisson mixture:
# log p(y) ~ log( mean_s dpois(y | mu_s) )
# Implemented with log-sum-exp for stability.

#' Monte Carlo Poisson log score using posterior predictive draws
#'
#' Approximates the expected log predictive density for Poisson counts by
#' averaging `dpois(y | mu_draw)` over posterior draws of `mu`.
#'
#' @param y Observed counts.
#' @param mu_draws Numeric vector/matrix of posterior draws for the mean.
#' @return Numeric scalar log score (higher is better).
poisson_log_score_mc <- function(y, mu_draws) {
  y <- as.numeric(y)
  mu_draws <- pmax(as.matrix(mu_draws), 1e-12)
  stopifnot(nrow(mu_draws) == length(y))
  
  lpd_i <- vapply(seq_along(y), function(i) {
    yi <- y[i]
    if (is.na(yi)) return(NA_real_)
    log_p <- stats::dpois(yi, lambda = mu_draws[i, ], log = TRUE)
    m <- max(log_p)
    m + log(mean(exp(log_p - m)))
  }, numeric(1))
  
  mean(lpd_i, na.rm = TRUE)
}

#' Score cross-validation predictions with multiple metrics
#'
#' Convenience wrapper that expects a data.frame of CV predictions and
#' returns a one-row summary with metrics like RMSE, MAE, deviance, and
#' interval coverage if available.
#'
#' @param df Data frame containing observed and predicted columns.
#' @param y_col Name of observed outcome column.
#' @param mu_col Name of predicted mean column.
#' @return A one-row tibble of scoring metrics.
score_cv_predictions <- function(df, y_col = "count_obs", mu_col = "mu_pred") {
  stopifnot(all(c(y_col, mu_col) %in% names(df)))
  
  y <- df[[y_col]]
  mu <- df[[mu_col]]
  
  # correlations can fail if constant vectors; guard with tryCatch
  cor_spear <- tryCatch(cor(y, mu, use = "complete.obs", method = "spearman"), error = function(e) NA_real_)
  cor_pear  <- tryCatch(cor(y, mu, use = "complete.obs", method = "pearson"),  error = function(e) NA_real_)
  
  tibble(
    n = sum(!is.na(y) & !is.na(mu)),
    rmse = rmse(y, mu),
    mae  = mae(y, mu),
    spearman = cor_spear,
    pearson  = cor_pear,
    poisson_dev = poisson_deviance(y, mu)
  )
}

# Read a saved CV plan if it exists; otherwise create it once and save.
# If a plan exists but the user requests additional repeats, generate ONLY the new repeats
# and save them to a separate file (do not overwrite the original plan file).

#' Read an existing CV plan, or create and save a new one
#'
#' In long workflows you want folds to be reproducible across runs. This
#' helper loads a CV plan from disk if it exists; otherwise it creates one
#' (typically from spatial blocks), saves it, and returns it.
#'
#' @param cv_plan_path Path to saved CV plan `.rds`.
#' @param blocks_sf Spatial blocks as `sf` with `block_id`.
#' @param n_folds Number of folds.
#' @param n_repeats Number of repeats.
#' @param seed Seed for reproducibility.
#' @return A CV plan object (data.frame/tibble).
read_cv_plan_or_create <- function(cv_plan_path,
                                   all_surveys,
                                   block_size_km,
                                   n_folds,
                                   n_repeats,
                                   cv_seed,
                                   out_dir) {
  
  if (file.exists(cv_plan_path)) {
    meta <- readRDS(cv_plan_path)
    
    # Basic sanity checks to avoid accidental mismatch across runs
    stopifnot(meta$block_size_km == block_size_km)
    stopifnot(meta$n_folds == n_folds)
    stopifnot(meta$seed == cv_seed)
    
    block_ids   <- meta$block_ids
    block_polys <- meta$block_polys
    cv_plan     <- meta$cv_plan
    
    stopifnot(all(c("rep","block_id","fold") %in% names(cv_plan)))
    
    saved_reps <- sort(unique(cv_plan$rep))
    max_saved <- if (length(saved_reps) == 0) 0 else max(saved_reps)
    
    if (n_repeats > max_saved) {
      extra_reps <- (max_saved + 1):n_repeats
      
      message("CV plan exists; generating additional repeats: ", paste(extra_reps, collapse = ", "))
      
      # Generate a deterministic plan up to n_repeats, then keep only the new reps
      cv_plan_full <- make_block_cv_plan(
        block_ids,
        n_folds = n_folds,
        n_repeats = n_repeats,
        seed = cv_seed
      )
      
      cv_plan_extra <- cv_plan_full %>% dplyr::filter(rep %in% extra_reps)
      
      # Save extras without overwriting the original plan
      extra_path <- file.path(out_dir, "meta",
                              paste0("cv_plan_blocks_", block_size_km, "km_extra_reps_",
                                     min(extra_reps), "_to_", max(extra_reps), ".rds"))
      save_atomic(cv_plan_extra, extra_path)
      message("Saved extra repeats plan: ", extra_path)
      
      cv_plan <- dplyr::bind_rows(cv_plan, cv_plan_extra)
    }
    
    return(list(block_ids = block_ids, block_polys = block_polys, cv_plan = cv_plan))
  }
  
  # No plan exists yet: create and save ONCE
  message("No existing CV plan found; creating and saving: ", cv_plan_path)
  
  stopifnot(inherits(all_surveys, "sf"))
  
  block_ids <- make_spatial_blocks_grid(all_surveys, block_size_km = block_size_km)
  
  cv_plan <- make_block_cv_plan(
    block_ids,
    n_folds = n_folds,
    n_repeats = n_repeats,
    seed = cv_seed
  )
  
  stopifnot(all(c("rep","block_id","fold") %in% names(cv_plan)))
  
  block_polys <- make_block_polygons_grid(
    block_ids = block_ids,
    block_size_km = block_size_km,
    crs = sf::st_crs(all_surveys)
  )
  
  save_atomic(list(
    block_size_km = block_size_km,
    n_folds = n_folds,
    n_repeats = max(cv_plan$rep),
    seed = cv_seed,
    block_ids = block_ids,
    cv_plan = cv_plan,
    block_polys = block_polys
  ), cv_plan_path)
  
  message("Saved CV plan: ", cv_plan_path)
  
  list(block_ids = block_ids, block_polys = block_polys, cv_plan = cv_plan)
}

# Return TRUE if this species/rep/fold has already been completed (i.e., all expected blocks exist in block summaries)

#' Check whether a CV fold output already exists on disk
#'
#' Used to skip expensive re-fitting for folds that have already been run.
#' Typically checks for expected output files under a block/fold directory.
#'
#' @param block_path Base output directory.
#' @param sp_code Species code used in filenames.
#' @param block_size_km Block size used (part of the directory structure).
#' @param rep Repeat index.
#' @param fold Fold index.
#' @param expected_pairs Character vector of expected (atlas) output labels.
#' @return TRUE if outputs exist and appear complete; FALSE otherwise.
fold_already_completed <- function(block_path, sp_code, block_size_km, rep, fold, expected_pairs) {
  if (!file.exists(block_path)) return(FALSE)
  
  x <- try(readRDS(block_path), silent = TRUE)
  if (inherits(x, "try-error") || is.null(x)) return(FALSE)
  if (inherits(x, "sf")) x <- sf::st_drop_geometry(x)
  
  # expected_pairs must have columns block_id and Atlas
  stopifnot(all(c("block_id","Atlas") %in% names(expected_pairs)))
  
  x2 <- x %>%
    dplyr::filter(
      sp_code == !!sp_code,
      block_size_km == !!block_size_km,
      rep == !!rep,
      fold == !!fold
    )
  
  if (nrow(x2) == 0) return(FALSE)
  
  done_pairs <- x2 %>% dplyr::distinct(block_id, Atlas)
  
  expected_pairs <- expected_pairs %>%
    dplyr::mutate(block_id = as.character(block_id), Atlas = as.character(Atlas)) %>%
    dplyr::distinct(block_id, Atlas)
  
  done_pairs <- done_pairs %>%
    dplyr::mutate(block_id = as.character(block_id), Atlas = as.character(Atlas)) %>%
    dplyr::distinct(block_id, Atlas)
  
  # Check whether every expected (block_id, Atlas) is present in done_pairs
  nrow(dplyr::anti_join(expected_pairs, done_pairs, by = c("block_id","Atlas"))) == 0
}

# ------------------------------------------------------------
# Spatial cross-validation helpers
# ------------------------------------------------------------

# Grid-based spatial blocks in projected CRS with km units.
# Returns a character vector block_id of length nrow(surveys_sf).

#' Build spatial xval blocks (regular grid) over survey locations
#'
#' Creates a regular grid of square blocks (in the workflow projection),
#' then assigns each survey point to a `block_id`. These blocks support
#' spatially-structured cross-validation (hold out whole blocks).
#'
#' @param surveys_sf Survey points as an `sf` object.
#' @param block_size_km Block width/height in kilometres (in projected units).
#' @return An `sf` polygon grid with a `block_id` column.
make_spatial_blocks_grid <- function(surveys_sf, block_size_km) {
  stopifnot(inherits(surveys_sf, "sf"))
  stopifnot(is.numeric(block_size_km), length(block_size_km) == 1, block_size_km > 0)
  
  if (sf::st_is_longlat(surveys_sf)) {
    stop("surveys_sf appears to be lon/lat. Transform to a projected CRS (km units expected) first.")
  }
  
  xy <- sf::st_coordinates(surveys_sf)
  if (nrow(xy) != nrow(surveys_sf)) {
    stop("st_coordinates() did not return one row per survey (possible MULTIPOINT/GEOMETRYCOLLECTION?).")
  }
  
  bx <- floor(xy[, 1] / block_size_km)
  by <- floor(xy[, 2] / block_size_km)
  
  paste0("B_", bx, "_", by)
}

# Given block_ids, create a fold plan (repeat x block -> fold).
# Returns data.frame with columns: repeat, block_id, fold.

#' Create a repeated k-fold xval plan over spatial blocks
#'
#' Given a vector of unique `block_id`s, randomly partitions blocks into
#' folds for each repeat. The returned object is typically saved and then
#' used to drive model fitting/prediction loops.
#'
#' @param block_ids Vector of block identifiers.
#' @param n_folds Number of folds per repeat.
#' @param n_repeats Number of repeats.
#' @param seed Random seed for reproducibility.
#' @return A data.frame/tibble describing repeat/fold membership per block.
make_block_cv_plan <- function(block_ids, n_folds = 5, n_repeats = 5, seed = 1) {
  stopifnot(is.vector(block_ids), length(block_ids) > 0)
  stopifnot(is.numeric(n_folds), length(n_folds) == 1, n_folds >= 2)
  stopifnot(is.numeric(n_repeats), length(n_repeats) == 1, n_repeats >= 1)
  stopifnot(is.numeric(seed), length(seed) == 1)
  
  blocks <- sort(unique(as.character(block_ids)))
  n_blocks <- length(blocks)
  if (n_blocks < n_folds) {
    stop("Number of unique blocks (", n_blocks, ") < n_folds (", n_folds, "). Reduce n_folds or block size.")
  }
  
  out <- vector("list", n_repeats)
  
  for (r in seq_len(n_repeats)) {
    set.seed(seed + r - 1)
    
    perm <- sample(blocks, size = n_blocks, replace = FALSE)
    fold_id <- rep(seq_len(n_folds), length.out = n_blocks)
    
    out[[r]] <- data.frame(
      repeat_id = r,          # <-- renamed (fix)
      block_id  = perm,
      fold      = fold_id,
      stringsAsFactors = FALSE
    )
  }
  
  do.call(rbind, out)
}

# ------------------------------------------------------------
# Convert block_id strings ("B_bx_by") to square polygons
# ------------------------------------------------------------

# Returns an sf with columns: block_id, geometry (POLYGON)
# Assumes survey coordinates are in km units.

#' Reconstruct block polygons from block IDs
#'
#' Utility to rebuild the geometry for a block grid when you only have the
#' identifiers (e.g., stored CV plans). It parses block IDs (encoding x/y)
#' to produce an `sf` polygon grid.
#'
#' @param block_ids Vector of block identifiers as produced by `make_spatial_blocks_grid()`.
#' @param block_size_km Block size in km.
#' @param crs CRS for the returned polygons (typically `get_aea_km_crs()`).
#' @return `sf` polygons for each block_id.
make_block_polygons_grid <- function(block_ids, block_size_km, crs) {
  stopifnot(is.vector(block_ids), length(block_ids) > 0)
  stopifnot(is.numeric(block_size_km), length(block_size_km) == 1, block_size_km > 0)
  stopifnot(!is.null(crs))
  
  ids <- sort(unique(as.character(block_ids)))
  
  # Parse bx/by from "B_<bx>_<by>" (supports negatives: B_-12_7)
  
  #' Parse a single block ID into grid coordinates
  #'
  #' Internal helper used by `make_block_polygons_grid()` to convert
  #' encoded IDs (e.g., 'x###_y###') into numeric indices.
  #'
  #' @param id Single block identifier string.
  #' @return A named list or vector with parsed components (x/y indices).
  parse_one <- function(id) {
    s <- sub("^B_", "", id)
    parts <- strsplit(s, "_", fixed = TRUE)[[1]]
    if (length(parts) != 2) stop("Unexpected block_id format: ", id)
    
    bx <- suppressWarnings(as.integer(parts[1]))
    by <- suppressWarnings(as.integer(parts[2]))
    if (is.na(bx) || is.na(by)) stop("Could not parse bx/by from block_id: ", id)
    
    c(bx = bx, by = by)
  }
  
  parsed <- t(vapply(ids, parse_one, numeric(2)))
  bx <- parsed[, "bx"]
  by <- parsed[, "by"]
  
  x0 <- bx * block_size_km
  x1 <- (bx + 1) * block_size_km
  y0 <- by * block_size_km
  y1 <- (by + 1) * block_size_km
  
  polys <- lapply(seq_along(ids), function(i) {
    coords <- matrix(
      c(
        x0[i], y0[i],
        x1[i], y0[i],
        x1[i], y1[i],
        x0[i], y1[i],
        x0[i], y0[i]
      ),
      ncol = 2,
      byrow = TRUE
    )
    sf::st_polygon(list(coords))
  })
  
  sf::st_as_sf(
    data.frame(block_id = ids, stringsAsFactors = FALSE),
    geometry = sf::st_sfc(polys, crs = crs)
  )
}