# ==============================================================================
# Species Distribution Model Validation: Observed vs. Predicted
# ==============================================================================
# Compares survey-based observations (honeycomb) against model predictions
# (raster) at multiple spatial scales (provincial, BCR, regional, local).
#
# Public API (in calling order):
#
#   GOAL 1 – Spatial grid
#     make_hex_grid()                 Create hexagon layer over study area
#
#   GOAL 2 – Per-hexagon summaries
#     extract_hex_predictions()       Extract raster predictions into each hex
#     summarize_surveys_by_hex()      Summarise observed counts per hex
#     summarize_hex()                 Convenience wrapper: runs both of the above
#                                     and merges the result
#
#   GOAL 3 – Region-wide statistics
#     compute_region_stats()          Correlation, P/A agreement, detection rate
#                                     and counts across all surveyed hexagons
#
#   GOAL 4 – Visualisation + output wrapper
#     assess_region()                 Master wrapper: accepts a pre-built hex
#                                     grid (or builds one), then summaries ->
#                                     stats -> three-panel figure + tidy output
#
#   Supporting plot helpers (used internally or standalone)
#     plot_raster_gg()
#     plot_hex_pred_obs()
#     plot_honeycomb()
#     process_water()                 Pre-process water polygons for mapping
#     resolve_water_layer()           Pick fine/medium/coarse water detail
# ==============================================================================

# Dependencies: sf, terra, ggplot2, ggforce, ggspatial, patchwork, dplyr, scales.
# Functions use explicit namespace calls where practical so this file can be
# sourced directly or moved into an R package.

# ==============================================================================
# Raster preprocessing
# ==============================================================================

#' Set raster values below a lower bound to NA
#'
#' A thin convenience wrapper for zeroing out low-abundance predictions before
#' passing a raster to `assess_region()` or the individual plot helpers.
#' Values strictly less than `lower_bound` are replaced with NA; values greater
#' than or equal to `lower_bound` are left unchanged.
#'
#' @param rast        A `SpatRaster`.
#' @param lower_bound Numeric scalar. Pixel values below this threshold are set
#'   to NA. Pass `0` to convert only negative values (e.g. to clean up a raster
#'   that should be non-negative but contains small floating-point negatives).
#' @return A `SpatRaster` with the same extent, resolution, and CRS as `rast`,
#'   with sub-threshold pixels set to NA.
#' @export
apply_raster_lower_bound <- function(rast, lower_bound) {
  
  stopifnot(inherits(rast, "SpatRaster"))
  stopifnot(is.numeric(lower_bound), length(lower_bound) == 1, is.finite(lower_bound))
  
  rast[rast < lower_bound] <- NA
  rast
}


# ==============================================================================
# GOAL 1 – Hexagon grid
# ==============================================================================

#' Create a hexagon grid over a study area
#'
#' Tiles the study area with regular hexagons so that the total number of
#' hexagons that *intersect* the boundary is approximately `n_hexagons`.
#' Internally uses EPSG:3978 (Canada Albers) when the input is geographic, then
#' reprojects back to the original CRS before returning.
#'
#' @param study_area  An `sf` object defining the region of interest.
#' @param n_hexagons  Target number of hexagons inside the study area (default 200).
#' @return An `sf` object with columns `hex_id` (integer) and `geometry`.
#' @export
make_hex_grid <- function(study_area, n_hexagons = 200) {
  
  stopifnot(inherits(study_area, "sf"))
  
  old_s2 <- sf::sf_use_s2()
  on.exit(sf::sf_use_s2(old_s2), add = TRUE)
  sf::sf_use_s2(FALSE)
  
  original_crs <- sf::st_crs(study_area)
  
  if (sf::st_is_longlat(study_area)) {
    study_area_work <- sf::st_transform(study_area, 3978)
  } else {
    study_area_work <- study_area
  }
  
  study_area_work <- study_area_work |>
    sf::st_make_valid() |>
    sf::st_collection_extract("POLYGON") |>
    sf::st_buffer(0)
  
  study_area_union <- sf::st_union(study_area_work)
  
  target_area  <- as.numeric(sf::st_area(study_area_union)) / n_hexagons
  hex_cellsize <- sqrt(target_area / (sqrt(3) / 2))
  
  hex_grid <- sf::st_make_grid(
    study_area_union,
    cellsize = hex_cellsize,
    square   = FALSE,
    what     = "polygons"
  )
  
  hex_sf <- sf::st_sf(
    hex_id   = seq_along(hex_grid),
    geometry = hex_grid,
    crs      = sf::st_crs(study_area_work)
  )
  
  keep   <- lengths(sf::st_intersects(hex_sf, study_area_union)) > 0
  hex_sf <- hex_sf[keep, ]
  hex_sf$hex_id <- seq_len(nrow(hex_sf))
  
  if (!is.na(original_crs)) {
    hex_sf <- sf::st_transform(hex_sf, original_crs)
  }
  
  hex_sf
}


# ==============================================================================
# GOAL 2 – Per-hexagon summaries
# ==============================================================================

#' Summarise point-count / ARU surveys by hexagon
#'
#' Spatially joins survey points to hexagons and computes, for each hexagon:
#' the number of surveys, the mean count-per-effort, and a binary
#' presence/absence flag.
#'
#' @param dat      An `sf` object with a `count_per_effort` numeric column.
#' @param hex_grid An `sf` hexagon grid (output of `make_hex_grid()`).
#' @return The full hex grid with added columns:
#'   \describe{
#'     \item{`n_surveys`}{Number of survey visits (0 for unsurveyed hexagons).}
#'     \item{`mean_count_per_effort`}{Mean count per effort (NA for unsurveyed).}
#'     \item{`obs_detected`}{Logical: TRUE if any survey detected the species.}
#'   }
#' @export
summarize_surveys_by_hex <- function(dat, hex_grid) {
  
  stopifnot(inherits(dat, "sf"))
  stopifnot(inherits(hex_grid, "sf"))
  stopifnot("count_per_effort" %in% names(dat))
  
  if (sf::st_crs(dat) != sf::st_crs(hex_grid)) {
    dat <- sf::st_transform(dat, sf::st_crs(hex_grid))
  }
  
  if (!"hex_id" %in% names(hex_grid)) {
    hex_grid$hex_id <- seq_len(nrow(hex_grid))
  }
  
  dat_hex <- sf::st_join(
    dat,
    hex_grid["hex_id"],
    join = sf::st_within,
    left = FALSE
  )
  
  hex_obs <- dat_hex |>
    sf::st_drop_geometry() |>
    dplyr::group_by(hex_id) |>
    dplyr::summarise(
      n_surveys             = dplyr::n(),
      mean_count_per_effort = mean(count_per_effort, na.rm = TRUE),
      obs_detected          = any(count_per_effort > 0, na.rm = TRUE),
      .groups = "drop"
    )
  
  hex_grid |>
    dplyr::left_join(hex_obs, by = "hex_id") |>
    dplyr::mutate(
      n_surveys             = dplyr::coalesce(n_surveys, 0L),
      mean_count_per_effort = dplyr::if_else(n_surveys > 0,
                                             mean_count_per_effort, NA_real_),
      obs_detected          = dplyr::if_else(n_surveys > 0,
                                             dplyr::coalesce(obs_detected, FALSE),
                                             NA)
    )
}


#' Extract raster model predictions into each hexagon
#'
#' Computes the mean predicted abundance within each hexagon, treating absent
#' pixels (0 or NA in the input raster) as zero so that partial occupancy is
#' correctly reflected in the hex mean. For example, if 50 % of pixels in a
#' hexagon are absent (0/NA) and the remaining 50 % have a predicted value of
#' 1, the returned `pred_mean` will be 0.5.
#'
#' The upper-quantile cap (`rast_max_q`) is applied to non-zero pixels before
#' extraction so that extreme outliers do not distort hex means, but absent
#' pixels are never excluded from the denominator.
#'
#' @param hex_grid     An `sf` hexagon grid (output of `make_hex_grid()`).
#' @param rast         A `SpatRaster` of model predictions. Absent cells should
#'   be encoded as 0 or NA; both are treated as zero during averaging.
#' @param rast_max_q   Upper quantile used to cap non-zero pixel values before
#'   extraction. Default `0.99`.
#' @return The hex grid with an added column `pred_mean` (mean predicted
#'   relative abundance within the hexagon, inclusive of absent pixels; NA only
#'   where the hexagon falls entirely outside the raster extent).
#' @export
extract_hex_predictions <- function(hex_grid,
                                    rast,
                                    rast_max_q = 0.99) {
  
  stopifnot(inherits(hex_grid, "sf"))
  stopifnot(inherits(rast, "SpatRaster"))
  
  # Align CRS
  hex_proj  <- sf::st_transform(hex_grid, terra::crs(rast))
  rast_crop <- terra::crop(rast, terra::vect(hex_proj))
  
  # Cap extreme non-zero values (outlier suppression only; does not affect
  # which pixels count as absent)
  r_max <- as.numeric(stats::quantile(
    terra::values(rast_crop)[terra::values(rast_crop) > 0],
    rast_max_q, na.rm = TRUE
  ))
  if (is.finite(r_max) && r_max > 0) {
    rast_crop <- terra::clamp(rast_crop, upper = r_max, values = TRUE)
  }
  
  # Convert NA to 0 so absent pixels are included in the denominator.
  # We work on a copy to avoid modifying the caller's raster object.
  rast_for_extract <- terra::subst(rast_crop, from = NA, to = 0)
  
  # Extract: na.rm = FALSE is intentional — there are no NAs left, and we
  # want every pixel (including the zeros) in the mean.
  pred_ex  <- terra::extract(rast_for_extract, terra::vect(hex_proj),
                             fun = mean, na.rm = FALSE)
  pred_col <- names(pred_ex)[2]
  
  hex_grid |>
    dplyr::mutate(
      # Hexagons with no raster coverage at all are left as NA
      pred_mean = dplyr::if_else(
        is.finite(pred_ex[[pred_col]]),
        pred_ex[[pred_col]],
        NA_real_
      )
    )
}


#' Combine observed-survey and model-prediction summaries at the hexagon level
#'
#' A convenience wrapper that runs `summarize_surveys_by_hex()` and
#' `extract_hex_predictions()` in sequence and merges both sets of columns onto
#' the hex grid. Also adds a `pred_detected` binary column (TRUE when
#' `pred_mean` strictly exceeds `pred_presence_threshold`).
#'
#' Absent pixels (0 or NA) in the raster are included in the hex mean
#' denominator, so `pred_mean` reflects true spatial average abundance across
#' the whole hexagon area (see `extract_hex_predictions()` for details).
#'
#' @param dat                    An `sf` object with a `count_per_effort` column.
#' @param hex_grid               Output of `make_hex_grid()`.
#' @param rast                   A `SpatRaster` of model predictions.
#' @param rast_max_q             Upper quantile cap passed to
#'   `extract_hex_predictions()`. Default `0.99`.
#' @param pred_presence_threshold `pred_mean` values strictly above this are
#'   classified as model-predicted present. Default `0` (any positive mean
#'   abundance counts as presence).
#' @return An `sf` object (one row per hexagon) with columns:
#'   `hex_id`, `n_surveys`, `mean_count_per_effort`, `obs_detected`,
#'   `pred_mean`, `pred_detected`, `geometry`.
#' @export
summarize_hex <- function(dat,
                          hex_grid,
                          rast,
                          rast_max_q              = 0.99,
                          pred_presence_threshold = 0) {
  
  hex_obs  <- summarize_surveys_by_hex(dat, hex_grid)
  
  hex_full <- extract_hex_predictions(hex_obs, rast,
                                      rast_max_q = rast_max_q)
  
  hex_full |>
    dplyr::mutate(
      pred_detected = dplyr::if_else(
        is.finite(pred_mean) & pred_mean > pred_presence_threshold,
        TRUE, FALSE
      )
    )
}


# ==============================================================================
# GOAL 3 – Region-wide statistics
# ==============================================================================

#' Compute region-wide model-validation statistics
#'
#' Takes a hex-level summary (output of `summarize_hex()`) and computes
#' summary statistics across all hexagons that were actually surveyed.
#'
#' **Statistics returned**
#' | Name | Description |
#' |------|-------------|
#' | `n_hex_surveyed` | Hexagons with >=1 survey |
#' | `n_hex_detected` | Hexagons where species was detected |
#' | `prop_hex_detected` | Detection rate (detected / surveyed) |
#' | `cor_obs_pred` | Pearson r between mean observed count and mean predicted abundance (unweighted) |
#' | `cor_obs_pred_wtd` | Pearson r weighted by number of surveys per hexagon |
#' | `pa_agreement` | Proportion of surveyed hexagons where observed P/A matches predicted P/A |
#' | `pa_sensitivity` | True-positive rate: model predicted present when observed present |
#' | `pa_specificity` | True-negative rate: model predicted absent when observed absent |
#' | `n_surveys_total` | Total individual survey visits across region |
#' | `mean_surveys_per_hex` | Average surveys per surveyed hexagon |
#' | `prop_det_in_pred_absent` | Of hexagons where species was detected, proportion where mean model prediction was at or below `absence_threshold` |
#' | `absence_threshold` | The threshold value used to define model-predicted absence (echoed from input) |
#' | `prop_absent_in_pred_high` | Of hexagons where model predicted high abundance (>= `high_abundance_threshold`), proportion where species was not detected |
#' | `high_abundance_threshold` | The threshold used to define high predicted abundance; either user-supplied or the 75th percentile of non-zero `pred_mean` values (echoed from input) |
#'
#' @param hex_summary       Output of `summarize_hex()` (or equivalent `sf` with
#'   columns `n_surveys`, `mean_count_per_effort`, `obs_detected`,
#'   `pred_mean`, `pred_detected`).
#' @param absence_threshold      Numeric scalar. Hexagons with `pred_mean` at
#'   or below this value (or with no raster coverage, i.e. `pred_mean` is NA)
#'   are considered model-predicted absent when computing
#'   `prop_det_in_pred_absent`. Default `0`.
#' @param high_abundance_threshold  Numeric scalar or `NULL`. Hexagons with
#'   `pred_mean` at or above this value are considered high-abundance when
#'   computing `prop_absent_in_pred_high`. When `NULL` (default), the threshold
#'   is set automatically to the 75th percentile of non-zero `pred_mean` values
#'   across surveyed hexagons.
#' @return A one-row `data.frame` of region-wide statistics.
#' @export
compute_region_stats <- function(hex_summary,
                                 absence_threshold       = 1e-4,
                                 high_abundance_threshold = NULL) {
  
  required_cols <- c("n_surveys", "mean_count_per_effort",
                     "obs_detected", "pred_mean", "pred_detected")
  missing <- setdiff(required_cols, names(hex_summary))
  if (length(missing) > 0) {
    stop("hex_summary is missing columns: ", paste(missing, collapse = ", "),
         "\nRun summarize_hex() first.")
  }
  
  surveyed <- hex_summary |>
    sf::st_drop_geometry() |>
    dplyr::filter(n_surveys > 0)
  
  n_hex_surveyed  <- nrow(surveyed)
  n_hex_detected  <- sum(surveyed$obs_detected, na.rm = TRUE)
  prop_detected   <- if (n_hex_surveyed > 0) n_hex_detected / n_hex_surveyed
  else NA_real_
  n_surveys_total <- sum(surveyed$n_surveys, na.rm = TRUE)
  mean_surveys_ph <- if (n_hex_surveyed > 0) n_surveys_total / n_hex_surveyed
  else NA_real_
  
  # -- Correlation (requires both obs and pred to be finite) ------------------
  cor_dat <- surveyed |>
    dplyr::filter(is.finite(mean_count_per_effort),
                  is.finite(pred_mean))
  
  cor_obs_pred <- if (nrow(cor_dat) >= 3) {
    stats::cor(cor_dat$pred_mean, cor_dat$mean_count_per_effort,
               use = "complete.obs", method = "pearson")
  } else NA_real_
  
  cor_obs_pred_wtd <- if (nrow(cor_dat) >= 3 &&
                          sum(cor_dat$n_surveys, na.rm = TRUE) > 0) {
    x  <- cor_dat$pred_mean
    y  <- cor_dat$mean_count_per_effort
    w  <- cor_dat$n_surveys
    mx <- sum(w * x) / sum(w)
    my <- sum(w * y) / sum(w)
    cov_xy <- sum(w * (x - mx) * (y - my)) / sum(w)
    var_x  <- sum(w * (x - mx)^2) / sum(w)
    var_y  <- sum(w * (y - my)^2) / sum(w)
    cov_xy / sqrt(var_x * var_y)
  } else NA_real_
  
  # -- Presence / absence agreement -------------------------------------------
  pa_dat <- surveyed |>
    dplyr::filter(!is.na(obs_detected), !is.na(pred_detected))
  
  n_pa   <- nrow(pa_dat)
  tp     <- sum( pa_dat$obs_detected &  pa_dat$pred_detected, na.rm = TRUE)
  tn     <- sum(!pa_dat$obs_detected & !pa_dat$pred_detected, na.rm = TRUE)
  fp     <- sum(!pa_dat$obs_detected &  pa_dat$pred_detected, na.rm = TRUE)
  fn     <- sum( pa_dat$obs_detected & !pa_dat$pred_detected, na.rm = TRUE)
  
  pa_agreement   <- if (n_pa > 0) (tp + tn) / n_pa       else NA_real_
  pa_sensitivity <- if ((tp + fn) > 0) tp / (tp + fn)    else NA_real_
  pa_specificity <- if ((tn + fp) > 0) tn / (tn + fp)    else NA_real_
  
  # -- Detections in model-predicted absent hexagons --------------------------
  # Numerator:   surveyed hexagons where the species was detected AND
  #              pred_mean is at or below absence_threshold (or NA, meaning
  #              no raster coverage -- the model implicitly predicts absence).
  # Denominator: all surveyed hexagons where the species was detected.
  # Interpretation: a high value flags a systematic model blind spot --
  # the species is frequently found in places the model considers absent.
  pred_absent <- is.na(surveyed$pred_mean) |
    surveyed$pred_mean <= absence_threshold
  
  n_det_in_pred_absent    <- sum(surveyed$obs_detected & pred_absent,
                                 na.rm = TRUE)
  prop_det_in_pred_absent <- if (n_hex_detected > 0)
    n_det_in_pred_absent / n_hex_detected
  else NA_real_
  
  # -- Non-detections in model-predicted high-abundance hexagons --------------
  # Threshold: default to the 75th percentile of non-zero pred_mean values
  # across surveyed hexagons; user may supply an explicit override.
  # Numerator:   surveyed hexagons where the species was NOT detected AND
  #              pred_mean >= high_abundance_threshold.
  # Denominator: all surveyed hexagons where pred_mean >= high_abundance_threshold.
  # Interpretation: a high value means the model is confidently predicting
  # presence in places where observers are not finding the species.
  nonzero_preds <- surveyed$pred_mean[is.finite(surveyed$pred_mean) &
                                        surveyed$pred_mean > 0]
  
  if (is.null(high_abundance_threshold)) {
    high_abundance_threshold <- if (length(nonzero_preds) >= 4)
      as.numeric(stats::quantile(nonzero_preds, probs = 0.75, na.rm = TRUE))
    else NA_real_
  }
  
  if (is.finite(high_abundance_threshold)) {
    pred_high      <- is.finite(surveyed$pred_mean) &
      surveyed$pred_mean >= high_abundance_threshold
    n_pred_high    <- sum(pred_high, na.rm = TRUE)
    n_absent_in_pred_high    <- sum(!surveyed$obs_detected & pred_high,
                                    na.rm = TRUE)
    prop_absent_in_pred_high <- if (n_pred_high > 0)
      n_absent_in_pred_high / n_pred_high
    else NA_real_
  } else {
    high_abundance_threshold  <- NA_real_
    prop_absent_in_pred_high  <- NA_real_
  }
  
  data.frame(
    n_hex_surveyed           = n_hex_surveyed,
    n_hex_detected           = n_hex_detected,
    prop_hex_detected        = prop_detected,
    n_surveys_total          = n_surveys_total,
    mean_surveys_per_hex     = mean_surveys_ph,
    cor_obs_pred             = cor_obs_pred,
    cor_obs_pred_wtd         = cor_obs_pred_wtd,
    pa_agreement             = pa_agreement,
    pa_sensitivity           = pa_sensitivity,
    pa_specificity           = pa_specificity,
    prop_det_in_pred_absent  = prop_det_in_pred_absent,
    absence_threshold        = absence_threshold,
    prop_absent_in_pred_high = prop_absent_in_pred_high,
    high_abundance_threshold = high_abundance_threshold
  )
}

# ==============================================================================
# GOAL 4 – Master wrapper
# ==============================================================================

#' Assess model predictions against survey observations at any spatial scale
#'
#' Orchestrates the full workflow:
#'   1. Accept a pre-built hexagon grid (or build one if `n_hexagons` is
#'      supplied and `hex_grid` is NULL).
#'   2. Summarise survey observations and raster predictions per hexagon
#'      (`summarize_hex`).
#'   3. Compute region-wide validation statistics (`compute_region_stats`).
#'   4. Build a three-panel patchwork figure:
#'      - Panel A: raw raster of model predictions (`plot_raster_gg`).
#'      - Panel B: hex-filled by mean predicted abundance, circles by observed
#'        mean count (`plot_hex_pred_obs`).
#'      - Panel C: honeycomb effort/detection map (`plot_honeycomb`).
#'
#' When looping over many species for the same region, build the hex grid once
#' with `make_hex_grid()` and pass it via `hex_grid` to avoid rebuilding it on
#' every iteration.
#'
#' @param region          `sf` polygon defining the study area extent. Used for
#'   plotting and, when `hex_grid` is NULL, for building the grid.
#' @param sp_dat          `sf` data frame with a `count_per_effort` column.
#' @param rast            `SpatRaster` of model predictions. Absent cells should
#'   be 0 or NA; both are treated as zero when computing hex means, and rendered
#'   as transparent in all map panels.
#' @param hex_grid        Pre-built hexagon grid (`sf` output of
#'   `make_hex_grid()`). When supplied, `n_hexagons` is ignored and no grid is
#'   built internally. Pass `NULL` (default) to build the grid from `region`
#'   and `n_hexagons`.
#' @param n_hexagons      Target number of hexagons. Only used when `hex_grid`
#'   is NULL. Default `1000`.
#' @param water           Optional `processed_water` object (from
#'                        `process_water()`). Pass `NULL` to omit water.
#' @param rast_max_q      Quantile (0-1) used to cap the colour scale upper end
#'   and clip outlier pixel values before extraction. Default `0.99`.
#' @param transform       Colour scale transform: `"identity"`, `"sqrt"`,
#'                        `"log"`, etc.
#' @param pred_presence_threshold  `pred_mean` values strictly above this are
#'   classified as model-predicted present in `hex_summary` and used for P/A
#'   statistics. Default `0`.
#' @param title           Optional character string added as a patchwork
#'                        annotation title above the three panels.
#' @param model_source    Optional character string identifying the model source
#'                        (e.g. `"BAM"`, `"Atlas"`). Displayed in the subtitle
#'                        when either `model_source` or `data_source` is
#'                        non-NULL.
#' @param data_source     Optional character string identifying the survey data
#'                        source (e.g. `"BAM"`, `"Atlas"`). Displayed in the
#'                        subtitle alongside `model_source`.
#'
#' @return A named list:
#'   \describe{
#'     \item{`plot_combined`}{A `patchwork` ggplot (three panels).}
#'     \item{`hex_summary`}{`sf` data frame with one row per hexagon
#'       containing columns from `summarize_hex()`.}
#'     \item{`region_stats`}{One-row `data.frame` from `compute_region_stats()`.}
#'   }
#' @export
assess_region <- function(region,
                          sp_dat,
                          rast,
                          hex_grid                = NULL,
                          n_hexagons              = 1000,
                          water                   = NULL,
                          rast_max_q              = 0.99,
                          transform               = "identity",
                          pred_presence_threshold = 0,
                          title                   = NULL,
                          model_source            = NULL,
                          data_source             = NULL) {
  
  stopifnot(inherits(region, "sf"))
  stopifnot(inherits(sp_dat, "sf"))
  stopifnot(inherits(rast, "SpatRaster"))
  if (!is.null(hex_grid)) stopifnot(inherits(hex_grid, "sf"))
  
  region <- sf::st_transform(region, sf::st_crs(sp_dat))
  sp_dat <- sf::st_filter(sp_dat, region, .predicate = sf::st_intersects)
  
  # -- Trim raster to region ---------------------------------------------------
  
  region_vect <- region |>
    sf::st_transform(terra::crs(rast)) |>
    terra::vect()
  
  rast_cropped <- rast |>
    terra::crop(region_vect) |>
    terra::mask(region_vect)
  
  # -- GOAL 1: Hexagon grid ---------------------------------------------------
  # Use the supplied grid if provided; build one otherwise.
  if (is.null(hex_grid)) {
    hex_grid <- make_hex_grid(region, n_hexagons = n_hexagons)
  }
  
  # -- GOAL 2: Per-hexagon summaries -----------------------------------------
  hex_summary <- summarize_hex(
    dat                     = sp_dat,
    hex_grid                = hex_grid,
    rast                    = rast_cropped,
    rast_max_q              = rast_max_q,
    pred_presence_threshold = pred_presence_threshold
  )
  
  # -- GOAL 3: Region-wide statistics ----------------------------------------
  region_stats <- compute_region_stats(hex_summary)
  
  # -- GOAL 4a: Panel A – raw raster -----------------------------------------
  p_rast <- plot_raster_gg(
    rast       = rast_cropped,
    study_area = region,
    water      = NULL,        # Water is hidden on the raster
    rast_max_q = rast_max_q,
    transform  = transform
  )
  
  # -- GOAL 4b: Panel B – hex predicted vs. observed -------------------------
  p_hex <- plot_hex_pred_obs(
    hex_summary = hex_summary,
    rast        = rast_cropped,
    study_area  = region,
    water       = water,
    rast_max_q  = rast_max_q,
    transform   = transform
  )
  
  # -- GOAL 4c: Panel C – honeycomb effort/detection -------------------------
  p_honey <- plot_honeycomb(
    hex_summary = hex_summary,
    study_area  = region,
    water       = water
  )
  
  # -- Compose figure --------------------------------------------------------
  margin_theme <- ggplot2::theme(plot.margin = ggplot2::margin(5, 10, 5, 10))
  
  combined <-
    (p_rast           + margin_theme) +
    (p_hex$plot       + margin_theme) +
    (p_honey          + margin_theme) +
    patchwork::plot_layout(ncol = 3)
  
  # Build annotation: title and optional subtitle from model/data source
  if (!is.null(title) || !is.null(model_source) || !is.null(data_source)) {
    
    subtitle <- NULL
    if (!is.null(model_source) || !is.null(data_source)) {
      parts <- c(
        if (!is.null(model_source)) paste("Model:", model_source),
        if (!is.null(data_source))  paste("Data:",  data_source)
      )
      subtitle <- paste(parts, collapse = "   |   ")
    }
    
    combined <- combined +
      patchwork::plot_annotation(
        title    = title,
        subtitle = subtitle,
        theme    = ggplot2::theme(
          plot.title    = ggplot2::element_text(size = 18, face = "bold",
                                                margin = ggplot2::margin(b = 2)),
          plot.subtitle = ggplot2::element_text(size = 12, colour = "grey30",
                                                margin = ggplot2::margin(b = 6))
        )
      )
  }
  
  list(
    plot_combined = combined,
    hex_summary   = hex_summary,
    region_stats  = region_stats
  )
}


# ==============================================================================
# Supporting functions – water pre-processing
# ==============================================================================

#' Pre-process water polygons for plotting
#'
#' Reprojects, crops, and creates three simplification levels so that the
#' appropriate detail is selected automatically at plot time.
#'
#' @param study_area             An `sf` object.
#' @param water                  An `sf` object of water polygons.
#' @param target_crs             CRS string for internal projection (LCC km).
#' @param area_breaks_km2        Two-element numeric: area thresholds (km²) that
#'   separate fine / medium / coarse simplification levels.
#' @param simplify_tolerances_km Named vector of `dTolerance` values (km).
#' @return A list of class `"processed_water"`.
#' @export
process_water <- function(study_area,
                          water,
                          target_crs        = "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=49 +lon_0=-95 +datum=NAD83 +units=km +no_defs",
                          area_breaks_km2   = c(25e3, 250e3),
                          simplify_tolerances_km = c(
                            fine   = 0.05,
                            medium = 1,
                            coarse = 10
                          )) {
  
  stopifnot(inherits(study_area, "sf"))
  stopifnot(inherits(water, "sf"))
  
  study_area <- sf::st_transform(study_area, target_crs)
  water      <- sf::st_transform(water, target_crs)
  
  study_area_union <- sf::st_union(study_area)
  
  water_sub <- water |>
    sf::st_crop(sf::st_bbox(study_area)) |>
    sf::st_filter(study_area_union, .predicate = sf::st_intersects)
  
  water_layers <- list(
    fine   = sf::st_simplify(water_sub, dTolerance = simplify_tolerances_km[["fine"]],
                             preserveTopology = TRUE),
    medium = sf::st_simplify(water_sub, dTolerance = simplify_tolerances_km[["medium"]],
                             preserveTopology = TRUE),
    coarse = sf::st_simplify(water_sub, dTolerance = simplify_tolerances_km[["coarse"]],
                             preserveTopology = TRUE)
  )
  
  structure(
    list(
      water_layers    = water_layers,
      study_area      = study_area,
      region_area_km2 = as.numeric(sf::st_area(study_area_union)),
      area_breaks_km2 = area_breaks_km2,
      crs             = sf::st_crs(study_area)
    ),
    class = c("processed_water", "list")
  )
}


#' Resolve the water layer to plot given the current region size
#'
#' @param water       A `processed_water` object.
#' @param study_area  An `sf` object (already in the plot CRS).
#' @return An `sf` polygon layer cropped and filtered to `study_area`.
resolve_water_layer <- function(water, study_area) {
  
  region_area_km2 <- as.numeric(sf::st_area(sf::st_union(study_area)))
  
  detail <- dplyr::case_when(
    region_area_km2 < water$area_breaks_km2[1] ~ "fine",
    region_area_km2 < water$area_breaks_km2[2] ~ "medium",
    TRUE                                        ~ "coarse"
  )
  
  water$water_layers[[detail]] |>
    sf::st_crop(study_area) |>
    sf::st_filter(study_area, .predicate = sf::st_intersects) |>
    sf::st_intersection(study_area)
}


# ==============================================================================
# Panel plot helpers
# ==============================================================================

# -- Internal: effort-class breaks and alpha scale ----------------------------

#' Build effort-class breaks and alpha scale for the honeycomb plot
#' @keywords internal
make_effort_classes <- function(max_surveys, alpha_min, alpha_max) {
  
  max_surveys <- ceiling(max_surveys)
  if (!is.finite(max_surveys) || max_surveys < 1) max_surveys <- 1
  
  n_positive_classes <- min(4, max_surveys)
  
  positive_breaks <- unique(ceiling(
    (seq(0, 1, length.out = n_positive_classes + 1)^2) * max_surveys
  ))
  positive_breaks <- positive_breaks[positive_breaks > 0]
  positive_breaks[length(positive_breaks)] <- max_surveys
  
  starts <- c(1, head(positive_breaks, -1) + 1)
  ends   <- positive_breaks
  valid  <- starts <= ends
  starts <- starts[valid]
  ends   <- ends[valid]
  
  positive_labels <- ifelse(
    starts == ends,
    as.character(starts),
    paste0(starts, "-", ends)
  )
  
  effort_levels <- c("0", positive_labels, paste0(">", max_surveys))
  
  alpha_values <- c(
    "0" = 0,
    setNames(
      seq(alpha_min, alpha_max, length.out = length(positive_labels) + 1),
      c(positive_labels, paste0(">", max_surveys))
    )
  )
  
  list(
    max_surveys     = max_surveys,
    starts          = starts,
    ends            = ends,
    positive_labels = positive_labels,
    effort_levels   = effort_levels,
    alpha_values    = alpha_values
  )
}


# -- Panel A: raw raster -------------------------------------------------------

#' Plot model-predicted relative abundance as a raster layer
#'
#' Pixels with a value of 0 or NA are rendered as transparent (absent).
#' All non-zero pixels are mapped onto the colour ramp, with the upper end
#' capped at the `rast_max_q` quantile to prevent outliers from compressing
#' the scale.
#'
#' @param rast         A `SpatRaster`. Absent cells should be 0 or NA.
#' @param study_area   `sf` polygon defining the region boundary.
#' @param water        Optional `processed_water` object.
#' @param rast_max_q   Quantile (0-1) used to cap the upper end of the colour
#'   scale. Default `0.99`.
#' @param palette      Character vector of colours (auto if NULL).
#' @param water_fill   Fill colour for water polygons.
#' @param transform    Scale transform: `"identity"`, `"sqrt"`, `"log"`, etc.
#' @return A `ggplot` object.
#' @export
plot_raster_gg <- function(rast,
                           study_area,
                           water      = NULL,
                           rast_max_q = 0.99,
                           palette    = NULL,
                           water_fill = "#b8dceb",
                           transform  = "identity") {
  
  stopifnot(inherits(rast, "SpatRaster"))
  stopifnot(inherits(study_area, "sf"))
  
  if (is.null(palette)) {
    palette <- grDevices::colorRampPalette(c(
      "#FBF7E2", "#FCF8D0", "#EEF7C2", "#CEF2B0",
      "#94E5A0", "#51C987", "#18A065", "#008C59",
      "#007F53", "#006344"
    ))(100)
  }
  
  water_to_plot <- NULL
  if (!is.null(water)) {
    stopifnot(inherits(water, "processed_water"))
    plot_crs   <- water$crs
    study_area <- sf::st_transform(study_area, plot_crs)
    rast       <- terra::project(rast, plot_crs$wkt)
    water_to_plot <- resolve_water_layer(water, study_area)
  } else {
    plot_crs <- sf::st_crs(study_area)
    rast     <- terra::project(rast, plot_crs$wkt)
  }
  
  rast <- terra::crop(rast, terra::vect(study_area))
  rast <- terra::mask(rast, terra::vect(study_area))
  
  # Blank out absent pixels (0 or NA to NA so they render as transparent)
  rast[rast <= 0] <- NA
  
  # Cap the upper end of the scale; compute from non-zero values only
  r_max <- as.numeric(stats::quantile(
    terra::values(rast), rast_max_q, na.rm = TRUE
  ))
  if (!is.finite(r_max) || r_max <= 0) r_max <- 1
  rast <- terra::clamp(rast, upper = r_max, values = TRUE)
  
  rast_df           <- as.data.frame(rast, xy = TRUE, na.rm = FALSE)
  names(rast_df)[3] <- "value"
  
  p <- ggplot2::ggplot() +
    ggplot2::geom_raster(data = rast_df,
                         ggplot2::aes(x = x, y = y, fill = value)) +
    ggplot2::geom_sf(data = study_area, fill = NA, colour = "gray30")
  
  if (!is.null(water_to_plot)) {
    p <- p + ggplot2::geom_sf(data = water_to_plot,
                              fill = water_fill, colour = water_fill,
                              linewidth = 0.1)
  }
  
  p +
    ggspatial::annotation_scale(location = "br", width_hint = 0.25) +
    ggspatial::annotation_north_arrow(
      location    = "tr",
      which_north = "true",
      style       = ggspatial::north_arrow_fancy_orienteering
    ) +
    ggplot2::scale_fill_gradientn(
      colours  = palette,
      limits   = c(NA, r_max),   # lower bound auto-detected from non-zero data
      oob      = scales::squish,
      na.value = "transparent",  # absent pixels show panel background through
      trans    = transform,
      name     = "Relative abundance"
    ) +
    ggplot2::coord_sf() +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position      = c(0.03, 0.03),
      legend.justification = c(0, 0),
      legend.background    = ggplot2::element_rect(fill = "white", colour = NA),
      legend.margin        = ggplot2::margin(4, 4, 4, 4),
      legend.box.margin    = ggplot2::margin(6, 6, 6, 6),
      plot.background      = ggplot2::element_rect(fill = "white", colour = "black",
                                                   linewidth = 0.5),
      panel.background     = ggplot2::element_rect(fill = "white", colour = NA),
      plot.title    = ggplot2::element_text(size = 14, face = "bold", hjust = 0),
      plot.subtitle = ggplot2::element_text(size = 12, hjust = 0,
                                            margin = ggplot2::margin(b = 6))
    ) +
    ggplot2::labs(title = "Model predictions")
}


# -- Panel B: hex predicted vs. observed ---------------------------------------

#' Plot observed mean counts over hexagon-averaged model predictions
#'
#' Fills each surveyed hexagon by mean predicted relative abundance and overlays
#' circles scaled by observed mean count per effort. Hexagons with a
#' `pred_mean` of 0 or NA are rendered as transparent (absent); all others are
#' mapped onto the colour ramp. The raster is used only to derive the upper cap
#' for the shared colour scale (`rast_max_q`).
#'
#' @param hex_summary          Output of `summarize_hex()`. Must contain
#'   `pred_mean`, `n_surveys`, and `mean_count_per_effort`.
#' @param rast                 A `SpatRaster` used to compute the colour scale
#'   upper limit. Not re-extracted here.
#' @param study_area           `sf` polygon.
#' @param water                Optional `processed_water` object.
#' @param rast_max_q           Quantile (0-1) used to cap the colour scale
#'   upper end. Default `0.99`.
#' @param max_count_per_effort Upper bound for observed circle scaling (auto).
#' @param palette              Fill palette; default matches `plot_raster_gg()`.
#' @param circle_fill          Fill colour for observed-count circles.
#' @param water_fill           Fill colour for water polygons.
#' @param transform            Fill scale transform.
#' @return A list: `plot` (ggplot), `hex_summary` (sf), `corr_unweighted` (dbl).
#' @export
plot_hex_pred_obs <- function(hex_summary,
                              rast,
                              study_area,
                              water                = NULL,
                              rast_max_q           = 0.99,
                              max_count_per_effort = NULL,
                              palette              = NULL,
                              circle_fill          = "black",
                              water_fill           = "#b8dceb",
                              transform            = "identity") {
  
  stopifnot(inherits(hex_summary, "sf"))
  stopifnot(inherits(rast, "SpatRaster"))
  stopifnot(inherits(study_area, "sf"))
  stopifnot("mean_count_per_effort" %in% names(hex_summary))
  stopifnot("n_surveys" %in% names(hex_summary))
  
  # Ensure pred_mean is present (re-extract if missing).
  if (!"pred_mean" %in% names(hex_summary)) {
    hex_summary <- extract_hex_predictions(hex_summary, rast,
                                           rast_max_q = rast_max_q)
  }
  
  if (is.null(palette)) {
    palette <- grDevices::colorRampPalette(c(
      "#FBF7E2", "#FCF8D0", "#EEF7C2", "#CEF2B0",
      "#94E5A0", "#51C987", "#18A065", "#008C59",
      "#007F53", "#006344"
    ))(100)
  }
  
  water_to_plot <- NULL
  
  if (!is.null(water)) {
    stopifnot(inherits(water, "processed_water"))
    plot_crs    <- water$crs
    study_area  <- sf::st_transform(study_area, plot_crs)
    hex_summary <- sf::st_transform(hex_summary, plot_crs)
    rast        <- terra::project(rast, plot_crs$wkt)
    water_to_plot <- resolve_water_layer(water, study_area)
  } else {
    plot_crs    <- sf::st_crs(hex_summary)
    study_area  <- sf::st_transform(study_area, plot_crs)
    hex_summary <- sf::st_transform(hex_summary, plot_crs)
    rast        <- terra::project(rast, plot_crs$wkt)
  }
  
  # Derive colour scale upper cap from the raster (display copy only)
  rast_display <- terra::crop(rast, terra::vect(study_area))
  rast_display <- terra::mask(rast_display, terra::vect(study_area))
  rast_display[rast_display <= 0] <- NA
  r_max <- as.numeric(stats::quantile(
    terra::values(rast_display), rast_max_q, na.rm = TRUE
  ))
  if (!is.finite(r_max) || r_max <= 0) r_max <- 1
  
  # Hexagons with pred_mean of 0 or NA map to NA to transparent
  hex_summary <- hex_summary |>
    dplyr::mutate(
      pred_mean_display = dplyr::if_else(
        is.finite(pred_mean) & pred_mean > 0,
        pred_mean,
        NA_real_
      )
    )
  
  if (is.null(max_count_per_effort)) {
    detected <- hex_summary$mean_count_per_effort[
      !is.na(hex_summary$mean_count_per_effort) &
        hex_summary$mean_count_per_effort > 0
    ]
    max_count_per_effort <- if (length(detected) == 0) 1 else
      as.numeric(stats::quantile(detected, probs = 0.99, na.rm = TRUE))
  }
  if (!is.finite(max_count_per_effort) || max_count_per_effort <= 0) {
    max_count_per_effort <- 1
  }
  
  hex_area       <- as.numeric(sf::st_area(hex_summary))
  inner_diameter <- 2 * sqrt(hex_area / (2 * sqrt(3)))
  
  circle_sf <- hex_summary |>
    dplyr::mutate(
      mean_count_per_effort = dplyr::coalesce(mean_count_per_effort, 0),
      n_surveys             = dplyr::coalesce(n_surveys, 0L),
      count_scaled  = pmin(mean_count_per_effort, max_count_per_effort) /
        max_count_per_effort,
      circle_diameter = dplyr::if_else(
        mean_count_per_effort > 0 & n_surveys > 0,
        0.1 * inner_diameter + sqrt(count_scaled) * (0.6 - 0.1) * inner_diameter,
        0
      ),
      circle_radius = circle_diameter / 2
    ) |>
    sf::st_point_on_surface()
  
  xy        <- sf::st_coordinates(circle_sf)
  circle_df <- circle_sf |>
    sf::st_drop_geometry() |>
    dplyr::mutate(x = xy[, 1], y = xy[, 2])
  
  corr_dat <- hex_summary |>
    sf::st_drop_geometry() |>
    dplyr::filter(n_surveys > 0,
                  is.finite(mean_count_per_effort),
                  is.finite(pred_mean))
  
  corr_unweighted <- if (nrow(corr_dat) >= 3) {
    stats::cor(corr_dat$pred_mean, corr_dat$mean_count_per_effort,
               use = "complete.obs", method = "pearson")
  } else NA_real_
  
  corr_label <- paste0(
    "Cor(Obs, Pred) = ",
    ifelse(is.finite(corr_unweighted), round(corr_unweighted, 2), "NA")
  )
  
  bbox    <- sf::st_bbox(study_area)
  corr_df <- data.frame(
    x     = bbox["xmax"] - 0.0001 * (bbox["xmax"] - bbox["xmin"]),
    y     = bbox["ymax"] - 0.0001 * (bbox["ymax"] - bbox["ymin"]),
    label = corr_label
  )
  
  p <- ggplot2::ggplot()
  
  if (!is.null(water_to_plot)) {
    p <- p + ggplot2::geom_sf(data = water_to_plot,
                              fill = water_fill, colour = water_fill,
                              linewidth = 0.1)
  }
  
  p <- p +
    ggplot2::geom_sf(
      data   = subset(hex_summary, n_surveys > 0),
      ggplot2::aes(fill = pred_mean_display),
      colour = "gray50", linewidth = 0.1, alpha = 1
    ) +
    ggplot2::geom_sf(data = study_area, fill = NA, colour = "gray30") +
    ggforce::geom_circle(
      data = dplyr::filter(circle_df, circle_radius > 0),
      ggplot2::aes(x0 = x, y0 = y, r = circle_radius),
      alpha = 1, fill = circle_fill, colour = NA
    ) +
    ggspatial::annotation_scale(location = "br", width_hint = 0.25) +
    ggplot2::scale_fill_gradientn(
      colours  = palette,
      limits   = c(NA, r_max),   # lower bound auto-detected from non-zero data
      oob      = scales::squish,
      na.value = "transparent",  # absent hexagons show panel background through
      trans    = transform,
      name     = "Mean predicted\nrelative abundance"
    ) +
    ggplot2::geom_label(
      data = corr_df,
      ggplot2::aes(x = x, y = y, label = label),
      hjust = 1, vjust = 1,
      fill = "white", colour = "black",
      label.size = 0, size = 4, family = "mono"
    ) +
    ggplot2::coord_sf() +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position      = c(0.03, 0.03),
      legend.justification = c(0, 0),
      legend.background    = ggplot2::element_rect(fill = "white", colour = NA),
      legend.margin        = ggplot2::margin(4, 4, 4, 4),
      legend.box.margin    = ggplot2::margin(6, 6, 6, 6),
      plot.background      = ggplot2::element_rect(fill = "white", colour = "black",
                                                   linewidth = 0.5),
      panel.background     = ggplot2::element_rect(fill = "white", colour = NA),
      plot.title    = ggplot2::element_text(size = 14, face = "bold", hjust = 0),
      plot.subtitle = ggplot2::element_text(size = 12, hjust = 0,
                                            margin = ggplot2::margin(b = 6))
    ) +
    ggplot2::labs(
      title    = "Predicted vs observed",
      subtitle = paste(
        "Hex fill = mean prediction within hexagon",
        "Circle size = observed mean count in hexagon",
        sep = "\n"
      )
    )
  
  list(plot = p, hex_summary = hex_summary, corr_unweighted = corr_unweighted)
}


# -- Panel C: honeycomb effort/detection map -----------------------------------

#' Plot a honeycomb survey-effort / detection map
#'
#' Fills hexagons by survey effort (number of visits) and overlays circles
#' scaled by observed mean count per effort.
#'
#' @param hex_summary          Output of `summarize_hex()` (or
#'   `summarize_surveys_by_hex()`).
#' @param study_area           `sf` polygon.
#' @param water                Optional `processed_water` object.
#' @param max_surveys          Upper bound for effort colour classes (auto).
#' @param max_count_per_effort Upper bound for circle scaling (auto).
#' @param alpha_min,alpha_max  Range of hex fill transparency.
#' @param hex_fill             Fill colour for hexagons.
#' @param circle_fill          Fill colour for detection circles.
#' @param water_fill           Fill colour for water polygons.
#' @return A `ggplot` object.
#' @export
plot_honeycomb <- function(hex_summary,
                           study_area,
                           water                = NULL,
                           max_surveys          = NULL,
                           max_count_per_effort = NULL,
                           alpha_min            = 0.15,
                           alpha_max            = 1,
                           hex_fill             = "#B38F47", 
                           circle_fill          = "black",
                           water_fill           = "#b8dceb") {
  
  stopifnot(inherits(hex_summary, "sf"))
  stopifnot(inherits(study_area, "sf"))
  stopifnot("n_surveys" %in% names(hex_summary))
  stopifnot("mean_count_per_effort" %in% names(hex_summary))
  
  water_to_plot <- NULL
  
  if (!is.null(water)) {
    stopifnot(inherits(water, "processed_water"))
    plot_crs    <- water$crs
    study_area  <- sf::st_transform(study_area, plot_crs)
    hex_summary <- sf::st_transform(hex_summary, plot_crs)
    water_to_plot <- resolve_water_layer(water, study_area)
  } else {
    plot_crs   <- sf::st_crs(hex_summary)
    study_area <- sf::st_transform(study_area, plot_crs)
  }
  
  if (is.null(max_surveys)) {
    surveyed    <- hex_summary$n_surveys[
      !is.na(hex_summary$n_surveys) & hex_summary$n_surveys >= 1
    ]
    max_surveys <- if (length(surveyed) == 0) 1 else
      as.numeric(stats::quantile(surveyed, probs = 0.8, na.rm = TRUE))
  }
  max_surveys <- ceiling(max_surveys)
  if (!is.finite(max_surveys) || max_surveys < 1) max_surveys <- 1
  
  if (is.null(max_count_per_effort)) {
    detected <- hex_summary$mean_count_per_effort[
      !is.na(hex_summary$mean_count_per_effort) &
        hex_summary$mean_count_per_effort > 0
    ]
    max_count_per_effort <- if (length(detected) == 0) 1 else
      as.numeric(stats::quantile(detected, probs = 0.99, na.rm = TRUE))
  }
  if (!is.finite(max_count_per_effort) || max_count_per_effort <= 0) {
    max_count_per_effort <- 1
  }
  
  effort_info   <- make_effort_classes(max_surveys, alpha_min, alpha_max)
  max_surveys   <- effort_info$max_surveys
  effort_levels <- effort_info$effort_levels
  alpha_values  <- effort_info$alpha_values
  
  hex_area       <- as.numeric(sf::st_area(hex_summary))
  inner_diameter <- 2 * sqrt(hex_area / (2 * sqrt(3)))
  
  hex_summary_plot <- hex_summary |>
    dplyr::mutate(
      n_surveys    = dplyr::coalesce(n_surveys, 0),
      effort_class = dplyr::case_when(
        n_surveys <= 0          ~ "0",
        n_surveys > max_surveys ~ paste0(">", max_surveys),
        TRUE ~ effort_info$positive_labels[
          findInterval(n_surveys,
                       vec = c(0, effort_info$ends),
                       rightmost.closed = TRUE)
        ]
      ),
      effort_class = factor(effort_class, levels = effort_levels)
    )
  
  circle_sf <- hex_summary_plot |>
    dplyr::mutate(
      mean_count_per_effort = dplyr::coalesce(mean_count_per_effort, 0),
      count_scaled  = pmin(mean_count_per_effort, max_count_per_effort) /
        max_count_per_effort,
      circle_diameter = dplyr::if_else(
        mean_count_per_effort > 0 & n_surveys > 0,
        0.1 * inner_diameter + sqrt(count_scaled) * (0.6 - 0.1) * inner_diameter,
        0
      ),
      circle_radius = circle_diameter / 2
    ) |>
    sf::st_point_on_surface()
  
  xy        <- sf::st_coordinates(circle_sf)
  circle_df <- circle_sf |>
    sf::st_drop_geometry() |>
    dplyr::mutate(x = xy[, 1], y = xy[, 2])
  
  p <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = study_area, fill = NA, colour = "gray80")
  
  if (!is.null(water_to_plot)) {
    p <- p + ggplot2::geom_sf(data = water_to_plot,
                              fill = water_fill, colour = water_fill,
                              linewidth = 0.1)
  }
  
  p +
    ggplot2::geom_sf(
      data = hex_summary_plot,
      ggplot2::aes(alpha = effort_class),
      fill = hex_fill, colour = NA
    ) +
    ggforce::geom_circle(
      data = dplyr::filter(circle_df, circle_radius > 0),
      ggplot2::aes(x0 = x, y0 = y, r = circle_radius),
      alpha = 1, fill = circle_fill, colour = NA
    ) +
    ggplot2::geom_point(
      data = dplyr::filter(circle_df, circle_radius > 0),
      ggplot2::aes(x = x, y = y, size = mean_count_per_effort),
      shape = 21, fill = circle_fill, colour = NA, alpha = 0
    ) +
    ggspatial::annotation_scale(location = "br", width_hint = 0.25) +
    ggplot2::scale_alpha_manual(
      name   = "Survey effort",
      values = alpha_values,
      drop   = FALSE,
      guide  = ggplot2::guide_legend(override.aes = list(fill = hex_fill,
                                                         colour = NA))
    ) +
    ggplot2::scale_size_continuous(name = "Mean count per effort", guide = "none") +
    ggplot2::guides(
      alpha = ggplot2::guide_legend(order = 1,
                                    override.aes = list(fill = hex_fill,
                                                        colour = NA))
    ) +
    ggplot2::coord_sf() +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position      = c(0.03, 0.03),
      legend.justification = c(0, 0),
      legend.background    = ggplot2::element_rect(fill = "white", colour = NA),
      legend.margin        = ggplot2::margin(4, 4, 4, 4),
      legend.box.margin    = ggplot2::margin(6, 6, 6, 6),
      plot.background      = ggplot2::element_rect(fill = "white", colour = "black",
                                                   linewidth = 0.5),
      panel.background     = ggplot2::element_rect(fill = "white", colour = NA),
      plot.title    = ggplot2::element_text(size = 14, face = "bold", hjust = 0),
      plot.subtitle = ggplot2::element_text(size = 12, hjust = 0,
                                            margin = ggplot2::margin(b = 6))
    ) +
    ggplot2::labs(
      title    = "Relative survey effort",
      subtitle = paste(
        "Hex fill = Number of surveys in hexagon",
        "Circle size = observed mean count in hexagon",
        sep = "\n"
      )
    )
}
