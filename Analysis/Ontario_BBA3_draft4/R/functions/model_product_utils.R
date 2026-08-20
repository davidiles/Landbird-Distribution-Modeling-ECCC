# ============================================================
# model_product_utils.R
#
# Purpose
#   Utility functions for generating model products from fitted
#   OBBA abundance/change models. Functions cover rasterization,
#   relative-abundance maps, hex-level change summaries, and
#   model-assessment panels comparing predictions with observations.
#
# Documentation style
#   Function comments are intentionally plain comments rather than
#   ROxygen comments, because this script is sourced directly rather
#   than built as an R package.
#
# Contents
#   Rasterization & scaling
#     rasterize_sf, resolve_positive_max, resolve_rast_max,
#     compute_assessment_scale_limits, prepare_relative_abundance_rasters
#   Maps
#     make_map, plot_raster_gg
#   Hex-level change summaries & maps
#     summarize_hex_draw_change, compute_hex_polygon_weights,
#     summarize_polygon_hex_draw_change,
#     classify_min_supported_change, get_min_supported_change_scale,
#     make_min_supported_change_legend, make_hex_abs_change_map
#   Region assessment core
#     assess_region, make_hex_grid, summarize_hex, summarize_surveys_by_hex,
#     extract_hex_predictions, compute_region_stats
#   Panel plot helpers
#     plot_honeycomb, make_effort_classes, make_change_comparison_plot,
#     make_paired_summary_table
#   Per-region / per-page assembly (used by 08_generate_model_products.R)
#     make_region_geom, load_species_surveys, require_bcr_key,
#     build_period_assessments, build_period_honeycombs, save_page
# ============================================================


# ------------------------------------------------------------
# Extract a pixel_id -> Biol_Region lookup from one prediction grid.
#
# grid        dat$grid_OBBA2 or dat$grid_OBBA3 (an sf carrying pixel_id, hex_id
#             and Biol_Region, all written by 06).
# region_col  Region-label column (default "Biol_Region").
#
# Returns a plain data.frame(pixel_id, hex_id, Biol_Region) -- no geometry, no
# spatial ops. Fails loudly if the grid predates 06's Biol_Region tagging.
# ------------------------------------------------------------
grid_region_lookup <- function(grid, region_col = "Biol_Region") {
  
  stopifnot(inherits(grid, "sf") || is.data.frame(grid))
  
  needed  <- c("pixel_id", "hex_id", region_col)
  missing <- setdiff(needed, names(grid))
  if (length(missing) > 0) {
    stop("grid_region_lookup(): grid is missing ", paste(missing, collapse = ", "),
         ". Re-run 06_filter_and_finalize_surveys.R so the saved grids carry ",
         "pixel_id, hex_id and ", region_col, ".", call. = FALSE)
  }
  
  g <- sf::st_drop_geometry(grid)
  data.frame(
    pixel_id    = g$pixel_id,
    hex_id      = g$hex_id,
    Biol_Region = as.character(g[[region_col]]),
    stringsAsFactors = FALSE
  )
}


# ------------------------------------------------------------
# Zero every prediction product for regions the species has no safe dates in.
#
# preds                  One prediction record (readRDS of a 07 output file).
# lookup_obba2           grid_region_lookup(dat$grid_OBBA2).
# lookup_obba3           grid_region_lookup(dat$grid_OBBA3).
# safe_regions           Regions WITH safe dates. NULL -> taken from
#                        preds$sp_safe_dates$Biol_Region.
# hex_no_safe_threshold  A hexagon is zeroed when this FRACTION (or more) of its
#                        pixels fall in no-safe-date regions. 0.5 = majority
#                        (default). Set 0 to zero a hex if ANY pixel is no-safe;
#                        set ~1 to zero only hexes lying entirely in no-safe
#                        regions.
# trim_unclassified      TRUE also zeros pixels/hexes with region = NA. FALSE
#                        (default) keeps them. Moot in practice, since 06 tags
#                        the grid with nearest_fallback = TRUE (no NA regions);
#                        retained as a guard against a stale grid.
#
# Returns preds with OBBA2/OBBA3/abs_change pixel summaries and hex_draws (both
# plain and _Corrected_for_Water) zeroed in place.
# ------------------------------------------------------------
trim_predictions_to_safe_dates <- function(preds,
                                           lookup_obba2,
                                           lookup_obba3,
                                           safe_regions          = NULL,
                                           hex_no_safe_threshold = 0.5,
                                           trim_unclassified     = FALSE,
                                           verbose               = TRUE) {
  
  # ---- Resolve the species' safe-date regions ----
  if (is.null(safe_regions)) {
    if (is.null(preds$sp_safe_dates) ||
        !"Biol_Region" %in% names(preds$sp_safe_dates)) {
      stop("trim_predictions_to_safe_dates(): safe_regions not supplied and ",
           "preds$sp_safe_dates$Biol_Region is unavailable.", call. = FALSE)
    }
    safe_regions <- preds$sp_safe_dates$Biol_Region
  }
  safe_regions <- unique(as.character(safe_regions))
  safe_regions <- safe_regions[!is.na(safe_regions)]
  
  # No safe-date table at all -> leave predictions untouched.
  if (length(safe_regions) == 0) {
    if (verbose) message("  safe-date trim: no safe regions listed; ",
                         "predictions left unchanged.")
    return(preds)
  }
  
  # ---- Guard against a Biol_Region label mismatch ----
  # Grid and safe-date labels both descend from Biol_Regions_ON, so this should
  # never fire; but if NONE of the safe regions match the grid labels, a naive
  # %in% would flag EVERY pixel as no-safe and erase the whole surface. Stop.
  grid_regions <- unique(c(lookup_obba2$Biol_Region, lookup_obba3$Biol_Region))
  grid_regions <- grid_regions[!is.na(grid_regions)]
  if (length(intersect(safe_regions, grid_regions)) == 0) {
    stop("trim_predictions_to_safe_dates(): none of the species' safe-date ",
         "regions match the Biol_Region labels on the prediction grid.\n",
         "  safe_regions: ", paste(safe_regions, collapse = ", "), "\n",
         "  grid regions: ", paste(sort(grid_regions), collapse = ", "), "\n",
         "This is almost certainly a Biol_Region label mismatch; refusing to ",
         "zero the whole surface.", call. = FALSE)
  }
  
  # ---- Per-pixel keep flag for a summary (matched by pixel_id, order-safe) ----
  pixel_keep <- function(pixel_ids, lookup) {
    region       <- lookup$Biol_Region[match(pixel_ids, lookup$pixel_id)]
    keep         <- region %in% safe_regions
    unclassified <- is.na(region)
    if (any(unclassified)) keep[unclassified] <- !trim_unclassified
    keep
  }
  
  # ---- Zero a pixel-summary data.frame in place ----
  zero_pixel_summary <- function(df, prefix, lookup) {
    if (is.null(df)) return(df)
    if (!"pixel_id" %in% names(df)) {
      stop("trim: summary for prefix '", prefix, "' has no pixel_id column.",
           call. = FALSE)
    }
    drop <- !pixel_keep(df$pixel_id, lookup)
    if (!any(drop)) return(df)
    
    val_cols  <- grep(paste0("^", prefix, "_"), names(df), value = TRUE)
    cv_cols   <- grep(paste0("^", prefix, "_cv"), names(df), value = TRUE)
    zero_cols <- setdiff(val_cols, cv_cols)
    
    for (cc in zero_cols) df[[cc]][drop] <- 0            # abundance / change -> 0
    for (cc in cv_cols)   df[[cc]][drop] <- NA_real_     # CV undefined at 0
    df
  }
  
  # ---- Per-hex no-safe fraction, from the OBBA2 grid (hex_draws' provenance) ----
  reg   <- lookup_obba2$Biol_Region
  nsafe <- !(reg %in% safe_regions)
  nsafe[is.na(reg)] <- isTRUE(trim_unclassified)
  hex_frac <- tapply(nsafe, lookup_obba2$hex_id, function(z) mean(z, na.rm = TRUE))
  
  # ---- Zero the draw vectors of no-safe hexes in place (length preserved) ----
  zero_hex_draws <- function(hd) {
    if (is.null(hd) || nrow(hd) == 0) return(hd)
    f <- hex_frac[as.character(hd$hex_id)]
    f[is.na(f)] <- 0                       # hex absent from lookup -> keep
    drop <- f >= hex_no_safe_threshold
    if (!any(drop)) return(hd)
    
    for (j in which(drop)) {
      hd$mu_OBBA2[[j]] <- numeric(length(hd$mu_OBBA2[[j]]))
      hd$mu_OBBA3[[j]] <- numeric(length(hd$mu_OBBA3[[j]]))
      if ("abs_change" %in% names(hd)) {
        hd$abs_change[[j]] <- numeric(length(hd$abs_change[[j]]))
      }
    }
    hd
  }
  
  # ---- Apply to every product variant ----
  # OBBA2 and abs_change carry the OBBA2 grid's pixel_id; OBBA3 the OBBA3 grid's.
  preds$OBBA2      <- zero_pixel_summary(preds$OBBA2,      "OBBA2",      lookup_obba2)
  preds$OBBA3      <- zero_pixel_summary(preds$OBBA3,      "OBBA3",      lookup_obba3)
  preds$abs_change <- zero_pixel_summary(preds$abs_change, "abs_change", lookup_obba2)
  
  preds$OBBA2_Corrected_for_Water <-
    zero_pixel_summary(preds$OBBA2_Corrected_for_Water, "OBBA2", lookup_obba2)
  preds$OBBA3_Corrected_for_Water <-
    zero_pixel_summary(preds$OBBA3_Corrected_for_Water, "OBBA3", lookup_obba3)
  preds$abs_change_Corrected_for_Water <-
    zero_pixel_summary(preds$abs_change_Corrected_for_Water, "abs_change", lookup_obba2)
  
  preds$hex_draws <- zero_hex_draws(preds$hex_draws)
  preds$hex_draws_Corrected_for_Water <-
    zero_hex_draws(preds$hex_draws_Corrected_for_Water)
  
  if (verbose) {
    n_reg_trim <- length(setdiff(grid_regions, safe_regions))
    n_hex_trim <- sum(hex_frac >= hex_no_safe_threshold, na.rm = TRUE)
    message("  safe-date trim: ", n_reg_trim, " region(s) without safe dates -> ",
            n_hex_trim, " hexagon(s) zeroed (threshold ",
            hex_no_safe_threshold, ").")
  }
  
  preds
}


# ------------------------------------------------------------
# Shared constants
# ------------------------------------------------------------

# Per-panel plot margin reused across composed assessment figures.
margin_theme <- ggplot2::theme(plot.margin = ggplot2::margin(5, 10, 5, 10))

# Rasterize fields from an sf grid to a terra SpatRaster.
#
# Purpose:
#   Converts one or more numeric columns in an sf prediction grid into raster
#   layers at a specified resolution. This is used to create map-ready raster
#   products from 1-km prediction grids. Metadata tags are attached to the
#   returned raster.
#
# Arguments:
#   grid_sf  - sf object containing polygon or point/grid-cell geometries.
#   field    - character vector of column names to rasterize.
#   res      - raster resolution in the units of grid_sf.
#   metadata - named character vector of metadata tags to attach.
#
# Returns:
#   A terra SpatRaster with one layer per requested field.
rasterize_sf <- function(grid_sf, field, res,metadata) {
  stopifnot(inherits(grid_sf, "sf"))
  stopifnot(field %in% names(grid_sf))
  
  v <- terra::vect(grid_sf)
  r_template <- terra::rast(v, res = res)
  r <- terra::rasterize(v, r_template, field = field, fun = mean, na.rm = TRUE)
  names(r) <- field
  metags(r) <- metadata
  
  r
}


# ------------------------------------------------------------
# Shared plotting-scale helpers
# ------------------------------------------------------------

# Resolve a positive upper plotting/extraction limit from either an explicit
# maximum or a quantile of positive finite values.
#
# Purpose:
#   Several model-assessment panels need to use identical scaling across atlas
#   periods. This helper centralizes the logic so a user can either supply an
#   explicit maximum, such as 1.5 expected birds per 5-min point count, or ask
#   for a quantile-based maximum, such as the 0.99 quantile.
#
# Arguments:
#   x              Numeric vector of values used to derive the limit.
#   max_value      Optional explicit positive maximum. If supplied, this takes
#                  precedence over max_q.
#   max_q          Quantile used when max_value is NULL.
#   default        Fallback value when x has no positive finite values.
#   argument_name  Name used in error messages.
#
# Returns:
#   A single positive numeric value.
resolve_positive_max <- function(x,
                                 max_value = NULL,
                                 max_q = 0.99,
                                 default = 1,
                                 argument_name = "max_value") {
  
  if (!is.null(max_value)) {
    if (!is.numeric(max_value) || length(max_value) != 1 ||
        !is.finite(max_value) || max_value <= 0) {
      stop(argument_name, " must be NULL or a single positive finite number.")
    }
    return(as.numeric(max_value))
  }
  
  if (!is.numeric(max_q) || length(max_q) != 1 ||
      !is.finite(max_q) || max_q <= 0 || max_q > 1) {
    stop("max_q must be a single finite number in (0, 1].")
  }
  
  vals <- x[is.finite(x) & x > 0]
  if (length(vals) == 0) return(default)
  
  out <- as.numeric(stats::quantile(vals, probs = max_q, na.rm = TRUE))
  if (!is.finite(out) || out <= 0) default else out
}

# Resolve an upper relative-abundance limit from a terra SpatRaster.
#
# Argument: rast       SpatRaster of model predictions.
# Argument: rast_max   Optional explicit upper limit for prediction scaling.
# Argument: rast_max_q Quantile used when rast_max is NULL.
# Returns: A single positive numeric value.
resolve_rast_max <- function(rast,
                             rast_max = NULL,
                             rast_max_q = 0.99) {
  stopifnot(inherits(rast, "SpatRaster"))
  
  resolve_positive_max(
    x             = terra::values(rast),
    max_value     = rast_max,
    max_q         = rast_max_q,
    default       = 1,
    argument_name = "rast_max"
  )
}

# Compute shared scale limits for model-assessment panels.
#
# Purpose:
#   When producing assessment panels for multiple atlas periods, this helper can
#   be run after creating `hex_summary` objects. The returned values can then
#   be supplied to
#   each call to `assess_region()` so all panels use identical scaling.
#
# Arguments:
#   hex_summaries        One or more `sf` objects returned by `summarize_hex()`
#                        or `assess_region()$hex_summary`.
#   rasts                Optional list of SpatRaster objects used to calculate a
#                        shared prediction limit directly from raster values.
#   rast_max             Optional explicit shared prediction maximum.
#   rast_max_q           Quantile used when rast_max is NULL.
#   max_count_per_effort Optional explicit shared observed-count maximum.
#   count_max_q          Quantile used when max_count_per_effort is NULL.
#   max_surveys          Optional explicit shared survey-effort maximum.
#   max_surveys_q        Quantile used when max_surveys is NULL.
#
# Returns:
#   A named list with `rast_max`, `max_count_per_effort`, and `max_surveys`.
compute_assessment_scale_limits <- function(...,
                                            hex_summaries = list(...),
                                            rasts = NULL,
                                            rast_max = NULL,
                                            rast_max_q = 0.99,
                                            max_count_per_effort = NULL,
                                            count_max_q = 0.99,
                                            max_surveys = NULL,
                                            max_surveys_q = 0.80) {
  
  if (length(hex_summaries) == 1 && is.list(hex_summaries[[1]]) &&
      !inherits(hex_summaries[[1]], "sf")) {
    hex_summaries <- hex_summaries[[1]]
  }
  
  if (length(hex_summaries) == 0) {
    stop("At least one hex_summary object must be supplied.")
  }
  
  if (!all(vapply(hex_summaries, inherits, logical(1), what = "sf"))) {
    stop("All hex_summaries must be sf objects.")
  }
  
  count_values <- unlist(lapply(hex_summaries, function(x) {
    if (!"mean_count_per_effort" %in% names(x)) {
      stop("Each hex_summary must contain mean_count_per_effort.")
    }
    x$mean_count_per_effort
  }), use.names = FALSE)
  
  survey_values <- unlist(lapply(hex_summaries, function(x) {
    if (!"n_surveys" %in% names(x)) {
      stop("Each hex_summary must contain n_surveys.")
    }
    x$n_surveys
  }), use.names = FALSE)
  
  if (!is.null(rasts)) {
    if (inherits(rasts, "SpatRaster")) rasts <- list(rasts)
    if (!all(vapply(rasts, inherits, logical(1), what = "SpatRaster"))) {
      stop("rasts must be NULL, a SpatRaster, or a list of SpatRaster objects.")
    }
    rast_values <- unlist(lapply(rasts, terra::values), use.names = FALSE)
  } else {
    rast_values <- unlist(lapply(hex_summaries, function(x) {
      if (!"pred_mean" %in% names(x)) {
        stop("Each hex_summary must contain pred_mean when rasts is NULL.")
      }
      x$pred_mean
    }), use.names = FALSE)
  }
  
  list(
    rast_max = resolve_positive_max(
      x             = rast_values,
      max_value     = rast_max,
      max_q         = rast_max_q,
      default       = 1,
      argument_name = "rast_max"
    ),
    max_count_per_effort = resolve_positive_max(
      x             = count_values,
      max_value     = max_count_per_effort,
      max_q         = count_max_q,
      default       = 1,
      argument_name = "max_count_per_effort"
    ),
    max_surveys = ceiling(resolve_positive_max(
      x             = survey_values,
      max_value     = max_surveys,
      max_q         = max_surveys_q,
      default       = 1,
      argument_name = "max_surveys"
    ))
  )
}

# ------------------------------------------------------------
# Relative abundance maps
# ------------------------------------------------------------

# # Prepare multiple relative-abundance rasters for side-by-side plotting.
# #
# # Purpose:
# #   Applies a shared absence threshold and a shared upper colour-scale limit
# #   across one or more rasters, so Atlas 2 and Atlas 3 maps use comparable
# #   legends. Low values are treated as absent, and high values are clamped at
# #   a common quantile-based maximum.
# #
# # Arguments:
# #   ...                 - named terra SpatRaster objects.
# #   rast_absent_limit   - values at or below this threshold are set to NA.
# #   rast_max_quantile   - quantile used to estimate the upper plotting limit.
# #
# # Returns:
# #   A list containing processed rasters, shared z limits, legend breaks, and
# #   original per-raster upper quantiles.
# prepare_relative_abundance_rasters <- function(...,
#                                                rast_absent_limit = 1/250,
#                                                rast_max_quantile = 0.99,
#                                                rast_max = NULL) {
#   
#   rasts <- list(...)
#   
#   if (length(rasts) == 0) {
#     stop("At least one raster must be supplied.")
#   }
#   
#   if (!all(vapply(rasts, inherits, logical(1), what = "SpatRaster"))) {
#     stop("All inputs passed through `...` must be terra SpatRaster objects.")
#   }
#   
#   process_raster <- function(rast) {
#     
#     r <- rast
#     
#     # Set low values to NA
#     r[r <= rast_absent_limit] <- NA
#     
#     # Calculate upper limit after low values are removed. If rast_max is
#     # supplied, it is used directly; otherwise use rast_max_quantile.
#     zmax <- resolve_rast_max(
#       rast       = r,
#       rast_max   = rast_max,
#       rast_max_q = rast_max_quantile
#     )
#     
#     list(
#       rast = r,
#       zmax = zmax
#     )
#   }
#   
#   processed <- lapply(rasts, process_raster)
#   
#   zmax_values <- vapply(processed, `[[`, numeric(1), "zmax")
#   
#   # Shared upper limit across all rasters. Use the maximum atlas-specific
#   # limit so all maps are displayed on the same absolute scale without
#   # saturating the higher-abundance atlas by default.
#   zmax_shared <- max(zmax_values, na.rm = TRUE)
#   
#   if (!is.finite(zmax_shared) || zmax_shared <= 0) {
#     zmax_shared <- 1
#   }
#   
#   # Clamp all rasters to the same upper limit
#   rasts_clamped <- lapply(processed, function(x) {
#     terra::clamp(
#       x$rast,
#       upper = zmax_shared,
#       values = TRUE
#     )
#   })
#   
#   names(rasts_clamped) <- names(rasts)
#   
#   zlim <- c(rast_absent_limit, zmax_shared)
#   
#   # Fixed legend breaks shared by all maps
#   zbreaks <- seq(zlim[1], zlim[2], length.out = 5)
#   
#   list(
#     rasters = rasts_clamped,
#     zlim = zlim,
#     zbreaks = zbreaks,
#     zmax_original = zmax_values
#   )
# }


prepare_relative_abundance_rasters <- function(...,
                                               absent_method = c("fixed", "cumulative", "hybrid"),

                                               # --- cumulative-population floor (absent_method = "cumulative") ---
                                               # White out the low-density pixels that TOGETHER hold only the
                                               # bottom (1 - coverage) of the predicted population, so the coloured
                                               # area retains at least `coverage` of the birds. Population per pixel
                                               # is density * cell area, so the floor is projection-robust rather
                                               # than assuming an equal-area grid.
                                               coverage         = 0.99,
                                               coverage_basis   = c("each", "pooled"),
                                               min_absent_limit = 0,     # absolute lower guard on the floor (0 = none)

                                               # --- fixed / hybrid floor (legacy) ---
                                               rast_absent_limit    = 1 / 250,
                                               absent_quantile      = 0.01,
                                               relative_to_quantile = 0.99,
                                               relative_fraction    = 0.01,

                                               # --- shared upper colour limit ---
                                               rast_max_quantile = 0.99,
                                               rast_max = NULL) {
  
  rasts <- list(...)
  absent_method  <- match.arg(absent_method)
  coverage_basis <- match.arg(coverage_basis)
  
  if (length(rasts) == 0) {
    stop("At least one raster must be supplied.")
  }
  
  if (!all(vapply(rasts, inherits, logical(1), what = "SpatRaster"))) {
    stop("All inputs passed through `...` must be terra SpatRaster objects.")
  }
  
  # ------------------------------------------------------------
  # Resolve absence floor
  # ------------------------------------------------------------
  # cumulative (recommended for abundance surfaces):
  #   Rank pixels by density and white out the smallest-density pixels whose
  #   summed POPULATION (density * cell area) is at most (1 - coverage) of the
  #   total. The coloured area then holds >= coverage of the predicted birds --
  #   the statement made in the map caption. Area-weighting makes the floor
  #   correct on any projection, not only an equal-area grid.
  #
  #   coverage_basis sets how the single, cross-period floor is formed:
  #     "each"   -> resolve the floor separately for each period and take the
  #                 MOST INCLUSIVE (smallest), so >= coverage is retained in BOTH
  #                 maps. A period that is nearly empty therefore pulls the shared
  #                 floor down -- deliberate: it never erases a declining or
  #                 colonising period's signal. Default.
  #     "pooled" -> one floor on the concatenated pixels of all periods (>= coverage
  #                 of the COMBINED population). Tighter maps, but the smaller
  #                 period can drop below coverage if its density is lower.
  #
  # fixed:
  #   Uses `rast_absent_limit` directly (a fixed detection-referenced density).
  #   This is the default, so callers that threshold non-abundance surfaces
  #   (e.g. probability-of-observation) are unaffected by the cumulative option.
  #
  # hybrid (legacy):
  #   max(rast_absent_limit, low prediction quantile, small fraction of a high
  #   prediction quantile). Retained for backward compatibility.
  # ------------------------------------------------------------
  
  # Extract finite, positive (density, cell-area) pairs from one raster.
  dens_area <- function(r) {
    v <- terra::values(r, mat = FALSE)
    a <- terra::values(terra::cellSize(r, unit = "km"), mat = FALSE)
    keep <- is.finite(v) & is.finite(a) & v > 0
    list(v = v[keep], a = a[keep])
  }
  
  # Population-weighted floor for a set of (density, area) pairs: the density at
  # or below which the cumulative population is <= (1 - coverage) of the total.
  cumulative_floor <- function(v, a) {
    if (length(v) == 0) return(0)
    o    <- order(v)                     # ascending density
    v    <- v[o]; a <- a[o]
    frac <- cumsum(v * a) / sum(v * a)   # cumulative population fraction from the bottom
    below <- which(frac <= (1 - coverage))  # pixels shed without eating into coverage
    if (length(below) == 0) 0 else v[max(below)]
  }
  
  resolve_absent_limit <- function(rasts) {
    
    if (absent_method == "fixed") {
      return(rast_absent_limit)
    }
    
    if (absent_method == "cumulative") {
      stopifnot(coverage > 0, coverage < 1)
      
      if (coverage_basis == "each") {
        floors <- vapply(rasts, function(r) {
          da <- dens_area(r)
          cumulative_floor(da$v, da$a)
        }, numeric(1))
        floor_val <- min(floors, na.rm = TRUE)      # most inclusive -> >= coverage in every period
      } else {
        pooled <- lapply(rasts, dens_area)
        v_all  <- unlist(lapply(pooled, `[[`, "v"))
        a_all  <- unlist(lapply(pooled, `[[`, "a"))
        floor_val <- cumulative_floor(v_all, a_all)
      }
      
      return(max(min_absent_limit, floor_val, na.rm = TRUE))
    }
    
    # hybrid (legacy)
    vals <- unlist(lapply(rasts, function(r) {
      terra::values(r, mat = FALSE, na.rm = TRUE)
    }))
    vals <- vals[is.finite(vals)]
    
    if (length(vals) == 0) {
      return(rast_absent_limit)
    }
    
    q_absent <- as.numeric(stats::quantile(vals, absent_quantile,      na.rm = TRUE))
    q_high   <- as.numeric(stats::quantile(vals, relative_to_quantile, na.rm = TRUE))
    relative_limit <- relative_fraction * q_high
    
    max(rast_absent_limit, q_absent, relative_limit, na.rm = TRUE)
  }
  
  resolved_absent_limit <- resolve_absent_limit(rasts)
  
  # ------------------------------------------------------------
  # Process each raster (threshold + per-raster upper limit)
  # ------------------------------------------------------------
  
  process_raster <- function(rast) {
    r <- rast
    r[r <= resolved_absent_limit] <- NA
    zmax <- resolve_rast_max(
      rast       = r,
      rast_max   = rast_max,
      rast_max_q = rast_max_quantile
    )
    list(rast = r, zmax = zmax)
  }
  
  processed   <- lapply(rasts, process_raster)
  zmax_values <- vapply(processed, `[[`, numeric(1), "zmax")
  
  zmax_shared <- max(zmax_values, na.rm = TRUE)
  if (!is.finite(zmax_shared) || zmax_shared <= 0) {
    zmax_shared <- 1
  }
  
  rasts_clamped <- lapply(processed, function(x) {
    terra::clamp(x$rast, upper = zmax_shared, values = TRUE)
  })
  names(rasts_clamped) <- names(rasts)
  
  zlim    <- c(resolved_absent_limit, zmax_shared)
  zbreaks <- seq(zlim[1], zlim[2], length.out = 5)
  
  list(
    rasters           = rasts_clamped,
    zlim              = zlim,
    zbreaks           = zbreaks,
    zmax_original     = zmax_values,
    rast_absent_limit = resolved_absent_limit,
    absent_method     = absent_method,
    coverage          = if (absent_method == "cumulative") coverage       else NA_real_,
    coverage_basis    = if (absent_method == "cumulative") coverage_basis else NA_character_
  )
}


# Create map of posterior mean relative abundance (or PObs) for one atlas period
make_map <- function(species_name,
                     subtitle,
                     legend_title = "Expected count\nper 5-min\npoint count",
                     rast,
                     region,
                     water      = NULL,
                     colpal     = NULL,
                     water_fill = "#b8dceb",
                     transform  = "identity",
                     zlim       = NULL,
                     zbreaks    = NULL,
                     legend_position = c(0.97, 0.97),
                     legend_justification = c(1, 1)) {
  
  region <- region |>
    sf::st_transform(terra::crs(rast))
  
  if (!is.null(water)) {
    water <- water |>
      sf::st_transform(sf::st_crs(region))
  }
  
  region_vect <- terra::vect(region)
  
  rast <- rast |>
    terra::crop(region_vect) |>
    terra::mask(region_vect)
  
  if (is.null(colpal)) {
    colpal <- grDevices::colorRampPalette(c(
      "#FBF7E2", "#FCF8D0", "#EEF7C2", "#CEF2B0",
      "#94E5A0", "#51C987", "#18A065", "#008C59",
      "#007F53", "#006344"
    ))(100)
  }
  
  rast_df <- as.data.frame(rast, xy = TRUE, na.rm = FALSE)
  names(rast_df)[3] <- "value"
  
  title_text <- paste0(
    "<span style='font-size:18pt;'><b>", species_name, "</b></span><br><br>",
    "<span style='font-size:14pt;'>",
    subtitle,
    "</span>"
  )
  
  if (is.null(zbreaks) && !is.null(zlim)) {
    zbreaks <- seq(zlim[1], zlim[2], length.out = 5)
  }
  
  ggplot2::ggplot() +
    ggplot2::geom_sf(
      data = region,
      fill = "white",
      colour = "black",
      linewidth = 0.3
    ) +
    ggplot2::geom_raster(
      data = rast_df,
      ggplot2::aes(x = x, y = y, fill = value)
    ) +
    {
      if (!is.null(water)) {
        ggplot2::geom_sf(
          data = water,
          fill = water_fill,
          colour = "transparent"
        )
      }
    } +
    ggplot2::geom_sf(
      data = region,
      fill = "transparent",
      colour = "black",
      linewidth = 0.3
    ) +
    ggspatial::annotation_scale(
      location = "bl",
      width_hint = 0.25
    ) +
    ggplot2::scale_fill_gradientn(
      limits   = zlim,
      breaks   = zbreaks,
      labels   = scales::label_number(accuracy = 0.01)(zbreaks),
      colours  = colpal,
      na.value = "transparent",
      trans    = "identity",
      name     =  legend_title,
      oob      = scales::squish,
      guide    = ggplot2::guide_colourbar(
        ticks.linewidth = 0.5,
        frame.linewidth = 0.5,
        frame.colour    = "black",
        ticks.colour    = "black"
      )
    )+
    ggplot2::coord_sf() +
    ggplot2::ggtitle(title_text) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggtext::element_markdown(lineheight = 1.1),
      
      legend.position      = legend_position,
      legend.justification = legend_justification,
      legend.background    = ggplot2::element_rect(fill = "white", colour = NA),
      
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.background  = ggplot2::element_rect(fill = "white", colour = "white"),
      
      axis.title = ggplot2::element_blank()
    )
}

# ------------------------------------------------------------
# Hex-level absolute population change maps
# ------------------------------------------------------------

summarize_hex_draw_change <- function(hex_grid,
                                      hex_draws,
                                      ci_level = 0.90,
                                      baseline_min = 1e-5,
                                      min_n_pixels = NULL,
                                      min_n_pixels_prop = 0.25) {
  
  stopifnot(inherits(hex_grid, "sf"))
  stopifnot("hex_id" %in% names(hex_grid))
  
  required_cols <- c("hex_id", "mu_OBBA2", "mu_OBBA3", "n_pixels")
  missing_cols <- setdiff(required_cols, names(hex_draws))
  
  if (length(missing_cols) > 0) {
    stop(
      "hex_draws is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  if (!is.numeric(ci_level) || length(ci_level) != 1 ||
      ci_level <= 0 || ci_level >= 1) {
    stop("ci_level must be a single number strictly between 0 and 1.")
  }
  
  if (is.null(min_n_pixels)) {
    min_n_pixels <- min_n_pixels_prop * max(hex_draws$n_pixels, na.rm = TRUE)
  }
  
  hex_draws <- hex_draws |>
    dplyr::filter(.data$n_pixels >= min_n_pixels)
  
  alpha <- (1 - ci_level) / 2
  
  n_hex <- length(hex_draws$hex_id)
  
  if (length(hex_draws$mu_OBBA2) != n_hex ||
      length(hex_draws$mu_OBBA3) != n_hex) {
    stop("hex_id, mu_OBBA2, and mu_OBBA3 must have the same length.")
  }
  
  summarize_one_hex <- function(mu2, mu3) {
    
    if (!is.numeric(mu2) || !is.numeric(mu3)) {
      stop("Each element of mu_OBBA2 and mu_OBBA3 must be numeric.")
    }
    
    if (length(mu2) != length(mu3)) {
      stop("Each paired mu_OBBA2 / mu_OBBA3 draw vector must have equal length.")
    }
    
    abs_change <- mu3 - mu2
    
    prop_change <- abs_change / mu2
    prop_change[mu2 <= baseline_min] <- NA_real_
    
    sym_change <- log(mu3 / mu2)
    sym_change[mu2 <= baseline_min] <- NA_real_
    
    tibble::tibble(
      mu2_mean = mean(mu2, na.rm = TRUE),
      mu2_median = stats::median(mu2, na.rm = TRUE),
      mu2_qlow = unname(stats::quantile(mu2, probs = alpha, na.rm = TRUE)),
      mu2_qhigh = unname(stats::quantile(mu2, probs = 1 - alpha, na.rm = TRUE)),
      
      mu3_mean = mean(mu3, na.rm = TRUE),
      mu3_median = stats::median(mu3, na.rm = TRUE),
      mu3_qlow = unname(stats::quantile(mu3, probs = alpha, na.rm = TRUE)),
      mu3_qhigh = unname(stats::quantile(mu3, probs = 1 - alpha, na.rm = TRUE)),
      
      abs_change_mean = mean(abs_change, na.rm = TRUE),
      abs_change_median = stats::median(abs_change, na.rm = TRUE),
      abs_change_qlow = unname(stats::quantile(abs_change, probs = alpha, na.rm = TRUE)),
      abs_change_qhigh = unname(stats::quantile(abs_change, probs = 1 - alpha, na.rm = TRUE)),
      
      prop_change_mean = mean(prop_change, na.rm = TRUE),
      prop_change_median = stats::median(prop_change, na.rm = TRUE),
      prop_change_qlow = unname(stats::quantile(prop_change, probs = alpha, na.rm = TRUE)),
      prop_change_qhigh = unname(stats::quantile(prop_change, probs = 1 - alpha, na.rm = TRUE)),
      
      sym_change_mean = mean(sym_change, na.rm = TRUE),
      sym_change_median = stats::median(sym_change, na.rm = TRUE),
      sym_change_qlow = unname(stats::quantile(sym_change, probs = alpha, na.rm = TRUE)),
      sym_change_qhigh = unname(stats::quantile(sym_change, probs = 1 - alpha, na.rm = TRUE)),
      
      p_increase = mean(abs_change > 0, na.rm = TRUE),
      p_decrease = mean(abs_change < 0, na.rm = TRUE)
    )
  }
  
  hex_summary <- purrr::map2_dfr(
    hex_draws$mu_OBBA2,
    hex_draws$mu_OBBA3,
    summarize_one_hex
  ) |>
    dplyr::mutate(
      hex_id = hex_draws$hex_id,
      n_pixels = hex_draws$n_pixels,
      min_n_pixels_used = min_n_pixels,
      direction = dplyr::case_when(
        p_increase > 0.975 ~ "increase",
        p_decrease > 0.975 ~ "decrease",
        TRUE ~ "uncertain"
      ),
      .before = 1
    )
  
  hex_grid |>
    dplyr::left_join(hex_summary, by = "hex_id")
}

# Area-weight of each hexagon inside a polygon (hex_id -> area_weight in [0, 1]).
#
# The weight is purely geometric: it depends only on `hex_grid` and `polygon`,
# NOT on any species' posterior draws. It is therefore identical for every
# species, so callers that summarise many species over the same set of regions
# (scripts 08/09) should compute it ONCE per region and pass it to
# summarize_polygon_hex_draw_change() via `hex_weights`, instead of paying for
# the st_intersection() below on every species.
compute_hex_polygon_weights <- function(hex_grid, polygon) {
  stopifnot(inherits(hex_grid, "sf"), inherits(polygon, "sf"))
  stopifnot("hex_id" %in% names(hex_grid))
  
  polygon <- sf::st_transform(polygon, sf::st_crs(hex_grid))
  polygon_union <- sf::st_union(sf::st_make_valid(polygon))
  
  hex_area <- sf::st_area(hex_grid)
  
  hex_intersection <- suppressWarnings(
    sf::st_intersection(
      hex_grid |> dplyr::select(hex_id),
      sf::st_sf(geometry = polygon_union)
    )
  )
  
  if (nrow(hex_intersection) == 0) {
    stop("No hexagons intersect the supplied polygon.")
  }
  
  hex_intersection |>
    dplyr::mutate(intersection_area = sf::st_area(geometry)) |>
    sf::st_drop_geometry() |>
    dplyr::group_by(hex_id) |>
    dplyr::summarise(
      intersection_area = sum(intersection_area),
      .groups = "drop"
    ) |>
    dplyr::left_join(
      tibble::tibble(hex_id = hex_grid$hex_id, hex_area = hex_area),
      by = "hex_id"
    ) |>
    dplyr::mutate(
      area_weight = as.numeric(intersection_area / hex_area),
      area_weight = pmin(pmax(area_weight, 0), 1)
    ) |>
    dplyr::select(hex_id, area_weight)
}

# Area-weighted posterior summary of OBBA2->OBBA3 change for one polygon.
#
# `hex_weights` (a tibble of hex_id + area_weight from compute_hex_polygon_weights())
# may be supplied to skip the polygon geometry entirely; this is the fast path
# used by scripts 08/09, which precompute one weight table per region and reuse
# it across every species. When `hex_weights` is NULL, `polygon` is required and
# the weights are computed on the fly (original behaviour, unchanged).
summarize_polygon_hex_draw_change <- function(hex_draws,
                                              hex_grid,
                                              polygon = NULL,
                                              ci_level = 0.90,
                                              baseline_min = 1e-5,
                                              return_draws = FALSE,
                                              hex_weights = NULL) {
  
  required_cols <- c("hex_id", "mu_OBBA2", "mu_OBBA3")
  missing_cols <- setdiff(required_cols, names(hex_draws))
  
  if (length(missing_cols) > 0) {
    stop(
      "hex_draws is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  stopifnot(inherits(hex_grid, "sf"))
  stopifnot("hex_id" %in% names(hex_grid))
  
  if (!is.numeric(ci_level) || length(ci_level) != 1 ||
      ci_level <= 0 || ci_level >= 1) {
    stop("ci_level must be a single number strictly between 0 and 1.")
  }
  
  alpha <- (1 - ci_level) / 2
  
  # ----------------------------------------------------------
  # Fraction of each hex inside the polygon (geometric; species-independent).
  # Reuse precomputed weights when supplied, else derive them from `polygon`.
  # ----------------------------------------------------------
  
  if (is.null(hex_weights)) {
    stopifnot(inherits(polygon, "sf"))
    hex_weights <- compute_hex_polygon_weights(hex_grid, polygon)
  } else {
    stopifnot(all(c("hex_id", "area_weight") %in% names(hex_weights)))
  }
  
  # ----------------------------------------------------------
  # Join weights to draws
  # ----------------------------------------------------------
  
  hex_draws_weighted <- hex_draws |>
    dplyr::inner_join(hex_weights, by = "hex_id")
  
  mu2_mat <- do.call(rbind, hex_draws_weighted$mu_OBBA2)
  mu3_mat <- do.call(rbind, hex_draws_weighted$mu_OBBA3)
  
  if (!identical(dim(mu2_mat), dim(mu3_mat))) {
    stop("mu_OBBA2 and mu_OBBA3 draw matrices do not have identical dimensions.")
  }
  
  weights <- hex_draws_weighted$area_weight
  
  # Apply area weights to each hex before summing across hexagons
  mu2_mat_weighted <- mu2_mat * weights
  mu3_mat_weighted <- mu3_mat * weights
  
  region_mu2_draws <- colSums(mu2_mat_weighted, na.rm = TRUE)
  region_mu3_draws <- colSums(mu3_mat_weighted, na.rm = TRUE)
  
  abs_change_draws <- region_mu3_draws - region_mu2_draws
  
  prop_change_draws <- abs_change_draws / region_mu2_draws
  sym_change_draws  <- log(region_mu3_draws / region_mu2_draws)
  
  prop_change_draws[region_mu2_draws <= baseline_min] <- NA_real_
  sym_change_draws[region_mu2_draws <= baseline_min]  <- NA_real_
  
  out <- tibble::tibble(
    n_hex_used = nrow(hex_draws_weighted),
    sum_area_weight = sum(weights, na.rm = TRUE),
    
    mu2_mean = mean(region_mu2_draws, na.rm = TRUE),
    mu2_median = stats::median(region_mu2_draws, na.rm = TRUE),
    mu2_qlow = unname(stats::quantile(region_mu2_draws, probs = alpha, na.rm = TRUE)),
    mu2_qhigh = unname(stats::quantile(region_mu2_draws, probs = 1 - alpha, na.rm = TRUE)),
    
    mu3_mean = mean(region_mu3_draws, na.rm = TRUE),
    mu3_median = stats::median(region_mu3_draws, na.rm = TRUE),
    mu3_qlow = unname(stats::quantile(region_mu3_draws, probs = alpha, na.rm = TRUE)),
    mu3_qhigh = unname(stats::quantile(region_mu3_draws, probs = 1 - alpha, na.rm = TRUE)),
    
    abs_change_mean = mean(abs_change_draws, na.rm = TRUE),
    abs_change_median = stats::median(abs_change_draws, na.rm = TRUE),
    abs_change_qlow = unname(stats::quantile(abs_change_draws, probs = alpha, na.rm = TRUE)),
    abs_change_qhigh = unname(stats::quantile(abs_change_draws, probs = 1 - alpha, na.rm = TRUE)),
    
    prop_change_mean = mean(prop_change_draws, na.rm = TRUE),
    prop_change_median = stats::median(prop_change_draws, na.rm = TRUE),
    prop_change_qlow = unname(stats::quantile(prop_change_draws, probs = alpha, na.rm = TRUE)),
    prop_change_qhigh = unname(stats::quantile(prop_change_draws, probs = 1 - alpha, na.rm = TRUE)),
    
    sym_change_mean = mean(sym_change_draws, na.rm = TRUE),
    sym_change_median = stats::median(sym_change_draws, na.rm = TRUE),
    sym_change_qlow = unname(stats::quantile(sym_change_draws, probs = alpha, na.rm = TRUE)),
    sym_change_qhigh = unname(stats::quantile(sym_change_draws, probs = 1 - alpha, na.rm = TRUE)),
    
    p_increase = mean(abs_change_draws > 0, na.rm = TRUE),
    p_decrease = mean(abs_change_draws < 0, na.rm = TRUE),
    
    direction = dplyr::case_when(
      p_increase > 0.975 ~ "increase",
      p_decrease > 0.975 ~ "decrease",
      TRUE ~ "uncertain"
    )
  )
  
  if (return_draws) {
    return(list(
      summary = out,
      draws = list(
        mu2 = region_mu2_draws,
        mu3 = region_mu3_draws,
        abs_change = abs_change_draws,
        prop_change = prop_change_draws,
        sym_change = sym_change_draws
      ),
      hex_weights = hex_weights
    ))
  }
  
  out
}

# ------------------------------------------------------------
# Change-map colours and legend
# ------------------------------------------------------------

classify_min_supported_change <- function(hex_sf) {
  
  required_cols <- c("sym_change_qlow", "sym_change_qhigh")
  missing_cols <- setdiff(required_cols, names(hex_sf))
  
  if (length(missing_cols) > 0) {
    stop(
      "hex_sf is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  hex_sf |>
    dplyr::mutate(
      min_supported_change = dplyr::case_when(
        
        is.na(.data$sym_change_qlow) | is.na(.data$sym_change_qhigh) ~ NA_character_,
        
        # Increases: lower CI bound above threshold
        .data$sym_change_qlow >= log(2.0) ~ "> 2x",
        .data$sym_change_qlow >= log(1.5) ~ "1.5x to 2x",
        .data$sym_change_qlow >= log(1.1) ~ "1.1x to 1.5x",
        
        # Declines: upper CI bound below threshold
        .data$sym_change_qhigh <= log(1 / 2.0) ~ "< 0.5x",
        .data$sym_change_qhigh <= log(1 / 1.5) ~ "0.5x to 0.67x",
        .data$sym_change_qhigh <= log(1 / 1.1) ~ "0.67x to 0.9x",
        
        TRUE ~ "0.9x to 1.1x"
      ),
      
      min_supported_change = factor(
        min_supported_change,
        levels = c(
          "< 0.5x",
          "0.5x to 0.67x",
          "0.67x to 0.9x",
          "0.9x to 1.1x",
          "1.1x to 1.5x",
          "1.5x to 2x",
          "> 2x"
        ),
        ordered = TRUE
      )
    )
}

get_min_supported_change_scale <- function() {
  
  change_levels <- c(
    "< 0.5x", "0.5x to 0.67x", "0.67x to 0.9x",
    "0.9x to 1.1x",
    "1.1x to 1.5x", "1.5x to 2x", "> 2x"
  )
  
  change_labels <- c(
    "-50%", "-33%", "-10%",
    "Weak or\nuncertain",
    "+10%", "+50%", "+100%"
  )
  
  tweak_colour <- function(col, lighten_amount = -0.15, saturate_amount = 0.2) {
    hsl <- as(colorspace::hex2RGB(col), "HLS")
    hsl@coords[, "L"] <- pmin(1, pmax(0, hsl@coords[, "L"] + lighten_amount))
    hsl@coords[, "S"] <- pmin(1, pmax(0, hsl@coords[, "S"] + saturate_amount))
    colorspace::hex(hsl)
  }
  
  change_colours <- RColorBrewer::brewer.pal(7, "RdBu")
  change_colours[4] <- "white"
  change_colours[-4] <- sapply(
    change_colours[-4],
    tweak_colour,
    lighten_amount = -0.1,
    saturate_amount = 0.1
  )
  
  list(
    levels = change_levels,
    labels = change_labels,
    colours = change_colours
  )
}


make_min_supported_change_legend <- function(ci_level = 0.90) {
  
  change_scale <- get_min_supported_change_scale()
  
  legend_df <- tibble::tibble(
    level = change_scale$levels,
    label = change_scale$labels,
    y = seq_along(change_scale$levels)
  )
  
  legend_title <- paste0(
    "Population changed\nby at least this much\n(",
    round(ci_level * 100),
    "% confident)"
  )
  
  title_y <- max(legend_df$y) + 1.6
  
  ggplot2::ggplot(legend_df, ggplot2::aes(x = 1, y = y)) +
    ggplot2::geom_tile(
      ggplot2::aes(fill = level),
      colour = "grey40",
      width = 0.6,
      height = 0.8,
      linewidth = 0.3
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = 1.55, label = label),
      hjust = 0,
      size = 3
    ) +
    ggplot2::annotate(
      "text",
      x = 0.6,
      y = title_y,
      label = legend_title,
      hjust = 0,
      vjust = 1,
      size = 3.5,
      fontface = "plain",
      lineheight = 0.9
    ) +
    ggplot2::scale_fill_manual(
      values = stats::setNames(change_scale$colours, change_scale$levels),
      guide = "none"
    ) +
    ggplot2::coord_cartesian(
      xlim = c(0.5, 6.5),
      ylim = c(0.5, title_y + 0.2),
      clip = "off"
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.margin = ggplot2::margin(2, 8, 2, 2)
    )
}

# ------------------------------------------------------------
# Make change map
# ------------------------------------------------------------

make_hex_abs_change_map <- function(species_name,
                                    hex_change_sf,
                                    region,
                                    region_boundaries = NULL,
                                    prov_change = NULL,
                                    ci_level = 0.90,
                                    size_var = "mu_mid_median",
                                    zlim = NULL,
                                    filter_low_abundance = TRUE,
                                    water = NULL,
                                    water_fill = "#b8dceb",
                                    legend_inset = c(
                                      left = 0.8,
                                      bottom = 0.5,
                                      right = 1.0,
                                      top = 0.99
                                    )) {
  
  stopifnot(inherits(hex_change_sf, "sf"))
  stopifnot(inherits(region, "sf"))
  stopifnot("min_supported_change" %in% names(hex_change_sf))
  
  region <- sf::st_transform(region, sf::st_crs(hex_change_sf))
  
  if (!is.null(water)) {
    water <- sf::st_transform(water, sf::st_crs(hex_change_sf))
  }
  
  if (!is.null(region_boundaries)) {
    region_boundaries <- sf::st_transform(region_boundaries, sf::st_crs(hex_change_sf)) %>%
      st_intersection(region)
  }
  
  # Show hexes where the species was present in either atlas
  hex_change_sf <- hex_change_sf |>
    dplyr::mutate(
      mu_max_median = pmax(mu2_median, mu3_median, na.rm = TRUE),
      mu_mid_median = (mu2_median + mu3_median) / 2
    )
  
  stopifnot(size_var %in% names(hex_change_sf))
  
  if (filter_low_abundance) {
    
    if (is.null(zlim) || length(zlim) < 1 || !is.finite(zlim[1])) {
      stop("zlim must be supplied when filter_low_abundance = TRUE")
    }
    
    absence_threshold <- zlim[1]
    
    hex_change_sf <- hex_change_sf |>
      dplyr::filter(mu_max_median > absence_threshold)
  }
  
  change_scale <- get_min_supported_change_scale()
  
  # Match old title formatting
  fmt_pct <- function(x, digits = 1) {
    ifelse(
      is.na(x),
      NA_character_,
      sprintf("%+.1f", round(100 * x, digits))
    )
  }
  
  if (!is.null(prov_change)) {
    title_text <- paste0(
      "<span style='font-size:18pt;'><b>", species_name, "</b></span><br><br>",
      "<span style='font-size:14pt;'>",
      "Overall change = ",
      fmt_pct(prov_change$prop_change_median), "% [",
      fmt_pct(prov_change$prop_change_qlow), "% to ",
      fmt_pct(prov_change$prop_change_qhigh), "%]",
      "</span>"
    )
  } else {
    title_text <- paste0(
      "<span style='font-size:18pt;'><b>", species_name, "</b></span>"
    )
  }
  
  legend_plot <- make_min_supported_change_legend(ci_level = ci_level)
  
  chg_plot <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = region, fill = "gray90")
  
  if (!is.null(water)) {
    chg_plot <- chg_plot +
      ggplot2::geom_sf(
        data = water,
        fill = water_fill,
        colour = "transparent"
      )
  }
  
  if (!is.null(region_boundaries)) {
    chg_plot <- chg_plot +
      ggplot2::geom_sf(
        data = region_boundaries,
        fill = "transparent",
        colour = "gray30"
      )
  }
  
  chg_plot <- chg_plot +
    ggplot2::geom_sf(
      data = sf::st_centroid(hex_change_sf),
      ggplot2::aes(
        colour = min_supported_change,
        size = .data[[size_var]]
      )
    ) +
    ggplot2::geom_sf(
      data = region,
      fill = "transparent",
      colour = "black",
      linewidth = 0.3
    ) +
    ggspatial::annotation_scale(
      location = "bl",
      width_hint = 0.25
    ) +
    ggplot2::scale_colour_manual(
      values = stats::setNames(change_scale$colours, change_scale$levels),
      limits = change_scale$levels,
      drop = FALSE,
      guide = "none"
    ) +
    ggplot2::scale_size_continuous(
      range = c(0, 2),
      guide = "none"
    ) +
    ggplot2::theme_bw() +
    ggplot2::ggtitle(title_text) +
    ggplot2::theme(
      plot.title = ggtext::element_markdown(lineheight = 1.1),
      panel.background = ggplot2::element_rect(fill = "transparent", colour = NA),
      plot.background  = ggplot2::element_rect(fill = "transparent", colour = NA),
      axis.title = ggplot2::element_blank()
    )
  
  chg_plot +
    patchwork::inset_element(
      legend_plot,
      left = legend_inset[["left"]],
      bottom = legend_inset[["bottom"]],
      right = legend_inset[["right"]],
      top = legend_inset[["top"]],
      on_top = TRUE
    )
}


# Assess model predictions against survey observations at any spatial scale
#
# Orchestrates the full workflow:
#   1. Accept a pre-built hexagon grid (or build one if `n_hexagons` is
#      supplied and `hex_grid` is NULL).
#   2. Summarise survey observations and raster predictions per hexagon
#      (`summarize_hex`).
#   3. Compute region-wide validation statistics (`compute_region_stats`).
#   4. Build a two-panel patchwork figure:
#      - Panel A: raw raster of model predictions (`plot_raster_gg`).
#      - Panel B: honeycomb effort/detection map (`plot_honeycomb`).
#
# When looping over many species for the same region, build the hex grid once
# with `make_hex_grid()` and pass it via `hex_grid` to avoid rebuilding it on
# every iteration.
#
# Argument: region          `sf` polygon defining the study area extent. Used for
#   plotting and, when `hex_grid` is NULL, for building the grid.
# Argument: sp_dat          `sf` data frame with a `count_per_effort` column.
# Argument: rast            `SpatRaster` of model predictions. Absent cells should
#   be 0 or NA; both are treated as zero when computing hex means, and rendered
#   as transparent in all map panels.
# Argument: hex_grid        Pre-built hexagon grid (`sf` output of
#   `make_hex_grid()`). When supplied, `n_hexagons` is ignored and no grid is
#   built internally. Pass `NULL` (default) to build the grid from `region`
#   and `n_hexagons`.
# Argument: n_hexagons      Target number of hexagons. Only used when `hex_grid`
#   is NULL. Default `1000`.
# Argument: water           Optional `processed_water` object (from
#                        `process_water()`). Pass `NULL` to omit water.
# Argument: rast_max_q      Quantile (0-1) used to cap the colour scale upper end
#   and clip outlier pixel values before extraction. Default `0.99`.
# Argument: transform       Colour scale transform: `"identity"`, `"sqrt"`,
#                        `"log"`, etc.
# Argument: pred_presence_threshold  `pred_mean` values strictly above this are
#   classified as model-predicted present in `hex_summary` and used for P/A
#   statistics. Default `0`.
# Argument: title           Optional character string added as a patchwork
#                        annotation title above the three panels.
# Argument: model_source    Optional character string identifying the model source
#                        (e.g. `"BAM"`, `"Atlas"`). Displayed in the subtitle
#                        when either `model_source` or `data_source` is
#                        non-NULL.
# Argument: data_source     Optional character string identifying the survey data
#                        source (e.g. `"BAM"`, `"Atlas"`). Displayed in the
#                        subtitle alongside `model_source`.
#
# Returns: A named list:
#   \describe{
#     \item{`plot_combined`}{A `patchwork` ggplot (two panels).}
#     \item{`hex_summary`}{`sf` data frame with one row per hexagon
#       containing columns from `summarize_hex()`.}
#     \item{`region_stats`}{One-row `data.frame` from `compute_region_stats()`.}
#   }
assess_region <- function(region,
                          region_boundaries       = NULL,
                          sp_dat,
                          rast,
                          hex_grid                = NULL,
                          n_hexagons              = 1000,
                          water                   = NULL,
                          rast_max_q              = 0.99,
                          rast_max                = NULL,
                          rast_min                = NULL,
                          max_count_per_effort    = NULL,
                          count_max_q             = 0.99,
                          max_surveys             = NULL,
                          max_surveys_q           = 0.80,
                          transform               = "identity",
                          pred_presence_threshold = 0,
                          title                   = NULL,
                          model_source            = NULL,
                          data_source             = NULL,
                          color_palette           = NULL,
                          count_size_label        = "observed mean count in hexagon"
                          ) {
  
  stopifnot(inherits(region, "sf"))
  stopifnot(inherits(sp_dat, "sf"))
  stopifnot(inherits(rast, "SpatRaster"))
  if (!is.null(hex_grid)) stopifnot(inherits(hex_grid, "sf"))
  
  region <- sf::st_transform(region, sf::st_crs(sp_dat))
  sp_dat <- sf::st_filter(sp_dat, region, .predicate = sf::st_intersects)
  
  if (!is.null(region_boundaries)) {
    region_boundaries <- region_boundaries |>
      sf::st_transform(sf::st_crs(region)) |>
      sf::st_intersection(region)
  }
  
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
    rast_max                = rast_max,
    pred_presence_threshold = pred_presence_threshold
  )
  
  # -- GOAL 3: Region-wide statistics ----------------------------------------
  region_stats <- compute_region_stats(hex_summary)
  
  # -- GOAL 3b: Resolve display-scale limits ---------------------------------
  # These values are passed to all relevant panels so scale limits are not
  # recomputed independently within each plotting function. To force identical
  # scaling across atlas periods, calculate shared values outside this function
  # and pass them through rast_max, max_count_per_effort, and max_surveys.
  scale_limits <- compute_assessment_scale_limits(
    hex_summaries        = list(hex_summary),
    rasts                = list(rast_cropped),
    rast_max             = rast_max,
    rast_max_q           = rast_max_q,
    max_count_per_effort = max_count_per_effort,
    count_max_q          = count_max_q,
    max_surveys          = max_surveys,
    max_surveys_q        = max_surveys_q
  )
  
  # -- GOAL 4a: Panel A – raw raster -----------------------------------------
  p_rast <- plot_raster_gg(
    rast       = rast_cropped,
    study_area = region,
    region_boundaries = region_boundaries,
    water      = NULL,        # Water is hidden on the raster
    rast_max_q = rast_max_q,
    rast_max   = scale_limits$rast_max,
    rast_min   = rast_min,
    transform  = transform,
    palette = color_palette
  )
  
  # -- GOAL 4b: Panel C – honeycomb effort/detection -------------------------
  p_honey <- plot_honeycomb(
    hex_summary          = hex_summary,
    study_area           = region,
    region_boundaries = region_boundaries,
    water                = water,
    max_surveys          = scale_limits$max_surveys,
    max_surveys_q        = max_surveys_q,
    max_count_per_effort   = scale_limits$max_count_per_effort,
    max_count_per_effort_q = count_max_q,
    count_size_label       = count_size_label
  )
  
  # -- Compose figure --------------------------------------------------------
  
  combined <-
    (p_rast           + margin_theme) +
    (p_honey          + margin_theme) +
    patchwork::plot_layout(ncol = 2)
  
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
    p_rast = p_rast,
    p_honey = p_honey,
    plot_combined = combined,
    hex_summary   = hex_summary,
    region_stats  = region_stats,
    scale_limits  = scale_limits
  )
}

# Create a hexagon grid over a study area
#
# Tiles the study area with regular hexagons so that the total number of
# hexagons that *intersect* the boundary is approximately `n_hexagons`.
# Internally uses EPSG:3978 (Canada Albers) when the input is geographic, then
# reprojects back to the original CRS before returning.
#
# Argument: study_area  An `sf` object defining the region of interest.
# Argument: n_hexagons  Target number of hexagons inside the study area (default 200).
# Returns: An `sf` object with columns `hex_id` (integer) and `geometry`.
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

# Combine observed-survey and model-prediction summaries at the hexagon level
#
# A convenience wrapper that runs `summarize_surveys_by_hex()` and
# `extract_hex_predictions()` in sequence and merges both sets of columns onto
# the hex grid. Also adds a `pred_detected` binary column (TRUE when
# `pred_mean` strictly exceeds `pred_presence_threshold`).
#
# Absent pixels (0 or NA) in the raster are included in the hex mean
# denominator, so `pred_mean` reflects true spatial average abundance across
# the whole hexagon area (see `extract_hex_predictions()` for details).
#
# Argument: dat                    An `sf` object with a `count_per_effort` column.
# Argument: hex_grid               Output of `make_hex_grid()`.
# Argument: rast                   A `SpatRaster` of model predictions.
# Argument: rast_max_q             Upper quantile cap passed to
#   `extract_hex_predictions()`. Default `0.99`.
# Argument: pred_presence_threshold `pred_mean` values strictly above this are
#   classified as model-predicted present. Default `0` (any positive mean
#   abundance counts as presence).
# Returns: An `sf` object (one row per hexagon) with columns:
#   `hex_id`, `n_surveys`, `mean_count_per_effort`, `obs_detected`,
#   `pred_mean`, `pred_detected`, `geometry`.
summarize_hex <- function(dat,
                          hex_grid,
                          rast,
                          rast_max_q              = 1,
                          rast_max                = NULL,
                          pred_presence_threshold = 0) {
  
  hex_obs  <- summarize_surveys_by_hex(dat, hex_grid)
  
  hex_full <- extract_hex_predictions(
    hex_obs,
    rast,
    rast_max_q = rast_max_q,
    rast_max   = rast_max
  )
  
  hex_full |>
    dplyr::mutate(
      pred_detected = dplyr::if_else(
        is.finite(pred_mean) & pred_mean > pred_presence_threshold,
        TRUE, FALSE
      )
    )
}

# Summarise point-count / ARU surveys by hexagon
#
# Spatially joins survey points to hexagons and computes, for each hexagon:
# the number of surveys, the mean count-per-effort, and a binary
# presence/absence flag.
#
# Argument: dat      An `sf` object with a `count_per_effort` numeric column.
# Argument: hex_grid An `sf` hexagon grid (output of `make_hex_grid()`).
# Returns: The full hex grid with added columns:
#   \describe{
#     \item{`n_surveys`}{Number of survey visits (0 for unsurveyed hexagons).}
#     \item{`mean_count_per_effort`}{Mean count per effort (NA for unsurveyed).}
#     \item{`obs_detected`}{Logical: TRUE if any survey detected the species.}
#   }
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

# Extract raster model predictions into each hexagon
#
# Computes the mean predicted abundance within each hexagon, treating absent
# pixels (0 or NA in the input raster) as zero so that partial occupancy is
# correctly reflected in the hex mean. For example, if 50 % of pixels in a
# hexagon are absent (0/NA) and the remaining 50 % have a predicted value of
# 1, the returned `pred_mean` will be 0.5.
#
# The upper-quantile cap (`rast_max_q`) is applied to non-zero pixels before
# extraction so that extreme outliers do not distort hex means, but absent
# pixels are never excluded from the denominator.
#
# Argument: hex_grid     An `sf` hexagon grid (output of `make_hex_grid()`).
# Argument: rast         A `SpatRaster` of model predictions. Absent cells should
#   be encoded as 0 or NA; both are treated as zero during averaging.
# Argument: rast_max_q   Upper quantile used to cap non-zero pixel values before
#   extraction. Default `0.99`.
# Returns: The hex grid with an added column `pred_mean` (mean predicted
#   relative abundance within the hexagon, inclusive of absent pixels; NA only
#   where the hexagon falls entirely outside the raster extent).
extract_hex_predictions <- function(hex_grid,
                                    rast,
                                    rast_max_q = 1,
                                    rast_max = NULL) {
  
  stopifnot(inherits(hex_grid, "sf"))
  stopifnot(inherits(rast, "SpatRaster"))
  
  # Align CRS
  hex_proj  <- sf::st_transform(hex_grid, terra::crs(rast))
  rast_crop <- terra::crop(rast, terra::vect(hex_proj))
  
  # Cap extreme non-zero values (outlier suppression only; does not affect
  # which pixels count as absent). The cap can be either an explicit maximum
  # or a quantile-derived maximum.
  r_max <- resolve_rast_max(
    rast       = rast_crop,
    rast_max   = rast_max,
    rast_max_q = rast_max_q
  )
  rast_crop <- terra::clamp(rast_crop, upper = r_max, values = TRUE)
  
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

# Compute region-wide model-validation statistics
#
# Takes a hex-level summary (output of `summarize_hex()`) and computes
# summary statistics across all hexagons that were actually surveyed.
#
# **Statistics returned**
# | Name | Description |
# |------|-------------|
# | `n_hex_surveyed` | Hexagons with >=1 survey |
# | `n_hex_detected` | Hexagons where species was detected |
# | `prop_hex_detected` | Detection rate (detected / surveyed) |
# | `cor_obs_pred` | Pearson r between mean observed count and mean predicted abundance (unweighted) |
# | `cor_obs_pred_wtd` | Pearson r weighted by number of surveys per hexagon |
# | `pa_agreement` | Proportion of surveyed hexagons where observed P/A matches predicted P/A |
# | `pa_sensitivity` | True-positive rate: model predicted present when observed present |
# | `pa_specificity` | True-negative rate: model predicted absent when observed absent |
# | `n_surveys_total` | Total individual survey visits across region |
# | `mean_surveys_per_hex` | Average surveys per surveyed hexagon |
# | `prop_det_in_pred_absent` | Of hexagons where species was detected, proportion where mean model prediction was at or below `absence_threshold` |
# | `absence_threshold` | The threshold value used to define model-predicted absence (echoed from input) |
# | `prop_absent_in_pred_high` | Of hexagons where model predicted high abundance (>= `high_abundance_threshold`), proportion where species was not detected |
# | `high_abundance_threshold` | The threshold used to define high predicted abundance; either user-supplied or the 75th percentile of non-zero `pred_mean` values (echoed from input) |
#
# Argument: hex_summary       Output of `summarize_hex()` (or equivalent `sf` with
#   columns `n_surveys`, `mean_count_per_effort`, `obs_detected`,
#   `pred_mean`, `pred_detected`).
# Argument: absence_threshold      Numeric scalar. Hexagons with `pred_mean` at
#   or below this value (or with no raster coverage, i.e. `pred_mean` is NA)
#   are considered model-predicted absent when computing
#   `prop_det_in_pred_absent`. Default `0`.
# Argument: high_abundance_threshold  Numeric scalar or `NULL`. Hexagons with
#   `pred_mean` at or above this value are considered high-abundance when
#   computing `prop_absent_in_pred_high`. When `NULL` (default), the threshold
#   is set automatically to the 75th percentile of non-zero `pred_mean` values
#   across surveyed hexagons.
# Returns: A one-row `data.frame` of region-wide statistics.
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

# -- Panel A: raw raster -------------------------------------------------------

# Plot model-predicted relative abundance as a raster layer
#
# Pixels with a value of 0 or NA are rendered as transparent (absent).
# All non-zero pixels are mapped onto the colour ramp, with the upper end
# capped at the `rast_max_q` quantile to prevent outliers from compressing
# the scale.
#
# Argument: rast         A `SpatRaster`. Absent cells should be 0 or NA.
# Argument: study_area   `sf` polygon defining the region boundary.
# Argument: water        Optional `processed_water` object.
# Argument: rast_max_q   Quantile (0-1) used to cap the upper end of the colour
#   scale. Default `0.99`.
# Argument: palette      Character vector of colours (auto if NULL).
# Argument: water_fill   Fill colour for water polygons.
# Argument: transform    Scale transform: `"identity"`, `"sqrt"`, `"log"`, etc.
# Returns: A `ggplot` object.
plot_raster_gg <- function(rast,
                           study_area,
                           region_boundaries = NULL,
                           water      = NULL,
                           rast_max_q = 1,
                           rast_max   = NULL,
                           rast_min   = NULL,
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
    plot_crs   <- water$crs
    study_area <- sf::st_transform(study_area, plot_crs)
    rast       <- terra::project(rast, plot_crs$wkt)
    water_to_plot <- water %>% st_intersection(study_area)
  } else {
    plot_crs <- sf::st_crs(study_area)
    rast     <- terra::project(rast, plot_crs$wkt)
  }
  
  rast <- terra::crop(rast, terra::vect(study_area))
  rast <- terra::mask(rast, terra::vect(study_area))
  
  # Blank out absent pixels (0 or NA to NA so they render as transparent)
  rast[rast <= 0] <- NA
  
  # Cap the upper end of the scale. The cap can be either an explicit
  # maximum or a quantile-derived maximum.
  r_max <- resolve_rast_max(
    rast       = rast,
    rast_max   = rast_max,
    rast_max_q = rast_max_q
  )
  rast <- terra::clamp(rast, upper = r_max, values = TRUE)
  
  # Lower end of the colour scale. Absent pixels were set to NA upstream
  # (values at or below the absence floor), so the faintest palette colour is
  # anchored to the lowest *present* density rather than to 0. This uses the
  # full ramp -- including the pale low end -- instead of compressing every
  # present pixel into the dark upper colours. A shared value can be supplied
  # via `rast_min` so paired (e.g. two-period) panels use an identical scale;
  # otherwise it is derived from this raster's surviving values.
  if (!is.null(rast_min) && is.finite(rast_min)) {
    r_min <- rast_min
  } else {
    present_vals <- terra::values(rast, mat = FALSE)
    present_vals <- present_vals[is.finite(present_vals) & present_vals > 0]
    r_min <- if (length(present_vals)) min(present_vals) else 0
  }
  if (!is.finite(r_min) || r_min >= r_max) r_min <- 0
  
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
  
  if (!is.null(region_boundaries)) {
    p <- p +
      ggplot2::geom_sf(
        data = region_boundaries,
        fill = "transparent",
        colour = "gray30"
      )
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
      limits   = c(r_min, r_max),
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


# -- Panel C: honeycomb effort/detection map -----------------------------------

# Plot a honeycomb survey-effort / detection map
#
# Fills hexagons by survey effort (number of visits) and overlays circles
# scaled by observed mean count per effort.
#
# Argument: hex_summary          Output of `summarize_hex()` (or
#   `summarize_surveys_by_hex()`).
# Argument: study_area           `sf` polygon.
# Argument: water                Optional `processed_water` object.
# Argument: max_surveys          Upper bound for effort colour classes (auto).
# Argument: max_count_per_effort Upper bound for circle scaling (auto).
# Argument: alpha_min,alpha_max  Range of hex fill transparency.
# Argument: hex_fill             Fill colour for hexagons.
# Argument: circle_fill          Fill colour for detection circles.
# Argument: water_fill           Fill colour for water polygons.
# Returns: A `ggplot` object.
plot_honeycomb <- function(hex_summary,
                           study_area,
                           region_boundaries    = NULL,
                           water                = NULL,
                           max_surveys          = NULL,
                           max_surveys_q        = 0.80,
                           max_count_per_effort = NULL,
                           max_count_per_effort_q = 0.99,
                           alpha_min            = 0.15,
                           alpha_max            = 1,
                           hex_fill             = "#B38F47", 
                           circle_fill          = "black",
                           water_fill           = "#b8dceb",
                           # Text used in the panel subtitle to describe what the
                           # circle size encodes. Defaults to raw mean count; set
                           # to a "per survey-minute" wording when count_per_effort
                           # is a per-minute rate (see prepare_honeycomb_surveys()).
                           count_size_label     = "observed mean count in hexagon") {
  
  stopifnot(inherits(hex_summary, "sf"))
  stopifnot(inherits(study_area, "sf"))
  stopifnot("n_surveys" %in% names(hex_summary))
  stopifnot("mean_count_per_effort" %in% names(hex_summary))
  
  # water_to_plot <- NULL
  # 
  # if (!is.null(water)) {
  #   stopifnot(inherits(water, "processed_water"))
  #   plot_crs    <- water$crs
  #   study_area  <- sf::st_transform(study_area, plot_crs)
  #   hex_summary <- sf::st_transform(hex_summary, plot_crs)
  #   water_to_plot <- resolve_water_layer(water, study_area)
  # } else {
  #   plot_crs   <- sf::st_crs(hex_summary)
  #   study_area <- sf::st_transform(study_area, plot_crs)
  # }
  
  water_to_plot <- NULL
  
  if (!is.null(water)) {
    
    if (inherits(water, "sf")) {
      plot_crs <- sf::st_crs(study_area)
      water_to_plot <- water |>
        sf::st_transform(plot_crs) |>
        sf::st_intersection(study_area)
      
    } else if (inherits(water, "processed_water")) {
      plot_crs <- water$crs
      study_area <- sf::st_transform(study_area, plot_crs)
      hex_summary <- sf::st_transform(hex_summary, plot_crs)
      water_to_plot <- resolve_water_layer(water, study_area)
      
    } else {
      stop("water must be NULL, an sf object, or a processed_water object.")
    }
    
  } else {
    plot_crs <- sf::st_crs(hex_summary)
    study_area <- sf::st_transform(study_area, plot_crs)
  }
  
  max_surveys <- ceiling(resolve_positive_max(
    x             = hex_summary$n_surveys,
    max_value     = max_surveys,
    max_q         = max_surveys_q,
    default       = 1,
    argument_name = "max_surveys"
  ))
  if (!is.finite(max_surveys) || max_surveys < 1) max_surveys <- 1
  
  max_count_per_effort <- resolve_positive_max(
    x             = hex_summary$mean_count_per_effort,
    max_value     = max_count_per_effort,
    max_q         = max_count_per_effort_q,
    default       = 1,
    argument_name = "max_count_per_effort"
  )
  
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
    ggplot2::geom_sf(data = study_area, fill = NA, colour = "gray30")
  
  if (!is.null(water_to_plot)) {
    p <- p + ggplot2::geom_sf(data = water_to_plot,
                              fill = water_fill, colour = water_fill,
                              linewidth = 0.1)
  }
  
  
  
  p <- p +
    ggplot2::geom_sf(
      data = hex_summary_plot,
      ggplot2::aes(alpha = effort_class),
      fill = hex_fill, colour = NA
    )
  
  if (!is.null(region_boundaries)) {
    p <- p +
      ggplot2::geom_sf(
        data = region_boundaries,
        fill = "transparent",
        colour = "gray30"
      )
  }
  
  p <- p +  
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
        paste0("Circle size = ", count_size_label),
        sep = "\n"
      )
    )
}

# ==============================================================================
# Panel plot helpers
# ==============================================================================

# -- Internal: effort-class breaks and alpha scale ----------------------------

# Build effort-class breaks and alpha scale for the honeycomb plot
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


# Plot region-specific estimates of population change
make_change_comparison_plot <- function(
    regional_estimates_FullModel,
    paired_change_summary,
    # South -> north up the y axis. Must contain every BCR_Label that can appear
    # in the inputs; 08 passes its own copy so the plot and the repeated-survey
    # table cannot drift apart.
    region_order = c(
      "Lower Great Lakes / St. Lawrence Plain",
      "Temperate Mixed",
      "Boreal Shield West",
      "Boreal Shield East",
      "Hudson Plains"),
    title = "Comparison of population change estimates",
    subtitle = "Median and 95% credible intervals"
) {
  
  pct_to_logratio <- function(x) log1p(x / 100)
  
  pct_label <- function(x) {
    pct <- 100 * (exp(x) - 1)
    dplyr::case_when(
      pct > 0  ~ paste0("+", round(pct), "%"),
      pct == 0 ~ "0%",
      TRUE     ~ paste0(round(pct), "%")
    )
  }
  
  x_breaks_pct <- c(-80, -50, -33, 0, 50, 100, 400)
  x_breaks <- pct_to_logratio(x_breaks_pct)
  
  region_lookup <- regional_estimates_FullModel %>%
    dplyr::select(
      BCR = Region_Number,
      Region_Name
    )
  
  full_model_plot_df <- regional_estimates_FullModel %>%
    dplyr::transmute(
      BCR = Region_Number,
      Region_Name,
      analysis = "Full spatial model",
      pct_change_median = pct_change_median,
      pct_change_qlow   = pct_change_qlow,
      pct_change_qhigh  = pct_change_qhigh
    )
  
  paired_plot_df <- paired_change_summary %>%
    dplyr::left_join(region_lookup, by = "BCR") %>%
    dplyr::mutate(
      analysis = dplyr::case_when(
        shared_radius_km == 0.1 ~ "Repeated (100 m buffer)",
        TRUE ~ paste0("Shared footprint = ", shared_radius_km * 1000, " m")
      )
    ) %>%
    dplyr::transmute(
      BCR,
      Region_Name,
      analysis,
      pct_change_median = pct_change_median,
      pct_change_qlow   = pct_change_qlow,
      pct_change_qhigh  = pct_change_qhigh
    )
  
  plot_df <- dplyr::bind_rows(
    full_model_plot_df,
    paired_plot_df
  )
  
  # factor() would turn an unlisted region into NA and drop it from the panel
  # without warning. Fail loudly instead, in the same spirit as
  # require_bcr_key(): a region silently missing from the figure is much worse
  # than a stopped run. This is what catches a stale region_order after a BCR
  # definition changes (e.g. the 76_77 split).
  unknown_regions <- setdiff(
    unique(stats::na.omit(plot_df$Region_Name)),
    region_order
  )
  
  if (length(unknown_regions) > 0) {
    stop(
      "Region label(s) not present in region_order: ",
      paste(unknown_regions, collapse = ", "),
      ". Update region_order (and the BCR_Label levels in ",
      "08_generate_model_products.R) to match the current BCR definitions.",
      call. = FALSE
    )
  }
  
  plot_df <- plot_df %>%
    dplyr::mutate(
      Region_Name = factor(Region_Name, levels = region_order),
      region_y = as.numeric(Region_Name),
      analysis = factor(
        analysis,
        levels = c(
          "Full spatial model",
          "Repeated (100 m buffer)"
        )
      ),
      y_pos = dplyr::case_when(
        analysis == "Full spatial model"    ~ region_y + 0.1,
        analysis == "Repeated (100 m buffer)"  ~ region_y - 0.1
      ),
      
      # Transform percent change to symmetric log-ratio scale
      x_med  = pct_to_logratio(pct_change_median),
      x_low  = pct_to_logratio(pct_change_qlow),
      x_high = pct_to_logratio(pct_change_qhigh)
    )
  
  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      y = y_pos,
      x = x_med,
      xmin = x_low,
      xmax = x_high,
      colour = analysis
    )
  ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linewidth = 0.4,
      linetype = "dashed",
      colour = "grey40"
    ) +
    ggplot2::geom_pointrange(
      linewidth = 1,
      fatten = 2.5
    ) +
    ggplot2::scale_x_continuous(
      breaks = x_breaks,
      labels = pct_label
    ) +
    ggplot2::scale_y_continuous(
      breaks = seq_along(region_order),
      labels = region_order
    ) +
    ggplot2::scale_colour_manual(
      values = c(
        "Full spatial model"    = "black",
        "Repeated (100 m buffer)"  = "#1FAA59"
      )
    ) +
    ggplot2::labs(
      x = "Estimated population change (%)",
      y = NULL,
      colour = "",
      title = title,
      subtitle = subtitle
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "top",
      plot.title = ggplot2::element_text(face = "bold")
    )
}

make_paired_summary_table <- function(paired_data_summary,
                                      title = "Repeated survey data summary") {
  
  table_df <- paired_data_summary %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      Atlas = dplyr::recode(
        Atlas,
        OBBA2 = "Atlas 2",
        OBBA3 = "Atlas 3"
      ),
      line = paste0(
        Atlas,
        ": mean count = ", sprintf("%.3f", mean_count),
        "   |   PObs = ", sprintf("%.3f", PObs),
        "   |   n_svy = ", scales::comma(n_surveys),
        "   |   n_sq = ", scales::comma(n_squares)
      )
    ) %>%
    dplyr::arrange(dplyr::desc(BCR_Label), Atlas) %>%
    dplyr::group_by(BCR_Label) %>%
    dplyr::summarise(
      text = paste0(
        "<b>", unique(BCR_Label), "</b>",
        "<br>",
        paste(line, collapse = "<br>")
      ),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(BCR_Label)) %>%
    dplyr::mutate(
      y = dplyr::row_number()
    )
  
  ggplot2::ggplot(table_df, ggplot2::aes(x = 0, y = -y, label = text)) +
    ggtext::geom_richtext(
      hjust = 0,
      vjust = 1,
      lineheight = 1.05,
      size = 3.2,
      fill = NA,
      label.color = NA
    ) +
    ggplot2::labs(title = title) +
    ggplot2::coord_cartesian(
      xlim = c(0, 1),
      ylim = c(-max(table_df$y) - 0.2, -0.5),
      expand = FALSE,
      clip = "off"
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        size = 13,
        hjust = 0,
        margin = ggplot2::margin(b = 6)
      ),
      plot.margin = ggplot2::margin(5, 5, 5, 5)
    )
}


# ==============================================================================
# Per-region / per-page assembly helpers
#
# These orchestrate the region- and page-level products consumed by
# 08_generate_model_products.R. They were previously defined inline in that
# script; moved here so 08 holds only configuration and top-level flow.
# They take all inputs as explicit arguments (no dependence on script globals
# or the `cfg` list); numeric defaults below mirror the cfg values 08 passes.
# ==============================================================================

# Region geometries + their assessment hex grids. These depend only on the
# study area and n_hexagons, so they are constant across all species.
make_region_geom <- function(bcr_codes, bcr_regions, study_boundary) {
  if (is.null(bcr_codes)) {
    return(study_boundary)
  }
  missing_codes <- setdiff(bcr_codes, bcr_regions$BCR)
  if (length(missing_codes) > 0) {
    stop("BCR code(s) not found in bcr_regions: ",
         paste(missing_codes, collapse = ", "), call. = FALSE)
  }
  subset(bcr_regions, BCR %in% bcr_codes) %>%
    sf::st_union() %>%
    sf::st_intersection(study_boundary) %>%
    sf::st_transform(sf::st_crs(study_boundary)) %>%
    sf::st_as_sf()
}


# Reconstruct a species' survey sf from its lightweight data_used record.
# The data_used file stores only survey_id + count; everything else (geometry,
# Survey_Type, Atlas, BCR, square_id, ...) is rejoined from all_surveys.
load_species_surveys <- function(dat_used, all_surveys) {
  all_surveys %>%
    dplyr::inner_join(dat_used$survey_counts, by = "survey_id") %>%
    dplyr::mutate(count_per_effort = count)
}


# Survey_Type values, matching the labels finalized in 06_filter_and_finalize
# _surveys.R. Structured surveys have (approximately) fixed effort; checklist
# surveys have variable effort and are OBBA3-only. Checklist surveys are present
# in a species' reconstructed data ONLY when checklists were fit for it, i.e.
# script 07 was run with survey_set = "PC_ARU_CL".
honeycomb_structured_types <- c("Point_Count", "ARU")
honeycomb_checklist_types  <- c("Breeding Bird Atlas", "Linear transect")


# Prepare a species' survey sf for the honeycomb effort/count panels.
#
# The honeycomb sizes its centroid circles by the per-hexagon mean of the
# `count_per_effort` column (via summarize_surveys_by_hex() -> plot_honeycomb()).
# This helper sets that column to a count-per-survey-minute rate for EVERY
# retained survey:
#
#     count_per_effort = count / Survey_Duration_Minutes
#
# This is done unconditionally (in both PC_ARU and PC_ARU_CL runs) so that:
#   * fixed-effort point counts / ARUs and variable-effort checklists share one
#     comparable scale (e.g. a count of 1 from a 5-min point count -> 1/5; a
#     count of 1 from a 28-min checklist -> 1/28); and
#   * the honeycomb circles are directly comparable between the two model
#     versions when reviewing PC_ARU vs PC_ARU_CL products side by side.
# The per-hexagon mean of these rates is then taken exactly as before.
#
# Checklist surveys ("Breeding Bird Atlas", "Linear transect") are retained
# whenever they are present. They appear in sp_surveys only when checklists were
# fit for the species (a PC_ARU_CL run leaves them in dat_used$survey_counts);
# for a PC_ARU run there are no checklist rows, so the type filter reduces to
# Point_Count + ARU. Checklists are OBBA3-only, but both atlas periods are
# normalised together (the caller splits by Atlas afterwards) so the shared
# Atlas 2 / Atlas 3 circle scale stays comparable.
#
# NOTE: this rescaling affects ONLY the honeycomb circles. The relative-
# abundance raster panels are untouched and still represent expected count per
# 5-minute point count.
#
# Argument: sp_surveys  sf from load_species_surveys(); must carry Survey_Type,
#   count, and Survey_Duration_Minutes.
# Returns: list(
#   data       = filtered sf with a count_per_effort column (count per minute),
#   per_minute = logical, always TRUE (count_per_effort is a per-minute rate)
# )
prepare_honeycomb_surveys <- function(sp_surveys) {
  stopifnot(inherits(sp_surveys, "sf"))
  stopifnot(all(c("Survey_Type", "count") %in% names(sp_surveys)))

  # Keep every survey type the model used. For a PC_ARU run sp_surveys holds no
  # checklist rows, so this reduces to Point_Count + ARU.
  keep_types <- c(honeycomb_structured_types, honeycomb_checklist_types)

  out <- sp_surveys %>%
    dplyr::filter(Survey_Type %in% keep_types)

  if (!"Survey_Duration_Minutes" %in% names(out)) {
    stop("prepare_honeycomb_surveys(): Survey_Duration_Minutes is missing; ",
         "cannot compute count per survey-minute.", call. = FALSE)
  }
  if (any(!is.finite(out$Survey_Duration_Minutes) |
          out$Survey_Duration_Minutes <= 0)) {
    stop("prepare_honeycomb_surveys(): non-positive or missing ",
         "Survey_Duration_Minutes; cannot compute count per survey-minute.",
         call. = FALSE)
  }

  out <- out %>%
    dplyr::mutate(count_per_effort = count / Survey_Duration_Minutes)

  list(data = out, per_minute = TRUE)
}


# Panel subtitle wording for the honeycomb circle-size legend. count_per_effort
# from prepare_honeycomb_surveys() is always a per-minute rate; the per_minute
# argument is retained so callers that size circles by raw counts elsewhere can
# still request the plain-count wording.
honeycomb_size_label <- function(per_minute = TRUE) {
  if (isTRUE(per_minute)) {
    "observed mean count per survey-minute in hexagon"
  } else {
    "observed mean count in hexagon"
  }
}


# Paired summaries produced by 07b must be keyed on BCR. If a legacy file keyed
# on Biol_Region is found, fail loudly rather than silently mis-joining: the
# posterior change summaries are computed per region upstream and cannot be
# retro-fitted spatially here.
require_bcr_key <- function(x, what, sp_english) {
  if ("BCR" %in% names(x)) {
    x$BCR <- as.character(x$BCR)
    return(x)
  }
  if ("Biol_Region" %in% names(x)) {
    stop(
      "Paired summary '", what, "' for ", sp_english, " is keyed on ",
      "Biol_Region (atlas biological regions). Rerun ",
      "07b_fit_models_to_shared_footprint.R with BCR regions so that the ",
      "shared change summaries are computed per BCR.",
      call. = FALSE
    )
  }
  stop("Paired summary '", what, "' for ", sp_english,
       " has no BCR column.", call. = FALSE)
}


# Build the Atlas 2 + Atlas 3 assessment panels for one region.
#
# Rasters are pre-cropped/masked to the region and survey data are pre-filtered
# so that the shared display limits are computed on exactly the pixels and
# points each panel will show. For the province the region is the full study
# boundary, so cropping is a no-op on values.
#
# Argument: transform  Colour-scale transform for the relative-abundance
#                      (raster) panel, passed through to assess_region() and on
#                      to plot_raster_gg(). One of `"identity"`, `"sqrt"`,
#                      `"log"`, etc. Applies to both Atlas 2 and Atlas 3 panels
#                      so their scales stay comparable.
#
# Returns: list(A2 = <assess_region output>, A3 = <...>, scale_limits = <...>)
build_period_assessments <- function(region_geom,
                                     region_hex,
                                     r2, r3,
                                     A2_dat, A3_dat,
                                     region_boundaries,
                                     sp_english,
                                     title_suffix = "",
                                     colpal,
                                     rast_max_q    = 0.99,
                                     count_max_q   = 0.99,
                                     max_surveys_q = 0.80,
                                     transform     = "identity",
                                     count_size_label = "observed mean count in hexagon") {
  
  region_vect <- terra::vect(sf::st_transform(region_geom, terra::crs(r2)))
  
  r2_reg <- r2 |> terra::crop(region_vect) |> terra::mask(region_vect)
  r3_reg <- r3 |> terra::crop(region_vect) |> terra::mask(region_vect)
  
  # st_filter requires matching CRS. The precomputed region geometry is in the
  # study_boundary / prediction-grid CRS, but the per-species survey data may be
  # stored in a different CRS, so align the region to each data layer's CRS.
  A2_reg <- sf::st_filter(
    A2_dat, sf::st_transform(region_geom, sf::st_crs(A2_dat)),
    .predicate = sf::st_intersects
  )
  A3_reg <- sf::st_filter(
    A3_dat, sf::st_transform(region_geom, sf::st_crs(A3_dat)),
    .predicate = sf::st_intersects
  )
  
  # Hex summaries first so display limits reflect both atlas periods together.
  A2_hex <- summarize_hex(dat = A2_reg, hex_grid = region_hex,
                          rast = r2_reg, rast_max_q = rast_max_q)
  A3_hex <- summarize_hex(dat = A3_reg, hex_grid = region_hex,
                          rast = r3_reg, rast_max_q = rast_max_q)
  
  lims <- compute_assessment_scale_limits(
    hex_summaries = list(A2_hex, A3_hex),
    rasts         = list(r2_reg, r3_reg),
    rast_max_q    = rast_max_q,
    count_max_q   = count_max_q,
    max_surveys_q = max_surveys_q
  )
  
  # Shared lower bound for the relative-abundance colour ramp: the faintest
  # present (non-absent) density across BOTH periods. Absent pixels are already
  # NA (values at or below the absence floor set upstream), so this is the lower
  # edge of 'present'. Sharing it keeps the two period panels on an identical
  # scale and lets the pale low end of the palette be used.
  present_min <- function(r) {
    v <- terra::values(r, mat = FALSE)
    v <- v[is.finite(v) & v > 0]
    if (length(v)) min(v) else NA_real_
  }
  rast_min_shared <- suppressWarnings(
    min(c(present_min(r2_reg), present_min(r3_reg)), na.rm = TRUE)
  )
  if (!is.finite(rast_min_shared)) rast_min_shared <- NULL
  
  title_for <- function(period) {
    if (nzchar(title_suffix)) {
      paste0(sp_english, " - ", period, " - ", title_suffix)
    } else {
      paste0(sp_english, " - ", period)
    }
  }
  
  assess_one <- function(period, dat, rast) {
    assess_region(
      region               = region_geom,
      region_boundaries    = region_boundaries,
      sp_dat               = dat,
      rast                 = rast,
      hex_grid             = region_hex,
      water                = NULL,
      rast_max             = lims$rast_max,
      rast_min             = rast_min_shared,
      max_count_per_effort = lims$max_count_per_effort,
      max_surveys          = lims$max_surveys,
      transform            = transform,
      title                = title_for(period),
      model_source         = NULL,
      data_source          = NULL,
      color_palette        = colpal,
      count_size_label     = count_size_label
    )
  }
  
  list(
    A2           = assess_one("Atlas 2", A2_reg, r2_reg),
    A3           = assess_one("Atlas 3", A3_reg, r3_reg),
    scale_limits = lims
  )
}


# Build the Atlas 2 + Atlas 3 *honeycomb-only* panels for one region.
#
# Used for species without a fitted model: there is no prediction raster, so the
# observed survey hex summary is rendered on its own (effort + observed mean
# count) via plot_honeycomb(). Atlas 2 and Atlas 3 share effort/count scaling so
# the two periods are directly comparable.
#
# Returns: list(A2 = <ggplot>, A3 = <ggplot>)
build_period_honeycombs <- function(region_geom,
                                    region_hex,
                                    region_bounds,
                                    A2_dat, A3_dat,
                                    sp_english,
                                    title_suffix  = "",
                                    max_surveys_q = 0.80,
                                    count_max_q   = 0.99,
                                    count_size_label = "observed mean count in hexagon") {
  
  A2_reg <- sf::st_filter(
    A2_dat, sf::st_transform(region_geom, sf::st_crs(A2_dat)),
    .predicate = sf::st_intersects
  )
  A3_reg <- sf::st_filter(
    A3_dat, sf::st_transform(region_geom, sf::st_crs(A3_dat)),
    .predicate = sf::st_intersects
  )
  
  # Survey-only hex summaries (no raster involved).
  A2_hex <- summarize_surveys_by_hex(dat = A2_reg, hex_grid = region_hex)
  A3_hex <- summarize_surveys_by_hex(dat = A3_reg, hex_grid = region_hex)
  
  # Shared scaling across both atlas periods. compute_assessment_scale_limits()
  # cannot be used here (it requires a pred_mean column when no raster is given),
  # so the two honeycomb limits are resolved directly with the same helper and
  # quantiles that plot_honeycomb / assess_region use.
  shared_max_surveys <- ceiling(resolve_positive_max(
    x = c(A2_hex$n_surveys, A3_hex$n_surveys),
    max_q = max_surveys_q, default = 1, argument_name = "max_surveys"
  ))
  shared_max_count <- resolve_positive_max(
    x = c(A2_hex$mean_count_per_effort, A3_hex$mean_count_per_effort),
    max_q = count_max_q, default = 1, argument_name = "max_count_per_effort"
  )
  
  honey_one <- function(hex_obs, period) {
    plot_honeycomb(
      hex_summary            = hex_obs,
      study_area             = region_geom,
      region_boundaries      = region_bounds,
      water                  = NULL,
      max_surveys            = shared_max_surveys,
      max_surveys_q          = max_surveys_q,
      max_count_per_effort   = shared_max_count,
      max_count_per_effort_q = count_max_q,
      count_size_label       = count_size_label
    ) +
      # Replace plot_honeycomb's generic title with the atlas period; its
      # descriptive subtitle (hex fill / circle size) is retained.
      ggplot2::labs(title = period)
  }
  
  list(
    A2 = honey_one(A2_hex, "Atlas 2"),
    A3 = honey_one(A3_hex, "Atlas 3")
  )
}


# Render a single ggplot/patchwork page to a PNG.
save_page <- function(plot, path,
                      width  = 20,
                      height = 10,
                      res    = 300) {
  ragg::agg_png(
    filename = path, width = width, height = height,
    units = "in", res = res, background = "white"
  )
  on.exit(grDevices::dev.off(), add = TRUE)
  print(plot)
  invisible(path)
}
