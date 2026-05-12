# ============================================================
# spatial_utils.R
#
# Spatial utilities used across the BBA workflow:
#   - CRS definition (projected; intended km units for INLA)
#   - raster crop/mask helpers (terra)
#   - spatial CV block utilities (grid blocks)
#   - covariate support helper (quantile-based, with optional zero spike)
#   - plotting helper for survey vs landscape hist overlays
# ============================================================

# ---------------------------
# CRS helper
# ---------------------------

# Define a projected CRS intended to have kilometre units to keep coordinate magnitudes small (INLA).
# Note: units=km is not universally supported in every GIS stack; we keep this because the pipeline
# is explicitly built around it. Downstream code should assume "distance units = km" when using this CRS.

#' Get the default projected CRS used throughout the workflow
#'
#' This helper returns a Lambert Conformal Conic (LCC) CRS intended for
#' Ontario-scale spatial modelling. The PROJ string specifies `+units=km`
#' so coordinates are in kilometres (helpful for INLA/SPDE meshes where
#' very large coordinate magnitudes can hurt numerical stability).
#'
#' @return An `sf::crs` object.
get_aea_km_crs <- function() {
  # Use a single-line PROJ string to avoid newline/whitespace surprises.
  proj <- paste(
    "+proj=lcc",
    "+lat_1=49",
    "+lat_2=77",
    "+lat_0=49",
    "+lon_0=-95",
    "+datum=NAD83",
    "+units=km",
    "+no_defs"
  )
  sf::st_crs(proj)
}

# ---------------------------
# Raster helpers (terra)
# ---------------------------

# Crop and mask a SpatRaster to an sf boundary.

#' Crop and mask a raster to a study boundary
#'
#' Convenience wrapper around `terra::crop()` and `terra::mask()` that:
#' 1) checks CRS existence on the raster,
#' 2) transforms the boundary to the raster CRS, and
#' 3) crops then masks the raster so cells outside the boundary become NA.
#'
#' @param r A `terra::SpatRaster`.
#' @param boundary_sf An `sf` or `sfc` boundary polygon.
#' @return A `terra::SpatRaster` cropped and masked to `boundary_sf`.
crop_mask_to_boundary <- function(r, boundary_sf) {
  stopifnot(inherits(r, "SpatRaster"))
  stopifnot(inherits(boundary_sf, "sf") || inherits(boundary_sf, "sfc"))
  
  # Get CRS as a single string
  r_crs <- terra::crs(r)
  
  if (is.null(r_crs) || is.na(r_crs) || !nzchar(r_crs)) {
    stop("Raster CRS is missing; cannot crop/mask safely.")
  }
  
  # Transform boundary to raster CRS
  boundary_sf <- sf::st_transform(boundary_sf, sf::st_crs(r_crs))
  b <- terra::vect(boundary_sf)
  
  r <- terra::crop(r, b)
  terra::mask(r, b)
}

# Make a template raster (resolution in metres) and mask it to a boundary.

#' Create a template raster at a target resolution (metres) for rasterization/warping
#'
#' Builds an empty `terra::SpatRaster` covering the bounding box of the
#' boundary (in `crs_wkt`) at resolution `res_m`, then masks it to the
#' boundary polygon. This template is used as a common grid for aligning
#' covariate rasters (road buffers, insects, water, etc.).
#'
#' @param boundary_sf `sf`/`sfc` boundary polygon.
#' @param res_m Numeric scalar; target resolution in metres.
#' @param crs_wkt Character WKT/PROJ string for the template CRS.
#' @return A masked `terra::SpatRaster` template.
make_template_raster_m <- function(boundary_sf, res_m, crs_wkt) {
  stopifnot(inherits(boundary_sf, "sf") || inherits(boundary_sf, "sfc"))
  stopifnot(is.numeric(res_m), length(res_m) == 1, res_m > 0)
  stopifnot(is.character(crs_wkt), length(crs_wkt) == 1, nzchar(crs_wkt))
  
  b <- sf::st_transform(boundary_sf, sf::st_crs(crs_wkt))
  bb <- sf::st_bbox(b)
  
  r <- terra::rast(
    xmin = bb["xmin"], xmax = bb["xmax"],
    ymin = bb["ymin"], ymax = bb["ymax"],
    resolution = res_m,
    crs = crs_wkt
  )
  
  # IMPORTANT: allocate values so terra::mask() works
  terra::values(r) <- 0
  
  crop_mask_to_boundary(r, b)
}

# ------------------------------------------------------------
# Covariate support by BCR
# ------------------------------------------------------------

# Quantile-based "support" score: support = 1 - 2*abs(q - 0.5), where q is the
# mid-rank ECDF quantile of the grid value within survey values.
#
# If the survey distribution has a notable point-mass at zero (>= zero_spike_min),
# treat zeros specially so a grid value of zero gets a quantile near the middle of
# the zero mass (p0/2), and nonzeros are ranked within the nonzero distribution then
# shifted by p0.

#' Quantile-based covariate support diagnostic by BCR
#'
#' Computes how well survey covariate values cover the landscape covariate
#' distribution, stratified by Bird Conservation Region (BCR). This is used
#' to identify extrapolation risk: covariate values common on the landscape
#' but rarely sampled by surveys.
#'
#' The implementation allows a 'zero spike' option for covariates with many
#' exact zeros (e.g., proportion of habitat types), where the zero mass is
#' treated separately from the continuous part.
#'
#' @param survey_df Data frame of survey covariates with a BCR column.
#' @param grid_df Data frame of landscape covariates (prediction grid) with BCR.
#' @param covars Character vector of covariate names to evaluate.
#' @param bcr_col Name of the BCR column.
#' @param probs Quantiles used to define support bins.
#' @param zero_spike_covars Optional vector of covariates to treat as zero-inflated.
#' @return A tidy table with support metrics per covariate × BCR.
covariate_support_by_bcr <- function(
    grid_sf,
    survey_sf,
    covariate,
    bcr_col = "BCR",
    zero_value = 0,
    zero_spike_min = 0.05,
    return_q = TRUE,
    return_n = TRUE
) {
  stopifnot(inherits(grid_sf, "sf"), inherits(survey_sf, "sf"))
  stopifnot(is.character(covariate), length(covariate) == 1)
  stopifnot(covariate %in% names(grid_sf), covariate %in% names(survey_sf))
  stopifnot(bcr_col %in% names(grid_sf), bcr_col %in% names(survey_sf))
  
  grid_vals <- grid_sf[[covariate]]
  grid_bcr  <- grid_sf[[bcr_col]]
  surv_vals <- survey_sf[[covariate]]
  surv_bcr  <- survey_sf[[bcr_col]]
  
  support <- rep(NA_real_, length(grid_vals))
  q_all   <- rep(NA_real_, length(grid_vals))
  n_used  <- rep(NA_integer_, length(grid_vals))
  
  for (bcr in sort(unique(stats::na.omit(grid_bcr)))) {
    
    grid_idx <- which(grid_bcr == bcr & !is.na(grid_vals))
    if (!length(grid_idx)) next
    
    s_vals <- surv_vals[surv_bcr == bcr]
    s_vals <- s_vals[!is.na(s_vals)]
    n <- length(s_vals)
    if (!n) next
    
    n_used[grid_idx] <- n
    g_vals <- grid_vals[grid_idx]
    
    p0 <- mean(s_vals == zero_value)
    
    if (p0 >= zero_spike_min) {
      s_nonzero <- s_vals[s_vals != zero_value]
      q <- rep(NA_real_, length(g_vals))
      
      is_zero <- (g_vals == zero_value)
      q[is_zero] <- p0 / 2
      
      is_nz <- !is_zero
      if (any(is_nz)) {
        v <- g_vals[is_nz]
        
        if (!length(s_nonzero)) {
          q[is_nz] <- 1
        } else {
          u <- sort(unique(s_nonzero))
          counts <- tabulate(match(s_nonzero, u))
          cum <- cumsum(counts)
          nn <- length(s_nonzero)
          
          pos_le <- findInterval(v, u, rightmost.closed = TRUE)
          pos_lt <- findInterval(v, u, left.open = TRUE, rightmost.closed = TRUE)
          
          F_le <- ifelse(pos_le <= 0, 0, cum[pos_le] / nn)
          F_lt <- ifelse(pos_lt <= 0, 0, cum[pos_lt] / nn)
          
          F_mid <- (F_lt + F_le) / 2
          q[is_nz] <- p0 + (1 - p0) * F_mid
        }
      }
      
    } else {
      u <- sort(unique(s_vals))
      counts <- tabulate(match(s_vals, u))
      cum <- cumsum(counts)
      nn <- length(s_vals)
      
      pos_le <- findInterval(g_vals, u, rightmost.closed = TRUE)
      pos_lt <- findInterval(g_vals, u, left.open = TRUE, rightmost.closed = TRUE)
      
      F_le <- ifelse(pos_le <= 0, 0, cum[pos_le] / nn)
      F_lt <- ifelse(pos_lt <= 0, 0, cum[pos_lt] / nn)
      
      q <- (F_lt + F_le) / 2
    }
    
    q_all[grid_idx] <- q
    support[grid_idx] <- 1 - 2 * abs(q - 0.5)
  }
  
  if (return_q || return_n) {
    out <- list(support = support)
    if (return_q) out$q <- q_all
    if (return_n) out$n_surveys <- n_used
    return(out)
  }
  
  support
}

# ------------------------------------------------------------
# Plotting helper
# ------------------------------------------------------------

#' Plot survey vs landscape covariate distributions (overlayed histograms)
#'
#' Convenience plotting function used in QC reports to compare:
#' - survey covariate distribution (by atlas / survey type), and
#' - landscape covariate distribution (prediction grid)
#' optionally stratified by BCR.
#'
#' @param survey_df Survey covariate table.
#' @param grid_df Landscape covariate table.
#' @param covar Name of the covariate to plot.
#' @param group_cols Columns defining survey subgroups to overlay (e.g., Atlas, Survey_Type).
#' @param bcr Optional BCR level or column name to facet by.
#' @return A `ggplot2` object.
plot_covariate_hist_overlay <- function(
    grid_sf,
    survey_sf,
    covariate,
    region_col = "region",
    bins = 30,
    binwidth = NULL,
    origin = NULL,
    survey_width_frac = 0.5,
    landscape_fill = "black",
    survey_fill = "red",
    facet_scales = "free_y",
    title = NULL
) {
  # Requires: dplyr, ggplot2, sf
  stopifnot(
    covariate %in% names(grid_sf),
    covariate %in% names(survey_sf),
    region_col %in% names(grid_sf),
    region_col %in% names(survey_sf)
  )
  
  grid_df <- sf::st_drop_geometry(grid_sf)[, c(region_col, covariate), drop = FALSE]
  names(grid_df) <- c("region", "value")
  grid_df$source <- "Landscape"
  
  svy_df <- sf::st_drop_geometry(survey_sf)[, c(region_col, covariate), drop = FALSE]
  names(svy_df) <- c("region", "value")
  svy_df$source <- "Survey"
  
  df <- dplyr::bind_rows(grid_df, svy_df) |>
    dplyr::filter(!is.na(.data$region), !is.na(.data$value))
  
  ann_svy_n <- df |>
    dplyr::filter(.data$source == "Survey") |>
    dplyr::count(.data$region, name = "n_surveys") |>
    dplyr::mutate(label = paste0("n = ", .data$n_surveys))
  
  # weights so bars sum to 1 within each (region, source)
  df <- df |>
    dplyr::group_by(.data$region, .data$source) |>
    dplyr::mutate(w = 1 / dplyr::n()) |>
    dplyr::ungroup()
  
  vmin <- min(df$value, na.rm = TRUE)
  vmax <- max(df$value, na.rm = TRUE)
  
  if (is.null(binwidth)) {
    rng <- vmax - vmin
    binwidth <- if (is.finite(rng) && rng > 0) rng / bins else 1
  }
  
  if (is.null(origin)) {
    origin <- floor(vmin / binwidth) * binwidth
  }
  
  breaks_seq <- seq(origin, vmax + binwidth, by = binwidth)
  
  df_hist <- df |>
    dplyr::mutate(
      bin_id = findInterval(.data$value, breaks_seq, rightmost.closed = TRUE),
      bin_id = pmax(1L, .data$bin_id),
      bin_left = breaks_seq[.data$bin_id],
      bin_center = .data$bin_left + binwidth / 2
    ) |>
    dplyr::group_by(.data$region, .data$source, .data$bin_center) |>
    dplyr::summarise(prop = sum(.data$w), .groups = "drop")
  
  df_land <- df_hist |> dplyr::filter(.data$source == "Landscape")
  df_svy  <- df_hist |> dplyr::filter(.data$source == "Survey")
  
  if (is.null(title)) {
    title <- paste0("Survey vs Landscape distribution: ", covariate)
  }
  
  ggplot2::ggplot() +
    ggplot2::geom_col(
      data = df_land,
      ggplot2::aes(x = .data$bin_center, y = .data$prop, fill = "Landscape"),
      width = binwidth
    ) +
    ggplot2::geom_col(
      data = df_svy,
      ggplot2::aes(x = .data$bin_center, y = .data$prop, fill = "Survey"),
      width = binwidth * survey_width_frac
    ) +
    ggplot2::facet_wrap(stats::as.formula("~region"), scales = facet_scales) +
    ggplot2::scale_fill_manual(
      values = c("Landscape" = landscape_fill, "Survey" = survey_fill),
      name = NULL
    ) +
    ggplot2::labs(
      title = title,
      x = covariate,
      y = "Proportion of surveys or pixels"
    ) +
    ggplot2::theme_bw() +
    ggplot2::geom_text(
      data = ann_svy_n,
      ggplot2::aes(x = -Inf, y = Inf, label = .data$label),
      inherit.aes = FALSE,
      hjust = -0.1,
      vjust = 1.2,
      size = 4
    )
}
