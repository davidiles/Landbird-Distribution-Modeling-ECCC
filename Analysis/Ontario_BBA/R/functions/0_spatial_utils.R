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
# Note: units=km is not universally supported in every GIS stack; we keep this because your pipeline
# is explicitly built around it. Downstream code should assume "distance units = km" when using this CRS.
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
crop_mask_to_boundary <- function(r, boundary_sf) {
  stopifnot(inherits(r, "SpatRaster"))
  stopifnot(inherits(boundary_sf, "sf") || inherits(boundary_sf, "sfc"))
  
  r_crs <- terra::crs(r, describe = TRUE)
  if (is.na(r_crs) || !nzchar(r_crs)) stop("Raster CRS is missing; cannot crop/mask safely.")
  
  # Transform boundary to raster CRS
  boundary_sf <- sf::st_transform(boundary_sf, sf::st_crs(terra::crs(r)))
  b <- terra::vect(boundary_sf)
  
  r <- terra::crop(r, b)
  terra::mask(r, b)
}

# Make a template raster (resolution in metres) and mask it to a boundary.
# NOTE: despite the file comment, this function is CRS-agnostic: you supply crs_wkt and res_m.
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
# Spatial CV helpers
# ------------------------------------------------------------

# Grid-based spatial blocks in projected CRS with km units.
# Returns a character vector block_id of length nrow(surveys_sf).
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
      repeat = r,
      block_id = perm,
      fold = fold_id,
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
make_block_polygons_grid <- function(block_ids, block_size_km, crs) {
  stopifnot(is.vector(block_ids), length(block_ids) > 0)
  stopifnot(is.numeric(block_size_km), length(block_size_km) == 1, block_size_km > 0)
  stopifnot(!is.null(crs))
  
  ids <- sort(unique(as.character(block_ids)))
  
  # Parse bx/by from "B_<bx>_<by>" (supports negatives: B_-12_7)
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