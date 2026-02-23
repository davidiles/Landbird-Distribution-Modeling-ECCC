
# Define the CRS; note kilometre units to ensure smaller coordinate values for INLA
get_aea_km_crs <- function() {
  sf::st_crs(
    "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=49 +lon_0=-95 
     +datum=NAD83 +units=km +no_defs"
  )
}

crop_mask_to_boundary <- function(r, boundary_sf) {
  # Transform boundary to raster CRS
  boundary_sf <- sf::st_transform(boundary_sf, sf::st_crs(terra::crs(r)))
  b <- terra::vect(boundary_sf)
  r <- terra::crop(r, b)
  terra::mask(r, b)
}

# Make a template raster in metre units
make_template_raster_m <- function(boundary_sf, res_m, crs_wkt) {
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

# Grid-based spatial blocks in projected CRS with km units
# Returns a character vector block_id of length nrow(surveys_sf)
make_spatial_blocks_grid <- function(surveys_sf, block_size_km) {
  stopifnot(inherits(surveys_sf, "sf"))
  stopifnot(is.numeric(block_size_km), length(block_size_km) == 1, block_size_km > 0)
  
  xy <- sf::st_coordinates(surveys_sf)
  if (nrow(xy) != nrow(surveys_sf)) stop("st_coordinates() did not return one row per survey.")
  
  bx <- floor(xy[, 1] / block_size_km)
  by <- floor(xy[, 2] / block_size_km)
  
  # stable, readable id
  paste0("B_", bx, "_", by)
}

# Given block_ids, create a fold plan (repeat x block -> fold)
# Returns a data.frame with columns: repeat, block_id, fold
make_block_cv_plan <- function(block_ids, n_folds = 5, n_repeats = 5, seed = 1) {
  stopifnot(is.vector(block_ids), length(block_ids) > 0)
  stopifnot(is.numeric(n_folds), length(n_folds) == 1, n_folds >= 2)
  stopifnot(is.numeric(n_repeats), length(n_repeats) == 1, n_repeats >= 1)
  
  blocks <- sort(unique(as.character(block_ids)))
  n_blocks <- length(blocks)
  if (n_blocks < n_folds) stop("Number of unique blocks (", n_blocks, ") < n_folds (", n_folds, ").")
  
  out <- vector("list", n_repeats)
  
  for (r in seq_len(n_repeats)) {
    set.seed(seed + r - 1)
    perm <- sample(blocks, size = n_blocks, replace = FALSE)
    
    fold_id <- rep(seq_len(n_folds), length.out = n_blocks)
    
    df_r <- data.frame(
      rep = r,
      block_id = perm,
      fold = fold_id,
      stringsAsFactors = FALSE
    )
    
    out[[r]] <- df_r
  }
  
  do.call(rbind, out)
}

# ------------------------------------------------------------
# Convert block_id strings ("B_bx_by") to square polygons
# ------------------------------------------------------------

# Returns an sf with columns: block_id, geometry (POLYGON)
# Assumes your survey coordinates are in km units (as your CRS helper suggests).
make_block_polygons_grid <- function(block_ids, block_size_km, crs) {
  stopifnot(is.vector(block_ids), length(block_ids) > 0)
  stopifnot(is.numeric(block_size_km), length(block_size_km) == 1, block_size_km > 0)
  stopifnot(!is.null(crs))
  
  ids <- sort(unique(as.character(block_ids)))
  
  # Parse bx/by from "B_<bx>_<by>"
  # Works with negative indices too (e.g., B_-12_7)
  parse_one <- function(id) {
    s <- sub("^B_", "", id)
    parts <- strsplit(s, "_", fixed = TRUE)[[1]]
    if (length(parts) != 2) stop("Unexpected block_id format: ", id)
    bx <- as.integer(parts[1])
    by <- as.integer(parts[2])
    c(bx = bx, by = by)
  }
  
  parsed <- t(vapply(ids, parse_one, numeric(2)))
  bx <- parsed[, "bx"]
  by <- parsed[, "by"]
  
  # Build square polygons: [x0,x1] x [y0,y1]
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
  stopifnot(covariate %in% names(grid_sf),
            covariate %in% names(survey_sf),
            region_col %in% names(grid_sf),
            region_col %in% names(survey_sf))
  
  # Pull values (no geometry)
  grid_df <- sf::st_drop_geometry(grid_sf)[, c(region_col, covariate), drop = FALSE]
  names(grid_df) <- c("region", "value")
  grid_df$source <- "Landscape"
  
  svy_df <- sf::st_drop_geometry(survey_sf)[, c(region_col, covariate), drop = FALSE]
  names(svy_df) <- c("region", "value")
  svy_df$source <- "Survey"
  
  df <- dplyr::bind_rows(grid_df, svy_df) |>
    dplyr::filter(!is.na(.data$region), !is.na(.data$value))
  
  # --- NEW: region-specific survey counts (after filtering NAs in region/value) ---
  ann_svy_n <- df |>
    dplyr::filter(.data$source == "Survey") |>
    dplyr::count(.data$region, name = "n_surveys") |>
    dplyr::mutate(label = paste0("n = ", .data$n_surveys))
  # ------------------------------------------------------------------------------
  
  # Weights so bars sum to 1 within each (region, source)
  df <- df |>
    dplyr::group_by(.data$region, .data$source) |>
    dplyr::mutate(w = 1 / dplyr::n()) |>
    dplyr::ungroup()
  
  # Choose binning
  vmin <- min(df$value, na.rm = TRUE)
  vmax <- max(df$value, na.rm = TRUE)
  
  if (is.null(binwidth)) {
    rng <- vmax - vmin
    binwidth <- if (is.finite(rng) && rng > 0) rng / bins else 1
  }
  
  if (is.null(origin)) {
    origin <- floor(vmin / binwidth) * binwidth
  }
  
  # Build breaks
  breaks_seq <- seq(origin, vmax + binwidth, by = binwidth)
  
  # Bin
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
    
    # Landscape (wide background bars)
    ggplot2::geom_col(
      data = df_land,
      ggplot2::aes(x = .data$bin_center, y = .data$prop),
      width = binwidth,
      fill = landscape_fill
    ) +
    # Survey (narrow foreground bars, centered)
    ggplot2::geom_col(
      data = df_svy,
      ggplot2::aes(x = .data$bin_center, y = .data$prop),
      width = binwidth * survey_width_frac,
      fill = survey_fill
    ) +
    ggplot2::facet_wrap(stats::as.formula("~region"), scales = facet_scales) +
    ggplot2::labs(
      title = title,
      x = covariate,
      y = "Proportion of surveys or pixels"
    ) +
    ggplot2::theme_bw() +
    
    # dummy fill legend entry that is invisible, but reserves space
    ggplot2::geom_rect(
      data = data.frame(xmin=0, xmax=0, ymin=0, ymax=0, dummy=" "),
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = dummy),
      inherit.aes = FALSE,
      alpha = 0,
      show.legend = TRUE
    ) +
    ggplot2::scale_fill_manual(
      values = c("Landscape" = "black", "Survey" = "red", " " = "white"),
      breaks = c("Landscape", "Survey", " "),
      name = NULL
    ) +
    ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(alpha = 1))) +
    
    # --- NEW: per-region n annotation (matches facets by 'region') ---
    ggplot2::geom_text(
      data = ann_svy_n,
      ggplot2::aes(x = -Inf, y = Inf, label = .data$label),
      inherit.aes = FALSE,
      hjust = -0.1,
      vjust = 1.2,
      size = 4
    )
}
