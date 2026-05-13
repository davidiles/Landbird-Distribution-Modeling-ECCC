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
# ============================================================

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
# Relative abundance maps
# ------------------------------------------------------------

# Prepare multiple relative-abundance rasters for side-by-side plotting.
#
# Purpose:
#   Applies a shared absence threshold and a shared upper colour-scale limit
#   across one or more rasters, so Atlas 2 and Atlas 3 maps use comparable
#   legends. Low values are treated as absent, and high values are clamped at
#   a common quantile-based maximum.
#
# Arguments:
#   ...                 - named terra SpatRaster objects.
#   rast_absent_limit   - values at or below this threshold are set to NA.
#   rast_max_quantile   - quantile used to estimate the upper plotting limit.
#
# Returns:
#   A list containing processed rasters, shared z limits, legend breaks, and
#   original per-raster upper quantiles.
prepare_relative_abundance_rasters <- function(...,
                                               rast_absent_limit = 1/250,
                                               rast_max_quantile = 0.99) {
  
  rasts <- list(...)
  
  if (length(rasts) == 0) {
    stop("At least one raster must be supplied.")
  }
  
  if (!all(vapply(rasts, inherits, logical(1), what = "SpatRaster"))) {
    stop("All inputs passed through `...` must be terra SpatRaster objects.")
  }
  
  process_raster <- function(rast) {
    
    r <- rast
    
    # Set low values to NA
    r[r <= rast_absent_limit] <- NA
    
    # Calculate upper quantile after low values are removed
    zmax <- as.numeric(stats::quantile(
      terra::values(r),
      probs = rast_max_quantile,
      na.rm = TRUE
    ))
    
    list(
      rast = r,
      zmax = zmax
    )
  }
  
  processed <- lapply(rasts, process_raster)
  
  zmax_values <- vapply(processed, `[[`, numeric(1), "zmax")
  
  # Shared upper limit across all rasters
  zmax_shared <- min(zmax_values, na.rm = TRUE)
  
  if (!is.finite(zmax_shared) || zmax_shared <= 0) {
    zmax_shared <- 1
  }
  
  # Clamp all rasters to the same upper limit
  rasts_clamped <- lapply(processed, function(x) {
    terra::clamp(
      x$rast,
      upper = zmax_shared,
      values = TRUE
    )
  })
  
  names(rasts_clamped) <- names(rasts)
  
  zlim <- c(rast_absent_limit, zmax_shared)
  
  # Fixed legend breaks shared by all maps
  zbreaks <- seq(zlim[1], zlim[2], length.out = 5)
  
  list(
    rasters = rasts_clamped,
    zlim = zlim,
    zbreaks = zbreaks,
    zmax_original = zmax_values
  )
}


# Create map of posterior mean relative abundance for one atlas period
make_relabund_map <- function(species_name,
                              atlas_label,
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
    "Relative Abundance - ", atlas_label,
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
      name     = "Expected count\nper 5-min\npoint count",
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

summarize_polygon_hex_draw_change <- function(hex_draws,
                                              hex_grid,
                                              polygon,
                                              ci_level = 0.90,
                                              baseline_min = 1e-5,
                                              return_draws = FALSE) {
  
  required_cols <- c("hex_id", "mu_OBBA2", "mu_OBBA3")
  missing_cols <- setdiff(required_cols, names(hex_draws))
  
  if (length(missing_cols) > 0) {
    stop(
      "hex_draws is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  stopifnot(inherits(hex_grid, "sf"))
  stopifnot(inherits(polygon, "sf"))
  stopifnot("hex_id" %in% names(hex_grid))
  
  if (!is.numeric(ci_level) || length(ci_level) != 1 ||
      ci_level <= 0 || ci_level >= 1) {
    stop("ci_level must be a single number strictly between 0 and 1.")
  }
  
  alpha <- (1 - ci_level) / 2
  
  # ----------------------------------------------------------
  # Calculate fraction of each hex inside polygon
  # ----------------------------------------------------------
  
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
  
  hex_weights <- hex_intersection |>
    dplyr::mutate(intersection_area = sf::st_area(geometry)) |>
    sf::st_drop_geometry() |>
    dplyr::group_by(hex_id) |>
    dplyr::summarise(
      intersection_area = sum(intersection_area),
      .groups = "drop"
    ) |>
    dplyr::left_join(
      tibble::tibble(
        hex_id = hex_grid$hex_id,
        hex_area = hex_area
      ),
      by = "hex_id"
    ) |>
    dplyr::mutate(
      area_weight = as.numeric(intersection_area / hex_area),
      area_weight = pmin(pmax(area_weight, 0), 1)
    ) |>
    dplyr::select(hex_id, area_weight)
  
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
#   4. Build a three-panel patchwork figure:
#      - Panel A: raw raster of model predictions (`plot_raster_gg`).
#      - Panel B: hex-filled by mean predicted abundance, circles by observed
#        mean count (`plot_hex_pred_obs`).
#      - Panel C: honeycomb effort/detection map (`plot_honeycomb`).
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
#     \item{`plot_combined`}{A `patchwork` ggplot (three panels).}
#     \item{`hex_summary`}{`sf` data frame with one row per hexagon
#       containing columns from `summarize_hex()`.}
#     \item{`region_stats`}{One-row `data.frame` from `compute_region_stats()`.}
#   }
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

# Plot observed mean counts over hexagon-averaged model predictions
#
# Fills each surveyed hexagon by mean predicted relative abundance and overlays
# circles scaled by observed mean count per effort. Hexagons with a
# `pred_mean` of 0 or NA are rendered as transparent (absent); all others are
# mapped onto the colour ramp. The raster is used only to derive the upper cap
# for the shared colour scale (`rast_max_q`).
#
# Argument: hex_summary          Output of `summarize_hex()`. Must contain
#   `pred_mean`, `n_surveys`, and `mean_count_per_effort`.
# Argument: rast                 A `SpatRaster` used to compute the colour scale
#   upper limit. Not re-extracted here.
# Argument: study_area           `sf` polygon.
# Argument: water                Optional `processed_water` object.
# Argument: rast_max_q           Quantile (0-1) used to cap the colour scale
#   upper end. Default `0.99`.
# Argument: max_count_per_effort Upper bound for observed circle scaling (auto).
# Argument: palette              Fill palette; default matches `plot_raster_gg()`.
# Argument: circle_fill          Fill colour for observed-count circles.
# Argument: water_fill           Fill colour for water polygons.
# Argument: transform            Fill scale transform.
# Returns: A list: `plot` (ggplot), `hex_summary` (sf), `corr_unweighted` (dbl).
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
