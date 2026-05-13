# ============================================================
# model_product_utils.R
# ============================================================

library(sf)
library(dplyr)
library(ggplot2)
library(terra)
library(patchwork)
library(magick)

# Rasterize an sf grid
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


#' Create map of posterior mean relative abundance for one atlas period
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


#' # ------------------------------------------------------------
#' # Small helpers
#' # ------------------------------------------------------------
#' 
#' # compute abundance limits for plotting
#' compute_shared_zmax <- function(preds, prob = 0.99, fallback = 1) {
#'   zmax2 <- as.numeric(stats::quantile(preds$OBBA2$OBBA2_q50, prob, na.rm = TRUE))
#'   zmax3 <- as.numeric(stats::quantile(preds$OBBA3$OBBA3_q50, prob, na.rm = TRUE))
#'   zmax <- max(zmax2, zmax3, na.rm = TRUE)
#'   if (!is.finite(zmax) || zmax <= 0) zmax <- fallback
#'   zmax
#' }
#' 
#' # build square overlays 
#' build_atlas_square_overlays <- function(sp_square_summary, atlas_sq_centroids_all) {
#'   if (is.null(sp_square_summary)) {
#'     return(list(OBBA2 = NULL, OBBA3 = NULL))
#'   }
#'   
#'   sq2 <- dplyr::filter(sp_square_summary, Atlas == "OBBA2")
#'   sq3 <- dplyr::filter(sp_square_summary, Atlas == "OBBA3")
#'   
#'   list(
#'     OBBA2 = dplyr::left_join(atlas_sq_centroids_all, sq2, by = "square_id"),
#'     OBBA3 = dplyr::left_join(atlas_sq_centroids_all, sq3, by = "square_id")
#'   )
#' }
#' 
#' 
#' #' Wrap a species label for nicer plot titles
#' #'
#' #' @param label Character species name.
#' #' @param max_length Maximum characters before wrapping.
#' #' @return Character string, optionally with "<br>" inserted.
#' wrap_species_label <- function(label, max_length = 18) {
#'   if (is.na(label) || nchar(label) <= max_length) return(label)
#'   
#'   words <- strsplit(label, " ")[[1]]
#'   if (length(words) <= 1) return(label)
#'   
#'   paste0(
#'     paste(words[-length(words)], collapse = " "),
#'     "<br>",
#'     words[length(words)]
#'   )
#' }
#' 
#' #' Expand a bounding box by a fraction
#' #'
#' #' @param bbox An sf bbox object.
#' #' @param frac Fractional expansion.
#' #' @return Named numeric vector xmin/ymin/xmax/ymax.
#' expand_bbox <- function(bbox, frac = 0.10) {
#'   bbox <- as.numeric(bbox)
#'   names(bbox) <- c("xmin", "ymin", "xmax", "ymax")
#'   
#'   xrange <- bbox["xmax"] - bbox["xmin"]
#'   yrange <- bbox["ymax"] - bbox["ymin"]
#'   
#'   bbox["xmin"] <- bbox["xmin"] - frac * xrange
#'   bbox["xmax"] <- bbox["xmax"] + frac * xrange
#'   bbox["ymin"] <- bbox["ymin"] - frac * yrange
#'   bbox["ymax"] <- bbox["ymax"] + frac * yrange
#'   
#'   bbox
#' }
#' 
#' #' Get a padded plotting bbox from a study boundary
#' #'
#' #' @param study_boundary sf polygon.
#' #' @param frac Fractional expansion.
#' #' @return Named numeric vector xmin/ymin/xmax/ymax.
#' map_bbox <- function(study_boundary, frac = 0.10) {
#'   expand_bbox(sf::st_bbox(study_boundary), frac = frac)
#' }
#' 
#' #' Rasterize an sf grid to a stars object
#' #'
#' #' @param grid_sf sf object with geometry and a value field.
#' #' @param field Name of field to rasterize.
#' #' @param res Raster resolution in map units.
#' #' @return stars object.
#' rasterize_sf_to_stars <- function(grid_sf, field, res) {
#'   stopifnot(inherits(grid_sf, "sf"))
#'   stopifnot(field %in% names(grid_sf))
#'   
#'   v <- terra::vect(grid_sf)
#'   r_template <- terra::rast(v, res = res)
#'   r <- terra::rasterize(v, r_template, field = field, fun = mean, na.rm = TRUE)
#'   
#'   stars::st_as_stars(r)
#' }
#' 

#' 
#' #' Standard project map theme
#' #'
#' #' @return ggplot theme object.
#' theme_map <- function() {
#'   theme_void() +
#'     theme(
#'       panel.background = element_rect(fill = "#F5F5F5", color = NA),
#'       plot.margin = unit(c(0, 0, 0, 0), "pt"),
#'       legend.title = ggtext::element_markdown(lineheight = 0.95),
#'       legend.position = c(0.02, 0.02),
#'       legend.justification = c(0, 0),
#'       legend.background = element_rect(fill = "white", color = "black"),
#'       legend.margin = margin(5, 5, 5, 5),
#'       panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
#'     )
#' }
#' 
#' #' Build a formatted legend title
#' #'
#' #' @param species_name Character species common name.
#' #' @param title Main panel title.
#' #' @param subtitle Secondary title.
#' #' @param stat_label Statistic label.
#' #' @return Character string containing markdown/html for ggtext legend title.
#' legend_title_map <- function(species_name, title, subtitle, stat_label) {
#'   paste0(
#'     "<span style='font-size:20pt; font-weight:bold'>",
#'     wrap_species_label(species_name),
#'     "</span><br><br>",
#'     "<span style='font-size:14pt'>", title, "</span><br>",
#'     "<span style='font-size:7pt'>", subtitle, "</span><br>",
#'     "<span style='font-size:7pt'>", stat_label, "</span>"
#'   )
#' }
#' 
#' #' Add common overlays to a map
#' #'
#' #' @param p ggplot object.
#' #' @param study_boundary sf boundary polygon.
#' #' @param bcr_sf Optional sf BCR polygons.
#' #' @param water_sf Optional sf water polygons.
#' #' @param atlas_squares_centroids Optional sf square centroids.
#' #' @param water_fill Fill for water layer.
#' #'   "plain" for black points, or "none" for no square overlay.
#' #' @return ggplot object.
#' add_base_map_layers <- function(p,
#'                                 study_boundary,
#'                                 bcr_sf = NULL,
#'                                 water_sf = NULL,
#'                                 atlas_squares_centroids = NULL,
#'                                 water_fill = "#EDF7FB") {
#'   
#'   if (!is.null(bcr_sf)) {
#'     p <- p + geom_sf(
#'       data = bcr_sf,
#'       colour = "gray80",
#'       fill = NA,
#'       linewidth = 0.15
#'     )
#'   }
#'   
#'   
#'   if (!is.null(water_sf)) {
#'     p <- p + geom_sf(
#'       data = water_sf,
#'       fill = water_fill,
#'       col = "transparent"
#'     )
#'   }
#'   
#'   p <- p + geom_sf(
#'     data = study_boundary,
#'     colour = "gray80",
#'     fill = NA,
#'     linewidth = 0.15,
#'     show.legend = FALSE
#'   )
#'   
#'   if (!is.null(atlas_squares_centroids)) {
#'     
#'     sq <- atlas_squares_centroids %>%
#'       dplyr::filter(!is.na(n_detections))
#'     
#'     if (nrow(sq) > 0) {
#'       sq <- sq %>%
#'         dplyr::mutate(detected = n_detections > 0)
#'       
#'       # surveyed, not detected
#'       p <- p + geom_sf(
#'         data = sq, # dplyr::filter(sq, !detected),
#'         colour = "black",
#'         shape = 1,
#'         size = 0.5,
#'         alpha = 1,
#'         stroke = 0.1
#'       )
#'       
#'       # surveyed, detected
#'       p <- p + geom_sf(
#'         data = dplyr::filter(sq, detected),
#'         colour = "black",
#'         shape = 16,   #46,
#'         size = 0.5, # 0.75,
#'         stroke = 0.1,
#'         alpha = 1
#'       )
#'     }
#'   }
#'   p
#' }
#' 
#' #' Add common annotations and map extent
#' #'
#' #' @param p ggplot object.
#' #' @param study_boundary sf boundary polygon.
#' #' @param bb Optional bbox. If NULL, computed from study_boundary.
#' #' @return ggplot object.
#' add_map_annotations <- function(p, study_boundary, bb = NULL) {
#'   if (is.null(bb)) {
#'     bb <- map_bbox(study_boundary)
#'   }
#'   
#'   p +
#'     theme_map() +
#'     ggspatial::annotation_scale(location = "br", width_hint = 0.3) +
#'     ggspatial::annotation_north_arrow(
#'       location = "tr",
#'       which_north = "true",
#'       pad_x = unit(0.2, "in"),
#'       pad_y = unit(0.2, "in"),
#'       style = ggspatial::north_arrow_fancy_orienteering()
#'     ) +
#'     coord_sf(
#'       crs = sf::st_crs(study_boundary),
#'       xlim = c(bb["xmin"], bb["xmax"]),
#'       ylim = c(bb["ymin"], bb["ymax"]),
#'       expand = FALSE
#'     )
#' }
#' 
#' #' Compute a shared upper plotting limit for OBBA2 and OBBA3 q50 surfaces
#' #'
#' #' @param preds Prediction object containing OBBA2 and OBBA3 summaries.
#' #' @param prob Quantile used for upper cap.
#' #' @param fallback Fallback value if cap is invalid.
#' #' @return Positive numeric scalar.
#' compute_shared_zmax <- function(preds, prob = 0.99, fallback = 1) {
#'   zmax2 <- as.numeric(stats::quantile(preds$OBBA2$OBBA2_q50, prob, na.rm = TRUE))
#'   zmax3 <- as.numeric(stats::quantile(preds$OBBA3$OBBA3_q50, prob, na.rm = TRUE))
#'   zmax <- max(zmax2, zmax3, na.rm = TRUE)
#'   
#'   if (!is.finite(zmax) || zmax <= 0) zmax <- fallback
#'   zmax
#' }
#' 
#' #' Build Atlas 2 / Atlas 3 atlas-square centroid overlays
#' #'
#' #' @param sp_square_summary Data frame with square summaries and Atlas column.
#' #' @param atlas_sq_centroids_all sf centroid layer with square_id.
#' #' @return Named list with OBBA2 and OBBA3 sf objects (or NULL values).
#' build_atlas_square_overlays <- function(sp_square_summary, atlas_sq_centroids_all) {
#'   if (is.null(sp_square_summary)) {
#'     return(list(OBBA2 = NULL, OBBA3 = NULL))
#'   }
#'   
#'   sq2 <- dplyr::filter(sp_square_summary, Atlas == "OBBA2")
#'   sq3 <- dplyr::filter(sp_square_summary, Atlas == "OBBA3")
#'   
#'   list(
#'     OBBA2 = dplyr::left_join(atlas_sq_centroids_all, sq2, by = "square_id"),
#'     OBBA3 = dplyr::left_join(atlas_sq_centroids_all, sq3, by = "square_id")
#'   )
#' }
#' 
#' # ------------------------------------------------------------
#' # Relative abundance maps
#' # ------------------------------------------------------------
#' 
#' #' Create relative-abundance maps for one atlas period
#' #'
#' #' Returns a posterior-median map and a CV map.
#' #'
#' #' @param species_name Common name for legend title.
#' #' @param grid_sf sf prediction grid.
#' #' @param pred_summary Data frame aligned to grid_sf.
#' #' @param prefix Either "OBBA2" or "OBBA3".
#' #' @param study_boundary sf boundary polygon.
#' #' @param bcr_sf Optional sf BCR overlay.
#' #' @param water_sf Optional sf water overlay.
#' #' @param atlas_squares_centroids Optional sf square-centroid overlay.
#' #' @param title Panel title. If NULL, inferred from prefix.
#' #' @param subtitle Subtitle for legend title.
#' #' @param res Raster resolution in map units.
#' #' @param bounds Optional list with lower and upper abundance bounds.
#' #' @return List with q50_plot, cv_plot, and bounds.
#' make_relabund_maps <- function(species_name,
#'                                grid_sf,
#'                                pred_summary,
#'                                prefix = c("OBBA2", "OBBA3"),
#'                                study_boundary,
#'                                bcr_sf = NULL,
#'                                water_sf = NULL,
#'                                atlas_squares_centroids = NULL,
#'                                title = NULL,
#'                                subtitle = "Relative abundance",
#'                                res = 1000,
#'                                bounds = NULL) {
#'   prefix <- match.arg(prefix)
#'   
#'   stopifnot(inherits(grid_sf, "sf"))
#'   stopifnot(inherits(study_boundary, "sf"))
#'   stopifnot(is.data.frame(pred_summary))
#'   
#'   q50_col <- paste0(prefix, "_q50")
#'   cv_col  <- paste0(prefix, "_cv_median")
#'   sd_col  <- paste0(prefix, "_sd")
#'   
#'   if (!(q50_col %in% names(pred_summary))) {
#'     stop("pred_summary is missing required column: ", q50_col)
#'   }
#'   
#'   grid <- grid_sf
#'   grid$pred_q50 <- pred_summary[[q50_col]]
#'   
#'   if (cv_col %in% names(pred_summary)) {
#'     grid$pred_cv <- pred_summary[[cv_col]]
#'   } else if (sd_col %in% names(pred_summary)) {
#'     grid$pred_cv <- pred_summary[[sd_col]] / pmax(pred_summary[[q50_col]], 1e-12)
#'   } else {
#'     stop("pred_summary must contain either ", cv_col, " or ", sd_col)
#'   }
#'   
#'   if (is.null(bounds)) {
#'     upper <- as.numeric(stats::quantile(grid$pred_q50, 0.99, na.rm = TRUE))
#'     if (!is.finite(upper) || upper <= 0) upper <- 1
#'     bounds <- list(lower = 0, upper = upper)
#'   } else {
#'     if (is.null(bounds$lower)) bounds$lower <- 0
#'     if (is.null(bounds$upper)) {
#'       bounds$upper <- as.numeric(stats::quantile(grid$pred_q50, 0.99, na.rm = TRUE))
#'       if (!is.finite(bounds$upper) || bounds$upper <= 0) bounds$upper <- 1
#'     }
#'   }
#'   
#'   grid$pred_capped_q50 <- pmax(pmin(grid$pred_q50, bounds$upper), bounds$lower)
#'   grid$pred_capped_cv  <- pmax(pmin(grid$pred_cv, 1), 0)
#'   
#'   q50_stars <- rasterize_sf_to_stars(grid, "pred_capped_q50", res = res)
#'   cv_stars  <- rasterize_sf_to_stars(grid, "pred_capped_cv",  res = res)
#'   
#'   if (is.null(title)) {
#'     title <- ifelse(prefix == "OBBA2", "Atlas 2", "Atlas 3")
#'   }
#'   
#'   bb <- map_bbox(study_boundary)
#'   
#'   q50_cols <- colorRampPalette(c(
#'     "#FEFEFE", "#FBF7E2", "#FCF8D0", "#EEF7C2", "#CEF2B0",
#'     "#94E5A0", "#51C987", "#18A065", "#008C59", "#007F53", "#006344"
#'   ))(10)
#'   
#'   cv_cols <- colorRampPalette(c(
#'     "#FEFEFE", "#FFF4B3", "#F5D271", "#F2B647", "#EC8E00", "#CA302A"
#'   ))(11)
#'   
#'   q50_plot <- ggplot() +
#'     stars::geom_stars(data = q50_stars)
#'   
#'   q50_plot <- add_base_map_layers(
#'     p = q50_plot,
#'     study_boundary = study_boundary,
#'     bcr_sf = bcr_sf,
#'     water_sf = water_sf,
#'     atlas_squares_centroids = atlas_squares_centroids,
#'     water_fill = "#EDF7FB"
#'   )
#'   
#'   q50_plot <- q50_plot +
#'     scale_fill_gradientn(
#'       name = legend_title_map(species_name, title, subtitle, "Posterior median"),
#'       colors = q50_cols,
#'       na.value = "transparent",
#'       limits = c(bounds$lower, bounds$upper)
#'     )
#'   
#'   q50_plot <- add_map_annotations(q50_plot, study_boundary, bb)
#'   
#'   cv_plot <- ggplot() +
#'     stars::geom_stars(data = cv_stars)
#'   
#'   cv_plot <- add_base_map_layers(
#'     p = cv_plot,
#'     study_boundary = study_boundary,
#'     bcr_sf = bcr_sf,
#'     water_sf = water_sf,
#'     atlas_squares_centroids = atlas_squares_centroids,
#'     water_fill = "#EDF7FB"
#'   )
#'   
#'   cv_plot <- cv_plot +
#'     scale_fill_gradientn(
#'       name = legend_title_map(species_name, title, subtitle, "Prediction CV"),
#'       colors = cv_cols,
#'       na.value = "transparent",
#'       breaks = seq(0, 1, length.out = 5),
#'       labels = c("0", "0.25", "0.5", "0.75", ">1"),
#'       limits = c(0, 1)
#'     )
#'   
#'   cv_plot <- add_map_annotations(cv_plot, study_boundary, bb)
#'   
#'   list(
#'     q50_plot = q50_plot,
#'     cv_plot = cv_plot,
#'     bounds = bounds
#'   )
#' }
#' 
#' # ------------------------------------------------------------
#' # Absolute change maps
#' # ------------------------------------------------------------
#' 
#' #' Create absolute-change maps between atlas periods
#' #'
#' #' Returns a posterior-median absolute-change map and a CI-width map.
#' #'
#' #' @param species_name Common name for legend title.
#' #' @param grid_sf sf prediction grid (typically the OBBA3 grid geometry).
#' #' @param abs_change_summary Data frame aligned to grid_sf.
#' #' @param study_boundary sf boundary polygon.
#' #' @param bcr_sf Optional sf BCR overlay.
#' #' @param water_sf Optional sf water overlay.
#' #' @param signif_poly Optional sf polygons classed as Increase/Decrease.
#' #' @param res Raster resolution in map units.
#' #' @param max_abs Optional symmetric upper plotting bound for change.
#' #' @param title Panel title.
#' #' @param subtitle Subtitle for legend title.
#' #' @return List with chg_plot, ciw_plot, and max_abs.
#' make_abs_change_maps <- function(species_name,
#'                                  grid_sf,
#'                                  abs_change_summary,
#'                                  study_boundary,
#'                                  bcr_sf = NULL,
#'                                  water_sf = NULL,
#'                                  signif_poly = NULL,
#'                                  res = 1000,
#'                                  max_abs = NULL,
#'                                  title = "Absolute change",
#'                                  subtitle = "OBBA2 to OBBA3") {
#'   stopifnot(inherits(grid_sf, "sf"))
#'   stopifnot(inherits(study_boundary, "sf"))
#'   stopifnot(is.data.frame(abs_change_summary))
#'   
#'   required_cols <- c("abs_change_q50", "abs_change_lower", "abs_change_upper")
#'   if (!all(required_cols %in% names(abs_change_summary))) {
#'     stop(
#'       "abs_change_summary must contain: ",
#'       paste(required_cols, collapse = ", ")
#'     )
#'   }
#'   
#'   grid <- grid_sf
#'   grid$chg_q50 <- abs_change_summary$abs_change_q50
#'   grid$chg_ciw <- abs_change_summary$abs_change_upper - abs_change_summary$abs_change_lower
#'   
#'   if (is.null(max_abs)) {
#'     max_abs <- as.numeric(stats::quantile(abs(grid$chg_q50), 0.99, na.rm = TRUE))
#'     
#'     if (!is.finite(max_abs) || max_abs <= 0) {
#'       max_abs <- max(abs(grid$chg_q50), na.rm = TRUE)
#'     }
#'     
#'     if (!is.finite(max_abs) || max_abs <= 0) {
#'       max_abs <- 1
#'     }
#'   } else {
#'     max_abs <- as.numeric(max_abs)
#'     if (!is.finite(max_abs) || max_abs <= 0) max_abs <- 1
#'   }
#'   
#'   grid$chg_capped <- pmax(pmin(grid$chg_q50, max_abs), -max_abs)
#'   grid$ciw_capped <- pmax(pmin(grid$chg_ciw, max_abs), 0)
#'   
#'   chg_stars <- rasterize_sf_to_stars(grid, "chg_capped", res = res)
#'   ciw_stars <- rasterize_sf_to_stars(grid, "ciw_capped", res = res)
#'   
#'   bb <- map_bbox(study_boundary)
#'   
#'   chg_cols <- colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(11)
#'   unc_cols <- colorRampPalette(c(
#'     "#FEFEFE", "#FFF4B3", "#F5D271", "#F2B647", "#EC8E00", "#CA302A"
#'   ))(11)
#'   
#'   chg_breaks <- seq(-max_abs, max_abs, length.out = 7)
#'   chg_breaks[4] <- 0
#'   ciw_breaks <- seq(0, max_abs, length.out = 7)
#'   
#'   chg_plot <- ggplot() +
#'     stars::geom_stars(data = chg_stars)
#'   
#'   chg_plot <- add_base_map_layers(
#'     p = chg_plot,
#'     study_boundary = study_boundary,
#'     bcr_sf = bcr_sf,
#'     water_sf = water_sf,
#'     atlas_squares_centroids = NULL,
#'     water_fill = "#F5F5F5"
#'   )
#'   
#'   chg_plot <- chg_plot +
#'     scale_fill_gradientn(
#'       name = legend_title_map(species_name, title, subtitle, "Posterior median"),
#'       colors = chg_cols,
#'       na.value = "transparent",
#'       breaks = chg_breaks,
#'       labels = signif(chg_breaks, 2),
#'       limits = c(-max_abs, max_abs)
#'     )
#'   
#'   if (!is.null(signif_poly) && nrow(signif_poly) > 0) {
#'     if ("classification" %in% names(signif_poly)) {
#'       if (any(signif_poly$classification == "Increase")) {
#'         chg_plot <- chg_plot +
#'           geom_sf(
#'             data = subset(signif_poly, classification == "Increase"),
#'             fill = "transparent",
#'             col = "dodgerblue",
#'             linewidth = 1
#'           )
#'       }
#'       
#'       if (any(signif_poly$classification == "Decrease")) {
#'         chg_plot <- chg_plot +
#'           geom_sf(
#'             data = subset(signif_poly, classification == "Decrease"),
#'             fill = "transparent",
#'             col = "orangered",
#'             linewidth = 1
#'           )
#'       }
#'     }
#'   }
#'   
#'   chg_plot <- add_map_annotations(chg_plot, study_boundary, bb)
#'   
#'   ciw_plot <- ggplot() +
#'     stars::geom_stars(data = ciw_stars)
#'   
#'   ciw_plot <- add_base_map_layers(
#'     p = ciw_plot,
#'     study_boundary = study_boundary,
#'     bcr_sf = bcr_sf,
#'     water_sf = water_sf,
#'     atlas_squares_centroids = NULL,
#'     water_fill = "#F5F5F5"
#'   )
#'   
#'   ciw_plot <- ciw_plot +
#'     scale_fill_gradientn(
#'       name = legend_title_map(species_name, title, subtitle, "90% CI width"),
#'       colors = unc_cols,
#'       na.value = "transparent",
#'       breaks = ciw_breaks,
#'       labels = signif(ciw_breaks, 2),
#'       limits = c(0, max_abs)
#'     )
#'   
#'   ciw_plot <- add_map_annotations(ciw_plot, study_boundary, bb)
#'   
#'   list(
#'     chg_plot = chg_plot,
#'     ciw_plot = ciw_plot,
#'     max_abs = max_abs
#'   )
#' }
#' 
#' # ------------------------------------------------------------
#' # Meaningful change polygons
#' # ------------------------------------------------------------
#' 
#' #' Flag hexagons with strong posterior support for change
#' #'
#' #' @param eta_draws_per_hex Posterior draws per polygon (hex), in the format
#' #'   expected by `summarize_polygon_hypothesis()`.
#' #' @param hexagon_sf `sf` object of hex polygons. Must include `hex_id`.
#' #' @param param Parameter name to test inside `summarize_polygon_hypothesis()`
#' #'   (e.g., "abs_change").
#' #' @param threshold Numeric threshold for change hypothesis.
#' #' @param prob_level Posterior probability threshold (e.g., 0.975).
#' #' @param direction Directions to compute in `summarize_polygon_hypothesis()`.
#' #' @param ci_probs Credible interval probs passed through.
#' #' @param include_summary Passed through to `summarize_polygon_hypothesis()`.
#' #'
#' #' @return An `sf` object of flagged hexagons joined to summaries, filtered to
#' #'   `classification != "None"`. If none flagged, returns 0-row `sf`.
#' #' @export
#' flag_hexagons_for_change <- function(eta_draws_per_hex,
#'                                      hexagon_sf,
#'                                      param = "abs_change",
#'                                      threshold = 0,
#'                                      prob_level = 0.975,
#'                                      direction = c("two_sided", "increase", "decrease"),
#'                                      ci_probs = c(0.05, 0.95),
#'                                      include_summary = TRUE) {
#'   hexagon_sf <- hexagon_sf %>% dplyr::mutate(hex_id = as.character(hex_id))
#'   
#'   flagged <- summarize_polygon_hypothesis(
#'     eta_draws_per_hex,
#'     param = param,
#'     threshold = threshold,
#'     prob_level = prob_level,
#'     direction = direction,
#'     ci_probs = ci_probs,
#'     include_summary = include_summary
#'   ) %>%
#'     dplyr::left_join(hexagon_sf, ., by = c("hex_id" = "poly_id")) %>%
#'     dplyr::filter(classification != "None") %>%
#'     sf::st_as_sf()
#'   
#'   flagged
#' }
#' 
#' 
#' #' Keep only spatially coherent patches above a minimum area threshold
#' #'
#' #' For each change class (e.g., Increase/Decrease), this:
#' #'   - builds a touches-based adjacency graph of hexagons,
#' #'   - finds connected components (patches),
#' #'   - computes patch areas after dissolving,
#' #'   - retains hexagons belonging to patches with area >= min_area_km2.
#' #'
#' #' @param flagged_hex_sf `sf` of flagged hexagons. Must include `classification`
#' #'   and geometry.
#' #' @param min_area_km2 Minimum patch area (km^2) to retain.
#' #'
#' #' @return `sf` of hexagons with an added `patch_id`, filtered to retained patches.
#' #'   If input has 0 rows, returns input unchanged.
#' #' @export
#' filter_flagged_to_large_patches <- function(flagged_hex_sf, min_area_km2) {
#'   if (nrow(flagged_hex_sf) == 0) return(flagged_hex_sf)
#'   
#'   flagged_patched <- flagged_hex_sf %>%
#'     dplyr::group_by(classification) %>%
#'     dplyr::group_modify(~{
#'       x <- .x
#'       
#'       # Neighbor list of touching polygons
#'       nb <- sf::st_touches(x)
#'       
#'       # Graph + connected components
#'       g <- igraph::graph_from_adj_list(nb, mode = "all")
#'       g <- igraph::as.undirected(g, mode = "collapse")
#'       x$patch_id <- igraph::components(g)$membership
#'       
#'       # Dissolve by patch and compute patch areas (km^2)
#'       patch_polys <- x %>%
#'         dplyr::group_by(patch_id) %>%
#'         # small "buffer in/out" trick to help remove tiny holes / slivers
#'         sf::st_buffer(sqrt(min_area_km2)) %>%
#'         dplyr::summarise(geometry = sf::st_union(geometry), .groups = "drop") %>%
#'         sf::st_buffer(-sqrt(min_area_km2)) %>%
#'         dplyr::mutate(
#'           patch_area_km2 =
#'             units::set_units(sf::st_area(geometry), "km^2") %>%
#'             units::drop_units()
#'         )
#'       
#'       keep_ids <- patch_polys$patch_id[patch_polys$patch_area_km2 >= min_area_km2]
#'       x %>% dplyr::filter(patch_id %in% keep_ids)
#'     }) %>%
#'     dplyr::ungroup() %>%
#'     sf::st_as_sf()
#'   
#'   flagged_patched
#' }
#' 
#' 
#' #' Dissolve hexagons into polygons by change classification
#' #'
#' #' @param flagged_hex_sf `sf` of flagged hexagons with a `classification` column.
#' #'
#' #' @return `sf` with one (multi)polygon per classification.
#' #' @export
#' dissolve_change_polygons <- function(flagged_hex_sf) {
#'   flagged_hex_sf %>%
#'     dplyr::group_by(classification) %>%
#'     dplyr::summarise(geometry = sf::st_union(geometry), .groups = "drop") %>%
#'     sf::st_as_sf()
#' }
#' 
#' 
#' #' Smooth change polygon boundaries and clip to the study boundary
#' #'
#' #' @param change_polys_sf `sf` polygons (usually output of `dissolve_change_polygons()`).
#' #' @param study_boundary `sf` polygon used for CRS target and clipping.
#' #' @param smoothing_bandwidth_m Numeric bandwidth (meters) for ksmooth.
#' #'
#' #' @return Smoothed, valid, boundary-clipped `sf` polygons.
#' #' @export
#' smooth_and_clip_change_polys <- function(change_polys_sf,
#'                                          study_boundary,
#'                                          smoothing_bandwidth_m) {
#'   change_polys_sf %>%
#'     sf::st_transform(sf::st_crs(study_boundary)) %>%
#'     smoothr::smooth(method = "ksmooth", bandwidth = smoothing_bandwidth_m) %>%
#'     sf::st_make_valid() %>%
#'     sf::st_intersection(study_boundary) %>%
#'     sf::st_as_sf()
#' }
#' 
#' 
#' #' Drop holes smaller than a threshold area (km^2) from polygons/multipolygons
#' #'
#' #' @param x An `sf` object or `sfc` geometry with POLYGON/MULTIPOLYGON features.
#' #' @param min_hole_area_km2 Minimum hole area to keep (km^2). Holes smaller than
#' #'   this are removed.
#' #'
#' #' @return Object of same class as input with holes removed where applicable.
#' #' @export
#' drop_small_holes <- function(x, min_hole_area_km2) {
#'   
#'   stopifnot(inherits(x, c("sf", "sfc")))
#'   
#'   geom <- sf::st_geometry(x)
#'   crs  <- sf::st_crs(geom)
#'   
#'   if (sf::st_is_longlat(geom)) {
#'     warning("Geometry is lon/lat. Areas will be computed using s2 (geodesic).")
#'   }
#'   
#'   drop_in_poly <- function(poly) {
#'     
#'     rings <- unclass(poly)
#'     outer <- rings[[1]]
#'     holes <- rings[-1]
#'     
#'     if (length(holes) == 0) return(poly)
#'     
#'     hole_areas_km2 <- vapply(
#'       holes,
#'       function(h) {
#'         a <- sf::st_area(sf::st_sfc(sf::st_polygon(list(h)), crs = crs))
#'         as.numeric(units::set_units(a, km^2))
#'       },
#'       numeric(1)
#'     )
#'     
#'     keep <- hole_areas_km2 >= min_hole_area_km2
#'     sf::st_polygon(c(list(outer), holes[keep]))
#'   }
#'   
#'   drop_in_mpoly <- function(mpoly) {
#'     
#'     polys <- unclass(mpoly)
#'     
#'     new_polys <- lapply(polys, function(rings) {
#'       drop_in_poly(sf::st_polygon(rings))
#'     })
#'     
#'     sf::st_multipolygon(lapply(new_polys, unclass))
#'   }
#'   
#'   new_geom <- lapply(geom, function(g) {
#'     if (inherits(g, "POLYGON")) {
#'       drop_in_poly(g)
#'     } else if (inherits(g, "MULTIPOLYGON")) {
#'       drop_in_mpoly(g)
#'     } else {
#'       g
#'     }
#'   })
#'   
#'   if (inherits(x, "sf")) {
#'     sf::st_geometry(x) <- sf::st_sfc(new_geom, crs = crs)
#'     return(x)
#'   } else {
#'     return(sf::st_sfc(new_geom, crs = crs))
#'   }
#' }
#' 
#' 
#' #' End-to-end builder for "meaningful change" polygons (flag → patch → dissolve → smooth)
#' #'
#' #' This is a convenience wrapper that matches your current inline workflow:
#' #' - flag hexagons based on posterior support,
#' #' - retain only patches exceeding `min_area_km2`,
#' #' - dissolve by Increase/Decrease,
#' #' - smooth and clip to study boundary,
#' #' - optionally remove small holes.
#' #'
#' #' @param eta_draws_per_hex Posterior draws per hex (see `flag_hexagons_for_change()`).
#' #' @param hexagon_sf `sf` hex grid with `hex_id`.
#' #' @param study_boundary `sf` boundary for CRS + clipping.
#' #' @param min_area_km2 Minimum patch size (km^2).
#' #' @param smoothing_bandwidth_m Smoothing bandwidth (meters).
#' #' @param ... Additional args forwarded to `flag_hexagons_for_change()`
#' #'   (e.g., `param`, `threshold`, `prob_level`, etc.).
#' #' @param drop_holes Logical; if TRUE, remove holes smaller than `min_area_km2`.
#' #'
#' #' @return Smoothed `sf` polygons by classification, or `NULL` if none flagged.
#' #' @export
#' build_meaningful_change_polys <- function(eta_draws_per_hex,
#'                                           hexagon_sf,
#'                                           study_boundary,
#'                                           min_area_km2 = 5000,
#'                                           smoothing_bandwidth_m = 75000,
#'                                           ...,
#'                                           drop_holes = TRUE) {
#'   
#'   flagged <- flag_hexagons_for_change(
#'     eta_draws_per_hex = eta_draws_per_hex,
#'     hexagon_sf = hexagon_sf,
#'     ...
#'   )
#'   
#'   if (nrow(flagged) == 0) return(NULL)
#'   
#'   flagged_patched <- filter_flagged_to_large_patches(
#'     flagged_hex_sf = flagged,
#'     min_area_km2 = min_area_km2
#'   )
#'   
#'   if (nrow(flagged_patched) == 0) return(NULL)
#'   
#'   change_polys <- dissolve_change_polygons(flagged_patched)
#'   
#'   change_polys_smoothed <- smooth_and_clip_change_polys(
#'     change_polys_sf = change_polys,
#'     study_boundary = study_boundary,
#'     smoothing_bandwidth_m = smoothing_bandwidth_m
#'   )
#'   
#'   if (drop_holes) {
#'     change_polys_smoothed <- drop_small_holes(change_polys_smoothed, min_area_km2)
#'   }
#'   
#'   change_polys_smoothed
#' }
#' 
#' 
#' 
#' 
#' #' Summarize posterior change within each hexagon and join to hex sf
#' #'
#' #' For each hexagon, compute posterior summaries for a parameter such as
#' #' absolute change, then join those summaries onto the hex geometry.
#' #'
#' #' @param eta_draws_per_hex Named list of posterior draws by hexagon.
#' #''   Expected structure:
#' #'   eta_draws_per_hex[[hex_id]][[param]] = numeric vector of draws
#' #' @param hexagon_sf sf object of hexagons. Must include `hex_id`.
#' #' @param param Name of parameter to summarize, default "abs_change".
#' #' @param threshold Threshold for evaluating probability of meaningful change.
#' #' @param ci_probs Length-2 numeric vector giving lower/upper credible interval
#' #'   probabilities, e.g. c(0.05, 0.95).
#' #' @param include_median Logical; if TRUE, also compute posterior median.
#' #'
#' #' @return sf object with one row per hexagon and posterior summary columns joined.
#' #' @export
#' summarize_change_by_hex <- function(eta_draws_per_hex,
#'                                     hexagon_sf,
#'                                     param = "abs_change",
#'                                     threshold = 0,
#'                                     ci_probs = c(0.05, 0.95),
#'                                     include_median = TRUE) {
#'   
#'   stopifnot(inherits(hexagon_sf, "sf"))
#'   stopifnot("hex_id" %in% names(hexagon_sf))
#'   stopifnot(length(ci_probs) == 2)
#'   
#'   hex_sf <- hexagon_sf %>%
#'     dplyr::mutate(hex_id = as.character(hex_id))
#'   
#'   # helper to extract a numeric draw vector for one hex
#'   extract_draws <- function(x, param) {
#'     if (is.null(x)) return(NULL)
#'     
#'     # most likely case: list element contains named components
#'     if (is.list(x) && !is.null(x[[param]])) {
#'       return(as.numeric(x[[param]]))
#'     }
#'     
#'     # alternative: data.frame/tibble with a column named param
#'     if (is.data.frame(x) && param %in% names(x)) {
#'       return(as.numeric(x[[param]]))
#'     }
#'     
#'     # fallback: if x itself is the numeric draw vector and param is implicit
#'     if (is.numeric(x)) {
#'       return(as.numeric(x))
#'     }
#'     
#'     NULL
#'   }
#'   
#'   # if list is unnamed, try to align by order
#'   if (is.null(names(eta_draws_per_hex))) {
#'     if (length(eta_draws_per_hex) != nrow(hex_sf)) {
#'       stop("eta_draws_per_hex is unnamed and does not match nrow(hexagon_sf).")
#'     }
#'     names(eta_draws_per_hex) <- hex_sf$hex_id
#'   }
#'   
#'   hex_summaries <- lapply(names(eta_draws_per_hex), function(id) {
#'     draws <- extract_draws(eta_draws_per_hex[[id]], param = param)
#'     
#'     if (is.null(draws) || length(draws) == 0 || all(is.na(draws))) {
#'       return(data.frame(
#'         hex_id = as.character(id),
#'         n_draws = 0,
#'         mean_change = NA_real_,
#'         median_change = NA_real_,
#'         lower = NA_real_,
#'         upper = NA_real_,
#'         prob_gt_threshold = NA_real_,
#'         prob_lt_neg_threshold = NA_real_,
#'         prob_abs_gt_threshold = NA_real_
#'       ))
#'     }
#'     
#'     draws <- draws[is.finite(draws)]
#'     
#'     if (length(draws) == 0) {
#'       return(data.frame(
#'         hex_id = as.character(id),
#'         n_draws = 0,
#'         mean_change = NA_real_,
#'         median_change = NA_real_,
#'         lower = NA_real_,
#'         upper = NA_real_,
#'         prob_gt_threshold = NA_real_,
#'         prob_lt_neg_threshold = NA_real_,
#'         prob_abs_gt_threshold = NA_real_
#'       ))
#'     }
#'     
#'     out <- data.frame(
#'       hex_id = as.character(id),
#'       n_draws = length(draws),
#'       mean_change = mean(draws),
#'       lower = as.numeric(stats::quantile(draws, ci_probs[1], na.rm = TRUE)),
#'       upper = as.numeric(stats::quantile(draws, ci_probs[2], na.rm = TRUE)),
#'       prob_gt_threshold = mean(draws > threshold),
#'       prob_lt_neg_threshold = mean(draws < -threshold),
#'       prob_abs_gt_threshold = mean(abs(draws) > threshold)
#'     )
#'     
#'     if (include_median) {
#'       out$median_change <- stats::median(draws, na.rm = TRUE)
#'       out <- out[, c(
#'         "hex_id", "n_draws", "mean_change", "median_change", "lower", "upper",
#'         "prob_gt_threshold", "prob_lt_neg_threshold", "prob_abs_gt_threshold"
#'       )]
#'     } else {
#'       out$median_change <- NULL
#'     }
#'     
#'     out
#'   })
#'   
#'   hex_summaries <- dplyr::bind_rows(hex_summaries)
#'   
#'   hex_sf %>%
#'     dplyr::left_join(hex_summaries, by = "hex_id") %>%
#'     sf::st_as_sf()
#' }
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' make_relabund_map <- function(species_name,
#'                               relabund_rast = NULL,
#'                               zlim = NULL,
#'                               layer = "mu_q50",
#'                               prefix = c("OBBA2", "OBBA3"),
#'                               study_boundary,
#'                               bcr_sf = NULL,
#'                               water_sf = NULL,
#'                               water_fill = "#EDF7FB",
#'                               sp_hex_grid = NULL,
#'                               survey_colour = "#bdbdbd",
#'                               detection_colour = "black",
#'                               max_circle_diameter_km = 10,
#'                               n_surveys_size_quantile = 0.90,
#'                               title = NULL,
#'                               subtitle = "Relative abundance") {
#'   
#'   prefix <- match.arg(prefix)
#'   
#'   stopifnot(inherits(study_boundary, "sf"))
#'   
#'   if (is.null(title)) {
#'     title <- ifelse(prefix == "OBBA2", "Atlas 2", "Atlas 3")
#'   }
#'   
#'   # ------------------------------------------------------------
#'   # Determine CRS
#'   # ------------------------------------------------------------
#'   if (!is.null(relabund_rast)) {
#'     plot_crs <- sf::st_crs(terra::crs(relabund_rast))
#'   } else {
#'     plot_crs <- sf::st_crs(study_boundary)
#'   }
#'   
#'   if (is.na(plot_crs)) {
#'     stop("study_boundary must have valid CRS if raster is NULL")
#'   }
#'   
#'   study_boundary <- sf::st_transform(study_boundary, plot_crs)
#'   if (!is.null(bcr_sf))   bcr_sf   <- sf::st_transform(bcr_sf, plot_crs)
#'   if (!is.null(water_sf)) water_sf <- sf::st_transform(water_sf, plot_crs)
#'   
#'   # ------------------------------------------------------------
#'   # Raster (optional)
#'   # ------------------------------------------------------------
#'   r_stars <- NULL
#'   
#'   if (!is.null(relabund_rast)) {
#'     if (is.null(zlim)) stop("zlim required when raster present")
#'     
#'     if (length(zlim) == 1) zlim <- c(0, zlim)
#'     
#'     r_plot <- relabund_rast[[layer]]
#'     names(r_plot) <- "relabund"
#'     r_stars <- stars::st_as_stars(r_plot)
#'   }
#'   
#'   q50_cols <- grDevices::colorRampPalette(c(
#'     "#FEFEFE", "#FBF7E2", "#FCF8D0", "#EEF7C2", "#CEF2B0",
#'     "#94E5A0", "#51C987", "#18A065", "#008C59", "#007F53", "#006344"
#'   ))(10)
#'   
#'   # ------------------------------------------------------------
#'   # Hex data → circle radii
#'   # ------------------------------------------------------------
#'   circle_df <- NULL
#'   
#'   if (!is.null(sp_hex_grid)) {
#'     
#'     required_cols <- c("Atlas", "n_surveys", "n_surveys_det")
#'     missing_cols <- setdiff(required_cols, names(sp_hex_grid))
#'     
#'     if (length(missing_cols) > 0) {
#'       stop("Missing columns: ", paste(missing_cols, collapse = ", "))
#'     }
#'     
#'     sp_hex_grid <- sp_hex_grid %>%
#'       sf::st_transform(plot_crs) %>%
#'       dplyr::filter(
#'         Atlas == prefix,
#'         !is.na(n_surveys),
#'         n_surveys > 0
#'       ) %>%
#'       dplyr::mutate(
#'         n_surveys_det = dplyr::coalesce(n_surveys_det, 0),
#'         n_surveys_det = pmin(n_surveys_det, n_surveys)
#'       )
#'     
#'     if (nrow(sp_hex_grid) > 0) {
#'       
#'       n_max <- stats::quantile(
#'         sp_hex_grid$n_surveys,
#'         probs = n_surveys_size_quantile,
#'         na.rm = TRUE
#'       )
#'       
#'       if (!is.finite(n_max) || n_max <= 0) {
#'         n_max <- max(sp_hex_grid$n_surveys, na.rm = TRUE)
#'       }
#'       
#'       # convert km → CRS units
#'       crs_units <- sf::st_crs(plot_crs)$units_gdal
#'       km_to_units <- ifelse(grepl("metre", crs_units), 1000, 1)
#'       
#'       max_radius <- (max_circle_diameter_km / 2) * km_to_units
#'       
#'       circle_df <- sp_hex_grid %>%
#'         dplyr::mutate(
#'           n_surveys_plot = pmin(n_surveys, n_max),
#'           n_det_plot     = pmin(n_surveys_det, n_max)
#'         ) %>%
#'         sf::st_centroid() %>%
#'         dplyr::mutate(
#'           x = sf::st_coordinates(.)[,1],
#'           y = sf::st_coordinates(.)[,2],
#'           
#'           # area ∝ n  → radius ∝ sqrt(n)
#'           r_gray  = sqrt(n_surveys_plot / n_max) * max_radius,
#'           r_black = sqrt(n_det_plot     / n_max) * max_radius
#'         ) %>%
#'         sf::st_drop_geometry()
#'     }
#'   }
#'   
#'   bb <- map_bbox(study_boundary)
#'   
#'   # ------------------------------------------------------------
#'   # Build plot
#'   # ------------------------------------------------------------
#'   p <- ggplot2::ggplot()
#'   
#'   if (!is.null(r_stars)) {
#'     p <- p +
#'       stars::geom_stars(data = r_stars) +
#'       ggplot2::scale_fill_gradientn(
#'         name = legend_title_map(species_name, title, subtitle, "Posterior median"),
#'         colors = q50_cols,
#'         na.value = "transparent",
#'         limits = zlim
#'       )
#'   }
#'   
#'   p <- add_base_map_layers(
#'     p = p,
#'     study_boundary = study_boundary,
#'     bcr_sf = bcr_sf,
#'     water_sf = water_sf,
#'     atlas_squares_centroids = NULL,
#'     water_fill = water_fill
#'   )
#'   
#'   # ------------------------------------------------------------
#'   # Draw circles (key change)
#'   # ------------------------------------------------------------
#'   if (!is.null(circle_df) && nrow(circle_df) > 0) {
#'     
#'     p <- p +
#'       ggforce::geom_circle(
#'         data = circle_df,
#'         ggplot2::aes(x0 = x, y0 = y, r = r_gray),
#'         fill = survey_colour,
#'         color = NA
#'       ) +
#'       ggforce::geom_circle(
#'         data = dplyr::filter(circle_df, r_black > 0),
#'         ggplot2::aes(x0 = x, y0 = y, r = r_black),
#'         fill = detection_colour,
#'         color = NA
#'       )
#'   }
#'   
#'   p <- add_map_annotations(p, study_boundary, bb)
#'   
#'   list(
#'     sp_plot = p,
#'     zlim = zlim
#'   )
#' }
#' 
#' 
#' 
#' 
#' make_relabund_map <- function(species_name,
#'                               relabund_rast = NULL,
#'                               zlim = NULL,
#'                               layer = "mu_q50",
#'                               prefix = c("OBBA2", "OBBA3"),
#'                               study_boundary,
#'                               bcr_sf = NULL,
#'                               water_sf = NULL,
#'                               water_fill = "#EDF7FB",
#'                               sp_hex_grid = NULL,
#'                               survey_colour = "#bdbdbd",
#'                               detection_colour = "black",
#'                               survey_circle_diameter_km = 20,
#'                               min_detection_diameter_km = 1,
#'                               max_detection_diameter_km = 15,
#'                               mean_count_scaling_max = NULL,
#'                               survey_scaling_max = 20,
#'                               alpha_min = 0.2,
#'                               alpha_max = 1,
#'                               title = NULL,
#'                               subtitle = "Relative abundance") {
#'   
#'   prefix <- match.arg(prefix)
#'   stopifnot(inherits(study_boundary, "sf"))
#'   
#'   if (survey_scaling_max <= 1) {
#'     stop("survey_scaling_max must be greater than 1.")
#'   }
#'   
#'   if (min_detection_diameter_km < 0) {
#'     stop("min_detection_diameter_km must be >= 0.")
#'   }
#'   
#'   if (max_detection_diameter_km <= 0) {
#'     stop("max_detection_diameter_km must be > 0.")
#'   }
#'   
#'   if (min_detection_diameter_km > max_detection_diameter_km) {
#'     stop("min_detection_diameter_km must be <= max_detection_diameter_km.")
#'   }
#'   
#'   if (alpha_min < 0 || alpha_min > 1 || alpha_max < 0 || alpha_max > 1) {
#'     stop("alpha_min and alpha_max must both be between 0 and 1.")
#'   }
#'   
#'   if (alpha_min > alpha_max) {
#'     stop("alpha_min must be <= alpha_max.")
#'   }
#'   
#'   if (is.null(title)) {
#'     title <- ifelse(prefix == "OBBA2", "Atlas 2", "Atlas 3")
#'   }
#'   
#'   # ------------------------------------------------------------
#'   # Determine plotting CRS
#'   # ------------------------------------------------------------
#'   if (!is.null(relabund_rast)) {
#'     stopifnot(inherits(relabund_rast, "SpatRaster"))
#'     plot_crs <- sf::st_crs(terra::crs(relabund_rast))
#'   } else {
#'     plot_crs <- sf::st_crs(study_boundary)
#'   }
#'   
#'   if (is.na(plot_crs)) {
#'     stop("study_boundary must have a valid CRS if relabund_rast is NULL.")
#'   }
#'   
#'   study_boundary <- sf::st_transform(study_boundary, plot_crs)
#'   if (!is.null(bcr_sf))   bcr_sf   <- sf::st_transform(bcr_sf, plot_crs)
#'   if (!is.null(water_sf)) water_sf <- sf::st_transform(water_sf, plot_crs)
#'   
#'   # ------------------------------------------------------------
#'   # Optional raster
#'   # ------------------------------------------------------------
#'   r_stars <- NULL
#'   
#'   if (!is.null(relabund_rast)) {
#'     if (is.null(zlim)) {
#'       stop("zlim is required when relabund_rast is not NULL.")
#'     }
#'     
#'     if (length(zlim) == 1) zlim <- c(0, zlim)
#'     
#'     if (!layer %in% names(relabund_rast)) {
#'       stop("relabund_rast is missing layer: ", layer)
#'     }
#'     
#'     r_plot <- relabund_rast[[layer]]
#'     names(r_plot) <- "relabund"
#'     r_stars <- stars::st_as_stars(r_plot)
#'   }
#'   
#'   q50_cols <- grDevices::colorRampPalette(c(
#'     "#FEFEFE", "#FBF7E2", "#FCF8D0", "#EEF7C2", "#CEF2B0",
#'     "#94E5A0", "#51C987", "#18A065", "#008C59", "#007F53", "#006344"
#'   ))(10)
#'   
#'   # ------------------------------------------------------------
#'   # Optional raw-data circles
#'   # ------------------------------------------------------------
#'   circle_df <- NULL
#'   
#'   if (!is.null(sp_hex_grid)) {
#'     stopifnot(inherits(sp_hex_grid, "sf"))
#'     
#'     required_cols <- c("Atlas", "n_surveys", "mean_count")
#'     missing_cols <- setdiff(required_cols, names(sp_hex_grid))
#'     
#'     if (length(missing_cols) > 0) {
#'       stop("sp_hex_grid is missing required column(s): ",
#'            paste(missing_cols, collapse = ", "))
#'     }
#'     
#'     sp_hex_grid <- sp_hex_grid %>%
#'       sf::st_transform(plot_crs) %>%
#'       dplyr::filter(
#'         Atlas == prefix,
#'         !is.na(n_surveys),
#'         n_surveys > 0
#'       ) %>%
#'       dplyr::mutate(
#'         mean_count = dplyr::coalesce(mean_count, 0)
#'       )
#'     
#'     if (nrow(sp_hex_grid) > 0) {
#'       crs_units <- sf::st_crs(plot_crs)$units_gdal
#'       
#'       km_to_units <- dplyr::case_when(
#'         grepl("metre|meter", crs_units, ignore.case = TRUE) ~ 1000,
#'         grepl("kilometre|kilometer", crs_units, ignore.case = TRUE) ~ 1,
#'         TRUE ~ 1000
#'       )
#'       
#'       survey_radius <- (survey_circle_diameter_km / 2) * km_to_units
#'       min_detection_radius <- (min_detection_diameter_km / 2) * km_to_units
#'       max_detection_radius <- (max_detection_diameter_km / 2) * km_to_units
#'       
#'       if (is.null(mean_count_scaling_max)) {
#'         mean_count_scaling_max <- max(sp_hex_grid$mean_count, na.rm = TRUE)
#'       }
#'       
#'       if (!is.finite(mean_count_scaling_max) || mean_count_scaling_max <= 0) {
#'         mean_count_scaling_max <- 1
#'       }
#'       
#'       circle_df <- sp_hex_grid %>%
#'         sf::st_centroid() %>%
#'         dplyr::mutate(
#'           x = sf::st_coordinates(.)[, 1],
#'           y = sf::st_coordinates(.)[, 2],
#'           
#'           survey_alpha = alpha_min +
#'             (pmin(n_surveys, survey_scaling_max) / survey_scaling_max) *
#'             (alpha_max - alpha_min),
#'           
#'           detection_alpha = pmax(survey_alpha, 0.4),
#'           
#'           r_survey = survey_radius,
#'           
#'           mean_count_scaled = pmin(mean_count, mean_count_scaling_max) /
#'             mean_count_scaling_max,
#'           
#'           r_detection = dplyr::if_else(
#'             mean_count > 0,
#'             min_detection_radius +
#'               sqrt(mean_count_scaled) *
#'               (max_detection_radius - min_detection_radius),
#'             0
#'           )
#'         ) %>%
#'         sf::st_drop_geometry()
#'     }
#'   }
#'   
#'   # ------------------------------------------------------------
#'   # Build plot
#'   # ------------------------------------------------------------
#'   bb <- map_bbox(study_boundary)
#'   
#'   p <- ggplot2::ggplot()
#'   
#'   if (!is.null(r_stars)) {
#'     p <- p +
#'       stars::geom_stars(data = r_stars) +
#'       ggplot2::scale_fill_gradientn(
#'         name = legend_title_map(species_name, title, subtitle, "Posterior median"),
#'         colors = q50_cols,
#'         na.value = "transparent",
#'         limits = zlim
#'       )
#'   }
#'   
#'   p <- add_base_map_layers(
#'     p = p,
#'     study_boundary = study_boundary,
#'     bcr_sf = bcr_sf,
#'     water_sf = water_sf,
#'     atlas_squares_centroids = NULL,
#'     water_fill = water_fill
#'   )
#'   
#'   if (!is.null(circle_df) && nrow(circle_df) > 0) {
#'     p <- p +
#'       ggforce::geom_circle(
#'         data = circle_df,
#'         ggplot2::aes(
#'           x0 = x,
#'           y0 = y,
#'           r = r_survey,
#'           alpha = survey_alpha
#'         ),
#'         fill = survey_colour,
#'         color = NA
#'       ) +
#'       ggforce::geom_circle(
#'         data = dplyr::filter(circle_df, mean_count > 0),
#'         ggplot2::aes(
#'           x0 = x,
#'           y0 = y,
#'           r = r_detection,
#'           alpha = detection_alpha
#'         ),
#'         fill = detection_colour,
#'         color = NA
#'       ) +
#'       ggplot2::scale_alpha_identity()
#'   }
#'   
#'   p <- add_map_annotations(p, study_boundary, bb)
#'   
#'   p <- p +
#'     ggplot2::theme(
#'       panel.background = ggplot2::element_rect(fill = "white", colour = NA),
#'       plot.background  = ggplot2::element_rect(fill = "white", colour = NA),
#'       panel.grid       = ggplot2::element_blank()
#'     )
#'   
#'   list(
#'     sp_plot = p,
#'     zlim = zlim
#'   )
#' }