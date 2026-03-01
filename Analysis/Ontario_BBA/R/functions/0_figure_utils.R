# ============================================================
# figure_utils.R
#
# Purpose:
#   Utilities for making publication-quality maps from
#   per-pixel prediction summaries.
#
# Philosophy:
#   - Functions RETURN ggplot objects (do not ggsave()).
#   - Caller script controls filenames, dpi, sizes.
#   - Minimal, consistent styling + a few sensible defaults.
# ============================================================

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(ggplot2)
  library(stars)
  library(terra)
  library(ggtext)
  library(ggspatial)
  library(RColorBrewer)
})

# ---------------------------
# Small helpers
# ---------------------------

wrap_species_label <- function(label, max_length = 18) {
  if (nchar(label) <= max_length) return(label)
  words <- strsplit(label, " ")[[1]]
  if (length(words) == 1) return(label)
  paste0(paste(words[-length(words)], collapse = " "), "<br>", words[length(words)])
}

expand_bbox <- function(bbox, frac = 0.10) {
  bbox <- as.numeric(bbox)
  names(bbox) <- c("xmin","ymin","xmax","ymax")  # ensure names exist
  
  xrange <- bbox["xmax"] - bbox["xmin"]
  yrange <- bbox["ymax"] - bbox["ymin"]
  
  bbox["xmin"] <- bbox["xmin"] - xrange * frac
  bbox["xmax"] <- bbox["xmax"] + xrange * frac
  bbox["ymin"] <- bbox["ymin"] - yrange * frac
  bbox["ymax"] <- bbox["ymax"] + yrange * frac
  
  bbox
}

# Rasterize an sf grid using a column name; returns stars object
rasterize_sf_to_stars <- function(grid_sf, field, res) {
  stopifnot(inherits(grid_sf, "sf"))
  stopifnot(field %in% names(grid_sf))
  
  v <- terra::vect(grid_sf)
  r_template <- terra::rast(v, res = res)
  r <- terra::rasterize(v, r_template, field = field, fun = mean, na.rm = TRUE)
  stars::st_as_stars(r)
}

# Build a consistent "map theme"
theme_map <- function() {
  theme_void() +
    theme(
      panel.background = element_rect(fill = "#F5F5F5", color = NA),
      plot.margin = unit(c(0, 0, 0, 0), "pt"),
      legend.title = ggtext::element_markdown(lineheight = 0.95),
      legend.position = c(0.02, 0.02),
      legend.justification = c(0, 0),
      legend.background = element_rect(fill = "white", color = "black"),
      legend.margin = margin(5, 5, 5, 5),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
    )
}

# ---------------------------
# Relative abundance maps
# ---------------------------

#' Build relabund maps (median + CV)
#'
#' @param grid_sf sf polygons for the atlas period (OBBA2 or OBBA3)
#' @param pred_summary data.frame with columns:
#'   - <prefix>_q50 (required)
#'   - <prefix>_cv_median (preferred) OR <prefix>_sd + <prefix>_q50 (fallback)
#' @param prefix "OBBA2" or "OBBA3"
#' @param res raster resolution for plotting (in map units)
#' @param bounds list with optional upper/lower; if NULL it is computed from q50
#' @return named list: list(q50_plot = <ggplot>, cv_plot = <ggplot>)
make_relabund_maps <- function(species_name,
                               grid_sf,
                               pred_summary,
                               prefix = c("OBBA2", "OBBA3"),
                               study_boundary,
                               bcr_sf = NULL,
                               water_sf = NULL,
                               atlas_squares_centroids = NULL,
                               title = NULL,
                               subtitle = "Relative Abundance",
                               res = 1000,
                               bounds = NULL) {
  
  prefix <- match.arg(prefix)
  stopifnot(inherits(grid_sf, "sf"))
  stopifnot(inherits(study_boundary, "sf"))
  
  q50_col <- paste0(prefix, "_q50")
  if (!(q50_col %in% names(pred_summary))) {
    stop("pred_summary missing required column: ", q50_col)
  }
  
  # attach q50
  grid <- grid_sf
  grid$pred_q50 <- pred_summary[[q50_col]]
  
  # attach CV (prefer cv_median; else compute sd/q50)
  cv_col <- paste0(prefix, "_cv_median")
  if (cv_col %in% names(pred_summary)) {
    grid$pred_cv <- pred_summary[[cv_col]]
  } else {
    sd_col <- paste0(prefix, "_sd")
    if (!(sd_col %in% names(pred_summary))) {
      stop("Need either ", cv_col, " or ", sd_col, " to make CV map.")
    }
    grid$pred_cv <- pred_summary[[sd_col]] / pmax(pred_summary[[q50_col]], 1e-12)
  }
  
  # bounds
  # If caller supplies bounds, we use them as-is (so OBBA2 and OBBA3 can share
  # identical colour limits for side-by-side comparison). Otherwise compute a
  # sensible default from the 99th percentile of the median surface.
  if (is.null(bounds)) {
    upper <- as.numeric(stats::quantile(grid$pred_q50, 0.99, na.rm = TRUE))
    if (!is.finite(upper) || upper <= 0) upper <- 1
    bounds <- list(lower = 0, upper = upper)
  } else {
    if (is.null(bounds$lower)) bounds$lower <- 0
    if (is.null(bounds$upper)) {
      bounds$upper <- as.numeric(stats::quantile(grid$pred_q50, 0.99, na.rm = TRUE))
    }
  }
  
  # cap for plotting (keeps visual scale stable)
  grid$pred_capped_q50 <- as.numeric(pmax(pmin(grid$pred_q50, bounds$upper), bounds$lower))
  grid$pred_capped_cv  <- as.numeric(pmax(pmin(grid$pred_cv, 1), 0))
  
  # rasterize to stars
  q50_stars <- rasterize_sf_to_stars(grid, "pred_capped_q50", res = res)
  cv_stars  <- rasterize_sf_to_stars(grid, "pred_capped_cv",  res = res)
  
  # bbox
  bb <- expand_bbox(st_bbox(study_boundary), 0.10)
  
  # palette
  colscale_q50 <- c(
    "#FEFEFE", "#FBF7E2", "#FCF8D0", "#EEF7C2", "#CEF2B0",
    "#94E5A0", "#51C987", "#18A065", "#008C59", "#007F53", "#006344"
  )
  colpal_q50 <- colorRampPalette(colscale_q50)
  
  colscale_cv <- c("#FEFEFE", "#FFF4B3", "#F5D271", "#F2B647", "#EC8E00", "#CA302A")
  colpal_cv <- colorRampPalette(colscale_cv)
  
  if (is.null(title)) title <- ifelse(prefix == "OBBA2", "Atlas 2", "Atlas 3")
  
  legend_title_q50 <- paste0(
    "<span style='font-size:20pt; font-weight:bold'>", wrap_species_label(species_name), "</span><br><br>",
    "<span style='font-size:14pt'>", title, "</span><br>",
    "<span style='font-size:7pt'>", subtitle, "</span><br>",
    "<span style='font-size:7pt'>Posterior Median</span>"
  )
  
  legend_title_cv <- paste0(
    "<span style='font-size:20pt; font-weight:bold'>", wrap_species_label(species_name), "</span><br><br>",
    "<span style='font-size:14pt'>", title, "</span><br>",
    "<span style='font-size:7pt'>", subtitle, "</span><br>",
    "<span style='font-size:7pt'>Prediction CV</span>"
  )
  
  base_layers <- list()
  if (!is.null(bcr_sf))   base_layers <- c(base_layers, list(geom_sf(data = bcr_sf, colour = "gray80", fill = NA, linewidth = 0.3)))
  base_layers <- c(base_layers, list(stars::geom_stars(data = q50_stars)))
  if (!is.null(water_sf)) base_layers <- c(base_layers, list(geom_sf(data = water_sf, fill = "#EDF7FB", col = "transparent")))
  if (!is.null(atlas_squares_centroids)) {
    # If the caller doesn't supply a detected flag, assume all TRUE.
    if (!("detected" %in% names(atlas_squares_centroids))) atlas_squares_centroids$detected <- TRUE
    base_layers <- c(base_layers, list(
      geom_sf(
        data = atlas_squares_centroids,
        aes(colour = detected),
        shape = 46,
        size = 0.2,
        alpha = 0.7
      )
    ))
  }
  base_layers <- c(base_layers, list(geom_sf(data = study_boundary, colour = "black", fill = NA, linewidth = 0.5, show.legend = FALSE)))
  
  q50_plot <- ggplot() +
    base_layers +
    scale_fill_gradientn(
      name = legend_title_q50,
      colors = colpal_q50(10),
      na.value = "transparent",
      limits = c(bounds$lower, bounds$upper)
    ) +
    {if (!is.null(atlas_squares_centroids)) scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "gray60"), guide = "none")} +
    theme_map() +
    ggspatial::annotation_scale(location = "br", width_hint = 0.3) +
    ggspatial::annotation_north_arrow(
      location = "tr",
      which_north = "true",
      pad_x = unit(0.2, "in"),
      pad_y = unit(0.2, "in"),
      style = ggspatial::north_arrow_fancy_orienteering()
    ) +
    coord_sf(
      crs = sf::st_crs(study_boundary),
      xlim = c(bb["xmin"], bb["xmax"]),
      ylim = c(bb["ymin"], bb["ymax"]),
      expand = FALSE
    )
  
  cv_plot <- ggplot() +
    {if (!is.null(bcr_sf)) geom_sf(data = bcr_sf, colour = "gray80", fill = NA, linewidth = 0.3)} +
    stars::geom_stars(data = cv_stars) +
    scale_fill_gradientn(
      name = legend_title_cv,
      colors = colpal_cv(11),
      na.value = "transparent",
      breaks = seq(0, 1, length.out = 5),
      labels = c("0", "0.25", "0.5", "0.75", ">1"),
      limits = c(0, 1)
    ) +
    {if (!is.null(water_sf)) geom_sf(data = water_sf, fill = "#EDF7FB", col = "transparent")} +
    {if (!is.null(atlas_squares_centroids)) geom_sf(data = atlas_squares_centroids, colour = "black", size = 0.2, shape = 46, stroke = 0.1, alpha = 0.5)} +
    theme_map() +
    ggspatial::annotation_scale(location = "br", width_hint = 0.3) +
    ggspatial::annotation_north_arrow(
      location = "tr",
      which_north = "true",
      pad_x = unit(0.2, "in"),
      pad_y = unit(0.2, "in"),
      style = ggspatial::north_arrow_fancy_orienteering()
    )+
    coord_sf(
      crs = sf::st_crs(study_boundary),
      xlim = c(bb["xmin"], bb["xmax"]),
      ylim = c(bb["ymin"], bb["ymax"]),
      expand = FALSE
    )
  
  list(q50_plot = q50_plot, cv_plot = cv_plot, bounds = bounds)
}

# ---------------------------
# Change maps (absolute change + CI width)
# ---------------------------

make_abs_change_maps <- function(species_name,
                                 grid_sf,               # usually OBBA3 grid (for geometry)
                                 abs_change_summary,    # from predictions/<sp>.rds
                                 study_boundary,
                                 bcr_sf = NULL,
                                 water_sf = NULL,
                                 signif_poly = NULL,
                                 res = 1000,
                                 # Optional: force symmetric bounds for change map.
                                 # If provided, colour limits will be -max_abs..+max_abs.
                                 max_abs = NULL,
                                 title = "Absolute change",
                                 subtitle = "OBBA2 to OBBA3") {
  
  stopifnot(inherits(grid_sf, "sf"))
  stopifnot(is.data.frame(abs_change_summary))
  
  q50_col <- "abs_change_q50"
  lcl_col <- "abs_change_lower"
  ucl_col <- "abs_change_upper"
  
  if (!all(c(q50_col, lcl_col, ucl_col) %in% names(abs_change_summary))) {
    stop("abs_change_summary must contain: abs_change_q50, abs_change_lower, abs_change_upper")
  }
  
  grid <- grid_sf
  grid$chg_q50 <- abs_change_summary[[q50_col]]
  grid$chg_ciw <- abs_change_summary[[ucl_col]] - abs_change_summary[[lcl_col]]
  
  # symmetric bounds around 0
  # If caller supplies max_abs, we use it (e.g., to match relabund colour scale).
  if (is.null(max_abs)) {
    max_abs <- as.numeric(stats::quantile(abs(grid$chg_q50), 0.99, na.rm = TRUE))
    if (!is.finite(max_abs) || max_abs == 0) max_abs <- max(abs(grid$chg_q50), na.rm = TRUE)
    if (!is.finite(max_abs) || max_abs == 0) max_abs <- 1
  } else {
    max_abs <- as.numeric(max_abs)
    if (!is.finite(max_abs) || max_abs <= 0) max_abs <- 1
  }
  
  grid$chg_capped <- pmax(pmin(grid$chg_q50,  max_abs), -max_abs)
  grid$ciw_capped <- pmax(pmin(grid$chg_ciw,  max_abs), 0)
  
  chg_stars <- rasterize_sf_to_stars(grid, "chg_capped", res = res)
  ciw_stars <- rasterize_sf_to_stars(grid, "ciw_capped", res = res)
  
  bb <- expand_bbox(st_bbox(study_boundary), 0.10)
  
  colscale_chg <- RColorBrewer::brewer.pal(11, "RdBu")
  colpal_chg <- colorRampPalette(colscale_chg)
  
  colscale_unc <- c("#FEFEFE", "#FFF4B3", "#F5D271", "#F2B647", "#EC8E00", "#CA302A")
  colpal_unc <- colorRampPalette(colscale_unc)
  
  legend_title_chg <- paste0(
    "<span style='font-size:20pt; font-weight:bold'>", wrap_species_label(species_name), "</span><br><br>",
    "<span style='font-size:14pt'>", title, "</span><br>",
    "<span style='font-size:7pt'>", subtitle, "</span><br>",
    "<span style='font-size:7pt'>Posterior Median</span>"
  )
  
  legend_title_ciw <- paste0(
    "<span style='font-size:20pt; font-weight:bold'>", wrap_species_label(species_name), "</span><br><br>",
    "<span style='font-size:14pt'>", title, "</span><br>",
    "<span style='font-size:7pt'>", subtitle, "</span><br>",
    "<span style='font-size:7pt'>90% CI width</span>"
  )
  
  breaks <- seq(-max_abs, max_abs, length.out = 7)
  breaks[4] <- 0
  
  chg_plot <- ggplot() +
    (if (!is.null(bcr_sf)) geom_sf(data = bcr_sf, colour = "gray80", fill = NA, linewidth = 0.3) else NULL) +
    stars::geom_stars(data = chg_stars) +
    scale_fill_gradientn(
      name = legend_title_chg,
      colors = colpal_chg(11),
      na.value = "transparent",
      breaks = breaks,
      labels = signif(breaks, 2),
      limits = c(-max_abs, max_abs)
    ) +
    (if (!is.null(water_sf)) geom_sf(data = water_sf, fill = "#F5F5F5", col = "transparent") else NULL) +
    (if (!is.null(signif_poly)) list(
      geom_sf(data = subset(signif_poly, signif_change == "Increase"),
              fill = "transparent", col = "dodgerblue", linewidth = 1),
      geom_sf(data = subset(signif_poly, signif_change == "Decrease"),
              fill = "transparent", col = "orangered", linewidth = 1)
    ) else NULL) +
    geom_sf(data = study_boundary, colour = "black", fill = NA, linewidth = 0.5, show.legend = FALSE) +
    theme_map() +
    ggspatial::annotation_scale(location = "br", width_hint = 0.3) +
    ggspatial::annotation_north_arrow(
      location = "tr",
      which_north = "true",
      pad_x = unit(0.2, "in"),
      pad_y = unit(0.2, "in"),
      style = ggspatial::north_arrow_fancy_orienteering()
    ) +
    coord_sf(
      crs = sf::st_crs(study_boundary),
      xlim = c(bb["xmin"], bb["xmax"]),
      ylim = c(bb["ymin"], bb["ymax"]),
      expand = FALSE
    )
  
  ciw_breaks <- seq(0, max_abs, length.out = 7)
  
  ciw_plot <- ggplot() +
    {if (!is.null(bcr_sf)) geom_sf(data = bcr_sf, colour = "gray80", fill = NA, linewidth = 0.3)} +
    stars::geom_stars(data = ciw_stars) +
    scale_fill_gradientn(
      name = legend_title_ciw,
      colors = colpal_unc(11),
      na.value = "transparent",
      breaks = ciw_breaks,
      labels = signif(ciw_breaks, 2),
      limits = c(0, max_abs)
    ) +
    {if (!is.null(water_sf)) geom_sf(data = water_sf, fill = "#F5F5F5", col = "transparent")} +
    geom_sf(data = study_boundary, colour = "black", fill = NA, linewidth = 0.5, show.legend = FALSE) +
    theme_map() +
    ggspatial::annotation_scale(location = "br", width_hint = 0.3) +
    ggspatial::annotation_north_arrow(
      location = "tr",
      which_north = "true",
      pad_x = unit(0.2, "in"),
      pad_y = unit(0.2, "in"),
      style = ggspatial::north_arrow_fancy_orienteering()
    ) +
    coord_sf(
      crs = sf::st_crs(study_boundary),
      xlim = c(bb["xmin"], bb["xmax"]),
      ylim = c(bb["ymin"], bb["ymax"]),
      expand = FALSE
    )
  
  list(chg_plot = chg_plot, ciw_plot = ciw_plot, max_abs = max_abs)
}
