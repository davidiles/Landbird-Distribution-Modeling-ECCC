# ============================================================
# figure_utils.R
#
# Utilities for making maps from pixel-level prediction summaries.
# Functions return ggplot objects; saving is handled elsewhere.
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
  library(igraph)
  library(smoothr)
  library(units)
})

# ------------------------------------------------------------
# Small helpers
# ------------------------------------------------------------

# compute abundance limits for plotting
compute_shared_zmax <- function(preds, prob = 0.99, fallback = 1) {
  zmax2 <- as.numeric(stats::quantile(preds$OBBA2$OBBA2_q50, prob, na.rm = TRUE))
  zmax3 <- as.numeric(stats::quantile(preds$OBBA3$OBBA3_q50, prob, na.rm = TRUE))
  zmax <- max(zmax2, zmax3, na.rm = TRUE)
  if (!is.finite(zmax) || zmax <= 0) zmax <- fallback
  zmax
}

# build square overlays 
build_atlas_square_overlays <- function(sp_square_summary, atlas_sq_centroids_all) {
  if (is.null(sp_square_summary)) {
    return(list(OBBA2 = NULL, OBBA3 = NULL))
  }
  
  sq2 <- dplyr::filter(sp_square_summary, Atlas == "OBBA2")
  sq3 <- dplyr::filter(sp_square_summary, Atlas == "OBBA3")
  
  list(
    OBBA2 = dplyr::left_join(atlas_sq_centroids_all, sq2, by = "square_id"),
    OBBA3 = dplyr::left_join(atlas_sq_centroids_all, sq3, by = "square_id")
  )
}


#' Wrap a species label for nicer plot titles
#'
#' @param label Character species name.
#' @param max_length Maximum characters before wrapping.
#' @return Character string, optionally with "<br>" inserted.
wrap_species_label <- function(label, max_length = 18) {
  if (is.na(label) || nchar(label) <= max_length) return(label)
  
  words <- strsplit(label, " ")[[1]]
  if (length(words) <= 1) return(label)
  
  paste0(
    paste(words[-length(words)], collapse = " "),
    "<br>",
    words[length(words)]
  )
}

#' Expand a bounding box by a fraction
#'
#' @param bbox An sf bbox object.
#' @param frac Fractional expansion.
#' @return Named numeric vector xmin/ymin/xmax/ymax.
expand_bbox <- function(bbox, frac = 0.10) {
  bbox <- as.numeric(bbox)
  names(bbox) <- c("xmin", "ymin", "xmax", "ymax")
  
  xrange <- bbox["xmax"] - bbox["xmin"]
  yrange <- bbox["ymax"] - bbox["ymin"]
  
  bbox["xmin"] <- bbox["xmin"] - frac * xrange
  bbox["xmax"] <- bbox["xmax"] + frac * xrange
  bbox["ymin"] <- bbox["ymin"] - frac * yrange
  bbox["ymax"] <- bbox["ymax"] + frac * yrange
  
  bbox
}

#' Get a padded plotting bbox from a study boundary
#'
#' @param study_boundary sf polygon.
#' @param frac Fractional expansion.
#' @return Named numeric vector xmin/ymin/xmax/ymax.
map_bbox <- function(study_boundary, frac = 0.10) {
  expand_bbox(sf::st_bbox(study_boundary), frac = frac)
}

#' Rasterize an sf grid to a stars object
#'
#' @param grid_sf sf object with geometry and a value field.
#' @param field Name of field to rasterize.
#' @param res Raster resolution in map units.
#' @return stars object.
rasterize_sf_to_stars <- function(grid_sf, field, res) {
  stopifnot(inherits(grid_sf, "sf"))
  stopifnot(field %in% names(grid_sf))
  
  v <- terra::vect(grid_sf)
  r_template <- terra::rast(v, res = res)
  r <- terra::rasterize(v, r_template, field = field, fun = mean, na.rm = TRUE)
  
  stars::st_as_stars(r)
}

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

#' Standard project map theme
#'
#' @return ggplot theme object.
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

#' Build a formatted legend title
#'
#' @param species_name Character species common name.
#' @param title Main panel title.
#' @param subtitle Secondary title.
#' @param stat_label Statistic label.
#' @return Character string containing markdown/html for ggtext legend title.
legend_title_map <- function(species_name, title, subtitle, stat_label) {
  paste0(
    "<span style='font-size:20pt; font-weight:bold'>",
    wrap_species_label(species_name),
    "</span><br><br>",
    "<span style='font-size:14pt'>", title, "</span><br>",
    "<span style='font-size:7pt'>", subtitle, "</span><br>",
    "<span style='font-size:7pt'>", stat_label, "</span>"
  )
}

#' Add common overlays to a map
#'
#' @param p ggplot object.
#' @param study_boundary sf boundary polygon.
#' @param bcr_sf Optional sf BCR polygons.
#' @param water_sf Optional sf water polygons.
#' @param atlas_squares_centroids Optional sf square centroids.
#' @param water_fill Fill for water layer.
#'   "plain" for black points, or "none" for no square overlay.
#' @return ggplot object.
add_base_map_layers <- function(p,
                                study_boundary,
                                bcr_sf = NULL,
                                water_sf = NULL,
                                atlas_squares_centroids = NULL,
                                water_fill = "#EDF7FB") {
  
  if (!is.null(bcr_sf)) {
    p <- p + geom_sf(
      data = bcr_sf,
      colour = "gray80",
      fill = NA,
      linewidth = 0.15
    )
  }
  
  
  if (!is.null(water_sf)) {
    p <- p + geom_sf(
      data = water_sf,
      fill = water_fill,
      col = "transparent"
    )
  }
  
  p <- p + geom_sf(
    data = study_boundary,
    colour = "gray80",
    fill = NA,
    linewidth = 0.15,
    show.legend = FALSE
  )
  
  if (!is.null(atlas_squares_centroids)) {
    
    sq <- atlas_squares_centroids %>%
      dplyr::filter(!is.na(n_detections))
    
    if (nrow(sq) > 0) {
      sq <- sq %>%
        dplyr::mutate(detected = n_detections > 0)
      
      # surveyed, not detected
      p <- p + geom_sf(
        data = sq, # dplyr::filter(sq, !detected),
        colour = "black",
        shape = 1,
        size = 0.5,
        alpha = 1,
        stroke = 0.1
      )
      
      # surveyed, detected
      p <- p + geom_sf(
        data = dplyr::filter(sq, detected),
        colour = "black",
        shape = 16,   #46,
        size = 0.5, # 0.75,
        stroke = 0.1,
        alpha = 1
      )
    }
  }
  p
}

#' Add common annotations and map extent
#'
#' @param p ggplot object.
#' @param study_boundary sf boundary polygon.
#' @param bb Optional bbox. If NULL, computed from study_boundary.
#' @return ggplot object.
add_map_annotations <- function(p, study_boundary, bb = NULL) {
  if (is.null(bb)) {
    bb <- map_bbox(study_boundary)
  }
  
  p +
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
}

#' Compute a shared upper plotting limit for OBBA2 and OBBA3 q50 surfaces
#'
#' @param preds Prediction object containing OBBA2 and OBBA3 summaries.
#' @param prob Quantile used for upper cap.
#' @param fallback Fallback value if cap is invalid.
#' @return Positive numeric scalar.
compute_shared_zmax <- function(preds, prob = 0.99, fallback = 1) {
  zmax2 <- as.numeric(stats::quantile(preds$OBBA2$OBBA2_q50, prob, na.rm = TRUE))
  zmax3 <- as.numeric(stats::quantile(preds$OBBA3$OBBA3_q50, prob, na.rm = TRUE))
  zmax <- max(zmax2, zmax3, na.rm = TRUE)
  
  if (!is.finite(zmax) || zmax <= 0) zmax <- fallback
  zmax
}

#' Build Atlas 2 / Atlas 3 atlas-square centroid overlays
#'
#' @param sp_square_summary Data frame with square summaries and Atlas column.
#' @param atlas_sq_centroids_all sf centroid layer with square_id.
#' @return Named list with OBBA2 and OBBA3 sf objects (or NULL values).
build_atlas_square_overlays <- function(sp_square_summary, atlas_sq_centroids_all) {
  if (is.null(sp_square_summary)) {
    return(list(OBBA2 = NULL, OBBA3 = NULL))
  }
  
  sq2 <- dplyr::filter(sp_square_summary, Atlas == "OBBA2")
  sq3 <- dplyr::filter(sp_square_summary, Atlas == "OBBA3")
  
  list(
    OBBA2 = dplyr::left_join(atlas_sq_centroids_all, sq2, by = "square_id"),
    OBBA3 = dplyr::left_join(atlas_sq_centroids_all, sq3, by = "square_id")
  )
}

# ------------------------------------------------------------
# Relative abundance maps
# ------------------------------------------------------------

#' Create relative-abundance maps for one atlas period
#'
#' Returns a posterior-median map and a CV map.
#'
#' @param species_name Common name for legend title.
#' @param grid_sf sf prediction grid.
#' @param pred_summary Data frame aligned to grid_sf.
#' @param prefix Either "OBBA2" or "OBBA3".
#' @param study_boundary sf boundary polygon.
#' @param bcr_sf Optional sf BCR overlay.
#' @param water_sf Optional sf water overlay.
#' @param atlas_squares_centroids Optional sf square-centroid overlay.
#' @param title Panel title. If NULL, inferred from prefix.
#' @param subtitle Subtitle for legend title.
#' @param res Raster resolution in map units.
#' @param bounds Optional list with lower and upper abundance bounds.
#' @return List with q50_plot, cv_plot, and bounds.
make_relabund_maps <- function(species_name,
                               grid_sf,
                               pred_summary,
                               prefix = c("OBBA2", "OBBA3"),
                               study_boundary,
                               bcr_sf = NULL,
                               water_sf = NULL,
                               atlas_squares_centroids = NULL,
                               title = NULL,
                               subtitle = "Relative abundance",
                               res = 1000,
                               bounds = NULL) {
  prefix <- match.arg(prefix)
  
  stopifnot(inherits(grid_sf, "sf"))
  stopifnot(inherits(study_boundary, "sf"))
  stopifnot(is.data.frame(pred_summary))
  
  q50_col <- paste0(prefix, "_q50")
  cv_col  <- paste0(prefix, "_cv_median")
  sd_col  <- paste0(prefix, "_sd")
  
  if (!(q50_col %in% names(pred_summary))) {
    stop("pred_summary is missing required column: ", q50_col)
  }
  
  grid <- grid_sf
  grid$pred_q50 <- pred_summary[[q50_col]]
  
  if (cv_col %in% names(pred_summary)) {
    grid$pred_cv <- pred_summary[[cv_col]]
  } else if (sd_col %in% names(pred_summary)) {
    grid$pred_cv <- pred_summary[[sd_col]] / pmax(pred_summary[[q50_col]], 1e-12)
  } else {
    stop("pred_summary must contain either ", cv_col, " or ", sd_col)
  }
  
  if (is.null(bounds)) {
    upper <- as.numeric(stats::quantile(grid$pred_q50, 0.99, na.rm = TRUE))
    if (!is.finite(upper) || upper <= 0) upper <- 1
    bounds <- list(lower = 0, upper = upper)
  } else {
    if (is.null(bounds$lower)) bounds$lower <- 0
    if (is.null(bounds$upper)) {
      bounds$upper <- as.numeric(stats::quantile(grid$pred_q50, 0.99, na.rm = TRUE))
      if (!is.finite(bounds$upper) || bounds$upper <= 0) bounds$upper <- 1
    }
  }
  
  grid$pred_capped_q50 <- pmax(pmin(grid$pred_q50, bounds$upper), bounds$lower)
  grid$pred_capped_cv  <- pmax(pmin(grid$pred_cv, 1), 0)
  
  q50_stars <- rasterize_sf_to_stars(grid, "pred_capped_q50", res = res)
  cv_stars  <- rasterize_sf_to_stars(grid, "pred_capped_cv",  res = res)
  
  if (is.null(title)) {
    title <- ifelse(prefix == "OBBA2", "Atlas 2", "Atlas 3")
  }
  
  bb <- map_bbox(study_boundary)
  
  q50_cols <- colorRampPalette(c(
    "#FEFEFE", "#FBF7E2", "#FCF8D0", "#EEF7C2", "#CEF2B0",
    "#94E5A0", "#51C987", "#18A065", "#008C59", "#007F53", "#006344"
  ))(10)
  
  cv_cols <- colorRampPalette(c(
    "#FEFEFE", "#FFF4B3", "#F5D271", "#F2B647", "#EC8E00", "#CA302A"
  ))(11)
  
  q50_plot <- ggplot() +
    stars::geom_stars(data = q50_stars)
  
  q50_plot <- add_base_map_layers(
    p = q50_plot,
    study_boundary = study_boundary,
    bcr_sf = bcr_sf,
    water_sf = water_sf,
    atlas_squares_centroids = atlas_squares_centroids,
    water_fill = "#EDF7FB"
  )
  
  q50_plot <- q50_plot +
    scale_fill_gradientn(
      name = legend_title_map(species_name, title, subtitle, "Posterior median"),
      colors = q50_cols,
      na.value = "transparent",
      limits = c(bounds$lower, bounds$upper)
    )
  
  q50_plot <- add_map_annotations(q50_plot, study_boundary, bb)
  
  cv_plot <- ggplot() +
    stars::geom_stars(data = cv_stars)
  
  cv_plot <- add_base_map_layers(
    p = cv_plot,
    study_boundary = study_boundary,
    bcr_sf = bcr_sf,
    water_sf = water_sf,
    atlas_squares_centroids = atlas_squares_centroids,
    water_fill = "#EDF7FB"
  )
  
  cv_plot <- cv_plot +
    scale_fill_gradientn(
      name = legend_title_map(species_name, title, subtitle, "Prediction CV"),
      colors = cv_cols,
      na.value = "transparent",
      breaks = seq(0, 1, length.out = 5),
      labels = c("0", "0.25", "0.5", "0.75", ">1"),
      limits = c(0, 1)
    )
  
  cv_plot <- add_map_annotations(cv_plot, study_boundary, bb)
  
  list(
    q50_plot = q50_plot,
    cv_plot = cv_plot,
    bounds = bounds
  )
}

# ------------------------------------------------------------
# Absolute change maps
# ------------------------------------------------------------

#' Create absolute-change maps between atlas periods
#'
#' Returns a posterior-median absolute-change map and a CI-width map.
#'
#' @param species_name Common name for legend title.
#' @param grid_sf sf prediction grid (typically the OBBA3 grid geometry).
#' @param abs_change_summary Data frame aligned to grid_sf.
#' @param study_boundary sf boundary polygon.
#' @param bcr_sf Optional sf BCR overlay.
#' @param water_sf Optional sf water overlay.
#' @param signif_poly Optional sf polygons classed as Increase/Decrease.
#' @param res Raster resolution in map units.
#' @param max_abs Optional symmetric upper plotting bound for change.
#' @param title Panel title.
#' @param subtitle Subtitle for legend title.
#' @return List with chg_plot, ciw_plot, and max_abs.
make_abs_change_maps <- function(species_name,
                                 grid_sf,
                                 abs_change_summary,
                                 study_boundary,
                                 bcr_sf = NULL,
                                 water_sf = NULL,
                                 signif_poly = NULL,
                                 res = 1000,
                                 max_abs = NULL,
                                 title = "Absolute change",
                                 subtitle = "OBBA2 to OBBA3") {
  stopifnot(inherits(grid_sf, "sf"))
  stopifnot(inherits(study_boundary, "sf"))
  stopifnot(is.data.frame(abs_change_summary))
  
  required_cols <- c("abs_change_q50", "abs_change_lower", "abs_change_upper")
  if (!all(required_cols %in% names(abs_change_summary))) {
    stop(
      "abs_change_summary must contain: ",
      paste(required_cols, collapse = ", ")
    )
  }
  
  grid <- grid_sf
  grid$chg_q50 <- abs_change_summary$abs_change_q50
  grid$chg_ciw <- abs_change_summary$abs_change_upper - abs_change_summary$abs_change_lower
  
  if (is.null(max_abs)) {
    max_abs <- as.numeric(stats::quantile(abs(grid$chg_q50), 0.99, na.rm = TRUE))
    
    if (!is.finite(max_abs) || max_abs <= 0) {
      max_abs <- max(abs(grid$chg_q50), na.rm = TRUE)
    }
    
    if (!is.finite(max_abs) || max_abs <= 0) {
      max_abs <- 1
    }
  } else {
    max_abs <- as.numeric(max_abs)
    if (!is.finite(max_abs) || max_abs <= 0) max_abs <- 1
  }
  
  grid$chg_capped <- pmax(pmin(grid$chg_q50, max_abs), -max_abs)
  grid$ciw_capped <- pmax(pmin(grid$chg_ciw, max_abs), 0)
  
  chg_stars <- rasterize_sf_to_stars(grid, "chg_capped", res = res)
  ciw_stars <- rasterize_sf_to_stars(grid, "ciw_capped", res = res)
  
  bb <- map_bbox(study_boundary)
  
  chg_cols <- colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(11)
  unc_cols <- colorRampPalette(c(
    "#FEFEFE", "#FFF4B3", "#F5D271", "#F2B647", "#EC8E00", "#CA302A"
  ))(11)
  
  chg_breaks <- seq(-max_abs, max_abs, length.out = 7)
  chg_breaks[4] <- 0
  ciw_breaks <- seq(0, max_abs, length.out = 7)
  
  chg_plot <- ggplot() +
    stars::geom_stars(data = chg_stars)
  
  chg_plot <- add_base_map_layers(
    p = chg_plot,
    study_boundary = study_boundary,
    bcr_sf = bcr_sf,
    water_sf = water_sf,
    atlas_squares_centroids = NULL,
    water_fill = "#F5F5F5"
  )
  
  chg_plot <- chg_plot +
    scale_fill_gradientn(
      name = legend_title_map(species_name, title, subtitle, "Posterior median"),
      colors = chg_cols,
      na.value = "transparent",
      breaks = chg_breaks,
      labels = signif(chg_breaks, 2),
      limits = c(-max_abs, max_abs)
    )
  
  if (!is.null(signif_poly) && nrow(signif_poly) > 0) {
    if ("classification" %in% names(signif_poly)) {
      if (any(signif_poly$classification == "Increase")) {
        chg_plot <- chg_plot +
          geom_sf(
            data = subset(signif_poly, classification == "Increase"),
            fill = "transparent",
            col = "dodgerblue",
            linewidth = 1
          )
      }
      
      if (any(signif_poly$classification == "Decrease")) {
        chg_plot <- chg_plot +
          geom_sf(
            data = subset(signif_poly, classification == "Decrease"),
            fill = "transparent",
            col = "orangered",
            linewidth = 1
          )
      }
    }
  }
  
  chg_plot <- add_map_annotations(chg_plot, study_boundary, bb)
  
  ciw_plot <- ggplot() +
    stars::geom_stars(data = ciw_stars)
  
  ciw_plot <- add_base_map_layers(
    p = ciw_plot,
    study_boundary = study_boundary,
    bcr_sf = bcr_sf,
    water_sf = water_sf,
    atlas_squares_centroids = NULL,
    water_fill = "#F5F5F5"
  )
  
  ciw_plot <- ciw_plot +
    scale_fill_gradientn(
      name = legend_title_map(species_name, title, subtitle, "90% CI width"),
      colors = unc_cols,
      na.value = "transparent",
      breaks = ciw_breaks,
      labels = signif(ciw_breaks, 2),
      limits = c(0, max_abs)
    )
  
  ciw_plot <- add_map_annotations(ciw_plot, study_boundary, bb)
  
  list(
    chg_plot = chg_plot,
    ciw_plot = ciw_plot,
    max_abs = max_abs
  )
}

# ------------------------------------------------------------
# Meaningful change polygons
# ------------------------------------------------------------

#' Flag hexagons with strong posterior support for change
#'
#' @param eta_draws_per_hex Posterior draws per polygon (hex), in the format
#'   expected by `summarize_polygon_hypothesis()`.
#' @param hexagon_sf `sf` object of hex polygons. Must include `hex_id`.
#' @param param Parameter name to test inside `summarize_polygon_hypothesis()`
#'   (e.g., "abs_change").
#' @param threshold Numeric threshold for change hypothesis.
#' @param prob_level Posterior probability threshold (e.g., 0.975).
#' @param direction Directions to compute in `summarize_polygon_hypothesis()`.
#' @param ci_probs Credible interval probs passed through.
#' @param include_summary Passed through to `summarize_polygon_hypothesis()`.
#'
#' @return An `sf` object of flagged hexagons joined to summaries, filtered to
#'   `classification != "None"`. If none flagged, returns 0-row `sf`.
#' @export
flag_hexagons_for_change <- function(eta_draws_per_hex,
                                     hexagon_sf,
                                     param = "abs_change",
                                     threshold = 0,
                                     prob_level = 0.975,
                                     direction = c("two_sided", "increase", "decrease"),
                                     ci_probs = c(0.05, 0.95),
                                     include_summary = TRUE) {
  hexagon_sf <- hexagon_sf %>% dplyr::mutate(hex_id = as.character(hex_id))
  
  flagged <- summarize_polygon_hypothesis(
    eta_draws_per_hex,
    param = param,
    threshold = threshold,
    prob_level = prob_level,
    direction = direction,
    ci_probs = ci_probs,
    include_summary = include_summary
  ) %>%
    dplyr::left_join(hexagon_sf, ., by = c("hex_id" = "poly_id")) %>%
    dplyr::filter(classification != "None") %>%
    sf::st_as_sf()
  
  flagged
}


#' Keep only spatially coherent patches above a minimum area threshold
#'
#' For each change class (e.g., Increase/Decrease), this:
#'   - builds a touches-based adjacency graph of hexagons,
#'   - finds connected components (patches),
#'   - computes patch areas after dissolving,
#'   - retains hexagons belonging to patches with area >= min_area_km2.
#'
#' @param flagged_hex_sf `sf` of flagged hexagons. Must include `classification`
#'   and geometry.
#' @param min_area_km2 Minimum patch area (km^2) to retain.
#'
#' @return `sf` of hexagons with an added `patch_id`, filtered to retained patches.
#'   If input has 0 rows, returns input unchanged.
#' @export
filter_flagged_to_large_patches <- function(flagged_hex_sf, min_area_km2) {
  if (nrow(flagged_hex_sf) == 0) return(flagged_hex_sf)
  
  flagged_patched <- flagged_hex_sf %>%
    dplyr::group_by(classification) %>%
    dplyr::group_modify(~{
      x <- .x
      
      # Neighbor list of touching polygons
      nb <- sf::st_touches(x)
      
      # Graph + connected components
      g <- igraph::graph_from_adj_list(nb, mode = "all")
      g <- igraph::as.undirected(g, mode = "collapse")
      x$patch_id <- igraph::components(g)$membership
      
      # Dissolve by patch and compute patch areas (km^2)
      patch_polys <- x %>%
        dplyr::group_by(patch_id) %>%
        # small "buffer in/out" trick to help remove tiny holes / slivers
        sf::st_buffer(sqrt(min_area_km2)) %>%
        dplyr::summarise(geometry = sf::st_union(geometry), .groups = "drop") %>%
        sf::st_buffer(-sqrt(min_area_km2)) %>%
        dplyr::mutate(
          patch_area_km2 =
            units::set_units(sf::st_area(geometry), "km^2") %>%
            units::drop_units()
        )
      
      keep_ids <- patch_polys$patch_id[patch_polys$patch_area_km2 >= min_area_km2]
      x %>% dplyr::filter(patch_id %in% keep_ids)
    }) %>%
    dplyr::ungroup() %>%
    sf::st_as_sf()
  
  flagged_patched
}


#' Dissolve hexagons into polygons by change classification
#'
#' @param flagged_hex_sf `sf` of flagged hexagons with a `classification` column.
#'
#' @return `sf` with one (multi)polygon per classification.
#' @export
dissolve_change_polygons <- function(flagged_hex_sf) {
  flagged_hex_sf %>%
    dplyr::group_by(classification) %>%
    dplyr::summarise(geometry = sf::st_union(geometry), .groups = "drop") %>%
    sf::st_as_sf()
}


#' Smooth change polygon boundaries and clip to the study boundary
#'
#' @param change_polys_sf `sf` polygons (usually output of `dissolve_change_polygons()`).
#' @param study_boundary `sf` polygon used for CRS target and clipping.
#' @param smoothing_bandwidth_m Numeric bandwidth (meters) for ksmooth.
#'
#' @return Smoothed, valid, boundary-clipped `sf` polygons.
#' @export
smooth_and_clip_change_polys <- function(change_polys_sf,
                                         study_boundary,
                                         smoothing_bandwidth_m) {
  change_polys_sf %>%
    sf::st_transform(sf::st_crs(study_boundary)) %>%
    smoothr::smooth(method = "ksmooth", bandwidth = smoothing_bandwidth_m) %>%
    sf::st_make_valid() %>%
    sf::st_intersection(study_boundary) %>%
    sf::st_as_sf()
}


#' Drop holes smaller than a threshold area (km^2) from polygons/multipolygons
#'
#' @param x An `sf` object or `sfc` geometry with POLYGON/MULTIPOLYGON features.
#' @param min_hole_area_km2 Minimum hole area to keep (km^2). Holes smaller than
#'   this are removed.
#'
#' @return Object of same class as input with holes removed where applicable.
#' @export
drop_small_holes <- function(x, min_hole_area_km2) {
  
  stopifnot(inherits(x, c("sf", "sfc")))
  
  geom <- sf::st_geometry(x)
  crs  <- sf::st_crs(geom)
  
  if (sf::st_is_longlat(geom)) {
    warning("Geometry is lon/lat. Areas will be computed using s2 (geodesic).")
  }
  
  drop_in_poly <- function(poly) {
    
    rings <- unclass(poly)
    outer <- rings[[1]]
    holes <- rings[-1]
    
    if (length(holes) == 0) return(poly)
    
    hole_areas_km2 <- vapply(
      holes,
      function(h) {
        a <- sf::st_area(sf::st_sfc(sf::st_polygon(list(h)), crs = crs))
        as.numeric(units::set_units(a, km^2))
      },
      numeric(1)
    )
    
    keep <- hole_areas_km2 >= min_hole_area_km2
    sf::st_polygon(c(list(outer), holes[keep]))
  }
  
  drop_in_mpoly <- function(mpoly) {
    
    polys <- unclass(mpoly)
    
    new_polys <- lapply(polys, function(rings) {
      drop_in_poly(sf::st_polygon(rings))
    })
    
    sf::st_multipolygon(lapply(new_polys, unclass))
  }
  
  new_geom <- lapply(geom, function(g) {
    if (inherits(g, "POLYGON")) {
      drop_in_poly(g)
    } else if (inherits(g, "MULTIPOLYGON")) {
      drop_in_mpoly(g)
    } else {
      g
    }
  })
  
  if (inherits(x, "sf")) {
    sf::st_geometry(x) <- sf::st_sfc(new_geom, crs = crs)
    return(x)
  } else {
    return(sf::st_sfc(new_geom, crs = crs))
  }
}


#' End-to-end builder for "meaningful change" polygons (flag → patch → dissolve → smooth)
#'
#' This is a convenience wrapper that matches your current inline workflow:
#' - flag hexagons based on posterior support,
#' - retain only patches exceeding `min_area_km2`,
#' - dissolve by Increase/Decrease,
#' - smooth and clip to study boundary,
#' - optionally remove small holes.
#'
#' @param eta_draws_per_hex Posterior draws per hex (see `flag_hexagons_for_change()`).
#' @param hexagon_sf `sf` hex grid with `hex_id`.
#' @param study_boundary `sf` boundary for CRS + clipping.
#' @param min_area_km2 Minimum patch size (km^2).
#' @param smoothing_bandwidth_m Smoothing bandwidth (meters).
#' @param ... Additional args forwarded to `flag_hexagons_for_change()`
#'   (e.g., `param`, `threshold`, `prob_level`, etc.).
#' @param drop_holes Logical; if TRUE, remove holes smaller than `min_area_km2`.
#'
#' @return Smoothed `sf` polygons by classification, or `NULL` if none flagged.
#' @export
build_meaningful_change_polys <- function(eta_draws_per_hex,
                                          hexagon_sf,
                                          study_boundary,
                                          min_area_km2 = 5000,
                                          smoothing_bandwidth_m = 75000,
                                          ...,
                                          drop_holes = TRUE) {
  
  flagged <- flag_hexagons_for_change(
    eta_draws_per_hex = eta_draws_per_hex,
    hexagon_sf = hexagon_sf,
    ...
  )
  
  if (nrow(flagged) == 0) return(NULL)
  
  flagged_patched <- filter_flagged_to_large_patches(
    flagged_hex_sf = flagged,
    min_area_km2 = min_area_km2
  )
  
  if (nrow(flagged_patched) == 0) return(NULL)
  
  change_polys <- dissolve_change_polygons(flagged_patched)
  
  change_polys_smoothed <- smooth_and_clip_change_polys(
    change_polys_sf = change_polys,
    study_boundary = study_boundary,
    smoothing_bandwidth_m = smoothing_bandwidth_m
  )
  
  if (drop_holes) {
    change_polys_smoothed <- drop_small_holes(change_polys_smoothed, min_area_km2)
  }
  
  change_polys_smoothed
}




#' Summarize posterior change within each hexagon and join to hex sf
#'
#' For each hexagon, compute posterior summaries for a parameter such as
#' absolute change, then join those summaries onto the hex geometry.
#'
#' @param eta_draws_per_hex Named list of posterior draws by hexagon.
#''   Expected structure:
#'   eta_draws_per_hex[[hex_id]][[param]] = numeric vector of draws
#' @param hexagon_sf sf object of hexagons. Must include `hex_id`.
#' @param param Name of parameter to summarize, default "abs_change".
#' @param threshold Threshold for evaluating probability of meaningful change.
#' @param ci_probs Length-2 numeric vector giving lower/upper credible interval
#'   probabilities, e.g. c(0.05, 0.95).
#' @param include_median Logical; if TRUE, also compute posterior median.
#'
#' @return sf object with one row per hexagon and posterior summary columns joined.
#' @export
summarize_change_by_hex <- function(eta_draws_per_hex,
                                    hexagon_sf,
                                    param = "abs_change",
                                    threshold = 0,
                                    ci_probs = c(0.05, 0.95),
                                    include_median = TRUE) {
  
  stopifnot(inherits(hexagon_sf, "sf"))
  stopifnot("hex_id" %in% names(hexagon_sf))
  stopifnot(length(ci_probs) == 2)
  
  hex_sf <- hexagon_sf %>%
    dplyr::mutate(hex_id = as.character(hex_id))
  
  # helper to extract a numeric draw vector for one hex
  extract_draws <- function(x, param) {
    if (is.null(x)) return(NULL)
    
    # most likely case: list element contains named components
    if (is.list(x) && !is.null(x[[param]])) {
      return(as.numeric(x[[param]]))
    }
    
    # alternative: data.frame/tibble with a column named param
    if (is.data.frame(x) && param %in% names(x)) {
      return(as.numeric(x[[param]]))
    }
    
    # fallback: if x itself is the numeric draw vector and param is implicit
    if (is.numeric(x)) {
      return(as.numeric(x))
    }
    
    NULL
  }
  
  # if list is unnamed, try to align by order
  if (is.null(names(eta_draws_per_hex))) {
    if (length(eta_draws_per_hex) != nrow(hex_sf)) {
      stop("eta_draws_per_hex is unnamed and does not match nrow(hexagon_sf).")
    }
    names(eta_draws_per_hex) <- hex_sf$hex_id
  }
  
  hex_summaries <- lapply(names(eta_draws_per_hex), function(id) {
    draws <- extract_draws(eta_draws_per_hex[[id]], param = param)
    
    if (is.null(draws) || length(draws) == 0 || all(is.na(draws))) {
      return(data.frame(
        hex_id = as.character(id),
        n_draws = 0,
        mean_change = NA_real_,
        median_change = NA_real_,
        lower = NA_real_,
        upper = NA_real_,
        prob_gt_threshold = NA_real_,
        prob_lt_neg_threshold = NA_real_,
        prob_abs_gt_threshold = NA_real_
      ))
    }
    
    draws <- draws[is.finite(draws)]
    
    if (length(draws) == 0) {
      return(data.frame(
        hex_id = as.character(id),
        n_draws = 0,
        mean_change = NA_real_,
        median_change = NA_real_,
        lower = NA_real_,
        upper = NA_real_,
        prob_gt_threshold = NA_real_,
        prob_lt_neg_threshold = NA_real_,
        prob_abs_gt_threshold = NA_real_
      ))
    }
    
    out <- data.frame(
      hex_id = as.character(id),
      n_draws = length(draws),
      mean_change = mean(draws),
      lower = as.numeric(stats::quantile(draws, ci_probs[1], na.rm = TRUE)),
      upper = as.numeric(stats::quantile(draws, ci_probs[2], na.rm = TRUE)),
      prob_gt_threshold = mean(draws > threshold),
      prob_lt_neg_threshold = mean(draws < -threshold),
      prob_abs_gt_threshold = mean(abs(draws) > threshold)
    )
    
    if (include_median) {
      out$median_change <- stats::median(draws, na.rm = TRUE)
      out <- out[, c(
        "hex_id", "n_draws", "mean_change", "median_change", "lower", "upper",
        "prob_gt_threshold", "prob_lt_neg_threshold", "prob_abs_gt_threshold"
      )]
    } else {
      out$median_change <- NULL
    }
    
    out
  })
  
  hex_summaries <- dplyr::bind_rows(hex_summaries)
  
  hex_sf %>%
    dplyr::left_join(hex_summaries, by = "hex_id") %>%
    sf::st_as_sf()
}


















make_relabund_map <- function(species_name,
                              relabund_rast = NULL,
                              zlim = NULL,
                              layer = "mu_q50",
                              prefix = c("OBBA2", "OBBA3"),
                              study_boundary,
                              bcr_sf = NULL,
                              water_sf = NULL,
                              water_fill = "#EDF7FB",
                              sp_hex_grid = NULL,
                              survey_colour = "#bdbdbd",
                              detection_colour = "black",
                              max_circle_diameter_km = 10,
                              n_surveys_size_quantile = 0.90,
                              title = NULL,
                              subtitle = "Relative abundance") {
  
  prefix <- match.arg(prefix)
  
  stopifnot(inherits(study_boundary, "sf"))
  
  if (is.null(title)) {
    title <- ifelse(prefix == "OBBA2", "Atlas 2", "Atlas 3")
  }
  
  # ------------------------------------------------------------
  # Determine CRS
  # ------------------------------------------------------------
  if (!is.null(relabund_rast)) {
    plot_crs <- sf::st_crs(terra::crs(relabund_rast))
  } else {
    plot_crs <- sf::st_crs(study_boundary)
  }
  
  if (is.na(plot_crs)) {
    stop("study_boundary must have valid CRS if raster is NULL")
  }
  
  study_boundary <- sf::st_transform(study_boundary, plot_crs)
  if (!is.null(bcr_sf))   bcr_sf   <- sf::st_transform(bcr_sf, plot_crs)
  if (!is.null(water_sf)) water_sf <- sf::st_transform(water_sf, plot_crs)
  
  # ------------------------------------------------------------
  # Raster (optional)
  # ------------------------------------------------------------
  r_stars <- NULL
  
  if (!is.null(relabund_rast)) {
    if (is.null(zlim)) stop("zlim required when raster present")
    
    if (length(zlim) == 1) zlim <- c(0, zlim)
    
    r_plot <- relabund_rast[[layer]]
    names(r_plot) <- "relabund"
    r_stars <- stars::st_as_stars(r_plot)
  }
  
  q50_cols <- grDevices::colorRampPalette(c(
    "#FEFEFE", "#FBF7E2", "#FCF8D0", "#EEF7C2", "#CEF2B0",
    "#94E5A0", "#51C987", "#18A065", "#008C59", "#007F53", "#006344"
  ))(10)
  
  # ------------------------------------------------------------
  # Hex data → circle radii
  # ------------------------------------------------------------
  circle_df <- NULL
  
  if (!is.null(sp_hex_grid)) {
    
    required_cols <- c("Atlas", "n_surveys", "n_surveys_det")
    missing_cols <- setdiff(required_cols, names(sp_hex_grid))
    
    if (length(missing_cols) > 0) {
      stop("Missing columns: ", paste(missing_cols, collapse = ", "))
    }
    
    sp_hex_grid <- sp_hex_grid %>%
      sf::st_transform(plot_crs) %>%
      dplyr::filter(
        Atlas == prefix,
        !is.na(n_surveys),
        n_surveys > 0
      ) %>%
      dplyr::mutate(
        n_surveys_det = dplyr::coalesce(n_surveys_det, 0),
        n_surveys_det = pmin(n_surveys_det, n_surveys)
      )
    
    if (nrow(sp_hex_grid) > 0) {
      
      n_max <- stats::quantile(
        sp_hex_grid$n_surveys,
        probs = n_surveys_size_quantile,
        na.rm = TRUE
      )
      
      if (!is.finite(n_max) || n_max <= 0) {
        n_max <- max(sp_hex_grid$n_surveys, na.rm = TRUE)
      }
      
      # convert km → CRS units
      crs_units <- sf::st_crs(plot_crs)$units_gdal
      km_to_units <- ifelse(grepl("metre", crs_units), 1000, 1)
      
      max_radius <- (max_circle_diameter_km / 2) * km_to_units
      
      circle_df <- sp_hex_grid %>%
        dplyr::mutate(
          n_surveys_plot = pmin(n_surveys, n_max),
          n_det_plot     = pmin(n_surveys_det, n_max)
        ) %>%
        sf::st_centroid() %>%
        dplyr::mutate(
          x = sf::st_coordinates(.)[,1],
          y = sf::st_coordinates(.)[,2],
          
          # area ∝ n  → radius ∝ sqrt(n)
          r_gray  = sqrt(n_surveys_plot / n_max) * max_radius,
          r_black = sqrt(n_det_plot     / n_max) * max_radius
        ) %>%
        sf::st_drop_geometry()
    }
  }
  
  bb <- map_bbox(study_boundary)
  
  # ------------------------------------------------------------
  # Build plot
  # ------------------------------------------------------------
  p <- ggplot2::ggplot()
  
  if (!is.null(r_stars)) {
    p <- p +
      stars::geom_stars(data = r_stars) +
      ggplot2::scale_fill_gradientn(
        name = legend_title_map(species_name, title, subtitle, "Posterior median"),
        colors = q50_cols,
        na.value = "transparent",
        limits = zlim
      )
  }
  
  p <- add_base_map_layers(
    p = p,
    study_boundary = study_boundary,
    bcr_sf = bcr_sf,
    water_sf = water_sf,
    atlas_squares_centroids = NULL,
    water_fill = water_fill
  )
  
  # ------------------------------------------------------------
  # Draw circles (key change)
  # ------------------------------------------------------------
  if (!is.null(circle_df) && nrow(circle_df) > 0) {
    
    p <- p +
      ggforce::geom_circle(
        data = circle_df,
        ggplot2::aes(x0 = x, y0 = y, r = r_gray),
        fill = survey_colour,
        color = NA
      ) +
      ggforce::geom_circle(
        data = dplyr::filter(circle_df, r_black > 0),
        ggplot2::aes(x0 = x, y0 = y, r = r_black),
        fill = detection_colour,
        color = NA
      )
  }
  
  p <- add_map_annotations(p, study_boundary, bb)
  
  list(
    sp_plot = p,
    zlim = zlim
  )
}




make_relabund_map <- function(species_name,
                              relabund_rast = NULL,
                              zlim = NULL,
                              layer = "mu_q50",
                              prefix = c("OBBA2", "OBBA3"),
                              study_boundary,
                              bcr_sf = NULL,
                              water_sf = NULL,
                              water_fill = "#EDF7FB",
                              sp_hex_grid = NULL,
                              survey_colour = "#bdbdbd",
                              detection_colour = "black",
                              survey_circle_diameter_km = 20,
                              min_detection_diameter_km = 1,
                              max_detection_diameter_km = 15,
                              mean_count_scaling_max = NULL,
                              survey_scaling_max = 20,
                              alpha_min = 0.2,
                              alpha_max = 1,
                              title = NULL,
                              subtitle = "Relative abundance") {
  
  prefix <- match.arg(prefix)
  stopifnot(inherits(study_boundary, "sf"))
  
  if (survey_scaling_max <= 1) {
    stop("survey_scaling_max must be greater than 1.")
  }
  
  if (min_detection_diameter_km < 0) {
    stop("min_detection_diameter_km must be >= 0.")
  }
  
  if (max_detection_diameter_km <= 0) {
    stop("max_detection_diameter_km must be > 0.")
  }
  
  if (min_detection_diameter_km > max_detection_diameter_km) {
    stop("min_detection_diameter_km must be <= max_detection_diameter_km.")
  }
  
  if (alpha_min < 0 || alpha_min > 1 || alpha_max < 0 || alpha_max > 1) {
    stop("alpha_min and alpha_max must both be between 0 and 1.")
  }
  
  if (alpha_min > alpha_max) {
    stop("alpha_min must be <= alpha_max.")
  }
  
  if (is.null(title)) {
    title <- ifelse(prefix == "OBBA2", "Atlas 2", "Atlas 3")
  }
  
  # ------------------------------------------------------------
  # Determine plotting CRS
  # ------------------------------------------------------------
  if (!is.null(relabund_rast)) {
    stopifnot(inherits(relabund_rast, "SpatRaster"))
    plot_crs <- sf::st_crs(terra::crs(relabund_rast))
  } else {
    plot_crs <- sf::st_crs(study_boundary)
  }
  
  if (is.na(plot_crs)) {
    stop("study_boundary must have a valid CRS if relabund_rast is NULL.")
  }
  
  study_boundary <- sf::st_transform(study_boundary, plot_crs)
  if (!is.null(bcr_sf))   bcr_sf   <- sf::st_transform(bcr_sf, plot_crs)
  if (!is.null(water_sf)) water_sf <- sf::st_transform(water_sf, plot_crs)
  
  # ------------------------------------------------------------
  # Optional raster
  # ------------------------------------------------------------
  r_stars <- NULL
  
  if (!is.null(relabund_rast)) {
    if (is.null(zlim)) {
      stop("zlim is required when relabund_rast is not NULL.")
    }
    
    if (length(zlim) == 1) zlim <- c(0, zlim)
    
    if (!layer %in% names(relabund_rast)) {
      stop("relabund_rast is missing layer: ", layer)
    }
    
    r_plot <- relabund_rast[[layer]]
    names(r_plot) <- "relabund"
    r_stars <- stars::st_as_stars(r_plot)
  }
  
  q50_cols <- grDevices::colorRampPalette(c(
    "#FEFEFE", "#FBF7E2", "#FCF8D0", "#EEF7C2", "#CEF2B0",
    "#94E5A0", "#51C987", "#18A065", "#008C59", "#007F53", "#006344"
  ))(10)
  
  # ------------------------------------------------------------
  # Optional raw-data circles
  # ------------------------------------------------------------
  circle_df <- NULL
  
  if (!is.null(sp_hex_grid)) {
    stopifnot(inherits(sp_hex_grid, "sf"))
    
    required_cols <- c("Atlas", "n_surveys", "mean_count")
    missing_cols <- setdiff(required_cols, names(sp_hex_grid))
    
    if (length(missing_cols) > 0) {
      stop("sp_hex_grid is missing required column(s): ",
           paste(missing_cols, collapse = ", "))
    }
    
    sp_hex_grid <- sp_hex_grid %>%
      sf::st_transform(plot_crs) %>%
      dplyr::filter(
        Atlas == prefix,
        !is.na(n_surveys),
        n_surveys > 0
      ) %>%
      dplyr::mutate(
        mean_count = dplyr::coalesce(mean_count, 0)
      )
    
    if (nrow(sp_hex_grid) > 0) {
      crs_units <- sf::st_crs(plot_crs)$units_gdal
      
      km_to_units <- dplyr::case_when(
        grepl("metre|meter", crs_units, ignore.case = TRUE) ~ 1000,
        grepl("kilometre|kilometer", crs_units, ignore.case = TRUE) ~ 1,
        TRUE ~ 1000
      )
      
      survey_radius <- (survey_circle_diameter_km / 2) * km_to_units
      min_detection_radius <- (min_detection_diameter_km / 2) * km_to_units
      max_detection_radius <- (max_detection_diameter_km / 2) * km_to_units
      
      if (is.null(mean_count_scaling_max)) {
        mean_count_scaling_max <- max(sp_hex_grid$mean_count, na.rm = TRUE)
      }
      
      if (!is.finite(mean_count_scaling_max) || mean_count_scaling_max <= 0) {
        mean_count_scaling_max <- 1
      }
      
      circle_df <- sp_hex_grid %>%
        sf::st_centroid() %>%
        dplyr::mutate(
          x = sf::st_coordinates(.)[, 1],
          y = sf::st_coordinates(.)[, 2],
          
          survey_alpha = alpha_min +
            (pmin(n_surveys, survey_scaling_max) / survey_scaling_max) *
            (alpha_max - alpha_min),
          
          detection_alpha = pmax(survey_alpha, 0.4),
          
          r_survey = survey_radius,
          
          mean_count_scaled = pmin(mean_count, mean_count_scaling_max) /
            mean_count_scaling_max,
          
          r_detection = dplyr::if_else(
            mean_count > 0,
            min_detection_radius +
              sqrt(mean_count_scaled) *
              (max_detection_radius - min_detection_radius),
            0
          )
        ) %>%
        sf::st_drop_geometry()
    }
  }
  
  # ------------------------------------------------------------
  # Build plot
  # ------------------------------------------------------------
  bb <- map_bbox(study_boundary)
  
  p <- ggplot2::ggplot()
  
  if (!is.null(r_stars)) {
    p <- p +
      stars::geom_stars(data = r_stars) +
      ggplot2::scale_fill_gradientn(
        name = legend_title_map(species_name, title, subtitle, "Posterior median"),
        colors = q50_cols,
        na.value = "transparent",
        limits = zlim
      )
  }
  
  p <- add_base_map_layers(
    p = p,
    study_boundary = study_boundary,
    bcr_sf = bcr_sf,
    water_sf = water_sf,
    atlas_squares_centroids = NULL,
    water_fill = water_fill
  )
  
  if (!is.null(circle_df) && nrow(circle_df) > 0) {
    p <- p +
      ggforce::geom_circle(
        data = circle_df,
        ggplot2::aes(
          x0 = x,
          y0 = y,
          r = r_survey,
          alpha = survey_alpha
        ),
        fill = survey_colour,
        color = NA
      ) +
      ggforce::geom_circle(
        data = dplyr::filter(circle_df, mean_count > 0),
        ggplot2::aes(
          x0 = x,
          y0 = y,
          r = r_detection,
          alpha = detection_alpha
        ),
        fill = detection_colour,
        color = NA
      ) +
      ggplot2::scale_alpha_identity()
  }
  
  p <- add_map_annotations(p, study_boundary, bb)
  
  p <- p +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.background  = ggplot2::element_rect(fill = "white", colour = NA),
      panel.grid       = ggplot2::element_blank()
    )
  
  list(
    sp_plot = p,
    zlim = zlim
  )
}