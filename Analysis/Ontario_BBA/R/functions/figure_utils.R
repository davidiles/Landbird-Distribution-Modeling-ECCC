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

#' Wrap a species label for nicer plot titles
#'
#' Long common names can overflow plot margins. This helper inserts line
#' breaks so labels stay readable in map facets and figure titles.
#'
#' @param label Character; species name/label.
#' @param max_length Maximum characters per line before wrapping.
#' @return A wrapped character string.
wrap_species_label <- function(label, max_length = 18) {
  if (nchar(label) <= max_length) return(label)
  words <- strsplit(label, " ")[[1]]
  if (length(words) == 1) return(label)
  paste0(paste(words[-length(words)], collapse = " "), "<br>", words[length(words)])
}

#' Expand a bounding box by a fraction
#'
#' Adds padding around an `sf` bbox so maps have breathing room around the
#' study area.
#'
#' @param bbox An `sf::bbox` object.
#' @param frac Fractional expansion (e.g., 0.10 adds 10% padding).
#' @return Expanded bbox (same class).
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

#' Rasterize an sf grid to a `stars` object for fast plotting
#'
#' Converts a regular prediction grid (sf points/polygons with values) into
#' a `stars` raster at resolution `res`. This can be faster than plotting
#' huge point layers directly with ggplot.
#'
#' @param grid_sf `sf` with geometry and a value column.
#' @param field Name of the value column to rasterize.
#' @param res Resolution (in map units) for the output raster.
#' @return A `stars` object.
rasterize_sf_to_stars <- function(grid_sf, field, res) {
  stopifnot(inherits(grid_sf, "sf"))
  stopifnot(field %in% names(grid_sf))
  
  v <- terra::vect(grid_sf)
  r_template <- terra::rast(v, res = res)
  r <- terra::rasterize(v, r_template, field = field, fun = mean, na.rm = TRUE)
  stars::st_as_stars(r)
}

# Build a consistent "map theme"

#' Standard ggplot theme for maps in this project
#'
#' Removes axes, sets consistent fonts, and standardizes legend placement.
#' Used across figure-building functions so outputs look uniform.
#'
#' @return A `ggplot2` theme object.
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

#' Create relative abundance maps for Atlas 2 and Atlas 3
#'
#' Produces publication-ready maps showing posterior summaries of predicted
#' relative abundance (or intensity) on the prediction grid for each atlas.
#' Designed for the reporting pipeline: handles common overlays (BCR,
#' water, boundary), legends, and optional file output.
#'
#' @param species_name Common name used for titles/filenames.
#' @param grid_sf Prediction grid as `sf`.
#' @param pred_summary Data frame of posterior summaries aligned to `grid_sf`.
#' @param prefix Character vector of atlas labels (default c("OBBA2","OBBA3")).
#' @param study_boundary Boundary polygon for outline.
#' @param bcr_sf Optional BCR polygons for overlay.
#' @param water_sf Optional water polygons/lines for context.
#' @param atlas_squares_centroids Optional centroids for atlas squares.
#' @param template_bbox Optional bbox to enforce identical map extents.
#' @param save_dir Optional directory; if provided, figures are saved there.
#' @return A named list of `ggplot` objects (and/or file paths).
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
    
    # If the caller doesn't supply a detected flag, set to NA.
    if (!("detected" %in% names(atlas_squares_centroids))) atlas_squares_centroids$detected <- NA
    
    base_layers <- c(base_layers, list(
      geom_sf(
        data = atlas_squares_centroids,
        aes(colour = detected),
        shape = 46,
        size = 0.5,
        alpha = 1
      )
    ))
  }
  base_layers <- c(base_layers, list(geom_sf(data = study_boundary, colour = "gray80", fill = NA, linewidth = 0.5, show.legend = FALSE)))
  
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

#' Create absolute change maps between atlases
#'
#' Plots pixel-wise differences in predicted abundance/intensity between
#' Atlas 3 and Atlas 2 (or other configured contrasts), using posterior
#' summaries (median, credible intervals, significance masks).
#'
#' @param species_name Common name used for titles/filenames.
#' @param grid_sf Prediction grid (`sf`).
#' @param change_summary Data frame of change summaries aligned to `grid_sf`.
#' @param study_boundary Boundary polygon.
#' @param bcr_sf Optional BCR overlay.
#' @param water_sf Optional water overlay.
#' @param template_bbox Optional bbox to match extents across figures.
#' @param save_dir Optional output directory to save maps.
#' @return A named list of `ggplot` objects (and/or file paths).
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
      geom_sf(data = subset(signif_poly, classification == "Increase"),
              fill = "transparent", col = "dodgerblue", linewidth = 1),
      geom_sf(data = subset(signif_poly, classification == "Decrease"),
              fill = "transparent", col = "orangered", linewidth = 1)
    ) else NULL) +
    geom_sf(data = study_boundary, colour = "gray80", fill = NA, linewidth = 0.5, show.legend = FALSE) +
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