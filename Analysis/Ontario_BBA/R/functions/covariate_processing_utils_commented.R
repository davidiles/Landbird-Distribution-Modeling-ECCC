#' Download an Earth Engine image as a GeoTIFF
#'
#' Wrapper around `rgee::ee_as_rast()` that exports an `ee$Image` to a local
#' GeoTIFF using the chosen export pathway (default: Google Drive via rgee).
#' This standardizes download arguments used in the covariate pipeline.
#'
#' @param ee_image An Earth Engine image (`ee$Image`).
#' @param region Earth Engine geometry defining the export region.
#' @param out_file Output GeoTIFF path.
#' @param scale_m Pixel scale in metres.
#' @param crs CRS for export (string).
#' @param via Export method passed to rgee (e.g., "drive").
#' @return The output filename (invisibly).
ee_download_tif <- function(ee_image, region, out_file, scale_m, crs, via = "drive") {
  if (!requireNamespace("rgee", quietly = TRUE)) {
    stop("Package 'rgee' is required for Earth Engine downloads.")
  }
  r <- rgee::ee_as_rast(
    image  = ee_image,
    region = region,
    scale  = scale_m,
    crs    = crs,
    via    = via
  )
  terra::writeRaster(r, out_file, overwrite = TRUE)
  invisible(out_file)
}

# Reclassify GHSL into different categories

#' Reclassify GHSL built-up / settlement classes into simplified bins
#'
#' Takes a GHSL categorical raster and collapses original class codes into
#' three broad urbanization categories used in modelling:
#'   1 = mostly uninhabited, 2 = suburban/peri-urban, 3 = dense urban.
#'
#' @param x A `terra::SpatRaster` (or compatible) of GHSL class codes.
#' @return A reclassified raster with values {1,2,3}. 
reclass_ghsl <- function(x) {
  x <- terra::ifel(x %in% c(10, 11), 1, x)                # mostly uninhabited
  x <- terra::ifel(x == 12,          2, x)                # suburban / peri-urban
  x <- terra::ifel(x %in% c(13,21,22,23,30), 3, x)        # dense urban
  x
}

# Calculate mean climate

#' Compute mean WorldClim climate over selected months/years
#'
#' Loads monthly WorldClim GeoTIFFs for a variable, crops/masks to the study
#' boundary in source CRS for speed, then projects each layer to the target
#' CRS and returns the mean across all selected layers.
#'
#' @param variable WorldClim variable name used in filenames (e.g., "tmin").
#' @param years Integer vector of years to include.
#' @param months Integer vector of months (1-12) to include.
#' @param base_dir Base directory containing per-variable subfolders.
#' @param boundary Study boundary `sf` polygon.
#' @param target_crs_wkt Target projection (WKT/PROJ string).
#' @param source_crs CRS of the source rasters (default EPSG:4326).
#' @param method Resampling method for `terra::project()` (e.g., "bilinear").
#' @return A `terra::SpatRaster` of the mean climate surface.
mean_climate <- function(variable,
                         years,
                         months,
                         base_dir,
                         boundary,
                         target_crs_wkt,
                         source_crs = "EPSG:4326",
                         method = "bilinear") {
  
  folder <- file.path(base_dir, variable)
  
  files <- c()
  for (yr in years) {
    for (mo in months) {
      f <- file.path(folder,
                     sprintf("wc2.1_cruts4.09_2.5m_%s_%d-%02d.tif", variable, yr, mo))
      if (!file.exists(f)) stop("Missing WorldClim file: ", f)
      files <- c(files, f)
    }
  }
  
  # boundary in source CRS for fast crop/mask
  boundary_src <- sf::st_transform(boundary, sf::st_crs(source_crs))
  boundary_v   <- terra::vect(boundary_src)
  
  projected_layers <- vector("list", length(files))
  
  for (i in seq_along(files)) {
    r <- terra::rast(files[i])
    
    # Crop/mask in source CRS first
    r <- terra::crop(r, boundary_v)
    r <- terra::mask(r, boundary_v)
    
    # Project to target CRS (metres-based recommended)
    r <- tryCatch(
      terra::project(r, target_crs_wkt, method = method),
      error = function(e) {
        stop("Warp failed for file: ", files[i], "\n", conditionMessage(e))
      }
    )
    
    projected_layers[[i]] <- r
  }
  
  r_stack <- terra::rast(projected_layers)
  terra::mean(r_stack, na.rm = TRUE)
}

# Convert road shapefiles to rasters (with 100 m resolution)

#' Rasterize road presence within a buffer distance
#'
#' Reads a road vector layer, buffers it by `buffer_m`, intersects with the
#' boundary, and rasterizes onto `template` so cells are 1 if they fall
#' within the road buffer (0/NA otherwise).
#'
#' @param roads_path Path to road vector data (any format `sf::st_read()` supports).
#' @param boundary_sf Study boundary polygon (`sf`).
#' @param template `terra::SpatRaster` template defining resolution/extent.
#' @param buffer_m Buffer distance in metres.
#' @return A `terra::SpatRaster` with 0/1 road-buffer presence.
rasterize_roads_buffer_presence <- function(roads_path, boundary_sf, template, buffer_m = 100) {
  
  # Read and transform to raster CRS
  roads <- sf::st_read(roads_path, quiet = TRUE) %>%
    sf::st_make_valid()
  
  roads <- sf::st_transform(roads, sf::st_crs(terra::crs(template)))
  
  # Crop to boundary (avoid buffering the whole province if inputs are huge)
  boundary_tr <- sf::st_transform(boundary_sf, sf::st_crs(terra::crs(template)))
  roads <- suppressWarnings(sf::st_intersection(roads, boundary_tr))
  if (nrow(roads) == 0) {
    # Return all zeros if nothing intersects
    out <- template
    terra::values(out) <- 0
    return(out)
  }
  
  # Buffer and dissolve
  roads_buf <- sf::st_buffer(roads, dist = buffer_m)
  roads_buf <- sf::st_union(roads_buf)
  roads_buf <- sf::st_as_sf(roads_buf)
  
  # Rasterize: 1 inside buffer, 0 outside
  out <- terra::rasterize(terra::vect(roads_buf), template, field = 1, background = 0, touches = TRUE)
  out <- crop_mask_to_boundary(out, boundary_tr)
  out
}

# Rasterize insect layers

#' Rasterize insect outbreak presence by outbreak type
#'
#' Reads an insect disturbance layer, subsets to a set of outbreak types,
#' optionally buffers, and rasterizes each type onto a common template.
#' Output is a multi-layer SpatRaster with one layer per outbreak type.
#'
#' @param insect_path Path to insect polygons/lines.
#' @param boundary_sf Study boundary (`sf`).
#' @param template `terra::SpatRaster` template.
#' @param type_col Column that encodes outbreak type.
#' @param keep_types Character vector of types to keep.
#' @param buffer_m Optional buffer distance (metres) before rasterization.
#' @return A multi-layer `terra::SpatRaster` (one layer per outbreak type).
rasterize_insect_presence_by_type <- function(insect_path,
                                              boundary_sf,
                                              template,
                                              years,
                                              insects_keep,
                                              ranking_exclude = "Light",
                                              simplify_tolerance_m = 100,
                                              touches = TRUE) {
  
  # ---- helpers ----

  #' Create an empty multi-layer raster stack matching a template
  #'
  #' Used when a vector source is missing or has no features after filtering.
  #' Returns a raster (or stack) of all zeros/NA with the same geometry as
  #' `template`, ensuring downstream code can proceed without special-casing.
  #'
  #' @param template `terra::SpatRaster` template.
  #' @return A `terra::SpatRaster` with the same geometry as `template`. 
  empty_stack <- function(template) {
    out <- template
    terra::values(out) <- 0
    names(out) <- "none"
    out
  }
  
  template_crs <- sf::st_crs(terra::crs(template))
  boundary_tr  <- sf::st_transform(boundary_sf, template_crs)
  
  # ---- read + harmonize CRS ----
  insect <- sf::st_read(insect_path, quiet = TRUE) |>
    sf::st_transform(template_crs)
  
  # ---- geometry hygiene (avoid GEOS TopologyException) ----
  # 1) repair + keep polygons only
  insect <- insect |>
    sf::st_make_valid() |>
    suppressWarnings(sf::st_collection_extract("POLYGON"))
  
  # 2) common extra repair for lingering side-location conflicts
  insect <- suppressWarnings(sf::st_buffer(insect, 0))
  
  # 3) drop empties early
  insect <- insect[!sf::st_is_empty(insect), , drop = FALSE]
  if (nrow(insect) == 0) return(empty_stack(template))
  
  # ---- attribute filters ----
  # (assumes fields are named INSECT, EVENT_YEAR, RANKING)
  insect <- insect |>
    dplyr::filter(.data$INSECT %in% insects_keep) |>
    dplyr::filter(.data$EVENT_YEAR %in% years)
  
  if (!is.null(ranking_exclude) && "RANKING" %in% names(insect)) {
    insect <- insect |> dplyr::filter(.data$RANKING != ranking_exclude)
  }
  
  if (nrow(insect) == 0) return(empty_stack(template))
  
  # ---- spatial subset WITHOUT computing intersections ----
  # st_filter + st_intersects only selects features; it does not alter geometry.
  insect <- sf::st_filter(insect, boundary_tr, .predicate = sf::st_intersects)
  if (nrow(insect) == 0) return(empty_stack(template))
  
  # ---- simplify (after filtering reduces complexity) ----
  if (!is.null(simplify_tolerance_m) && simplify_tolerance_m > 0) {
    insect <- suppressWarnings(sf::st_simplify(insect, dTolerance = simplify_tolerance_m))
  }
  
  # ---- union by insect type (wrap for clearer error messages) ----
  insect_union <- tryCatch(
    insect |>
      dplyr::group_by(.data$INSECT) |>
      dplyr::summarise(geometry = sf::st_union(geometry), .groups = "drop") |>
      sf::st_make_valid() |>
      suppressWarnings(sf::st_buffer(., 0)),
    error = function(e) {
      stop(
        "Union failed while dissolving insect polygons by INSECT.\n",
        "Try increasing simplify_tolerance_m (e.g., 250 or 500).\n",
        "Original error: ", conditionMessage(e)
      )
    }
  )
  
  if (nrow(insect_union) == 0) return(empty_stack(template))
  
  # ---- rasterize each insect type to a layer ----
  layers <- vector("list", nrow(insect_union))
  
  for (i in seq_len(nrow(insect_union))) {
    nm   <- insect_union$INSECT[i]
    poly <- insect_union[i, , drop = FALSE]
    
    r_i <- terra::rasterize(
      x = terra::vect(poly),
      y = template,
      field = 1,
      background = 0,
      touches = touches
    )
    names(r_i) <- make.names(nm)
    layers[[i]] <- r_i
  }
  
  out <- terra::rast(layers)
  
  # final clip/mask to boundary (safe; raster operation)
  out <- crop_mask_to_boundary(out, boundary_tr)
  
  out
}



# Rasterize water (rivers and waterbodies)

#' Process OHN water features into raster covariates
#'
#' Reads Ontario Hydro Network (OHN) water layers, optionally buffers/filters
#' features (e.g., lakes, rivers), rasterizes onto `template`, and writes
#' outputs to disk. Designed to be called from the covariate build script.
#'
#' @param water_path Path to OHN water vector data.
#' @param boundary_sf Study boundary (`sf`).
#' @param template Template raster for alignment.
#' @param out_dir Output directory for derived rasters.
#' @param buffer_m Buffer distance in metres (if applicable).
#' @return A named list of output file paths (or rasters), depending on implementation.
process_ohn_water <- function(water_path,
                              boundary_sf,
                              crs_wkt,
                              raster_res_m = 30,
                              river_buffer_m = 100,
                              simplify_tolerance_m = 10,
                              min_area_m2_for_filtered = 1e6,
                              touches = TRUE) {
  
  # ---- prep CRS ----
  crs_target <- sf::st_crs(crs_wkt)
  boundary_tr <- sf::st_transform(boundary_sf, crs_target)
  
  # ---- read + geometry hygiene ----
  water <- sf::st_read(water_path, quiet = TRUE) %>%
    sf::st_transform(crs_target) %>%
    sf::st_make_valid()
  
  # keep polygons only; OHN can have mixed geometry artifacts
  water <- suppressWarnings(sf::st_collection_extract(water, "POLYGON"))
  water <- suppressWarnings(sf::st_buffer(water, 0))   # extra repair step
  water <- water[!sf::st_is_empty(water), , drop = FALSE]
  
  if (nrow(water) == 0) {
    template <- make_template_raster_m(boundary_tr, raster_res_m, crs_wkt)
    terra::values(template) <- 0
    return(list(
      water_open = template,
      water_river = template,
      water_filtered = boundary_tr[0, ]
    ))
  }
  
  # ---- spatial subset WITHOUT intersection ----
  # Keeps only features that touch the boundary, without computing new geometry.
  water <- sf::st_filter(water, boundary_tr, .predicate = sf::st_intersects)
  if (nrow(water) == 0) {
    template <- make_template_raster_m(boundary_tr, raster_res_m, crs_wkt)
    terra::values(template) <- 0
    return(list(
      water_open = template,
      water_river = template,
      water_filtered = boundary_tr[0, ]
    ))
  }
  
  # ---- classify waterbodies ----
  water <- water %>%
    dplyr::mutate(
      water_class = dplyr::case_when(
        .data$WATERBODY_ %in% c("Lake", "Kettle Lake", "Pond", "Reservoir", "Beaver Pond") ~ "Open Water",
        .data$WATERBODY_ %in% c("Canal", "River") ~ "River",
        TRUE ~ NA_character_
      )
    )
  
  water_open  <- water %>% dplyr::filter(.data$water_class == "Open Water")
  water_river <- water %>% dplyr::filter(.data$water_class == "River")
  

  # ---- simplify (optional) ----
  if (!is.null(simplify_tolerance_m) && simplify_tolerance_m > 0) {
    if (nrow(water_open)  > 0) water_open  <- suppressWarnings(sf::st_simplify(water_open,  dTolerance = simplify_tolerance_m))
    if (nrow(water_river) > 0) water_river <- suppressWarnings(sf::st_simplify(water_river, dTolerance = simplify_tolerance_m))
  }
  
  # Buffer river layer (to incorporate some notion of riparian habitat)
  water_river <- water_river %>% st_buffer(river_buffer_m)
  
  # ---- template raster ----
  template <- make_template_raster_m(
    boundary_sf = boundary_tr,
    res_m = raster_res_m,
    crs_wkt = crs_wkt
  )
  
  # Ensure template has values (prevents "SpatRaster has no values" in mask)
  if (!terra::hasValues(template)) terra::values(template) <- 0
  
  # ---- rasterize presence ----
  r_open <- template; terra::values(r_open) <- 0
  r_riv  <- template; terra::values(r_riv)  <- 0
  
  if (nrow(water_open) > 0) {
    r_open <- terra::rasterize(terra::vect(water_open), template, field = 1, background = 0, touches = touches)
  }
  if (nrow(water_river) > 0) {
    r_riv  <- terra::rasterize(terra::vect(water_river), template, field = 1, background = 0, touches = touches)
  }
  
  # ---- large waterbodies layer for masking/plotting ----
  # Avoid intersection; just subset to those intersecting boundary
  water_filtered <- sf::st_read(water_path, quiet = TRUE) %>%
    sf::st_transform(crs_target) %>%
    sf::st_make_valid() %>%
    suppressWarnings(sf::st_collection_extract("POLYGON")) %>%
    suppressWarnings(sf::st_buffer(., 0)) %>%
    dplyr::mutate(area_m2 = as.numeric(sf::st_area(geometry))) %>%
    dplyr::filter(.data$area_m2 >= min_area_m2_for_filtered)
  
  water_filtered <- water_filtered[!sf::st_is_empty(water_filtered), , drop = FALSE]
  if (nrow(water_filtered) > 0) {
    water_filtered <- sf::st_filter(water_filtered, boundary_tr, .predicate = sf::st_intersects)
  }
  
  # Final clip via raster mask (safe)
  list(
    water_open   = crop_mask_to_boundary(r_open, boundary_tr),
    water_river  = crop_mask_to_boundary(r_riv,  boundary_tr),
    water_filtered = water_filtered
  )
}

# Safe raster load: returns NULL if missing

#' Safe raster loader with informative errors
#'
#' Attempts to read a raster from `path` using `terra::rast()`. If loading
#' fails (corrupt file, missing file, driver issues), throws a clearer error
#' message so pipeline failures are easier to diagnose.
#'
#' @param path Path to a raster file.
#' @return A `terra::SpatRaster`.
try_rast <- function(path) {
  if (!file.exists(path)) return(NULL)
  terra::rast(path)
}

# Collapse insect multi-layer raster into two layers:
# broadleaf = max(FTC, Gypsy Moth); needleleaf = max(JPB, Spruce Budworm)

#' Collapse insect layers into broader outbreak indicators
#'
#' Many insect datasets have multiple layers by type or year. This helper
#' collapses them into a smaller set of indicators (e.g., any-insect, or
#' broadleaf vs needleleaf), typically using max/any logic.
#'
#' @param insect_rast Multi-layer `terra::SpatRaster` from insect processing.
#' @return A simplified `terra::SpatRaster` with fewer layers.
collapse_insect_layers <- function(insect_rast) {
  stopifnot(inherits(insect_rast, "SpatRaster"))
  
  nms <- names(insect_rast)
  
  # Helper: find a layer by any of several possible names

  #' Pick the first existing file from a set of candidate paths
  #'
  #' Convenience helper for handling alternative filenames/locations.
  #'
  #' @param candidates Character vector of file paths.
  #' @return The first path that exists; errors if none exist.
  pick <- function(candidates) {
    idx <- match(candidates, nms)
    idx <- idx[!is.na(idx)]
    if (length(idx) == 0) return(NULL)
    insect_rast[[idx[1]]]
  }
  
  # Support raw names + make.names variants
  ftc <- pick(c("Forest Tent Caterpillar", "Forest.Tent.Caterpillar"))
  gm  <- pick(c("Gypsy Moth", "Gypsy.Moth"))
  jpb <- pick(c("Jack Pine Budworm", "Jack.Pine.Budworm"))
  sbw <- pick(c("Spruce Budworm", "Spruce.Budworm"))
  
  # Broadleaf = max(FTC, GM); Needleleaf = max(JPB, SBW)
  broadleaf <- NULL
  needleleaf <- NULL
  
  if (!is.null(ftc) && !is.null(gm)) {
    broadleaf <- max(ftc, gm)  # <-- use base max() (dispatches correctly)
  } else if (!is.null(ftc)) {
    broadleaf <- ftc
  } else if (!is.null(gm)) {
    broadleaf <- gm
  }
  
  if (!is.null(jpb) && !is.null(sbw)) {
    needleleaf <- max(jpb, sbw)  # <-- base max()
  } else if (!is.null(jpb)) {
    needleleaf <- jpb
  } else if (!is.null(sbw)) {
    needleleaf <- sbw
  }
  
  # If one category is missing entirely, create a zero raster so downstream code is stable
  if (is.null(broadleaf)) {
    broadleaf <- insect_rast[[1]]
    terra::values(broadleaf) <- 0
  }
  if (is.null(needleleaf)) {
    needleleaf <- insect_rast[[1]]
    terra::values(needleleaf) <- 0
  }
  
  out <- c(broadleaf, needleleaf)
  names(out) <- c("insect_broadleaf", "insect_needleleaf")
  out
}

# Drop frac_0.* columns that correspond to absence (0 class)

#' Drop covariate columns that are all (or nearly all) zero
#'
#' Habitat proportion covariates can be entirely absent in a study area.
#' This helper removes fraction/proportion columns that are effectively
#' constant zero to avoid singularities and reduce clutter.
#'
#' @param df Data frame containing covariate columns.
#' @return `df` with zero-only fraction columns removed.
drop_frac0_cols <- function(df) {
  df %>% dplyr::select(-dplyr::matches("^frac_0\\.?"))
}

# Clean exactextractr "frac_1.<name>" names into just "<name>"

#' Standardize names of fractional habitat covariates
#'
#' Renames fraction/proportion columns into a consistent naming scheme
#' expected by modelling scripts (e.g., ensuring prefixes like `lc_` and
#' removing awkward characters).
#'
#' @param df Data frame with habitat fraction columns.
#' @return A renamed data frame.
clean_frac_names <- function(df) {
  names(df) <- gsub("^frac_1\\.", "", names(df))
  names(df) <- gsub("^frac1\\.", "", names(df))
  df
}
