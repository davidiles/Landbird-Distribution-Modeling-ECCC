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



# Rasterize rivers
process_river_proximity <- function(x,
                                    boundary_sf,
                                    crs_wkt,
                                    raster_res_m = 100,
                                    buffer_m = 250,
                                    simplify_tolerance_m = 10,
                                    touches = TRUE,
                                    use_raster_distance_buffer = TRUE) {
  
  crs_target  <- sf::st_crs(crs_wkt)
  boundary_tr <- sf::st_transform(boundary_sf, crs_target)
  
  template <- make_template_raster_m(
    boundary_sf = boundary_tr,
    res_m = raster_res_m,
    crs_wkt = crs_wkt
  )
  terra::values(template) <- 0
  
  empty_raster <- function() {
    r <- template
    terra::values(r) <- 0
    r
  }
  
  if (nrow(x) == 0) return(empty_raster())
  
  old_s2 <- sf::sf_use_s2()
  on.exit(sf::sf_use_s2(old_s2), add = TRUE)
  sf::sf_use_s2(FALSE)
  
  # Transform + clean
  x <- sf::st_transform(x, crs_target)
  x <- x[!sf::st_is_empty(x), , drop = FALSE]
  
  # Simplify if needed
  if (!is.null(simplify_tolerance_m) && simplify_tolerance_m > 0) {
    x <- suppressWarnings(sf::st_simplify(
      x,
      dTolerance = simplify_tolerance_m,
      preserveTopology = TRUE
    ))
    x <- x[!sf::st_is_empty(x), , drop = FALSE]
  }
  
  if (nrow(x) == 0) return(empty_raster())
  
  # Rasterize / buffer
  if (use_raster_distance_buffer && buffer_m > 0) {
    
    r0 <- terra::rasterize(
      x = terra::vect(x),
      y = template,
      field = 1,
      background = NA,
      touches = touches,
      fun = "max"
    )
    
    d <- terra::distance(r0)
    r <- terra::ifel(d <= buffer_m, 1, 0)
    
  } else {
    
    if (!is.null(buffer_m) && buffer_m > 0) {
      x <- suppressWarnings(sf::st_buffer(x, dist = buffer_m))
      x <- x[!sf::st_is_empty(x), , drop = FALSE]
      if (nrow(x) == 0) return(empty_raster())
    }
    
    r <- terra::rasterize(
      x = terra::vect(x),
      y = template,
      field = 1,
      background = 0,
      touches = touches,
      fun = "max"
    )
  }
  
  r <- crop_mask_to_boundary(r, boundary_tr)
  names(r) <- "river_proximity"
  
  r
}

# Rasterize lakes
process_lake_proximity <- function(x,
                                   boundary_sf,
                                   crs_wkt,
                                   raster_res_m = 100,
                                   buffer_m = 250,
                                   simplify_tolerance_m = 10,
                                   touches = TRUE,
                                   use_raster_distance_buffer = TRUE,
                                   out_name = "lake_proximity") {
  
  crs_target  <- sf::st_crs(crs_wkt)
  boundary_tr <- sf::st_transform(boundary_sf, crs_target)
  
  template <- make_template_raster_m(
    boundary_sf = boundary_tr,
    res_m = raster_res_m,
    crs_wkt = crs_wkt
  )
  terra::values(template) <- 0
  
  empty_raster <- function() {
    r <- template
    terra::values(r) <- 0
    names(r) <- out_name
    r
  }
  
  if (is.null(x) || nrow(x) == 0) {
    return(empty_raster())
  }
  
  old_s2 <- sf::sf_use_s2()
  on.exit(sf::sf_use_s2(old_s2), add = TRUE)
  sf::sf_use_s2(FALSE)
  
  x <- sf::st_transform(x, crs_target)
  x <- suppressWarnings(sf::st_make_valid(x))
  x <- suppressWarnings(sf::st_collection_extract(x, "POLYGON"))
  x <- x[!sf::st_is_empty(x), , drop = FALSE]
  
  if (nrow(x) == 0) {
    return(empty_raster())
  }
  
  x <- suppressWarnings(sf::st_crop(x, sf::st_bbox(boundary_tr)))
  x <- x[!sf::st_is_empty(x), , drop = FALSE]
  
  if (nrow(x) == 0) {
    return(empty_raster())
  }
  
  if (!is.null(simplify_tolerance_m) && simplify_tolerance_m > 0) {
    x <- suppressWarnings(sf::st_simplify(
      x,
      dTolerance = simplify_tolerance_m,
      preserveTopology = TRUE
    ))
    x <- x[!sf::st_is_empty(x), , drop = FALSE]
  }
  
  if (nrow(x) == 0) {
    return(empty_raster())
  }
  
  if (use_raster_distance_buffer && buffer_m > 0) {
    
    r0 <- terra::rasterize(
      x = terra::vect(x),
      y = template,
      field = 1,
      background = NA,
      touches = touches,
      fun = "max"
    )
    
    d <- terra::distance(r0)
    r <- terra::ifel(d <= buffer_m, 1, 0)
    
  } else {
    
    if (!is.null(buffer_m) && buffer_m > 0) {
      x <- suppressWarnings(sf::st_buffer(x, dist = buffer_m))
      x <- x[!sf::st_is_empty(x), , drop = FALSE]
      
      if (nrow(x) == 0) {
        return(empty_raster())
      }
    }
    
    r <- terra::rasterize(
      x = terra::vect(x),
      y = template,
      field = 1,
      background = 0,
      touches = touches,
      fun = "max"
    )
  }
  
  r <- crop_mask_to_boundary(r, boundary_tr)
  names(r) <- out_name
  
  r
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

#' Assert that required input files exist
#'
#' Convenience validator to ensure that all required input files are present
#' before running downstream processing steps. If any files are missing, the
#' function stops with a readable error listing all missing paths.
#'
#' @param x Character vector of file paths to check.
#' @param label Character scalar used in the error message to describe what the
#'   files are for (e.g., `"landcover rasters"`).
#'
#' @return Invisibly returns `NULL` if all files exist.
#' @export
need_files <- function(x, label) {
  missing <- x[!file.exists(x)]
  if (length(missing) > 0) {
    stop(
      "Missing required input(s) for ", label, ":\n  - ",
      paste(missing, collapse = "\n  - ")
    )
  }
}

#' Construct a path inside the covariate directory
#'
#' Helper for building paths relative to a covariate directory (`cov_dir`).
#' This is convenient for consistently referencing covariate inputs.
#'
#' @param x Character scalar (or vector) giving a filename or relative path
#'   within `cov_dir`.
#'
#' @details
#' This function assumes that `cov_dir` exists in the calling environment.
#' It does not validate that the resulting path exists.
#'
#' @return A character vector of full file paths.
#' @export
r_path <- function(x) file.path(cov_dir, x)

#' Crop and mask a raster to a boundary (NULL-safe)
#'
#' Crops and masks a `terra` SpatRaster to a study boundary polygon for speed
#' and memory efficiency. Designed to be pipe-friendly: if `r` is `NULL`, the
#' function returns `NULL` without error.
#'
#' @param r A `terra::SpatRaster`, or `NULL`.
#' @param boundary_sf An `sf` polygon object defining the boundary used for
#'   cropping/masking.
#'
#' @return A cropped and masked `terra::SpatRaster`, or `NULL` if `r` is `NULL`.
#' @export
crop_if <- function(r, boundary_sf) {
  if (is.null(r)) return(NULL)
  v <- terra::vect(boundary_sf)
  r <- terra::crop(r, v)
  terra::mask(r, v)
}

#' Transform grid polygons to the raster CRS
#'
#' Ensures that polygon geometries used for extraction match the CRS of a
#' `terra` raster. This is particularly important for `exactextractr`, which
#' expects CRS agreement.
#'
#' @param grid_cells An `sf` object containing polygon geometries (grid cells).
#' @param r A `terra::SpatRaster` whose CRS should be matched.
#'
#' @return `grid_cells` transformed to the CRS of `r` (an `sf` object).
#' @export
match_grid_to_raster_crs <- function(grid_cells, r) {
  st_transform(grid_cells, st_crs(terra::crs(r)))
}

#' Assign polygon IDs to point features
#'
#' Assigns a polygon identifier to each point. Points are first matched using
#' `sf::st_within()`. Optionally, points that do not fall within any polygon can
#' be assigned to the nearest polygon using `sf::st_nearest_feature()`.
#'
#' @param pts An `sf` object with point geometries.
#' @param polys An `sf` object with polygon geometries.
#' @param id_col Character scalar naming the ID column in `polys` to return.
#' @param nearest_fallback Logical; if `TRUE`, points with no containing polygon
#'   will be assigned the nearest polygon ID.
#'
#' @details
#' `polys` is transformed to the CRS of `pts` before spatial operations.
#' If a point falls within multiple polygons, the first match is used.
#'
#' @return A vector of polygon IDs (from `polys[[id_col]]`) aligned to rows of
#'   `pts`. If `nearest_fallback = FALSE`, unmatched points will be `NA`.
#' @export
assign_poly_id <- function(pts, polys, id_col, nearest_fallback = TRUE) {
  stopifnot(inherits(pts, "sf"), inherits(polys, "sf"))
  stopifnot(id_col %in% names(polys))
  
  polys <- st_transform(polys, st_crs(pts))
  
  hit <- st_within(pts, polys)               # list column of indices
  idx <- rep(NA_integer_, nrow(pts))
  
  has_hit <- lengths(hit) > 0L
  if (any(has_hit)) {
    idx[has_hit] <- vapply(hit[has_hit], function(x) x[1], integer(1))
  }
  
  if (nearest_fallback) {
    missing <- is.na(idx)
    if (any(missing)) {
      idx[missing] <- st_nearest_feature(pts[missing, , drop = FALSE], polys)
    }
  }
  
  polys[[id_col]][idx]
}

#' Extract polygon-level mean from a raster
#'
#' Computes the mean raster value within each grid polygon using exact,
#' area-weighted extraction via `exactextractr::exact_extract()`. The raster is
#' cropped/masked to the study boundary first for efficiency.
#'
#' @param r A `terra::SpatRaster`, or `NULL`.
#' @param grid_cells An `sf` object of polygons over which to extract values.
#' @param boundary_sf An `sf` polygon defining the study boundary used to crop
#'   and mask the raster.
#'
#' @return A numeric vector of length `nrow(grid_cells)` giving the mean value
#'   for each polygon. If `r` is `NULL`, returns an `NA_real_` vector of that
#'   length.
#' @export
extract_mean <- function(r, grid_cells, boundary_sf) {
  if (is.null(r)) return(rep(NA_real_, nrow(grid_cells)))
  r <- crop_if(r, boundary_sf)
  exact_extract(r, grid_cells, "mean")
}

#' Extract fractional coverage for categorical rasters with consistent naming
#'
#' Extracts fractional coverage (proportion of polygon area) for each raster
#' category within each grid polygon using `exactextractr::exact_extract(..., "frac")`.
#' Includes optional column cleanup steps: dropping all-zero categories, adding
#' a prefix to category columns, and normalizing names returned by exactextractr
#' for named raster layers.
#'
#' @param r A categorical `terra::SpatRaster`, or `NULL`.
#' @param grid_cells An `sf` object of polygons over which to extract fractions.
#' @param boundary_sf An `sf` polygon defining the study boundary used to crop
#'   and mask the raster.
#' @param prefix Optional character scalar. If provided, renames columns like
#'   `frac_1`/`frac1` to `<prefix>_1`, `<prefix>_2`, etc.
#' @param drop0 Logical; if `TRUE`, drop fraction columns that are entirely zero.
#' @param clean_names Logical; if `TRUE`, clean column names in cases where
#'   exactextractr returns layer-qualified names (e.g., `frac_1.<layername>`).
#'
#' @details
#' The grid polygons are transformed to the raster CRS before extraction (via
#' `match_grid_to_raster_crs()`).
#'
#' This function expects helper utilities `drop_frac0_cols()` and
#' `clean_frac_names()` to be defined elsewhere in your codebase.
#'
#' @return A data frame of fractional coverage columns, or `NULL` if `r` is
#'   `NULL` or extraction returns no rows.
#' @export
extract_frac_clean <- function(r, grid_cells, boundary_sf, prefix = NULL, drop0 = TRUE, clean_names = FALSE) {
  if (is.null(r)) return(NULL)
  r <- crop_if(r, boundary_sf)
  gc <- match_grid_to_raster_crs(grid_cells, r)
  
  df <- suppressWarnings(exact_extract(r, gc, "frac"))
  if (is.null(df) || nrow(df) == 0) return(NULL)
  
  if (drop0) df <- drop_frac0_cols(df)
  
  if (!is.null(prefix)) {
    # frac_1, frac_2 -> <prefix>_1, <prefix>_2 (works for both frac_# and frac# variants)
    names(df) <- gsub("^frac_?", paste0(prefix, "_"), names(df))
  }
  
  if (clean_names) {
    # For named layers (e.g., insect_broadleaf), exactextractr may return frac_1.<layername>
    df <- clean_frac_names(df)
  }
  
  df
}
