# ============================================================
# covariate_processing_utils.R
#
# Utility functions for rasterizing vector covariates, matching sf/raster
# coordinate reference systems, and extracting raster summaries to grid cells.
# CRS transformations are intentionally handled inside helper functions; do not
# simplify or replace CRS definitions without checking downstream assumptions.
# ============================================================

# Simplify, buffer, rasterize, and optionally save an sf layer as a binary raster.
buf_rasterize_save <- function(sf_obj, buffer_dist_km, out_path, label,
                               template     = template_30m,
                               simplify_tol = SIMPLIFY_KM) {
  
  message(sprintf("[%s] %d features | simplify %.0f m | buffer %.0f m | rasterizing...",
                  label, nrow(sf_obj), simplify_tol * 1000, buffer_dist_km * 1000))
  
  v <- sf_obj %>%
    sf::st_simplify(preserveTopology = TRUE, dTolerance = simplify_tol) %>%
    sf::st_make_valid() %>%
    sf::st_buffer(dist = buffer_dist_km) %>%
    sf::st_make_valid() %>%
    terra::vect()
  
  r <- terra::rasterize(v, template, field = 1L,
                        background = 0L, touches = TRUE)
  
  if (!is.null(out_path)) {
    terra::writeRaster(r, out_path, datatype = "INT1U",
                       overwrite = TRUE, gdal = GTIFF_OPTS)
    message(sprintf("[%s] Saved -> %s", label, basename(out_path)))
  }
  
  invisible(r)
}

# Load a raster when the file exists; return NULL otherwise.
try_rast <- function(path) {
  if (!file.exists(path)) return(NULL)
  terra::rast(path)
}

# Create square polygons centred on point geometries.
make_square_buffer <- function(points, half_width = 500) {
  coords <- st_coordinates(points)
  
  squares <- lapply(seq_len(nrow(coords)), function(i) {
    x <- coords[i, "X"]
    y <- coords[i, "Y"]
    st_polygon(list(matrix(c(
      x - half_width, y - half_width,
      x + half_width, y - half_width,
      x + half_width, y + half_width,
      x - half_width, y + half_width,
      x - half_width, y - half_width   # close the polygon
    ), ncol = 2, byrow = TRUE)))
  })
  
  st_sfc(squares, crs = st_crs(points))
}

# Assign each point the ID of its containing polygon, optionally using nearest fallback.
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

# Extract mean raster values within grid cells, optionally treating raster NA as zero.
extract_mean <- function(r, grid_cells, boundary_sf, na_as_zero = FALSE) {
  if (is.null(r)) return(rep(NA_real_, nrow(grid_cells)))
  
  boundary_sf <- match_sf_to_raster_crs(boundary_sf, r)
  grid_cells  <- match_sf_to_raster_crs(grid_cells, r)
  
  r <- crop_if(r, boundary_sf)
  
  if (na_as_zero) {
    r[is.na(r)] <- 0
  }
  
  exactextractr::exact_extract(r, grid_cells, "mean")
}

# Extract fractional coverage of categorical raster classes within grid cells.
extract_frac <- function(r, grid_cells, boundary_sf,
                         prefix = NULL,
                         drop0 = TRUE,
                         clean_names = FALSE) {
  if (is.null(r)) return(NULL)
  
  boundary_sf <- match_sf_to_raster_crs(boundary_sf, r)
  grid_cells  <- match_sf_to_raster_crs(grid_cells, r)
  
  r <- crop_if(r, boundary_sf)
  
  df <- suppressWarnings(exactextractr::exact_extract(r, grid_cells, "frac"))
  if (is.null(df) || nrow(df) == 0) return(NULL)
  
  if (drop0) df <- drop_frac0_cols(df)
  
  if (!is.null(prefix)) {
    names(df) <- gsub("^frac_?", paste0(prefix, "_"), names(df))
  }
  
  if (clean_names) {
    df <- clean_frac_names(df)
  }
  
  df
}

# Transform sf data to match a raster CRS after checking that both CRS values exist.
match_sf_to_raster_crs <- function(x, r) {
  r_crs <- sf::st_crs(terra::crs(r, proj = TRUE))
  
  if (is.na(r_crs)) {
    stop("Raster has missing CRS.")
  }
  if (is.na(sf::st_crs(x))) {
    stop("sf object has missing CRS.")
  }
  
  if (sf::st_crs(x) != r_crs) {
    x <- sf::st_transform(x, r_crs)
  }
  
  x
}

# Clean class-fraction column names returned by exactextractr.
clean_frac_names <- function(df) {
  names(df) <- gsub("^frac_1\\.", "", names(df))
  names(df) <- gsub("^frac1\\.", "", names(df))
  df
}

# Drop fraction columns for the zero/background class.
drop_frac0_cols <- function(df) {
  df %>% dplyr::select(-dplyr::matches("^frac_0\\.?"))
}
