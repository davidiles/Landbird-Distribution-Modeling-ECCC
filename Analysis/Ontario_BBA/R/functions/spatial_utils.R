
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