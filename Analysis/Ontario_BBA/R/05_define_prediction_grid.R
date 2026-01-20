# ============================================================
# 05_define_prediction_grid.R
#
# Purpose:
#   Create a regular prediction grid over the Ontario study area
#   for mapping and model prediction. Produces both grid-cell polygons
#   (for landscape summaries) and grid centroids (for prediction points).
#
# Key choices:
#   - Grid is constructed in EPSG:3978 (meters) to align with processed rasters.
#   - Centroids are used as the prediction locations; polygons are retained for
#     buffer-based covariate extraction (exactextractr) in the next step.
#
# Output:
#   data_clean/spatial/prediction_grid.rds
# ============================================================

library(sf)
library(dplyr)

dir.create("data_clean/spatial", recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Config
# ------------------------------------------------------------

grid_cellsize_m <- 1000       # in metres
grid_crs <- sf::st_crs(3978)  # EPSG:3978 (meters)

# ------------------------------------------------------------
# Inputs
# ------------------------------------------------------------

study_area <- readRDS("data_clean/spatial/study_area.rds")

# Use the unbuffered boundary for defining prediction domain
boundary <- study_area$boundary %>%
  st_transform(grid_crs) %>%
  st_make_valid()

# ------------------------------------------------------------
# Build grid polygons
# ------------------------------------------------------------

grid_polys <- st_make_grid(
  boundary,
  cellsize = grid_cellsize_m,
  what = "polygons",
  square = TRUE
) %>%
  st_as_sf() %>%
  rename(geometry = x) %>%
  mutate(pixel_id = row_number())

# Compute centroids
grid_pts <- st_centroid(grid_polys) %>%
  dplyr::select(pixel_id)

# Keep cells whose centroid is inside the boundary
inside <- st_within(grid_pts, boundary, sparse = FALSE)[, 1]

grid_polys <- grid_polys[inside, ]
grid_pts   <- grid_pts[inside, ]

# Redefine pixel ids
grid_polys$pixel_id <- 1:nrow(grid_polys)
grid_pts$pixel_id <- 1:nrow(grid_pts)

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

saveRDS(
  list(
    grid_cells = grid_polys,
    grid_centroids = grid_pts,
    cellsize_m = grid_cellsize_m,
    crs = grid_crs,
    date_created = Sys.time()
  ),
  file = "data_clean/spatial/prediction_grid.rds"
)
