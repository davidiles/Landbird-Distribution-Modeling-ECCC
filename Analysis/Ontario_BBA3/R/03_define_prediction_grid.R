# ============================================================
# 03_define_prediction_grid.R
#
# Purpose:
#   Create a regular prediction grid over the Ontario study area
#   for mapping and model prediction. Produces both grid-cell polygons
#   (for landscape summaries) and grid centroids (for prediction points).
#
# Key choices:
#   - Grid is constructed in EPSG:3978 (meters) to align with processed rasters.
#   - Centroids are used as prediction locations; polygons are retained for
#     area-based covariate extraction in the next step.
#
# Output:
#   data_clean/spatial/prediction_grid.rds
# ============================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(here)
})

# Centralized paths
source(here::here("R", "00_config_paths.R"))

# Ensure output directory exists
out_dir <- file.path(paths$data_clean, "spatial")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Config
# ------------------------------------------------------------

grid_cellsize_m <- 1000       # metres
grid_crs <- sf::st_crs(3978)  # EPSG:3978 (Canada Albers, metres)

# ------------------------------------------------------------
# Inputs
# ------------------------------------------------------------

study_area_path <- file.path(paths$data_clean, "spatial", "study_area.rds")
if (!file.exists(study_area_path)) {
  stop("Cannot find study area at: ", study_area_path,
       "\nHave you run 01_define_study_area.R?")
}
study_area <- readRDS(study_area_path)

# Use the buffered boundary for prediction domain (will be trimmed later)
boundary <- study_area$boundary_buffer_5km %>%
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
  rename(geometry = x)

# Compute centroids (safe: grid is projected + cells are regular squares)
grid_pts <- st_centroid(grid_polys)

# Keep cells whose centroid intersects the boundary
inside <- st_intersects(grid_pts, boundary, sparse = FALSE)[, 1]

grid_polys <- grid_polys[inside, , drop = FALSE]
grid_pts   <- grid_pts[inside, , drop = FALSE]

# Assign pixel ids AFTER filtering (stable within a given run)
grid_polys <- grid_polys %>% mutate(pixel_id = row_number()) %>% relocate(pixel_id)
grid_pts   <- grid_pts   %>% mutate(pixel_id = row_number()) %>% select(pixel_id, geometry)

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

out_path <- file.path(out_dir, paste0("prediction_grid_",grid_cellsize_m,"_m.rds"))

saveRDS(
  list(
    grid_cells     = grid_polys,
    grid_centroids = grid_pts,
    cellsize_m     = grid_cellsize_m,
    crs            = grid_crs,
    bbox           = st_bbox(boundary),
    n_cells        = nrow(grid_polys),
    date_created   = Sys.time()
  ),
  file = out_path
)

message("03_define_prediction_grid.R complete.")
message("Saved to: ", out_path)
