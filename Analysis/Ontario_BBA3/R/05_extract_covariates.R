# ============================================================
# 05_extract_covariates.R
#
# Purpose:
#   1) Compute covariates on the prediction grid for OBBA2 and OBBA3
#      from processed rasters (EPSG:3978), summarizing native-resolution
#      rasters onto the harmonized grid (step 05).
#   2) Attach covariates to survey points by nearest grid centroid
#      (OBBA2 surveys -> OBBA2 grid covs; OBBA3 surveys -> OBBA3 grid covs).
#
# Inputs:
#   - data_clean/spatial/study_area.rds
#   - data_clean/surveys/surveys_raw.rds
#   - data_clean/surveys/count_matrix_raw.rds
#   - data_clean/metadata/species_list.rds
#   - data_clean/spatial/*.tif (from 04_*)
#   - data_clean/spatial/prediction_grid.rds
#
# Outputs:
#   - data_clean/birds/analysis_data_covariates.rds
# ============================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(terra)
  library(exactextractr)
  library(here)
})

# Centralized paths
source(here::here("R", "00_config_paths.R"))

# ------------------------------------------------------------
# Load utilities
# ------------------------------------------------------------

spatial_utils_path <- file.path(paths$functions, "spatial_utils.R")
covariate_processing_utils_path <- file.path(paths$functions, "covariate_processing_utils.R")

source(spatial_utils_path)
source(covariate_processing_utils_path)

# Ensure output directory exists
birds_dir <- file.path(paths$data_clean, "birds")
dir.create(birds_dir, recursive = TRUE, showWarnings = FALSE)

spatial_dir <- file.path(paths$data_clean, "spatial")

# ------------------------------------------------------------
# Config
# ------------------------------------------------------------

crs_rast <- sf::st_crs(3978)   # EPSG:3978

cov_dir  <- file.path(spatial_dir)
out_file <- file.path(paths$data_clean, "birds", "analysis_data_covariates.rds")

in_grid         <- file.path(spatial_dir,  "prediction_grid_1000_m.rds")

in_surveys      <- file.path(paths$data_clean, "surveys",  "surveys_raw.rds")
in_count_matrix <- file.path(paths$data_clean, "surveys",  "count_matrix_raw.rds")
in_species      <- file.path(paths$data_clean, "metadata", "species_list.rds")

in_study_area   <- file.path(spatial_dir,  "study_area.rds")

# Optional national layers (shared Data/)
in_atlas_squares <- file.path(paths$data, "Spatial", "National", "AtlasSquares", "NationalSquares_FINAL.shp")

# ------------------------------------------------------------
# Load inputs
# ------------------------------------------------------------

all_surveys  <- readRDS(in_surveys)
count_matrix <- readRDS(in_count_matrix)
all_species  <- readRDS(in_species)

grid_obj <- readRDS(in_grid)
grid_cells     <- grid_obj$grid_cells
grid_centroids <- grid_obj$grid_centroids

study_area <- readRDS(in_study_area)
boundary_buf <- study_area$boundary_buffer_5km %>%
  st_transform(crs_rast) %>%
  st_make_valid()

# Ensure grid is in extraction CRS
grid_cells     <- st_transform(grid_cells, crs_rast)
grid_centroids <- st_transform(grid_centroids, crs_rast)

# Stable ordering field
if (!("obs_idx" %in% names(all_surveys))) {
  all_surveys$obs_idx <- seq_len(nrow(all_surveys))
}
all_surveys <- all_surveys %>% relocate(obs_idx)

# Surveys in extraction CRS
all_surveys_m <- st_transform(all_surveys, crs_rast)

# Create 1-km square polygons centered on each survey
survey_pixels <- make_square_buffer(all_surveys_m, half_width = 500)
survey_pixels <- st_set_geometry(all_surveys_m, survey_pixels)

# ------------------------------------------------------------
# Load covariate files
# ------------------------------------------------------------

# ---- Vector files
roads_2005  <- try_rast(file.path(spatial_dir,"roads_2005_buf_125m.tif"))
roads_2025  <- try_rast(file.path(spatial_dir,"roads_2025_buf_125m.tif"))

rivers_large <- try_rast(file.path(spatial_dir,"rivers_large_buf_250m.tif"))
rivers_small <- try_rast(file.path(spatial_dir,"rivers_small_buf_125m.tif"))

lakes_large <- try_rast(file.path(spatial_dir,"lakes_large_buf_250m.tif"))
lakes_small <- try_rast(file.path(spatial_dir,"lakes_small_buf_250m.tif"))

coastline   <- try_rast(file.path(spatial_dir,"coastline_buf_500m.tif"))
great_lakes <- try_rast(file.path(spatial_dir,"great_lakes_buf_500m.tif"))
open_water  <- try_rast(file.path(spatial_dir,"waterbodies.tif"))

# Land cover
# NOTE: THIS SHOULD GET SUPPLEMENTED WITH ADDITIONAL RASTERS (e.g., SCANFI, CANLABS, Far North Land Cover, etc)
LCC2020 <- try_rast(file.path(paths$data,"Spatial","National","LandCoverCanada2020","landcover-2020-classification.tif"))

# ------------------------------------------------------------
# Build grid covariate tables (OBBA2 / OBBA3)
# ------------------------------------------------------------

# Surveys in each atlas (for model fitting)
surveys2 <- survey_pixels %>% filter(Atlas == "OBBA2")
surveys3 <- survey_pixels %>% filter(Atlas == "OBBA3")

# Grid cells for each atlas (for prediction)
grid_OBBA2 <- grid_centroids
grid_OBBA3 <- grid_centroids

# LCC2020 fractions -> LCC2020_#
LCC2020_frac <- extract_frac(LCC2020, grid_cells, boundary_buf, prefix = "LCC2020", drop0 = TRUE, clean_names = FALSE)
grid_OBBA2 <- bind_cols(grid_OBBA2, LCC2020_frac)
grid_OBBA3 <- bind_cols(grid_OBBA3, LCC2020_frac)

LCC2020_frac_surveys2 <- extract_frac(LCC2020, surveys2, boundary_buf, prefix = "LCC2020", drop0 = TRUE, clean_names = FALSE)
surveys2   <- bind_cols(surveys2, LCC2020_frac_surveys2)
LCC2020_frac_surveys3 <- extract_frac(LCC2020, surveys3, boundary_buf, prefix = "LCC2020", drop0 = TRUE, clean_names = FALSE)
surveys3   <- bind_cols(surveys3, LCC2020_frac_surveys3)

# Roads
grid_OBBA2$road <- extract_mean(roads_2005, grid_cells, boundary_buf)
surveys2$road <- extract_mean(roads_2005, surveys2, boundary_buf)
grid_OBBA3$road <- extract_mean(roads_2025, grid_cells, boundary_buf)
surveys3$road <- extract_mean(roads_2025, surveys3, boundary_buf)

# Large rivers
grid_OBBA2$rivers_large <- extract_mean(rivers_large, grid_cells, boundary_buf)
grid_OBBA3$rivers_large <- grid_OBBA2$rivers_large
surveys2$rivers_large <- extract_mean(rivers_large, surveys2, boundary_buf)
surveys3$rivers_large <- extract_mean(rivers_large, surveys3, boundary_buf)

# Small rivers
grid_OBBA2$rivers_small <- extract_mean(rivers_small, grid_cells, boundary_buf)
grid_OBBA3$rivers_small <- grid_OBBA2$rivers_small
surveys2$rivers_small <- extract_mean(rivers_small, surveys2, boundary_buf)
surveys3$rivers_small <- extract_mean(rivers_small, surveys3, boundary_buf)

# Large lakes
grid_OBBA2$lakes_large <- extract_mean(lakes_large, grid_cells, boundary_buf)
grid_OBBA3$lakes_large <- grid_OBBA2$lakes_large
surveys2$lakes_large <- extract_mean(lakes_large, surveys2, boundary_buf)
surveys3$lakes_large <- extract_mean(lakes_large, surveys3, boundary_buf)

# Small lakes
grid_OBBA2$lakes_small <- extract_mean(lakes_small, grid_cells, boundary_buf)
grid_OBBA3$lakes_small <- grid_OBBA2$lakes_small
surveys2$lakes_small <- extract_mean(lakes_small, surveys2, boundary_buf)
surveys3$lakes_small <- extract_mean(lakes_small, surveys3, boundary_buf)

# Great lakes
grid_OBBA2$great_lakes <- extract_mean(great_lakes, grid_cells, boundary_buf)
grid_OBBA3$great_lakes <- grid_OBBA2$great_lakes
surveys2$great_lakes <- extract_mean(great_lakes, surveys2, boundary_buf)
surveys3$great_lakes <- extract_mean(great_lakes, surveys3, boundary_buf)

# Coastline
grid_OBBA2$coastline <- extract_mean(coastline, grid_cells, boundary_buf)
grid_OBBA3$coastline <- grid_OBBA2$coastline
surveys2$coastline <- extract_mean(coastline, surveys2, boundary_buf)
surveys3$coastline <- extract_mean(coastline, surveys3, boundary_buf)

# Coastline
grid_OBBA2$open_water <- extract_mean(open_water, grid_cells, boundary_buf)
grid_OBBA3$open_water <- grid_OBBA2$open_water
surveys2$open_water <- extract_mean(open_water, surveys2, boundary_buf)
surveys3$open_water <- extract_mean(open_water, surveys3, boundary_buf)

# Keep geometry last
grid_OBBA2 <- grid_OBBA2 %>% relocate(geometry, .after = last_col())
grid_OBBA3 <- grid_OBBA3 %>% relocate(geometry, .after = last_col())

# Recombine surveys2 and surveys3 into sf point object
all_surveys_with_covs <- bind_rows(surveys2,surveys3) %>%
  arrange(obs_idx) %>%
  st_centroid()

# ------------------------------------------------------------
# Assign each survey to an atlas square
# ------------------------------------------------------------

if (file.exists(in_atlas_squares)) {
  atlas_squares <- st_read(in_atlas_squares, quiet = TRUE) %>%
    st_make_valid() %>%
    st_transform(st_crs(all_surveys_with_covs)) %>%
    select(square_id)
  
  hit <- st_within(all_surveys_with_covs, atlas_squares)
  
  square_id <- rep(NA_character_, nrow(all_surveys_with_covs))
  has_hit <- lengths(hit) > 0L
  if (any(has_hit)) {
    idx <- vapply(hit[has_hit], function(x) x[1], integer(1))
    square_id[has_hit] <- atlas_squares$square_id[idx]
  }
  all_surveys_with_covs$square_id <- square_id
} else {
  message("Atlas squares file not found; skipping square_id assignment.")
  all_surveys_with_covs$square_id <- NA_character_
}

# ------------------------------------------------------------
# Crop prediction grids to exact study area boundary
# ------------------------------------------------------------

st_crs(study_area$boundary) == st_crs(grid_OBBA2)
bndry <- study_area$boundary %>% st_transform(st_crs(grid_OBBA2))
grid_OBBA2 <- st_filter(
  grid_OBBA2,
  bndry,
  .predicate = st_intersects
)

grid_OBBA3 <- st_filter(
  grid_OBBA3,
  bndry,
  .predicate = st_intersects
)


# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

# Return surveys to original CRS (km units) for downstream scripts
all_surveys_with_covs <- st_transform(all_surveys_with_covs, st_crs(all_surveys)) %>%
  relocate(geometry, .after = last_col())

saveRDS(
  list(
    all_surveys_with_covs = all_surveys_with_covs,
    count_matrix = count_matrix,
    grid_OBBA2 = grid_OBBA2,   # EPSG:3978
    grid_OBBA3 = grid_OBBA3,   # EPSG:3978
    grid_centroids = grid_centroids,
    grid_cells = grid_cells,
    all_species = all_species,
    boundary = study_area$boundary,
    boundary_buffer_5km = study_area$boundary_buffer_5km,
    date_created = Sys.time()
  ),
  file = out_file
)

message("05_extract_covariates.R complete.")