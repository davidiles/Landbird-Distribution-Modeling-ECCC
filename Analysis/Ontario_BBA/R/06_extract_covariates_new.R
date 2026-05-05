# ============================================================
# 06_extract_covariates.R
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

# ------------------------------------------------------------
# Config
# ------------------------------------------------------------

crs_rast <- sf::st_crs(3978)   # EPSG:3978

cov_dir  <- file.path(paths$data_clean, "spatial")
out_file <- file.path(paths$data_clean, "birds", "analysis_data_covariates.rds")

in_grid         <- file.path(paths$data_clean, "spatial",  "prediction_grid_1000_m.rds")

in_surveys      <- file.path(paths$data_clean, "surveys",  "surveys_raw.rds")
in_count_matrix <- file.path(paths$data_clean, "surveys",  "count_matrix_raw.rds")
in_species      <- file.path(paths$data_clean, "metadata", "species_list.rds")

in_study_area   <- file.path(paths$data_clean, "spatial",  "study_area.rds")

# Optional national layers (shared Data/)
in_atlas_squares <- file.path(paths$data, "Spatial", "National", "AtlasSquares", "NationalSquares_FINAL.shp")
in_bcr           <- file.path(paths$data, "Spatial", "National", "BCR", "BCR_Terrestrial_master.shp")

# ------------------------------------------------------------
# Load inputs
# ------------------------------------------------------------

need_files(c(in_surveys, in_count_matrix, in_species, in_grid, in_study_area), "core inputs")

all_surveys  <- readRDS(in_surveys)
count_matrix <- readRDS(in_count_matrix)
all_species  <- readRDS(in_species)

grid_obj <- readRDS(in_grid)
grid_cells     <- grid_obj$grid_cells
grid_centroids <- grid_obj$grid_centroids

study_area <- readRDS(in_study_area)
boundary_buf_m <- study_area$boundary_buffer_25km %>%
  st_transform(crs_rast) %>%
  st_make_valid()

# Ensure grid is in extraction CRS
grid_cells     <- st_transform(grid_cells, crs_rast)
grid_centroids <- st_transform(grid_centroids, crs_rast)

# Stable ordering field
if (!("obs_idx" %in% names(all_surveys))) {
  all_surveys$obs_idx <- seq_len(nrow(all_surveys))
}

# Surveys in extraction CRS
all_surveys_m <- st_transform(all_surveys, crs_rast)

# ------------------------------------------------------------
# Load rasters
# ------------------------------------------------------------

# Required processed rasters (from 04_process_covariates.R)
req_obba2 <- c(
  r_path("prec_OBBA2.tif"),
  r_path("tmax_OBBA2.tif"),
  r_path("Urban_2000.tif"),
  r_path("MODIS_LC_OBBA2_mode.tif")
)

req_obba3 <- c(
  r_path("prec_OBBA3.tif"),
  r_path("tmax_OBBA3.tif"),
  r_path("Urban_2020.tif"),
  r_path("MODIS_LC_OBBA3_mode.tif")
)

need_files(req_obba2, "required OBBA2 rasters (run 04_process_covariates.R)")
need_files(req_obba3, "required OBBA3 rasters (run 04_process_covariates.R)")

prec2  <- rast(r_path("prec_OBBA2.tif"))
tmax2  <- rast(r_path("tmax_OBBA2.tif"))
urban2 <- rast(r_path("Urban_2000.tif"))
MODIS2    <- rast(r_path("MODIS_LC_OBBA2_mode.tif"))
road2  <- try_rast(r_path("roadside_2005_100m.tif"))
insect2 <- try_rast(r_path("insect_OBBA2.tif"))

prec3  <- rast(r_path("prec_OBBA3.tif"))
tmax3  <- rast(r_path("tmax_OBBA3.tif"))
urban3 <- rast(r_path("Urban_2020.tif"))
MODIS3    <- rast(r_path("MODIS_LC_OBBA3_mode.tif"))
road3  <- try_rast(r_path("roadside_2025_100m.tif"))
insect3 <- try_rast(r_path("insect_OBBA3.tif"))

# New layers
LCC2020 <- try_rast(file.path(paths$data,"Spatial","National","LandCoverCanada2020","landcover-2020-classification.tif"))
water_river <- rast(r_path("water_river.tif"))
water_lake_large <- rast(r_path("water_lake_large.tif"))
water_lake_small <- rast(r_path("water_lake_small.tif"))


# ------------------------------------------------------------
# Build grid covariate tables (OBBA2 / OBBA3)
# ------------------------------------------------------------

grid_OBBA2 <- grid_centroids
grid_OBBA3 <- grid_centroids

# Continuous means
grid_OBBA2$prec <- extract_mean(prec2, grid_cells, boundary_buf_m)
grid_OBBA2$tmax <- extract_mean(tmax2, grid_cells, boundary_buf_m)
grid_OBBA3$prec <- extract_mean(prec3, grid_cells, boundary_buf_m)
grid_OBBA3$tmax <- extract_mean(tmax3, grid_cells, boundary_buf_m)

# Roads (0/1) as mean = fraction of pixel within 100 m of a road
if (!is.null(road2)) grid_OBBA2$road <- extract_mean(road2, grid_cells, boundary_buf_m)
if (!is.null(road3)) grid_OBBA3$road <- extract_mean(road3, grid_cells, boundary_buf_m)

# MODIS LC fractions -> MODIS_#
MODIS_frac2 <- extract_frac_clean(MODIS2, grid_cells, boundary_buf_m, prefix = "MODIS", drop0 = TRUE, clean_names = FALSE)
MODIS_frac3 <- extract_frac_clean(MODIS3, grid_cells, boundary_buf_m, prefix = "MODIS", drop0 = TRUE, clean_names = FALSE)
if (!is.null(MODIS_frac2)) grid_OBBA2 <- bind_cols(grid_OBBA2, MODIS_frac2)
if (!is.null(MODIS_frac3)) grid_OBBA3 <- bind_cols(grid_OBBA3, MODIS_frac3)

# LCC2020 fractions -> LCC2020_#
LCC2020_frac <- extract_frac_clean(LCC2020, grid_cells, boundary_buf_m, prefix = "LCC2020", drop0 = TRUE, clean_names = FALSE)
if (!is.null(LCC2020_frac)){
  grid_OBBA2 <- bind_cols(grid_OBBA2, LCC2020_frac)
  grid_OBBA3 <- bind_cols(grid_OBBA3, LCC2020_frac)
}

# Urban fractions -> urban_#
urb_frac2 <- extract_frac_clean(urban2, grid_cells, boundary_buf_m, prefix = "urban", drop0 = TRUE, clean_names = FALSE)
urb_frac3 <- extract_frac_clean(urban3, grid_cells, boundary_buf_m, prefix = "urban", drop0 = TRUE, clean_names = FALSE)
if (!is.null(urb_frac2)) grid_OBBA2 <- bind_cols(grid_OBBA2, urb_frac2)
if (!is.null(urb_frac3)) grid_OBBA3 <- bind_cols(grid_OBBA3, urb_frac3)

# Insects -> collapse to broadleaf/needleleaf then fractions, cleaned names
if (!is.null(insect2)) {
  insect2 <- crop_if(insect2, boundary_buf_m)
  insects2 <- collapse_insect_layers(insect2)
  frac_in2 <- suppressWarnings(exact_extract(insects2, match_grid_to_raster_crs(grid_cells, insects2), "frac"))
  frac_in2 <- drop_frac0_cols(frac_in2)
  frac_in2 <- clean_frac_names(frac_in2)
  grid_OBBA2 <- bind_cols(grid_OBBA2, frac_in2)
}
if (!is.null(insect3)) {
  insect3 <- crop_if(insect3, boundary_buf_m)
  insects3 <- collapse_insect_layers(insect3)
  frac_in3 <- suppressWarnings(exact_extract(insects3, match_grid_to_raster_crs(grid_cells, insects3), "frac"))
  frac_in3 <- drop_frac0_cols(frac_in3)
  frac_in3 <- clean_frac_names(frac_in3)
  grid_OBBA3 <- bind_cols(grid_OBBA3, frac_in3)
}

# # Water 

# proximity to small and large rivers
  river_frac <- extract_frac_clean(water_river, grid_cells, boundary_buf_m, prefix = "river", drop0 = TRUE, clean_names = FALSE)
  if (!is.null(river_frac)){
    grid_OBBA2 <- bind_cols(grid_OBBA2,river_frac)
    grid_OBBA3 <- bind_cols(grid_OBBA3,river_frac)
  }
  
# proximity to small and large lakes
lake_lg <- extract_frac_clean(water_lake_large, grid_cells, boundary_buf_m, prefix = "lake_lg", drop0 = TRUE, clean_names = FALSE)
grid_OBBA2 <- bind_cols(grid_OBBA2, lake_lg)
grid_OBBA3 <- bind_cols(grid_OBBA3, lake_lg)

lake_sm <- extract_frac_clean(water_lake_small, grid_cells, boundary_buf_m, prefix = "lake_sm", drop0 = TRUE, clean_names = FALSE)
grid_OBBA2 <- bind_cols(grid_OBBA2, lake_sm)
grid_OBBA3 <- bind_cols(grid_OBBA3, lake_sm)

# -----------------------
# ***********************
# To do: add proximity to great lakes and James/Hudson Bay
# -----------------------
# ***********************

# Keep geometry last
grid_OBBA2 <- grid_OBBA2 %>% relocate(geometry, .after = last_col())
grid_OBBA3 <- grid_OBBA3 %>% relocate(geometry, .after = last_col())

# ------------------------------------------------------------
# Assign surveys covariates based on their pixel membership (i.e., nearest grid centroid)
# ------------------------------------------------------------

if (!("Atlas" %in% names(all_surveys_m))) {
  stop("Survey object must contain an 'Atlas' column with values like 'OBBA2'/'OBBA3'.")
}

surveys2 <- all_surveys_m %>% filter(Atlas == "OBBA2")
surveys3 <- all_surveys_m %>% filter(Atlas == "OBBA3")

idx2 <- st_nearest_feature(surveys2, grid_OBBA2)
idx3 <- st_nearest_feature(surveys3, grid_OBBA3)

# Avoid duplicating pixel_id when binding
cov2 <- st_drop_geometry(grid_OBBA2[idx2, , drop = FALSE])
cov3 <- st_drop_geometry(grid_OBBA3[idx3, , drop = FALSE])

surveys2_covs <- bind_cols(surveys2, cov2)
surveys3_covs <- bind_cols(surveys3, cov3)

all_surveys_with_covs_m <- bind_rows(surveys2_covs, surveys3_covs) %>%
  arrange(obs_idx)

# ------------------------------------------------------------
# Assign each survey to an atlas square
# ------------------------------------------------------------

if (file.exists(in_atlas_squares)) {
  atlas_squares <- st_read(in_atlas_squares, quiet = TRUE) %>%
    st_make_valid() %>%
    st_transform(st_crs(all_surveys_with_covs_m)) %>%
    select(square_id)
  
  hit <- st_within(all_surveys_with_covs_m, atlas_squares)
  
  square_id <- rep(NA_character_, nrow(all_surveys_with_covs_m))
  has_hit <- lengths(hit) > 0L
  if (any(has_hit)) {
    idx <- vapply(hit[has_hit], function(x) x[1], integer(1))
    square_id[has_hit] <- atlas_squares$square_id[idx]
  }
  all_surveys_with_covs_m$square_id <- square_id
} else {
  message("Atlas squares file not found; skipping square_id assignment.")
  all_surveys_with_covs_m$square_id <- NA_character_
}

# ------------------------------------------------------------
# Assign each survey and grid cell to a BCR within Ontario
# ------------------------------------------------------------

if (file.exists(in_bcr)) {
  bcr_on <- st_read(in_bcr, quiet = TRUE) %>%
    st_make_valid() %>%
    select(BCR, PROVINCE_S) %>%
    filter(PROVINCE_S == "ONTARIO", BCR %in% c(7, 8, 12, 13)) %>%
    group_by(BCR) %>%
    summarise(geometry = st_union(geometry), .groups = "drop") %>%
    st_make_valid()
  
  all_surveys_with_covs_m$BCR <- assign_poly_id(all_surveys_with_covs_m, bcr_on, id_col = "BCR", nearest_fallback = TRUE)
  
  grid_bcr <- assign_poly_id(grid_centroids, bcr_on %>% st_transform(st_crs(grid_centroids)),
                             id_col = "BCR", nearest_fallback = TRUE)
  grid_OBBA2$BCR <- grid_bcr
  grid_OBBA3$BCR <- grid_bcr
} else {
  message("BCR shapefile not found; skipping BCR assignment.")
  all_surveys_with_covs_m$BCR <- NA_integer_
  grid_OBBA2$BCR <- NA_integer_
  grid_OBBA3$BCR <- NA_integer_
}

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

# Return surveys to original CRS (km units) for downstream scripts
all_surveys_with_covs <- st_transform(all_surveys_with_covs_m, st_crs(all_surveys)) %>%
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
    boundary_buffer_25km = study_area$boundary_buffer_25km,
    bcr_sf = bcr_on,
    date_created = Sys.time()
  ),
  file = out_file
)

message("06_extract_covariates.R complete.")
