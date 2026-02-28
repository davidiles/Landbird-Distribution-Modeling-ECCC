# ============================================================
# 04_process_covariates.R
#
# Purpose:
#   Convert raw raster downloads (Earth Engine / WorldClim) and local vector layers
#   (roads, insect damage, hydrography) into analysis-ready covariates clipped to the
#   Ontario study area (with buffer) and stored in a single raster CRS for consistent
#   extraction and prediction.
#
# Raster CRS used in this script:
#   EPSG:3978 (Canada Albers, metres)
#
# Notes:
#   - In this step, covariates are reprojected + clipped, but kept at their native resolution.
#     Harmonization to a common grid and covariate extraction happens in step 06.
#   - The km-units CRS is used later for INLA/SPDE mesh coordinates and sf objects.
# ============================================================

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(terra)
})

source("R/functions/spatial_utils.R")
source("R/functions/covariate_processing_utils.R")

dir.create("data_clean/spatial", recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Config (user-editable)
# ------------------------------------------------------------

# MODIS settings
years_modis_obba2 <- 2001:2005
years_modis_obba3 <- 2020:2024

# Climate settings
years_prec_obba2 <- 2000:2004
years_prec_obba3 <- 2020:2024
months_prec <- 3:7

years_tmax_obba2 <- 2000:2004
years_tmax_obba3 <- 2020:2024
months_tmax <- 3:7

# Input dirs
raw_modis_dir <- "../../Data/Spatial/MODIS"
raw_ghsl_dir  <- "../../Data/Spatial/GHSL"
worldclim_dir <- "../../Data/Spatial/WorldClim"  # expects /prec and /tmax subfolders

# Extra layers
roads_2005_src <- "../../Data/Spatial/RoadNetwork/2005/grnf035r05a_e.shp"
roads_2025_src <- "../../Data/Spatial/RoadNetwork/2025/lrnf000r25a_e.shp"
insect_src     <- "../../Data/Spatial/Insect_Damage/FOREST_INSECT_DAMAGE_EVENT.shp"
water_src      <- "../../Data/Spatial/Ontario_Hydro_Network_(OHN)_-_Waterbody/Ontario_Hydro_Network_(OHN)_-_Waterbody.shp"

# Canonical raster CRS (metres)
crs_raster_m <- sf::st_crs(3978)$wkt

# ------------------------------------------------------------
# Study area (buffered boundary)
# ------------------------------------------------------------

study_area   <- readRDS("data_clean/spatial/study_area.rds")
boundary_buf <- study_area$boundary_buffer_25km

# Work in EPSG:3978 inside this script
boundary_buf_m <- st_transform(boundary_buf, st_crs(3978))

# A reusable 100 m template for vector-derived layers (roads/insects)
template_100m <- make_template_raster_m(
  boundary_sf = boundary_buf_m,
  res_m = 100,
  crs_wkt = crs_raster_m
)

# ------------------------------------------------------------
# 1) MODIS composites (mode across years)
# ------------------------------------------------------------

need_files <- function(paths, label) {
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0) {
    stop(
      "Missing required files for ", label, ":\n  - ",
      paste(missing, collapse = "\n  - "),
      "\n\nDid you run 03_download_raw_covariates.R first?"
    )
  }
}

modis_files_obba2 <- file.path(raw_modis_dir, sprintf("modis_lc_%d.tif", years_modis_obba2))
modis_files_obba3 <- file.path(raw_modis_dir, sprintf("modis_lc_%d.tif", years_modis_obba3))

need_files(modis_files_obba2, "MODIS OBBA2")
need_files(modis_files_obba3, "MODIS OBBA3")

rasters_obba2 <- rast(modis_files_obba2)
rasters_obba3 <- rast(modis_files_obba3)

# Reproject to canonical raster CRS (categorical -> nearest)
rasters_obba2 <- project(rasters_obba2, crs_raster_m, method = "near")
rasters_obba3 <- project(rasters_obba3, crs_raster_m, method = "near")

# Crop/mask to buffered boundary
rasters_obba2 <- crop_mask_to_boundary(rasters_obba2, boundary_buf_m)
rasters_obba3 <- crop_mask_to_boundary(rasters_obba3, boundary_buf_m)

mode_obba2 <- terra::modal(rasters_obba2, ties = "first", na.rm = TRUE)
mode_obba3 <- terra::modal(rasters_obba3, ties = "first", na.rm = TRUE)

writeRaster(mode_obba2, "data_clean/spatial/MODIS_LC_OBBA2_mode.tif", overwrite = TRUE)
writeRaster(mode_obba3, "data_clean/spatial/MODIS_LC_OBBA3_mode.tif", overwrite = TRUE)

# ------------------------------------------------------------
# 2) GHSL reclass (Urban_2000 / Urban_2020)
# ------------------------------------------------------------

ghsl_2000_path <- file.path(raw_ghsl_dir, "ghsl_urban_2000.tif")
ghsl_2020_path <- file.path(raw_ghsl_dir, "ghsl_urban_2020.tif")
need_files(c(ghsl_2000_path, ghsl_2020_path), "GHSL")

ghsl_2000 <- rast(ghsl_2000_path)
ghsl_2020 <- rast(ghsl_2020_path)

# IMPORTANT: explicit crop_mask_to_boundary call (avoids pipe-arg confusion)
ghsl_2000 <- crop_mask_to_boundary(project(ghsl_2000, crs_raster_m, method = "near"), boundary_buf_m)
ghsl_2020 <- crop_mask_to_boundary(project(ghsl_2020, crs_raster_m, method = "near"), boundary_buf_m)

ghsl_2000_reclass <- reclass_ghsl(ghsl_2000)
ghsl_2020_reclass <- reclass_ghsl(ghsl_2020)

writeRaster(ghsl_2000_reclass, "data_clean/spatial/Urban_2000.tif", overwrite = TRUE)
writeRaster(ghsl_2020_reclass, "data_clean/spatial/Urban_2020.tif", overwrite = TRUE)

# ------------------------------------------------------------
# 3) Climate means (WorldClim 2.1)
# ------------------------------------------------------------

mean_prec_obba2 <- mean_climate("prec", years_prec_obba2, months_prec, worldclim_dir, boundary_buf_m, crs_raster_m)
mean_prec_obba3 <- mean_climate("prec", years_prec_obba3, months_prec, worldclim_dir, boundary_buf_m, crs_raster_m)

mean_tmax_obba2 <- mean_climate("tmax", years_tmax_obba2, months_tmax, worldclim_dir, boundary_buf_m, crs_raster_m)
mean_tmax_obba3 <- mean_climate("tmax", years_tmax_obba3, months_tmax, worldclim_dir, boundary_buf_m, crs_raster_m)

writeRaster(mean_prec_obba2, "data_clean/spatial/prec_OBBA2.tif", overwrite = TRUE)
writeRaster(mean_prec_obba3, "data_clean/spatial/prec_OBBA3.tif", overwrite = TRUE)
writeRaster(mean_tmax_obba2, "data_clean/spatial/tmax_OBBA2.tif", overwrite = TRUE)
writeRaster(mean_tmax_obba3, "data_clean/spatial/tmax_OBBA3.tif", overwrite = TRUE)

# ------------------------------------------------------------
# 4) Additional layers (roads, insects, water) in EPSG:3978
# ------------------------------------------------------------

# 4a) Roads: within 100 m
if (file.exists(roads_2005_src) && file.exists(roads_2025_src)) {
  roadside_2005 <- rasterize_roads_buffer_presence(roads_2005_src, boundary_buf_m, template_100m, buffer_m = 100)
  roadside_2025 <- rasterize_roads_buffer_presence(roads_2025_src, boundary_buf_m, template_100m, buffer_m = 100)
  
  writeRaster(roadside_2005, "data_clean/spatial/roadside_2005_100m.tif", overwrite = TRUE)
  writeRaster(roadside_2025, "data_clean/spatial/roadside_2025_100m.tif", overwrite = TRUE)
} else {
  message("Road layers not found; skipping roads covariates.")
}

# 4b) Insects: areas defoliated by insect type
if (file.exists(insect_src)) {
  insect_OBBA2 <- rasterize_insect_presence_by_type(
    insect_path = insect_src,
    boundary_sf = boundary_buf_m,
    template = template_100m,
    years = 2000:2004,
    insects_keep = c("Gypsy Moth", "Spruce Budworm", "Forest Tent Caterpillar", "Jack Pine Budworm"),
    ranking_exclude = "Light",
    simplify_tolerance_m = 100
  )
  
  insect_OBBA3 <- rasterize_insect_presence_by_type(
    insect_path = insect_src,
    boundary_sf = boundary_buf_m,
    template = template_100m,
    years = 2020:2024,
    insects_keep = c("Gypsy Moth", "Spruce Budworm", "Forest Tent Caterpillar", "Jack Pine Budworm"),
    ranking_exclude = "Light",
    simplify_tolerance_m = 100
  )
  
  writeRaster(insect_OBBA2, "data_clean/spatial/insect_OBBA2.tif", overwrite = TRUE)
  writeRaster(insect_OBBA3, "data_clean/spatial/insect_OBBA3.tif", overwrite = TRUE)
} else {
  message("Insect damage layer not found; skipping insect covariates.")
}

# 4c) Water: open water vs river + large-water polygons for masking/plotting
if (file.exists(water_src)) {
  water_layers <- process_ohn_water(
    water_path = water_src,
    boundary_sf = boundary_buf_m,
    crs_wkt = crs_raster_m,
    raster_res_m = 30,
    river_buffer_m = 100,
    simplify_tolerance_m = 10,
    min_area_m2_for_filtered = 1e6
  )
  
  writeRaster(water_layers$water_open,  "data_clean/spatial/water_open.tif",  overwrite = TRUE)
  writeRaster(water_layers$water_river, "data_clean/spatial/water_river.tif", overwrite = TRUE)
  
  st_write(
    water_layers$water_filtered,
    "data_clean/spatial/water_filtered.shp",
    delete_layer = TRUE,
    quiet = TRUE
  )
} else {
  message("OHN water layer not found; skipping water covariates.")
}

message("04_process_covariates.R complete.")