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
#   - The km-units Lambert CRS is used later for INLA/SPDE mesh coordinates and sf objects,
#     but rasters are kept in metres (EPSG:3978) for robust GDAL/terra reprojection.
# ============================================================

library(sf)
library(dplyr)
library(terra)

source("R/functions/spatial_utils.R")
source("R/functions/covariate_processing_utils.R")

dir.create("data_clean/spatial", recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Config
# ------------------------------------------------------------

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
worldclim_dir <- "../../Data/Spatial/WorldClim"  # expects /prec and /tmax subfolders as before

# ------------------------------------------------------------
# Study area + raster processing CRS (EPSG:3978)
# ------------------------------------------------------------

study_area <- readRDS("data_clean/spatial/study_area.rds")
boundary_buf <- study_area$boundary_buffer_25km
crs_aea_km <- study_area$crs

# Terra wants WKT; sf::st_crs works fine too
terra_crs_raster_m <- sf::st_crs(3978)$wkt

# ------------------------------------------------------------
# 1) MODIS composites (mode across years)
# ------------------------------------------------------------

modis_files_obba2 <- file.path(raw_modis_dir, sprintf("modis_lc_%d.tif", years_modis_obba2))
modis_files_obba3 <- file.path(raw_modis_dir, sprintf("modis_lc_%d.tif", years_modis_obba3))

stopifnot(all(file.exists(modis_files_obba2)))
stopifnot(all(file.exists(modis_files_obba3)))

rasters_obba2 <- rast(modis_files_obba2)
rasters_obba3 <- rast(modis_files_obba3)

# Reproject to canonical CRS (categorical -> nearest)
rasters_obba2 <- project(rasters_obba2, terra_crs_raster_m, method = "near")
rasters_obba3 <- project(rasters_obba3, terra_crs_raster_m, method = "near")

# Crop/mask to buffered boundary
rasters_obba2 <- crop_mask_to_boundary(rasters_obba2, boundary_buf)
rasters_obba3 <- crop_mask_to_boundary(rasters_obba3, boundary_buf)

mode_obba2 <- terra::modal(rasters_obba2, ties = "first", na.rm = TRUE)
mode_obba3 <- terra::modal(rasters_obba3, ties = "first", na.rm = TRUE)

writeRaster(mode_obba2, "data_clean/spatial/MCD12Q1_OBBA2_mode.tif", overwrite = TRUE)
writeRaster(mode_obba3, "data_clean/spatial/MCD12Q1_OBBA3_mode.tif", overwrite = TRUE)

# ------------------------------------------------------------
# 2) GHSL reclass (Urban_2000 / Urban_2020)
# ------------------------------------------------------------

ghsl_2000 <- rast(file.path(raw_ghsl_dir, "ghsl_urban_2000.tif"))
ghsl_2020 <- rast(file.path(raw_ghsl_dir, "ghsl_urban_2020.tif"))

ghsl_2000 <- project(ghsl_2000, terra_crs_raster_m, method = "near") |> crop_mask_to_boundary(boundary_buf)
ghsl_2020 <- project(ghsl_2020, terra_crs_raster_m, method = "near") |> crop_mask_to_boundary(boundary_buf)

ghsl_2000_reclass <- reclass_ghsl(ghsl_2000)
ghsl_2020_reclass <- reclass_ghsl(ghsl_2020)

writeRaster(ghsl_2000_reclass, "data_clean/spatial/Urban_2000.tif", overwrite = TRUE)
writeRaster(ghsl_2020_reclass, "data_clean/spatial/Urban_2020.tif", overwrite = TRUE)

# ------------------------------------------------------------
# 3) Climate means (WorldClim 2.1)
# ------------------------------------------------------------

mean_prec_obba2 <- mean_climate(
  variable = "prec",
  years = years_prec_obba2,
  months = months_prec,
  base_dir = worldclim_dir,
  boundary = boundary_buf,
  target_crs_wkt = terra_crs_raster_m
)

mean_prec_obba3 <- mean_climate(
  variable = "prec",
  years = years_prec_obba3,
  months = months_prec,
  base_dir = worldclim_dir,
  boundary = boundary_buf,
  target_crs_wkt = terra_crs_raster_m
)

mean_tmax_obba2 <- mean_climate(
  variable = "tmax",
  years = years_tmax_obba2,
  months = months_tmax,
  base_dir = worldclim_dir,
  boundary = boundary_buf,
  target_crs_wkt = terra_crs_raster_m
)

mean_tmax_obba3 <- mean_climate(
  variable = "tmax",
  years = years_tmax_obba3,
  months = months_tmax,
  base_dir = worldclim_dir,
  boundary = boundary_buf,
  target_crs_wkt = terra_crs_raster_m
)

writeRaster(mean_prec_obba2, "data_clean/spatial/prec_OBBA2.tif", overwrite = TRUE)
writeRaster(mean_prec_obba3, "data_clean/spatial/prec_OBBA3.tif", overwrite = TRUE)
writeRaster(mean_tmax_obba2, "data_clean/spatial/tmax_OBBA2.tif", overwrite = TRUE)
writeRaster(mean_tmax_obba3, "data_clean/spatial/tmax_OBBA3.tif", overwrite = TRUE)

# ------------------------------------------------------------
# 4) Additional layers (roads, insects, water)
#    Outputs are in EPSG:3978 (meters)
# ------------------------------------------------------------

# -------------------------
# 4a) Roads: within 100 m
# -------------------------

roads_2005_src <- "../../Data/Spatial/RoadNetwork/2005/grnf035r05a_e.shp"
roads_2025_src <- "../../Data/Spatial/RoadNetwork/2025/lrnf000r25a_e.shp"

if (file.exists(roads_2005_src) && file.exists(roads_2025_src)) {
  
  template_100m <- make_template_raster_m(
    boundary_sf = boundary_buf,
    res_m = 100,
    crs_wkt = terra_crs_raster_m
  )
  
  roadside_2005 <- rasterize_roads_buffer_presence(
    roads_path = roads_2005_src,
    boundary_sf = boundary_buf,
    template = template_100m,
    buffer_m = 100
  )
  
  roadside_2025 <- rasterize_roads_buffer_presence(
    roads_path = roads_2025_src,
    boundary_sf = boundary_buf,
    template = template_100m,
    buffer_m = 100
  )
  
  terra::writeRaster(roadside_2005, "data_clean/spatial/roadside_2005_100m.tif", overwrite = TRUE)
  terra::writeRaster(roadside_2025, "data_clean/spatial/roadside_2025_100m.tif", overwrite = TRUE)
}

# -------------------------
# 4b) Insects: areas defoliated by insect type
# -------------------------

insect_src <- "../../Data/Spatial/Insect_Damage/FOREST_INSECT_DAMAGE_EVENT.shp"

if (file.exists(insect_src)) {
  
  template_100m <- make_template_raster_m(
    boundary_sf = boundary_buf,
    res_m = 100,
    crs_wkt = terra_crs_raster_m
  )
  
  insect_OBBA2 <- rasterize_insect_presence_by_type(
    insect_path = insect_src,
    boundary_sf = boundary_buf,
    template = template_100m,
    years = 2000:2004,
    insects_keep = c("Gypsy Moth","Spruce Budworm","Forest Tent Caterpillar","Jack Pine Budworm"),
    ranking_exclude = "Light",
    simplify_tolerance_m = 100
  )
  
  insect_OBBA3 <- rasterize_insect_presence_by_type(
    insect_path = insect_src,
    boundary_sf = boundary_buf,
    template = template_100m,
    years = 2020:2024,
    insects_keep = c("Gypsy Moth","Spruce Budworm","Forest Tent Caterpillar","Jack Pine Budworm"),
    ranking_exclude = "Light",
    simplify_tolerance_m = 100
  )
  
  terra::writeRaster(insect_OBBA2, "data_clean/spatial/insect_OBBA2.tif", overwrite = TRUE)
  terra::writeRaster(insect_OBBA3, "data_clean/spatial/insect_OBBA3.tif", overwrite = TRUE)
}

# -------------------------
# 4c) Water: open water vs river + large-water polygons for masking/plotting
# -------------------------

water_src <- "../../Data/Spatial/Ontario_Hydro_Network_(OHN)_-_Waterbody/Ontario_Hydro_Network_(OHN)_-_Waterbody.shp"

if (file.exists(water_src)) {
  
  water_layers <- process_ohn_water(
    water_path = water_src,
    boundary_sf = boundary_buf,
    crs_wkt = terra_crs_raster_m,
    raster_res_m = 30,
    simplify_tolerance_m = 10,
    min_area_m2_for_filtered = 1e6
  )
  
  terra::writeRaster(water_layers$water_open,  "data_clean/spatial/water_open.tif",  overwrite = TRUE)
  terra::writeRaster(water_layers$water_river, "data_clean/spatial/water_river.tif", overwrite = TRUE)
  
  sf::st_write(
    water_layers$water_filtered,
    "data_clean/spatial/water_filtered.shp",
    delete_layer = TRUE,
    quiet = TRUE
  )
}
