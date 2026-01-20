# ============================================================
# 03_download_raw_covariates.R
#
# Downloads RAW covariate rasters (cache) from Google Earth Engine.
# No processing/compositing here.
#
# Outputs:
#   data/spatial/MODIS/modis_lc_<year>.tif
#   data/spatial/GHSL/ghsl_urban_<year>.tif
# ============================================================

library(sf)
library(dplyr)
library(terra)

# Earth Engine via rgee
library(rgee)

source("R/functions/spatial_utils.R")
source("R/functions/covariate_processing_utils.R")

# ------------------------------------------------------------
# Config
# ------------------------------------------------------------

years_modis_obba2 <- 2001:2005
years_modis_obba3 <- 2020:2024
years_ghsl <- c(2000, 2020)

out_dir_modis <- "../../Data/Spatial/MODIS"
out_dir_ghsl  <- "../../Data/Spatial/GHSL"

dir.create(out_dir_modis, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_ghsl,  recursive = TRUE, showWarnings = FALSE)

# Export CRS (meters) for EE downloads (stable and standard)
# Processing to km CRS happens in step 04.
export_crs <- "EPSG:3978"

# Typical native-ish scales
scale_modis_m <- 500   # MCD12Q1 nominal 500 m
scale_ghsl_m  <- 1000  # GHSL SMOD nominal 1 km

# ------------------------------------------------------------
# Load study area (use buffered boundary for downloads)
# ------------------------------------------------------------

study_area <- readRDS("data_clean/spatial/study_area.rds")
download_boundary <- study_area$boundary_buffer_25km %>% st_transform(4326)

# Convert to EE geometry
# Note: rgee expects sf in lon/lat for conversion.
ee_Initialize()  # assumes you’ve already authenticated in your environment
download_region_ee <- rgee::sf_as_ee(download_boundary)

# ------------------------------------------------------------
# MODIS MCD12Q1 (Land Cover Type 1: IGBP)
# Dataset: MODIS/061/MCD12Q1
# Band: LC_Type1
# ------------------------------------------------------------

modis_ic <- ee$ImageCollection("MODIS/061/MCD12Q1")

for (yr in c(years_modis_obba2, years_modis_obba3)) {
  out_file <- file.path(out_dir_modis, sprintf("modis_lc_%d.tif", yr))
  if (file.exists(out_file)) next
  
  message("Downloading MODIS land cover year: ", yr)
  
  img <- modis_ic$
    filterDate(paste0(yr, "-01-01"), paste0(yr, "-12-31"))$
    first()$
    select("LC_Type1")$
    clip(download_region_ee)$
    reproject(crs = export_crs, scale = scale_modis_m)
  
  ee_download_tif(
    ee_image = img,
    region   = download_region_ee,
    out_file = out_file,
    scale_m  = scale_modis_m,
    crs      = export_crs,
    via      = "drive"
  )
}

# ------------------------------------------------------------
# GHSL SMOD (Degree of Urbanization)
# Dataset: JRC/GHSL/P2023A/GHS_SMOD_V2-0
# Band: smod_code
# ------------------------------------------------------------

ghsl_ic <- ee$ImageCollection("JRC/GHSL/P2023A/GHS_SMOD_V2-0")

for (yr in years_ghsl) {
  out_file <- file.path(out_dir_ghsl, sprintf("ghsl_urban_%d.tif", yr))
  if (file.exists(out_file)) next
  
  message("Downloading GHSL SMOD year: ", yr)
  
  img <- ghsl_ic$
    filterDate(paste0(yr, "-01-01"), paste0(yr, "-12-31"))$
    first()$
    select("smod_code")$
    clip(download_region_ee)$
    reproject(crs = export_crs, scale = scale_ghsl_m)
  
  ee_download_tif(
    ee_image = img,
    region   = download_region_ee,
    out_file = out_file,
    scale_m  = scale_ghsl_m,
    crs      = export_crs,
    via      = "drive"
  )
}

message("03_download_raw_covariates.R complete.")