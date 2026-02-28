# ============================================================
# 03_download_raw_covariates.R
#
# Downloads RAW covariate rasters (cache) from Google Earth Engine.
# No processing/compositing here.
#
# Outputs:
#   ../../Data/Spatial/MODIS/modis_lc_<year>.tif
#   ../../Data/Spatial/GHSL/ghsl_urban_<year>.tif
#   ../../Data/Spatial/metadata/ee_download_receipt_03.rds
# ============================================================

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(terra)
  library(rgee)
})

source("R/functions/spatial_utils.R")
source("R/functions/covariate_processing_utils.R") # contains ee_download_tif()

# ------------------------------------------------------------
# Config (user-editable)
# ------------------------------------------------------------

years_modis_obba2 <- 2001:2005
years_modis_obba3 <- 2020:2024
years_ghsl <- c(2000, 2020)

years_modis <- sort(unique(c(years_modis_obba2, years_modis_obba3)))

out_dir_modis <- "../../Data/Spatial/MODIS"
out_dir_ghsl  <- "../../Data/Spatial/GHSL"
out_dir_meta  <- "../../Data/Spatial/metadata"

dir.create(out_dir_modis, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_ghsl,  recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_meta,  recursive = TRUE, showWarnings = FALSE)

# Export CRS (meters) for EE downloads (stable and standard).
# Processing to km CRS happens in step 04.
export_crs <- "EPSG:3978"

# Typical nominal scales
scale_modis_m <- 500   # MCD12Q1 nominal 500 m
scale_ghsl_m  <- 1000  # GHSL SMOD nominal 1 km

# ------------------------------------------------------------
# Load study area (use buffered boundary for downloads)
# ------------------------------------------------------------

study_area <- readRDS("data_clean/spatial/study_area.rds")

# rgee needs lon/lat geometry for sf_as_ee()
download_boundary_ll <- study_area$boundary_buffer_25km %>%
  st_transform(4326) %>%
  st_make_valid()

# ------------------------------------------------------------
# Initialize Earth Engine
# ------------------------------------------------------------

tryCatch(
  ee_Initialize(),
  error = function(e) {
    stop(
      "Earth Engine initialization failed.\n",
      "Common fixes:\n",
      "  - Run ee_Initialize() interactively once to authenticate\n",
      "  - Ensure you have Google Drive enabled/authorized if using via='drive'\n",
      "Original error:\n", conditionMessage(e)
    )
  }
)

download_region_ee <- rgee::sf_as_ee(download_boundary_ll)

# ------------------------------------------------------------
# MODIS MCD12Q1 (Land Cover Type 1: IGBP)
# Collection: MODIS/061/MCD12Q1
# Band: LC_Type1
# ------------------------------------------------------------

modis_collection_id <- "MODIS/061/MCD12Q1"
modis_band <- "LC_Type1"
modis_ic <- ee$ImageCollection(modis_collection_id)

for (yr in years_modis) {
  out_file <- file.path(out_dir_modis, sprintf("modis_lc_%d.tif", yr))
  if (file.exists(out_file)) next
  
  message("Downloading MODIS land cover year: ", yr)
  
  img_year <- modis_ic$
    filterDate(paste0(yr, "-01-01"), paste0(yr, "-12-31"))$
    select(modis_band)
  
  # Ensure the year exists in the collection
  n_img <- img_year$size()$getInfo()
  if (is.null(n_img) || n_img < 1) {
    warning("No MODIS image found for year ", yr, " (skipping).")
    next
  }
  
  img <- img_year$
    first()$
    clip(download_region_ee)
  
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
# Collection: JRC/GHSL/P2023A/GHS_SMOD_V2-0
# Band: smod_code
# ------------------------------------------------------------

ghsl_collection_id <- "JRC/GHSL/P2023A/GHS_SMOD_V2-0"
ghsl_band <- "smod_code"
ghsl_ic <- ee$ImageCollection(ghsl_collection_id)

for (yr in years_ghsl) {
  out_file <- file.path(out_dir_ghsl, sprintf("ghsl_urban_%d.tif", yr))
  if (file.exists(out_file)) next
  
  message("Downloading GHSL SMOD year: ", yr)
  
  img_year <- ghsl_ic$
    filterDate(paste0(yr, "-01-01"), paste0(yr, "-12-31"))$
    select(ghsl_band)
  
  n_img <- img_year$size()$getInfo()
  if (is.null(n_img) || n_img < 1) {
    warning("No GHSL SMOD image found for year ", yr, " (skipping).")
    next
  }
  
  img <- img_year$
    first()$
    clip(download_region_ee)
  
  ee_download_tif(
    ee_image = img,
    region   = download_region_ee,
    out_file = out_file,
    scale_m  = scale_ghsl_m,
    crs      = export_crs,
    via      = "drive"
  )
}

# ------------------------------------------------------------
# Save a small receipt (reproducibility)
# ------------------------------------------------------------

receipt <- list(
  script = "03_download_raw_covariates.R",
  run_time_utc = format(Sys.time(), tz = "UTC"),
  export_crs = export_crs,
  years_modis = years_modis,
  years_ghsl = years_ghsl,
  modis = list(collection = modis_collection_id, band = modis_band, scale_m = scale_modis_m),
  ghsl  = list(collection = ghsl_collection_id, band = ghsl_band,  scale_m = scale_ghsl_m),
  boundary_bbox_ll = sf::st_bbox(download_boundary_ll)
)

saveRDS(receipt, file.path(out_dir_meta, "ee_download_receipt_03.rds"))

message("03_download_raw_covariates.R complete.")