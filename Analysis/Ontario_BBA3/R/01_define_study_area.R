# ============================================================
# 01_define_study_area.R
#
# Purpose:
#   - Define the spatial domain of the Ontario BBA analysis
#   - Establish the canonical CRS (AEA, intended km units)
#   - Create study boundary and buffers for edge-safe processing
#
# Inputs:
#   - ../../Data/Spatial/BCR/BCR_Terrestrial_master.shp
#   - functions: R/functions/spatial_utils.R (or update path if different)
#
# Outputs:
#   - data_clean/spatial/study_area.rds
#     A list with elements:
#       $crs
#       $boundary
#       $boundary_buffer_25km
# ============================================================

rm(list=ls())

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(here)
})

# Centralized paths; loads an object called "paths" that has paths to each subfolder
source(here::here("R", "00_config_paths.R"))

# ------------------------------------------------------------
# Load spatial utility functions
# ------------------------------------------------------------

source(file.path(paths$functions, "spatial_utils.R"))
crs_aea_km <- get_aea_km_crs()

# ------------------------------------------------------------
# Load raw boundary data
# ------------------------------------------------------------

bcr_path <- file.path(paths$data, "Spatial", "BCR", "BCR_Terrestrial_master.shp")
bcr <- st_read(bcr_path, quiet = TRUE)

required_cols <- c("PROVINCE_S", "BCR")
missing_cols <- setdiff(required_cols, names(bcr))
if (length(missing_cols) > 0) {
  stop("Boundary file is missing required columns: ",
       paste(missing_cols, collapse = ", "))
}

# Ontario boundary (single polygon geometry)
ontario <- bcr %>%
  filter(PROVINCE_S == "ONTARIO") %>%
  st_make_valid() %>%
  select(BCR, PROVINCE_S) %>%
  st_union() %>%                 # dissolve into a single geometry
  st_transform(crs_aea_km)

# ------------------------------------------------------------
# Create buffers
# ------------------------------------------------------------

# Buffer distance is in the units of the CRS. This script assumes the CRS is
# projected and uses km units (so dist = 5 means 5 km).
ontario_buffer_5km <- st_buffer(ontario, dist = 5) %>%
  st_make_valid()

# ------------------------------------------------------------
# Sanity checks
# ------------------------------------------------------------

if (sf::st_is_longlat(ontario)) {
  stop("Ontario geometry is in lon/lat after transform; expected a projected CRS.")
}

stopifnot(st_is_valid(ontario))
stopifnot(st_is_valid(ontario_buffer_5km))

# Be explicit about CRS equivalence
if (is.na(st_crs(ontario)$wkt) || is.na(crs_aea_km$wkt)) {
  stop("CRS WKT missing; cannot verify CRS equivalence reliably.")
}
stopifnot(st_crs(ontario) == crs_aea_km)

# ------------------------------------------------------------
# Save study area object
# ------------------------------------------------------------

study_area <- list(
  crs = crs_aea_km,
  boundary = ontario,
  boundary_buffer_5km = ontario_buffer_5km
)

out_dir <- file.path(paths$data_clean, "spatial")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

saveRDS(study_area, file = file.path(out_dir, "study_area.rds"))

message("01_define_study_area.R complete.")