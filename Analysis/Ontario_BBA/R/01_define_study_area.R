# ============================================================
# 01_define_study_area.R
#
# Purpose:
#   - Define the spatial domain of the Ontario BBA analysis
#   - Establish the canonical CRS (AEA, km units)
#   - Create study boundary and buffers
#
# Outputs:
#   data_clean/spatial/study_area.rds
# ============================================================

library(sf)
library(dplyr)

# ------------------------------------------------------------
# Load spatial utility functions
# ------------------------------------------------------------

source("R/functions/spatial_utils.R")

crs_aea_km <- get_aea_km_crs()

# ------------------------------------------------------------
# Load raw boundary data
# ------------------------------------------------------------

# Ontario boundary
ontario <- st_read("../../Data/Spatial/BCR/BCR_Terrestrial_master.shp") %>%
  subset(PROVINCE_S == "ONTARIO") %>%
  st_make_valid() %>%
  dplyr::select(BCR, PROVINCE_S) %>%
  st_union() %>%
  st_transform(crs_aea_km)

# ------------------------------------------------------------
# Create buffers
# ------------------------------------------------------------

# 25 km buffer for covariate extraction near edges
ontario_buffer_25km <- st_buffer(ontario, dist = 25)
  
# union geometry to avoid multipart issues downstream
ontario <- st_union(ontario)
ontario_buffer_25km <- st_union(ontario_buffer_25km)

# ------------------------------------------------------------
# Sanity checks
# ------------------------------------------------------------

stopifnot(st_is_valid(ontario))
stopifnot(st_is_valid(ontario_buffer_25km))
stopifnot(st_crs(ontario) == crs_aea_km)

# ------------------------------------------------------------
# Save study area object
# ------------------------------------------------------------

study_area <- list(
  crs = crs_aea_km,
  boundary = ontario,
  boundary_buffer_25km = ontario_buffer_25km
)

dir.create("data_clean/spatial", recursive = TRUE, showWarnings = FALSE)

saveRDS(
  study_area,
  file = "data_clean/spatial/study_area.rds"
)
