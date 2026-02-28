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

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
})

# ------------------------------------------------------------
# Load spatial utility functions
# ------------------------------------------------------------

spatial_utils_path <- "R/functions/spatial_utils.R"
if (!file.exists(spatial_utils_path)) {
  stop("Cannot find spatial utils at: ", spatial_utils_path,
       "\nUpdate 'spatial_utils_path' to the correct location.")
}
source(spatial_utils_path)

crs_aea_km <- get_aea_km_crs()

# ------------------------------------------------------------
# Load raw boundary data
# ------------------------------------------------------------

bcr_path <- "../../Data/Spatial/BCR/BCR_Terrestrial_master.shp"
if (!file.exists(bcr_path)) {
  stop("Cannot find boundary shapefile at: ", bcr_path)
}

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
# projected and uses km units (so dist = 25 means 25 km).
ontario_buffer_25km <- st_buffer(ontario, dist = 25) %>%
  st_make_valid()

# ------------------------------------------------------------
# Sanity checks
# ------------------------------------------------------------

if (sf::st_is_longlat(ontario)) {
  stop("Ontario geometry is in lon/lat after transform; expected a projected CRS.")
}

stopifnot(st_is_valid(ontario))
stopifnot(st_is_valid(ontario_buffer_25km))

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
  boundary_buffer_25km = ontario_buffer_25km
)

out_dir <- "data_clean/spatial"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

saveRDS(study_area, file = file.path(out_dir, "study_area.rds"))

message("01_define_study_area.R complete.")