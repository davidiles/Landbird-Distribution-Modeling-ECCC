# ============================================================
# 06_extract_covariates.R
#
# Purpose:
#   1) Compute covariates on the prediction grid for OBBA2 and OBBA3
#      from processed rasters (EPSG:3978).
#   2) Attach covariates to survey points by nearest grid cell
#      (OBBA2 surveys -> OBBA2 grid covs; OBBA3 surveys -> OBBA3 grid covs).
#
# Inputs:
#   - data_clean/spatial/study_area.rds                (from 01_*)
#   - data_clean/birds/point_counts_clean.rds          (from 02_*)
#   - data_clean/spatial/*.tif                         (from 04_*)
#   - data_clean/spatial/prediction_grid.rds           (from 05_*)
#
# Outputs:
#   - data_clean/birds/analysis_data_covariates.rds
#     containing:
#       * all_surveys_with_covs (sf)
#       * count_matrix (matrix)
#       * grid_OBBA2 (sf points)
#       * grid_OBBA3 (sf points)
#       * study boundary objects
# ============================================================

library(sf)
library(dplyr)
library(terra)
library(exactextractr)

source("R/functions/spatial_utils.R")
source("R/functions/covariate_processing_utils.R")

dir.create("data_clean/birds",  recursive = TRUE, showWarnings = FALSE)
dir.create("data_clean/spatial", recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Config
# ------------------------------------------------------------

# Raster CRS (processing/extraction CRS)
crs_rast <- sf::st_crs(3978)     # EPSG:3978
crs_rast_wkt <- crs_rast$wkt

# Paths (adjust to match what 02_* writes)
in_surveys       <- "data_clean/surveys/surveys_raw.rds"
in_count_matrix  <- "data_clean/surveys/count_matrix_raw.rds"
in_species       <- "data_clean/metadata/species_list.rds"
in_grid          <- "data_clean/spatial/prediction_grid.rds"
in_study_area    <- "data_clean/spatial/study_area.rds"
in_atlas_squares <- "../../Data/Spatial/National/AtlasSquares/NationalSquares_FINAL.shp"
in_bcr           <- "../../Data/Spatial/National/BCR/BCR_Terrestrial_master.shp"

cov_dir          <- "data_clean/spatial"
out_file         <- "data_clean/birds/analysis_data_covariates.rds"

# ------------------------------------------------------------
# Load inputs
# ------------------------------------------------------------

stopifnot(file.exists(in_surveys))
stopifnot(file.exists(in_count_matrix))
stopifnot(file.exists(in_grid))
stopifnot(file.exists(in_study_area))

all_surveys  <- readRDS(in_surveys)
count_matrix <- readRDS(in_count_matrix)    
all_species  <- readRDS(in_species)      

grid_obj <- readRDS(in_grid)
grid_cells     <- grid_obj$grid_cells      # sf polygons (EPSG:3978)
grid_centroids <- grid_obj$grid_centroids  # sf points (EPSG:3978)

study_area <- readRDS(in_study_area)
boundary_buf <- study_area$boundary_buffer_25km %>% st_transform(crs_rast)

# Ensure grid is in raster CRS
grid_cells     <- st_transform(grid_cells, crs_rast)
grid_centroids <- st_transform(grid_centroids, crs_rast)

# Ensure surveys have obs_idx for later ordering
if (!("obs_idx" %in% names(all_surveys))) {
  all_surveys$obs_idx <- seq_len(nrow(all_surveys))
}

# Work in raster CRS for joins
all_surveys_m <- st_transform(all_surveys, crs_rast)

# ------------------------------------------------------------
# Load covariate rasters (processed outputs from 04_*)
# ------------------------------------------------------------

r_path <- function(x) file.path(cov_dir, x)

# Time-invariant
water_open  <- try_rast(r_path("water_open.tif"))
water_river <- try_rast(r_path("water_river.tif"))

# OBBA2
prec2  <- rast(r_path("prec_OBBA2.tif"))
tmax2  <- rast(r_path("tmax_OBBA2.tif"))
urban2 <- rast(r_path("Urban_2000.tif"))
lc2    <- rast(r_path("MCD12Q1_OBBA2_mode.tif"))
road2  <- try_rast(r_path("roadside_2005_100m.tif"))
insect2 <- try_rast(r_path("insect_OBBA2.tif"))

# OBBA3
prec3  <- rast(r_path("prec_OBBA3.tif"))
tmax3  <- rast(r_path("tmax_OBBA3.tif"))
urban3 <- rast(r_path("Urban_2020.tif"))
lc3    <- rast(r_path("MCD12Q1_OBBA3_mode.tif"))
road3  <- try_rast(r_path("roadside_2025_100m.tif"))
insect3 <- try_rast(r_path("insect_OBBA3.tif"))

# Optionally crop rasters to buffered boundary for speed
crop_to <- terra::vect(boundary_buf)

crop_if <- function(r) {
  if (is.null(r)) return(NULL)
  r <- terra::crop(r, crop_to)
  terra::mask(r, crop_to)
}

water_open  <- crop_if(water_open)
water_river <- crop_if(water_river)

prec2 <- crop_if(prec2); tmax2 <- crop_if(tmax2); urban2 <- crop_if(urban2); lc2 <- crop_if(lc2)
prec3 <- crop_if(prec3); tmax3 <- crop_if(tmax3); urban3 <- crop_if(urban3); lc3 <- crop_if(lc3)

road2 <- crop_if(road2); road3 <- crop_if(road3)
insect2 <- crop_if(insect2); insect3 <- crop_if(insect3)

# ------------------------------------------------------------
# Extract grid covariates (OBBA2 / OBBA3)
# ------------------------------------------------------------

grid_base <- grid_centroids %>%
  st_drop_geometry() %>%
  select(pixel_id)

grid_OBBA2 <- grid_centroids
grid_OBBA3 <- grid_centroids

# Continuous means
grid_OBBA2$prec <- exact_extract(prec2, grid_cells, "mean")
grid_OBBA2$tmax <- exact_extract(tmax2, grid_cells, "mean")
grid_OBBA3$prec <- exact_extract(prec3, grid_cells, "mean")
grid_OBBA3$tmax <- exact_extract(tmax3, grid_cells, "mean")

# Roads - fraction within footprint
if (!is.null(road2)) grid_OBBA2$road <- exact_extract(road2, grid_cells, "mean")
if (!is.null(road3)) grid_OBBA3$road <- exact_extract(road3, grid_cells, "mean")

# Categorical fractions: MODIS LC -> lc_*
lc_frac2 <- exact_extract(lc2, st_transform(grid_cells, st_crs(lc2)), "frac") %>% suppressWarnings()
names(lc_frac2) <- gsub("^frac", "lc", names(lc_frac2))
grid_OBBA2 <- bind_cols(grid_OBBA2, lc_frac2)

lc_frac3 <- exact_extract(lc3, st_transform(grid_cells, st_crs(lc3)), "frac") %>% suppressWarnings()
names(lc_frac3) <- gsub("^frac", "lc", names(lc_frac3))
grid_OBBA3 <- bind_cols(grid_OBBA3, lc_frac3)

# Categorical fractions: Urban -> urban_*
urb_frac2 <- exact_extract(urban2, st_transform(grid_cells, st_crs(urban2)), "frac") %>% suppressWarnings()
names(urb_frac2) <- gsub("^frac", "urban", names(urb_frac2))
grid_OBBA2 <- bind_cols(grid_OBBA2, urb_frac2)

urb_frac3 <- exact_extract(urban3, st_transform(grid_cells, st_crs(urban3)), "frac") %>% suppressWarnings()
names(urb_frac3) <- gsub("^frac", "urban", names(urb_frac3))
grid_OBBA3 <- bind_cols(grid_OBBA3, urb_frac3)

# Insects -> collapse to broadleaf/needleleaf then fractions
if (!is.null(insect2)) {
  insects2 <- collapse_insect_layers(insect2)
  frac_in2 <- exact_extract(insects2, st_transform(grid_cells, st_crs(insects2)), "frac") %>% suppressWarnings()
  frac_in2 <- drop_frac0_cols(frac_in2)
  grid_OBBA2 <- bind_cols(grid_OBBA2, clean_frac_names(frac_in2))
}

if (!is.null(insect3)) {
  insects3 <- collapse_insect_layers(insect3)
  frac_in3 <- exact_extract(insects3, st_transform(grid_cells, st_crs(insects3)), "frac") %>% suppressWarnings()
  frac_in3 <- drop_frac0_cols(frac_in3)
  grid_OBBA3 <- bind_cols(grid_OBBA3, clean_frac_names(frac_in3))
}

# Water –  interprets as fraction; is static through time so values are reused for atlas 2 and atlas 3
if (!is.null(water_open))  grid_OBBA2$water_open  <- exact_extract(water_open,  grid_cells, "mean")
if (!is.null(water_river)) grid_OBBA2$water_river <- exact_extract(water_river, grid_cells, "mean")
grid_OBBA3$water_open  <- grid_OBBA2$water_open
grid_OBBA3$water_river  <- grid_OBBA2$water_river

# Ensure no geometry duplication issues
grid_OBBA2 <- grid_OBBA2 %>% mutate(pixel_id = grid_centroids$pixel_id)
grid_OBBA3 <- grid_OBBA3 %>% mutate(pixel_id = grid_centroids$pixel_id)

# Ensure geometry is in final column position
grid_OBBA2 <- grid_OBBA2 %>% dplyr::relocate(geometry, .after = dplyr::last_col())
grid_OBBA3 <- grid_OBBA3 %>% dplyr::relocate(geometry, .after = dplyr::last_col())

# ------------------------------------------------------------
# Assign grid covariates to surveys (nearest grid centroid)
# ------------------------------------------------------------

if (!("Atlas" %in% names(all_surveys_m))) {
  stop("Survey object must contain an 'Atlas' column with values like 'OBBA2'/'OBBA3'.")
}

surveys2 <- all_surveys_m %>% filter(Atlas == "OBBA2")
surveys3 <- all_surveys_m %>% filter(Atlas == "OBBA3")

idx2 <- st_nearest_feature(surveys2, grid_OBBA2)
idx3 <- st_nearest_feature(surveys3, grid_OBBA3)

surveys2_covs <- surveys2 %>%
  bind_cols(st_drop_geometry(grid_OBBA2[idx2, ]))

surveys3_covs <- surveys3 %>%
  bind_cols(st_drop_geometry(grid_OBBA3[idx3, ]))

all_surveys_with_covs_m <- bind_rows(surveys2_covs, surveys3_covs) %>%
  arrange(obs_idx)


# ------------------------------------------------------------
# Assign each survey to an atlas square (typically 10 x 10 km)
# ------------------------------------------------------------

atlas_squares <- st_read(in_atlas_squares, quiet = TRUE) %>%
  st_make_valid() %>%
  st_transform(st_crs(all_surveys_with_covs_m)) %>%
  dplyr::select(square_id)

# index of the polygon each point falls in (NA if none)
hit <- st_within(all_surveys_with_covs_m, atlas_squares)

# stable: we assign in-place, no row reordering
all_surveys_with_covs_m$square_id <- dplyr::if_else(
  lengths(hit) == 0L,
  NA_character_,
  atlas_squares$square_id[vapply(hit, `[[`, integer(1), 1)]
)

# ------------------------------------------------------------
# Assign each survey and grid cell to a BCR within Ontario
# ------------------------------------------------------------

# Helper function
assign_poly_id <- function(pts, polys, id_col = "BCR", nearest_fallback = TRUE) {
  stopifnot(inherits(pts, "sf"), inherits(polys, "sf"))
  stopifnot(id_col %in% names(polys))
  
  # CRS match
  polys <- st_transform(polys, st_crs(pts))
  
  # within-index: list of polygon row indices per point
  hit <- st_within(pts, polys)
  
  # first hit (NA if none)
  idx <- ifelse(
    lengths(hit) == 0L,
    NA_integer_,
    vapply(hit, function(x) x[1], integer(1))
  )
  
  # optional nearest fallback for NAs
  if (nearest_fallback) {
    missing <- is.na(idx)
    if (any(missing)) {
      idx[missing] <- st_nearest_feature(pts[missing, ], polys)
    }
  }
  
  polys[[id_col]][idx]
}

# bcrs in Ontario
bcr_on <- st_read(in_bcr, quiet = TRUE) %>%
  st_make_valid() %>%
  dplyr::select(BCR) %>%
  dplyr::filter(BCR %in% c(7, 8, 12, 13)) %>%
  dplyr::group_by(BCR) %>%
  dplyr::summarise(geometry = st_union(geometry), .groups = "drop") %>%
  st_make_valid()

# Associate each empirical survey with a BCR
all_surveys_with_covs_m$BCR <- assign_poly_id(all_surveys_with_covs_m, bcr_on, "BCR", nearest_fallback = TRUE)

# Associate each prediction grid cell with a BCR
bcr_for_grid <- st_transform(bcr_on, st_crs(grid_centroids))
grid_bcr <- assign_poly_id(grid_centroids, bcr_for_grid, "BCR", nearest_fallback = TRUE)

grid_OBBA2$BCR <- grid_bcr
grid_OBBA3$BCR <- grid_bcr

# Sanity checks
table(is.na(all_surveys_with_covs_m$BCR))
table(is.na(grid_OBBA2$BCR))

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

# Return to original survey CRS (in km units)
all_surveys_with_covs <- st_transform(all_surveys_with_covs_m, st_crs(all_surveys)) %>% 
  dplyr::relocate(geometry, .after = dplyr::last_col())

grid_OBBA2 <- grid_OBBA2 %>% dplyr::relocate(geometry, .after = dplyr::last_col())
grid_OBBA3 <- grid_OBBA3 %>% dplyr::relocate(geometry, .after = dplyr::last_col())

saveRDS(
  list(
    all_surveys_with_covs = all_surveys_with_covs,
    count_matrix = count_matrix,
    grid_OBBA2 = grid_OBBA2,   # stays in EPSG:3978
    grid_OBBA3 = grid_OBBA3,   # stays in EPSG:3978
    grid_centroids = grid_centroids,
    grid_cells = grid_cells,
    all_species = all_species,
    boundary = study_area$boundary,
    boundary_buffer_25km = study_area$boundary_buffer_25km,
    date_created = Sys.time()
  ),
  file = out_file
)

message("06_extract_covariates.R complete.")
