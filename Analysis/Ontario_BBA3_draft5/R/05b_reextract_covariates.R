# ============================================================
# 05b_reextract_covariates.R
#
# Re-extract a subset of covariates whose source rasters changed and splice
# the new values into the existing covariate file in place. Everything else
# (LCC fractions, roads, square_id, boundary crop) is left untouched.
# ============================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(terra)
  library(exactextractr)
  library(here)
})

source(here::here("R", "00_config_paths.R"))
source(file.path(paths$functions, "covariate_processing_utils.R"))

# ------------------------------------------------------------
# Config
# ------------------------------------------------------------

crs_rast       <- sf::st_crs(3978)                 # extraction CRS (matches rasters)
birds_dir      <- file.path(paths$data_clean, "birds")
spatial_dir    <- file.path(paths$data_clean, "spatial")
covariate_file <- file.path(birds_dir, "analysis_data_covariates.rds")
covariate_file_new <- file.path(birds_dir, "analysis_data_covariates_new.rds")

# Covariates to re-extract -> their (updated) rasters. Edit this list whenever a
# different set of rasters changes. All are extract_mean() layers and
# time-invariant (same raster for both atlases), so each is extracted once and
# written to OBBA2 + OBBA3 + surveys.
reextract_rasters <- list(
  rivers_large = file.path(spatial_dir, "rivers_large_buf_500m.tif"),
  rivers_small = file.path(spatial_dir, "rivers_small_buf_250m.tif"),
  lakes_large  = file.path(spatial_dir, "lakes_large_buf_500m.tif"),
  lakes_small  = file.path(spatial_dir, "lakes_small_buf_250m.tif"),
  coastline    = file.path(spatial_dir, "coastline_buf_3km.tif"),
  great_lakes  = file.path(spatial_dir, "great_lakes_buf_3km.tif")
)

stopifnot(all(file.exists(unlist(reextract_rasters))))

# ------------------------------------------------------------
# Load the previously prepared covariate file
# ------------------------------------------------------------

dat <- readRDS(covariate_file)

# Columns being replaced must already exist (guards against a typo silently
# appending a new column instead of overwriting).
stopifnot(
  all(names(reextract_rasters) %in% names(dat$grid_OBBA2)),
  all(names(reextract_rasters) %in% names(dat$grid_OBBA3)),
  all(names(reextract_rasters) %in% names(dat$all_surveys_with_covs))
)

boundary_buf <- dat$boundary_buffer_5km %>%
  st_transform(crs_rast) %>%
  st_make_valid()

# ------------------------------------------------------------
# Rebuild grid extraction geometry (retained cells only)
#
# grid_OBBA2/OBBA3 hold boundary-cropped *centroids*; area-weighted extraction
# needs the *cell polygons*. grid_cells and grid_centroids are the full,
# same-order grid, so recover each retained row's cell by a geometry match.
# ------------------------------------------------------------

grid_cells     <- dat$grid_cells
grid_centroids <- dat$grid_centroids
stopifnot(nrow(grid_cells) == nrow(grid_centroids))

match2 <- st_equals(dat$grid_OBBA2, grid_centroids)
match3 <- st_equals(dat$grid_OBBA3, grid_centroids)
stopifnot(all(lengths(match2) == 1L), all(lengths(match3) == 1L))
cell_idx  <- vapply(match2, `[`, integer(1), 1L)
cell_idx3 <- vapply(match3, `[`, integer(1), 1L)
stopifnot(identical(cell_idx, cell_idx3))   # both grids: same cells, same order

cells_keep <- grid_cells[cell_idx, ]
stopifnot(nrow(cells_keep) == nrow(dat$grid_OBBA2))

# ------------------------------------------------------------
# Rebuild survey extraction geometry (1-km square footprints)
# ------------------------------------------------------------

surveys_m  <- st_transform(dat$all_surveys_with_covs, crs_rast)
survey_pix <- make_square_buffer(surveys_m, half_width = 500)
survey_pix <- st_set_geometry(surveys_m, survey_pix)
stopifnot(nrow(survey_pix) == nrow(dat$all_surveys_with_covs))

# ------------------------------------------------------------
# Re-extract and splice in (each raster read once)
# ------------------------------------------------------------

for (nm in names(reextract_rasters)) {
  r <- try_rast(reextract_rasters[[nm]])
  
  grid_vals <- extract_mean(r, cells_keep, boundary_buf)
  dat$grid_OBBA2[[nm]] <- grid_vals
  dat$grid_OBBA3[[nm]] <- grid_vals   # hydro is time-invariant: same for both
  
  dat$all_surveys_with_covs[[nm]] <- extract_mean(r, survey_pix, boundary_buf)
}

# ------------------------------------------------------------
# Record the re-extraction and overwrite the file in place
# ------------------------------------------------------------

dat$reextraction_log <- c(
  dat$reextraction_log,
  list(list(
    covariates = names(reextract_rasters),
    rasters    = unname(unlist(reextract_rasters)),
    date       = Sys.time()
  ))
)

saveRDS(dat, covariate_file_new)
message("05b_reextract_covariates.R complete: re-extracted ",
        paste(names(reextract_rasters), collapse = ", "))