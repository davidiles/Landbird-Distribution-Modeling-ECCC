# ============================================================
# 05b_refresh_survey_covariates.R
#
# Purpose
#   Survey-only refresh of analysis_data_covariates.rds. Re-extracts covariates
#   for the (updated) survey locations and patches all_surveys_with_covs +
#   count_matrix into the existing output, WITHOUT re-running the expensive
#   province-wide grid extraction. grid_OBBA2/3, grid geometry, species list,
#   and boundaries are carried over unchanged from the previous 05 run.
#
#   Run the full 05_extract_covariates.R whenever the grid, rasters, boundary,
#   or land-cover inputs change. Run THIS when only the survey inputs changed:
#   the carried-over grid covariates are only valid if their inputs are stale.
#
# Output (patched in place)
#   data_clean/birds/analysis_data_covariates.rds
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
# Config and inputs
# ------------------------------------------------------------

crs_rast    <- sf::st_crs(3978)   # extraction CRS (matches processed rasters)
BACKUP_PREV <- TRUE               # copy existing rds to a timestamped .bak first

birds_dir   <- file.path(paths$data_clean, "birds")
spatial_dir <- file.path(paths$data_clean, "spatial")

out_file         <- file.path(birds_dir, "analysis_data_covariates.rds")
in_surveys       <- file.path(paths$data_clean, "surveys", "surveys_raw.rds")
in_count_matrix  <- file.path(paths$data_clean, "surveys", "count_matrix_raw.rds")
in_study_area    <- file.path(spatial_dir, "study_area.rds")
in_atlas_squares <- file.path(paths$data, "Spatial", "National", "AtlasSquares",
                              "NationalSquares_FINAL.shp")

stopifnot(file.exists(out_file))  # this refresh patches a prior full-05 output

# ------------------------------------------------------------
# Load updated survey inputs
# ------------------------------------------------------------

# Original files
all_surveys  <- readRDS(in_surveys)
count_matrix <- readRDS(in_count_matrix)

# Re-assert survey <-> count row alignment BY KEY (see 05 for full rationale):
# obs_idx is a positional index that 06 uses to pull count rows, valid only if
# count_matrix shares all_surveys' order.
stopifnot(
  setequal(rownames(count_matrix), all_surveys$survey_id),
  anyDuplicated(all_surveys$survey_id) == 0
)


count_matrix <- count_matrix[all_surveys$survey_id, , drop = FALSE]
stopifnot(identical(rownames(count_matrix), all_surveys$survey_id))

study_area   <- readRDS(in_study_area)
boundary_buf <- study_area$boundary_buffer_5km %>%
  st_transform(crs_rast) %>%
  st_make_valid()

# Stable ordering field for surveys.
if (!("obs_idx" %in% names(all_surveys))) {
  all_surveys$obs_idx <- seq_len(nrow(all_surveys))
}
all_surveys   <- all_surveys %>% relocate(obs_idx)
all_surveys_m <- st_transform(all_surveys, crs_rast)

# 1-km square footprints centred on each survey (for area-based extraction).
survey_pixels <- make_square_buffer(all_surveys_m, half_width = 500)
survey_pixels <- st_set_geometry(all_surveys_m, survey_pixels)

# ------------------------------------------------------------
# Load covariate rasters (lazy; only survey footprints are extracted here)
# ------------------------------------------------------------

roads_2005   <- try_rast(file.path(spatial_dir, "roads_2005_buf_125m.tif"))
roads_2025   <- try_rast(file.path(spatial_dir, "roads_2025_buf_125m.tif"))
rivers_large <- try_rast(file.path(spatial_dir, "rivers_large_buf_500m.tif"))
rivers_small <- try_rast(file.path(spatial_dir, "rivers_small_buf_250m.tif"))
lakes_large  <- try_rast(file.path(spatial_dir, "lakes_large_buf_500m.tif"))
lakes_small  <- try_rast(file.path(spatial_dir, "lakes_small_buf_250m.tif"))
coastline    <- try_rast(file.path(spatial_dir, "coastline_buf_3km.tif"))
great_lakes  <- try_rast(file.path(spatial_dir, "great_lakes_buf_3km.tif"))
open_water   <- try_rast(file.path(spatial_dir, "waterbodies.tif"))

LCC2020 <- try_rast(file.path(paths$data, "Spatial", "National",
                              "LandCoverCanada2020", "landcover-2020-classification.tif"))
LCC2010 <- try_rast(file.path(paths$data, "Spatial", "National",
                              "LandCoverCanada2010", "landcover-2010-classification.tif"))

# ------------------------------------------------------------
# Extract survey covariates (OBBA2 = 2010 LCC, OBBA3 = 2020 LCC)
#
# Same functions / rasters / footprints / boundary_buf as 05, so results are
# identical to what a full 05 run would produce for the survey side.
# ------------------------------------------------------------

extract_survey_covs <- function(surveys, lcc, roads) {
  surveys <- bind_cols(
    surveys,
    extract_frac(lcc, surveys, boundary_buf, prefix = "LCC", drop0 = TRUE, clean_names = FALSE)
  )
  surveys$road         <- extract_mean(roads, surveys, boundary_buf)
  surveys$rivers_large <- extract_mean(rivers_large, surveys, boundary_buf)
  surveys$rivers_small <- extract_mean(rivers_small, surveys, boundary_buf)
  surveys$lakes_large  <- extract_mean(lakes_large, surveys, boundary_buf)
  surveys$lakes_small  <- extract_mean(lakes_small, surveys, boundary_buf)
  surveys$great_lakes  <- extract_mean(great_lakes, surveys, boundary_buf)
  surveys$coastline    <- extract_mean(coastline, surveys, boundary_buf)
  surveys$open_water   <- extract_mean(open_water, surveys, boundary_buf)
  surveys
}

surveys2 <- extract_survey_covs(surveys = survey_pixels %>% filter(Atlas == "OBBA2"),
                                lcc = LCC2010, roads = roads_2005)
surveys3 <- extract_survey_covs(surveys = survey_pixels %>% filter(Atlas == "OBBA3"),
                                lcc = LCC2020, roads = roads_2025)

# Recombine into a single point sf, ordered by obs_idx.
all_surveys_with_covs <- bind_rows(surveys2, surveys3) %>%
  arrange(obs_idx) %>%
  st_centroid()

# ------------------------------------------------------------
# Assign each survey to an atlas square
# ------------------------------------------------------------

if (file.exists(in_atlas_squares)) {
  atlas_squares <- st_read(in_atlas_squares, quiet = TRUE) %>%
    st_make_valid() %>%
    st_transform(st_crs(all_surveys_with_covs)) %>%
    select(square_id)

  hit       <- st_within(all_surveys_with_covs, atlas_squares)
  square_id <- rep(NA_character_, nrow(all_surveys_with_covs))
  has_hit   <- lengths(hit) > 0L
  if (any(has_hit)) {
    idx <- vapply(hit[has_hit], function(x) x[1], integer(1))
    square_id[has_hit] <- atlas_squares$square_id[idx]
  }
  all_surveys_with_covs$square_id <- square_id
} else {
  message("Atlas squares file not found; skipping square_id assignment.")
  all_surveys_with_covs$square_id <- NA_character_
}

# Surveys returned to km-unit CRS for downstream scripts.
all_surveys_with_covs <- st_transform(all_surveys_with_covs, st_crs(all_surveys)) %>%
  relocate(geometry, .after = last_col())

# ------------------------------------------------------------
# Patch the existing output in place
# (grid_OBBA2/3, grid geometry, species, boundaries carried over unchanged)
# ------------------------------------------------------------

prev <- readRDS(out_file)

# Guard against a stale carried-over species list after adding surveys.
new_species <- setdiff(colnames(count_matrix), colnames(prev$count_matrix))
if (length(new_species)) {
  warning("count_matrix gained ", length(new_species),
          " species column(s) absent from the previous run; carried-over ",
          "all_species may be stale: ", paste(head(new_species, 10), collapse = ", "))
}

if (BACKUP_PREV) {
  bak <- sub("\\.rds$", paste0("_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".bak.rds"), out_file)
  file.copy(out_file, bak, overwrite = FALSE)
  message("Backed up previous output -> ", basename(bak))
}

prev$all_surveys_with_covs <- all_surveys_with_covs
prev$count_matrix          <- count_matrix
prev$date_created          <- Sys.time()
prev$refreshed_surveys_only <- Sys.time()  # record that grids were NOT re-extracted

saveRDS(prev, file = out_file)

message("05b_refresh_survey_covariates.R complete (survey-only refresh).")
