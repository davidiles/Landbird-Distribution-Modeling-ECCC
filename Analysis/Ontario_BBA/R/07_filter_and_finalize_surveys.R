# ============================================================
# 07_filter_and_finalize_surveys.R
#
# Purpose
#   Create the final analysis-ready survey dataset used by downstream
#   modeling scripts.
#
# What this script does
#   1) Loads survey, count, grid, and boundary objects produced by
#      06_extract_covariates.R.
#   2) Filters surveys to a common analysis window:
#        - 5-minute duration
#        - Hours_Since_Sunrise in [-1, 6]
#        - DayOfYear in [135, 196]
#        - within the study boundary
#   3) Removes duplicate survey locations, retaining one survey per
#      rounded x/y location.
#   4) Subsets the survey × species count matrix to match the filtered
#      surveys, preserving alignment through obs_idx.
#   5) Adds survey-level helper variables used in modeling.
#   6) Flags prediction-grid points that intersect optional large-water
#      polygons, if available.
#   7) Standardizes selected continuous covariates using statistics
#      calculated from the filtered surveys.
#   8) Adds derived north/south BCR covariates and simple road/river
#      indicators.
#   9) Builds species summary tables:
#        - detections and occupied squares by atlas/BCR
#        - raw-count dispersion summaries for Poisson vs nbinomial triage
#        - safe-date summaries by atlas/BCR
#   10) Builds a hex grid for later regional summaries.
#
# Inputs
#   - data_clean/birds/analysis_data_covariates.rds
#       must contain:
#         * all_surveys_with_covs (sf)
#         * count_matrix
#         * grid_OBBA2, grid_OBBA3 (sf)
#         * boundary (sf)
#         * bcr_sf (sf)
#         * all_species
#   - R/functions/survey_processing_utils.R
#       used functions:
#         * as_counts_tbl()
#         * dedupe_by_location()
#         * standardize_covars()
#   - R/functions/inla_model_utils.R
#       used functions:
#         * make_hex_grid()
#   - Optional:
#       data_clean/spatial/water_filtered.shp
#       data_clean/metadata/safe_dates_OBBA.xlsx
#
# Output
#   - data_clean/birds/data_ready_for_analysis.rds
#       list with:
#         * study_boundary
#         * bcr_sf
#         * all_surveys
#         * counts
#         * grid_OBBA2, grid_OBBA3
#         * species_to_model
#         * outside_safe_dates
#         * missing_safe_dates
#         * dispersion_screening
#         * hex_grid
#         * safe_dates
#         * date_created
#
# Notes
#   - obs_idx is treated as the row-position key linking filtered surveys
#     to the saved count matrix.
#   - dedupe_by_location() uses random sampling; a fixed seed is set here
#     for reproducibility.
# ============================================================

rm(list=ls())

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(here)
  library(readxl)
})

# ------------------------------------------------------------
# Paths and utilities
# ------------------------------------------------------------

source(here::here("R", "00_config_paths.R"))
source(file.path(paths$functions, "survey_processing_utils.R"))
source(file.path(paths$functions, "inla_model_utils.R"))

dir.create(file.path(paths$data_clean, "birds"), recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Config
# ------------------------------------------------------------

in_file  <- file.path(paths$data_clean, "birds", "analysis_data_covariates.rds")
out_file <- file.path(paths$data_clean, "birds", "data_ready_for_analysis.rds")

keep_duration_min <- 5
hss_min <- -1
hss_max <- 6
doy_min <- 135
doy_max <- 196

dup_coord_digits <- 2
dedupe_seed <- 1

covars_for_stats <- c(
  "insect_broadleaf", "insect_needleleaf", "road", "water_river", "prec", "tmax",
  "urban_2", "urban_3",
  "lc_1", "lc_4", "lc_5", "lc_8", "lc_9", "lc_10", "lc_11", "lc_12", "lc_14", "lc_17"
)
covars_to_standardize <- c("prec", "tmax")

south_bcr <- c(12, 13)
north_bcr <- c(7, 8)

water_filtered_path <- file.path(paths$data_clean, "spatial", "water_filtered.shp")
safe_date_path <- file.path(paths$data_clean, "metadata", "safe_dates_OBBA.xlsx")

# ------------------------------------------------------------
# Load input object
# ------------------------------------------------------------

dat <- readRDS(in_file)

all_surveys    <- dat$all_surveys_with_covs
count_matrix   <- dat$count_matrix
grid_OBBA2     <- dat$grid_OBBA2
grid_OBBA3     <- dat$grid_OBBA3
study_boundary <- dat$boundary

stopifnot(inherits(all_surveys, "sf"))
stopifnot(inherits(grid_OBBA2, "sf"))
stopifnot(inherits(grid_OBBA3, "sf"))

water_filtered <- if (file.exists(water_filtered_path)) {
  st_read(water_filtered_path, quiet = TRUE)
} else {
  NULL
}

all_species_unique <- dat$all_species %>%
  distinct(species_id, .keep_all = TRUE) %>%
  transmute(
    species_id = as.character(species_id),
    english_name
  )

# ------------------------------------------------------------
# Basic data checks and setup
# ------------------------------------------------------------

if (!("obs_idx" %in% names(all_surveys))) {
  all_surveys <- all_surveys %>%
    mutate(obs_idx = row_number(), .before = 1)
} else {
  all_surveys <- all_surveys %>%
    arrange(obs_idx)
}

counts <- as_counts_tbl(count_matrix)

if (!("DayOfYear" %in% names(all_surveys))) {
  all_surveys <- all_surveys %>%
    mutate(DayOfYear = yday(Date_Time))
}

required_cols <- c("Survey_Duration_Minutes", "Hours_Since_Sunrise", "DayOfYear", "Date_Time", "Atlas")
missing_cols <- setdiff(required_cols, names(all_surveys))
if (length(missing_cols) > 0) {
  stop("Missing required survey columns: ", paste(missing_cols, collapse = ", "))
}

# ------------------------------------------------------------
# 1) Filter surveys to the analysis window
# ------------------------------------------------------------

surveys_f <- all_surveys %>%
  filter(
    round(Survey_Duration_Minutes) == keep_duration_min,
    Hours_Since_Sunrise >= hss_min,
    Hours_Since_Sunrise <= hss_max,
    DayOfYear >= doy_min,
    DayOfYear <= doy_max
  ) %>%
  arrange(obs_idx)

# ------------------------------------------------------------
# 2) Keep only surveys inside the study boundary
# ------------------------------------------------------------

study_boundary <- st_transform(study_boundary, st_crs(surveys_f))
inside <- st_within(surveys_f, study_boundary, sparse = FALSE)[, 1]

surveys_f <- surveys_f[inside, , drop = FALSE] %>%
  arrange(obs_idx)

# ------------------------------------------------------------
# 3) Remove duplicate survey locations
# ------------------------------------------------------------

set.seed(dedupe_seed)
surveys_f <- dedupe_by_location(surveys_f, digits = dup_coord_digits) %>%
  arrange(obs_idx)

# ------------------------------------------------------------
# 4) Subset counts to the retained surveys
# ------------------------------------------------------------

stopifnot(max(surveys_f$obs_idx) <= nrow(counts))

counts_f <- counts %>%
  slice(surveys_f$obs_idx) %>%
  mutate(obs_idx = surveys_f$obs_idx)

stopifnot(nrow(counts_f) == nrow(surveys_f))

# ------------------------------------------------------------
# 5) Add helper columns used downstream
# ------------------------------------------------------------

surveys_f <- surveys_f %>%
  mutate(
    Atlas3   = if_else(Atlas == "OBBA3", 1L, 0L),
    Atlas3_c = Atlas3 - 0.5,
    days_rescaled = DayOfYear - 166,
    BCR_idx = as.integer(factor(BCR))
  )

if ("Survey_Type" %in% names(surveys_f)) {
  surveys_f <- surveys_f %>%
    mutate(ARU = if_else(Survey_Type == "ARU", 1L, 0L))
} else {
  surveys_f <- surveys_f %>%
    mutate(ARU = 0L)
}

if ("square_id" %in% names(surveys_f)) {
  surveys_f <- surveys_f %>%
    mutate(
      square_atlas = as.integer(factor(paste0(square_id, "-", Atlas))),
      square_year  = as.integer(factor(paste0(square_id, "-", year(Date_Time))))
    )
} else {
  surveys_f <- surveys_f %>%
    mutate(
      square_atlas = NA_integer_,
      square_year  = NA_integer_
    )
}

# ------------------------------------------------------------
# 6) Prediction grids: add on_water and pixel_id
# ------------------------------------------------------------

if (!is.null(water_filtered)) {
  water_filtered <- water_filtered %>%
    st_make_valid() %>%
    st_transform(st_crs(grid_OBBA2))
  
  grid_OBBA2$on_water <- lengths(st_intersects(grid_OBBA2, water_filtered)) > 0
  grid_OBBA3$on_water <- grid_OBBA2$on_water
} else {
  grid_OBBA2$on_water <- NA
  grid_OBBA3$on_water <- NA
}

core_grid_vars <- intersect(c("prec", "tmax"), names(grid_OBBA2))
if (length(core_grid_vars) > 0) {
  grid_OBBA2 <- grid_OBBA2 %>%
    filter(if_all(all_of(core_grid_vars), ~ !is.na(.x)))
  
  grid_OBBA3 <- grid_OBBA3 %>%
    filter(if_all(all_of(core_grid_vars), ~ !is.na(.x)))
}

if (!("pixel_id" %in% names(grid_OBBA2))) {
  grid_OBBA2$pixel_id <- seq_len(nrow(grid_OBBA2))
}
if (!("pixel_id" %in% names(grid_OBBA3))) {
  grid_OBBA3$pixel_id <- seq_len(nrow(grid_OBBA3))
}

# ------------------------------------------------------------
# 7) Standardize selected continuous covariates
# ------------------------------------------------------------

std <- standardize_covars(
  surveys_sf = surveys_f,
  grid2 = grid_OBBA2,
  grid3 = grid_OBBA3,
  covars_for_stats = covars_for_stats,
  covars_to_standardize = covars_to_standardize
)

surveys_f   <- std$surveys
grid_OBBA2  <- std$grid2
grid_OBBA3  <- std$grid3
covar_stats <- std$stats

# ------------------------------------------------------------
# 8) Add derived north/south BCR covariates
# ------------------------------------------------------------

add_derived_covariates <- function(df, south_bcr = c(12, 13), north_bcr = c(7, 8)) {
  if (!("BCR" %in% names(df))) {
    stop("BCR missing; required for derived covariates.")
  }
  
  get_or_zero <- function(nm) {
    if (nm %in% names(df)) df[[nm]] else rep(0, nrow(df))
  }
  
  df <- df %>%
    mutate(
      on_river = as.integer(get_or_zero("water_river") > 0),
      on_road  = as.integer(get_or_zero("road") > 0)
    )
  
  for (nm in c("lc_1", "lc_4", "lc_5", "lc_8", "lc_9", "lc_10")) {
    df[[paste0(nm, "S")]] <- if_else(df$BCR %in% south_bcr, get_or_zero(nm), 0)
    df[[paste0(nm, "N")]] <- if_else(df$BCR %in% north_bcr, get_or_zero(nm), 0)
  }
  
  df %>%
    mutate(
      onriver_S = if_else(BCR %in% south_bcr, on_river, 0L),
      onriver_N = if_else(BCR %in% north_bcr, on_river, 0L)
    )
}

surveys_f  <- add_derived_covariates(surveys_f, south_bcr, north_bcr)
grid_OBBA2 <- add_derived_covariates(grid_OBBA2, south_bcr, north_bcr)
grid_OBBA3 <- add_derived_covariates(grid_OBBA3, south_bcr, north_bcr)

# ------------------------------------------------------------
# 9) Build survey metadata and long-format count table once
# ------------------------------------------------------------

survey_meta <- surveys_f %>%
  st_drop_geometry() %>%
  mutate(
    ecoregion = case_when(
      BCR == 7 ~ "Hudson Plains",
      BCR %in% c(8, 12) ~ "Boreal Shield",
      BCR == 13 ~ "Mixedwood Plains",
      TRUE ~ NA_character_
    ),
    survey_doy = yday(as.Date(Date_Time))
  )

counts_long <- counts_f %>%
  pivot_longer(
    cols = -obs_idx,
    names_to = "species_id",
    values_to = "count"
  ) %>%
  mutate(species_id = as.character(species_id)) %>%
  left_join(survey_meta, by = "obs_idx") %>%
  left_join(all_species_unique, by = "species_id")

# ------------------------------------------------------------
# 10) Build hex grid for later regional summaries
# ------------------------------------------------------------

hex_grid <- make_hex_grid(study_boundary, width_km = 25)

# ------------------------------------------------------------
# 11) Load and apply safe dates
# ------------------------------------------------------------

species_manual_safedates <- c("Canada Jay", "Common Raven", "White-winged Crossbill", "Red Crossbill")

if (!file.exists(safe_date_path)) {
  stop("Cannot find safe date file at: ", safe_date_path)
}

safe_dates <- read_xlsx(safe_date_path) %>%
  mutate(
    start_doy = yday(as.Date(start)) - 4,
    end_doy   = yday(as.Date(end)) + 4
  ) %>%
  select(sp_english, ecoregion, start_doy, end_doy)

median_start <- median(safe_dates$start_doy, na.rm = TRUE)
median_end   <- median(safe_dates$end_doy, na.rm = TRUE)

safe_dates <- safe_dates %>%
  mutate(
    start_doy = if_else(sp_english %in% species_manual_safedates, median_start, start_doy),
    end_doy   = if_else(sp_english %in% species_manual_safedates, median_end, end_doy)
  )

counts_long <- counts_long %>%
  left_join(
    safe_dates,
    by = c("english_name" = "sp_english", "ecoregion" = "ecoregion")
  ) %>%
  mutate(
    safe_date_status = case_when(
      is.na(start_doy) | is.na(end_doy) ~ "no_safe_dates",
      survey_doy >= start_doy & survey_doy <= end_doy ~ "inside",
      TRUE ~ "outside"
    )
  )

species_bcr_summary <- counts_long %>%
  group_by(species_id, english_name, Atlas, BCR, ecoregion) %>%
  summarise(
    detections_total = sum(count > 0, na.rm = TRUE),
    detections_safe = sum(count > 0 & safe_date_status == "inside", na.rm = TRUE),
    detections_outside = sum(count > 0 & safe_date_status == "outside", na.rm = TRUE),
    detections_no_safe = sum(count > 0 & safe_date_status == "no_safe_dates", na.rm = TRUE),
    prop_outside = if_else(
      detections_total > 0,
      detections_outside / detections_total,
      NA_real_
    ),
    prop_no_safe = if_else(
      detections_total > 0,
      detections_no_safe / detections_total,
      NA_real_
    ),
    .groups = "drop"
  )

outside_safe_dates <- species_bcr_summary %>%
  filter(prop_outside >= 0.1, detections_total >= 100) %>%
  arrange(english_name, BCR, Atlas)

missing_safe_dates <- species_bcr_summary %>%
  filter(detections_total > 10, prop_no_safe > 0) %>%
  arrange(desc(detections_total))

species_safecounts_per_atlas <- counts_long %>%
  group_by(species_id, english_name, Atlas) %>%
  summarise(
    detections_safe = sum(count > 0 & safe_date_status == "inside", na.rm = TRUE),
    n_squares_safe = n_distinct(square_id[count > 0 & safe_date_status == "inside"]),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Atlas,
    values_from = c(detections_safe, n_squares_safe),
    names_glue = "{.value}_{Atlas}",
    values_fill = 0
  )

species_to_model <- species_safecounts_per_atlas %>%
  arrange(desc(n_squares_safe_OBBA3)) %>%
  relocate(species_id, english_name)


# ------------------------------------------------------------
# Save covariate rasters
# ------------------------------------------------------------
rast_dir  <- file.path(paths$data_clean,"spatial")

rast_path_a2 <- file.path(rast_dir, "atlas2_cov_rasterstack.tif")
rast_path_a3 <- file.path(rast_dir, "atlas3_cov_rasterstack.tif")

# contains function for rasterizing
figure_utils_path <- file.path(paths$functions, "figure_utils.R")
source(figure_utils_path)

vars_to_rasterize <- c("prec",
                       "tmax",
                       "on_road",
                       "on_river",
                       "lc_1",
                       "lc_4",
                       "lc_5",
                       "lc_8",
                       "lc_9",
                       "lc_10",
                       "lc_11",
                       "lc_12",
                       "lc_14",
                       "lc_17",
                       "insect_broadleaf",
                       "insect_needleleaf")

# Covariate raster stack for atlas 2
r2 = rasterize_sf(grid_OBBA2,
                  vars_to_rasterize,
                  res = 1001,
                  metadata = c(
                    description = c("Covariates at 1 km resolution")
                  ))
writeRaster(r2,filename = rast_path_a2,overwrite = TRUE)


# Covariate raster stack for atlas 2
r3 = rasterize_sf(grid_OBBA3,
                  vars_to_rasterize,
                  res = 1001,
                  metadata = c(
                    description = c("Covariates at 1 km resolution")
                  ))
writeRaster(r3,filename = rast_path_a3,overwrite = TRUE)

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

saveRDS(
  list(
    
    study_boundary = study_boundary,
    bcr_sf = dat$bcr_sf,
    
    all_surveys = surveys_f,
    counts = counts_f,
    grid_OBBA2 = grid_OBBA2,
    grid_OBBA3 = grid_OBBA3,
    
    species_to_model = species_to_model,
    species_bcr_summary = species_bcr_summary,
    
    outside_safe_dates = outside_safe_dates,
    missing_safe_dates = missing_safe_dates,
    
    hex_grid = hex_grid,
    safe_dates = safe_dates,
    date_created = Sys.time()
  ),
  file = out_file
)

message("07_filter_and_finalize_surveys.R complete")