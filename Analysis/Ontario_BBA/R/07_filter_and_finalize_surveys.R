# ============================================================
# 07_filter_and_finalize_surveys.R
#
# Purpose:
#   Finalize the analysis-ready dataset by applying all survey filters
#   and final bookkeeping in ONE place.
#
#   This script:
#     1) Loads data_clean/birds/analysis_data_covariates.rds (from script 06)
#     2) Filters surveys to analysis window:
#          - 5-minute surveys
#          - Hours since sunrise in [-1, 6]
#          - Day-of-year in [135, 196]
#          - inside study boundary
#     3) Removes duplicate locations (keep one survey per rounded x/y)
#     4) Subsets count matrix to remaining surveys
#     5) Adds helper columns: Atlas3, ARU, square_atlas, square_year
#     6) Flags grid points overlapping large water polygons (optional)
#     7) Standardizes prec and tmax only (using survey means/sds)
#     8) Ensures grid has no NA rows; assigns pixel_id
#     9) Builds species summary tables for modeling triage
#
# Inputs:
#   data_clean/birds/analysis_data_covariates.rds
#
# Outputs:
#   data_clean/birds/data_ready_for_analysis.rds
#
# Notes:
#   - This script assumes surveys already have:
#       * Survey_Duration_Minutes
#       * Hours_Since_Sunrise
#       * Date_Time
#       * Atlas (OBBA2/OBBA3)
#       * square_id
#       * BCR
#   - Grids are expected to be in EPSG:3978 (meters).
# ============================================================

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(tidyr)
  library(lubridate)
})

# helper-functions
source("R/functions/survey_processing_utils.R")  

# ------------------------------------------------------------
# Config
# ------------------------------------------------------------

in_file  <- "data_clean/birds/analysis_data_covariates.rds"
out_file <- "data_clean/birds/data_ready_for_analysis.rds"

# Filtering windows
keep_duration_min <- 5
hss_min <- -1
hss_max <- 6
doy_min <- 135
doy_max <- 196

# Duplicate location handling (rounding in projected coordinates)
dup_coord_digits <- 5

# Covariates used for computing means/sds (some may be absent; that's OK)
covars_for_stats <- c(
  "insect_broadleaf","insect_needleleaf","road","water_river","prec","tmax",
  "urban_2","urban_3",
  "lc_1","lc_4","lc_5","lc_8","lc_9","lc_10","lc_11","lc_12","lc_14","lc_17"
)
covars_to_standardize <- c("prec","tmax")

# ------------------------------------------------------------
# Load input object
# ------------------------------------------------------------

stopifnot(file.exists(in_file))
dat <- readRDS(in_file)

all_surveys <- dat$all_surveys_with_covs
count_matrix <- dat$count_matrix
grid_OBBA2 <- dat$grid_OBBA2
grid_OBBA3 <- dat$grid_OBBA3
study_boundary <- dat$boundary
water_filtered <- sf::st_read("data_clean/spatial/water_filtered.shp")

stopifnot(inherits(all_surveys, "sf"))
stopifnot(inherits(grid_OBBA2, "sf"))
stopifnot(inherits(grid_OBBA3, "sf"))

# Ensure obs_idx exists and is stable for subsetting
if (!("obs_idx" %in% names(all_surveys))) {
  all_surveys <- all_surveys %>% mutate(obs_idx = row_number(), .before = 1)
} else {
  all_surveys <- all_surveys %>% arrange(obs_idx)
}

counts <- as_counts_tbl(count_matrix)

# Day-of-year (if missing)
if (!("DayOfYear" %in% names(all_surveys))) {
  all_surveys <- all_surveys %>% mutate(DayOfYear = yday(Date_Time))
}

# ------------------------------------------------------------
# 1) Filter surveys (duration, sunrise window, season)
# ------------------------------------------------------------

req <- c("Survey_Duration_Minutes", "Hours_Since_Sunrise", "DayOfYear", "Date_Time")
missing <- setdiff(req, names(all_surveys))
if (length(missing) > 0) stop("Missing required survey columns: ", paste(missing, collapse = ", "))

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
# 2) Keep only surveys inside study boundary
# ------------------------------------------------------------

# Ensure CRS compatibility
study_boundary <- st_transform(study_boundary, st_crs(surveys_f))

inside <- st_within(surveys_f, study_boundary, sparse = FALSE)[, 1]
surveys_f <- surveys_f[inside, ] %>% arrange(obs_idx)

# ------------------------------------------------------------
# 3) Retain a single survey per repeated location
# ------------------------------------------------------------

surveys_f <- dedupe_by_location(surveys_f, digits = dup_coord_digits) %>%
  arrange(obs_idx)

# ------------------------------------------------------------
# 4) Subset counts matrix to remaining surveys
# ------------------------------------------------------------

if (!all(surveys_f$obs_idx %in% counts$obs_idx)) {
  stop("Some filtered surveys have obs_idx not present in counts. Check that count_matrix row order matches all_surveys.")
}

counts_f <- counts[surveys_f$obs_idx,]

# ------------------------------------------------------------
# 5) Add helper columns (Atlas3, ARU, square_* ids)
# ------------------------------------------------------------

if (!("Atlas" %in% names(surveys_f))) stop("Missing 'Atlas' column in surveys (expected OBBA2/OBBA3).")

surveys_f <- surveys_f %>%
  mutate(
    Atlas3 = ifelse(Atlas == "OBBA2", 0, 1),
    ARU = ifelse("Survey_Type" %in% names(surveys_f) & Survey_Type == "ARU", 1, 0)
  )

if ("square_id" %in% names(surveys_f)) {
  surveys_f <- surveys_f %>%
    mutate(
      square_atlas = as.numeric(factor(paste0(square_id, "-", Atlas))),
      square_year  = as.numeric(factor(paste0(square_id, "-", year(Date_Time))))
    )
}

# ------------------------------------------------------------
# 6) Grid: on_water flag (optional) + NA removal + pixel_id
# ------------------------------------------------------------

# Use water polygon layer if available, for masking out predictions
if (!is.null(water_filtered)) {
  water_filtered <- st_transform(water_filtered, st_crs(grid_OBBA2))
  overlaps2 <- st_intersects(st_centroid(grid_OBBA2), water_filtered)
  overlaps3 <- st_intersects(st_centroid(grid_OBBA3), water_filtered)
  grid_OBBA2$on_water <- lengths(overlaps2) > 0
  grid_OBBA3$on_water <- lengths(overlaps3) > 0
}

grid_OBBA2 <- na.omit(grid_OBBA2)
grid_OBBA3 <- na.omit(grid_OBBA3)

# Keep a consistent pixel_id across periods (assumes same row order/geometry)
grid_OBBA2$pixel_id <- seq_len(nrow(grid_OBBA2))
grid_OBBA3$pixel_id <- seq_len(nrow(grid_OBBA3))

# ------------------------------------------------------------
# 7) Standardize prec and tmax only (based on surveys)
# ------------------------------------------------------------

std <- standardize_covars(
  surveys_sf = surveys_f,
  grid2 = grid_OBBA2,
  grid3 = grid_OBBA3,
  covars_for_stats = covars_for_stats,
  covars_to_standardize = covars_to_standardize
)

surveys_f <- std$surveys
grid_OBBA2 <- std$grid2
grid_OBBA3 <- std$grid3
covar_stats <- std$stats

# ------------------------------------------------------------
# 8) Create several "derived" covariates, to separate their effects into
#    north and south components
# ------------------------------------------------------------

# Create “derived” covariates consistently across surveys + grids
add_derived_covariates <- function(surveys, grid2, grid3,
                                   south_bcr = c(12, 13), north_bcr = c(7, 8)) {
  
  make_bits <- function(df) {
    if (!("water_river" %in% names(df))) df$water_river <- NA_real_
    if (!("BCR" %in% names(df))) stop("BCR missing; required for lc_* split.")
    
    df %>%
      mutate(
        on_river = as.integer(!is.na(water_river) & water_river > 0),
        on_road = as.integer(!is.na(road) & road > 0),
        
        lc_1S  = ifelse(BCR %in% south_bcr, lc_1, 0),
        lc_1N  = ifelse(BCR %in% north_bcr, lc_1, 0),
        
        lc_4S  = ifelse(BCR %in% south_bcr, lc_4, 0),
        lc_4N  = ifelse(BCR %in% north_bcr, lc_4, 0),
        
        lc_5S  = ifelse(BCR %in% south_bcr, lc_5, 0),
        lc_5N  = ifelse(BCR %in% north_bcr, lc_5, 0),
        
        lc_8S  = ifelse(BCR %in% south_bcr, lc_8, 0),
        lc_8N  = ifelse(BCR %in% north_bcr, lc_8, 0),
        
        lc_9S  = ifelse(BCR %in% south_bcr, lc_9, 0),
        lc_9N  = ifelse(BCR %in% north_bcr, lc_9, 0),
        
        lc_10S = ifelse(BCR %in% south_bcr, lc_10, 0),
        lc_10N = ifelse(BCR %in% north_bcr, lc_10, 0)
      ) %>%
      mutate(onriver_S  = ifelse(BCR %in% south_bcr, on_river, 0),
             onriver_N  = ifelse(BCR %in% north_bcr, on_river, 0))
  }
  
  surveys <- make_bits(surveys)
  grid2   <- make_bits(grid2)
  grid3   <- make_bits(grid3)
  
  list(surveys = surveys, grid2 = grid2, grid3 = grid3)
}

# Several covariates will be divied into north and south components
south_bcr <- c(12, 13)
north_bcr <- c(7, 8)

tmp <- add_derived_covariates(surveys_f, grid_OBBA2, grid_OBBA3, south_bcr, north_bcr)
surveys_f <- tmp$surveys
grid_OBBA2  <- tmp$grid2
grid_OBBA3  <- tmp$grid3

# ------------------------------------------------------------
# 9) Species summary tables (for triage / reporting)
# ------------------------------------------------------------

all_species_unique <- dat$all_species %>%
  group_by(species_id) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(species_id = as.character(species_id)) %>%
  dplyr::select(species_id, english_name)

# counts_f has obs_idx + species columns
species_square_summary <- counts_f %>%
  pivot_longer(cols = -obs_idx, names_to = "species_id", values_to = "count") %>%
  left_join(st_drop_geometry(surveys_f), by = "obs_idx") %>%
  group_by(species_id, square_id, Atlas, BCR) %>%
  summarise(n_surveys_det = sum(count > 0, na.rm = TRUE), .groups = "drop") %>%
  group_by(species_id, Atlas, BCR) %>%
  summarise(
    n_det = sum(n_surveys_det, na.rm = TRUE),
    n_sq  = sum(n_surveys_det > 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(all_species_unique)

species_total_squares <- species_square_summary %>%
  group_by(species_id, Atlas) %>%
  summarise(
    total_detections = sum(n_det, na.rm = TRUE),
    total_squares    = sum(n_sq,  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = Atlas,
              values_from = c(total_detections, total_squares),
              values_fill = 0) %>%
  left_join(all_species_unique)

species_to_model <- species_total_squares %>%
  arrange(desc(total_squares_OBBA3)) %>%
  relocate(species_id,english_name)

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

saveRDS(
  list(
    study_boundary = study_boundary,
    all_surveys = surveys_f,
    counts = counts_f,
    grid_OBBA2 = grid_OBBA2,
    grid_OBBA3 = grid_OBBA3,
    species_square_summary = species_square_summary,
    species_total_squares = species_total_squares,
    species_to_model = species_to_model,
    covar_stats = covar_stats,
    date_created = Sys.time()
  ),
  file = out_file
)

message("07_filter_and_finalize_surveys.R complete: ", out_file)
