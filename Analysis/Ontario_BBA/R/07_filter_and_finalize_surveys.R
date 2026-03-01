# ============================================================
# 07_filter_and_finalize_surveys.R
#
# Purpose
#   Create the final, analysis-ready dataset for modeling by applying *all*
#   survey filters and final bookkeeping in one place. This is the last step
#   before model fitting and cross-validation scripts.
#
# What this script does
#   1) Loads the covariate-augmented surveys and grids from:
#        data_clean/birds/analysis_data_covariates.rds   (script 06)
#   2) Filters surveys to a consistent “analysis window”:
#        - duration == 5 minutes
#        - Hours_Since_Sunrise in [-1, 6]
#        - DayOfYear in [135, 196]
#        - within the Ontario study boundary
#   3) Removes duplicate sampling locations (keeps one survey per rounded x/y).
#      NOTE: dedupe_by_location() uses random sampling; we set a seed here so
#      results are reproducible.
#   4) Subsets the survey × species count matrix to match the filtered surveys.
#      IMPORTANT: in the current workflow, the link between surveys and the
#      count matrix is row order (via obs_idx). This script preserves and uses
#      obs_idx consistently to maintain alignment.
#   5) Adds helper variables used by multiple downstream scripts:
#        - Atlas3 indicator (0/1 for OBBA2/OBBA3)
#        - ARU indicator (0/1)
#        - square_atlas and square_year (factor IDs for random effects)
#   6) (Optional) Flags prediction-grid points overlapping large-water polygons
#      if data_clean/spatial/water_filtered.shp is available.
#   7) Standardizes selected continuous covariates (currently prec and tmax)
#      using means/sds computed from the filtered surveys (not the full grid).
#   8) Creates derived “north vs south BCR” covariate splits used in some model
#      formulations (e.g., lc_1S vs lc_1N), plus simple binary on_river/on_road
#      indicators.
#   9) Builds species summary tables (detections and occupied squares) to help
#      triage which species to model and to support reporting/QA.
#
# Inputs
#   - data_clean/birds/analysis_data_covariates.rds   (from script 06)
#       must contain:
#         * all_surveys_with_covs (sf)
#         * count_matrix (matrix or df)
#         * grid_OBBA2, grid_OBBA3 (sf points; EPSG:3978)
#         * boundary (sf; Ontario boundary)
#   - R/functions/survey_processing_utils.R
#       used functions:
#         * as_counts_tbl()
#         * dedupe_by_location()
#         * standardize_covars()
#   - Optional:
#       data_clean/spatial/water_filtered.shp
#
# Outputs
#   - data_clean/birds/data_ready_for_analysis.rds
#       list with:
#         * study_boundary (sf)
#         * all_surveys (sf; filtered + finalized)
#         * counts (tibble; aligned with all_surveys by obs_idx)
#         * grid_OBBA2, grid_OBBA3 (sf; covariates + optional on_water flag)
#         * species_square_summary (tbl)
#         * species_total_squares (tbl)
#         * species_to_model (tbl)
#         * covar_stats (tbl; means/sds used for standardization)
#         * date_created
#
# Key design choices / conventions
#   - CRS: surveys may be stored in km-units CRS for modeling, but the prediction
#     grids are expected to be in EPSG:3978 (metres). This script does not
#     reproject grids; it assumes they are already in EPSG:3978 from script 06.
#   - Reproducibility: because dedupe_by_location() uses slice_sample(), set a
#     fixed seed before de-duplication to ensure deterministic outputs.
#   - Alignment: obs_idx is treated as the stable “row-position key” linking
#     surveys and the count matrix within this workflow. If you modify earlier
#     scripts, ensure this alignment remains valid.
#
# Common adaptations for other atlases / agencies
#   - Adjust the seasonal window (DayOfYear) to match local breeding phenology.
#   - Adjust the sunrise window (Hours_Since_Sunrise) for survey protocols.
#   - Modify duplicate-location handling (rounding scale and selection rule).
#   - Update BCR groupings (north/south) or remove these derived covariates if
#     covariate effects should not differ by region
#
# ============================================================

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(tidyr)
  library(lubridate)
})

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

# Duplicate location handling (rounding in projected coordinates; kiolmetre CRS)
# digits = 2 means nearest 10 m; 1 nearest 100 m.
dup_coord_digits <- 2

# Make dedupe reproducible (dedupe_by_location uses slice_sample)
dedupe_seed <- 1

covars_for_stats <- c(
  "insect_broadleaf","insect_needleleaf","road","water_river","prec","tmax",
  "urban_2","urban_3",
  "lc_1","lc_4","lc_5","lc_8","lc_9","lc_10","lc_11","lc_12","lc_14","lc_17"
)
covars_to_standardize <- c("prec","tmax")

# Optional water polygon layer path
water_filtered_path <- "data_clean/spatial/water_filtered.shp"

# ------------------------------------------------------------
# Load input object
# ------------------------------------------------------------

stopifnot(file.exists(in_file))
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
  sf::st_read(water_filtered_path, quiet = TRUE)
} else {
  NULL
}

# Ensure obs_idx exists and is stable
if (!("obs_idx" %in% names(all_surveys))) {
  all_surveys <- all_surveys %>% mutate(obs_idx = row_number(), .before = 1)
} else {
  all_surveys <- all_surveys %>% arrange(obs_idx)
}

# Convert count matrix to tibble with row-position obs_idx = 1..N
# IMPORTANT: this obs_idx reflects row order in the saved matrix,
# not any persistent survey identifier.
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
if (length(missing) > 0) {
  stop("Missing required survey columns: ", paste(missing, collapse = ", "))
}

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

study_boundary <- st_transform(study_boundary, st_crs(surveys_f))
inside <- st_within(surveys_f, study_boundary, sparse = FALSE)[, 1]
surveys_f <- surveys_f[inside, , drop = FALSE] %>% arrange(obs_idx)

# ------------------------------------------------------------
# 3) Retain a single survey per repeated location (reproducible)
# ------------------------------------------------------------

set.seed(dedupe_seed)
surveys_f <- dedupe_by_location(surveys_f, digits = dup_coord_digits) %>%
  arrange(obs_idx)

# ------------------------------------------------------------
# 4) Subset counts to remaining surveys
# ------------------------------------------------------------
# Because as_counts_tbl() defines obs_idx = row_number(), we subset by row position.
# This is valid IF the original count_matrix rows were in the same order as all_surveys
# when they were saved together in script 06 (which they should be).
stopifnot(max(surveys_f$obs_idx) <= nrow(counts))

counts_f <- counts %>%
  slice(surveys_f$obs_idx)

# Ensure alignment
stopifnot(nrow(counts_f) == nrow(surveys_f))

# Replace counts_f$obs_idx with the survey obs_idx for downstream joins clarity
counts_f$obs_idx <- surveys_f$obs_idx

# ------------------------------------------------------------
# 5) Add helper columns
# ------------------------------------------------------------

if (!("Atlas" %in% names(surveys_f))) stop("Missing 'Atlas' column in surveys (expected OBBA2/OBBA3).")

surveys_f <- surveys_f %>%
  mutate(
    Atlas3 = ifelse(Atlas == "OBBA2", 0L, 1L),
    ARU    = ifelse("Survey_Type" %in% names(surveys_f) & Survey_Type == "ARU", 1L, 0L)
  )

if ("square_id" %in% names(surveys_f)) {
  surveys_f <- surveys_f %>%
    mutate(
      square_atlas = as.numeric(factor(paste0(square_id, "-", Atlas))),
      square_year  = as.numeric(factor(paste0(square_id, "-", year(Date_Time))))
    )
} else {
  surveys_f$square_atlas <- NA_real_
  surveys_f$square_year  <- NA_real_
}

# ------------------------------------------------------------
# 6) Grid: on_water flag + NA handling + pixel_id
# ------------------------------------------------------------

if (!is.null(water_filtered)) {
  water_filtered <- st_make_valid(water_filtered)
  water_filtered <- st_transform(water_filtered, st_crs(grid_OBBA2))
  
  overlaps2 <- st_intersects(grid_OBBA2, water_filtered)
  overlaps3 <- st_intersects(grid_OBBA3, water_filtered)
  
  grid_OBBA2$on_water <- lengths(overlaps2) > 0
  grid_OBBA3$on_water <- lengths(overlaps3) > 0
} else {
  grid_OBBA2$on_water <- NA
  grid_OBBA3$on_water <- NA
}

# Avoid na.omit() removing rows due to optional covariates missing
core_grid_vars <- intersect(c("prec", "tmax"), names(grid_OBBA2))
if (length(core_grid_vars) > 0) {
  grid_OBBA2 <- grid_OBBA2 %>% filter(if_all(all_of(core_grid_vars), ~ !is.na(.x)))
  grid_OBBA3 <- grid_OBBA3 %>% filter(if_all(all_of(core_grid_vars), ~ !is.na(.x)))
}

# Preserve existing pixel_id if present
if (!("pixel_id" %in% names(grid_OBBA2))) grid_OBBA2$pixel_id <- seq_len(nrow(grid_OBBA2))
if (!("pixel_id" %in% names(grid_OBBA3))) grid_OBBA3$pixel_id <- seq_len(nrow(grid_OBBA3))

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

surveys_f   <- std$surveys
grid_OBBA2  <- std$grid2
grid_OBBA3  <- std$grid3
covar_stats <- std$stats

# ------------------------------------------------------------
# 8) Derived covariates (north/south components)
# ------------------------------------------------------------

add_derived_covariates <- function(df, south_bcr = c(12, 13), north_bcr = c(7, 8)) {
  if (!("BCR" %in% names(df))) stop("BCR missing; required for lc_* split.")
  get0 <- function(nm) if (nm %in% names(df)) df[[nm]] else 0
  
  df %>%
    mutate(
      on_river = as.integer(!is.na(get0("water_river")) & get0("water_river") > 0),
      on_road  = as.integer(!is.na(get0("road")) & get0("road") > 0),
      
      lc_1S  = ifelse(BCR %in% south_bcr, get0("lc_1"),  0),
      lc_1N  = ifelse(BCR %in% north_bcr, get0("lc_1"),  0),
      lc_4S  = ifelse(BCR %in% south_bcr, get0("lc_4"),  0),
      lc_4N  = ifelse(BCR %in% north_bcr, get0("lc_4"),  0),
      lc_5S  = ifelse(BCR %in% south_bcr, get0("lc_5"),  0),
      lc_5N  = ifelse(BCR %in% north_bcr, get0("lc_5"),  0),
      lc_8S  = ifelse(BCR %in% south_bcr, get0("lc_8"),  0),
      lc_8N  = ifelse(BCR %in% north_bcr, get0("lc_8"),  0),
      lc_9S  = ifelse(BCR %in% south_bcr, get0("lc_9"),  0),
      lc_9N  = ifelse(BCR %in% north_bcr, get0("lc_9"),  0),
      lc_10S = ifelse(BCR %in% south_bcr, get0("lc_10"), 0),
      lc_10N = ifelse(BCR %in% north_bcr, get0("lc_10"), 0),
      
      onriver_S = ifelse(BCR %in% south_bcr, on_river, 0),
      onriver_N = ifelse(BCR %in% north_bcr, on_river, 0)
    )
}

south_bcr <- c(12, 13)
north_bcr <- c(7, 8)

surveys_f  <- add_derived_covariates(surveys_f,  south_bcr, north_bcr)
grid_OBBA2 <- add_derived_covariates(grid_OBBA2, south_bcr, north_bcr)
grid_OBBA3 <- add_derived_covariates(grid_OBBA3, south_bcr, north_bcr)

# ------------------------------------------------------------
# 9) Species summary tables
# ------------------------------------------------------------

all_species_unique <- dat$all_species %>%
  group_by(species_id) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(species_id = as.character(species_id)) %>%
  select(species_id, english_name)

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
  left_join(all_species_unique, by = "species_id")

species_total_squares <- species_square_summary %>%
  group_by(species_id, Atlas) %>%
  summarise(
    total_detections = sum(n_det, na.rm = TRUE),
    total_squares    = sum(n_sq,  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Atlas,
    values_from = c(total_detections, total_squares),
    values_fill = 0
  ) %>%
  left_join(all_species_unique, by = "species_id")

species_to_model <- species_total_squares %>%
  arrange(desc(total_squares_OBBA3)) %>%
  relocate(species_id, english_name)

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