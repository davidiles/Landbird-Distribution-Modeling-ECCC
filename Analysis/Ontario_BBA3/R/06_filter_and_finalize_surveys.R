# ============================================================
# 06_filter_and_finalize_surveys.R
#
# Purpose
#   Create the final analysis-ready survey and prediction datasets used
#   by downstream modeling scripts.
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
#   6) Prepares prediction grids:
#        - flags optional open-water pixels
#        - drops cells missing core covariates
#        - ensures stable pixel_id values
#   8) Adds derived BCR-specific habitat covariates (e.g., effect of a land cover type differs between BCRs)
#   9) Builds survey metadata and a long-format count table used for
#      species summaries.
#  10) Loads breeding safe dates, expands them from ecoregions to BCRs,
#      and summarizes detections inside/outside safe dates.
#  11) Builds a 25-km hex grid for later regional summaries
#  12) Optionally writes covariate raster stacks for atlas-wide mapping.
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
#         * build_pixel_polygon_index()
#   - Optional:
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
#         * species_bcr_summary
#         * outside_safe_dates
#         * missing_safe_dates
#         * safe_dates
#         * safe_dates_bcr
#         * covar_stats
#         * hex_grid
#         * prediction_chunk_lookup
#         * date_created
#
# Notes
#   - obs_idx is treated as the row-position key linking filtered surveys
#     to the saved count matrix.
#   - dedupe_by_location() uses random sampling; a fixed seed is set here
#     for reproducibility.
#   - The same hex-based chunk lookup is reused later for chunked model
#     prediction, so it is created once here and saved.
# ============================================================


# NOTES:
# - update safe dates

rm(list = ls())

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
source(file.path(paths$functions, "inla_model_utils2.R"))
source(file.path(paths$functions, "spatial_utils.R"))
source(file.path(paths$functions, "covariate_processing_utils.R"))

dir.create(file.path(paths$data_clean, "birds"), recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Config
# ------------------------------------------------------------

in_file  <- file.path(paths$data_clean, "birds", "analysis_data_covariates.rds")
out_file <- file.path(paths$data_clean, "birds", "data_ready_for_analysis.rds")

min_duration <- 3
max_duration <- 10
hss_min <- -1
hss_max <- 6
doy_min <- 135
doy_max <- 196

dup_coord_digits <- 2 # 10 m accuracy
dedupe_seed <- 1

# NOTE: needs to be updated with new safe dates
safe_date_path <- file.path(paths$data_clean, "metadata", "safe_dates_OBBA.xlsx")

# Uses old BCR boundary definitions
bcr_path       <- file.path(paths$data, "Spatial", "National", "BCR", "BCR_Terrestrial_master.shp")

covars_to_rasterize <- c("ForestNeedleleaf","ForestBroadleaf","ForestMixed","Wetland","Cropland",
                         "Urban","On_Road",
                         "Grassland_BCR7_8","Grassland_BCR12_13",
                         "Shrubland_BCR13","Shrubland_BCR7_8_12",
                         "Lake_Lg","Lake_Sm_BCR7","Lake_Sm_BCR8","Lake_Sm_BCR12","Lake_Sm_BCR13",
                         "GreatLakes","HudsonBayCoast",
                         "River_Lg_BCR7","River_Lg_BCR13_12_8","River_Sm_BCR7","River_Sm_BCR13_12_8")

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

required_cols <- c("Survey_Duration_Minutes", "Hours_Since_Sunrise", "DayOfYear", "Date_Time", "Atlas","Survey_Type")
missing_cols <- setdiff(required_cols, names(all_surveys))
if (length(missing_cols) > 0) {
  stop("Missing required survey columns: ", paste(missing_cols, collapse = ", "))
}

# ------------------------------------------------------------
# 1) Filter surveys to the analysis window, and appropriate type
# ------------------------------------------------------------

surveys_f <- all_surveys %>%
  filter(
    round(Survey_Duration_Minutes) >= min_duration & round(Survey_Duration_Minutes) <= max_duration,
    Hours_Since_Sunrise >= hss_min,
    Hours_Since_Sunrise <= hss_max,
    DayOfYear >= doy_min,
    DayOfYear <= doy_max,
    Survey_Type %in% c("Point_Count","ARU","Breeding Bird Atlas","Linear transect"),
    
    # Remove "Breeding Bird Atlas" surveys with a non-zero distance (a tiny minority of them)
    !(Survey_Type == "Breeding Bird Atlas" & Distance_Traveled_m > 0)
  ) %>%
  arrange(obs_idx) %>%
  
  # Add indicator variables for stationary count (SC) and linear transect (LT) checklists
  mutate(SC = ifelse(Survey_Type == "Breeding Bird Atlas",1,0),
         LT = ifelse(Survey_Type == "Linear transect",1,0))
         

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

n_orig <- nrow(surveys_f)
set.seed(dedupe_seed)
surveys_f <- dedupe_by_location(surveys_f, digits = 2, target_doy = 160, date_col = "Date_Time") %>%
  arrange(obs_idx)
n_thinned <- nrow(surveys_f)
n_removed <- n_orig - n_thinned
print(n_removed)

# ------------------------------------------------------------
# 4) Subset counts to the retained surveys
# ------------------------------------------------------------

stopifnot(max(surveys_f$obs_idx) <= nrow(counts))

counts_f <- counts %>%
  slice(surveys_f$obs_idx) %>%
  mutate(obs_idx = surveys_f$obs_idx)

stopifnot(nrow(counts_f) == nrow(surveys_f))

# ------------------------------------------------------------
# 5) Create derived covariates, which include:
#          - categorical urban, road, river, and lake effects
#          - separate effects of grassland for BCRs 7/8 and 12/13
#          - separate effects of shrubland in BCR 13
#          - 
# ------------------------------------------------------------

if (file.exists(bcr_path)) {
  bcr_on <- st_read(bcr_path, quiet = TRUE) %>%
    st_make_valid() %>%
    select(BCR, PROVINCE_S) %>%
    filter(PROVINCE_S == "ONTARIO", BCR %in% c(7, 8, 12, 13)) %>%
    group_by(BCR) %>%
    summarise(geometry = st_union(geometry), .groups = "drop") %>%
    st_make_valid()

  surveys_f$BCR <- assign_poly_id(surveys_f, bcr_on, id_col = "BCR", nearest_fallback = TRUE)

  grid_bcr <- assign_poly_id(grid_OBBA2, bcr_on %>% st_transform(st_crs(grid_OBBA2)),
                             id_col = "BCR", nearest_fallback = TRUE)
  grid_OBBA2$BCR <- grid_bcr
  grid_OBBA3$BCR <- grid_bcr
} else {
  message("BCR shapefile not found; skipping BCR assignment.")
  all_surveys_with_covs$BCR <- NA_integer_
  grid_OBBA2$BCR <- NA_integer_
  grid_OBBA3$BCR <- NA_integer_
}

# Add some helper columns, including BCR membership
surveys_f <- surveys_f %>%
  mutate(
    Atlas3   = if_else(Atlas == "OBBA3", 1L, 0L),
    Atlas3_c = Atlas3 - 0.5,
    days_rescaled = DayOfYear - 166,
    BCR_idx = as.integer(factor(BCR))
  )

surveys_f <- surveys_f %>%
  mutate(ARU = if_else(Survey_Type == "ARU", 1L, 0L),
         Checklist = ifelse(Survey_Type %in% c("Breeding Bird Atlas","Linear transect"), 1L,0L))

# add pixel id and square id
surveys_f <- surveys_f %>%
  mutate(
    square_atlas = as.integer(factor(paste0(square_id, "-", Atlas)))
  )

add_derived_covariates <- function(df,
                                   urban_var = "LCC2020_17",
                                   grassland_var = "LCC2020_10",
                                   shrubland_var = "LCC2020_8") {
  
  if (!("BCR" %in% names(df))) {
    stop("BCR missing; required for derived covariates.")
  }
  
  get_or_zero <- function(nm) {
    if (nm %in% names(df)) df[[nm]] else rep(0, nrow(df))
  }
  
  df %>%
    mutate(
      
      # ------------------------------------------------------
      # Binary derived covariates
      # ------------------------------------------------------
      
      # If pixel has more than 50% Urban, call it an "Urban" pixel
      Urban = as.integer(
        get_or_zero(urban_var) > 0.5
      ),
      
      # Remove roads from Urban areas (Urban is assumed to include road effects)
      On_Road = as.integer(get_or_zero("road") > 0 & Urban == 0),
      
      # ------------------------------------------------------
      # Grassland split:
      # BCRs 7/8 (Likely represents burn scars / regen here) vs BCRs 12/13 (Likely represents grassland here)
      # ------------------------------------------------------
      
      Grassland_BCR7_8 = if_else(
        BCR %in% c(7, 8),
        get_or_zero(grassland_var),
        0
      ),
      
      Grassland_BCR12_13 = if_else(
        BCR %in% c(12, 13),
        get_or_zero(grassland_var),
        0
      ),
      
      # ------------------------------------------------------
      # Shrubland split:
      # BCR 13 (represents hedgerows/dogwood/etc) vs BCRs 7/8/12 (represents salix/betula/other dwarf shrubs)
      # ------------------------------------------------------
      
      Shrubland_BCR13 = if_else(
        BCR == 13,
        get_or_zero(shrubland_var),
        0
      ),
      
      Shrubland_BCR7_8_12 = if_else(
        BCR %in% c(7, 8, 12),
        get_or_zero(shrubland_var),
        0
      ),
      
      # ------------------------------------------------------
      # Seprate small lake effect for each BCR
      # ------------------------------------------------------
      
      Lake_Sm_BCR7 = if_else(
        BCR == 7,
        lakes_small,
        0L
      ),
      
      Lake_Sm_BCR8 = if_else(
        BCR == 8,
        lakes_small,
        0L
      ),
      
      Lake_Sm_BCR12 = if_else(
        BCR == 12,
        lakes_small,
        0L
      ),
      
      Lake_Sm_BCR13 = if_else(
        BCR == 13,
        lakes_small,
        0L
      ),
      
      # ------------------------------------------------------
      # Separate river effects in BCR 7 (Hudson Bay Lowlands)
      # Rivers represent a unique habitat type in HBL
      # ------------------------------------------------------
      
      River_Lg_BCR7 = if_else(
        BCR == 7,
        get_or_zero("rivers_large"),
        0
      ),
      
      River_Lg_BCR13_12_8 = if_else(
        BCR %in% c(13, 12, 8),
        get_or_zero("rivers_large"),
        0
      ),
      
      River_Sm_BCR7 = if_else(
        BCR == 7,
        get_or_zero("rivers_small"),
        0
      ),
      
      River_Sm_BCR13_12_8 = if_else(
        BCR %in% c(13, 12, 8),
        get_or_zero("rivers_small"),
        0
      ),
      
    ) %>%
    dplyr::rename(
      "GreatLakes" = "great_lakes",
      "HudsonBayCoast" = "coastline",
      "Lake_Lg" = "lakes_large",
      "River_Lg" = "rivers_large",
      "River_Sm" = "rivers_small",
      "ForestNeedleleaf" = "LCC2020_1",
      "ForestBroadleaf" = "LCC2020_5",  
      "ForestMixed" = "LCC2020_6", 
      "Wetland" = "LCC2020_14", 
      "Cropland" = "LCC2020_15"
    )
}

surveys_f  <- add_derived_covariates(surveys_f)
grid_OBBA2 <- add_derived_covariates(grid_OBBA2)
grid_OBBA3 <- add_derived_covariates(grid_OBBA3)

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
# 10) Load safe dates, expand them to BCR, and summarize species
# ------------------------------------------------------------

# These species are given custom safe dates (medians of other species)
species_manual_safedates <- c(
  "Canada Jay", "Common Raven", "White-winged Crossbill", "Red Crossbill"
)

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

safe_dates_bcr <- safe_dates %>%
  mutate(
    BCR = case_when(
      ecoregion == "Mixedwood Plains" ~ "13",
      ecoregion == "Hudson Plains"    ~ "7",
      ecoregion == "Boreal Shield"    ~ "8, 12",
      TRUE ~ NA_character_
    )
  ) %>%
  separate_rows(BCR, sep = ",") %>%
  mutate(BCR = as.integer(trimws(BCR))) %>%
  select(sp_english, BCR, start_doy, end_doy) %>%
  arrange(sp_english, BCR) %>%
  mutate(midpoint = (start_doy + end_doy) / 2)

counts_long <- counts_long %>%
  left_join(
    safe_dates_bcr,
    by = c("english_name" = "sp_english", "BCR" = "BCR")
  ) %>%
  mutate(
    safe_date_status = case_when(
      is.na(start_doy) | is.na(end_doy) ~ "no_safe_dates",
      survey_doy >= start_doy & survey_doy <= end_doy ~ "inside",
      TRUE ~ "outside"
    )
  )

# Summarize detections inside / outside safe date windows
species_bcr_summary <- counts_long %>%
  group_by(species_id, english_name, Atlas, BCR) %>%
  summarise(
    detections_total   = sum(count > 0, na.rm = TRUE),
    detections_safe    = sum(count > 0 & safe_date_status == "inside", na.rm = TRUE),
    detections_outside = sum(count > 0 & safe_date_status == "outside", na.rm = TRUE),
    detections_no_safe = sum(count > 0 & safe_date_status == "no_safe_dates", na.rm = TRUE),
    prop_outside = if_else(detections_total > 0, detections_outside / detections_total, NA_real_),
    prop_no_safe = if_else(detections_total > 0, detections_no_safe / detections_total, NA_real_),
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
    n_squares_safe  = n_distinct(square_id[count > 0 & safe_date_status == "inside"]),
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
# 11) Build coarse hex grid for storing individual posterior draws (and presenting change estimates)
#     Individual posterior draws are saved only at the 25-km hexagon level
#     (Note: can swap in atlas squares here)
# ------------------------------------------------------------

hex_grid_25km <- make_hex_grid(study_boundary, width_km = 25)

grid_OBBA2 <- grid_OBBA2 %>% mutate(Atlas = "OBBA2")
grid_OBBA3 <- grid_OBBA3 %>% mutate(Atlas = "OBBA3")

# Assign each prediction pixel to a hexagon
hex_idx <- build_pixel_polygon_index(
  grid_sf     = grid_OBBA2,
  polygons_sf = hex_grid_25km,
  poly_id_col = "hex_id",
  join        = "within"
)

prediction_hex_lookup <- tibble(
  pixel_id = grid_OBBA2$pixel_id,
  hex_id   = hex_idx$pix_poly_id
) %>%
  filter(!is.na(hex_id))

grid_OBBA2 <- grid_OBBA2 %>% left_join(.,prediction_hex_lookup)
grid_OBBA3 <- grid_OBBA3 %>% left_join(.,prediction_hex_lookup)

grid_OBBA2 <- grid_OBBA2 %>% na.omit()
grid_OBBA3 <- grid_OBBA3 %>% na.omit()

# ------------------------------------------------------------
# Save vector data
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
    
    safe_dates = safe_dates,
    safe_dates_bcr = safe_dates_bcr,
    
    hex_grid_25km = hex_grid_25km,
    
    date_created = Sys.time()
  ),
  file = out_file
)

# ------------------------------------------------------------
# Save covariate rasters
# ------------------------------------------------------------

rast_dir  <- file.path(paths$data_clean, "spatial")
rast_path_a2 <- file.path(rast_dir, "atlas2_cov_rasterstack.tif")
rast_path_a3 <- file.path(rast_dir, "atlas3_cov_rasterstack.tif")

figure_utils_path <- file.path(paths$functions, "figure_utils.R")
source(figure_utils_path)

# Covariate raster stack for atlas 2
r2 <- rasterize_sf(
  grid_OBBA2,
  covars_to_rasterize,
  res = 1001,
  metadata = c(description = c("Covariates at 1 km resolution"))
)
# Save as raster stack
writeRaster(r2, filename = rast_path_a2, overwrite = TRUE)

# Covariate raster stack for atlas 3
r3 <- rasterize_sf(
  grid_OBBA3,
  covars_to_rasterize,
  res = 1001,
  metadata = c(description = c("Covariates at 1 km resolution"))
)
# Save as raster stack
writeRaster(r3, filename = rast_path_a3, overwrite = TRUE)


# Also save as separate rasters (for import into GIS applications)
dir.create(file.path(rast_dir,"A2_rasters"), recursive = TRUE, showWarnings = FALSE)
terra::writeRaster(
  r2,
  filename = file.path(rast_dir,"A2_rasters", paste0(names(r2), ".tif")),
  overwrite = TRUE
)

dir.create(file.path(rast_dir,"A3_rasters"), recursive = TRUE, showWarnings = FALSE)
terra::writeRaster(
  r3,
  filename = file.path(rast_dir,"A3_rasters", paste0(names(r3), ".tif")),
  overwrite = TRUE
)


message("06_filter_and_finalize_surveys.R complete")