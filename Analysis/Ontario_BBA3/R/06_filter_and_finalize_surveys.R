# ============================================================
# 06_filter_and_finalize_surveys.R
#
# Build the analysis-ready survey, count, species-summary, prediction-grid,
# and covariate-raster products used by downstream atlas models.
#
# Important: obs_idx links filtered surveys back to rows of count_matrix.
# Keep that alignment unchanged when editing this script.
# ============================================================

# Clean session and load packages
rm(list = ls())

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(here)
  library(readxl)
})

# ------------------------------------------------------------------------------
# Paths and helper scripts
# ------------------------------------------------------------------------------

source(here::here("R", "00_config_paths.R"))
source(file.path(paths$functions, "survey_processing_utils.R"))
source(file.path(paths$functions, "covariate_processing_utils.R"))

dir.create(file.path(paths$data_clean, "birds"), recursive = TRUE, showWarnings = FALSE)

# Analysis settings
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

safe_date_path <- file.path(paths$data, "Bird_Data_Raw","OBBA3", "OnAtlasSafeDates_update_for_NatureCounts_2026-04-27.xlsx")

bcr_path       <- file.path(paths$data, "Spatial", "National", "BCR", "BCR_Terrestrial_master.shp")
Biol_Regions_path <- file.path(paths$data, "Spatial", "Ontario_Atlas_biol_regions", "Atlas_biol_regions.shp")

covars_to_rasterize <- c("ForestNeedleleaf","ForestBroadleaf","ForestMixed","Wetland","Cropland",
                         "Urban","On_Road",
                         "Grassland_BCR7_8","Grassland_BCR12_13",
                         "Shrubland_BCR13","Shrubland_BCR7_8_12",
                         "Lake_Lg","Lake_Sm_BCR7","Lake_Sm_BCR8","Lake_Sm_BCR12","Lake_Sm_BCR13",
                         "GreatLakes","HudsonBayCoast",
                         "River_Lg_BCR7","River_Lg_BCR13_12_8","River_Sm_BCR7","River_Sm_BCR13_12_8")

# ------------------------------------------------------------------------------
# Load covariate-enriched survey and prediction objects
# ------------------------------------------------------------------------------

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

# ------------------------------------------------------------------------------
# Validate required fields and create stable count table
# ------------------------------------------------------------------------------

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

# ------------------------------------------------------------------------------
# Filter surveys to the analysis window and supported survey types
# ------------------------------------------------------------------------------

surveys_f <- all_surveys %>%
  filter(
    round(Survey_Duration_Minutes) >= min_duration & round(Survey_Duration_Minutes) <= max_duration,
    Hours_Since_Sunrise >= hss_min,
    Hours_Since_Sunrise <= hss_max,
    DayOfYear >= doy_min,
    DayOfYear <= doy_max,
    Survey_Type %in% c("Point_Count","ARU","Breeding Bird Atlas","Linear transect"),
    
    !(Survey_Type == "Breeding Bird Atlas" & Distance_Traveled_m > 0)
  ) %>%
  arrange(obs_idx) %>%
  
  mutate(SC = ifelse(Survey_Type == "Breeding Bird Atlas",1,0),
         LT = ifelse(Survey_Type == "Linear transect",1,0))

# ------------------------------------------------------------------------------
# Keep surveys inside the study boundary
# ------------------------------------------------------------------------------

study_boundary <- st_transform(study_boundary, st_crs(surveys_f))
inside <- st_within(surveys_f, study_boundary, sparse = FALSE)[, 1]

surveys_f <- surveys_f[inside, , drop = FALSE] %>%
  arrange(obs_idx)

# ------------------------------------------------------------------------------
# Thin duplicate locations within year
# ------------------------------------------------------------------------------

n_orig <- nrow(surveys_f)
set.seed(dedupe_seed)
surveys_f <- dedupe_by_location(surveys_f, digits = dup_coord_digits, target_doy = 160, date_col = "Date_Time") %>%
  arrange(obs_idx)
n_thinned <- nrow(surveys_f)
n_removed <- n_orig - n_thinned
print(n_removed)

# ------------------------------------------------------------------------------
# Subset counts to the retained survey rows
# ------------------------------------------------------------------------------

stopifnot(max(surveys_f$obs_idx) <= nrow(counts))

counts_f <- counts %>%
  slice(surveys_f$obs_idx) %>%
  mutate(obs_idx = surveys_f$obs_idx)

stopifnot(nrow(counts_f) == nrow(surveys_f))

# ------------------------------------------------------------------------------
# Assign BCR membership used by derived habitat covariates
# ------------------------------------------------------------------------------

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
}

# ------------------------------------------------------------------------------
# Assign atlas biological regions used by safe-date summaries
# ------------------------------------------------------------------------------

if (file.exists(Biol_Regions_path)) {
  Biol_Regions_ON <- st_read(Biol_Regions_path, quiet = TRUE) %>%
    st_make_valid() %>%
    dplyr::rename(Biol_Region = Biol_Regio)
  
  surveys_f$Biol_Region <- assign_poly_id(surveys_f, Biol_Regions_ON, id_col = "Biol_Region", nearest_fallback = TRUE)
  
  grid_Biol_Region <- assign_poly_id(grid_OBBA2, Biol_Regions_ON %>% st_transform(st_crs(grid_OBBA2)),
                                     id_col = "Biol_Region", nearest_fallback = TRUE)
  grid_OBBA2$Biol_Region <- grid_Biol_Region
  grid_OBBA3$Biol_Region <- grid_Biol_Region
}

# Add survey-level model helper variables
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

surveys_f <- surveys_f %>%
  mutate(
    square_atlas = as.integer(factor(paste0(square_id, "-", Atlas)))
  )

# ------------------------------------------------------------------------------
# Create derived habitat covariates used by the models and maps
# ------------------------------------------------------------------------------

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
      
      Urban = as.integer(
        get_or_zero(urban_var) > 0.5
      ),
      
      On_Road = as.integer(get_or_zero("road") > 0 & Urban == 0),
      
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
      "Lake_Sm" = "lakes_small",
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

# ------------------------------------------------------------------------------
# Build survey metadata and long-format species counts
# ------------------------------------------------------------------------------

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

# ------------------------------------------------------------------------------
# Load and complete breeding safe-date windows
# Median level-1 and 2 safe dates within each Biol_Region
# NOTE: SOME SPECIES HAVE SAFE DATES THAT SPAN ACROSS JAN 1; THIS IS CURRENTLY NOT SUPPORTED
#       - by default it only uses safe dates up to Sep 1
# ------------------------------------------------------------------------------

# User controls
safe_levels <- c(1, 2)
min_safe_doy <- doy_min   
max_safe_doy <- doy_max 

Biol_Regions_df <- Biol_Regions_ON %>%
  as.data.frame() %>%
  dplyr::select(Biol_Region, ECOZONE_NA)

safe_dates <- readxl::read_xlsx(safe_date_path) %>%
  dplyr::rename(
    sp_english  = english_name,
    Biol_Region = biol_region,
    start_doy   = `start_dt (julian)`,
    end_doy     = `end_dt (julian)`
  )

# Collapse selected levels to one window per species x region.
safe_date_windows <- safe_dates %>%
  dplyr::filter(level %in% safe_levels) %>%
  dplyr::filter(start_doy < 240) %>%
  dplyr::group_by(sp_english, Biol_Region) %>%
  dplyr::summarise(
    start_doy = min(start_doy, na.rm = TRUE),
    end_doy   = max(end_doy, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    start_doy = pmax(start_doy, min_safe_doy),
    end_doy   = pmin(end_doy, max_safe_doy)
  ) %>%
  dplyr::filter(start_doy < end_doy)

species_fallbacks <- safe_date_windows %>%
  dplyr::group_by(sp_english) %>%
  dplyr::summarise(
    species_latest_start_doy = max(start_doy, na.rm = TRUE),
    species_earliest_end_doy = min(end_doy, na.rm = TRUE),
    .groups = "drop"
  )

region_medians <- safe_date_windows %>%
  dplyr::group_by(Biol_Region) %>%
  dplyr::summarise(
    median_start_doy = median(start_doy, na.rm = TRUE),
    median_end_doy   = median(end_doy, na.rm = TRUE),
    .groups = "drop"
  )

overall_median_start <- median(safe_date_windows$start_doy, na.rm = TRUE)
overall_median_end   <- median(safe_date_windows$end_doy, na.rm = TRUE)

all_sp_regions <- safe_dates %>%
  dplyr::distinct(sp_english, Biol_Region)

safe_dates_breeding <- all_sp_regions %>%
  dplyr::left_join(
    safe_date_windows,
    by = c("sp_english", "Biol_Region")
  ) %>%
  dplyr::left_join(species_fallbacks, by = "sp_english") %>%
  dplyr::left_join(region_medians, by = "Biol_Region") %>%
  dplyr::mutate(
    start_doy = dplyr::coalesce(
      start_doy,
      species_latest_start_doy,
      median_start_doy,
      overall_median_start
    ),
    end_doy = dplyr::coalesce(
      end_doy,
      species_earliest_end_doy,
      median_end_doy,
      overall_median_end
    ),
    midpoint = (start_doy + end_doy) / 2
  ) %>%
  dplyr::select(sp_english, Biol_Region, start_doy, end_doy, midpoint) %>%
  dplyr::full_join(Biol_Regions_df, by = "Biol_Region")

# ------------------------------------------------------------------------------
# Classify surveys by safe-date status
# ------------------------------------------------------------------------------

counts_long <- counts_long %>%
  left_join(
    safe_dates_breeding,
    by = c("english_name" = "sp_english", "Biol_Region" = "Biol_Region")
  ) %>%
  mutate(
    safe_date_status = case_when(
      is.na(start_doy) | is.na(end_doy) ~ "no_safe_dates",
      survey_doy >= start_doy & survey_doy <= end_doy ~ "inside",
      TRUE ~ "outside"
    )
  )

# ------------------------------------------------------------------------------
# Summarize detections by species, atlas, and biological region
# ------------------------------------------------------------------------------

species_Biol_Region_summary <- counts_long %>%
  group_by(species_id, english_name, Atlas, Biol_Region) %>%
  summarise(
    detections_total   = sum(count > 0, na.rm = TRUE),
    detections_safe    = sum(count > 0 & safe_date_status == "inside", na.rm = TRUE),
    detections_outside = sum(count > 0 & safe_date_status == "outside", na.rm = TRUE),
    detections_no_safe = sum(count > 0 & safe_date_status == "no_safe_dates", na.rm = TRUE),
    prop_outside = if_else(detections_total > 0, detections_outside / detections_total, NA_real_),
    prop_no_safe = if_else(detections_total > 0, detections_no_safe / detections_total, NA_real_),
    .groups = "drop"
  )

outside_safe_dates <- species_Biol_Region_summary %>%
  filter(prop_outside >= 0.1, detections_total >= 100) %>%
  arrange(english_name, Biol_Region, Atlas)

missing_safe_dates <- species_Biol_Region_summary %>%
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

# ------------------------------------------------------------------------------
# Build 25-km hex grid and pixel-to-hex lookup
# ------------------------------------------------------------------------------

hex_grid_25km <- make_hex_grid(study_boundary, width_km = 25)

grid_OBBA2 <- grid_OBBA2 %>% mutate(Atlas = "OBBA2")
grid_OBBA3 <- grid_OBBA3 %>% mutate(Atlas = "OBBA3")

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

# ------------------------------------------------------------------------------
# Save analysis-ready vector/list objects
# ------------------------------------------------------------------------------

saveRDS(
  list(
    study_boundary = study_boundary,
    bcr_sf = dat$bcr_sf,
    
    all_surveys = surveys_f,
    counts = counts_f,
    grid_OBBA2 = grid_OBBA2,
    grid_OBBA3 = grid_OBBA3,
    
    species_to_model = species_to_model,
    
    safe_dates_breeding = safe_dates_breeding,
    species_Biol_Region_summary = species_Biol_Region_summary,
    outside_safe_dates = outside_safe_dates,
    missing_safe_dates = missing_safe_dates,
    
    hex_grid_25km = hex_grid_25km,
    
    date_created = Sys.time()
  ),
  file = out_file
)

# ------------------------------------------------------------------------------
# Save covariate raster stacks later use if desired
# ------------------------------------------------------------------------------

# rast_dir  <- file.path(paths$data_clean, "spatial")
# rast_path_a2 <- file.path(rast_dir, "atlas2_cov_rasterstack.tif")
# rast_path_a3 <- file.path(rast_dir, "atlas3_cov_rasterstack.tif")
# 
# figure_utils_path <- file.path(paths$functions, "figure_utils.R")
# source(figure_utils_path)
# 
# r2 <- rasterize_sf(
#   grid_OBBA2,
#   covars_to_rasterize,
#   res = 1001,
#   metadata = c(description = c("Covariates at 1 km resolution"))
# )
# writeRaster(r2, filename = rast_path_a2, overwrite = TRUE)
# 
# r3 <- rasterize_sf(
#   grid_OBBA3,
#   covars_to_rasterize,
#   res = 1001,
#   metadata = c(description = c("Covariates at 1 km resolution"))
# )
# writeRaster(r3, filename = rast_path_a3, overwrite = TRUE)
# 
# dir.create(file.path(rast_dir,"A2_rasters"), recursive = TRUE, showWarnings = FALSE)
# terra::writeRaster(
#   r2,
#   filename = file.path(rast_dir,"A2_rasters", paste0(names(r2), ".tif")),
#   overwrite = TRUE
# )
# 
# dir.create(file.path(rast_dir,"A3_rasters"), recursive = TRUE, showWarnings = FALSE)
# terra::writeRaster(
#   r3,
#   filename = file.path(rast_dir,"A3_rasters", paste0(names(r3), ".tif")),
#   overwrite = TRUE
# )

message("06_filter_and_finalize_surveys.R complete")
