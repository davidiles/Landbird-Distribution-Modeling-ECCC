# ============================================================
# 02_load_and_format_point_counts.R
#
# Purpose:
#   - Load and harmonize NatureCounts OBBA2 and OBBA3 point count data
#   - Remove unusable surveys (missing date/time/lat/lon/species/count)
#   - Construct survey-level sf object
#   - Construct survey × species count matrix
#   - Compute local time zone and hours since sunrise
#
# Inputs:
#   - data_clean/spatial/study_area.rds
#   - ../../Data/Bird_Data_Raw/OBBA3/Point_Count/onatlas3pc_naturecounts_data.txt
#   - ../../Data/Bird_Data_Raw/OBBA2/Point_Count/obba2pc_naturecounts_data.txt
#   - R/functions/spatial_utils.R
#   - R/functions/survey_processing_utils.R
#
# Outputs:
#   - data_clean/surveys/surveys_raw.rds
#   - data_clean/surveys/count_matrix_raw.rds
#   - data_clean/metadata/species_list.rds
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(sf)
  library(lubridate)
  library(naturecounts)
  library(purrr)
  library(readr)
  library(tidyr)
})

# ------------------------------------------------------------
# Load utilities and study area
# ------------------------------------------------------------

source("R/functions/spatial_utils.R")
source("R/functions/survey_processing_utils.R")

study_area <- readRDS("data_clean/spatial/study_area.rds")
crs_aea_km <- study_area$crs

# ------------------------------------------------------------
# Helper: build POSIXct datetime from date + fractional hours
# ------------------------------------------------------------

make_datetime_from_frac_hours <- function(date_ymd, frac_hours, tz = "UTC") {
  # date_ymd: Date
  # frac_hours: numeric, e.g., 13.5 means 13:30
  stopifnot(inherits(date_ymd, "Date"))
  
  # Convert to seconds (avoid minute==60 edge cases)
  secs <- round(frac_hours * 3600)
  as.POSIXct(date_ymd, tz = tz) + secs
}

# ------------------------------------------------------------
# Species list (NatureCounts canonical)
# ------------------------------------------------------------

all_species <- search_species_code() %>%
  rename(
    sp_code = BSCDATA,
    species_scientific_name = scientific_name
  )

english_lookup <- all_species %>%
  group_by(species_id) %>%
  summarise(
    english_name = dplyr::first(english_name),
    .groups = "drop"
  )

# ------------------------------------------------------------
# Load and process OBBA3 data
# ------------------------------------------------------------

obba3_path <- "../../Data/Bird_Data_Raw/OBBA3/Point_Count/onatlas3pc_naturecounts_data.txt"
NC_OBBA3_raw <- read.table(
  obba3_path,
  sep = "\t", header = TRUE, fill = TRUE, comment.char = "", quote = ""
) %>%
  select(where(~ sum(!is.na(.x)) > 0)) %>%
  rename(species_id = SpeciesCode) %>%
  mutate(
    Date_Time = make_datetime_from_frac_hours(
      ymd(ObservationDate),
      TimeObservationsStarted,
      tz = "UTC"
    ),
    Project_Name = "OBBA3",
    Data_Source = "NatureCounts",
    Survey_Type = infer_survey_type_OBBA3(
      Remarks, Remarks2, EffortMeasurement1, SurveyAreaIdentifier
    ),
    Special_Survey = if_else(EffortMeasurement1 == "Special", "Special", "No")
  )

NC_OBBA3 <- NC_OBBA3_raw %>%
  select(
    Data_Source, Project_Name, Survey_Type, Special_Survey,
    DecimalLatitude, DecimalLongitude,
    Date_Time, DurationInHours,
    species_id, ObservationCount
  )

# ------------------------------------------------------------
# Load and process OBBA2 data
# ------------------------------------------------------------

obba2_path <- "../../Data/Bird_Data_Raw/OBBA2/Point_Count/obba2pc_naturecounts_data.txt"
NC_OBBA2_raw <- read.table(
  obba2_path,
  sep = "\t", header = TRUE, fill = TRUE, comment.char = "", quote = ""
) %>%
  select(where(~ sum(!is.na(.x)) > 0)) %>%
  left_join(
    all_species[, c("sp_code", "species_id")],
    by = c("SpeciesCode" = "sp_code"),
    multiple = "first"
  ) %>%
  rename(TimeObservationsStarted = TimeCollected) %>%
  mutate(
    Date_Time = make_datetime_from_frac_hours(
      ymd(paste(YearCollected, MonthCollected, DayCollected, sep = "-")),
      TimeObservationsStarted,
      tz = "UTC"
    ),
    Project_Name = "OBBA2",
    Data_Source = "NatureCounts",
    Survey_Type = "Point_Count",
    Special_Survey = "No"
  )

NC_OBBA2 <- NC_OBBA2_raw %>%
  select(
    Data_Source, Project_Name, Survey_Type, Special_Survey,
    DecimalLatitude, DecimalLongitude,
    Date_Time, DurationInHours,
    species_id, ObservationCount
  )

# ------------------------------------------------------------
# Combine atlases and basic cleaning
# ------------------------------------------------------------

# NOTE: the longitude filter removes obvious data errors (e.g., missing sign).
# A stricter spatial filter (within Ontario boundary) occurs later in the pipeline.
NC_long <- bind_rows(NC_OBBA3, NC_OBBA2) %>%
  filter(
    !is.na(DecimalLatitude),
    !is.na(DecimalLongitude),
    DecimalLongitude < -60,
    !is.na(species_id),
    !is.na(ObservationCount),
    !is.na(Date_Time)
  ) %>%
  rename(
    Latitude = DecimalLatitude,
    Longitude = DecimalLongitude
  ) %>%
  mutate(
    # enforce stable ID construction across reruns/adaptations
    Latitude  = round(Latitude, 5),
    Longitude = round(Longitude, 5),
    Survey_Duration_Minutes = round(DurationInHours * 60),
    Max_Distance_Metres = Inf,
    survey_id = make_survey_id(Project_Name, Latitude, Longitude, Date_Time, Survey_Type)
  )

# ------------------------------------------------------------
# Fix several species_id discrepancies within/between atlases
# ------------------------------------------------------------

# Common (id = 40919), Hoary (id = 20400), and undistinguished (id = 45264) Redpoll should get lumped
NC_long$species_id[NC_long$species_id %in% c(40919,20400,45264)] <- 45264 

# Northern flicker: Colaptes auratus auratus/luteus (id = 10470) should get grouped with Colaptes auratus (id = 48798)
NC_long$species_id[NC_long$species_id %in% c(10470,48798)] <- 48798 

# Common/Wilson's Snipe should get lumped: Gallinago gallinago (id = 4940) and Gallinago delicata (id = 4950)
NC_long$species_id[NC_long$species_id %in% c(4940,4950)] <- 4950

# Yellow-rumped warbler: Setophaga coronata coronata (id = 16620) should get grouped with YRWA (id = 16610)
NC_long$species_id[NC_long$species_id %in% c(16620,16610)] <- 16610

NC_long$species_id[NC_long$species_id %in% c(40876,19530)] <- 19530


# ------------------------------------------------------------
# Sort species by relative abundance;
# used to identify species that could be mis-labeled in one atlas
# ------------------------------------------------------------

species_wide <- NC_long %>%
  
  # Add species names
  left_join(.,unique(all_species[,c("species_id","species_scientific_name")])) %>%


  group_by(species_id, Project_Name) %>%
  summarize(sum_count = sum(ObservationCount, na.rm = TRUE), .groups = "drop") %>%
  # wide: one row per species, one column per atlas
  pivot_wider(
    names_from  = Project_Name,
    values_from = sum_count,
    values_fill = 0
  ) %>%
  # differences / effect sizes
  mutate(
    
    total = OBBA2+OBBA3,
    
    # absolute difference (for "largest raw change")
    abs_diff = abs(OBBA3 - OBBA2),
    
    # signed difference (positive = higher in OBBA3)
    diff = OBBA3 - OBBA2,
    
    # fold-change (add +1 to avoid division by zero)
    fold = (OBBA3 + 1) / (OBBA2 + 1),
    
    # log fold-change (symmetric, nicer for ranking big relative shifts)
    log_fold = log(fold)
  )

# Examine species with "large" differences (+/- 30%) and at least 100 detections
# extreme differences could imply differences in species_id between atlases

large_differences <- species_wide %>%
  left_join(unique(all_species[,c("species_id","species_scientific_name")]))%>%
  filter(total>=100 & abs(log_fold) >= abs(log(0.7))) %>%
  arrange(desc(total)) %>%
  
  # Attach english names
  left_join(english_lookup)

large_differences %>% as.data.frame()

# ------------------------------------------------------------
# Construct survey-level table
# ------------------------------------------------------------

survey_info <- NC_long %>%
  group_by(survey_id) %>%
  summarise(
    Project_Name = first(Project_Name),
    Data_Source  = first(Data_Source),
    Survey_Type  = first(Survey_Type),
    Special_Survey = first(Special_Survey),
    Latitude     = first(Latitude),
    Longitude    = first(Longitude),
    Date_Time    = first(Date_Time),
    DurationInHours = median(DurationInHours, na.rm = TRUE),
    Survey_Duration_Minutes = round(DurationInHours * 60),
    Max_Distance_Metres = Inf,
    .groups = "drop"
  ) %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE)

# ------------------------------------------------------------
# Construct survey × species count matrix
# ------------------------------------------------------------

count_matrix <- build_count_matrix(
  NC_long,
  survey_id_col = "survey_id",
  species_col   = "species_id",
  count_col     = "ObservationCount"
)

# ------------------------------------------------------------
# Time zone and hours since sunrise (can be slow)
# ------------------------------------------------------------

survey_info <- add_hours_since_sunrise(survey_info)

# ------------------------------------------------------------
# CRS transform
# ------------------------------------------------------------

survey_info <- survey_info %>%
  mutate(Atlas = Project_Name) %>%
  st_transform(crs_aea_km) %>%
  relocate(geometry, .after = last_col())

# ------------------------------------------------------------
# Sanity checks
# ------------------------------------------------------------

stopifnot(
  length(unique(NC_long$survey_id)) == nrow(survey_info),
  length(unique(rownames(count_matrix))) == nrow(survey_info)
)

# ------------------------------------------------------------
# Save outputs
# ------------------------------------------------------------

dir.create("data_clean/surveys", recursive = TRUE, showWarnings = FALSE)
dir.create("data_clean/metadata", recursive = TRUE, showWarnings = FALSE)

saveRDS(survey_info, "data_clean/surveys/surveys_raw.rds")
saveRDS(count_matrix, "data_clean/surveys/count_matrix_raw.rds")
saveRDS(all_species, "data_clean/metadata/species_list.rds")

message("02_load_format_point_counts.R complete.")