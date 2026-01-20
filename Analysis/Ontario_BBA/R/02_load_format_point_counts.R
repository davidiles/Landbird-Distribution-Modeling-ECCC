# ============================================================
# 02_load_and_format_point_counts.R
#
# Purpose:
#   - Load and harmonize NatureCounts OBBA2 and OBBA3 data
#   - Remove unusable surveys (those missing dates/times/lat/lon)
#   - Construct survey-level sf object
#   - Construct survey × species count matrix
#   - Compute local time zone and hours since sunrise
#
# Outputs:
#   data_clean/surveys/surveys_raw.rds
#   data_clean/surveys/count_matrix_raw.rds
#   data_clean/metadata/species_list.rds
# ============================================================

library(tidyverse)
library(sf)
library(lubridate)
library(lutz)
library(suncalc)
library(naturecounts)
library(purrr)

# ------------------------------------------------------------
# Load utilities and study area
# ------------------------------------------------------------

source("R/functions/spatial_utils.R")
source("R/functions/survey_processing_utils.R")

study_area <- readRDS("data_clean/spatial/study_area.rds")
crs_aea_km <- study_area$crs

# ------------------------------------------------------------
# Species list (NatureCounts canonical)
# ------------------------------------------------------------

all_species <- search_species_code() %>%
  rename(
    sp_code = BSCDATA,
    species_scientific_name = scientific_name
  )

# ------------------------------------------------------------
# Load and process OBBA3 data
# ------------------------------------------------------------

NC_OBBA3 <- read.table(
  "../../Data/Bird_Data_Raw/OBBA3/Point_Count/onatlas3pc_naturecounts_data.txt",
  sep = "\t", header = TRUE, fill = TRUE, comment.char = "", quote = ""
) %>%
  select(where(~ sum(!is.na(.x)) > 0)) %>%
  rename(species_id = SpeciesCode) %>%
  mutate(
    Date_Time = ymd(ObservationDate) +
      hours(floor(TimeObservationsStarted)) +
      minutes(round(60 * (TimeObservationsStarted %% 1))),
    Project_Name = "OBBA3",
    Data_Source = "NatureCounts",
    Survey_Type = infer_survey_type_OBBA3(
      Remarks, Remarks2, EffortMeasurement1, SurveyAreaIdentifier
    )
  ) %>%
  select(
    Data_Source, Project_Name, Survey_Type,
    DecimalLatitude, DecimalLongitude,
    Date_Time, DurationInHours,
    species_id, ObservationCount
  )

# ------------------------------------------------------------
# Load and process OBBA2 data
# ------------------------------------------------------------

NC_OBBA2 <- read.table(
  "../../Data/Bird_Data_Raw/OBBA2/Point_Count/obba2pc_naturecounts_data.txt",
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
    Date_Time = ymd(paste(YearCollected, MonthCollected, DayCollected, sep = "-")) +
      hours(floor(TimeObservationsStarted)) +
      minutes(round(60 * (TimeObservationsStarted %% 1))),
    Project_Name = "OBBA2",
    Data_Source = "NatureCounts",
    Survey_Type = "Point_Count"
  ) %>%
  select(
    Data_Source, Project_Name, Survey_Type,
    DecimalLatitude, DecimalLongitude,
    Date_Time, DurationInHours,
    species_id, ObservationCount
  )

# ------------------------------------------------------------
# Combine atlases and basic cleaning
# ------------------------------------------------------------

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
    Survey_Duration_Minutes = round(DurationInHours * 60),
    Max_Distance_Metres = Inf,
    survey_id = make_survey_id(Project_Name, Latitude, Longitude, Date_Time,Survey_Type)
  )

# ------------------------------------------------------------
# Construct survey-level table
# ------------------------------------------------------------

survey_info <- NC_long %>%
  group_by(survey_id) %>%
  summarise(
    Project_Name = first(Project_Name),
    Data_Source  = first(Data_Source),
    Survey_Type  = first(Survey_Type),
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
# Time zone and hours since sunrise (takes a while)
# ------------------------------------------------------------

survey_info <- add_hours_since_sunrise(survey_info)

# ------------------------------------------------------------
# Atlas indicator and CRS transform
# ------------------------------------------------------------

survey_info <- survey_info %>%
  mutate(
    Atlas = ifelse(year(Date_Time) > 2020, "OBBA3", "OBBA2")
  ) %>%
  st_transform(crs_aea_km) %>%
  dplyr::relocate(geometry, .after = dplyr::last_col())

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

saveRDS(
  survey_info,
  "data_clean/surveys/surveys_raw.rds"
)

saveRDS(
  count_matrix,
  "data_clean/surveys/count_matrix_raw.rds"
)

saveRDS(
  all_species,
  "data_clean/metadata/species_list.rds"
)
