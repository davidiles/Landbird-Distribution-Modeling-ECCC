# ============================================================
# 02_load_format_point_counts.R
#
# Purpose
#   - Load and harmonize NatureCounts OBBA2/OBBA3 point counts + OBBA3 checklists
#   - Remove unusable surveys (missing date/time/lat/lon/species/count)
#   - Lump a handful of species_id discrepancies within/between atlases
#   - Build a survey-level sf object and a survey x species count matrix
#   - Enforce identical row order between the two before saving
#
# Inputs
#   data_clean/spatial/study_area.rds
#   <Data>/Bird_Data_Raw/OBBA3/Point_Count/... , OBBA2/Point_Count/...  (via NatureCounts)
#   R/functions/survey_processing_utils.R
#
# Outputs
#   data_clean/surveys/surveys_raw.rds
#   data_clean/surveys/count_matrix_raw.rds
#   data_clean/metadata/species_list.rds
# ============================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(sf)
  library(lubridate)
  library(naturecounts)
  library(purrr)
  library(readr)
  library(tidyr)
  library(here)
})

source(here::here("R", "00_config_paths.R"))
source(file.path(paths$functions, "survey_processing_utils.R"))

study_area <- readRDS(file.path(paths$data_clean, "spatial", "study_area.rds"))
crs_aea_km <- study_area$crs

# ------------------------------------------------------------
# Species list (NatureCounts canonical)
# ------------------------------------------------------------

all_species <- search_species_code() %>%
  rename(sp_code = BSCDATA, species_scientific_name = scientific_name)

english_lookup <- all_species %>%
  group_by(species_id) %>%
  summarise(english_name = dplyr::first(english_name), .groups = "drop")

# ------------------------------------------------------------
# Atlas 3 point counts
# ------------------------------------------------------------

NC_OBBA3_PC_raw <- nc_data_dl(collections = "ONATLAS3PC", username = "ilesd",
                              fields_set = "core", info = "Ontario Atlas") %>%
  select(where(~ sum(!is.na(.x)) > 0)) %>%
  mutate(
    Date_Time = make_datetime_from_frac_hours(
      ymd(ObservationDate), as.numeric(TimeObservationsStarted), tz = "UTC"
    ),
    Project_Name   = "OBBA3",
    Data_Source    = "NatureCounts",
    Special_Survey = if_else(EffortMeasurement1 == "Special", "Special", "No")
  )

NC_OBBA3_PC <- NC_OBBA3_PC_raw %>%
  mutate(across(c(DecimalLatitude, DecimalLongitude, DurationInHours,
                  species_id, ObservationCount), as.numeric)) %>%
  mutate(DecimalLongitude = if_else(DecimalLongitude > 0, -DecimalLongitude, DecimalLongitude)) %>%
  mutate(
    # Infer survey type from NatureCounts free-text/effort fields.
    Survey_Type = infer_survey_type_OBBA3_PC(
      Remarks, Remarks2, EffortMeasurement1, SurveyAreaIdentifier
    )
  ) %>%
  dplyr::relocate(
    Data_Source, Project_Name, Survey_Type, Special_Survey,
    DecimalLatitude, DecimalLongitude, Date_Time, DurationInHours,
    species_id, ObservationCount
  )

# ------------------------------------------------------------
# Atlas 3 checklists
# ------------------------------------------------------------

NC_OBBA3_CL_raw <- nc_data_dl(collections = "ONATLAS3BE_DO", username = "ilesd",
                              fields_set = "extended", info = "Ontario Atlas") %>%
  select(where(~ sum(!is.na(.x)) > 0)) %>%
  mutate(across(c(species_id, DecimalLatitude, DecimalLongitude,
                  TimeObservationsStarted, DurationInHours, ObservationCount), as.numeric)) %>%
  mutate(
    Survey_Type         = ProtocolType,
    Distance_Traveled_m = as.numeric(EffortMeasurement1),
    Date_Time           = make_datetime_from_frac_hours(
      ymd(ObservationDate), TimeObservationsStarted, tz = "UTC"
    ),
    Project_Name   = "OBBA3",
    Data_Source    = "NatureCounts",
    Special_Survey = if_else(EffortMeasurement1 == "Special", "Special", "No")
  ) %>%
  # Drop mislabeled effort units.
  filter(EffortUnits1 == "distance_meters")

NC_OBBA3_CL <- NC_OBBA3_CL_raw %>%
  mutate(
    Survey_Type = infer_survey_type_OBBA3_CL(
      Survey_Type, Remarks, Remarks2, EffortMeasurement1, SurveyAreaIdentifier
    )
  ) %>%
  dplyr::relocate(
    Data_Source, Project_Name, Survey_Type,
    DecimalLatitude, DecimalLongitude, Date_Time, DurationInHours,
    species_id, ObservationCount
  ) %>%
  mutate(TimeObservationsStarted = as.character(TimeObservationsStarted))

# ------------------------------------------------------------
# Atlas 2 point counts
# ------------------------------------------------------------

NC_OBBA2_PC_raw <- nc_data_dl(collections = "OBBA2PC", username = "ilesd",
                              fields_set = "extended", info = "Ontario Atlas") %>%
  select(where(~ sum(!is.na(.x)) > 0)) %>%
  rename(TimeObservationsStarted = TimeCollected) %>%
  mutate(
    Date_Time = make_datetime_from_frac_hours(
      ymd(paste(YearCollected, MonthCollected, DayCollected, sep = "-")),
      as.numeric(TimeObservationsStarted), tz = "UTC"
    ),
    Project_Name   = "OBBA2",
    Data_Source    = "NatureCounts",
    Survey_Type    = "Point_Count",
    Special_Survey = "No"
  )

NC_OBBA2_PC <- NC_OBBA2_PC_raw %>%
  mutate(across(c(DecimalLatitude, DecimalLongitude, DurationInHours,
                  species_id, ObservationCount), as.numeric)) %>%
  dplyr::relocate(
    Data_Source, Project_Name, Survey_Type, Special_Survey,
    DecimalLatitude, DecimalLongitude, Date_Time, DurationInHours,
    species_id, ObservationCount
  )

# ------------------------------------------------------------
# Combine atlases and basic cleaning
# ------------------------------------------------------------

# NOTE: the longitude filter removes obvious data errors (e.g. missing sign).
# A stricter within-boundary spatial filter happens later in the pipeline.
NC_long <- bind_rows(NC_OBBA3_PC, NC_OBBA2_PC, NC_OBBA3_CL) %>%
  filter(
    !is.na(DecimalLatitude), !is.na(DecimalLongitude),
    DecimalLongitude < -60, DecimalLatitude < 60,
    !is.na(species_id), !is.na(ObservationCount), !is.na(Date_Time)
  ) %>%
  rename(Latitude = DecimalLatitude, Longitude = DecimalLongitude) %>%
  mutate(
    # Enforce stable ID construction across reruns/adaptations.
    Latitude  = round(Latitude, 6),
    Longitude = round(Longitude, 6),
    Survey_Duration_Minutes = round(DurationInHours * 60),
    Max_Distance_Metres     = Inf,
    survey_id = make_survey_id(Project_Name, Latitude, Longitude, Date_Time, Survey_Type)
  )

# ------------------------------------------------------------
# Fix species_id discrepancies within/between atlases (lump synonyms/subspecies)
# ------------------------------------------------------------

# Common (40919), Hoary (20400) and Redpolls (45264) -> lump to 45264.
NC_long$species_id[NC_long$species_id %in% c(40919, 20400, 45264)] <- 45264
# Northern Flicker auratus/luteus (10470) -> Colaptes auratus (48798).
NC_long$species_id[NC_long$species_id %in% c(10470, 48798)] <- 48798
# Common/Wilson's Snipe (4940, 4950) -> 4950.
NC_long$species_id[NC_long$species_id %in% c(4940, 4950)] <- 4950
# Yellow-rumped Warbler coronata (16620) -> YRWA (16610).
NC_long$species_id[NC_long$species_id %in% c(16620, 16610)] <- 16610
# Red-winged Blackbird listed under two ids (40876 not in the NC species table).
NC_long$species_id[NC_long$species_id %in% c(40876, 19530)] <- 19530
# Winter Wren (14830, 40702) -> 40702.
NC_long$species_id[NC_long$species_id %in% c(14830, 40702)] <- 40702

# ------------------------------------------------------------
# Optional QA: flag species with large cross-atlas differences
# ------------------------------------------------------------
# Extreme OBBA2-vs-OBBA3 differences can indicate species_id mislabelling
# between atlases; this print-out informs the manual lumping above. It does not
# feed any saved output, so it can be commented out for non-interactive runs.

large_differences <- NC_long %>%
  left_join(unique(all_species[, c("species_id", "species_scientific_name")]), by = "species_id") %>%
  group_by(species_id, Project_Name) %>%
  summarize(sum_count = sum(ObservationCount, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Project_Name, values_from = sum_count, values_fill = 0) %>%
  mutate(
    total    = OBBA2 + OBBA3,
    log_fold = log((OBBA3 + 1) / (OBBA2 + 1))
  ) %>%
  left_join(unique(all_species[, c("species_id", "species_scientific_name")]), by = "species_id") %>%
  filter(total >= 100, abs(log_fold) >= abs(log(0.7))) %>%
  arrange(desc(total)) %>%
  left_join(english_lookup, by = "species_id")

print(as.data.frame(large_differences))

# ------------------------------------------------------------
# Survey-level table
# ------------------------------------------------------------

survey_info <- NC_long %>%
  group_by(survey_id) %>%
  summarise(
    Project_Name  = first(Project_Name),
    Data_Source   = first(Data_Source),
    Survey_Type   = first(Survey_Type),
    Special_Survey = first(Special_Survey),
    Latitude      = first(Latitude),
    Longitude     = first(Longitude),
    Date_Time     = first(Date_Time),
    DurationInHours = median(DurationInHours, na.rm = TRUE),
    Survey_Duration_Minutes = round(DurationInHours * 60),
    Distance_Traveled_m     = first(Distance_Traveled_m),
    NumberOfObservers       = first(NumberOfObservers),
    EffortMeasurement1 = first(EffortMeasurement1),
    EffortUnits1       = first(EffortUnits1),
    EffortMeasurement4 = first(EffortMeasurement4),
    EffortUnits4       = first(EffortUnits4),
    EffortMeasurement5 = first(EffortMeasurement5),
    EffortUnits5       = first(EffortUnits5),
    ObservationDescriptor  = first(ObservationDescriptor),
    ObservationCount2      = first(ObservationCount2),
    ObservationDescriptor2 = first(ObservationDescriptor2),
    ObservationCount3      = first(ObservationCount3),
    ObservationDescriptor3 = first(ObservationDescriptor3),
    ObservationDate        = first(ObservationDate),
    DateUncertaintyInDays  = first(DateUncertaintyInDays),
    AllIndividualsReported = first(AllIndividualsReported),
    AllSpeciesReport       = first(AllSpeciesReported),
    SurveyAreaIdentifier   = first(SurveyAreaIdentifier),
    Remarks  = first(Remarks),
    Remarks2 = first(Remarks2),
    .groups = "drop"
  ) %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE)

# ------------------------------------------------------------
# Survey x species count matrix
# ------------------------------------------------------------

count_matrix <- build_count_matrix(
  NC_long,
  survey_id_col = "survey_id",
  species_col   = "species_id",
  count_col     = "ObservationCount"
)

# ------------------------------------------------------------
# CRS transform
# ------------------------------------------------------------

survey_info <- survey_info %>%
  mutate(Atlas = Project_Name) %>%
  st_transform(crs_aea_km) %>%
  relocate(geometry, .after = last_col())

# ------------------------------------------------------------
# CRITICAL: enforce identical row order between survey_info and count_matrix
# ------------------------------------------------------------
# survey_info rows come from group_by(survey_id) %>% summarise() (dplyr's group
# ordering, C locale in dplyr >= 1.1), while count_matrix rownames come from
# build_count_matrix()'s sort(unique(survey_id)) (base sort() under LC_COLLATE).
# Those two ordering engines DISAGREE on a non-C locale, so the same set of
# surveys can be held in DIFFERENT row orders. Downstream (05/06) matches counts
# to surveys positionally by obs_idx, so a disagreement silently attaches the
# wrong counts. Re-key count_matrix into survey_info's order and assert it.

stopifnot(
  length(unique(NC_long$survey_id)) == nrow(survey_info),
  length(unique(rownames(count_matrix))) == nrow(survey_info),
  setequal(rownames(count_matrix), survey_info$survey_id),
  anyDuplicated(survey_info$survey_id)  == 0,
  anyDuplicated(rownames(count_matrix)) == 0
)
count_matrix <- count_matrix[survey_info$survey_id, , drop = FALSE]
stopifnot(identical(rownames(count_matrix), survey_info$survey_id))

# ------------------------------------------------------------
# Save outputs
# ------------------------------------------------------------

dir.create(file.path(paths$data_clean, "surveys"),  recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(paths$data_clean, "metadata"), recursive = TRUE, showWarnings = FALSE)

saveRDS(survey_info,  file.path(paths$data_clean, "surveys",  "surveys_raw.rds"))
saveRDS(count_matrix, file.path(paths$data_clean, "surveys",  "count_matrix_raw.rds"))
saveRDS(all_species,  file.path(paths$data_clean, "metadata", "species_list.rds"))

message("02_load_format_point_counts.R complete.")
