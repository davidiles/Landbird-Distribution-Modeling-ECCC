# ============================================================================
# ONTARIO BREEDING BIRD ATLAS DATA PREPARATION
# Process NatureCounts Data
# ============================================================================

# ------------------------------------------------
# Load/install packages
# ------------------------------------------------

my_packs <- c('tidyverse', 'magrittr', 'sf', 'terra', 'wildrtrax', 'ggspatial', 
              'naturecounts', 'suncalc', 'lubridate', 'openxlsx', 'readxl', 'lutz','purrr')

if (any(!my_packs %in% installed.packages()[, 'Package'])) {
  install.packages(my_packs[which(!my_packs %in% installed.packages()[, 'Package'])], dependencies = TRUE)
}
lapply(my_packs, require, character.only = TRUE)

rm(list = ls())

setwd("C:/Users/IlesD/OneDrive - EC-EC/Iles/Projects/Landbirds/Landbird-Distribution-Modeling-ECCC/")

`%!in%` <- Negate(`%in%`)

# ============================================================================
# Initial prep
# ============================================================================

# Projections and study area
AEA_proj <- st_crs(3978)

# Study area boundary (Ontario + 100km buffer)
Study_Area <- st_read("Data/Spatial/BCR/BCR_Terrestrial_master.shp") %>%
  subset(PROVINCE_S == "ONTARIO") %>%
  st_make_valid() %>%
  dplyr::select(BCR, PROVINCE_S) %>%
  st_union() %>%
  st_transform(AEA_proj)

# Get species list
all_species <- search_species_code() %>% 
  rename(sp_code = BSCDATA, species_scientific_name = scientific_name)

# ============================================================================
# PART 2: NATURECOUNTS DATA PROCESSING
# ============================================================================

# List of unique species codes
sp_id_unique <- unique(all_species$species_id)

# ---- Process OBBA-3 data

NCData_OBBA3 <- read.table("Data/Bird_Data_Raw/OBBA3/Point_Count/onatlas3pc_naturecounts_data.txt", 
                           sep = "\t", header = TRUE, fill = TRUE, comment.char = "", quote = "") %>%
  dplyr::select(where(~sum(!is.na(.x)) > 0)) %>%
  rename(species_id = SpeciesCode) %>%
  mutate(
    Date_Time = ymd(ObservationDate) + hours(floor(TimeObservationsStarted)) + 
      minutes(round(60 * (TimeObservationsStarted - floor(TimeObservationsStarted)))),
    Data_Source = "NatureCounts",
    Project_Name = "OBBA3",
    
    # Assume ARU survey if any of the following are true
    Survey_Type = case_when(
      str_detect(Remarks2, "ARU") ~ "ARU",
      str_detect(EffortMeasurement1, "ARU") ~ "ARU",
      str_detect(Remarks2, "WildtraxID") ~ "ARU",
      str_detect(SurveyAreaIdentifier, "Wildtrax") ~ "ARU",
      str_detect(Remarks, "WildtraxID") ~ "ARU",
      TRUE ~ "Point_Count"  # fallback
    )
  ) %>%
  dplyr::select(Data_Source, Project_Name, Survey_Type, DecimalLatitude, DecimalLongitude,
                Date_Time, DurationInHours, species_id, ObservationCount)


# ---- Process OBBA-2 data
NCData_OBBA2 <- read.table("Data/Bird_Data_Raw/OBBA2/Point_Count/obba2pc_naturecounts_data.txt", 
                           sep = "\t", header = TRUE, fill = TRUE, comment.char = "", quote = "") %>%
  dplyr::select(where(~sum(!is.na(.x)) > 0)) %>%
  left_join(all_species[, c("sp_code", "species_id")], by = c("SpeciesCode" = "sp_code"), multiple = "first") %>%
  rename(TimeObservationsStarted = TimeCollected) %>%
  mutate(
    Date_Time = ymd(paste(YearCollected, MonthCollected, DayCollected, sep = "-")) + 
      hours(floor(TimeObservationsStarted)) + 
      minutes(round(60 * (TimeObservationsStarted - floor(TimeObservationsStarted)))),
    Data_Source = "NatureCounts",
    Project_Name = "OBBA2",
    Survey_Type = "Point_Count") %>%
  dplyr::select(Data_Source, Project_Name, Survey_Type, DecimalLatitude, DecimalLongitude,
                Date_Time, DurationInHours, species_id, ObservationCount)

# ---- Combine into single data object

NCData <- bind_rows(NCData_OBBA3, NCData_OBBA2) %>%
  mutate(survey_ID = as.numeric(as.factor(paste(Project_Name, DecimalLatitude, DecimalLongitude, Date_Time)))) %>%
  subset(!is.na(DecimalLatitude) & !is.na(DecimalLongitude) & DecimalLongitude < -60 & 
           !is.na(species_id) & !is.na(ObservationCount)) %>%
  dplyr::rename(Latitude = DecimalLatitude,
                Longitude = DecimalLongitude) %>%
  mutate(Survey_Duration_Minutes = round(DurationInHours*60),
         Max_Distance_Metres = Inf) %>%
  subset(Survey_Duration_Minutes <= 10)

# Create survey info
NC_surveyinfo <- NCData %>%
  dplyr::select(-species_id, -ObservationCount) %>%
  distinct() %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE) %>%
  subset(!is.na(Date_Time))

# Create count matrix
NC_matrix <- matrix(0, nrow = nrow(NC_surveyinfo), ncol = length(sp_id_unique),
                    dimnames = list(NC_surveyinfo$survey_ID, sp_id_unique))

# Populate matrix 
count_summary <- NCData %>%
  group_by(survey_ID, species_id) %>%
  summarize(total_count = sum(ObservationCount), .groups = 'drop')

# Precompute index lookups
row_indices <- match(count_summary$survey_ID, rownames(NC_matrix))
col_indices <- match(count_summary$species_id, colnames(NC_matrix))

# Filter valid indices (non-NA)
valid <- !is.na(row_indices) & !is.na(col_indices)

# Use matrix indexing to assign values
NC_matrix[cbind(row_indices[valid], col_indices[valid])] <- count_summary$total_count[valid]

# Remove species never observed
NC_matrix <- NC_matrix[, -which(colSums(NC_matrix, na.rm = TRUE) == 0)]

# ============================================================================
# PART 5: CALCULATE TIME ZONES AND TIME SINCE SUNRISE FOR EACH SURVEY LOCATION
# ============================================================================

# Identify the timezone of each survey location
NC_surveyinfo_wgs84 <- st_transform(NC_surveyinfo, crs = 4326)

# Get time zones using point geometries
NC_surveyinfo$timezone <- lutz::tz_lookup(NC_surveyinfo_wgs84, method = "fast")

# Extract key data from sf object (drop geometry for processing)
survey_data <- NC_surveyinfo %>%
  st_drop_geometry() %>%
  select(Latitude, Longitude, Date_Time, timezone)

# Create unique combinations of lat/lon/date/timezone for efficient sunrise lookup
unique_sun_query <- survey_data %>%
  mutate(date = as.Date(Date_Time)) %>%
  distinct(Latitude, Longitude, date, timezone)

# Use suncalc to get local sunrise times (takes a while)
sunrise_times <- pmap_dfr(unique_sun_query, function(Latitude, Longitude, date, timezone) {
  suncalc::getSunlightTimes(
    date = date,
    lat = Latitude,
    lon = Longitude,
    keep = "sunrise",
    tz = timezone
  ) %>%
    select(date, lat, lon, sunrise) %>%
    rename(
      Latitude = lat,
      Longitude = lon
    ) %>%
    mutate(timezone = timezone)
})

# Join sunrise times back to full data
survey_data_with_sunrise <- survey_data %>%
  mutate(date = as.Date(Date_Time)) %>%
  left_join(sunrise_times, by = c("date", "Latitude", "Longitude", "timezone")) %>%
  mutate(
    # Ensure both are in the same time zone before subtraction
    Hours_Since_Sunrise = as.numeric(difftime(
      force_tz(Date_Time, tzone = timezone),
      force_tz(sunrise, tzone = timezone),
      units = "hours"
    ))
  ) %>%
  select(Hours_Since_Sunrise)

# Add Hours_Since_Sunrise back to sf object
NC_surveyinfo$Hours_Since_Sunrise <- survey_data_with_sunrise$Hours_Since_Sunrise

# Add atlas cycle
NC_surveyinfo$Atlas <- ifelse(year(NC_surveyinfo$Date_Time) > 2020, "OBBA3", "OBBA2")

# Remove erroneous columns
NC_surveyinfo <- NC_surveyinfo %>% dplyr::select(-DurationInHours) %>%
  relocate(survey_ID,Data_Source,Project_Name,Atlas,Survey_Type,Latitude,Longitude,Date_Time,timezone,Hours_Since_Sunrise) %>%
  st_transform(AEA_proj)

# ============================================================================
# PART 7: SAVE DATA
# ============================================================================

# Save final combined dataset
analysis_data <- list(
  NC_surveyinfo = NC_surveyinfo,
  full_count_matrix = NC_matrix,
  all_species = all_species,
  date_created = Sys.time()
)

saveRDS(analysis_data, file = "Analysis/Ontario_BBA/Data_Cleaned/analysis_data.rds")

# # ============================================================================
# # QUICK SUMMARY DATA PLOTS
# # ============================================================================
# 
# # Load Ontario boundary for plotting
# ONBoundary <- st_read("data/spatial/BCR/BCR_Terrestrial_master.shp") %>%
#   subset(PROVINCE_S == "ONTARIO") %>%
#   st_make_valid() %>%
#   dplyr::select(BCR, PROVINCE_S) %>%
#   st_transform(AEA_proj)
# 
# # ---- Summary of survey locations in each atlas
# ggplot() +
#   geom_sf(data = NC_surveyinfo, aes(col = Survey_Type), size = 0.3) +
#   geom_sf(data = ONBoundary, fill = "transparent", col = "gray60") +
#   facet_grid(Survey_Type ~ Atlas) +
#   theme_bw() +
#   labs(title = "Ontario Breeding Bird Atlas - Data Coverage",
#        subtitle = sprintf("%d surveys across %d species", nrow(NC_surveyinfo), ncol(NC_matrix)-1)) +
#   theme(axis.text = element_blank(), axis.ticks = element_blank())
# 
# # ---- Summary of surveys by day of year
# ggplot(NC_surveyinfo, aes(x = yday(Date_Time), fill = Survey_Type)) +
#   geom_histogram(binwidth = 1) +
#   facet_grid(Atlas~.) +
#   labs(title = "Survey Distribution Across Days of Year", x = "Day of Year", y = "Survey Count")+
#   theme_bw()
# 
# # ---- Summary of surveys by time of day
# ggplot(NC_surveyinfo, aes(x = hour(Date_Time) +
#                           minute(Date_Time) / 60 +
#                           second(Date_Time) / 3600, 
#                         fill = Survey_Type)) +
#   geom_histogram(binwidth = 0.1) +
#   facet_grid(Atlas~.) +
#   labs(title = "Survey Distribution Across Hours of the Day", x = "Hour (local time)", y = "Survey Count")+
#   theme_bw()
# 
# # ---- Summary of surveys by Hours Since Sunrise
# ggplot(NC_surveyinfo, aes(x = Hours_Since_Sunrise, fill = Survey_Type)) +
#   geom_histogram(binwidth = 0.1) +
#   facet_grid(Atlas~.) +
#   labs(title = "Surveys Relative to Sunrise", x = "Hours Since Sunrise", y = "Survey Count")+
#   theme_bw()
