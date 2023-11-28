# ---------------------------------------------------------------------------
# ***************************************************************************
# Download WildTrax datasets within a study area boundary
# ***************************************************************************
# ---------------------------------------------------------------------------

# ------------------------------------------------
# Load/install packages and set graphical themes / working directory
# ------------------------------------------------
my_packs = c('tidyverse',
             'openxlsx',
             'RColorBrewer',
             'viridis',
             'ggrepel',
             'scales',
             'wildRtrax',
             'lubridate','sf')

if (any(!my_packs %in% installed.packages()[, 'Package'])) {install.packages(my_packs[which(!my_packs %in% installed.packages()[, 'Package'])],dependencies = TRUE)}
lapply(my_packs, require, character.only = TRUE)

rm(list=ls())

#remotes::install_github("ABbiodiversity/wildRtrax@development")

# ------------------------------------------------
# Set working directory
# ------------------------------------------------

stub <- function() {}
thisPath <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  if (length(grep("^-f$", cmdArgs)) > 0) {
    # R console option
    normalizePath(dirname(cmdArgs[grep("^-f", cmdArgs) + 1]))[1]
  } else if (length(grep("^--file=", cmdArgs)) > 0) {
    # Rscript/R console option
    scriptPath <- normalizePath(dirname(sub("^--file=", "", cmdArgs[grep("^--file=", cmdArgs)])))[1]
  } else if (Sys.getenv("RSTUDIO") == "1") {
    # RStudio
    dirname(rstudioapi::getSourceEditorContext()$path)
  } else if (is.null(attr(stub, "srcref")) == FALSE) {
    # 'source'd via R console
    dirname(normalizePath(attr(attr(stub, "srcref"), "srcfile")$filename))
  } else {
    stop("Cannot find file path")
  }
}

dirname <- thisPath()
setwd(dirname)

# -----------------------------------------------
# Useful functions
# -----------------------------------------------

`%!in%` <- Negate(`%in%`)

# -----------------------------------------------
# Load list of species prepared by script 0_species_list.R
# -----------------------------------------------

all_species <- readRDS("../Data_Cleaned/all_species.RDS")

# -----------------------------------------------
# Polygon delineating study area boundary
# -----------------------------------------------

Study_Area <- st_read("../../../Data/Spatial/National/BCR/BCR_Terrestrial_master.shp")  %>%
  subset(PROVINCE_S == "ONTARIO") %>%
  st_make_valid() %>%
  dplyr::select(BCR, PROVINCE_S)

# -----------------------------------------------
# WildTrax credentials
# -----------------------------------------------

wt_auth() # Need your wildtrax username

# ---------------------------------------------------------
# Identify projects that are part of the BBMP
# ---------------------------------------------------------

ARU_projects <- wt_get_download_summary(sensor_id = 'ARU') %>% subset(sensor == "ARU") %>% arrange(project)
PC_projects <- wt_get_download_summary(sensor_id = 'PC') %>% subset(sensor == "PC") %>% arrange(project)

# # ---------------------------------------------------------
# # Download projects with 'ARU' data
# # ---------------------------------------------------------
# 
# ARU_fulldat <- data.frame()
# 
# for (i in 1:nrow(ARU_projects)){
#   print(i)
# 
#   loc <- wt_download_report(project_id = ARU_projects$project_id[i], sensor_id = "ARU", reports = c("location"), weather_cols = FALSE)
# 
#   if (is.null(nrow(loc))) next
# 
#   loc <- loc %>%
#     select(location_id,latitude,longitude) %>%
#     unique()
# 
#   rec <- wt_download_report(project_id = ARU_projects$project_id[i], sensor_id = "ARU", reports = c("recording"), weather_cols = FALSE) %>%
#     mutate(year = as.numeric(year(recording_date_time))) %>%
#     select(project_id,location_id,year) %>%
#     unique()
# 
#   loc <- wt_download_report(project_id = ARU_projects$project_id[i], sensor_id = "ARU", reports = c("location"), weather_cols = FALSE) %>%
#     select(location_id,latitude,longitude) %>%
#     unique()
# 
#   dat <- full_join(loc,rec)
# 
#   ARU_fulldat <- rbind(ARU_fulldat, dat)
# }
# 
# write.csv(ARU_fulldat, file = "../Data_Cleaned/WildTrax/WildTrax_ARU_locations.csv", row.names = FALSE)
# 
# # ---------------------------------------------------------
# # Download projects with 'PC' data
# # ---------------------------------------------------------
# 
# PC_fulldat <- data.frame()
# 
# for (i in 1:nrow(PC_projects)){
#   print(i)
# 
#   loc <- wt_download_report(project_id = PC_projects$project_id[i], sensor_id = "PC", reports = c("location"), weather_cols = FALSE)
# 
#   if (is.null(nrow(loc))) next
# 
#   loc <- loc %>%
#     select(location_id,latitude,longitude) %>%
#     unique()
# 
#   rec <- wt_download_report(project_id = PC_projects$project_id[i], sensor_id = "PC", reports = c("point_count"), weather_cols = FALSE) %>%
#     mutate(year = as.numeric(year(survey_date))) %>%
#     select(project_id,location_id,year) %>%
#     unique()
# 
#   loc <- wt_download_report(project_id = PC_projects$project_id[i], sensor_id = "PC", reports = c("location"), weather_cols = FALSE) %>%
#     select(location_id,latitude,longitude) %>%
#     unique()
# 
#   dat <- full_join(loc,rec)
# 
#   PC_fulldat <- rbind(PC_fulldat, dat)
# }
# 
# write.csv(PC_fulldat, file = "../Data_Cleaned/WildTrax/WildTrax_PC_locations.csv", row.names = FALSE)
# 
# ---------------------------------------------------------
# Subset to data within Study Area Boundary
# ---------------------------------------------------------

ARU_fulldat <- read.csv(file = "../Data_Cleaned/WildTrax/WildTrax_ARU_locations.csv")
PC_fulldat <- read.csv(file = "../Data_Cleaned/WildTrax/WildTrax_PC_locations.csv")

# Remove sites outside study area
ARU_sf <- ARU_fulldat %>%
  na.omit() %>%
  st_as_sf(coords=c("longitude","latitude"),crs=4326, remove = FALSE) %>%
  st_transform(crs = st_crs(Study_Area)) %>%
  st_intersection(Study_Area)

PC_sf <- PC_fulldat %>%
  na.omit() %>%
  st_as_sf(coords=c("longitude","latitude"),crs=4326, remove = FALSE) %>%
  st_transform(crs = st_crs(Study_Area))%>%
  st_intersection(Study_Area)

# Summarize number of surveys from each project
ARU_summary <- ARU_sf %>%
  as.data.frame() %>%
  group_by(project_id) %>%
  summarize(n_surveys = n()) %>%
  subset(n_surveys > 5) %>%
  left_join(ARU_projects[,c("project_id","organization","project")]) %>%
  arrange(desc(n_surveys)) %>%
  select(n_surveys,organization,project)

PC_summary <- PC_sf %>%
  as.data.frame() %>%
  group_by(project_id) %>%
  summarize(n_surveys = n()) %>%
  subset(n_surveys > 5) %>%
  left_join(PC_projects[,c("project_id","organization","project")]) %>%
  arrange(desc(n_surveys)) %>%
  select(n_surveys,organization,project)

# Save summary of projects within the study area boundary
write.csv(ARU_summary,file = "../Data_Cleaned/WildTrax/WildTrax_ARU_summary.csv",row.names = FALSE)
write.csv(PC_summary,file = "../Data_Cleaned/WildTrax/WildTrax_PC_summary.csv",row.names = FALSE)



# *********************************************************
# *********************************************************
# DOWNLOAD / PROCESS ARU DATA
# *********************************************************
# *********************************************************

# ---------------------------------------------------------
# Download species records and recording info from each ARU survey
# ---------------------------------------------------------

# Ideally would download reports = "main" but this is not working (times out)
ARU_tags <- data.frame()
ARU_recordings <- data.frame()
for (i in 1:nrow(ARU_summary)){
  project_name = ARU_summary$project[i]
  PID = subset(ARU_projects, project == project_name)$project_id
  recs <- wt_download_report(project_id = PID, sensor_id = "ARU", reports = c("recording"), weather_cols = FALSE)
  tags <- wt_download_report(project_id = PID, sensor_id = "ARU", reports = c("tag"), weather_cols = FALSE)
  ARU_recordings <- rbind(ARU_recordings,recs)
  ARU_tags <- rbind(ARU_tags,tags)
}

# ---------------------------------------------------------
# Remove individual recordings outside the study area boundary
# ---------------------------------------------------------

ARU_recordings <- left_join(ARU_recordings, unique(ARU_fulldat[,c("location_id","latitude","longitude")])) %>%
  subset(!is.na(latitude) & !is.na(longitude)) %>%
  st_as_sf(coords=c("longitude","latitude"),crs=4326, remove = FALSE) %>%
  st_transform(crs = st_crs(Study_Area)) %>%
  st_intersection(Study_Area)

# ---------------------------------------------------------
# Only consider transcribed recordings
# ---------------------------------------------------------

ARU_recordings <- subset(ARU_recordings, aru_task_status == "Transcribed")

# ---------------------------------------------------------
# Date of recordings
# ---------------------------------------------------------

ARU_recordings$recording_date_time <- lubridate::ymd_hms(ARU_recordings$recording_date_time)

# ---------------------------------------------------------
# Extract duration and transcription method (SPT or SPM) for each recording
# ---------------------------------------------------------

ARU_recordings <- ARU_recordings %>%
  mutate(Duration_Seconds = str_split_i(method," ",1) %>% str_sub(start = 1,end = -2) %>% as.numeric(),
         Transcription_Method = str_split_i(method," ",2) %>% str_sub(start = 2,end = -1)) 

# ---------------------------------------------------------
# Convert TMTT observations in ARU_tags dataframe
# ---------------------------------------------------------

ARU_tags$individual_count <- ARU_tags$abundance
ARU_tags <- wt_replace_tmtt(ARU_tags)

# ---------------------------------------------------------
# Fix abundance 
# ---------------------------------------------------------

ARU_tags$individual_count <- as.numeric(ARU_tags$individual_count)
ARU_tags <- subset(ARU_tags, !is.na(individual_count))

# ---------------------------------------------------------
# Convert WildTrax species codes to Birds Canada species codes when they disagree
# ---------------------------------------------------------

species <- subset(all_species,WT_spcd %in% unique(ARU_tags$species_code))

species_to_fix <- subset(species,WT_spcd != BSC_spcd)
for (i in 1:nrow(species_to_fix)){
  ARU_tags$species_code[which(ARU_tags$species_code == species_to_fix$WT_spcd[i])] <- species_to_fix$BSC_spcd[i]
}

# Remove species with no BSC code
subset(ARU_tags, species_code %!in% all_species$BSC_spcd)$species_common_name %>% table() %>% sort(decreasing = TRUE)

ARU_tags <- subset(ARU_tags, species_code %in% all_species$BSC_spcd)

# *********************************************************
# *********************************************************
# PROCESS SPT transcriptions
# *********************************************************
# *********************************************************

ARU_recordings_SPT <- subset(ARU_recordings, Transcription_Method == "SPT")

# ---------------------------------------------------------
# Remove duplicate rows (caused by multiple observers - records will be averaged)
# ---------------------------------------------------------

ARU_recordings_SPT <- dplyr::select(ARU_recordings_SPT,
                                    recording_id,latitude,longitude,
                                    Transcription_Method,Duration_Seconds,
                                    recording_date_time,equipment,
                                    organization,project,project_id,location,location_id
) %>% unique()

# Duplicated recording_ids
duplicated_ids <- ARU_recordings_SPT$recording_id[duplicated(ARU_recordings_SPT$recording_id)] %>% unique()
subset(ARU_recordings_SPT,recording_id %in% duplicated_ids) %>% arrange(recording_id)

# Remove duplicated recording IDs 
# (!!! THIS NEEDS TO BE CAREFULLY INSPECTED IN THE FUTURE  !!!)
duplicated_rows <- which(duplicated(ARU_recordings_SPT$recording_id))
ARU_recordings_SPT <- ARU_recordings_SPT[-duplicated_rows,]

# ---------------------------------------------------------
# Process tags (count individuals)
# ---------------------------------------------------------

# Extract correct tags (based on recording id)
ARU_tags_SPT <- subset(ARU_tags, recording_id %in% ARU_recordings_SPT$recording_id)
ARU_tags_SPT$recording_id <- factor(ARU_tags_SPT$recording_id,levels = unique(ARU_recordings_SPT$recording_id))

# Total number of individuals detected per survey (sum of individual counts)
ARU_counts_SPT <- ARU_tags_SPT %>%
  
  # Task ID is different observers
  group_by(species_code,recording_id,task_id) %>%
  summarize(total_count = sum(individual_count)) %>%
  
  # Take mean if there are multiple tasks
  group_by(species_code,recording_id) %>%
  summarize(mean_count = round(mean(total_count))) %>%
  
  as.data.frame() %>%
  pivot_wider(names_from = species_code,
              values_from = mean_count,
              values_fill = 0,
              id_expand = TRUE,
              names_expand = TRUE) %>%
  arrange(recording_id)

# Same ordering as ARU_recordings_SPT
mean(ARU_counts_SPT$recording_id == ARU_recordings_SPT$recording_id) # should be 1

# *********************************************************
# *********************************************************
# PROCESS SPM TRANSCRIPTIONS (individuals can appear more than once)
# *********************************************************
# *********************************************************

ARU_recordings_SPM <- subset(ARU_recordings, Transcription_Method == "SPM")

# ---------------------------------------------------------
# Remove duplicate rows (caused by multiple observers - records will be averaged)
# ---------------------------------------------------------

ARU_recordings_SPM <- dplyr::select(ARU_recordings_SPM,
                                    recording_id,latitude,longitude,
                                    Transcription_Method,Duration_Seconds,
                                    recording_date_time,equipment,
                                    organization,project,project_id,location,location_id
) %>% unique()

# In cases where a recording was transcribed using multiple SPM methods (occurs once), only use the first one
ARU_recordings_SPM <- ARU_recordings_SPM[!duplicated(ARU_recordings_SPM$recording_id),]

# ---------------------------------------------------------
# Process tags (count individuals)
# ---------------------------------------------------------

# Extract correct tags (based on recording id)
ARU_tags_SPM <- subset(ARU_tags, recording_id %in% ARU_recordings_SPM$recording_id)

ARU_tags_SPM$recording_id <- factor(ARU_tags_SPM$recording_id,
                                    levels = ARU_recordings_SPM$recording_id)


# Total number of individuals detected per survey (sum of individual counts)
ARU_counts_SPM <- ARU_tags_SPM %>%
  
  # Task ID is different observers
  group_by(species_code,recording_id,task_id,individual_order) %>%
  
  # Remove repeated observations of individual birds (individual order)
  summarize(individual_count = round(mean(individual_count))) %>%
  
  # Task ID is different observers
  group_by(species_code,recording_id,task_id) %>%
  summarize(total_count = sum(individual_count)) %>%
  
  # Take mean if there are multiple tasks
  group_by(species_code,recording_id) %>%
  summarize(mean_count = round(mean(total_count))) %>%
  
  as.data.frame() %>%
  pivot_wider(names_from = species_code,
              values_from = mean_count,
              values_fill = 0,
              id_expand = TRUE,
              names_expand = TRUE) %>%
  arrange(recording_id)

# Same ordering as ARU_recordings_SPM
mean(ARU_counts_SPM$recording_id == ARU_recordings_SPM$recording_id) # should be 1


# *********************************************************
# *********************************************************
# PROCESS POINT COUNT DATA
# *********************************************************
# *********************************************************

PC_summary  <- read.csv(file = "../Data_Cleaned/WildTrax/WildTrax_PC_summary.csv")

# ---------------------------------------------------------
# Download species records and recording info from each Point Count survey
# ---------------------------------------------------------

# Note that all location information is stored in point count dataframe
PC_counts <- data.frame()
for (i in 1:nrow(PC_summary)){
  project_name = PC_summary$project[i]
  PID = subset(PC_projects, project == project_name)$project_id
  counts <- wt_download_report(project_id = PID, sensor_id = "PC", reports = c("main"), weather_cols = FALSE)
  PC_counts <- bind_rows(PC_counts,counts)
}

# ---------------------------------------------------------
# Convert WildTrax species codes to Birds Canada species codes when they disagree
# ---------------------------------------------------------

species <- subset(all_species,WT_spcd %in% unique(PC_counts$species_code))

species_to_fix <- subset(species,WT_spcd != BSC_spcd)
for (i in 1:nrow(species_to_fix)){
  PC_counts$species_code[which(PC_counts$species_code == species_to_fix$WT_spcd[i])] <- species_to_fix$BSC_spcd[i]
}

# Remove species with no BSC code
subset(PC_counts, species_code %!in% all_species$BSC_spcd)$species_common_name %>% table() %>% sort(decreasing = TRUE)

PC_counts <- subset(PC_counts, species_code %in% all_species$BSC_spcd)


# ---------------------------------------------------------
# Dataframe to track information associated with each survey
# ---------------------------------------------------------

# 1 row per survey
PC_surveys <- PC_counts %>%
  select(-detection_distance,-detection_time,-species_code,-species_common_name,-species_scientific_name,-individual_count,-detection_heard,-detection_seen,-detection_comments) %>%
  unique()

# ---------------------------------------------------------
# Date of each survey
# ---------------------------------------------------------

PC_surveys$survey_date <- lubridate::ymd_hms(PC_surveys$survey_date)

# ---------------------------------------------------------
# Append total duration of each survey
# ---------------------------------------------------------

table(PC_surveys$survey_duration_method)
method_definitions <- rbind(c(survey_duration_method = "0-1-2-3-4-5-6-7-8-9-10min", Survey_Duration_Minutes = 10),
                            c(survey_duration_method = "0-10min", Survey_Duration_Minutes = 10),
                            c(survey_duration_method = "0-3-10min", Survey_Duration_Minutes = 10),
                            c(survey_duration_method = "0-3-5-10min", Survey_Duration_Minutes = 10),
                            c(survey_duration_method = "0-3-5min", Survey_Duration_Minutes = 5),
                            c(survey_duration_method = "0-3min", Survey_Duration_Minutes = 3),
                            c(survey_duration_method = "0-5-10min", Survey_Duration_Minutes = 5)) %>%
  as.data.frame()

PC_surveys <- left_join(PC_surveys, method_definitions)
table(PC_surveys$Survey_Duration_Minutes, useNA = "always")


# ---------------------------------------------------------
# Append maximum distance of each survey
# ---------------------------------------------------------

table(PC_surveys$survey_distance_method)
method_definitions <- rbind(c(survey_distance_method = "0m-50m-100m", Max_Distance_Metres = 100),
                            c(survey_distance_method = "0m-50m-100m-150m-INF", Max_Distance_Metres = Inf),
                            c(survey_distance_method = "0m-50m-100m-INF", Max_Distance_Metres = Inf),
                            c(survey_distance_method = "0m-INF", Max_Distance_Metres = Inf)) %>%
  as.data.frame()

PC_surveys <- left_join(PC_surveys, method_definitions)

# ---------------------------------------------------------
# Remove records outside the study area boundary
# ---------------------------------------------------------

PC_surveys <- PC_surveys %>%
  subset(!is.na(latitude) & !is.na(longitude)) %>%
  st_as_sf(coords=c("longitude","latitude"),crs=4326, remove = FALSE) %>%
  st_transform(crs = st_crs(Study_Area)) %>%
  st_intersection(Study_Area)

# ---------------------------------------------------------
# Create a matrix of counts for each species
# ---------------------------------------------------------

PC_counts <- subset(PC_counts, survey_id %in% PC_surveys$survey_id)
PC_counts$survey_id <- factor(PC_counts$survey_id, levels = PC_surveys$survey_id)

# Total number of individuals detected per survey (sum of individual counts)
PC_counts <- PC_counts %>%
  
  group_by(species_code,survey_id) %>%
  summarize(total_count = sum(individual_count)) %>%
  
  as.data.frame() %>%
  pivot_wider(names_from = species_code,
              values_from = total_count,
              values_fill = 0,
              id_expand = TRUE,
              names_expand = TRUE) %>%
  arrange(survey_id)

# Same ordering as PC_surveys
mean(PC_counts$survey_id == PC_surveys$survey_id) # should be 1



# *********************************************************
# *********************************************************
# COMBINE ARU AND POINT COUNTS INTO SINGLE DATAFRAME
# RENAME AND SELECT RELEVANT COLUMNS
# *********************************************************
# *********************************************************

# Combine SPT and SPM into single dataframe
ARU_recordings_combined <- bind_rows(ARU_recordings_SPT,ARU_recordings_SPM) %>%
  
  rename(Project_Name = project, 
         survey_ID = recording_id, 
         Date_Time = recording_date_time, 
         Latitude = latitude, 
         Longitude = longitude) %>%
  
  mutate(Data_Source = "WildTrax",
         Survey_Type = paste0("ARU_",Transcription_Method),
         Survey_Duration_Minutes = Duration_Seconds/60,
         Max_Distance_Metres = Inf) %>%
  
  dplyr::select(Data_Source,
                Project_Name,
                Survey_Type,
                survey_ID,
                Latitude,
                Longitude,
                Date_Time,
                Survey_Duration_Minutes,
                Max_Distance_Metres)
  
# Point Counts
PC_surveys <- PC_surveys %>%
  
  rename(Project_Name = project, 
         survey_ID = survey_id, 
         Date_Time = survey_date, 
         Latitude = latitude, 
         Longitude = longitude) %>%
  
  mutate(Data_Source = "WildTrax",
         Survey_Type = "Point_Count",
         Survey_Duration_Minutes = as.numeric(Survey_Duration_Minutes),
         Max_Distance_Metres = as.numeric(Max_Distance_Metres)) %>%
  
  dplyr::select(Data_Source,
                Project_Name,
                Survey_Type,
                survey_ID,
                Latitude,
                Longitude,
                Date_Time,
                Survey_Duration_Minutes,
                Max_Distance_Metres)
  
# Combined
WT_surveyinfo <- bind_rows(ARU_recordings_combined, PC_surveys)


# *********************************************************
# *********************************************************
# Count matrix (combined)
# *********************************************************
# *********************************************************

# All species in WT database
WT_species <- wildRtrax::wt_get_species() %>% subset(species_class == "AVES")
WT_species_codes <- WT_species$species_code

# Fill in matrix
WT_matrix <- matrix(0,nrow=nrow(WT_surveyinfo), ncol = length(WT_species$species_code),
                            dimnames = list(NULL,WT_species$species_code))

for (spp in WT_species_codes){
  
  if (spp %in% colnames(ARU_counts_SPT)) WT_matrix[which(WT_surveyinfo$Survey_Type == "ARU_SPT"),which(colnames(WT_matrix) == spp)] <- as.data.frame(ARU_counts_SPT)[,spp]
  if (spp %in% colnames(ARU_counts_SPM)) WT_matrix[which(WT_surveyinfo$Survey_Type == "ARU_SPM"),which(colnames(WT_matrix) == spp)] <- as.data.frame(ARU_counts_SPM)[,spp]
  if (spp %in% colnames(PC_counts)) WT_matrix[which(WT_surveyinfo$Survey_Type == "Point_Count"),which(colnames(WT_matrix) == spp)] <- as.data.frame(PC_counts)[,spp]
  
}

# Remove species that were never detected
WT_matrix <- WT_matrix[,-which(colSums(WT_matrix)==0)]

# Fix particular species manually
#WT_matrix[,"YRWA"] <- WT_matrix[,"YRWA"] + WT_matrix[,"MYWA"]
#WT_matrix <- WT_matrix[,-which(colnames(WT_matrix) == "MYWA")]

# *********************************************************
# *********************************************************
# Save file
# *********************************************************
# *********************************************************

WT_dat <- list(WT_surveyinfo = WT_surveyinfo,
               WT_matrix = WT_matrix)
saveRDS(WT_dat, file = "../Data_Cleaned/WildTrax/WT_dat.rds")
