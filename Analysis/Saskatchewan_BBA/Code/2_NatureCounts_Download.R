# ---------------------------------------------------------------------------
# ***************************************************************************
# Download NatureCounts datasets within a study area boundary
#  - merge with WildTrax data, and remove duplicated records
# ***************************************************************************
# ---------------------------------------------------------------------------

# ------------------------------------------------
# Load/install packages and set graphical themes / working directory
# ------------------------------------------------

my_packs = c('tidyverse',
             'magrittr',
             'sf',
             'terra',
             'ggspatial',
             'naturecounts',
             'suntools',
             "readxl")

if (any(!my_packs %in% installed.packages()[, 'Package'])) {install.packages(my_packs[which(!my_packs %in% installed.packages()[, 'Package'])],dependencies = TRUE)}
lapply(my_packs, require, character.only = TRUE)

rm(list=ls())

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

# ------------------------------------------------
# ggplot theme
# ------------------------------------------------

CustomTheme <- theme_set(theme_bw())
CustomTheme <- theme_update(legend.key = element_rect(colour = NA), 
                            legend.key.height = unit(1.2, "line"),
                            panel.grid.major = element_line(colour = 'gray95'),
                            #panel.grid.minor = element_line(colour = 'transparent'),
                            panel.border = element_rect(linetype = "solid",
                                                        colour = "black",
                                                        linewidth = 1, fill = NA),
                            axis.line = element_line(colour = "black"),
                            strip.text = element_text(size = 12, colour = "black"),
                            strip.background = element_rect(colour = "black",
                                                            fill = "lightblue2",
                                                            linetype = "solid"),
                            axis.title.y = element_text(margin = margin(0,10,0,0)),
                            axis.title.x = element_text(margin = margin(10,0,0,0)),
                            panel.background = element_rect(fill = "white"))


# -----------------------------------------------
# Useful functions
# -----------------------------------------------

`%!in%` <- Negate(`%in%`)

# ************************************************************
# ************************************************************
# Prepare data - based on code from Dean Evans
# ************************************************************
# ************************************************************

###--------> Last downloaded/run on Nov 16 2023

# # Get Species information #
# species <- search_species_code()
# 
# duplicates <- species %>% 
#   group_by(species_id) %>% 
#   filter(n()>1)
# 
# species %<>% filter(BSCDATA %!in% c("MCLO","EWPE","CAGO","UNDD","GRPA","ASTK","MONP","COSN","RWBU","PEPL","UNSG","GUSP","RODO","ECDO","WPWI","TTWO","GRAJ","SIFL","HYBG","YWAR","NSTS","OTHR","DOWS","EWWR","UNWX","COMO","SBDU","JUNI","RIPH","BNOW","BDOW"))
# 
# ### Point Counts
# PCData <- read_xlsx(path = "../../../Data/Bird_Data_Raw/NatureCounts/SK_atlas/skatlas1pc_naturecounts_data.xlsx")
# 
# ### Checklists
# DOData <- read_xlsx(path = "../../../Data/Bird_Data_Raw/NatureCounts/SK_atlas/skatlas1be_do_naturecounts_data.xlsx")
# 
# ### Highest Breeding evidence
# BEData <- read_xlsx(path = "../../../Data/Bird_Data_Raw/NatureCounts/SK_atlas/skatlas1be_summ_naturecounts_data.xlsx")
# 
# #### Data Processing ####
# ### Remove empty columns
# PCData %<>% dplyr::select(where(~sum(!is.na(.x)) > 0))
# DOData %<>% dplyr::select(where(~sum(!is.na(.x)) > 0))
# BEData %<>% dplyr::select(where(~sum(!is.na(.x)) > 0))
# 
# ### Add Species codes
# PCData <- PCData %>% rename(species_id = SpeciesCode)
# DOData <- DOData %>% rename(species_id = SpeciesCode)
# BEData <- BEData %>% rename(species_id = SpeciesCode)
# 
# PCData %<>% left_join(species %>% dplyr::select(species_id,BSCDATA)) %>% rename(spcd=BSCDATA)
# DOData %<>% left_join(species %>% dplyr::select(species_id,BSCDATA)) %>% rename(spcd=BSCDATA)
# BEData %<>% left_join(species %>% dplyr::select(species_id,BSCDATA)) %>% rename(spcd=BSCDATA)
# 
# ###Check for missing spcd
# missing_spcd <- PCData %>% filter(is.na(spcd)) %>% dplyr::select(species_id,CommonName) %>% distinct()
# PCData %<>% mutate(spcd=case_when(species_id==32671 ~ "GHOW",
#                                   species_id==40183 ~ "RNEP",
#                                   species_id==46051 ~ "EUCD",
#                                   .default = spcd)) %>% filter(!is.na(species_id)) %>% filter(species_id!=12285)
# 
# missing_spcd <- DOData %>% filter(is.na(spcd)) %>% dplyr::select(species_id,CommonName) %>% distinct()
# DOData %<>% mutate(spcd=case_when(species_id==32671 ~ "GHOW",
#                                   species_id==40183 ~ "RNEP",
#                                   species_id==46051 ~ "EUCD",
#                                   species_id==40309 ~ "RTHA",
#                                   species_id==10485 ~ "NOFL",
#                                   species_id==40320 ~ "MERL",
#                                   species_id==45168 ~ "DEJU",
#                                   species_id==16563 ~ "YEWA",
#                                   species_id==2551 ~ "GBHE",
#                                   species_id==40667 ~ "BARS",
#                                   species_id==16801 ~ "PAWA",
#                                   species_id==40845 ~ "FOSP",
#                                   species_id==40182 ~ "COME",
#                                   .default = spcd)) %>% filter(!is.na(species_id)) %>% filter(species_id!=12285)
# 
# missing_spcd <- BEData %>% filter(is.na(spcd)) %>% dplyr::select(species_id,CommonName) %>% distinct()
# BEData %<>% mutate(spcd=case_when(species_id==32671 ~ "GHOW",
#                                   species_id==40183 ~ "RNEP",
#                                   species_id==46051 ~ "EUCD",
#                                   .default = spcd)) %>% filter(!is.na(species_id)) %>% filter(species_id!=12285)
# rm(missing_spcd)
# 
# # Remove breeding evidence of x or NA from all
# DOData %<>% filter(BreedingBirdAtlasCode!="X",!is.na(BreedingBirdAtlasCode))
# BEData %<>% filter(BreedingBirdAtlasCode!="X",!is.na(BreedingBirdAtlasCode))
# 
# saveRDS(PCData,"../Data_Cleaned/NatureCounts/PCData.rds")
# saveRDS(DOData,"../Data_Cleaned/NatureCounts/DOData.rds")


# ************************************************************
# ************************************************************
# Merge data with WildTrax
# ************************************************************
# ************************************************************

# -----------------------------------------------
# Data from WildTrax (prepared by previous script)
# -----------------------------------------------

AEA_proj <- "+proj=aea +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-106 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs "

WT_dat <- readRDS(file = "../Data_Cleaned/WildTrax/WT_dat.rds")
WT_surveyinfo <- WT_dat$WT_surveyinfo %>% st_transform(crs = AEA_proj)
WT_matrix <- WT_dat$WT_matrix

# -----------------------------------------------
# Load list of species prepared by script 0_species_list.R
# -----------------------------------------------

all_species <- readRDS("../Data_Cleaned/all_species.RDS")
BSC_species <- readRDS("../Data_Cleaned/BSC_species.RDS")

# ************************************************************
# ************************************************************
# PART 1: PREPARE POINT COUNT DATA
#
# This script creates two objects:
#   1) PC_surveyinfo - an sf object with "survey information" for each point count
#                    - e.g., Date, Location, Observer, Point count duration, etc
#   2) PC_matrix - a matrix object with counts for each species for each point count
#  
#   Note that PC_surveyinfo and PC_matrix should have same number of rows
#
# ************************************************************
# ************************************************************

# Prepared by a separate script; downloaded from NatureCounts
PCData <- readRDS("../Data_Cleaned/NatureCounts/PCData.rds")

# Create a "surveyID" column (unique to each point count)
PCData %<>% mutate(surveyID = as.numeric(factor(paste0(DecimalLatitude,DecimalLongitude,ObservationDate,TimeCollected))))

# Duration in minutes for each point count
PCData  %<>% mutate(DurationInMinutes=round(PCData$DurationInHours * 60))

# ---------------------------------------------------------
# QA/QC
# ---------------------------------------------------------

# Remove NA counts
PCData %<>%  filter(!is.na(ObservationCount)) 

# ---------------------------------------------------------
# Dataframe with one row per survey
# ---------------------------------------------------------

# Construct a dataset that contains JUST survey information (not counts)
PC_surveyinfo <- PCData %>% 
  as.data.frame() %>%
  dplyr::select(surveyID,Locality,DecimalLatitude,DecimalLongitude,ObservationDate,TimeCollected,DurationInMinutes,EffortMeasurement1) %>%
  distinct() %>%
  rename(sq_id = Locality) %>%
  st_as_sf(.,coords = c("DecimalLongitude","DecimalLatitude"), crs = st_crs(4326), remove = FALSE)

# Counts for each species on each survey are stored in a matrix
PC_matrix <- matrix(0,nrow = nrow(PC_surveyinfo),ncol = nrow(BSC_species),
                    dimnames = list(PC_surveyinfo$surveyID,BSC_species$BSC_spcd))

pb = txtProgressBar(min = 0, max = nrow(PCData), initial = 0, style = 3) 
for (i in 1:nrow(PCData)){
  PC_matrix[which(rownames(PC_matrix) == PCData$surveyID[i]),which(colnames(PC_matrix) == PCData$spcd[i])] <- PCData$ObservationCount[i]
  setTxtProgressBar(pb,i)
}
close(pb)

# Remove columns (species) that were never observed to save storage space
PC_matrix <- PC_matrix[,-which(colSums(PC_matrix,na.rm = T)==0)]

# ---------------------------------------------------------
# Ensure column names are consistent WildTrax dataset
# ---------------------------------------------------------

NC_PC_surveyinfo <- PC_surveyinfo %>%
  
  rename(survey_ID = surveyID, 
         Latitude = DecimalLatitude, 
         Longitude = DecimalLongitude,
         Survey_Duration_Minutes = DurationInMinutes,
         Survey_Type = EffortMeasurement1) %>%
  
  mutate(Max_Distance_Metres = Inf, 
         Date_Time = ymd(ObservationDate) + hours(floor(TimeCollected)) + minutes(round(60*(TimeCollected-floor(TimeCollected)))),
         Data_Source = "NatureCounts",
         Project_Name = "SK Breeding Bird Atlas") %>%
  
  dplyr::select(Data_Source,
                Project_Name,
                Survey_Type,
                survey_ID,
                Latitude,
                Longitude,
                Date_Time,
                Survey_Duration_Minutes,
                Max_Distance_Metres)

# ************************************************************
# ************************************************************
# PART 2: PREPARE CHECKLIST DATA
#
# This script creats two objects:
#   1) DO_surveyinfo - an sf object with "survey information" for each checklist
#                    - e.g., Date, Location, Observer, Checklist effort, etc
#   2) DO_matrix - a matrix object with counts for each species on each checklist
#  
#   Note that DO_surveyinfo and DO_matrix should have same number of rows
#
# ************************************************************
# ************************************************************

# Prepared by a separate script; downloaded from NatureCounts
DOData <- readRDS("../Data_Cleaned/NatureCounts/DOData.rds")

# ------------------------------------------------------------
# QA/QC
# ------------------------------------------------------------

# Only include complete checklists (AllSpeciesReported == "Y")
DOData %<>% mutate(AllSpeciesReported=case_when(AllSpeciesReported=="y" ~ "Y", .default = AllSpeciesReported))
table(DOData$AllSpeciesReported, useNA = "always")
DOData %<>% filter(AllSpeciesReported == "Y")

# ----------------------------------------------------------------
# ASSUME NA'S CORRESPOND TO NON-DETECTIONS (SET TO 0)
# ----------------------------------------------------------------

DOData$ObservationCount[which(is.na(DOData$ObservationCount))] <- 0

# ---------------------------------------------------------
# Dataframe with one row per checklist
# ---------------------------------------------------------

# Summarize remaining checklist information
DO_surveyinfo <- DOData %>% 
  dplyr::rename(TravelDistance_m = EffortMeasurement1,
                IncludesPointCounts = EffortMeasurement3) %>%
  group_by(SamplingEventIdentifier,
           Locality,
           ObservationDate,
           TimeCollected,
           CollectorNumber,
           ProtocolType,
           AllIndividualsReported,
           AllSpeciesReported,
           
           # Effort covariates
           DurationInHours,
           NumberOfObservers,
           
           # Travel distance
           TravelDistance_m,
           
           # Does checklist include point counts (EffortMeasurement3 == 1)
           IncludesPointCounts) %>%
  
  summarize(n_Species_Detected = length(unique(spcd)),
            n_Birds_Detected = sum(ObservationCount),
            
            # There is one checklist with a minor discrepancy in lat/lon coordinates among rows, so
            # calculate the mean lat/lon.  All other checklists are consistent among rows
            Latitude = mean(DecimalLatitude),
            Longitude = mean(DecimalLongitude)) %>%
  
  mutate(Date_Time = ymd(ObservationDate) + hours(floor(TimeCollected)) + minutes(round(60*(TimeCollected-floor(TimeCollected)))),
         Data_Source = "NatureCounts",
         Project_Name = "SK Breeding Bird Atlas",
         Survey_Duration_Minutes = DurationInHours*60) %>%
  
  rename(Atlas_Square_ID = Locality,
         survey_ID = SamplingEventIdentifier,
         Travel_Distance_Metres = TravelDistance_m,
         Survey_Type = ProtocolType) %>%
  ungroup() %>%
  st_as_sf(.,coords = c("Longitude","Latitude"), crs = st_crs(4326), remove = FALSE) %>%
  subset(!is.na(Date_Time)) %>% 
  dplyr::select(Data_Source,Project_Name,Survey_Type,survey_ID,Latitude,Longitude,
                Date_Time,Survey_Duration_Minutes,Travel_Distance_Metres,n_Species_Detected,n_Birds_Detected,IncludesPointCounts,geometry)

# ---------------------------------------------------------
# Create matrix of species counts (each row corresponds to one in DO_surveyinfo)
# ---------------------------------------------------------

DO_matrix <- matrix(0,nrow = nrow(DO_surveyinfo),ncol = nrow(BSC_species),
                    dimnames = list(DO_surveyinfo$survey_ID,BSC_species$BSC_spcd))

pb = txtProgressBar(min = 0, max = nrow(DOData), initial = 0, style = 3) 
for (i in 1:nrow(DOData)) {
  DO_matrix[which(rownames(DO_matrix) == DOData$SamplingEventIdentifier[i]),which(colnames(DO_matrix) == DOData$spcd[i])] <- DOData$ObservationCount[i]
  setTxtProgressBar(pb,i)
}
close(pb)

# Remove columns (species) that were never observed to save storage space
DO_matrix <- DO_matrix[,-which(colSums(DO_matrix,na.rm = TRUE)==0)]

# ************************************************************
# ************************************************************
# CHECK FOR OVERLAP BETWEEN NATURECOUNTS AND WILDTRAX POINT COUNT DATASETS
#
# - assume ARU data in WildTrax data is 'authoritative' because it can be reviewed/updated more often
#
# - assume Human Point counts in NatureCounts are authoritative because they have been thoroughly QA/QC'd
# ************************************************************
# ************************************************************

NC_PC_surveyinfo <- NC_PC_surveyinfo %>% st_transform(crs = AEA_proj)

# Identify NatureCounts point counts that are within 10 m of WildTrax locations
WT_buff <- st_buffer(WT_surveyinfo,10) %>% st_union()

NC_PC_surveyinfo$to_evaluate <- FALSE
NC_PC_surveyinfo$to_evaluate[which(as.matrix(st_intersects(NC_PC_surveyinfo,WT_buff)))] <- TRUE

distance_matrix <- st_distance(subset(NC_PC_surveyinfo, to_evaluate),WT_surveyinfo)

# Flags for which observations to remove
NC_PC_surveyinfo$Obs_Index <- 1:nrow(NC_PC_surveyinfo)
WT_surveyinfo$Obs_Index <- 1:nrow(WT_surveyinfo)
NC_PC_surveyinfo$to_remove <- FALSE
WT_surveyinfo$to_remove <- FALSE

# Check for and flag overlap
for (i in 1:nrow(subset(NC_PC_surveyinfo, to_evaluate))){
  
  obs_to_evaluate <- subset(NC_PC_surveyinfo, to_evaluate)[i,]
  dists <- distance_matrix[i,]
  WT_within_5m <- WT_surveyinfo[which(dists <= units::set_units(5,"m")),]
  time_diff <- abs(difftime(obs_to_evaluate$Date_Time, WT_within_5m$Date_Time, units='mins')) %>% as.numeric()
  
  # If the survey was ARU-based, use WildTrax as authoritative
  if (min(time_diff)<=10 & obs_to_evaluate$Survey_Type != "IN_PERSON") NC_PC_surveyinfo$to_remove[which(NC_PC_surveyinfo$Obs_Index == obs_to_evaluate$Obs_Index)] <- TRUE
  
  # If the survey was Human-based, use NatureCounts as authoritative
  if (min(time_diff)<=10 & obs_to_evaluate$Survey_Type == "IN_PERSON") WT_surveyinfo$to_remove[which(WT_surveyinfo$survey_ID %in% WT_within_5m[which(time_diff<=10),]$survey_ID)] <- TRUE
}

# Remove duplicated records from NatureCounts dataset
PC_removed <- subset(NC_PC_surveyinfo,to_remove)
NC_PC_surveyinfo <- subset(NC_PC_surveyinfo,!to_remove)
NC_PC_matrix <- PC_matrix[NC_PC_surveyinfo$Obs_Index,]
NC_PC_surveyinfo <- NC_PC_surveyinfo %>% dplyr::select(-Obs_Index,-to_evaluate, -to_remove)

# Remove duplicated records from WildTrax dataset
WT_removed <- subset(WT_surveyinfo,to_remove)
WT_surveyinfo <- subset(WT_surveyinfo,!to_remove)
WT_matrix <- WT_matrix[WT_surveyinfo$Obs_Index,]
WT_surveyinfo <- WT_surveyinfo %>% dplyr::select(-Obs_Index, -to_remove)

# ************************************************************
# ************************************************************
# COMBINE ALL DATASETS INTO A SINGLE FILE
# ************************************************************
# ************************************************************

WT_surveyinfo$survey_ID <- as.character(WT_surveyinfo$survey_ID)
NC_PC_surveyinfo$survey_ID <- as.character(NC_PC_surveyinfo$survey_ID)
DO_surveyinfo$survey_ID <- as.character(DO_surveyinfo$survey_ID)

all_surveys <- bind_rows(WT_surveyinfo %>% st_transform(crs = AEA_proj) %>% mutate(Survey_Class = "WT"),
                         NC_PC_surveyinfo %>% st_transform(crs = AEA_proj) %>% mutate(Survey_Class = "NC_PC"),
                         DO_surveyinfo %>% st_transform(crs = AEA_proj) %>% mutate(Survey_Class = "DO"))

# Set time zone
tz(all_surveys$Date_Time) <- "Canada/Central"

# Hours since sunrise
all_surveys$Sunrise <- suntools::sunriset(crds = st_transform(all_surveys,crs = 4326),
                                          dateTime = all_surveys$Date_Time,
                                          direction = c("sunrise"),
                                          POSIXct.out = TRUE)$time
all_surveys$Hours_Since_Sunrise <- with(all_surveys, difftime(Date_Time,Sunrise,units="hours") )


# -----------------------------------------------
# Combine species count matrices
# -----------------------------------------------

NC_species <- colnames(PC_matrix)
WT_species <- colnames(WT_matrix)
DO_species <- colnames(DO_matrix)

# NC_species only detected in NatureCounts point counts
NC_only <- NC_species[NC_species %!in% WT_species & NC_species %!in% DO_species]

# NC_species only detected in WildTrax
WT_only <- WT_species[WT_species %!in% NC_species & WT_species %!in% DO_species]

# NC_species only detected in checklists
DO_only <- DO_species[DO_species %!in% NC_species & DO_species %!in% WT_species]

subset(all_species,BSC_spcd %in% NC_only)
subset(all_species,BSC_spcd %in% WT_only)
subset(all_species,BSC_spcd %in% DO_only)

species_combined <- c(NC_species,WT_species,DO_species) %>% unique() %>% sort()

full_count_matrix <- matrix(0, nrow=nrow(all_surveys), ncol = length(species_combined),
                            dimnames = list(all_surveys$survey_ID,species_combined))

for (spp in species_combined){
  
  if (spp %in% colnames(WT_matrix)) full_count_matrix[which(all_surveys$Survey_Class == "WT"),spp] <- WT_matrix[,spp] 
  if (spp %in% colnames(PC_matrix)) full_count_matrix[which(all_surveys$Survey_Class == "NC_PC"),spp] <- NC_PC_matrix[,spp] 
  if (spp %in% colnames(DO_matrix)) full_count_matrix[which(all_surveys$Survey_Class == "DO"),spp] <- DO_matrix[,spp] 
  
}

# ************************************************************
# ************************************************************
# CHECK FOR OVERLAP BETWEEN STATIONARY COUNT CHECKLISTS AND POINT COUNTS
# - REMOVE DUPLICATES
# ************************************************************
# ************************************************************

# all_surveys$Obs_Index <- 1:nrow(all_surveys)
# 
# point_counts <- subset(all_surveys, Survey_Type %in% c("Point_Count","IN_PERSON",
#                                                        
#                                                        "ARU_SPT","ARU_SPM",
#                                                        "ARU_SM2","ARU_BAR_LT","ARU_SM_UNKN",
#                                                        "ARU_SM4","ZOOM_H2N","ARU_IRIVER_E",
#                                                        "ARU_MARANTZ","ARU_IRIVER","ARU_UNKNOWN"))
# 
# stationary_counts <- subset(all_surveys, Survey_Type %in% c("Breeding Bird Atlas")) 
# -----------------------------------------------
# Save 
# -----------------------------------------------

analysis_data <- list(all_surveys = all_surveys,
                      full_count_matrix = full_count_matrix)

saveRDS(analysis_data, file = "../Data_Cleaned/analysis_data.rds")

# ************************************************************
# ************************************************************
# Plot Data Availability
# ************************************************************
# ************************************************************

SK_BCR <- st_read("../../../Data/Spatial/National/BCR/BCR_Terrestrial_master.shp")  %>%
  subset(PROVINCE_S == "SASKATCHEWAN") %>%
  st_make_valid() %>%
  st_transform(st_crs(AEA_proj)) %>%
  dplyr::select(BCR, PROVINCE_S)

SaskBoundary <- SK_BCR %>% st_union()

all_surveys$Survey_Class <- factor(all_surveys$Survey_Class, levels = c("WT","NC_PC","DO"))
ggplot()+
  geom_sf(data = SK_BCR, fill = "gray95", col = "transparent") +
  geom_sf(data = subset(all_surveys, year(Date_Time) > 2010),aes(col = Survey_Class), size = 0.2)+
  geom_sf(data = SaskBoundary, fill = "transparent", col = "black") +
  scale_color_manual(values=c("orangered","black","dodgerblue"),name = "Data Source")+
  facet_grid(.~Survey_Class)+
  ggtitle("Data availability")
