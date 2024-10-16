# ---------------------------------------------------------------------------
# ***************************************************************************
# Download NatureCounts datasets within a study area boundary
#  - merge with WildTrax data, and remove duplicated records
# ***************************************************************************
# ---------------------------------------------------------------------------

# ------------------------------------------------
# Load/install packages and set graphical themes / working directory
# ------------------------------------------------

my_packs = c('tidyverse','magrittr','sf','terra','ggspatial','naturecounts','suntools','readxl')
if (any(!my_packs %in% installed.packages()[, 'Package'])) {install.packages(my_packs[which(!my_packs %in% installed.packages()[, 'Package'])],dependencies = TRUE)}
lapply(my_packs, require, character.only = TRUE)

rm(list=ls())

# Set working directory
setwd("C:/Users/IlesD/OneDrive - EC-EC/Iles/Projects/Landbirds/Landbird-Distribution-Modeling-ECCC/Analysis/Ontario_BBA/Code")

# Useful functions
`%!in%` <- Negate(`%in%`)

all_species <- readRDS("../Data_Cleaned/all_species.RDS")

# ************************************************************
# ************************************************************
# PROCESS OBBA-3 DATA
#  -> Last downloaded/run on 12 Sep 2024
# ************************************************************
# ************************************************************

# Point_Counts
PCData_OBBA3 <- read.table("../../../Data/Bird_Data_Raw/NatureCounts/OBBA3/onatlas3pc_naturecounts_data.txt", sep="\t",
                           header=TRUE, fill = TRUE,comment.char = "",quote = "") %>%
  dplyr::select(where(~sum(!is.na(.x)) > 0)) %>%                            # Remove empty columns
  rename(species_id = SpeciesCode) %>%                                      # Rename species code

  # Add relevant columns
  mutate(Date_Time = ymd(ObservationDate) + hours(floor(TimeObservationsStarted)) + minutes(round(60*(TimeObservationsStarted-floor(TimeObservationsStarted)))),
         Data_Source = "NatureCounts",
         Project_Name = "OBBA3",
         Survey_Type = "Point_Count") %>%

  # Only use relevant columns
  dplyr::select(Data_Source,
                Project_Name,
                Survey_Type,
                DecimalLatitude,
                DecimalLongitude,
                Date_Time,
                DurationInHours,
                species_id,
                ObservationCount)

saveRDS(PCData_OBBA3,"../Data_Cleaned/NatureCounts/PCData_OBBA3.rds")

# # Checklists
# DOData_OBBA3  <- read.table("../../../Data/Bird_Data_Raw/NatureCounts/OBBA3/onatlas3be_do_naturecounts_data.txt", sep="\t",
#                             header=TRUE, fill = TRUE,comment.char = "",quote = "") %>%
#   dplyr::select(where(~sum(!is.na(.x)) > 0)) %>%                            # Remove empty columns
#   rename(species_id = SpeciesCode) %>%                                      # Rename species code
#   left_join(., species %>% dplyr::select(species_id,species_id)) %>%              # Join species codes/names from naturecounts
#   filter(BreedingBirdAtlasCode!="X",!is.na(BreedingBirdAtlasCode))          # Remove records with 'X' or 'NA' for breeding code
# saveRDS(DOData_OBBA3,"../Data_Cleaned/NatureCounts/DOData_OBBA3.rds")

# ************************************************************
# ************************************************************
# PROCESS OBBA-2 DATA
# - Note that for OBBA-2 data, it appears that SpeciesCode is the 4 letter abbreviation (named 'species_id') rather than a numeric code (unlike atlas 3 dataset)
# ************************************************************
# ************************************************************

# Point_Counts
PCData_OBBA2 <- read.table("../../../Data/Bird_Data_Raw/NatureCounts/OBBA2/Point_Count/obba2pc_naturecounts_data.txt", sep="\t",
                           header=TRUE, fill = TRUE,comment.char = "",quote = "") %>%
  dplyr::select(where(~sum(!is.na(.x)) > 0)) %>% 
  left_join(.,all_species[,c("BSC_spcd","species_id")], by = c("SpeciesCode" = "BSC_spcd"), multiple = "first") %>%
  # Remove empty columns
  rename(TimeObservationsStarted = TimeCollected) %>%   # Rename species code

  # Add relevant columns
  mutate(Date_Time = ymd(paste(YearCollected,MonthCollected,DayCollected,sep="-")) + hours(floor(TimeObservationsStarted)) + minutes(round(60*(TimeObservationsStarted-floor(TimeObservationsStarted)))),
         Data_Source = "NatureCounts",
         Project_Name = "OBBA2",
         Survey_Type = "Point_Count") %>%

  # Only use relevant columns
  dplyr::select(Data_Source,
                Project_Name,
                Survey_Type,
                DecimalLatitude,
                DecimalLongitude,
                Date_Time,
                DurationInHours,
                species_id,
                ObservationCount)

saveRDS(PCData_OBBA2,"../Data_Cleaned/NatureCounts/PCData_OBBA2.rds")

# ### Breeding Evidence
# DOData_OBBA2 <- read.table("../../../Data/Bird_Data_Raw/NatureCounts/OBBA2/Raw_Breeding_Evidence/obba2be_raw_naturecounts_data.txt", sep="\t",
#                            header=TRUE, fill = TRUE,comment.char = "",quote = "") %>%
#   dplyr::select(where(~sum(!is.na(.x)) > 0)) %>%                            # Remove empty columns
#   rename(species_id = SpeciesCode) %>%                                            # Rename species code
#   left_join(., species %>% dplyr::select(species_id,species_id)) %>%              # Join species codes/names from naturecounts
#   filter(BreedingBirdAtlasCode!="X",!is.na(BreedingBirdAtlasCode))          # Remove records with 'X' or 'NA' for breeding code
#saveRDS(DOData_OBBA2,"../Data_Cleaned/NatureCounts/DOData_OBBA2.rds")

# ### Highest Breeding evidence
# BEData_OBBA2 <- read.table("../../../Data/Bird_Data_Raw/NatureCounts/OBBA2/Highest_Breeding_Evidence/obba2be_summ_naturecounts_data.txt", sep="\t",
#                            header=TRUE, fill = TRUE,comment.char = "",quote = "") %>%
#   dplyr::select(where(~sum(!is.na(.x)) > 0)) %>%                            # Remove empty columns
#   rename(species_id = SpeciesCode) %>%                                            # Rename species code
#   left_join(., species %>% dplyr::select(species_id,species_id)) %>%              # Join species codes/names from naturecounts
#   filter(BreedingBirdAtlasCode!="X",!is.na(BreedingBirdAtlasCode))          # Remove records with 'X' or 'NA' for breeding code
# saveRDS(BEData_OBBA2,"../Data_Cleaned/NatureCounts/BEData_OBBA2.rds")
# 
# 
# # ************************************************************
# # ************************************************************
# # Combine data from atlases, then separate survey information (lat, lon, date, etc) from counts
# # ************************************************************
# # ************************************************************
# 
# PCData <- bind_rows(PCData_OBBA3,PCData_OBBA2) %>%
#   mutate(survey_ID = as.numeric(as.factor(paste(Project_Name,DecimalLatitude,DecimalLongitude,Date_Time)))) %>%
#   subset(!is.na(DecimalLatitude) & !is.na(DecimalLongitude) & DecimalLongitude < -60 & !is.na(species_id) & !is.na(ObservationCount))
# 
# # Construct a dataset that contains JUST survey information (not counts)
# PC_surveyinfo <- PCData %>%
#   dplyr::select(-species_id,-ObservationCount) %>%
#   distinct() %>%
#   st_as_sf(.,coords = c("DecimalLongitude","DecimalLatitude"), crs = st_crs(4326), remove = FALSE) %>%
#   subset(!is.na(Date_Time))
# 
# 
# # Counts for each species on each survey are stored in a matrix
# PC_matrix <- matrix(0,nrow = nrow(PC_surveyinfo),ncol = nrow(all_species),
#                     dimnames = list(PC_surveyinfo$survey_ID,all_species$species_id))
# 
# # Populate matrix (takes a long time, but only needs to be done once)
# pb = txtProgressBar(min = 0, max = nrow(PCData), initial = 0, style = 3)
# for (i in 1:nrow(PCData)){
#   PC_matrix[which(rownames(PC_matrix) == PCData$survey_ID[i]),which(colnames(PC_matrix) == PCData$species_id[i])] <- PCData$ObservationCount[i]
#   setTxtProgressBar(pb,i)
# }
# close(pb)
# 
# # Remove columns (species) that were never observed to save storage space
# PC_matrix <- PC_matrix[,-which(colSums(PC_matrix,na.rm = T)==0)]
# 
# PC_dat <- list(PC_surveyinfo = PC_surveyinfo, PC_matrix = PC_matrix)
# saveRDS(PC_dat, file = "../Data_Cleaned/NatureCounts/PC_dat.rds")

PC_dat <- readRDS(file = "../Data_Cleaned/NatureCounts/PC_dat.rds")
PC_surveyinfo <- PC_dat$PC_surveyinfo
PC_matrix <- PC_dat$PC_matrix

# # ************************************************************
# # ************************************************************
# # Merge data with WildTrax, and check for duplicated data between the two sources
#       - assume ARU data in WildTrax data is 'authoritative' because it can be reviewed/updated more often
#       - assume Human Point_Counts in NatureCounts are authoritative because they have been thoroughly QA/QC'd
# # ************************************************************
# # ************************************************************

WT_dat <- readRDS(file = "../Data_Cleaned/WildTrax/WT_dat.rds")
WT_surveyinfo <- WT_dat$WT_surveyinfo
WT_matrix <- WT_dat$WT_matrix

# Remove any observations without date/time information
to_remove <- which(is.na(WT_surveyinfo$Date_Time))
WT_surveyinfo <- WT_surveyinfo[-to_remove,]
WT_matrix <- WT_matrix[-to_remove,]

# ************************************************************
# ************************************************************
# CHECK FOR OVERLAP BETWEEN NATURECOUNTS AND WILDTRAX Point_Count DATASETS
#
#  - assume ARU data in WildTrax data is 'authoritative' because it can be reviewed/updated more often
#  - assume Human Point_Counts in NatureCounts are authoritative because they have been thoroughly QA/QC'd
#
# - currently, deleting paired human/ARU records (only using human records in those cases - THIS NEEDS TO BE FIXED)
# ************************************************************
# ************************************************************

AEA_proj <- "+proj=aea +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-85 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs "

PC_surveyinfo <- st_transform(PC_surveyinfo, crs = AEA_proj)
WT_surveyinfo <- st_transform(WT_surveyinfo, crs = AEA_proj)

# Identify NatureCounts Point_Counts that are within 10 m of WildTrax locations
WT_buff <- st_buffer(WT_surveyinfo,10) %>% st_union()
PC_surveyinfo$to_evaluate <- FALSE
PC_surveyinfo$to_evaluate[which(as.matrix(st_intersects(PC_surveyinfo,WT_buff)))] <- TRUE
distance_matrix <- st_distance(subset(PC_surveyinfo, to_evaluate),WT_surveyinfo)

# Flags for which observations to remove
PC_surveyinfo$Obs_Index <- 1:nrow(PC_surveyinfo)
WT_surveyinfo$Obs_Index <- 1:nrow(WT_surveyinfo)
PC_surveyinfo$to_remove <- FALSE
WT_surveyinfo$to_remove <- FALSE

# Check for and flag overlap
for (i in 1:nrow(subset(PC_surveyinfo, to_evaluate))){
  
  obs_to_evaluate <- subset(PC_surveyinfo, to_evaluate)[i,]
  dists <- distance_matrix[i,]
  WT_within_10m <- WT_surveyinfo[which(dists <= units::set_units(10,"m")),]
  time_diff <- abs(difftime(obs_to_evaluate$Date_Time, WT_within_10m$Date_Time, units='mins')) %>% as.numeric()
  
  # Ifsurvey was ARU-based, use WildTrax as authoritative
  if (min(time_diff)<=10 & obs_to_evaluate$Survey_Type != "Point_Count") PC_surveyinfo$to_remove[which(PC_surveyinfo$Obs_Index == obs_to_evaluate$Obs_Index)] <- TRUE
  
  # If survey was Human-based, use NatureCounts as authoritative
  if (min(time_diff)<=10 & obs_to_evaluate$Survey_Type == "Point_Count") WT_surveyinfo$to_remove[which(WT_surveyinfo$survey_ID %in% WT_within_10m[which(time_diff<=10),]$survey_ID)] <- TRUE
  
  print(i)
}

# Remove duplicated records from NatureCounts dataset
PC_removed <- subset(PC_surveyinfo,to_remove)
PC_surveyinfo <- subset(PC_surveyinfo,!to_remove)
PC_matrix <- PC_matrix[PC_surveyinfo$Obs_Index,]
PC_surveyinfo <- PC_surveyinfo %>% dplyr::select(-Obs_Index,-to_evaluate, -to_remove)

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
PC_surveyinfo$survey_ID <- as.character(PC_surveyinfo$survey_ID)

all_surveys <- bind_rows(WT_surveyinfo %>% st_transform(crs = AEA_proj) %>% mutate(Survey_Class = "WT"),
                         PC_surveyinfo %>% st_transform(crs = AEA_proj) %>% mutate(Survey_Class = "PC"))

# Set time zone
tz(all_surveys$Date_Time) <- "Canada/Eastern"

# Hours since sunrise
all_surveys$Sunrise <- suntools::sunriset(crds = st_transform(all_surveys,crs = 4326),
                                          dateTime = all_surveys$Date_Time,
                                          direction = c("sunrise"),
                                          POSIXct.out = TRUE)$time
all_surveys$Hours_Since_Sunrise <- with(all_surveys, difftime(Date_Time,Sunrise,units="hours") )

# Fix survey duration when recorded in either hours or minutes
all_surveys$Survey_Duration_Minutes[is.na(all_surveys$Survey_Duration_Minutes)] <- all_surveys$DurationInHours[is.na(all_surveys$Survey_Duration_Minutes)]*60


# -----------------------------------------------
# Combine species count matrices
# -----------------------------------------------

NC_species <- colnames(PC_matrix)
WT_species <- colnames(WT_matrix)

# Species only detected in NatureCounts
NC_only <- NC_species[NC_species %!in% WT_species] 
NC_only <- colSums(PC_matrix[,NC_only]) %>% sort(decreasing = TRUE)
for (i in 1:length(NC_only)) names(NC_only)[i] <- subset(all_species, species_id == names(NC_only[i]))$english_name

# Species only detected in WildTrax
WT_only <- WT_species[WT_species %!in% NC_species ]
WT_only <- colSums(WT_matrix[,WT_only]) %>% sort(decreasing = TRUE)
for (i in 1:length(WT_only)) names(WT_only)[i] <- subset(all_species, species_id == names(WT_only[i]))$english_name


species_combined <- c(NC_species,WT_species) %>% unique() %>% sort()

full_count_matrix <- matrix(0, nrow=nrow(all_surveys), ncol = length(species_combined),
                            dimnames = list(all_surveys$survey_ID,species_combined))

for (spp in species_combined){

  if (spp %in% colnames(WT_matrix)) full_count_matrix[which(all_surveys$Survey_Class == "WT"),spp] <- WT_matrix[,spp]
  if (spp %in% colnames(PC_matrix)) full_count_matrix[which(all_surveys$Survey_Class == "PC"),spp] <- PC_matrix[,spp]
  #if (spp %in% colnames(DO_matrix)) full_count_matrix[which(all_surveys$Survey_Class == "DO"),spp] <- DO_matrix[,spp]

}

# -----------------------------------------------
# Combine species count matrices
# -----------------------------------------------


# # -----------------------------------------------
# # Save
# # -----------------------------------------------
# 
# analysis_data <- list(all_surveys = all_surveys,
#                       full_count_matrix = full_count_matrix)
# 
# saveRDS(analysis_data, file = "../Data_Cleaned/analysis_data.rds")


# ************************************************************
# ************************************************************
# Plot Data Availability
# ************************************************************
# ************************************************************
analysis_data <- readRDS(file = "../Data_Cleaned/analysis_data.rds")
all_surveys <- analysis_data$all_surveys
full_count_matrix <- analysis_data$full_count_matrix

all_surveys$Survey_Duration_Minutes <- round(all_surveys$Survey_Duration_Minutes)

all_surveys$sdur <- cut(all_surveys$Survey_Duration_Minutes,breaks = c(0,0.5,1.5,2.5,3.5,4.5,5.5,Inf), ordered_result = TRUE,right=FALSE)
table(all_surveys$sdur)

ON_BCR <- st_read("../../../Data/Spatial/National/BCR/BCR_Terrestrial_master.shp")  %>%
  subset(PROVINCE_S == "ONTARIO") %>%
  st_make_valid() %>%
  st_transform(st_crs(AEA_proj)) %>%
  dplyr::select(BCR, PROVINCE_S)

ONBoundary <- ON_BCR %>% st_union()

all_surveys$Survey_Class <- factor(all_surveys$Survey_Class, levels = c("WT","PC")) # "DO"
all_surveys$Atlas <- "OBBA-2"
all_surveys$Atlas[year(all_surveys$Date_Time) > 2020] <- "OBBA-3"

ggplot()+
  geom_sf(data = ON_BCR, fill = "gray95", col = "transparent") +
  geom_sf(data = all_surveys, aes(col = sdur), size = 0.2)+
  geom_sf(data = ONBoundary, fill = "transparent", col = "gray80") +
  facet_grid(Survey_Type~Atlas)+
  ggtitle("Data availability")
