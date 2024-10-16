# ************************************************************************
# Prepare data for analysis
#  - select date / time windows for survey inclusion
#  - merge in covariates
# ************************************************************************

# ------------------------------------------------
# Load/install packages and set graphical themes / working directory
# ------------------------------------------------

my_packs = c('tidyverse',
             'sf',
             'lubridate',
             'suntools',
             'factoextra',
             'viridis',
             'napops',
             'naturecounts',
             'terra',
             'exactextractr',
             'magrittr')

if (any(!my_packs %in% installed.packages()[, 'Package'])) {install.packages(my_packs[which(!my_packs %in% installed.packages()[, 'Package'])],dependencies = TRUE)}
lapply(my_packs, require, character.only = TRUE)


rm(list=ls())

# ------------------------------------------
# Set working directory
# ------------------------------------------

setwd("C:/Users/IlesD/OneDrive - EC-EC/Iles/Projects/Landbirds/Landbird-Distribution-Modeling-ECCC/Analysis/Ontario_BBA/Code")

`%!in%` <- Negate(`%in%`)

# ******************************************************************
# PART 1: Select data meeting criteria for inclusion (dates, time since sunrise, etc)
# ******************************************************************

AEA_proj <- "+proj=aea +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-87 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs "

analysis_data <- readRDS(file = "../Data_Cleaned/analysis_data.rds")

all_surveys <- analysis_data$all_surveys %>% mutate(Obs_Index = 1:nrow(.))
full_count_matrix <- analysis_data$full_count_matrix
all_surveys$Survey_Type[all_surveys$Survey_Type == "Point Count"] <- "Point_Count"

# ------------------------------------------
# Round survey duration
# ------------------------------------------
all_surveys$Survey_Duration_Minutes <- round(all_surveys$Survey_Duration_Minutes)
table(all_surveys$Survey_Duration_Minutes)

# ------------------------------------------
# Select Point Counts / ARUs to use
# ------------------------------------------

PC_to_use <- subset(all_surveys,
                    Survey_Type %in% c("Point_Count","ARU_SPT") &
                      
                      Survey_Duration_Minutes >= 1 &
                      Survey_Duration_Minutes <= 10 &
                      
                      Hours_Since_Sunrise >= -2 &
                      Hours_Since_Sunrise <= 4 &
                      
                      yday(Date_Time) >= yday(ymd("2022-05-15")) &
                      yday(Date_Time) <= yday(ymd("2022-07-15")) #&
                      
                      #year(Date_Time) >= 2021 &
                      #year(Date_Time) <= 2025
)
# 
# # ------------------------------------------
# # Select STATIONARY COUNT data to use
# # ------------------------------------------
# 
# # Select stationary counts to use
# SC_to_use <- subset(all_surveys,
#                     
#                     Survey_Type %in% c("Stationary Count") & 
#                       
#                       Hours_Since_Sunrise >= -2 &
#                       Hours_Since_Sunrise <= 4 &
#                       
#                       Survey_Duration_Minutes >= 1 &
#                       Survey_Duration_Minutes <= 120 &
#                       
#                       yday(Date_Time) >= yday(ymd("2022-05-15")) &
#                       yday(Date_Time) <= yday(ymd("2022-07-15")) & 
#                       
#                       year(Date_Time) >= 2021 &
#                       year(Date_Time) <= 2025)
# 
# # ------------------------------------------
# # Select LINEAR TRANSECT data to use
# # ------------------------------------------
# 
# LT_to_use <- subset(all_surveys,
#                     Survey_Type %in% c("Linear transect") & 
#                       
#                       Hours_Since_Sunrise >= -2 &
#                       Hours_Since_Sunrise <= 15 &
#                       
#                       Survey_Duration_Minutes > 10 &
#                       Survey_Duration_Minutes <= 120 &
#                       Travel_Distance_Metres <= 10000 &
#                       
#                       yday(Date_Time) >= yday(ymd("2022-05-28")) &
#                       yday(Date_Time) <= yday(ymd("2022-07-07")) &
#                       
#                       year(Date_Time) >= 2017 &
#                       year(Date_Time) <= 2021)

# ------------------------------------------
# Subset
# ------------------------------------------

surveys_to_use <- c(PC_to_use$Obs_Index) # , SC_to_use$Obs_Index , LT_to_use$Obs_Index
all_surveys <- subset(all_surveys, Obs_Index %in% surveys_to_use)
full_count_matrix <- full_count_matrix[surveys_to_use,]

# ------------------------------------------
# Atlas Squares in which each survey is located
# ------------------------------------------

# ONSquares <- st_read("../../../Data/Spatial/Ontario/ONSquares/ONSquares.shp") %>%
#   st_transform(st_crs(all_surveys)) %>%
#   dplyr::select(SQUARE_ID) %>%
#   rename(sq_id = SQUARE_ID)

all_surveys <- all_surveys %>% 
  mutate(Obs_Index = 1:nrow(.)) #%>%st_intersection(ONSquares)
full_count_matrix <- full_count_matrix[all_surveys$Obs_Index,]

all_surveys$Obs_Index <- 1:nrow(all_surveys)

# ------------------------------------------
# CONVERT TO APPROPRIATE CRS (of eventual map projection)
# ------------------------------------------

# Raster with target properties
target_crs <- AEA_proj

all_surveys <- st_transform(all_surveys,target_crs)


# ******************************************************************
# ******************************************************************
# PART 2: SELECT COVARIATES FOR ANALYSIS
# ******************************************************************
# ******************************************************************

# Dataframe to store covariates
all_surveys_covariates <- all_surveys %>% dplyr::select(Obs_Index)

# ---------------------------------------------------
# Prepare a grid on which maps/predictions will eventually be made
# ---------------------------------------------------

ONBoundary <- st_read("../../../Data/Spatial/National/BCR/BCR_Terrestrial_master.shp")  %>%
  subset(PROVINCE_S == "ONTARIO") %>%
  st_make_valid() %>%
  st_union() %>%
  st_transform(target_crs)

# # 1 km x 1 km grid
# ONGrid <- st_make_grid(
#   ONBoundary,
#   cellsize = units::set_units(1*1,km^2),
#   what = "polygons",
#   square = TRUE,
#   flat_topped = FALSE)%>%
#   st_as_sf() %>%
#   st_intersection(ONBoundary) %>%
#   na.omit()
# 
# ONGrid$point_id <- 1:nrow(ONGrid)
# ONGrid_centroid <- st_centroid(ONGrid)

# ---------------------------------------------------
# Load covariate layers
# ---------------------------------------------------

# Path to spatial covariates
covar_folder <- "../../../Data/Spatial/"

# Add 20 km buffer so that covariates can extend slightly outside province if possible
ONBoundary_buffer <- ONBoundary %>% st_buffer(20000)

# Annual mean temperature
AMT <- rast(paste0(covar_folder,"National/AnnualMeanTemperature/wc2.1_30s_bio_1.tif")) %>% 
  crop(st_transform(ONBoundary_buffer,crs(.))) %>%  
  project(st_as_sf(ONBoundary_buffer), res = 250)

# Land cover of Canada 2020 (not reprojected... faster to reproject grid)
lcc2020 <- rast(paste0(covar_folder,"National/LandCoverCanada2020/landcover-2020-classification.tif")) %>% 
  crop(st_transform(ONBoundary_buffer,crs(.))) 

# Stand canopy closure
SCC <- rast(paste0(covar_folder,"National/NationalForestInventory/NFI_MODIS250m_2011_kNN_Structure_Stand_CrownClosure_v1.tif")) %>% 
  crop(st_transform(ONBoundary_buffer,crs(.))) %>% 
  project(st_as_sf(ONBoundary_buffer), res = 250)

# ---------------------------------------------------
# For each survey, extract covariates within 1 km of location
# ---------------------------------------------------

# Continuous covariates
all_surveys_1km <- all_surveys %>%
  st_buffer(1000) %>%
  mutate(#elevation_1km = exact_extract(elevation, ., 'mean'),
         AMT_1km = exact_extract(AMT, ., 'mean'),
         SCC_1km = exact_extract(SCC, ., 'mean'))

# Proportion of each land cover class
prop_LCC_1km <- exact_extract(lcc2020,all_surveys_1km,"frac") %>% suppressWarnings()
names(prop_LCC_1km) <- paste0(str_replace(names(prop_LCC_1km),"frac","LCC"),"_1km")
prop_LCC_1km[setdiff(paste0("LCC_",seq(1,18),"_1km"),names(prop_LCC_1km))] <- 0
prop_LCC_1km %<>% dplyr::select(sort(names(.)))

# Fix land cover class names
prop_LCC_1km <- prop_LCC_1km %>%
  mutate(
    Needleleaf_forest_1km = LCC_1_1km + LCC_2_1km,
    Mixed_forest_1km = LCC_5_1km + LCC_6_1km,
    Grass_shrub_1km = LCC_8_1km + LCC_10_1km,
    Crop_1km = LCC_15_1km,
    Urban_1km = LCC_17_1km,
    Wetland_1km = LCC_14_1km,
    Water_1km = LCC_18_1km) %>%
  dplyr::select(Needleleaf_forest_1km:Water_1km)

prop_LCC_1km$Obs_Index <- all_surveys_1km$Obs_Index
all_surveys_1km <- left_join(all_surveys_1km,prop_LCC_1km)

# Select covariate columns
all_surveys_1km <- all_surveys_1km %>%
  as.data.frame() %>%
  dplyr::select(Obs_Index,colnames(all_surveys_1km)[colnames(all_surveys_1km) %!in% colnames(all_surveys)])

# Join with survey dataset
all_surveys_covariates <- all_surveys_covariates %>% left_join(all_surveys_1km)

# ---------------------------------------------------
# For each point in ONGrid, extract covariates within 1 km of location
# ---------------------------------------------------

# Continuous covariates
ONGrid_1km <- ONGrid_centroid %>%
  st_buffer(1000) %>%
  mutate(#elevation_1km = exact_extract(elevation, ., 'mean'),
         AMT_1km = exact_extract(AMT, ., 'mean'),
         SCC_1km = exact_extract(SCC, ., 'mean'))

# Proportion of each land cover class
prop_LCC_1km <- exact_extract(lcc2020,ONGrid_1km,"frac") %>% suppressWarnings()
names(prop_LCC_1km) <- paste0(str_replace(names(prop_LCC_1km),"frac","LCC"),"_1km")
prop_LCC_1km[setdiff(paste0("LCC_",seq(1,18),"_1km"),names(prop_LCC_1km))] <- 0
prop_LCC_1km %<>% dplyr::select(sort(names(.)))

# Fix land cover class names
prop_LCC_1km <- prop_LCC_1km %>%
  mutate(
    Needleleaf_forest_1km = LCC_1_1km + LCC_2_1km,
    Mixed_forest_1km = LCC_5_1km + LCC_6_1km,
    Grass_shrub_1km = LCC_8_1km + LCC_10_1km,
    Crop_1km = LCC_15_1km,
    Urban_1km = LCC_17_1km,
    Wetland_1km = LCC_14_1km,
    Water_1km = LCC_18_1km) %>%
  dplyr::select(Needleleaf_forest_1km:Water_1km)

prop_LCC_1km$point_id <- ONGrid_1km$point_id

ONGrid_1km <- left_join(ONGrid_1km,prop_LCC_1km)

# Join with ONGrid
ONGrid <- ONGrid %>% left_join(as.data.frame(ONGrid_1km) %>% dplyr::select(-x))

# --------------------------------
# Conduct Principal Components Analysis on covariates to identify axes of major variation in habitat
# --------------------------------

# Remove surveys with no covariate information
surveys_to_remove <- as.data.frame(all_surveys_covariates) %>% dplyr::select(-geometry)
surveys_to_remove <- which(is.na(rowSums(surveys_to_remove)))
all_surveys_covariates <- all_surveys_covariates[-surveys_to_remove,]
all_surveys <- all_surveys[-surveys_to_remove,]
full_count_matrix <- full_count_matrix[-surveys_to_remove,]

covars_for_PCA <- all_surveys_covariates %>%
  as.data.frame() %>%
  dplyr::select(AMT_1km:Water_1km)

pca <- prcomp(covars_for_PCA, scale = TRUE)

# ------------------------------------------
# Interpretation of specific axes (e.g., axes 1 and 2)
# ------------------------------------------

summary(pca)   # Proportion variance explaind by axes
fviz_eig(pca)  # Scree plot (first 5 axes explain 85% of variation in habitat between sites)
pca            # Variable loadings

fviz_pca_var(pca,
             axes = c(1,2),
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = viridis(10),
             repel = TRUE     # Avoid text overlapping
)


# ------------------------------------------
# Predict PCA values for each survey location and standardize (mean = 0, sd = 1)
# ------------------------------------------

all_surveys_PCA <- predict(pca, newdata = as.data.frame(all_surveys_covariates)[,names(pca$center)])
ONGrid_PCA <- predict(pca, newdata = as.data.frame(ONGrid)[,names(pca$center)])

for (covar in colnames(all_surveys_PCA)){

  covar_mean <- mean(all_surveys_PCA[,covar],na.rm = TRUE)
  covar_sd <- sd(all_surveys_PCA[,covar],na.rm = TRUE)

  all_surveys_PCA[,covar] <- (as.data.frame(all_surveys_PCA)[,covar] - covar_mean)/covar_sd
  ONGrid_PCA[,covar] <- (as.data.frame(ONGrid_PCA)[,covar] - covar_mean)/covar_sd
  
}

all_surveys_covariates <- all_surveys_covariates %>%
  as.data.frame() %>%
  dplyr::select(-geometry) %>%
  bind_cols(all_surveys_PCA)


all_surveys <- full_join(all_surveys,all_surveys_covariates)
ONGrid <- bind_cols(ONGrid,ONGrid_PCA)

# ONGrid_centroid <- st_centroid(ONGrid)
# ONGrid_list = list(ONGrid = ONGrid,
#                    ONGrid_centroid = ONGrid_centroid)
# saveRDS(ONGrid_list,
#         file = "../Data_Cleaned/Spatial/ONGrid.rds")

# ******************************************************************
# PART 3: Identify list of species to run analysis for
# ******************************************************************

# # ------------------------------------------
# # Process species names / labels (from Birds Canada)
# # ------------------------------------------
# 
# ON_spcd <- search_species_code() %>% rename(spcd = BSCDATA, CommonName = english_name)
# countSpaces <- function(s) { sapply(gregexpr(" ", s), function(p) { sum(p>=0) } ) }
# 
# # Process Species Names so they fit
# ON_spcd$Label <- NA
# ON_spcd$CommonName[ON_spcd$CommonName=="Rock Pigeon (Feral Pigeon)"] <- "Rock Pigeon"
# 
# for (i in 1:nrow(ON_spcd)) {
#   Name <- ON_spcd$CommonName[i]
#   if(nchar(Name) > 13){
#     if(countSpaces(Name)>0){
#       ON_spcd$Label[i] <- gsub(" "," \n",Name)
#     }
#     
#   }
#   else {
#     ON_spcd$Label[i] <- ON_spcd$CommonName[i]
#   }
#   
# }

# ------------------------------------------
# Only fit models for species detected in at least 50 atlas squares)
# ------------------------------------------


BSC_species <- readRDS("../Data_Cleaned/BSC_species.RDS")

n_detections <- full_count_matrix %>% 
  reshape2::melt() %>%
  dplyr::rename(sq_id = Var1, Species_Code_BSC = Var2, detected = value) %>%
  subset(detected>0) %>%
  group_by(Species_Code_BSC) %>%
  summarize(n_detections = sum(detected>0)) %>%
  arrange(desc(n_detections)) %>%
  left_join(.,BSC_species[,c("species_id","english_name")], by = c("Species_Code_BSC" = "species_id"))

species_to_model <- n_detections %>%
  
  # Remove erroneous observations
  subset(english_name %!in% c("passerine sp.",
                              "duck sp.",
                              "new world sparrow sp.",
                              "new world warbler sp.",
                              "blackbird sp.",
                              "gull sp.",
                              "Greater/Lesser Yellowlegs",
                              "Scolopacidae sp.",
                              "vireo sp.",
                              "Catharus sp.",
                              "Worthen's Sparrow",
                              "woodpecker sp.",
                              "Alder/Willow Flycatcher (Traill's Flycatcher)",
                              "Aythya sp.",
                              "Bat sp.",
                              "bird sp.",
                              "chickadee sp.",
                              "finch sp.",
                              "Golden-winged/Blue-winged Warbler",
                              "goose sp.",
                              "hawk sp.",
                              "new world flycatcher sp.",
                              "Northern Flicker (Yellow-shafted)",
                              "owl sp.",
                              "Philadelphia/Red-eyed Vireo",
                              "swallow sp.",
                              "swan sp.",
                              "tern sp.",
                              "Ula-ai-hawane",
                              "Yellow-rumped Warbler (Myrtle)",
                              "Brewster's Warbler (hybrid)",
                              "moorhen/coot/gallinule sp.",
                              "Lesser/Greater Yellowlegs",
                              "Coccyzus sp.",
                              "Bohemian/Cedar Waxwing",
                              "Chuck-will's-widow"))

dim(species_to_model) 
species_to_model$english_name %>% sort()

# ******************************************************************
# PART 4: IDENTIFY SPECIES WITH DETECTABILITY OFFSETS AVAILABLE
# ******************************************************************
# 
# napops_species <- list_species() %>% rename(Species_Code_NAPOPS = Species,
#                                             Common_Name_NAPOPS = Common_Name,
#                                             Scientific_Name_NAPOPS = Scientific_Name)
# 
# species_to_model <- left_join(species_to_model,napops_species[,c("Species_Code_NAPOPS","Common_Name_NAPOPS","Removal","Distance")],
#                               by = c("Species_Code_BSC" = "Species_Code_NAPOPS"))
# 
# 
# species_to_model$offset_exists <- FALSE
# species_to_model$EDR <- NA
# species_to_model$cue_rate <- NA
# species_to_model$log_offset_5min <- 0
# 
# # Extract QPAD offsets if available
# for (i in 1:nrow(species_to_model)){
#   
#   sp = species_to_model$Species_Code_BSC[i]
#   print(sp)
#   offset_exists <- FALSE
#   sp_napops <- subset(napops_species,Species_Code_NAPOPS == sp)
#   
#   if (nrow(sp_napops)>0){
#     if (sp_napops$Removal == 1 & sp_napops$Distance == 1){
#       
#       species_to_model$offset_exists[i] <- TRUE
#       species_to_model$cue_rate[i] <- cue_rate(species = sp,od = 153, tssr = 0, model = 1)[3] %>% as.numeric()
#       species_to_model$EDR[i] <- edr(species = sp,road = FALSE, forest = 0.5,model = 1)[3] %>% as.numeric()
#       
#       # Calculate A and p, which jointly determine offset
#       A_metres <- c(pi*species_to_model$EDR[i]^2)
#       p <- 1-exp(-5*species_to_model$cue_rate[i])
#       
#       species_to_model$log_offset_5min[i] <- log(A_metres * p)
#       
#     }
#   }
# }

# ******************************************************************
# PART 5: DOWNLOAD, TRIM, AND SAVE eBIRD RANGE LIMITS FOR SPECIES
# ******************************************************************

# ------------------------------------------------
# Load packages
# ------------------------------------------------

require(ebirdst)
require(terra)

# Need to set access key (only once) before downloading ranges
# ebirdst::set_ebirdst_access_key()

# ------------------------------------------------
# Ensure path is correctly set
# ------------------------------------------------
usethis::edit_r_environ()

# # Spatial layers for ON
# ON_BCR <- st_read("../../../Data/Spatial/National/BCR/BCR_Terrestrial_master.shp")  %>%
#   subset(PROVINCE_S == "ONTARIO") %>%
#   st_make_valid() %>%
#   dplyr::select(BCR, PROVINCE_S)
# 
# ONBoundary <- ON_BCR %>% st_union()

# List to contain sf objects for species ranges
species_ranges <- list()

for (species_name in species_to_model$english_name){
  print(species_name)
  if (file.exists("../Data_Cleaned/Spatial/eBird_ranges_ON.RDS")) species_ranges <- readRDS("../Data_Cleaned/Spatial/eBird_ranges_ON.RDS")
  
  # ------------------------------------------------
  # Check is species already in the list
  # ------------------------------------------------
  if (species_name %in% names(species_ranges)) next
  
  # ------------------------------------------------
  # Download and trim ebird range for this species
  # ------------------------------------------------
  
  # Check if species is available
  check <- get_species(species_name)
  if (length(check)==0){
    print(paste0(species_name, " not available in eBird"))
    next
  }
  
  if (length(check)>0 & is.na(check)){
    print(paste0(species_name, " not available in eBird"))
    next
  }
  
  ebirdst_download_status(species_name, 
                          download_abundance = FALSE,
                          download_ranges = TRUE,
                          pattern = "_smooth_27km_") 
  
  path <- get_species_path(species_name)
  
  range <- load_ranges(species_name, resolution = "27km",smoothed = TRUE)
  
  range <- range %>% subset(season %in% c("resident","breeding")) %>% 
    st_transform(.,crs = st_crs(ONBoundary)) %>% 
    st_union() %>%
    st_crop(ONBoundary)
  
  species_ranges[[species_name]] <- range
  saveRDS(species_ranges,file="../Data_Cleaned/Spatial/eBird_ranges_ON.RDS")
  
  
} # close species loop


# ******************************************************************
# ******************************************************************
# Save
# ******************************************************************
# ******************************************************************

analysis_data_package <- list(
  
  all_surveys = all_surveys, # Survey information (protocol, location, time, etc)
  full_count_matrix = full_count_matrix, # counts of each species for each survey
  
  pca = pca,
  
  # Contains covariates on ON-wide grid
  ONGrid = ONGrid,
  
  # Species to include in analysis (and number of detections in point counts)
  species_to_model = species_to_model,
  
  # Species codes and common names
  species_to_model = species_to_model,
  
  # eBird range limits
  species_ranges = species_ranges
  
)

saveRDS(analysis_data_package,"../Data_Cleaned/analysis_data_package.rds")
