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
             'naturecounts')

if (any(!my_packs %in% installed.packages()[, 'Package'])) {install.packages(my_packs[which(!my_packs %in% installed.packages()[, 'Package'])],dependencies = TRUE)}
lapply(my_packs, require, character.only = TRUE)

rm(list=ls())

# ------------------------------------------
# Set working directory
# ------------------------------------------

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

`%!in%` <- Negate(`%in%`)

# ******************************************************************
# PART 1: Select data meeting criteria for inclusion (dates, time since sunrise, etc)
# ******************************************************************

analysis_data <- readRDS(file = "../Data_Cleaned/analysis_data.rds")

all_surveys <- analysis_data$all_surveys %>% mutate(Obs_Index = 1:nrow(.))
full_count_matrix <- analysis_data$full_count_matrix

# ------------------------------------------
# Select Point Counts / ARUs to use
# ------------------------------------------

PC_to_use <- subset(all_surveys,
                    Survey_Type %in% c("Point_Count","ARU_SPT","ARU_SPM") &
                      Survey_Duration_Minutes > 1 &
                    Survey_Duration_Minutes <= 10 &
                      Hours_Since_Sunrise >= -2 &
                      Hours_Since_Sunrise <= 4 &
                      yday(Date_Time) >= yday(ymd("2022-05-31")) &
                      yday(Date_Time) <= yday(ymd("2022-07-07")))

# ------------------------------------------
# Select STATIONARY COUNT data to use (assuming these are 'breeding bird atlas' checklists; vast majority do not have distance)
# ------------------------------------------

# Select stationary counts to use
SC_to_use <- subset(all_surveys,
                    Survey_Type %in% c("Breeding Bird Atlas") & 
                      Hours_Since_Sunrise >= -2 &
                      Hours_Since_Sunrise <= 6 &
                      yday(Date_Time) >= yday(ymd("2022-05-31")) &
                      yday(Date_Time) <= yday(ymd("2022-07-07")) &
                      Survey_Duration_Minutes <= 60)


# ------------------------------------------
# Select LINEAR TRANSECT data to use
# ------------------------------------------

LT_to_use <- subset(all_surveys,
                    Survey_Type %in% c("Linear transect") & 
                      Hours_Since_Sunrise >= -2 &
                      Hours_Since_Sunrise <= 15 &
                      yday(Date_Time) >= yday(ymd("2022-05-31")) &
                      yday(Date_Time) <= yday(ymd("2022-07-07")) &
                      
                      Survey_Duration_Minutes > 10 &
                      Survey_Duration_Minutes <= 120 &
                      Travel_Distance_Metres <= 10000)

# ------------------------------------------
# Subset
# ------------------------------------------

surveys_to_use <- c(PC_to_use$Obs_Index, SC_to_use$Obs_Index, LT_to_use$Obs_Index)
all_surveys <- subset(all_surveys, Obs_Index %in% surveys_to_use)
full_count_matrix <- full_count_matrix[surveys_to_use,]

# ------------------------------------------
# Atlas Squares in which each survey is located
# ------------------------------------------

SaskSquares <- st_read("../../../Data/Spatial/Saskatchewan/SaskSquares/SaskSquares.shp") %>%
  st_transform(st_crs(all_surveys)) %>%
  dplyr::select(SQUARE_ID) %>%
  rename(sq_id = SQUARE_ID)

all_surveys <- all_surveys %>% 
  mutate(Obs_Index = 1:nrow(.)) %>%
  st_intersection(SaskSquares)
full_count_matrix <- full_count_matrix[all_surveys$Obs_Index,]

# ******************************************************************
# PART 2: PREPARE COVARIATES FOR ANALYSIS
# ******************************************************************

# note that covariate layers were prepared by script "2_Raster_Grid_Prep.R", on a 1km x 1km grid
SaskGrid <- read_rds("../Data_Cleaned/Spatial/SaskGrid_covariates.rds")

# ---------------------------------------------------
# Prepare covariates - combine several land cover classes
# ---------------------------------------------------

SaskGrid <- SaskGrid %>%
  mutate(
    
    Needleleaf_forest_5km = LCC_1_5km + LCC_2_5km,
    Mixed_forest_5km = LCC_5_5km + LCC_6_5km,
    Grass_shrub_5km = LCC_8_5km + LCC_10_5km,
    Crop_5km = LCC_15_5km,
    Urban_5km = LCC_17_5km,
    Wetland_5km = LCC_14_5km,
    Water_5km = LCC_18_5km,
    
    Needleleaf_forest_1km = LCC_1_1km + LCC_2_1km,
    Mixed_forest_1km = LCC_5_1km + LCC_6_1km,
    Grass_shrub_1km = LCC_8_1km + LCC_10_1km,
    Crop_1km = LCC_15_1km,
    Urban_1km = LCC_17_1km,
    Wetland_1km = LCC_14_1km,
    Water_1km = LCC_18_1km)

# --------------------------------
# Conduct Principal Components Analysis
# --------------------------------

covars <- as.data.frame(SaskGrid) %>%
  
  dplyr::select(
    
    Needleleaf_forest_5km,
    Mixed_forest_5km,
    Grass_shrub_5km,
    Crop_5km,
    Urban_5km,
    Wetland_5km,
    Water_5km,
    
    Needleleaf_forest_1km,
    Mixed_forest_1km,
    Grass_shrub_1km,
    Crop_1km,
    Urban_1km,
    Wetland_1km,
    Water_1km,
    
    # Other covariates
    elevation_5km, # elevation
    AMT_5km, # Annual mean temperature
    SCC_5km, # Stand canopy closure,
    elevation_1km, # elevation
    AMT_1km, # Annual mean temperature
    SCC_1km # Stand canopy closure,
    
  ) %>%
  na.omit()

pca <- prcomp(covars, scale = TRUE)

# ------------------------------------------
# Interpretation of specific axes (e.g., axes 1 and 2)
# ------------------------------------------

summary(pca)   # Proportion variance explaind by axes
fviz_eig(pca)  # Scree plot (first 7 axes explain most variation in habitat between sites)
pca            # Variable loadings

fviz_pca_var(pca,
             axes = c(1,2),
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = viridis(10),
             repel = TRUE     # Avoid text overlapping
)

# ------------------------------------------
# Save PCA values for each covariate grid pixel
# ------------------------------------------

PCA_values <- predict(pca,newdata = as.data.frame(SaskGrid)[,names(pca$center)]) %>% as.data.frame() 

# Join with SaskGrid
PCA_values$point_id <- SaskGrid$point_id
SaskGrid <- left_join(SaskGrid,PCA_values)

range(SaskGrid$PC1) # [1] -15.013819   5.020434

# ------------------------------------------
# Extract PCA values for each survey location
# ------------------------------------------

all_surveys$Obs_Index <- 1:nrow(all_surveys)

all_surveys <- all_surveys %>%
  st_transform(st_crs(SaskGrid)) %>%
  st_intersection(SaskGrid) %>%
  arrange(Obs_Index) %>% 
  relocate(geometry, .after = last_col())
full_count_matrix <- full_count_matrix[all_surveys$Obs_Index,]

# ------------------------------------------
# Identify covariates to include in 'omnibus' model
# ------------------------------------------

# Any variables with correlation higher than threshold are removed (reduce co-linearity)
threshold <- 0.5  # Note that higher thresholds will include more covariates

covar_df <- all_surveys %>%
  as.data.frame() %>%
  dplyr::select(elevation_1km:PC4) %>%
  relocate(PC1,PC2,PC3,PC4)

# Remove covariates that have no variation
covar_sd <- apply(covar_df,2,sd)
covar_df <- covar_df[,-which(covar_sd == 0)]

for (covar in colnames(covar_df)){
  
  if (!(covar %in% colnames(covar_df))) next
  
  # Remove covariates that have no variation
  if (sd(covar_df[,covar])==0){
    covar_df <- covar_df[,-which(colnames(covar_df) == covar)]
    next
  }
  
  # Check correlation with all other covariates
  cor_mat <- abs(cor(covar_df))
  cor_mat <- cor_mat[-which(rownames(cor_mat)==covar),covar]
  
  # Remove any covariates with correlation greater than threshold
  if (sum(abs(cor_mat)>threshold)>0){
    vars_to_remove <- names(cor_mat[which(cor_mat > threshold)])
    cols_to_remove <- which(colnames(covar_df) %in% vars_to_remove)
    covar_df <- covar_df[,-cols_to_remove]
  }
  
}

# Covariates to include in models
covariates_to_include <- colnames(covar_df)

covariates_to_include

# ******************************************************************
# Standardize covariates in analysis
# ******************************************************************

covar_df <- as.data.frame(all_surveys) %>% dplyr::select(covariates_to_include)

for (covar in covariates_to_include){
  
  covar_mean <- mean(covar_df[,covar],na.rm = TRUE)
  covar_sd <- sd(covar_df[,covar],na.rm = TRUE)
  
  SaskGrid[,covar] <- (as.data.frame(SaskGrid)[,covar] - covar_mean)/covar_sd
  all_surveys[,covar] <- (as.data.frame(all_surveys)[,covar] - covar_mean)/covar_sd
  
}

# ******************************************************************
# PART 3: Identify list of species to run analysis for
# ******************************************************************

# ------------------------------------------
# Process species names / labels (from Birds Canada)
# ------------------------------------------

Sask_spcd <- search_species_code() %>% rename(spcd = BSCDATA, CommonName = english_name)
countSpaces <- function(s) { sapply(gregexpr(" ", s), function(p) { sum(p>=0) } ) }

# Process Species Names so they fit
Sask_spcd$Label <- NA
Sask_spcd$CommonName[Sask_spcd$CommonName=="Rock Pigeon (Feral Pigeon)"] <- "Rock Pigeon"

for (i in 1:nrow(Sask_spcd)) {
  Name <- Sask_spcd$CommonName[i]
  if(nchar(Name) > 13){
    if(countSpaces(Name)>0){
      Sask_spcd$Label[i] <- gsub(" "," \n",Name)
    }
    
  }
  else {
    Sask_spcd$Label[i] <- Sask_spcd$CommonName[i]
  }
  
}

# ------------------------------------------
# Only fit models for species detected at least 50 times in point count & ARU surveys (this is an arbitrary lower bound)
# ------------------------------------------

n_detections <- (full_count_matrix[which(all_surveys$Survey_Type %in% c("Point_Count","ARU_SPT","ARU_SPM")),]>0) %>% 
  colSums() %>%
  sort(decreasing = TRUE)

species_to_model <- n_detections[n_detections >= 50] 

species_to_model <- data.frame(Species_Code_BSC = names(species_to_model),
                               Number_of_Detections = species_to_model)

BSC_species <- readRDS("../Data_Cleaned/BSC_species.RDS")

species_to_model <- left_join(species_to_model,BSC_species %>% dplyr::select(-index),
                              by = c("Species_Code_BSC" = "BSC_spcd")) %>%
  
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
                              "Worthen's Sparrow"))

dim(species_to_model) # 180 species

# ******************************************************************
# PART 4: IDENTIFY SPECIES WITH DETECTABILITY OFFSETS AVAILABLE
# ******************************************************************

napops_species <- list_species() %>% rename(Species_Code_NAPOPS = Species,
                                            Common_Name_NAPOPS = Common_Name,
                                            Scientific_Name_NAPOPS = Scientific_Name)

species_to_model <- left_join(species_to_model,napops_species[,c("Species_Code_NAPOPS","Common_Name_NAPOPS","Removal","Distance")],
                              by = c("Species_Code_BSC" = "Species_Code_NAPOPS"))


species_to_model$offset_exists <- FALSE
species_to_model$EDR <- NA
species_to_model$cue_rate <- NA
species_to_model$log_offset_5min <- 0

# Extract QPAD offsets if available
for (i in 1:nrow(species_to_model)){
  
  sp = species_to_model$Species_Code_BSC[i]
  print(sp)
  offset_exists <- FALSE
  sp_napops <- subset(napops_species,Species_Code_NAPOPS == sp)
  
  if (nrow(sp_napops)>0){
    if (sp_napops$Removal == 1 & sp_napops$Distance == 1){
      
      species_to_model$offset_exists[i] <- TRUE
      species_to_model$cue_rate[i] <- cue_rate(species = sp,od = 153, tssr = 0, model = 1)[3] %>% as.numeric()
      species_to_model$EDR[i] <- edr(species = sp,road = FALSE, forest = 0.5,model = 1)[3] %>% as.numeric()
      
      # Calculate A and p, which jointly determine offset
      A_metres <- c(pi*species_to_model$EDR[i]^2)
      p <- 1-exp(-5*species_to_model$cue_rate[i])
      
      species_to_model$log_offset_5min[i] <- log(A_metres * p)
      
    }
  }
}

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

# This should read:
# EBIRDST_KEY='ntqm1ha68fov'
# EBIRDST_DATA_DIR='D:/Working_Files/1_Projects/Landbirds/SK_BBA_analysis/SaskAtlas/!Data/!Spatial/eBird/'

# Spatial layers for SK
SK_BCR <- st_read("../../../Data/Spatial/National/BCR/BCR_Terrestrial_master.shp")  %>%
  subset(PROVINCE_S == "SASKATCHEWAN") %>%
  st_make_valid() %>%
  dplyr::select(BCR, PROVINCE_S)

SaskBoundary <- SK_BCR %>% st_union()

# List to contain sf objects for species ranges
species_ranges <- list()

for (sp_code in species_to_model$Species_Code_BSC){
  print(sp_code)
  if (file.exists("../Data_Cleaned/Spatial/eBird_ranges_SK.RDS")) species_ranges <- readRDS("../Data_Cleaned/Spatial/eBird_ranges_SK.RDS")
  
  # ------------------------------------------------
  # Check is species already in the list
  # ------------------------------------------------
  if (sp_code %in% names(species_ranges)) next
  
  # ------------------------------------------------
  # Download and trim ebird range for this species
  # ------------------------------------------------
  
  species_name = Sask_spcd$CommonName[which(Sask_spcd$spcd == sp_code)]
  species_label = Sask_spcd$Label[which(Sask_spcd$spcd == sp_code)]
  
  # Check if species is available
  check <- get_species(species_name)
  if (length(check)==0){
    print(paste0(sp_code, " not available in eBird"))
    next
  }
  
  if (length(check)>0 & is.na(check)){
    print(paste0(sp_code, " not available in eBird"))
    next
  }
  
  ebirdst_download(species_name, pattern = "range") 
  
  path <- get_species_path(species_name)
  
  range <- load_ranges(path, resolution = "lr")
  
  range <- range %>% subset(season %in% c("resident","breeding")) %>% 
    st_transform(.,crs = st_crs(SaskBoundary)) %>% 
    st_union() %>%
    st_crop(SaskBoundary)
  
  species_ranges[[sp_code]] <- range
  saveRDS(species_ranges,file="../Data_Cleaned/Spatial/eBird_ranges_SK.RDS")
  
  
} # close species loop


# ******************************************************************
# ******************************************************************
# Save
# ******************************************************************
# ******************************************************************

analysis_data_package <- list(
  
  all_surveys = all_surveys, # Survey information (protocol, location, time, etc)
  full_count_matrix = full_count_matrix, # counts of each species for each survey
  
  # Contains covariates on SK-wide grid
  SaskGrid = SaskGrid,
  
  # Species to include in analysis (and number of detections in point counts)
  species_to_model = species_to_model,
  
  # Species codes and common names
  Sask_spcd = Sask_spcd,
  
  # Covariates to include in analysis
  covariates_to_include = covariates_to_include,
  
  # eBird range limits
  species_ranges = species_ranges
  
)

saveRDS(analysis_data_package,"../Data_Cleaned/analysis_data_package.rds")
