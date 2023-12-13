# ************************************************************************
# Prepare spatial layers for analysis
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
             'magrittr',
             'circular',
             'stars',
             'corrplot',
             'ggpubr')

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

my.wt.circ.mean <- function(aspect, slope) {
  sinr <- sum(sin(rad(aspect)) * slope,na.rm=T)
  cosr <- sum(cos(rad(aspect)) * slope,na.rm=T)
  circmean <- atan2(sinr, cosr)
  return(deg(circmean) %% 360)
}


# ---------------------------------------------------
# CRS in which maps will be produced
# ---------------------------------------------------

target_crs <- "+proj=aea +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-106 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs "

# ---------------------------------------------------
# Boundary of Saskatchewan
# ---------------------------------------------------

SaskBoundary <- st_read("../../../Data/Spatial/National/BCR/BCR_Terrestrial_master.shp")  %>%
  subset(PROVINCE_S == "SASKATCHEWAN") %>%
  st_make_valid() %>%
  st_union() %>%
  st_transform(target_crs)

SaskBoundary_buffer <- st_buffer(SaskBoundary,10000)

# *******************************************************************
# *******************************************************************
# Relevant shapefiles
# *******************************************************************
# *******************************************************************

# Path to spatial covariates
covar_folder <- "../../../Data/Spatial/"

# Example raster with target properties
target_raster <- rast(paste0(covar_folder,"/National/AnnualMeanTemperature/wc2.1_30s_bio_1.tif")) %>% 
  crop(st_transform(SaskBoundary_buffer ,crs(.))) %>% 
  project(target_crs, res = 1000)

# ---------------------------------------------------
# Geographic covariates
# ---------------------------------------------------

# Elevation
elevation <- rast(paste0(covar_folder,"Saskatchewan/SaskElevation/Sask_Elevation.tif")) %>%  
  project(target_raster, align = TRUE, method = "bilinear") %>%
  resample(y = target_raster,method="bilinear")

# Slope
slope <- rast(paste0(covar_folder,"Saskatchewan/SaskSlope/SaskSlope.tif")) %>%  
  project(target_raster, align = TRUE, method = "bilinear")%>%
  resample(y = target_raster,method="bilinear")

# Aspect
aspect <- rast(paste0(covar_folder,"Saskatchewan/SaskAspect/SaskAspect.tif")) %>%  
  project(target_raster, align = TRUE, method = "bilinear")%>%
  resample(y = target_raster,method="bilinear")

# ---------------------------------------------------
# Vegetation / land cover / habitat
# ---------------------------------------------------

# LCC2020 - NOTE THAT THIS IS NOT SCALED UP TO 1KM X 1KM PIXELS
lcc2020 <- rast(paste0(covar_folder,"National/LandCoverCanada2020/landcover-2020-classification.tif")) %>% 
  crop(st_transform(SaskBoundary_buffer,crs(.)))

# Percent broadleaf
prcBroadLeaf <- rast(paste0(covar_folder,"Saskatchewan/SaskSCANFI/SK_and_prairies_SCANFI_sps_prcB_SW_2020_v1.tif")) %>%  
  project(target_raster, align = TRUE, method = "bilinear")%>%
  resample(y = target_raster,method="bilinear")

# Percent conifer
prcConifer <- rast(paste0(covar_folder,"Saskatchewan/SaskSCANFI/SK_and_prairies_SCANFI_sps_prcC_other_SW_2020_v1.tif")) %>%  
  project(target_raster, align = TRUE, method = "bilinear")%>%
  resample(y = target_raster,method="bilinear")

# Biomass
biomass <- rast(paste0(covar_folder,"Saskatchewan/SaskSCANFI/SK_and_prairies_SCANFI_att_biomass_SW_2020_v1.tif")) %>%  
  project(target_raster, align = TRUE, method = "bilinear")%>%
  resample(y = target_raster,method="bilinear")

# Canopy height
height <- rast(paste0(covar_folder,"Saskatchewan/SaskSCANFI/SK_and_prairies_SCANFI_att_height_SW_2020_v1.tif")) %>%  
  project(target_raster, align = TRUE, method = "bilinear")%>%
  resample(y = target_raster,method="bilinear")

# Stand canopy closure
closure <- rast(paste0(covar_folder,"Saskatchewan/SaskSCANFI/SK_and_prairies_SCANFI_att_closure_SW_2020_v1.tif")) %>%  
  project(target_raster, align = TRUE, method = "bilinear")%>%
  resample(y = target_raster,method="bilinear")

# Distance to H20
dist2water <- rast(paste0(covar_folder,"Saskatchewan/SaskDist2Water/Dist_to_H20resample.tif")) %>%  
  project(target_raster, align = TRUE, method = "bilinear")%>%
  resample(y = target_raster,method="bilinear")

# Wetland recurrence index
wetland_recurr <- rast(paste0(covar_folder,"Saskatchewan/SaskWetland/wetland_recurr_Sask.tif")) %>%  
  project(target_raster, align = TRUE, method = "bilinear")%>%
  resample(y = target_raster,method="bilinear")

# ---------------------------------------------------
# Climate
# ---------------------------------------------------

# Mean Annual Temperature
MAT <- rast(paste0(covar_folder,"National/Bioclimate/Normal_1991_2020_MAT.tif")) %>% 
  crop(st_transform(SaskBoundary_buffer,crs(.))) %>%
  project(target_raster, align = TRUE, method = "bilinear")%>%
  resample(y = target_raster,method="bilinear")

# Mean Annual Precipitation
MAP <- rast(paste0(covar_folder,"National/Bioclimate/Normal_1991_2020_MAP.tif")) %>% 
  crop(st_transform(SaskBoundary_buffer,crs(.))) %>%
  project(target_raster, align = TRUE, method = "bilinear")%>%
  resample(y = target_raster,method="bilinear")

# ---------------------------------------------------
# Combine into single raster stack
# ---------------------------------------------------

raster_list <- list(elevation,slope, aspect,
                    prcBroadLeaf,prcConifer,biomass,height,closure,dist2water,wetland_recurr,
                    MAT,MAP)
raster_stack <- rast(raster_list)
names(raster_stack) = c("elevation","slope","aspect","prcBroadLeaf","prcConifer","biomass","height","closure","dist2water","wetland_recurr",
                        "MAT","MAP")

plot(raster_stack, col = viridis(10))


# *******************************************************************
# *******************************************************************
# Create grid across study area
# *******************************************************************
# *******************************************************************

# Polygons for every 1km x 1km pixel
SaskGrid <- st_as_stars(raster_stack) %>% st_as_sf(., as_points = FALSE, merge = FALSE)

# Centroids of every pixel
SaskGrid_centroid <- st_as_stars(raster_stack) %>% st_as_sf(., as_points = TRUE, merge = FALSE)

# ---------------------------------------------------
# For each point in SaskGrid, extract covariates within 1 km of location
# ---------------------------------------------------

SaskGrid_1km <- SaskGrid_centroid %>% st_buffer(1000) 

# Proportion of each land cover class from lcc 2020
prop_LCC_1km <- exact_extract(lcc2020,SaskGrid_1km,"frac") %>% suppressWarnings()
names(prop_LCC_1km) <- paste0(str_replace(names(prop_LCC_1km),"frac","LCC"),"_1km")
prop_LCC_1km[setdiff(paste0("LCC_",seq(1,18),"_1km"),names(prop_LCC_1km))] <- 0
prop_LCC_1km %<>% dplyr::select(sort(names(.)))

# Combine land cover class names
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

SaskGrid <- bind_cols(SaskGrid,prop_LCC_1km)

saveRDS(SaskGrid,file="../Data_Cleaned/Spatial/SaskGrid.RDS")


# *******************************************************************
# *******************************************************************
# Extract covariates at each survey location
# *******************************************************************
# *******************************************************************

analysis_data <- readRDS(file = "../Data_Cleaned/analysis_data.rds")
all_surveys <- analysis_data$all_surveys
full_count_matrix <- analysis_data$full_count_matrix

# ------------------------------------------
# Atlas Squares in which each survey is located
# ------------------------------------------

SaskSquares <- st_read("../../../Data/Spatial/Saskatchewan/SaskSquares/SaskSquares.shp") %>%
  st_transform(st_crs(all_surveys)) %>%
  dplyr::select(SQUARE_ID) %>%
  rename(sq_id = SQUARE_ID)

all_surveys <- all_surveys %>% 
  mutate(Obs_Index = 1:nrow(.)) %>%
  st_intersection(SaskSquares) %>%
  arrange(Obs_Index)
full_count_matrix <- full_count_matrix[all_surveys$Obs_Index,]

all_surveys$Obs_Index <- 1:nrow(all_surveys)

# ---------------------------------------------------
# Extract covariates for each survey from raster stack
# ---------------------------------------------------

all_surveys = terra::extract(raster_stack,vect(all_surveys) , bind = TRUE) %>% st_as_sf()

all_surveys_1km <- all_surveys %>% st_buffer(1000)

# Proportion of each land cover class from lcc2020
prop_LCC_1km <- exact_extract(lcc2020,all_surveys_1km,"frac") %>% suppressWarnings()
names(prop_LCC_1km) <- paste0(str_replace(names(prop_LCC_1km),"frac","LCC"),"_1km")
prop_LCC_1km[setdiff(paste0("LCC_",seq(1,18),"_1km"),names(prop_LCC_1km))] <- 0
prop_LCC_1km %<>% dplyr::select(sort(names(.)))

# Combine land cover classes
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

all_surveys <- bind_cols(all_surveys,prop_LCC_1km)

# Remove several surveys with missing covariate information
all_surveys <- all_surveys %>%
  mutate(Obs_Index = 1:nrow(all_surveys)) %>%
  drop_na(elevation:Water_1km)

full_count_matrix <- full_count_matrix[all_surveys$Obs_Index,]
all_surveys$Obs_Index <- 1:nrow(all_surveys)

# *******************************************************************
# *******************************************************************
# Evaluate correlations between covariates
# *******************************************************************
# *******************************************************************

M = as.data.frame(all_surveys) %>%
  select(elevation:Water_1km,-aspect) %>%
  cor()

corrplot(M, method = 'number', type = "upper", order = 'hclust',
         addrect = 3, rect.col = 'blue', rect.lwd = 3, tl.pos = 'd')

# height, biomass, and closure are extremely tightly related.  

covars_for_PCA <- c("Urban_1km",
                    "Crop_1km",
                    "elevation",
                    "MAT",
                    "MAP",
                    "prcConifer",
                    "prcBroadLeaf",
                    "closure",
                    "dist2water",
                    "Needleleaf_forest_1km",
                    "Wetland_1km",
                    "wetland_recurr",
                    "slope",
                    "Grass_shrub_1km")

M = as.data.frame(all_surveys) %>%
  select(covars_for_PCA) %>%
  cor()

corrplot(M, method = 'number', type = "upper", order = 'hclust',
         addrect = 3, rect.col = 'blue', rect.lwd = 3, tl.pos = 'd')


# *******************************************************************
# *******************************************************************
# Conduct principal components analysis
# *******************************************************************
# *******************************************************************

dat_for_PCA <- all_surveys %>%
  as.data.frame() %>%
  dplyr::select(covars_for_PCA)

pca <- prcomp(dat_for_PCA, scale = TRUE)

# ------------------------------------------
# Interpretation of specific axes (e.g., axes 1 and 2)
# ------------------------------------------

summary(pca)   # Proportion variance explaind by axes
fviz_eig(pca)  # Scree plot (first 5 axes explain 85% of variation in habitat between sites)
pca            # Variable loadings

fviz_pca_var(pca,
             axes = c(1,12),
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = viridis(10),
             repel = TRUE     # Avoid text overlapping
)

# ------------------------------------------
# Predict PCA values for each survey location and standardize (mean = 0, sd = 1)
# ------------------------------------------

all_surveys_PCA <- predict(pca, newdata = as.data.frame(all_surveys)[,names(pca$center)])
SaskGrid_PCA <- predict(pca, newdata = as.data.frame(SaskGrid)[,names(pca$center)])

for (covar in colnames(all_surveys_PCA)){

  covar_mean <- mean(all_surveys_PCA[,covar],na.rm = TRUE)
  covar_sd <- sd(all_surveys_PCA[,covar],na.rm = TRUE)

  all_surveys_PCA[,covar] <- (as.data.frame(all_surveys_PCA)[,covar] - covar_mean)/covar_sd
  SaskGrid_PCA[,covar] <- (as.data.frame(SaskGrid_PCA)[,covar] - covar_mean)/covar_sd
  
}

all_surveys <- all_surveys %>% bind_cols(all_surveys_PCA)
SaskGrid <- bind_cols(SaskGrid,SaskGrid_PCA)

# ------------------------------------------
# Plot maps of covariates
# ------------------------------------------

SaskRast <- SaskGrid %>% 
  dplyr::select(covars_for_PCA,PC1:PC14) %>%
  stars::st_rasterize()

covar_to_plot <- names(SaskRast)

covar_plotlist <- list()

for (covar in covar_to_plot){
  cplot <- ggplot() + geom_stars(data = SaskRast, aes(fill = !!sym(covar)))+
    scale_fill_gradientn(colours = viridis(10), name = covar)+
    geom_sf(data = SaskBoundary,colour="black",fill=NA,lwd=0.3,show.legend = F)+
    ggtitle(covar)+
    theme_bw()+
    xlab("")+ylab("")
  
  covar_plotlist[[covar]] <- cplot
}

covar_plots <- ggarrange(plotlist = covar_plotlist,nrow=2,ncol=length(covars_for_PCA))

png("../Output/Covariate_Maps/Covariate_Maps.png", width=60, height=10, units="in", res=300, type="cairo")
print(covar_plots)
dev.off()


# ******************************************************************
# Prepare species names
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
# Summary of number of detections (in point counts/ARUs) for each species
# ------------------------------------------

BSC_species <- readRDS("../Data_Cleaned/BSC_species.RDS")

PC_surveys <- which(all_surveys$Survey_Type %in% c("Point_Count","ARU_SPM","ARU_SPT"))
n_detections <- full_count_matrix
rownames(n_detections) <- all_surveys$sq_id
n_detections <- n_detections[PC_surveys,]

n_detections <- n_detections %>% 
  reshape2::melt() %>%
  dplyr::rename(sq_id = Var1, Species_Code_BSC = Var2, detected = value) %>%
  subset(detected>0) %>%
  group_by(Species_Code_BSC) %>%
  summarize(n_squares = length(unique(sq_id)),
            n_detections = sum(detected>0)) 

species_to_model <- left_join(n_detections,BSC_species %>% dplyr::select(-index),
                              by = c("Species_Code_BSC" = "BSC_spcd")) %>%
  
  # Remove erroneous observations
  subset(english_name %!in% c("Accipiter sp." ,
                              "passerine sp.",
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
                              "Buteo sp.",
                              "Bat sp.",
                              "woodpecker sp.",
                              "falcon sp.",
                              "finch sp.",
                              "owl sp.",
                              "peep sp.",
                              "swan sp.",
                              "tern sp.",
                              "wren sp.",
                              "new world flycatcher sp."))

dim(species_to_model)
species_to_model$english_name %>% sort()

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
# EBIRDST_DATA_DIR='C:/Users/IlesD/OneDrive - EC-EC/Iles/Projects/Landbirds/Landbird-Distribution-Modeling-ECCC/Data/Spatial/eBird/'

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
  
  pca = pca,
  
  # Contains covariates on SK-wide grid
  SaskGrid = SaskGrid,
  
  # Species to include in analysis (and number of detections in point counts)
  species_to_model = species_to_model,
  
  # Species codes and common names
  Sask_spcd = Sask_spcd,
  
  # eBird range limits
  species_ranges = species_ranges
  
)

saveRDS(analysis_data_package,"../Data_Cleaned/analysis_data_package.rds")

write.csv(analysis_data_package$species_to_model,file = "../Data_Cleaned/species_to_model.csv",row.names = FALSE)
