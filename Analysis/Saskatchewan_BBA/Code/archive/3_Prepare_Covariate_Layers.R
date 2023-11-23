# ------------------------------------------------
# Load/install packages and set graphical themes / working directory
# ------------------------------------------------

my_packs = c('tidyverse',
             'sf',
             'terra',
             'exactextractr',
             'magrittr')

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
# Load study area boundary
# ------------------------------------------------

AEA_proj <- "+proj=aea +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-106 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs "

SaskBoundary <- st_read("../../../Data/Spatial/National/BCR/BCR_Terrestrial_master.shp")  %>%
  subset(PROVINCE_S == "SASKATCHEWAN") %>%
  st_make_valid() %>%
  dplyr::select(BCR, PROVINCE_S) %>%
  st_transform(st_crs(AEA_proj))

# Add 20 km buffer so that covariates can extend slightly outside province
SaskBoundary_buffer <- SaskBoundary %>% st_buffer(20000)

# ------------------------------------------------
# Create fine-scale grid across study area (e.g., 1 km x 1 km pixels)
# ------------------------------------------------

SaskGrid <- st_make_grid(
  SaskBoundary,
  cellsize = units::set_units(1*1,km^2),
  what = "polygons",
  square = TRUE,
  flat_topped = FALSE)%>%
  st_as_sf() %>%
  st_intersection(SaskBoundary) %>%
  na.omit()

SaskGrid$point_id <- 1:nrow(SaskGrid)
SaskGrid_centroid <- st_centroid(SaskGrid)
  
# --------------------------------------------------------
# Load and reproject covariate rasters
#   - note that landcover raster is too large to reproject
# --------------------------------------------------------

# Path to spatial covariates
covar_folder <- "../../../Data/Spatial/"

# Elevation
elevation <- rast(paste0(covar_folder,"Saskatchewan/SaskElevation/Sask_Elevation.tif")) %>%  
  project(crs(SaskGrid), res = 250)

# Annual mean temperature
AMT <- rast(paste0(covar_folder,"National/AnnualMeanTemperature/wc2.1_30s_bio_1.tif")) %>% crop(st_transform(SaskBoundary_buffer,crs(.))) %>% project(crs(SaskGrid), res = 250)

# Stand canopy closure
SCC <- rast(paste0(covar_folder,"National/NationalForestInventory/NFI_MODIS250m_2011_kNN_Structure_Stand_CrownClosure_v1.tif")) %>% crop(st_transform(SaskBoundary_buffer,crs(.))) %>% project(crs(SaskGrid), res = 250)

# Land cover of Canada 2020 (not reprojected... faster to reproject grid)
lcc2020 <- rast(paste0(covar_folder,"National/LandCoverCanada2020/landcover-2020-classification.tif")) %>% crop(st_transform(SaskBoundary_buffer,crs(.))) 

# --------------------------------------------------------
# Extract covariate values within 1 km of each grid cell centroid
# --------------------------------------------------------

# 1 km circular buffer around each centroid
SaskGrid_1km <- st_buffer(SaskGrid_centroid, 1000)

# Continuous covariates
SaskGrid$elevation_1km <- exact_extract(elevation, SaskGrid_1km, 'mean')
SaskGrid$AMT_1km <- exact_extract(AMT, SaskGrid_1km, 'mean')
SaskGrid$SCC_1km <- exact_extract(SCC, SaskGrid_1km, 'mean')

# Proportion land cover within buffer
prop_LCC_1km <- exact_extract(lcc2020,SaskGrid_1km,"frac") %>% suppressWarnings()

names(prop_LCC_1km) <- paste0(str_replace(names(prop_LCC_1km),"frac","LCC"),"_1km")
prop_LCC_1km[setdiff(paste0("LCC_",seq(1,18),"_1km"),names(prop_LCC_1km))] <- 0
prop_LCC_1km %<>% 
  select(sort(names(.)))

# Join with SaskGrid
prop_LCC_1km$point_id <- SaskGrid$point_id
SaskGrid <- left_join(SaskGrid,prop_LCC_1km)

# Remove large objects to free up memory
rm(prop_LCC_1km,SaskGrid_1km)

# # --------------------------------------------------------
# # Extract covariate values within 5 km of each grid cell centroid
# # --------------------------------------------------------
# 
# # 5 km circular buffer around each centroid
# SaskGrid_5km <- st_buffer(SaskGrid_centroid, 5000)
# 
# # Continuous covariates
# SaskGrid$elevation_5km <- exact_extract(elevation, SaskGrid_5km, 'mean')
# SaskGrid$AMT_5km <- exact_extract(AMT, SaskGrid_5km, 'mean')
# SaskGrid$SCC_5km <- exact_extract(SCC, SaskGrid_5km, 'mean')
# 
# # Proportion land cover within buffer
# prop_LCC_matrix_5km <- matrix(NA,nrow = nrow(SaskGrid),ncol = 18) # Empty matrix to store results
# 
# # Modified for faster and less memory intensive way of getting proportion of LCC
# prop_LCC_5km <- exact_extract(lcc2020,SaskGrid_5km,"frac") %>% suppressWarnings()
# 
# names(prop_LCC_5km) <- paste0(str_replace(names(prop_LCC_5km),"frac","LCC"),"_5km")
# 
# prop_LCC_5km[setdiff(paste0("LCC_",seq(1,18),"_5km"),names(prop_LCC_5km))] <- 0
# prop_LCC_5km %<>% 
#   select(sort(names(.)))
# 
# # Join with SaskGrid
# prop_LCC_5km$point_id <- SaskGrid$point_id
# SaskGrid <- left_join(SaskGrid,prop_LCC_5km)
# 
# # Remove large objects to free up memory
# rm(prop_LCC_5km)

# --------------------------------------------------------
# Save
# --------------------------------------------------------

SaskGrid <- SaskGrid %>% relocate(point_id)

saveRDS(SaskGrid,"../Data_Cleaned/Spatial/SaskGrid_covariates.rds")
