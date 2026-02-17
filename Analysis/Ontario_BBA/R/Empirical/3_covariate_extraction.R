# ============================================================================
# DEFINE STUDY AREA
# EXTRACT COVARIATES FOR EACH SURVEY LOCATION AND GRID ACROSS STUDY AREA
# ============================================================================

# ============================================================================
# Load/install packages
# ============================================================================

my_packs <- c('tidyverse', 'sf','exactextractr','terra','factoextra','viridis','stars','ebirdst','readxl','dbscan')

if (any(!my_packs %in% installed.packages()[, 'Package'])){
  install.packages(my_packs[which(!my_packs %in% installed.packages()[, 'Package'])], dependencies = TRUE)
}
lapply(my_packs, require, character.only = TRUE)

rm(list = ls())

setwd("C:/Users/IlesD/OneDrive - EC-EC/Iles/Projects/Landbirds/Landbird-Distribution-Modeling-ECCC/")

`%!in%` <- Negate(`%in%`)

# ============================================================================
# Load data package prepared by script 1
# ============================================================================

analysis_data <- readRDS(file = "Analysis/Ontario_BBA/Data_Cleaned/analysis_data.rds")
all_surveys <- analysis_data$NC_surveyinfo
count_matrix <- analysis_data$full_count_matrix
all_species <- analysis_data$all_species
rm(analysis_data)

AEA_proj <- st_crs(all_surveys)

# Prepare study boundary, which includes BCRs 13, 12, and lower 100km of BCR 8
focal_area <- st_read("Data/Spatial/BCR/BCR_Terrestrial_master.shp") %>%
  subset(PROVINCE_S == "ONTARIO") %>%
  st_make_valid() %>%
  dplyr::select(BCR, PROVINCE_S) %>%
  st_transform(AEA_proj)%>%
  st_union() %>%
  st_set_precision(1) %>%
  st_buffer(100000)

study_boundary <- st_read("Data/Spatial/BCR/BCR_Terrestrial_master.shp") %>%
  subset(PROVINCE_S == "ONTARIO") %>%
  #st_make_valid() %>%
  dplyr::select(BCR, PROVINCE_S) %>%
  st_transform(AEA_proj)%>%
  st_union() %>%
  st_set_precision(1) %>%
  st_make_valid() 

# Subset to surveys occurring within the study boundary
inside <- st_within(all_surveys, study_boundary, sparse = FALSE) %>%
  as.logical()
all_surveys <- all_surveys[which(inside), ]
count_matrix <- count_matrix[which(inside), ]

# ============================================================================
# Load covariate layers
# ============================================================================

# Path to spatial covariates
covar_folder <- "Analysis/Ontario_BBA/Data_Cleaned/Spatial/"

# Add 10 km buffer so that covariates can extend outside of study area
study_boundary_buffer <- study_boundary %>% st_buffer(10000)

# ---- Time-invariant layers
water_open <- terra::rast(paste0(covar_folder,"water_open.tif")) %>% crop(st_transform(study_boundary_buffer,crs(.))) 
water_river <- terra::rast(paste0(covar_folder,"water_river.tif")) %>% crop(st_transform(study_boundary_buffer,crs(.))) 

# ---- Layers for OBBA2
prec_OBBA2 <- terra::rast(paste0(covar_folder,"prec_OBBA2.tif")) %>% crop(st_transform(study_boundary_buffer,crs(.)))
tmax_OBBA2 <- terra::rast(paste0(covar_folder,"tmax_OBBA2.tif")) %>% crop(st_transform(study_boundary_buffer,crs(.)))
urban_OBBA2 <- terra::rast(paste0(covar_folder,"Urban_2000.tif")) %>% crop(st_transform(study_boundary_buffer,crs(.)))
road_OBBA2 <- terra::rast(paste0(covar_folder,"roadside_2005_100m.tif")) %>% crop(st_transform(study_boundary_buffer,crs(.))) 
insect_OBBA2 <- terra::rast(paste0(covar_folder,"insect_OBBA2.tif")) %>% crop(st_transform(study_boundary_buffer,crs(.))) 
lc_OBBA2 <- terra::rast(paste0(covar_folder,"MCD12Q1_OBBA2_mode.tif")) %>% crop(st_transform(study_boundary_buffer,crs(.)))

# ---- Layers for OBBA3
prec_OBBA3 <- terra::rast(paste0(covar_folder,"prec_OBBA3.tif")) %>% crop(st_transform(study_boundary_buffer,crs(.)))
tmax_OBBA3 <- terra::rast(paste0(covar_folder,"tmax_OBBA3.tif")) %>% crop(st_transform(study_boundary_buffer,crs(.)))
urban_OBBA3 <- terra::rast(paste0(covar_folder,"Urban_2020.tif")) %>% crop(st_transform(study_boundary_buffer,crs(.)))
road_OBBA3 <- terra::rast(paste0(covar_folder,"roadside_2025_100m.tif")) %>% crop(st_transform(study_boundary_buffer,crs(.))) 
insect_OBBA3 <- terra::rast(paste0(covar_folder,"insect_OBBA3.tif")) %>% crop(st_transform(study_boundary_buffer,crs(.))) 
lc_OBBA3 <- terra::rast(paste0(covar_folder,"MCD12Q1_OBBA3_mode.tif")) %>% crop(st_transform(study_boundary_buffer,crs(.)))

# Names of land cover classes from modis
lc_classes <- readxl::read_xlsx("Data/Spatial/MODIS/modis_lc_classes.xlsx")

# ============================================================================
# Prepare a 1 km grid on which maps/predictions will eventually be made
# ============================================================================

set.seed(123)

grid_resolution_m <- 1000
# Create grid (cellsize is in km)
grid <- st_make_grid(
  study_boundary,
  cellsize = units::set_units(grid_resolution_m, m),
  what = "polygons",
  square = TRUE,
  flat_topped = FALSE
)

grid_sf <- st_as_sf(grid) %>% rename(geometry = x)

# Filter grid to those intersecting the boundary
grid_sf <- st_filter(grid_sf, study_boundary, .predicate = st_intersects)
grid_sf <- grid_sf %>% st_centroid()

# ============================================================================
# For each grid cell, extract covariates from OBBA2
# ============================================================================

# ---- Continuous covariates (calculate landscape mean values)
grid_OBBA2 <- grid_sf %>% st_buffer(grid_resolution_m/2)

grid_OBBA2$water_river <- exact_extract(water_river, grid_OBBA2, 'mean')
grid_OBBA2$prec <- exact_extract(prec_OBBA2, grid_OBBA2, 'mean')
grid_OBBA2$tmax <- exact_extract(tmax_OBBA2, grid_OBBA2, 'mean')
grid_OBBA2$road <- exact_extract(road_OBBA2, grid_OBBA2, 'mean')

# ---- Categorical covariates (calculate landscape fractions within buffer)
# Proportion of MODIS land cover classes
prop_lc_OBBA2 <- exact_extract(lc_OBBA2,st_transform(grid_OBBA2,st_crs(lc_OBBA2)),"frac") %>% suppressWarnings()
names(prop_lc_OBBA2) <- str_replace(names(prop_lc_OBBA2),"frac","lc")
grid_OBBA2 <- bind_cols(grid_OBBA2,prop_lc_OBBA2)

# Proportion of landscape comprised by each of 3 urban classes
prop_urban_OBBA2 <- exact_extract(urban_OBBA2,st_transform(grid_OBBA2,st_crs(urban_OBBA2)),"frac") %>% suppressWarnings()
names(prop_urban_OBBA2) <- str_replace(names(prop_urban_OBBA2),"frac","urban")
grid_OBBA2 <- bind_cols(grid_OBBA2,prop_urban_OBBA2)

# Proportion of landscape comprised by each insect class
r_broadleaf_OBBA2 <- max(insect_OBBA2[[c("Forest Tent Caterpillar", "Gypsy Moth")]])
r_needleleaf_OBBA2 <- max(insect_OBBA2[[c("Jack Pine Budworm", "Spruce Budworm")]])
insects_OBBA2 <- c(r_broadleaf_OBBA2, r_needleleaf_OBBA2)
names(insects_OBBA2) <- c("insect_broadleaf", "insect_needleleaf")
prop_insect_OBBA2 <- exact_extract(insects_OBBA2,st_transform(grid_OBBA2,st_crs(insect_OBBA2)),"frac") %>% suppressWarnings()
prop_insect_OBBA2 <- prop_insect_OBBA2 %>% dplyr::select(-starts_with("frac_0."))
names(prop_insect_OBBA2) <- str_replace(names(prop_insect_OBBA2),"frac_1.","")

# Combine insect columns
grid_OBBA2 <- bind_cols(grid_OBBA2,prop_insect_OBBA2)

# ============================================================================
# For each grid cell, extract covariates from OBBA3
# ============================================================================

# ---- Continuous covariates (calculate landscape mean values)
grid_OBBA3 <- grid_OBBA2 %>% dplyr::select(water_river)
grid_OBBA3$prec <- exact_extract(prec_OBBA3, grid_OBBA3, 'mean')
grid_OBBA3$tmax <- exact_extract(tmax_OBBA3, grid_OBBA3, 'mean')
grid_OBBA3$road <- exact_extract(road_OBBA3, grid_OBBA3, 'mean')

# ---- Categorical covariates (calculate landscape fractions within buffer)
# Proportion of MODIS land cover classes
prop_lc_OBBA3 <- exact_extract(lc_OBBA3,st_transform(grid_OBBA3,st_crs(lc_OBBA3)),"frac") %>% suppressWarnings()
names(prop_lc_OBBA3) <- str_replace(names(prop_lc_OBBA3),"frac","lc")
grid_OBBA3 <- bind_cols(grid_OBBA3,prop_lc_OBBA3)

# Proportion of landscape comprised by each of 3 urban classes
prop_urban_OBBA3 <- exact_extract(urban_OBBA3,st_transform(grid_OBBA3,st_crs(urban_OBBA3)),"frac") %>% suppressWarnings()
names(prop_urban_OBBA3) <- str_replace(names(prop_urban_OBBA3),"frac","urban")
grid_OBBA3 <- bind_cols(grid_OBBA3,prop_urban_OBBA3)

# Proportion of landscape comprised by each insect class
r_broadleaf_OBBA3 <- max(insect_OBBA3[[c("Forest Tent Caterpillar", "Gypsy Moth")]])
r_needleleaf_OBBA3 <- max(insect_OBBA3[[c("Jack Pine Budworm", "Spruce Budworm")]])
insects_OBBA3 <- c(r_broadleaf_OBBA3, r_needleleaf_OBBA3)
names(insects_OBBA3) <- c("insect_broadleaf", "insect_needleleaf")
prop_insect_OBBA3 <- exact_extract(insects_OBBA3,st_transform(grid_OBBA3,st_crs(insect_OBBA3)),"frac") %>% suppressWarnings()
prop_insect_OBBA3 <- prop_insect_OBBA3 %>% dplyr::select(-starts_with("frac_0."))
names(prop_insect_OBBA3) <- str_replace(names(prop_insect_OBBA3),"frac_1.","")

# Combine insect columns
grid_OBBA3 <- bind_cols(grid_OBBA3,prop_insect_OBBA3)

# ============================================================================
# Examine plots of covariates
# ============================================================================

# ---- Plot for OBBA2
df_OBBA2 <- grid_OBBA2 %>%
  st_centroid() %>%
  bind_cols(st_coordinates(.)) %>%
  as.data.frame() %>%
  dplyr::select(-geometry) %>%
  tidyr::pivot_longer(cols = -c(X,Y), names_to = "variable", values_to = "value") %>%
  group_by(variable) %>%
  mutate(value_scaled = (value - min(value, na.rm = TRUE)) /
           (max(value, na.rm = TRUE) - min(value, na.rm = TRUE))) %>%
  ungroup()

covar_plot_OBBA2 <- ggplot(df_OBBA2, aes(x = X, y = Y, fill = value_scaled))+
  geom_tile()+
  facet_wrap(~variable, scales = "free")+
  scale_fill_viridis_c() +
  theme_minimal()

# Save the plot as PDF
ggsave(
  filename = "Analysis/Ontario_BBA/Output/Covariate_Maps/covariate_maps_OBBA2.pdf",
  plot = covar_plot_OBBA2,
  width = 20,        # width in inches
  height = 10        # height in inches
)

# ---- Plot for OBBA3

df_OBBA3 <- grid_OBBA3 %>%
  st_centroid() %>%
  bind_cols(st_coordinates(.)) %>%
  as.data.frame() %>%
  dplyr::select(-geometry) %>%
  tidyr::pivot_longer(cols = -c(X,Y), names_to = "variable", values_to = "value") %>%
  group_by(variable) %>%
  mutate(value_scaled = (value - min(value, na.rm = TRUE)) /
           (max(value, na.rm = TRUE) - min(value, na.rm = TRUE))) %>%
  ungroup()

covar_plot_OBBA3 <- ggplot(df_OBBA3, aes(x = X, y = Y, fill = value_scaled))+
  geom_tile()+
  facet_wrap(~variable, scales = "free")+
  scale_fill_viridis_c() +
  theme_minimal()

# Save the plot as PDF
ggsave(
  filename = "Analysis/Ontario_BBA/Output/Covariate_Maps/covariate_maps_OBBA3.pdf",
  plot = covar_plot_OBBA3,
  width = 20,        # width in inches
  height = 10        # height in inches
)

# ============================================================================
# Examine correlations between covariates
# ============================================================================

# DECISION TO REMOVE THE FOLLOWING COVARIATES:
# "lc_13",  # MODIS-based classification of urban habitat, which is already captured by urban_3 from GHSL,
# "urban_1", # GHSL classified as "unpopulated", which is the converse of urban_3
# "water_open" # based on Ontario Hydrology layer; not time-varying, and highly correlated with MODIS lc_17

# Both atlases
pca_input <- bind_rows(grid_OBBA2,grid_OBBA3) %>%
  as.data.frame() %>%
  # Remove geometry and several covariates that are either rare or captured by other covariates
  dplyr::select(-geometry, -lc_3, -lc_7, -lc_16,-lc_13,-urban_1)

# Compute correlation matrix
cor_matrix <- cor(pca_input, use = "pairwise.complete.obs")

# Threshold for "high correlation"
thresh <- 0.5  
sum( abs(cor_matrix[upper.tri(cor_matrix)]) >= thresh )

# Get upper triangle only
cor_pairs <- which(abs(cor_matrix) > thresh & upper.tri(cor_matrix), arr.ind = TRUE)

# Make a tidy dataframe of correlated pairs
high_cor_df <- data.frame(
  var1 = rownames(cor_matrix)[cor_pairs[,1]],
  var2 = colnames(cor_matrix)[cor_pairs[,2]],
  correlation = cor_matrix[cor_pairs]
) %>%
  arrange(desc(abs(correlation)))

high_cor_df


# ============================================================================
# Assign covariates from grid_sf to surveys, referenced to the correct time point
# ============================================================================

all_surveys$obs_idx <- 1:nrow(all_surveys)

# Split surveys by Atlas
surveys_OBBA2 <- all_surveys %>% filter(Atlas == "OBBA2")
surveys_OBBA3 <- all_surveys %>% filter(Atlas == "OBBA3")

# Nearest join to corresponding grid
nearest_OBBA2 <- st_nearest_feature(surveys_OBBA2, grid_OBBA2)
nearest_OBBA3 <- st_nearest_feature(surveys_OBBA3, grid_OBBA3)

# Bind covariates
surveys_OBBA2_covs <- surveys_OBBA2 %>%
  bind_cols(st_drop_geometry(grid_OBBA2[nearest_OBBA2, ]))

surveys_OBBA3_covs <- surveys_OBBA3 %>%
  bind_cols(st_drop_geometry(grid_OBBA3[nearest_OBBA3, ]))

# Combine back
all_surveys_with_covs <- bind_rows(surveys_OBBA2_covs, surveys_OBBA3_covs) %>%
  arrange(obs_idx)
rownames(all_surveys_with_covs) <- NULL

# ============================================================================
# Download and store species range limits (from eBird)
# ============================================================================

# Limit to species actually encountered in the atlas dataset
count_matrix <- count_matrix[,colSums(count_matrix)!= 0]
species_encountered <- colnames(count_matrix)

all_species <- all_species %>%
  mutate(species_id = as.character(species_id)) %>%
  subset(species_id %in% species_encountered)

# List to contain sf objects for species ranges
species_ranges <- list()
if (file.exists("Analysis/Ontario_BBA/Data_Cleaned/Spatial/eBird_ranges_ON.rds")) species_ranges <- readRDS("Analysis/Ontario_BBA/Data_Cleaned/Spatial/eBird_ranges_ON.rds")

for (species_name in unique(all_species$english_name)){

  print(species_name)

  if (species_name %in% names(species_ranges)) next
  
  if (file.exists("Analysis/Ontario_BBA/Data_Cleaned/Spatial/eBird_ranges_ON.rds")) species_ranges <- readRDS("Analysis/Ontario_BBA/Data_Cleaned/Spatial/eBird_ranges_ON.rds")

  # ------------------------------------------------
  # Check is species already in the list
  # ------------------------------------------------
  
  if (species_name %in% names(species_ranges)) next

  # ------------------------------------------------
  # Download and trim ebird range for this species
  # ------------------------------------------------
  
  sp_code_ebird <- subset(ebirdst::ebirdst_runs, common_name == species_name)$species_code
  if (species_name == "Yellow-bellied Sapsucker") sp_code_ebird <- sp_code_ebird[2]

  # Check if species is available
  check <- get_species(sp_code_ebird)
  if (length(check)==0){
    print(paste0(species_name, " not available in eBird"))
    next
  }

  if (length(check)>0 & is.na(check)){
    print(paste0(species_name, " not available in eBird"))
    next
  }

  ebirdst_download_status(sp_code_ebird,
                          download_abundance = FALSE,
                          download_ranges = TRUE,
                          pattern = "_smooth_9km_")

  path <- get_species_path(sp_code_ebird)

  range <- load_ranges(sp_code_ebird, resolution = "9km",smoothed = TRUE)

  range <- range %>% subset(season %in% c("resident","breeding")) %>%
    st_make_valid() %>%
    st_union()

  species_ranges[[species_name]] <- range
  saveRDS(species_ranges,file="Analysis/Ontario_BBA/Data_Cleaned/Spatial/eBird_ranges_ON.rds")

} # close species loop

# ============================================================================
# Append BCR to each survey
# ============================================================================

# BCRs in Ontario
BCR <- st_read("Data/Spatial/BCR/BCR_Terrestrial_master.shp") %>%
  subset(PROVINCE_S == "ONTARIO") %>%
  st_make_valid() %>%
  group_by(BCR,BCRNAME) %>%
  summarize() %>%
  st_transform(AEA_proj)

# Determine BCR membership for each survey
all_surveys_with_covs <- all_surveys_with_covs %>%
  st_join(BCR %>% select(BCR), left = TRUE)  # keep only the BCR column

# Identify surveys that are still outside polygons and associate with nearest BCR
missing_BCR <- which(is.na(all_surveys_with_covs$BCR))

if(length(missing_BCR) > 0){
  # Find nearest BCR polygon for these points
  nearest_idx <- st_nearest_feature(all_surveys_with_covs[missing_BCR, ], BCR)
  
  # Assign BCR from nearest polygon
  all_surveys_with_covs$BCR[missing_BCR] <- BCR$BCR[nearest_idx]
}

# Determine BCR membership for each location in grid_sf
grid_sf <- grid_sf %>%
  st_join(BCR %>% select(BCR), left = TRUE)  # keep only the BCR column

# Identify surveys that are still outside polygons and associate with nearest BCR
missing_BCR <- which(is.na(grid_sf$BCR))

if(length(missing_BCR) > 0){
  # Find nearest BCR polygon for these points
  nearest_idx <- st_nearest_feature(grid_sf[missing_BCR, ], BCR)
  
  # Assign BCR from nearest polygon
  grid_sf$BCR[missing_BCR] <- BCR$BCR[nearest_idx]
}

# ============================================================================
# Calculate day of year for each survey
# ============================================================================

all_surveys_with_covs <- all_surveys_with_covs %>%
  mutate(DayOfYear = yday(Date_Time))

# ============================================================================
# Determine atlas square in which each survey is located
# ============================================================================

# Load shapefile containing atlas squares (10km x 10km)
atlas_squares <- st_read("Data/Spatial/National/AtlasSquares/NationalSquares_FINAL.shp")  %>%
  subset(prov == "ON") %>%
  st_transform(AEA_proj)

# Create square_id identifier
all_surveys_with_covs <- all_surveys_with_covs %>% 
  st_join(atlas_squares %>% select(square_id), left = TRUE) %>%
  relocate(geometry,.after = last_col())

# ============================================================================
# Save objects
# ============================================================================

analysis_data_covariates <- list(
  all_surveys_with_covs = all_surveys_with_covs,
  count_matrix = count_matrix,
  grid_OBBA2 = grid_OBBA2,
  grid_OBBA3 = grid_OBBA3,
  all_species = all_species,
  species_ranges = species_ranges,
  study_boundary = study_boundary,
  date_created = Sys.time()
)

saveRDS(analysis_data_covariates, file = "Analysis/Ontario_BBA/Data_Cleaned/analysis_data_covariates.rds")
