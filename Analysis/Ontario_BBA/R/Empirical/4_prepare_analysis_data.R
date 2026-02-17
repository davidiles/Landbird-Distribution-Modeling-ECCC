library(tidyverse)
library(sf)
library(readxl)
library(lubridate)

# Clear working memory
rm(list=ls())

set.seed(123)
setwd("C:/Users/IlesD/OneDrive - EC-EC/Iles/Projects/Landbirds/Landbird-Distribution-Modeling-ECCC/")

# ==========================
# Load raw/processed data
# ==========================

analysis_data_covariates <- readRDS("Analysis/Ontario_BBA/Data_Cleaned/analysis_data_covariates.rds")

all_species <- analysis_data_covariates$all_species
species_ranges <- analysis_data_covariates$species_ranges
all_surveys <- analysis_data_covariates$all_surveys_with_covs
counts <- analysis_data_covariates$count_matrix
grid_OBBA2 <- analysis_data_covariates$grid_OBBA2
grid_OBBA3 <- analysis_data_covariates$grid_OBBA3

# Note units in km, which helps INLA fit spatial fields more easily
AEA_proj  <- "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=49 +lon_0=-95 +datum=NAD83 +units=km +no_defs"

# ==========================
# Load shapefiles
# ==========================

# Study boundary
study_boundary <- analysis_data_covariates$study_boundary %>%
  st_transform(st_crs(AEA_proj))

# All BCRs
allBCR <- st_read("Data/Spatial/BCR/BCR_Terrestrial_master.shp") %>%
  st_transform(st_crs(study_boundary)) %>%
  st_make_valid() %>%
  st_buffer(0) %>%
  group_by(BCR,BCRNAME,PROVINCE_S) %>%
  summarise() %>%
  st_transform(st_crs(AEA_proj))

# Waterbodies
water <- st_read("Analysis/Ontario_BBA/Data_Cleaned/Spatial/water_filtered.shp") %>%
  select(WATERBODY_, geometry) %>%
  st_transform(st_crs(grid_OBBA2))%>%
  st_transform(st_crs(AEA_proj))

Great_Lakes <- allBCR %>% subset(BCRNAME == "GREAT LAKES") %>%
  st_transform(st_crs(water)) %>%
  mutate(WATERBODY_ = "Lake") %>%
  select(WATERBODY_, geometry)%>%
  st_transform(st_crs(AEA_proj))

water <- rbind(water, Great_Lakes)

# ==========================
# CRS transformations
# ==========================

all_surveys <- all_surveys %>% st_transform(st_crs(AEA_proj))
grid_OBBA2 <- grid_OBBA2 %>% st_centroid() %>% st_transform(st_crs(AEA_proj))
grid_OBBA3 <- grid_OBBA3 %>% st_centroid() %>% st_transform(st_crs(AEA_proj))
study_boundary <- study_boundary %>% st_transform(st_crs(AEA_proj))

# ==========================
# BCR membership
# ==========================

sf::sf_use_s2(FALSE)

ONBCR <- subset(allBCR, PROVINCE_S == "ONTARIO")

grid_OBBA2 <- st_join(
  grid_OBBA2,
  ONBCR[, "BCR"],
  join = st_within,
  left = TRUE
)

grid_OBBA3 <- st_join(
  grid_OBBA3,
  ONBCR[, "BCR"],
  join = st_within,
  left = TRUE
)
  
# ==========================
# Filter surveys 
# ==========================

all_surveys <- all_surveys %>%
  subset(round(Survey_Duration_Minutes) == 5 &
           Hours_Since_Sunrise >= -1 &
           Hours_Since_Sunrise <= 6 &
           DayOfYear >= 135 &
           DayOfYear <= 196) %>%
  arrange(obs_idx) %>%
  relocate(obs_idx)

# keep only points that fall inside study boundary
all_surveys <- all_surveys[st_within(all_surveys, study_boundary, sparse = FALSE), ]

# retain a single survey per repeated location
all_surveys <- all_surveys %>%
  mutate(x = round(st_coordinates(.)[,1], digits = 5),
         y = round(st_coordinates(.)[,2], digits = 5)) %>%
  group_by(x, y) %>%
  slice_sample(n = 1) %>%
  ungroup() %>%
  select(-x, -y)

# subset counts to remaining surveys
counts <- counts %>% as_tibble() %>% mutate(obs_idx = seq_len(n()))
counts <- counts[all_surveys$obs_idx,]

# Add some additional columns
all_surveys <- all_surveys %>%
  mutate(Atlas3 = ifelse(Atlas == "OBBA2",0,1),
         square_atlas = as.numeric(factor(paste0(square_id,"-",Atlas))),
         square_year = as.numeric(as.factor(paste0(square_id,year(Date_Time)))),
         ARU = ifelse(Survey_Type == "ARU",1,0))

# ==========================
# Water overlap
# ==========================

overlaps <- st_intersects(grid_OBBA2, water)
grid_OBBA2$on_water <- lengths(overlaps) > 0
grid_OBBA3$on_water <- lengths(overlaps) > 0

# ==========================
# Covariate standardization: only for temp and precip
# ==========================

 covars <- c("insect_broadleaf","insect_needleleaf","road","water_river","prec","tmax",
             "urban_2","urban_3","lc_1","lc_4","lc_5","lc_8","lc_9","lc_10","lc_11","lc_12",
             "lc_14","lc_17")


# compute mean and sd for each covariate
covar_stats <- all_surveys %>%
  as.data.frame() %>%
  dplyr::select(all_of(covars)) %>%
  summarise(across(
    everything(),
    list(mean = \(x) mean(x, na.rm = TRUE),
         sd   = \(x) sd(x, na.rm = TRUE)
    )
  ))

covars <- c("prec","tmax")

# standardize
for (col in covars) {
  mu <- covar_stats[[paste0(col, "_mean")]]
  sigma <- covar_stats[[paste0(col, "_sd")]]
  all_surveys[[col]] <- (all_surveys[[col]] - mu) / sigma
  grid_OBBA2[[col]] <- (grid_OBBA2[[col]] - mu) / sigma
  grid_OBBA3[[col]] <- (grid_OBBA3[[col]] - mu) / sigma
}

# ==========================
# Covariate standardization: only for temp and precip
# ==========================


# ==========================
# Ensure grid has no NA values
# ==========================

grid_OBBA2 <- na.omit(grid_OBBA2)
grid_OBBA3 <- na.omit(grid_OBBA3)

grid_OBBA2$pixel_id <- grid_OBBA3$pixel_id <- 1:nrow(grid_OBBA2)

# ==========================
# Summarize data availability for each species
# ==========================

species_list <- analysis_data_covariates$all_species %>%
  mutate(species_id = as.character(species_id)) %>% 
  distinct(species_id, .keep_all = TRUE)

species_square_summary <- counts %>%
  pivot_longer(
    cols = -obs_idx,            
    names_to = "species_id",
    values_to = "count"
  ) %>%
  left_join(all_surveys, by = "obs_idx") %>%
  group_by(species_id, square_id, Atlas,BCR) %>%
  
  # Number of surveys that detected the species within each square
  summarise(count = sum(count>0), .groups = "drop") %>%
  group_by(species_id, Atlas,BCR) %>%
  
  # Number of squares in which the species was detected
  summarise(n_det = sum(count),
            n_sq = sum(count > 0), .groups = "drop") %>%
  left_join(species_list[,c("species_id","english_name")], by = "species_id") %>%
  dplyr::select(english_name, species_id, BCR,Atlas, n_sq,n_det) %>%
  distinct()

species_total_squares <- species_square_summary %>%
  group_by(species_id, english_name, Atlas) %>%
  summarise(total_detections = sum(n_det),
            total_squares = sum(n_sq), 
            .groups = "drop") %>%
  pivot_wider(names_from = Atlas, values_from = c(total_detections,total_squares), values_fill = 0)

species_to_model <- species_total_squares %>% 
  arrange(desc(total_squares_OBBA3))

# ==========================
# Save as rds file
# ==========================

saveRDS(
  
  list(
    study_boundary = study_boundary,
    all_surveys = all_surveys,
    counts = counts,
    grid_OBBA2 = grid_OBBA2,
    grid_OBBA3 = grid_OBBA3,
    species_square_summary = species_square_summary,
    species_to_model = species_to_model,
    species_square_summary = species_square_summary,
    species_total_squares = species_total_squares,
    water = water), 
  
  file = "Analysis/Ontario_BBA/Data_Cleaned/data_ready_for_analysis.rds")
