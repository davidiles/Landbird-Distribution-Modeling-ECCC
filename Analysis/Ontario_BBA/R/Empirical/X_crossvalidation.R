# ===============================================================================
# BAYESIAN ANALYSIS / SPECIES DISTRIBUTION MODELS FOR ONTARIO BREEDING BIRD ATLAS
# ===============================================================================

library(tidyverse)
library(sf)
library(fmesher)
library(INLA)
library(blockCV)
library(scoringRules)

rm(list=ls())

# Load birdDistribution package
devtools::load_all("C:/Users/IlesD/OneDrive - EC-EC/Iles/Projects/Landbirds/birdDistribution")

# ===============================================================================
# Load processed data and covariates, created by previous scripts
# ===============================================================================

analysis_data_covariates <- readRDS(file = "data_clean/analysis_data_covariates.rds")
all_species <- analysis_data_covariates$all_species
species_ranges <- analysis_data_covariates$species_ranges

# Survey-level covariates
all_surveys <- analysis_data_covariates$all_surveys %>%
  mutate(Atlas = ifelse(year(Date_Time)>2020,"OBBA3","OBBA2"),
         obs_idx = seq_len(n()))

# Counts of each species at each survey
counts <- analysis_data_covariates$full_count_matrix  

# 1 km grid across Ontario with covariates
grid_sf <- analysis_data_covariates$grid_sf %>%
  # Only keep complete rows (no NAs in attribute columns) & valid geometry
  filter(complete.cases(st_drop_geometry(.)) &  !is.na(st_geometry(.)))

AEA_proj <- st_crs(all_surveys)

# Study area boundary
study_boundary <- st_read("data/spatial/BCR/BCR_Terrestrial_master.shp") %>%
  subset(PROVINCE_S == "ONTARIO") %>%
  st_make_valid() %>%
  dplyr::select(BCR, PROVINCE_S) %>%
  st_transform(AEA_proj)%>%
  st_union() %>%
  st_set_precision(1)

# BCRs in Ontario
BCR <- st_read("data/spatial/BCR/BCR_Terrestrial_master.shp") %>%
  subset(PROVINCE_S == "ONTARIO") %>%
  st_make_valid() %>%
  group_by(BCR,BCRNAME) %>%
  summarize() %>%
  st_transform(AEA_proj)

# ---- Determine atlas square in which each survey is located

# Load shapefile containing atlas squares (10km x 10km)
atlas_squares <- st_read("C:/Users/IlesD/OneDrive - EC-EC/Iles/Projects/Landbirds/Landbird-Distribution-Modeling-ECCC/Data/Spatial/National/AtlasSquares/NationalSquares_FINAL.shp")  %>%
  subset(prov == "ON") %>%
  st_transform(AEA_proj)

# Create square_id and square_atlas identifiers
all_surveys <- all_surveys %>% 
  st_intersection(atlas_squares %>% dplyr::select(square_id)) %>%
  mutate(square_id = as.numeric(as.factor(square_id)),
         square_atlas = as.numeric(as.factor(paste0(square_id,Atlas))),
         square_year = as.numeric(as.factor(paste0(square_id,year(Date_Time)))))%>%
  relocate(geometry,.after = last_col()) %>%
  arrange(obs_idx)


# ===============================================================================
# Select data to use in models
# ===============================================================================
all_surveys$Survey_Duration_Minutes <- round(all_surveys$Survey_Duration_Minutes)
all_surveys <- subset(all_surveys,
                      
                      Survey_Duration_Minutes == 5 &
                        
                        Hours_Since_Sunrise >= -1 &
                        Hours_Since_Sunrise <= 6  &
                        
                        yday(Date_Time) >= yday(ymd("2022-06-15")) &
                        yday(Date_Time) <= yday(ymd("2022-07-15")) &
                        Project_Name %in% c("OBBA2",
                                            "OBBA3",
                                            "Ontario Breeding Bird Atlas Digital Point Counts 2024",
                                            "Ontario Breeding Bird Atlas Digital Point Counts 2023",
                                            "Ontario Breeding Bird Atlas Digital Point Counts 2022",
                                            "Ontario Breeding Bird Atlas Digital Point Counts 2021")) %>%
  arrange(obs_idx)  %>%
  relocate(obs_idx)

counts <- counts %>% rename(obs_idx = survey_ID) %>% mutate(obs_idx = seq_len(n()))
counts <- counts[all_surveys$obs_idx,]

# Checks
sum(is.na(counts)) # should be 0
mean(counts$obs_idx == all_surveys$obs_idx) # should be 1
dim(all_surveys) # 119760 surveys
dim(counts) # 119760 surveys

# ===============================================================================
# Summarize data availability for each species
# ===============================================================================

# ---- Which projects are in the final dataset?
survey_summary <- all_surveys %>%
  as.data.frame() %>%
  group_by(Survey_Type,Project_Name) %>%
  summarize(n = n()) %>%
  pivot_wider(names_from = Survey_Type, values_from = n) %>%
  mutate(Total = rowSums(across(where(is.numeric)), na.rm = TRUE)) %>%
  relocate(Project_Name,Total,Point_Count,ARU) %>%
  arrange(desc(Total))

survey_summary

# ---- Summarize number and proportion of atlas squares that detected each species
species_list <- analysis_data_covariates$all_species %>%
  mutate(species_id = as.character(species_id))

species_square_summary <- counts %>%
  pivot_longer(
    cols = -obs_idx,            
    names_to = "species_id",
    values_to = "count"
  ) %>%
  left_join(all_surveys, by = "obs_idx") %>%
  group_by(species_id, square_id, Atlas) %>%
  summarise(count = sum(count), .groups = "drop") %>%
  group_by(species_id, Atlas) %>%
  summarise(n_sq = sum(count > 0), .groups = "drop") %>%
  left_join(species_list, by = "species_id") %>%
  select(english_name, BSC_spcd, species_id, Atlas, n_sq)

species_square_summary

# Model species that were detected in at least 100 squares in both atlases
species_to_model <- species_square_summary %>%
  group_by(english_name, BSC_spcd, species_id, Atlas) %>%
  summarise(n_sq = sum(n_sq), .groups = "drop") %>%   # collapse duplicates
  pivot_wider(
    names_from = Atlas,
    values_from = n_sq,
    values_fill = 0
  ) %>%
  mutate(Total = OBBA2 + OBBA3) %>%
  filter(OBBA2 > 100, OBBA3 > 100) %>%
  arrange(desc(Total))

species_to_model

# ===============================================================================
# Define covariates
# ===============================================================================

# Data frame of covariates with priors
sd_linear <- 0.5
cov_df <- data.frame(covariate = paste0("PC", 1:4), 
                     beta = 1,
                     model = "linear",
                     mean = 0, prec = 1 / sd_linear^2) 
cov_df <- cov_df %>% bind_rows(cov_df %>% mutate(beta = 2, prec = 1 / (sd_linear / 2)^2))

# For incorporation into prediction formulas below
covariates <- cov_df %>% mutate(formula = paste0("Beta", beta, "_", covariate, "*", covariate, "^", beta))
covariate_terms <- paste(covariates$formula, collapse = " + ")

# ===============================================================================
# Define spatial blocks based on survey locations
# ===============================================================================

set.seed(123)
cv <- cv_spatial(x = all_surveys, # sf or SpatialPoints of sample data (e.g. species data)
                 size = 50, # size of the blocks in km
                 k = 5,      # number of folds
                 hexagon = TRUE, 
                 selection = "random", 
                 iteration = 100,
                 plot = TRUE)

all_surveys$cv_fold <- factor(cv$folds_ids)
cb_palette <- c(
  "#4477AA", # blue
  "#CC6677", # red/pink
  "#DDCC77", # yellow
  "#117733", # green
  "#88CCEE"  # cyan
)

ggplot(all_surveys)+ 
  geom_sf(aes(col = cv_fold))+ 
  scale_color_manual(values = cb_palette)+
  facet_grid(Atlas~.)+
  theme_bw()

# ===============================================================================
# Loop through species and CV folds, and conduct cross-validation analysis
# ===============================================================================

CV_folds <- unique(all_surveys$cv_fold) %>% sort()

sp_list <- c("Yellow Warbler","Winter Wren","White-throated Sparrow","American Robin","Lincoln's Sparrow","Savannah Sparrow","Fox Sparrow","Yellow-rumped Warbler",
             "Red-eyed Vireo","Swainson's Thrush","Nashville Warbler","Magnolia Warbler","Ruby-crowned Kinglet","Ovenbird",
             "American Goldfinch","Eastern Wood-Pewee","Chestnut-sided Warbler","Black-capped Chickadee","Song Sparrow","Blue Jay")

for (sp_english in sp_list){
  
  sp_code <- subset(species_to_model, english_name == sp_english)$species_id[1] %>% as.character()
  sp_range_poly <- analysis_data_covariates$species_ranges[[sp_english]] %>% sf::st_transform(sf::st_crs(all_surveys))
  
  # Prepare species survey data
  sp_dat <- all_surveys %>% 
    mutate(count = counts[[sp_code]],
           distance_from_range = as.numeric(sf::st_distance(.,sp_range_poly)),
           CRPS_1 = NA,
           CRPS_2 = NA,
           CRPS_3 = NA)
  
  for (fold in CV_folds[1:5]){
    
    print(paste0(sp_english, " - fold ",fold))
    
    # ---------------------------------------------
    # ---- Separate into training and validation data
    # ---------------------------------------------
    
    sp_dat_val <- subset(sp_dat, Atlas == "OBBA3" & cv_fold == fold)
    sp_dat_train <- subset(sp_dat, !(Atlas == "OBBA3" & cv_fold == fold))
    
    # ---------------------------------------------
    # ---- MODEL 1: Analysis of single atlas cycle (OBBA3)
    # ---------------------------------------------
    
    # Fit model
    mod1 <- fit_inla1(sp_dat = subset(sp_dat_train, Atlas == "OBBA3"),  
                      study_boundary, 
                      covariates = cov_df,
                      error_type = "poisson",
                      prior_range_abund = c(500,0.5),  # 50% chance range is smaller than 500 km
                      prior_sigma_abund = c(1,0.1))  # 10% chance SD is larger than 0.5
    
    # Predictions to withheld data
    pred_formula <- paste(
      c("~ Intercept",
        "spde_abund",
        "HSS",
        "range_effect * distance_from_range",
        covariate_terms),
      collapse = " + "
    ) %>% as.formula()
    
    preds1 <- inlabru::generate(mod1,
                                sp_dat_val,
                                formula = pred_formula,
                                n.samples = 1000,
                                seed = 123)
    
    # ---------------------------------------------
    # ---- MODEL 2: Analysis of both atlas cycles (OBBA2 and OBBA3) with no change model
    # ---------------------------------------------
    
    # Fit model
    mod2 <- fit_inla1(sp_dat = sp_dat_train,  
                      study_boundary, 
                      covariates = cov_df,
                      error_type = "poisson",
                      prior_range_abund = c(500,0.5),  # 50% chance range is smaller than 500 km
                      prior_sigma_abund = c(1,0.1))  # 10% chance SD is larger than 0.5
    
    # Predictions to withheld data
    pred_formula <- paste(
      c("~ Intercept",
        "spde_abund",
        "HSS",
        "range_effect * distance_from_range",
        covariate_terms),
      collapse = " + "
    ) %>% as.formula()
    
    preds2 <- inlabru::generate(mod2,
                                sp_dat_val,
                                formula = pred_formula,
                                n.samples = 1000,
                                seed = 123)
    
    # ---------------------------------------------
    # ---- MODEL 3: Analysis of both atlas cycles (OBBA2 and OBBA3) with change model
    # ---------------------------------------------
    
    # Fit model
    mod3 <- fit_inla2(sp_dat = sp_dat_train,  
                      study_boundary, 
                      covariates = cov_df,
                      error_type = "poisson",
                      prior_range_abund = c(500,0.5),  # 50% chance range is smaller than 500 km
                      prior_sigma_abund = c(1,0.1),  # 10% chance SD is larger than 1
                      prior_range_change = c(500,0.5), # 50% chance range is smaller than 500 km
                      prior_sigma_change = c(1,0.1)) # 10% chance SD is larger than 1
    
    # Predictions to withheld data
    pred_formula <- paste(
      c("~ Intercept_OBBA2",
        "effect_OBBA3",
        "spde_abund",
        "spde_change",
        "HSS",
        "range_effect * distance_from_range",
        covariate_terms),
      collapse = " + "
    ) %>% 
      as.formula()
    
    preds3 <- inlabru::generate(mod3,
                                sp_dat_val,
                                formula = pred_formula,
                                n.samples = 1000,
                                seed = 123)
    
    # ---------------------------------------------
    # ---- Calculate crossvalidation scores
    # ---------------------------------------------
    
    obs <- sp_dat_val$count  
    
    preds1 <- exp(preds1)
    preds2 <- exp(preds2)
    preds3 <- exp(preds3)
    
    # preds1[i, ] is the vector of 1000 posterior predictive samples for datapoint i
    crps1 <- sapply(seq_len(nrow(preds1)), function(i) {
      scoringRules::crps_sample(y = obs[i], dat = preds1[i, ])
    })
    
    crps2 <- sapply(seq_len(nrow(preds2)), function(i) {
      scoringRules::crps_sample(y = obs[i], dat = preds2[i, ])
    })
    
    crps3 <- sapply(seq_len(nrow(preds3)), function(i) {
      scoringRules::crps_sample(y = obs[i], dat = preds3[i, ])
    })
    
    # get the row indices of the validation points in sp_dat
    val_idx <- which(sp_dat$Atlas == "OBBA3" & sp_dat$cv_fold == fold)
    
    # fill in CRPS for those rows
    sp_dat$CRPS_1[val_idx] <- crps1
    sp_dat$CRPS_2[val_idx] <- crps2
    sp_dat$CRPS_3[val_idx] <- crps3
    
    
  } # close fold loop
  
  sp_dat_OBBA3 <- subset(sp_dat, Atlas == "OBBA3") %>%
    as.data.frame() %>%
    mutate(sp_english = sp_english, sp_code = sp_code) %>%
    dplyr::select(sp_english, sp_code, obs_idx, cv_fold, CRPS_1, CRPS_2, CRPS_3)
  
  saveRDS(sp_dat_OBBA3, file = paste0("model_output/spatial_xval/",sp_english,"_xval.rds"))
  
} # close species loop
