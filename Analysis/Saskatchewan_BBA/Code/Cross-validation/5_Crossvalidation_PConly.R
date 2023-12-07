# ************************************************
# BAYESIAN ANALYSIS / SPECIES DISTRIBUTION MODELS FOR SASKATCHEWAN BREEDING BIRD ATLAS
# 
# 1) A 'processed data package' for analysis is prepared by previous script
# ************************************************

rm(list=ls())

# ------------------------------------------------
# Set working directory
# ------------------------------------------------

dirname <- "C:/Users/IlesD/OneDrive - EC-EC/Iles/Projects/Landbirds/Landbird-Distribution-Modeling-ECCC/Analysis/Saskatchewan_BBA/Code"
setwd(dirname)

`%!in%` <- Negate(`%in%`)

# ------------------------------------------------
# Load packages
# ------------------------------------------------

require(tidyverse)
require(sf)
require(INLA)     # INLA_22.05.07  
require(inlabru)  # inlabru_2.7.0
require(viridis)
require(terra)
require(stars)
require(ggtext)
library(pROC)

# Timeout INLA after 10 minutes (if it has not fit by then, it has likely stalled)
inla.setOption(inla.timeout = 60*10)

# ------------------------------------------------
# Function to rasterize a series of spatial predictions (needed for plotting)
# ------------------------------------------------

cut.fn <- function(df = NA, 
                   target_raster = NA, 
                   column_name = NA, 
                   lower_bound = NA, 
                   upper_bound = NA){
  
  max_val <- upper_bound
  max_val <- ifelse(is.na(max_val), 0, max_val)
  max_lev <- ifelse(max_val > 1.6, 4,ifelse(max_val > 0.8, 4, 3))
  
  cut_levs <- signif(max_val/(2^((max_lev-1):0)), 2)
  cut_levs <- unique(cut_levs)
  cut_levs <- ifelse(is.na(cut_levs), 0, cut_levs)
  
  if (lower_bound %in% cut_levs) cut_levs <- cut_levs[-which(cut_levs == lower_bound)]
  if (lower_bound > min(cut_levs)) cut_levs = cut_levs[-which(cut_levs < lower_bound)]
  
  max_lev <- length(cut_levs)
  
  cut_levs_labs <- c(paste0("0-",lower_bound),
                     paste(lower_bound, cut_levs[1], sep="-"),
                     paste(cut_levs[-max_lev], cut_levs[-1], sep="-"),
                     paste(cut_levs[max_lev], "+"))
  
  cut_levs <- c(-1, lower_bound, cut_levs, 1000) %>% unique()
  
  df$levs <- cut(as.data.frame(df)[,column_name], cut_levs, labels=cut_levs_labs)
  
  tmp = stars::st_rasterize(df %>% dplyr::select(levs, geometry))
  
  return(list(raster = tmp,cut_levs = cut_levs))
}

# ******************************************************************
# LOAD AND SUBSET BIRD DATA BASED ON SPECIFIC CRITERIA (DATE RANGES, ETC.)
# ******************************************************************

analysis_data <- readRDS("../Data_Cleaned/analysis_data_package.rds")
attach(analysis_data)

# Fix survey type labels
all_surveys$Survey_Type[all_surveys$Survey_Type == "IN_PERSON"] <- "Point_Count"
all_surveys$Survey_Type[all_surveys$Survey_Type %in% c("ARU_SM2","ARU_BAR_LT","ARU_SM_UNKN",
                                                       "ARU_SM4","ZOOM_H2N","ARU_IRIVER_E",
                                                       "ARU_MARANTZ","ARU_IRIVER","ARU_UNKNOWN")] <- "ARU_SPT"
AEA_proj <- "+proj=aea +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-106 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs "

SaskBoundary <- st_read("../../../Data/Spatial/National/BCR/BCR_Terrestrial_master.shp")  %>%
  subset(PROVINCE_S == "SASKATCHEWAN") %>%
  st_make_valid() %>%
  st_union() %>%
  st_transform(st_crs(AEA_proj))

SaskBoundary_buffer <- SaskBoundary %>% st_buffer(20000)

# Raster with target properties
target_raster <- rast("../../../Data/Spatial/National/AnnualMeanTemperature/wc2.1_30s_bio_1.tif") %>% 
  crop(st_transform(SaskBoundary,crs(.))) %>% 
  project(AEA_proj, res = 1000) %>%
  mask(vect(SaskBoundary))

SaskWater <- read_sf("../../../Data/Spatial/Saskatchewan/SaskWater/SaskWaterClip.shp") %>% 
  st_transform(st_crs(target_raster))
SaskWater$area <- st_area(SaskWater)
SaskWater$area <- as.numeric(SaskWater$area)
SaskWater <- SaskWater[SaskWater$area>2.580e+06 ,]

# Atlas Squares
SaskSquares <- st_read("../../../Data/Spatial/Saskatchewan/SaskSquares/SaskSquares.shp") %>%
  st_transform(st_crs(target_raster)) %>%
  dplyr::select(SQUARE_ID) %>%
  rename(sq_id = SQUARE_ID)

SaskGrid <- SaskGrid %>%
  st_intersection(SaskBoundary)

# ------------------------------------------
# Select Point Counts / ARUs to use
# ------------------------------------------

PC_to_use <- subset(all_surveys,
                    Survey_Type %in% c("Point_Count","ARU_SPT","ARU_SPM") &
                      
                      Survey_Duration_Minutes > 1 &
                      Survey_Duration_Minutes <= 10 &
                      
                      Hours_Since_Sunrise >= -2 &
                      Hours_Since_Sunrise <= 4 &
                      
                      yday(Date_Time) >= yday(ymd("2022-05-28")) &
                      yday(Date_Time) <= yday(ymd("2022-07-07")) &
                      
                      year(Date_Time) >= 2017 &
                      year(Date_Time) <= 2021
)

# ------------------------------------------
# Select STATIONARY COUNT data to use (assuming these are 'breeding bird atlas' checklists; vast majority do not have distance)
# ------------------------------------------

# Select stationary counts to use
SC_to_use <- subset(all_surveys,
                    
                    Survey_Type %in% c("Breeding Bird Atlas") &
                      
                      Hours_Since_Sunrise >= -2 &
                      Hours_Since_Sunrise <= 8 &
                      
                      Survey_Duration_Minutes >= 1 &
                      Survey_Duration_Minutes <= 120 &
                      
                      yday(Date_Time) >= yday(ymd("2022-05-28")) &
                      yday(Date_Time) <= yday(ymd("2022-07-07")) &
                      
                      year(Date_Time) >= 2017 &
                      year(Date_Time) <= 2021)

# ------------------------------------------
# Select LINEAR TRANSECT data to use
# ------------------------------------------

LT_to_use <- subset(all_surveys,
                    Survey_Type %in% c("Linear transect") &
                      
                      Hours_Since_Sunrise >= -2 &
                      Hours_Since_Sunrise <= 8 &
                      
                      Survey_Duration_Minutes > 10 &
                      Survey_Duration_Minutes <= 120 &
                      Travel_Distance_Metres <= 10000 &
                      
                      yday(Date_Time) >= yday(ymd("2022-05-28")) &
                      yday(Date_Time) <= yday(ymd("2022-07-07")) &
                      
                      year(Date_Time) >= 2017 &
                      year(Date_Time) <= 2021)

# ------------------------------------------
# Subset
# ------------------------------------------

surveys_to_use <- c(PC_to_use$Obs_Index, SC_to_use$Obs_Index, LT_to_use$Obs_Index)
all_surveys <- subset(all_surveys, Obs_Index %in% surveys_to_use) %>% arrange(Obs_Index)
full_count_matrix <- full_count_matrix[all_surveys$Obs_Index,]

# ------------------------------------------
# Calculate number of squares in which each species was detected, using each survey type
# ------------------------------------------

species_summary <- data.frame()

for (sp_code in colnames(full_count_matrix)){
  
  species_detections <- all_surveys %>% 
    as.data.frame() %>%
    mutate(count = full_count_matrix[,sp_code]) %>%
    mutate(Survey_Type = replace(Survey_Type, Survey_Type %in% c("Point_Count","ARU_SPT","ARU_SPM"), "PC/ARU")) %>%
    group_by(sq_id, Survey_Type) %>%
    summarize(n = sum(count>0,na.rm = TRUE))%>%
    group_by(Survey_Type) %>%
    subset(n>0)%>%
    summarize(n_squares = length(unique(sq_id))) %>%
    pivot_wider(names_from = Survey_Type, values_from = n_squares) %>%
    mutate(sp_code = sp_code) %>%
    relocate(sp_code)
  
  species_summary <- bind_rows(species_summary,species_detections)
  
}
species_summary[is.na(species_summary)] <- 0

# Number of squares in which each species was detected, with each survey type
species_summary <- species_summary %>% arrange(desc(`PC/ARU`))
head(species_summary)

species_summary <- species_summary %>% 
  subset(`PC/ARU` >=50) %>%
  left_join(.,species_to_model, by = c("sp_code" = "Species_Code_BSC")) %>%
  subset(!is.na(species_scientific_name))

# ******************************************************************
# Loop through species and conduct cross-validation
# ******************************************************************

# Create 10 x 10 km spatial folds
set.seed(999)
n_folds <- 5
Crossval_Grid <- st_make_grid(
  st_buffer(SaskBoundary,50000),
  cellsize = units::set_units(10*10,km^2),
  what = "polygons",
  square = TRUE,
  flat_topped = FALSE)%>%
  st_as_sf() %>%
  na.omit() %>%
  mutate(id = sample(1:nrow(.),nrow(.))) %>%
  mutate(Crossval_Fold = cut(id,breaks = seq(0,max(id)+1,length.out = n_folds+1)) %>% as.factor() %>% as.numeric())%>%
  dplyr::select(Crossval_Fold,x)

all_surveys$Obs_Index <- 1:nrow(all_surveys)
all_surveys <- all_surveys %>% st_intersection(Crossval_Grid) %>% arrange(Obs_Index)

results <- data.frame()
results_path <- "../Output/Crossvalidation/Crossval_results_PConly.rds"

set.seed(999)
species_summary <- sample_n(species_summary,25)

for (xval_fold in 1:n_folds){
#for (sp_code in species_summary$sp_code){
  for (sp_code in c("SPGR","PAWA","KILL","BLBW")){
    
    if (file.exists(results_path)){
      
      results <- readRDS(results_path)
      if (nrow(results[which(results$sp_code == sp_code & results$Crossval_Fold == xval_fold),])>0) next
      
    }
    
    print(paste0(sp_code," - crossval fold ", xval_fold))
    
    # ----------------------------------------------------
    # Extract counts/data for this species
    # ----------------------------------------------------
    
    # Prepare data for this species
    sp_dat <- all_surveys %>% 
      mutate(count = full_count_matrix[,sp_code],
             presence = as.numeric(full_count_matrix[,sp_code]>0)) %>%
      
      subset(Survey_Type %in% c("Point_Count","ARU_SPT","ARU_SPM",
                                "Breeding Bird Atlas",
                                "Linear transect"))
    
    # ----------------------------------------------------
    # Extract ebird range for this species (if it exists); prepared by 4_Prep_Analysis.R
    # ----------------------------------------------------
    
    if (sp_code %in% names(species_ranges)){
      
      range <- species_ranges[[sp_code]] %>% st_transform(st_crs(sp_dat))
      
      # Identify distance of each survey to the edge of species range (in km)
      sp_dat$distance_from_range <- ((st_distance(sp_dat, range)) %>% as.numeric())/1000
    } else{
      sp_dat$distance_from_range <- 0
    }
    
    # ----------------------------------------------------
    # Generate QPAD offsets for each survey (assumes unlimited distance point counts)
    # ----------------------------------------------------
    
    species_offsets <- subset(species_to_model, Species_Code_BSC == sp_code)
    if (nrow(species_offsets)==0) next
    
    if (species_offsets$offset_exists == FALSE) sp_dat$log_QPAD_offset <- 0
    
    if (species_offsets$offset_exists == TRUE){
      A_metres <- pi*species_offsets$EDR^2
      p <- 1-exp(-sp_dat$Survey_Duration_Minutes*species_offsets$cue_rate)
      sp_dat$log_QPAD_offset <- log(A_metres * p)
    }
    
    # ----------------------------------------------------
    # Separate Data types
    # ----------------------------------------------------
    sp_dat$presence <- as.numeric(sp_dat$count>0)
    
    validation_data <- subset(sp_dat, Crossval_Fold == xval_fold) %>%
      subset(Survey_Type %in% c("Point_Count","ARU_SPT","ARU_SPM")) %>% 
      as('Spatial')
    
    training_data <- subset(sp_dat, Crossval_Fold != xval_fold)
    PC_dat <- subset(training_data,Survey_Type %in% c("Point_Count","ARU_SPT","ARU_SPM")) %>% as('Spatial')
    SC_dat <- subset(training_data,Survey_Type %in% c("Breeding Bird Atlas")) %>% as('Spatial')
    LT_dat <- subset(training_data,Survey_Type %in% c("Linear transect")) %>% as('Spatial')
    
    # ----------------------------------------------------
    # Create a spatial mesh, which is used to fit the residual spatial field
    # NOTE: MAKING A COARSER MESH (LARGER VALUES OF MAX_EDGE) WILL SPEED UP COMPUTATIONS AND REDUCE RAM REQUIREMENTS, 
    # BUT WITH POTENTIAL LOSS OF ACCURACY/PRECISION
    # DEFAULT SHOULD CREATE APPROX 5000
    # ----------------------------------------------------
    
    # make a two extension hulls and mesh for spatial model
    hull <- fm_extensions(
      SaskBoundary,
      convex = c(50000, 200000),
      concave = c(350000, 500000)
    )
    
    # Setting max.edge = c(30000,100000) requires longer to fit (7 min), but may increase model accuracy
    # Setting max.edge = c(50000,100000) allows the model to fit much more quickly (3.5 min), though possibly with reduced accuracy
    mesh_spatial <- fm_mesh_2d_inla(
      boundary = hull, 
      max.edge = c(50000, 100000), # km inside and outside
      cutoff = 5000, 
      crs = fm_crs(sp_dat)
    )
    
    mesh_locs <- mesh_spatial$loc[,c(1,2)] %>% as.data.frame()
    
    # dim(mesh_locs)
    # plot(mesh_spatial)
    
    # Controls the 'residual spatial field'.  This can be adjusted to create smoother surfaces.
    prior_range <- c(300000,0.1) # 10% chance range is smaller than 300000
    prior_sigma <- c(0.5,0.1)    # 10% chance sd is larger than 0.5
    matern_coarse <- inla.spde2.pcmatern(mesh_spatial,
                                         prior.range = prior_range, 
                                         prior.sigma = prior_sigma,
                                         constr = TRUE
    )
    
    # ----------------------------------------------------
    # Create mesh to model effect of time since sunrise (TSS)
    # Hours since sunrise is fit as a GAM-type effect.  
    # Recommend leaving these settings at these defaults
    # ----------------------------------------------------
    
    TSS_range <- range(sp_dat$Hours_Since_Sunrise)
    TSS_meshpoints <- seq(TSS_range[1]-0.1,TSS_range[2]+0.1,length.out = 11)
    TSS_mesh1D = inla.mesh.1d(TSS_meshpoints,boundary="free")
    TSS_spde = inla.spde2.pcmatern(TSS_mesh1D,
                                   prior.range = c(6,0.1), # 10% range is smaller than 6
                                   prior.sigma = c(1,0.1)) # 10% chance sd is larger than 1
    
    # ------------------------------------------------
    # Create mesh to model effect of checklist duration
    # ------------------------------------------------
    
    SC_duration_meshpoints <- seq(min(SC_dat$Survey_Duration_Minutes)-0.1,max(SC_dat$Survey_Duration_Minutes)+0.1,length.out = 11)
    SC_duration_mesh1D = inla.mesh.1d(SC_duration_meshpoints,boundary="free")
    SC_duration_spde = inla.spde2.pcmatern(SC_duration_mesh1D,
                                           prior.range = c(60,0.1),
                                           prior.sigma = c(1,0.5))
    
    LT_duration_meshpoints <- seq(min(LT_dat$Survey_Duration_Minutes)-10,max(LT_dat$Survey_Duration_Minutes)+10,length.out = 11)
    LT_duration_mesh1D = inla.mesh.1d(LT_duration_meshpoints,boundary="free")
    LT_duration_spde = inla.spde2.pcmatern(LT_duration_mesh1D,
                                           prior.range = c(60,0.1),
                                           prior.sigma = c(1,0.5))
    
    LT_distance_meshpoints <- seq(min(LT_dat$Travel_Distance_Metres)-1000,max(LT_dat$Travel_Distance_Metres)+1000,length.out = 11)
    LT_distance_mesh1D = inla.mesh.1d(LT_distance_meshpoints,boundary="free")
    LT_distance_spde = inla.spde2.pcmatern(LT_distance_mesh1D,
                                           prior.range = c(5000,0.1),
                                           prior.sigma = c(1,0.5))
    
    # ***************************************************************
    # Model formulas
    # ***************************************************************
    
    # Names of covariates to include in model. Ideally, covariates should be uncorrelated with each other.
    # Each covariate will include a quadratic effect to allow for 'intermediate optimum' effects
    covariates_to_include <- paste0("PC",1:7) 
    
    # How much shrinkage should be applied to covariate effects?
    sd_linear <- 1 # Change to smaller value (e.g., 0.1), if you want to heavily shrink covariate effects and potentially create smoother surfaces
    prec_linear <-  c(1/sd_linear^2,1/(sd_linear/2)^2)
    
    #SC_duration(main = Survey_Duration_Minutes,model = SC_duration_spde) +
    #Intercept_SC(1)+
    model_components = as.formula(paste0('~
            Intercept_PC(1)+
            range_effect(1,model="linear", mean.linear = -0.046, prec.linear = 10000)+
            TSS(main = Hours_Since_Sunrise,model = TSS_spde) +
            spde_coarse(main = coordinates, model = matern_coarse) +',
                                         paste0("Beta1_",covariates_to_include,'(1,model="linear", mean.linear = 0, prec.linear = ', prec_linear[1],')', collapse = " + "),
                                         '+',
                                         paste0("Beta2_",covariates_to_include,'(1,model="linear", mean.linear = 0, prec.linear = ', prec_linear[2],')', collapse = " + ")))
    
    model_formula_PC = as.formula(paste0('count ~
                  Intercept_PC +
                  log_QPAD_offset +
                  TSS +
                  range_effect * distance_from_range +
                  spde_coarse +',
                                         paste0("Beta1_",covariates_to_include,'*',covariates_to_include, collapse = " + "),
                                         '+',
                                         paste0("Beta2_",covariates_to_include,'*',covariates_to_include, collapse = " + ")))
    
    model_formula_SC = as.formula(paste0('count ~
                  Intercept_SC +
                  log_QPAD_offset +
                  TSS +
                  SC_duration +
                  range_effect * distance_from_range +
                  spde_coarse +',
                                         paste0("Beta1_",covariates_to_include,'*',covariates_to_include, collapse = " + "),
                                         '+',
                                         paste0("Beta2_",covariates_to_include,'*',covariates_to_include, collapse = " + ")))
    
    model_formula_LT = as.formula(paste0('presence ~ log(1/exp(-exp(
                  Intercept_LT +
                  TSS +
                  spde_coarse +
                  LT_distance +
                  LT_duration +
                  range_effect * distance_from_range +',
                                         paste0("Beta1_",covariates_to_include,'*',covariates_to_include, collapse = " + "),
                                         '+',
                                         paste0("Beta2_",covariates_to_include,'*',covariates_to_include, collapse = " + "),"))-1)"))
    
    # ----------------------------------------------------------------------------
    # Fit model to point counts and ARUs with negative binomial error
    # ----------------------------------------------------------------------------
    
    start <- Sys.time()
    fit_INLA <- NULL
    while(is.null(fit_INLA)){
      fit_INLA <- bru(components = model_components,
                      
                      like(family = "nbinomial",
                           formula = model_formula_PC,
                           data = PC_dat),
                      # 
                      # like(family = "nbinomial",
                      #      formula = model_formula_SC,
                      #      data = SC_dat),
                      
                      # like(family = "binomial",
                      #      formula = model_formula_LT,
                      #      data = LT_dat),
                      
                      options = list(
                        control.compute = list(waic = FALSE, cpo = FALSE),
                        bru_verbose = 4))
      if ("try-error" %in% class(fit_INLA)) fit_INLA <- NULL
    }
    
    end <- Sys.time()
    runtime_INLA <- difftime( end,start, units="mins") %>% round(2)
    print(paste0(sp_code," - ",runtime_INLA," min to fit model")) # 2 min
    
    # ****************************************************************************
    # ****************************************************************************
    # Validation
    # ****************************************************************************
    # ****************************************************************************
    
    pred_formula_PC = as.formula(paste0(' ~
                  Intercept_PC +
                  log_QPAD_offset +
                  range_effect * distance_from_range +
                  spde_coarse +',
                                        paste0("Beta1_",covariates_to_include,'*',covariates_to_include, collapse = " + "),
                                        '+',
                                        paste0("Beta2_",covariates_to_include,'*',covariates_to_include, collapse = " + ")))
    
    
    pred <- NULL
    pred <- generate(fit_INLA,
                     validation_data,
                     formula =  pred_formula_PC,
                     n.samples = 1000)
    
    pred <- exp(pred)
    pred <- apply(pred,1,function(x) median(x))
    
    # Parameters defining the negative binomial distribution
    size <- fit_INLA$summary.hyperpar$'0.5quant'[1]
    prob <- size/(size + pred)
    
    # Probability of zero
    prob_zero <- dnbinom(0,mu=pred ,size=size)
    prob_nonzero <- 1-prob_zero
    
    # Cross-validation statistics
    AUC <- auc(response = validation_data$presence,predictor = prob_nonzero)
    lppd <- sum(dnbinom(validation_data$count,mu=pred,size=size,log = TRUE))
    
    # *********************************************************************
    # Save results
    # *********************************************************************
    
    run_results <- data.frame(sp_code = sp_code,
                              n_detections_val = sum(validation_data$count>0),
                              mean_count_val_obs = mean(validation_data$count),
                              mean_count_val_pred = mean(pred),
                              Crossval_Fold = xval_fold,
                              lppd_PConly = lppd,
                              AUC_PConly = AUC)
    
    if (file.exists(results_path)) results <- readRDS(results_path)
    
    results <- rbind(results, run_results)
    print(results)
    saveRDS(results,results_path)
    
  } # close species loop
} # close xval loop


