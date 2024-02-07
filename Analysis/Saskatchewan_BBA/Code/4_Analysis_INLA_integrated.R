# ************************************************
# BAYESIAN ANALYSIS / SPECIES DISTRIBUTION MODELS FOR SASKATCHEWAN BREEDING BIRD ATLAS
# 
# 1) A 'processed data package' for analysis is prepared by previous script
# ************************************************

rm(list=ls())

# ------------------------------------------------
# Set working directory
# ------------------------------------------------
# 
# stub <- function() {}
# thisPath <- function() {
#   cmdArgs <- commandArgs(trailingOnly = FALSE)
#   if (length(grep("^-f$", cmdArgs)) > 0) {
#     # R console option
#     normalizePath(dirname(cmdArgs[grep("^-f", cmdArgs) + 1]))[1]
#   } else if (length(grep("^--file=", cmdArgs)) > 0) {
#     # Rscript/R console option
#     scriptPath <- normalizePath(dirname(sub("^--file=", "", cmdArgs[grep("^--file=", cmdArgs)])))[1]
#   } else if (Sys.getenv("RSTUDIO") == "1") {
#     # RStudio
#     dirname(rstudioapi::getSourceEditorContext()$path)
#   } else if (is.null(attr(stub, "srcref")) == FALSE) {
#     # 'source'd via R console
#     dirname(normalizePath(attr(attr(stub, "srcref"), "srcfile")$filename))
#   } else {
#     stop("Cannot find file path")
#   }
# }
# 
# dirname <- thisPath()
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
require(exactextractr)

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

analysis_data <- readRDS("../Data_Cleaned/analysis_data_package_FromDean.rds")
attach(analysis_data)
all_surveys$Hours_Since_Sunrise <- as.numeric(all_surveys$Hours_Since_Sunrise)

table(as.data.frame(all_surveys)[,c("Survey_Type","Project_Name")])

survey_summary <- all_surveys %>%
  subset(year(Date_Time) >= 2017 &
           year(Date_Time) <= 2021) %>%
  as.data.frame() %>%
  group_by(Survey_Type,Project_Name) %>%
  summarize(n = n()) %>%
  pivot_wider(names_from = Survey_Type, values_from = n) %>%
  mutate(Total = rowSums(across(where(is.numeric)), na.rm = TRUE)) %>%
  relocate(Project_Name,Total,Point_Count,ARU_SPT) %>%
  arrange(desc(Total))

AEA_proj <- "+proj=aea +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-106 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs "

SaskBoundary <- st_read("../../../Data/Spatial/National/BCR/BCR_Terrestrial_master.shp")  %>%
  subset(PROVINCE_S == "SASKATCHEWAN") %>%
  st_make_valid() %>%
  st_union() %>%
  st_transform(st_crs(AEA_proj))
SaskGrid <- SaskGrid %>% st_intersection(SaskBoundary)

# Raster with target properties
target_raster <- rast("../../../Data/Spatial/National/AnnualMeanTemperature/wc2.1_30s_bio_1.tif") %>% 
  crop(st_transform(SaskBoundary,crs(.))) %>% 
  project(AEA_proj, res = 1000) %>%
  mask(vect(SaskBoundary))

SaskWater <- read_sf("../../../Data/Spatial/Saskatchewan/SaskWater/SaskWaterClip.shp") %>% 
  st_transform(st_crs(target_raster)) %>%
  mutate(area = as.numeric(st_area(.))) %>%
  subset(area > 2.580e+06)

# Atlas Squares
SaskSquares <- st_read("../../../Data/Spatial/Saskatchewan/SaskSquares/SaskSquares.shp") %>%
  st_transform(st_crs(target_raster)) %>%
  dplyr::select(SQUARE_ID) %>%
  rename(sq_id = SQUARE_ID)

# ------------------------------------------
# Select Point Counts / ARUs to use
# Decision to remove "SPM" aru transcriptions, as they are considered less reliable
# ------------------------------------------

PC_to_use <- subset(all_surveys,
                    Survey_Type %in% c("Point_Count","ARU_SPT") &
                      
                      Survey_Duration_Minutes > 1 &
                      Survey_Duration_Minutes <= 10 &
                      
                      Hours_Since_Sunrise >= -2 &
                      Hours_Since_Sunrise <= 4 &
                      
                      yday(Date_Time) >= yday(ymd("2022-05-28")) &
                      yday(Date_Time) <= yday(ymd("2022-07-07")) &
                      
                      year(Date_Time) >= 2017 &
                      year(Date_Time) <= 2021
)

dim(PC_to_use) # 26228

# ------------------------------------------
# Select STATIONARY COUNT data to use (assuming these are 'breeding bird atlas' checklists; vast majority do not have distance)
# ------------------------------------------

# Note that Dean Evans prepared a new data table describing overall checklist effort and species
# presence/absence records within each atlas square.  This is going to be used in place of observation-level
# checklists because it is not clear if stationary count checklists contain duplicate information from point counts.
# Using square-level summaries is anticipated to be less strong affected by duplicate information.

# Also note that covariate values associated with each square need to be re-calculated (the values in Dean's
# table are missing in many cases).  This is done below.

# # Select stationary counts to use
# SC_to_use <- subset(all_surveys,
#                     
#                     Survey_Type %in% c("Breeding Bird Atlas") &
#                       
#                       Hours_Since_Sunrise >= -2 &
#                       Hours_Since_Sunrise <= 8 &
#                       
#                       Survey_Duration_Minutes >= 10 &
#                       Survey_Duration_Minutes <= 120 &
#                       
#                       yday(Date_Time) >= yday(ymd("2022-05-28")) &
#                       yday(Date_Time) <= yday(ymd("2022-07-07")) &
#                       
#                       year(Date_Time) >= 2017 &
#                       year(Date_Time) <= 2021 &
#                       
#                       IncludesPointCounts == 0)
# 
# dim(SC_to_use) # 2026

# ------------------------------------------
# Calculate mean covariate valueswithin each Atlas square (covariates taken from SaskGrid)
# Note: these will be merged with new stationary checklist dataframe ("SaskSquareChecklist") within the species loop
# ------------------------------------------

# This vector should contain names of each covariate needed for analysis
covars <- paste0("PC",1:14)
cov_rast <- stars::st_rasterize(SaskGrid %>% dplyr::select(covars, geometry)) %>% rast()
SaskSquares <- cbind(SaskSquares,exact_extract(cov_rast,SaskSquares,"mean") %>% suppressWarnings())
colnames(SaskSquares)[2:(ncol(SaskSquares)-1)] <- covars
SaskSquares_centroids <- st_centroid(SaskSquares)

# ------------------------------------------
# Select LINEAR TRANSECT data to use
# ------------------------------------------

LT_to_use <- subset(all_surveys,
                    Survey_Type %in% c("Linear transect") &
                      
                      Hours_Since_Sunrise >= -2 &
                      Hours_Since_Sunrise <= 8 &
                      
                      Survey_Duration_Minutes > 10 &
                      Survey_Duration_Minutes <= 120 &
                      
                      Travel_Distance_Metres > 250 &
                      Travel_Distance_Metres <= 10000 &
                      
                      # Ensure speed is less than 3 metres per second
                      (Travel_Distance_Metres / (Survey_Duration_Minutes * 60)) <= 3 &
                      
                      yday(Date_Time) >= yday(ymd("2022-05-28")) &
                      yday(Date_Time) <= yday(ymd("2022-07-07")) &
                      
                      year(Date_Time) >= 2017 &
                      year(Date_Time) <= 2021 &
                      IncludesPointCounts == 0)

dim(LT_to_use) # 1748

# ------------------------------------------
# Subset
# ------------------------------------------

surveys_to_use <- c(PC_to_use$Obs_Index,  LT_to_use$Obs_Index) # SC_to_use$Obs_Index,
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
    mutate(Survey_Type = replace(Survey_Type, Survey_Type %in% c("Point_Count","ARU_SPT"), "PC/ARU")) %>%
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

# ******************************************************************
# Loop through species, fit models, generate maps
# ******************************************************************

species_summary <- species_summary %>% subset(`PC/ARU` >=20)

for (sp_code in rev(species_summary$sp_code)){
 
  # ----------------------------------------------------
  # Extract counts/data for this species
  # ----------------------------------------------------
  
  print(sp_code)
  
  # Prepare data for this species
  sp_dat <- all_surveys %>% 
    mutate(count = full_count_matrix[,sp_code],
           presence = as.numeric(full_count_matrix[,sp_code]>0)) %>%
    
    subset(Survey_Type %in% c("Point_Count","ARU_SPT",
                              #"Breeding Bird Atlas",
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
  # Prepare SC data (based on SaskSquareChecklist object created by Dean)
  # ----------------------------------------------------
  
  # Assign SC checklists (from Dean) to each crossval fold
  # Note that covariate information is missing/incorrect in the SaskSquareChecklist object, so need to replace it with SaskSquare_centroids
  SC_dat <- SaskSquareChecklist %>%
    as.data.frame() %>%
    dplyr::rename(species_code = sp_code) %>%
    subset(species_code == sp_code) %>%
    dplyr::select(sq_id,presence,total_checklists,total_effort,checklists_with_species,effort_to_species) %>% 
    left_join(SaskSquares_centroids,.) %>%
    subset(total_checklists > 0 & !is.na(total_effort)) 
  
  SC_dat$presence <- as.numeric(SC_dat$presence)
  
  # Extract ebird range for this species (if it exists); prepared by 4_Prep_Analysis.R
  if (sp_code %in% names(species_ranges)){
    
    range <- species_ranges[[sp_code]] %>% st_transform(st_crs(sp_dat))
    
    # Identify distance of each survey to the edge of species range (in km)
    SC_dat$distance_from_range <- ((st_distance(SC_dat, range)) %>% as.numeric())/1000
  } else{
    SC_dat$distance_from_range <- 0
  }
  SC_dat <- as(SC_dat, 'Spatial')
  
  # ----------------------------------------------------
  # Separate Data types
  # ----------------------------------------------------
  
  # Point count and ARU data treated as the same thing in this analysis
  PC_dat <- subset(sp_dat,Survey_Type %in% c("Point_Count","ARU_SPT")) %>% as('Spatial')

  # Stationary counts 
  #SC_dat <- subset(sp_dat, Survey_Type == "Breeding Bird Atlas") %>% mutate(presence = as.numeric(count>0)) %>% as('Spatial')
  
  # Linear transects
  LT_dat <- subset(sp_dat, Survey_Type == "Linear transect") %>% mutate(presence = as.numeric(count>0)) %>% as('Spatial')
  
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
  prior_range <- c(300000, 0.1) # 10% chance range is smaller than 300000
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
  # Create mesh to model effect of linear transect effort
  # ------------------------------------------------
  
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
  
  # ------------------------------------------------
  # Create mesh to model effect of stationary checklist effort
  # ------------------------------------------------
  
  SC_effort_meshpoints <- seq(min(SC_dat$total_effort)-1,max(SC_dat$total_effort)+1,length.out = 11)
  SC_effort_mesh1D = inla.mesh.1d(SC_effort_meshpoints,boundary="free")
  SC_effort_spde = inla.spde2.pcmatern(SC_effort_mesh1D,
                                       prior.range = c(30,0.1),
                                       prior.sigma = c(1,0.5))
  
  # ***************************************************************
  # Model formulas
  # ***************************************************************
  
  # Names of covariates to include in model. Ideally, covariates should be uncorrelated with each other.
  # Each covariate will include a quadratic effect to allow for 'intermediate optimum' effects
  covariates_to_include <- paste0("PC",1:7) 
  
  # How much shrinkage should be applied to covariate effects?
  sd_linear <- 0.1  # Change to smaller value (e.g., 0.1), if you want to heavily shrink covariate effects and potentially create smoother surfaces
  prec_linear <-  c(1/sd_linear^2,1/(sd_linear/2)^2)
  
  #Intercept_SC(1)+
  #  SC_effort(main = total_effort,model = SC_effort_spde) +
    
  model_components = as.formula(paste0('~
            Intercept_PC(1)+
            range_effect(1,model="linear", mean.linear = -0.046, prec.linear = 10000)+
            TSS(main = Hours_Since_Sunrise,model = TSS_spde) +
            
            Intercept_LT(1)+
            LT_duration(main = Survey_Duration_Minutes,model = LT_duration_spde) +
            LT_distance(main = Travel_Distance_Metres,model = LT_distance_spde) +
            
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
                                       paste0("Beta2_",covariates_to_include,'*',covariates_to_include,"^2", collapse = " + ")))
  
  model_formula_LT = as.formula(paste0('presence ~ log(1/exp(-exp(
                  Intercept_LT +
                  LT_distance +
                  LT_duration +
                  TSS +
                  spde_coarse +
                  range_effect * distance_from_range +',
                                       paste0("Beta1_",covariates_to_include,'*',covariates_to_include, collapse = " + "),
                                       '+',
                                       paste0("Beta2_",covariates_to_include,'*',covariates_to_include,"^2", collapse = " + "),"))-1)"))
  
  model_formula_SC = as.formula(paste0('presence ~ log(1/exp(-exp(
                  Intercept_SC +
                  SC_effort +
                  spde_coarse +
                  range_effect * distance_from_range +',
                                       paste0("Beta1_",covariates_to_include,'*',covariates_to_include, collapse = " + "),
                                       '+',
                                       paste0("Beta2_",covariates_to_include,'*',covariates_to_include,"^2", collapse = " + "),"))-1)"))
  
  # ----------------------------------------------------------------------------
  # Fit model to point counts and ARUs with negative binomial error
  # ----------------------------------------------------------------------------
  
  start <- Sys.time()
  fit_INLA <- NULL
  while(is.null(fit_INLA)){
    
    fit_model <- function(){
      tryCatch(expr = {bru(components = model_components,
                           
                           like(family = "nbinomial",
                                formula = model_formula_PC,
                                data = PC_dat),
                           
                           like(family = "binomial",
                                formula = model_formula_LT,
                                data = LT_dat),
                           
                            like(family = "binomial",
                                   formula = model_formula_SC,
                                   data = SC_dat),
                           
                           
                           options = list(bru_initial = list(lsig = 2,
                                                             Intercept_PC = -10,
                                                             Intercept_LT = 0,
                                                             Intercept_SC = 0),
                                          control.compute = list(waic = FALSE, cpo = FALSE),
                                          bru_verbose = 4))},
               error = function(e){NULL})
    }
    fit_INLA <- fit_model()
    
    if ("try-error" %in% class(fit_INLA)) fit_INLA <- NULL
  }
  
  end <- Sys.time()
  runtime_INLA <- difftime( end,start, units="mins") %>% round(2)
  print(paste0(sp_code," - ",runtime_INLA," min to fit model")) 
  
  # ****************************************************************************
  # ****************************************************************************
  # GENERATE MAPS
  # ****************************************************************************
  # ****************************************************************************
  
  # For every pixel on landscape, extract distance (in km) from eBird range limit
  SaskGrid_species <- SaskGrid %>%
    mutate(distance_from_range = (st_centroid(.) %>% 
                                    st_distance( . , range) %>% 
                                    as.numeric())/1000)
  
  # ----------------------------------------------------
  # QPAD offsets associated with a 5-minute unlimited distance survey
  # ----------------------------------------------------
  
  species_offsets <- subset(species_to_model, Species_Code_BSC == sp_code)
  log_offset_5min <- 0
  if (species_offsets$offset_exists == TRUE) log_offset_5min <- species_offsets$log_offset_5min
  
  # ----------------------------------------------------
  # Generate predictions on SaskGrid_species raster
  # Per 15 point counts; hence the log(15) in the prediction formula
  # ----------------------------------------------------
  
  pred_formula_PC = as.formula(paste0(' ~
                  Intercept_PC +
                  log_offset_5min +
                  log(15)+
                  range_effect * distance_from_range +
                  spde_coarse +',
                                      paste0("Beta1_",covariates_to_include,'*',covariates_to_include, collapse = " + "),
                                      '+',
                                      paste0("Beta2_",covariates_to_include,'*',covariates_to_include, collapse = " + ")))
  
  # Note that predictions are initially on log scale
  start2 <- Sys.time()
  pred <- NULL
  pred <- generate(fit_INLA,
                   as(SaskGrid_species,'Spatial'),
                   formula =  pred_formula_PC,
                   n.samples = 1000)
  
  pred <- exp(pred)
  
  # Median and upper/lower credible intervals (90% CRI)
  prediction_quantiles = apply(pred,1,function(x) quantile(x,c(0.05,0.5,0.95),na.rm = TRUE))
  SaskGrid_species$pred_q05 <- prediction_quantiles[1,]
  SaskGrid_species$pred_q50 <- prediction_quantiles[2,]
  SaskGrid_species$pred_q95 <- prediction_quantiles[3,]
  SaskGrid_species$pred_CI_width_90 <- prediction_quantiles[3,] - prediction_quantiles[1,]
  
  # Probability of observing species
  size <- fit_INLA$summary.hyperpar$'0.5quant'[1] # parameter of negative binomial
  prob_zero_PC <- dnbinom(0,mu=prediction_quantiles[2,],size=size)
  SaskGrid_species$pObs_5min <- 1-prob_zero_PC
  
  # Upper/lower 90% credible intervals on PObs
  SaskGrid_species$pObs_5min_q05 <- 1-dnbinom(0,mu=prediction_quantiles[1,],size=size)
  SaskGrid_species$pObs_5min_q95 <- 1-dnbinom(0,mu=prediction_quantiles[3,],size=size)
  
  
  # Convert to CRS of target raster
  SaskGrid_species <- SaskGrid_species %>% st_transform(st_crs(target_raster))
  
  end2 <- Sys.time() 
  runtime_pred <- difftime( end2,start2, units="mins") %>% round(2)
  print(paste0(sp_code," - ",runtime_pred," min to generate predictions")) # 4 min
  
  # ------------------------------------------------
  # Summarize SaskSquares where species was detected (for plotting)
  # ------------------------------------------------
  
  PC_detected <- sp_dat %>%
    subset(Survey_Type %in% c("Point_Count","ARU_SPT")) %>% 
    as.data.frame() %>%
    group_by(sq_id) %>%
    summarize(PC_detected = as.numeric(sum(count)>0),
              PC_mean_count = mean(count) %>% round(2))
  
  CL_detected <- sp_dat %>%
    subset(Survey_Type %in% c("Linear transect")) %>%
    as.data.frame() %>%
    group_by(sq_id) %>%
    summarize(CL_detected = as.numeric(sum(count)>0),
              CL_mean_count = mean(count))
  
  SaskSquares_species <- SaskSquares %>%
    relocate(geometry,.after = last_col()) %>%
    left_join(PC_detected) %>% 
    left_join(CL_detected)
  
  SaskSquares_centroids <- st_centroid(SaskSquares_species)
  
  # ------------------------------------------------
  # Label for figure
  # ------------------------------------------------
  
  species_name = Sask_spcd$CommonName[which(Sask_spcd$spcd == sp_code)]
  species_label = Sask_spcd$Label[which(Sask_spcd$spcd == sp_code)]
  
  # ------------------------------------------------
  # sf object for ebird range limit (optional - not needed for plotting, but potentially helpful)
  # ------------------------------------------------
  
  range <- NA
  if (sp_code %in% names(species_ranges)){
    range <- species_ranges[[sp_code]]  %>% 
      st_transform(st_crs(SaskBoundary)) %>%
      st_intersection(SaskBoundary)
  }
  
  # ****************************************************************************
  # FIGURE 1: RELATIVE ABUNDANCE
  # ****************************************************************************
  
  colscale_relabund <-c("#FEFEFE", "#FBF7E2", "#FCF8D0", "#EEF7C2", "#CEF2B0", "#94E5A0", "#51C987", "#18A065", "#008C59", "#007F53", "#006344")
  colpal_relabund <- colorRampPalette(colscale_relabund)

  lower_bound <- 0.15
  upper_bound <- quantile(SaskGrid_species$pred_q50,0.95,na.rm = TRUE) %>% signif(2)
  if (lower_bound >= (upper_bound/5)) lower_bound <- (upper_bound/5) %>% signif(2)
  
  sp_cut <- cut.fn(df = SaskGrid_species,
                   target_raster = target_raster,
                   column_name = "pred_q50",
                   lower_bound = lower_bound,
                   upper_bound = upper_bound)
  
  raster_q50 <- sp_cut$raster 
  
  plot_q50 <- ggplot() +
    
    geom_stars(data = raster_q50) +
    scale_fill_manual(name = "<span style='font-size:13pt'>Relative Abundance</span><br><span style='font-size:7pt'>Per 15 point counts</span>",
                      values = colpal_relabund(length(levels(raster_q50$levs))), drop=FALSE,na.translate=FALSE)+
    
    geom_sf(data = SaskWater,colour=NA,fill="#59F3F3",show.legend = F)+
    
    geom_sf(data = SaskBoundary,colour="black",fill=NA,lwd=0.3,show.legend = F) +
    
    #geom_sf(data = range,colour="darkred",fill=NA,show.legend = F) +
    
    # Surveyed squares
    geom_sf(data = subset(SaskSquares_centroids, !is.na(PC_detected) | !is.na(CL_detected)), col = "gray70", pch = 19, size = 0.1)+
    
    # Point count detections
    geom_sf(data = subset(SaskSquares_centroids, !is.na(PC_detected) & PC_detected > 0), col = "black", pch = 19, size = 0.2)+
    
    # Checklist detections
    geom_sf(data = subset(SaskSquares_species, !is.na(CL_detected) & CL_detected > 0), col = "black", size = 0.2, fill = "transparent")+
    
    
    coord_sf(clip = "off",xlim = range(as.data.frame(st_coordinates(SaskBoundary))$X))+
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())+
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    theme(legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(5,10,5,-20),
          legend.title.align=0.5,
          legend.title = element_markdown(lineheight=.9,hjust = "left"))+
    theme(legend.key = element_rect(fill = "transparent", colour = "transparent"))+
    annotate(geom="text",x=346000,y=2050000, label= paste0(species_label),lineheight = .85,hjust = 0,size=6,fontface =2) +
    annotate(geom="text",x=200000,y=970000, label= paste0("Prepared on ",Sys.Date()),size=1.5,lineheight = .75,hjust = 0,color="gray75")+
    guides(fill = guide_legend(order = 1), 
           size = guide_legend(order = 2))
  
  png(paste0("../Output/Prediction_Maps/Relative_Abundance/",sp_code,"_q50.png"), width=6.5, height=8, units="in", res=300, type="cairo")
  print(plot_q50)
  dev.off()
  
  # ****************************************************************************
  # FIGURE 2: UNCERTAINTY EXPRESSED AS WIDTH OF 90% CREDIBLE INTERVAL
  # ****************************************************************************
  
  colscale_uncertainty <- c("#FEFEFE", "#FFF4B3", "#F5D271", "#F2B647", "#EC8E00", "#CA302A")
  colpal_uncertainty <- colorRampPalette(colscale_uncertainty)
  
  lower_bound <- 0.15
  upper_bound <- quantile(SaskGrid_species$pred_CI_width_90,0.99,na.rm = TRUE) %>% signif(2)
  if (lower_bound >= (upper_bound/5)) lower_bound <- (upper_bound/5) %>% signif(2)
  
  raster_CI_width_90 <- cut.fn(df = SaskGrid_species,
                               target_raster = target_raster,
                               column_name = "pred_CI_width_90",
                               lower_bound = lower_bound,
                               upper_bound = upper_bound)$raster
  
  plot_CI_width_90 <- ggplot() +
    
    geom_stars(data = raster_CI_width_90) +
    scale_fill_manual(name = "<span style='font-size:13pt'>Uncertainty</span><br><span style='font-size:7pt'>Per 15 point counts</span><br><span style='font-size:7pt'>Width of 90% CI</span>",
                      values = colpal_uncertainty(length(levels(raster_CI_width_90$levs))), drop=FALSE,na.translate=FALSE)+
    
    geom_sf(data = SaskWater,colour=NA,fill="#59F3F3",show.legend = F)+
    
    geom_sf(data = SaskBoundary,colour="black",fill=NA,lwd=0.3,show.legend = F) +
    
    #geom_sf(data = range,colour="darkred",fill=NA,show.legend = F) +
    
    # Surveyed squares
    geom_sf(data = subset(SaskSquares_centroids, !is.na(PC_detected) | !is.na(CL_detected)), col = "gray70", pch = 19, size = 0.1)+
    
    # Point count detections
    geom_sf(data = subset(SaskSquares_centroids, !is.na(PC_detected) & PC_detected > 0), col = "black", pch = 19, size = 0.2)+
    
    # Checklist detections
    geom_sf(data = subset(SaskSquares_centroids, !is.na(CL_detected) & CL_detected > 0), col = "black", pch = 4, size = 0.2)+
    
    coord_sf(clip = "off",xlim = range(as.data.frame(st_coordinates(SaskBoundary))$X))+
    theme(panel.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    theme(legend.margin=margin(0,0,0,0),legend.box.margin=margin(5,10,5,-20),legend.title.align=0.5,
          legend.title = element_markdown(lineheight=.9,hjust = "left"))+
    annotate(geom="text",x=346000,y=1850000, label= paste0(species_label),lineheight = .85,hjust = 0,size=6,fontface =2) +
    annotate(geom="text",x=346000,y=1400000, label= paste0("Prepared on ",Sys.Date()),size=3,lineheight = .75,hjust = 0,color="#3b3b3b")
  print(plot_CI_width_90)
  
  png(paste0("../Output/Prediction_Maps/Uncertainty/",sp_code,"_CI_width_90.png"), width=6.5, height=8, units="in", res=300, type="cairo")
  print(plot_CI_width_90)
  dev.off()
  
  # ****************************************************************************
  # FIGURE 3: PROBABILITY OF OBSERVING SPECIES AFTER CONDUCTING 15 POINT COUNTS
  # ****************************************************************************
  
  colscale_pObs <- c("#FEFEFE",RColorBrewer::brewer.pal(5,"BuGn")[2:5])
  colpal_pObs <- colorRampPalette(colscale_pObs)
  
  cut_levs <- c(-0.1,0.05,0.125,0.25,0.5,1)
  cut_levs_labs <- c("0 to 0.05", 
                     "0.05 to 0.125", 
                     "0.125 to 0.25", 
                     "0.25 to 0.50", 
                     "0.50 to 1")
  
  SaskGrid_species$pObs_levs <- cut(as.data.frame(SaskGrid_species)[,"pObs_5min"], 
                                    cut_levs,labels=cut_levs_labs)
  
  raster_pObs = stars::st_rasterize(SaskGrid_species %>% dplyr::select(pObs_levs, geometry))
  
  plot_pObs <- ggplot() +
    
    geom_stars(data = raster_pObs) +
    scale_fill_manual(name = "<span style='font-size:13pt'>Prob. of Observation</span><br><span style='font-size:7pt'>Per 15 point counts </span><br><span style='font-size:7pt'>(Posterior Median)</span>",
                      values = colpal_pObs(length(levels(raster_pObs$pObs_levs))), 
                      drop=FALSE,na.translate=FALSE)+
    
    #geom_sf(data = SaskWater,colour=NA,fill="#59F3F3",show.legend = F)+
    
    geom_sf(data = SaskBoundary,colour="black",fill=NA,lwd=0.3,show.legend = F) +
    
    #geom_sf(data = range,colour="darkred",fill=NA,show.legend = F) +
    
    # Surveyed squares
    geom_sf(data = subset(SaskSquares_centroids, !is.na(PC_detected) | !is.na(CL_detected)), col = "gray70", pch = 19, size = 0.1)+
    
    # Point count detections
    geom_sf(data = subset(SaskSquares_centroids, !is.na(PC_detected) & PC_detected > 0), col = "black", pch = 19, size = 0.2)+
    
    # Checklist detections
    geom_sf(data = subset(SaskSquares_centroids, !is.na(CL_detected) & CL_detected > 0), col = "black", pch = 4, size = 0.2)+
    
    coord_sf(clip = "off",xlim = range(as.data.frame(st_coordinates(SaskBoundary))$X))+
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())+
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    theme(legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(5,10,5,-20),
          legend.title.align=0.5,
          legend.title = element_markdown(lineheight=.9,hjust = "left"))+
    annotate(geom="text",x=346000,y=1850000, label= paste0(species_label),lineheight = .85,hjust = 0,size=6,fontface =2) +
    annotate(geom="text",x=346000,y=1400000, label= paste0("Prepared on ",Sys.Date()),size=3,lineheight = .75,hjust = 0,color="#3b3b3b")
  
  png(paste0("../Output/Prediction_Maps/PObs/",sp_code,"_PObs.png"), width=6.5, height=8, units="in", res=300, type="cairo")
  print(plot_pObs)
  dev.off()

} # close species loop


