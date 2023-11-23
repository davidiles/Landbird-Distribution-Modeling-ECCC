# ************************************************
# BAYESIAN ANALYSIS / SPECIES DISTRIBUTION MODELS FOR SASKATCHEWAN BREEDING BIRD ATLAS
# 
# 1) A 'processed data package' for analysis is prepared by previous script
# ************************************************

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

`%!in%` <- Negate(`%in%`)

# ------------------------------------------------
# Load packages
# ------------------------------------------------

require(tidyverse)
require(sf)
require(INLA)     # INLA_22.05.07  
require(inlabru)  # inlabru_2.7.0

# ------------------------------------------------
# CONDUCT ANALYSIS
# ------------------------------------------------

analysis_data <- readRDS("../Data_Cleaned/analysis_data_package.rds")
attach(analysis_data)

for (sp_code in species_to_model$Species_Code_BSC){
  
  print(sp_code)
  
  model_file <- paste0("../Output/Fitted_Models/Model_",sp_code,"_squareday.rds")
  
  # Skip this species if already run
  if (file.exists(model_file)) next
  
  # Prepare data for this species
  sp_dat <- all_surveys %>% 
    mutate(count = full_count_matrix[,sp_code]) %>%
    
    # Only select point counts
    subset(Survey_Type %in% c("Point_Count","ARU_SPT","ARU_SPM"))
    
  
  # ------------------------------------------------
  # Extract ebird range for this species (if it exists); prepared by 4_Prep_Analysis.R
  # ------------------------------------------------
  
  if (sp_code %in% names(species_ranges)){
    
    range <- species_ranges[[sp_code]] %>% st_transform(st_crs(sp_dat))
    
    # Identify distance of each survey to the edge of species range (in km)
    sp_dat$distance_from_range <- ((st_distance(sp_dat, range)) %>% as.numeric())/1000
  } else{
    sp_dat$distance_from_range <- 0
  }
  
  # ------------------------------------------------
  # # Standardize 'time since sunrise' covariate (TSS)
  # ------------------------------------------------
  
  sp_dat$Hours_Since_Sunrise <- as.data.frame(sp_dat$Hours_Since_Sunrise)[,1] %>% scale()
  
  # ------------------------------------------------
  # Prepare iid random effects
  # ------------------------------------------------
  
  # Define 'square ID' covariate which will be treated as a random effect
  sp_dat$sq_idx <- factor(sp_dat$sq_id) %>% as.numeric()
  
  # Define 'square-day' covariate which will be treated as a random effect
  sp_dat$square_day <- paste0(sp_dat$sq_id,"-",yday(sp_dat$Date_Time)) %>% factor() %>% as.numeric()
  
  # Define unique survey location / year covariate
  coords <- st_coordinates(sp_dat) %>% round() %>% as.data.frame()
  X <- coords$X
  Y <- coords$Y
  Year <- year(sp_dat$Date_Time)
  sp_dat$loc_year <- paste0(Year,"-",X,"-",Y) %>% factor() %>% as.numeric()
  
  # --------------------------------
  # Generate QPAD offsets for each survey (assumes unlimited distance point counts)
  # --------------------------------
  
  species_offsets <- subset(species_to_model, Species_Code_BSC == sp_code)
  
  if (species_offsets$offset_exists == FALSE) sp_dat$log_QPAD_offset <- 0
  
  if (species_offsets$offset_exists == TRUE){
    A_metres <- pi*species_offsets$EDR^2
    p <- 1-exp(-sp_dat$Survey_Duration_Minutes*species_offsets$cue_rate)
    sp_dat$log_QPAD_offset <- log(A_metres * p)
    
  }
  
  # ------------------------------------------------
  # Separate Data types
  # ------------------------------------------------
  
  # Prepare point count data for this species
  PC_sf <- subset(sp_dat,Survey_Type == "Point_Count") 
  
  # Prepare ARU data for this species
  ARU_sf <- subset(sp_dat,Survey_Type %in% c("ARU_SPT","ARU_SPM")) 
  
  # Stationary counts
  #SC_sf <- subset(sp_dat, Survey_Type == "Breeding Bird Atlas") %>% mutate(presence = as.numeric(count>0))
  
  # Linear transects
  #LT_sf <- subset(sp_dat, Survey_Type == "Linear transect") %>% mutate(presence = as.numeric(count>0))
  
  # ------------------------------------------------
  # z-standardize effort covariates for checklists
  # ------------------------------------------------
  
  #SC_sf$SC_duration <- scale(SC_sf$Survey_Duration_Minutes)
  #LT_sf$LT_duration <- scale(LT_sf$Survey_Duration_Minutes)
  #LT_sf$LT_distance <- scale(LT_sf$Travel_Distance_Metres)
  
  # ------------------------------------------------
  # Convert to spatial objects
  # ------------------------------------------------
  
  PC_sp <- as(PC_sf,'Spatial')
  ARU_sp <- as(ARU_sf,'Spatial')
  #SC_sp <- as(SC_sf,'Spatial')
  #LT_sp <- as(LT_sf,'Spatial')
  
  # ------------------------------------------------
  # Create a spatial mesh, which is used to fit the residual spatial field
  # ------------------------------------------------
  
  # Note: mesh developed using tutorial at: https://rpubs.com/jafet089/886687
  max.edge = diff(range(st_coordinates(sp_dat)[,1]))/15
  bound.outer = diff(range(st_coordinates(sp_dat)[,1]))/3
  cutoff = max.edge/5
  bound.outer = diff(range(st_coordinates(sp_dat)[,1]))/3
  
  mesh_spatial <- inla.mesh.2d(loc = st_coordinates(sp_dat),
                               cutoff = max.edge/2, #max.edge/5,
                               max.edge = c(1,2)*max.edge,
                               offset=c(max.edge, bound.outer))
  mesh_locs <- mesh_spatial$loc[,c(1,2)] %>% as.data.frame()
  dim(mesh_locs)
  
  matern_coarse <- inla.spde2.pcmatern(mesh_spatial,
                                       prior.range = c(250000,0.1), # 10% chance range is smaller than 250000
                                       prior.sigma = c(2, 0.1), # 10% chance sd is larger than 2
                                       constr = TRUE # sum to 0 constraint
  )                
  
  # ------------------------------------------------
  # Define random effect prior
  # ------------------------------------------------
  
  kappa_prec <- list(prior = "pcprec", param = c(1,0.1))
  
  # ------------------------------------------------
  # Create mesh to model effect of time since sunrise (TSS)
  # ------------------------------------------------
  
  TSS_range <- range(sp_dat$Hours_Since_Sunrise)
  TSS_meshpoints <- seq(TSS_range[1]-0.1,TSS_range[2]+0.1,length.out = 11)
  TSS_mesh1D = inla.mesh.1d(TSS_meshpoints,boundary="free")
  TSS_spde = inla.spde2.pcmatern(TSS_mesh1D,
                                 prior.range = c(3,0.5), 
                                 prior.sigma = c(2,0.1)) # 10% chance sd is larger than 2
  # 
  # # ------------------------------------------------
  # # Create mesh to model effect of checklist duration
  # # ------------------------------------------------
  # 
  # SC_duration_meshpoints <- seq(min(SC_sf$SC_duration)-0.1,max(SC_sf$SC_duration)+0.1,length.out = 11)
  # SC_duration_mesh1D = inla.mesh.1d(SC_duration_meshpoints,boundary="free")
  # SC_duration_spde = inla.spde2.pcmatern(SC_duration_mesh1D,
  #                                        prior.range = c(1,0.5), 
  #                                        prior.sigma = c(2,0.1)) # 10% chance sd is larger than 2
  # 
  # LT_duration_meshpoints <- seq(min(LT_sf$LT_duration)-0.1,max(LT_sf$LT_duration)+0.1,length.out = 11)
  # LT_duration_mesh1D = inla.mesh.1d(LT_duration_meshpoints,boundary="free")
  # LT_duration_spde = inla.spde2.pcmatern(LT_duration_mesh1D,
  #                                        prior.range = c(1,0.5), 
  #                                        prior.sigma = c(2,0.1)) # 10% chance sd is larger than 2
  # 
  # LT_distance_meshpoints <- seq(min(LT_sf$LT_distance)-0.1,max(LT_sf$LT_distance)+0.1,length.out = 11)
  # LT_distance_mesh1D = inla.mesh.1d(LT_distance_meshpoints,boundary="free")
  # LT_distance_spde = inla.spde2.pcmatern(LT_distance_mesh1D,
  #                                        prior.range = c(1,0.5), 
  #                                        prior.sigma = c(2,0.1)) # 10% chance sd is larger than 2
  # 
  # ------------------------------------------------
  # Model formulas
  # ------------------------------------------------
  
  # kappa_loc_year + 
  
  model_formula_PC = as.formula(paste0('count ~
                  Intercept_PC +
                  log_QPAD_offset +
                  TSS +
                  kappa_squareday +
                  kappa_squareID +
                  spde_coarse +
                  range_effect * distance_from_range +
                   ',
                                       paste0("Beta1_",covariates_to_include,'*',covariates_to_include, collapse = " + ")))
  
  model_formula_ARU = as.formula(paste0('count ~
                  Intercept_ARU +
                  log_QPAD_offset +
                  TSS +
                  spde_coarse +
                  kappa_squareday +
                  kappa_squareID +
                  range_effect * distance_from_range +
                   ',
                                       paste0("Beta1_",covariates_to_include,'*',covariates_to_include, collapse = " + ")))
  
  # ------------------------------------------------
  # Fit model to point counts only
  # ------------------------------------------------
  
  # kappa_loc_year(loc_year, model = "iid", constr = TRUE, hyper = list(prec = kappa_prec))+
  bru_options_reset()
  
  inits <- c(-5,-5,rep(0,length(covariates_to_include))) %>% as.list()
  names(inits) <- c("Intercept_PC","Intercept_ARU",paste0("Beta1_",covariates_to_include))
  
  model_components = as.formula(paste0('~
  Intercept_PC(1)+
  Intercept_ARU(1)+
  range_effect(1,model="linear", mean.linear = -0.046, prec.linear = 10000)+
  TSS(main = Hours_Since_Sunrise,model = TSS_spde) +
  kappa_squareday(square_day, model = "iid", constr = TRUE, hyper = list(prec = kappa_prec))+
  kappa_squareID(sq_idx, model = "iid", constr = TRUE, hyper = list(prec = kappa_prec))+
  spde_coarse(main = coordinates, model = matern_coarse) +
  ',
                                       
                                       paste0("Beta1_",covariates_to_include,'(1,model="linear", mean.linear = 0, prec.linear = 100)', collapse = " + "))
  )
  
  start <- Sys.time()
  fit_PConly <- bru(components = model_components,
                    
                    like(family = "poisson",
                         formula = model_formula_PC,
                         data = PC_sp),
                    
                    like(family = "poisson",
                         formula = model_formula_ARU,
                         data = ARU_sp),
                    
                    options = list(
                      control.compute = list(waic = TRUE, cpo = TRUE),
                      bru_verbose = 4,
                      bru_initial = inits))
  end <- Sys.time() 
  
  runtime_PConly <- difftime( end,start, units="mins")
  print(paste0(sp_code," - ",runtime_PConly)) # ~10 min 
  
  species_fit <- list(sp_dat = sp_dat,
                      fit_PConly = fit_PConly)
  
  saveRDS(species_fit,model_file)
  
  gc(reset = TRUE) 
  
  # 
  # model_formula_SC = as.formula(paste0('presence ~ log(1/exp(-exp(
  # 
  #                 Intercept_SC +
  #                 TSS +
  #                 kappa_squareID +
  #                 kappa_squareday +
  #                 spde_coarse +
  #                 SC_duration +
  #                 range_effect * distance_from_range +
  #                                      ',
  #                                      paste0("Beta1_",covariates_to_include,'*',covariates_to_include, collapse = " + "),
  #                                      "))-1)"))
  # 
  # model_formula_LT = as.formula(paste0('presence ~ log(1/exp(-exp(
  # 
  #                 Intercept_LT +
  #                 TSS +
  #                 kappa_squareID +
  #                 kappa_squareday +
  #                 spde_coarse +
  #                 LT_distance +
  #                 LT_duration +
  #                 range_effect * distance_from_range +
  #                                      ',
  #                                      paste0("Beta1_",covariates_to_include,'*',covariates_to_include, collapse = " + "),
  #                                      "))-1)"))
} # close species loop
