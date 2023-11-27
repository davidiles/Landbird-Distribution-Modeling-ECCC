# ************************************************
# BAYESIAN ANALYSIS / SPECIES DISTRIBUTION MODELS FOR SASKATCHEWAN BREEDING BIRD ATLAS
# 
# 1) A 'processed data package' for analysis is prepared by previous script
# ************************************************

rm(list=ls())

# ------------------------------------------------
# Set working directory
# ------------------------------------------------

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
dirname <- "C:/Users/IlesD/OneDrive - EC-EC/Iles/Backup/2023-11-09/SDM_ECCC/Analysis/Saskatchewan/Code"
setwd(dirname)

`%!in%` <- Negate(`%in%`)

# ------------------------------------------------
# Load packages
# ------------------------------------------------

require(tidyverse)
require(sf)
require(INLA)     # INLA_22.05.07  
require(inlabru)  # inlabru_2.7.0
library(dismo)
library(terra)
library(exactextractr)
library(magrittr)
library(factoextra)
library(pROC)
library(viridis)
library(ggpubr)

# Timeout INLA after 4 minutes (if it has not fit by then, it has likely stalled)
inla.setOption(inla.timeout = 60*10)

# ----------------------------------------------------------------
# Function to fit BRT
# ----------------------------------------------------------------

fit_brt <- function(model_data,
                    response_column = NA,
                    covariate_columns = NA){
  mod_brt <- NULL
  
  ntrees <- 50
  tcomplexity <- 5
  lrate <- 0.01
  m <- 0
  
  while(is.null(mod_brt)){
    
    m <- m + 1
    if(m < 11){
      ntrees <- 50
      lrate <- 0.01
    } else if(m < 21){
      lrate <- 0.001
    } else if(m < 31){
      ntrees <- 25
      lrate <- 0.001
    } else if(m < 41){
      ntrees <- 25
      lrate <- 0.0001
    } else if(m < 51){
      ntrees <- 10
      lrate <- 0.00001
    }
    else{
      break
    }
    
    ptm <- proc.time()
    if(inherits(try(
      mod_brt <- dismo::gbm.step(data = model_data,
                                 gbm.x = covariate_columns,
                                 gbm.y = response_column,
                                 offset = model_data$log_QPAD_offset,
                                 family = "poisson",
                                 tree.complexity = tcomplexity,
                                 learning.rate = lrate,
                                 n.trees = ntrees,
                                 n.folds = 5,
                                 max.trees = 10000)
    ), "try-error")){
      cat("Couldn't fit model", n, "in the iteration", m, "\n")
    }
    t <- proc.time() - ptm
  }
  if(is.null(mod_brt)){
    next
  }
  return(mod_brt)
}

# ------------------------------------------------
# LOAD DATA
# ------------------------------------------------

AEA_proj <- "+proj=aea +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-106 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs "

analysis_data <- readRDS("../Data_Cleaned/analysis_data_package.rds")
all_surveys <- analysis_data$all_surveys
full_count_matrix <- analysis_data$full_count_matrix
species_to_model <- analysis_data$species_to_model
species_ranges <- analysis_data$species_ranges
pca <- analysis_data$pca

all_surveys_PC <- all_surveys %>%
  subset(Survey_Type %in% c("Point_Count","ARU_SPT","ARU_SPM"))
full_count_matrix_PC <- analysis_data$full_count_matrix[all_surveys_PC$Obs_Index,]
dim(all_surveys_PC)

all_surveys <- all_surveys %>%
  subset(Survey_Type %in% c("Point_Count","ARU_SPT","ARU_SPM","Breeding Bird Atlas", "Linear transect"))
full_count_matrix <- analysis_data$full_count_matrix[all_surveys$Obs_Index,]
dim(all_surveys)

# ------------------------------------------------
# Divide Study Area into 50 x 50 km blocks for cross-validation
# ------------------------------------------------

SaskBoundary <- st_read("../../../Data/Spatial/National/BCR/BCR_Terrestrial_master.shp")  %>%
  subset(PROVINCE_S == "SASKATCHEWAN") %>%
  st_make_valid() %>%
  st_union() %>%
  st_transform(st_crs(all_surveys))

# Create spatial folds
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

# *********************************************************************
# Fit data to each subset, predict into withheld blocks
# *********************************************************************

set.seed(222)

# Prepare data for this species
PC_surveys <- all_surveys %>% subset(Survey_Type %in% c("Point_Count","ARU_SPT","ARU_SPM"))
PC_counts <- full_count_matrix[PC_surveys$Obs_Index,]
rownames(PC_counts) <- PC_surveys$sq_id

n_detections <- PC_counts %>% 
  reshape2::melt() %>%
  dplyr::rename(sq_id = Var1, Species_Code_BSC = Var2, detected = value) %>%
  subset(detected>0) %>%
  group_by(Species_Code_BSC) %>%
  summarize(n_squares = length(unique(sq_id)),
            n_detections = sum(detected>0)) %>%
  subset(Species_Code_BSC %in% species_to_model$Species_Code_BSC)

species_for_crossval <- n_detections %>% 
  subset(n_squares >=50 & n_detections >= 100) %>%
  sample_n(50) %>%
  left_join(species_to_model %>% dplyr::select(-n_squares,-n_detections))%>%
  arrange(n_squares,n_detections) %>%
  as.data.frame()

results <- data.frame()
results_path <- "../Output/Crossvalidation/Crossval_results_Default.rds"

for (sp_code in species_for_crossval$Species_Code_BSC){
  for (xval_fold in 1:n_folds){
    
    print(sp_code)
    if (file.exists(results_path)){
      results <- readRDS(results_path)
      
      if (nrow(subset(results, Species == sp_code & Crossval_Fold == xval_fold))>0) next
    }
    
    # Prepare data for this species
    sp_dat <- all_surveys %>% mutate(count = full_count_matrix[,sp_code]) %>%
      st_intersection(Crossval_Grid)
    
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
    # Standardize 'time since sunrise' covariate (TSS)
    # ------------------------------------------------
    
    sp_dat$Hours_Since_Sunrise <- as.data.frame(sp_dat$Hours_Since_Sunrise)[,1] %>% scale()
    
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
    # Separate data types
    # ------------------------------------------------
    sp_dat$presence <- as.numeric(sp_dat$count>0)
    
    validation_data <- subset(sp_dat, Crossval_Fold == xval_fold) %>%
      subset(Survey_Type %in% c("Point_Count","ARU_SPT","ARU_SPM")) %>% 
      as('Spatial')
    
    training_data <- subset(sp_dat, Crossval_Fold != xval_fold)
    training_data_PC <- subset(training_data,Survey_Type %in% c("Point_Count","ARU_SPT","ARU_SPM")) %>% as('Spatial')
    training_data_SC <- subset(training_data,Survey_Type %in% c("Breeding Bird Atlas") & IncludesPointCounts == 0) %>% as('Spatial')
    training_data_LT <- subset(training_data,Survey_Type %in% c("Linear transect") & IncludesPointCounts == 0) %>%
      mutate(Travel_Speed_mps = Travel_Distance_Metres/(Survey_Duration_Minutes*60)) %>%
      subset(Travel_Distance_Metres <= 5000 & Travel_Speed_mps <= 1.5) %>%
      as('Spatial')
    
    (training_data_SC$count>0) %>% sum()
    (training_data_LT$count>0) %>% sum()
    
    # *********************************************************************
    # FIT WITH INLA
    # *********************************************************************
    
    # Use first 8 PCA axes
    covariates_to_include <- paste0("PC",1:8)
    
    # ------------------------------------------------
    # Create a spatial mesh, which is used to fit the residual spatial field
    # ------------------------------------------------
    
    # Note: mesh developed using tutorial at: https://rpubs.com/jafet089/886687
    max.edge = diff(range(st_coordinates(sp_dat)[,1]))/10
    bound.outer = diff(range(st_coordinates(sp_dat)[,1]))/3
    cutoff = max.edge/5
    bound.outer = diff(range(st_coordinates(sp_dat)[,1]))/3
    
    mesh_spatial <- inla.mesh.2d(loc = st_coordinates(sp_dat),
                                 cutoff = max.edge/2,
                                 max.edge = c(1,2)*max.edge,
                                 offset=c(max.edge, bound.outer))
    mesh_locs <- mesh_spatial$loc[,c(1,2)] %>% as.data.frame()
    
    prior_range <- c(300000,0.1)
    prior_sigma <- c(0.5,0.1)
    matern_coarse <- inla.spde2.pcmatern(mesh_spatial,
                                         prior.range = prior_range, # 10% chance range is smaller than 100000
                                         prior.sigma = prior_sigma # 50% chance sd is larger than 1
    )
    
    
    # ------------------------------------------------
    # Create mesh to model effect of time since sunrise (TSS)
    # ------------------------------------------------
    
    TSS_range <- range(sp_dat$Hours_Since_Sunrise)
    TSS_meshpoints <- seq(TSS_range[1]-0.1,TSS_range[2]+0.1,length.out = 11)
    TSS_mesh1D = inla.mesh.1d(TSS_meshpoints,boundary="free")
    TSS_spde = inla.spde2.pcmatern(TSS_mesh1D,
                                   prior.range = c(1,0.5),
                                   prior.sigma = c(1,0.1)) # 10% chance sd is larger than 1
    
    # ------------------------------------------------
    # Create mesh to model effect of checklist effort
    # ------------------------------------------------
    
    SC_duration_meshpoints <- seq(min(training_data_SC$Survey_Duration_Minutes)/1.2,max(training_data_SC$Survey_Duration_Minutes)*1.2,length.out = 11)
    SC_duration_mesh1D = inla.mesh.1d(SC_duration_meshpoints,boundary="free")
    SC_duration_spde = inla.spde2.pcmatern(SC_duration_mesh1D,
                                           prior.range = c(10,0.5),
                                           prior.sigma = c(1,0.1)) # 10% chance sd is larger than 1
    
    LT_speed_meshpoints <- seq(min(training_data_LT$Travel_Speed_mps)/1.2,max(training_data_LT$Travel_Speed_mps)*1.2,length.out = 11)
    LT_speed_mesh1D = inla.mesh.1d(LT_speed_meshpoints,boundary="free")
    LT_speed_spde = inla.spde2.pcmatern(LT_speed_mesh1D,
                                        prior.range = c(0.2,0.5),
                                        prior.sigma = c(1,0.1)) # 10% chance sd is larger than 1
    
    LT_distance_meshpoints <- seq(min(training_data_LT$Travel_Distance_Metres)/1.2,max(training_data_LT$Travel_Distance_Metres)*1.2,length.out = 11)
    LT_distance_mesh1D = inla.mesh.1d(LT_distance_meshpoints,boundary="free")
    LT_distance_spde = inla.spde2.pcmatern(LT_distance_mesh1D,
                                           prior.range = c(2000,0.5),
                                           prior.sigma = c(1,0.1)) # 10% chance sd is larger than 1
    
    # ------------------------------------------------
    # Fit to point counts + ARU only
    # ------------------------------------------------
    
    # model_components = as.formula(paste0('~
    #         Intercept_PC(1)+
    #         range_effect(1,model="linear", mean.linear = -0.046, prec.linear = 10000)+
    #         TSS(main = Hours_Since_Sunrise,model = TSS_spde) +
    #         spde_coarse(main = coordinates, model = matern_coarse) +',
    #                                      paste0("Beta1_",covariates_to_include,'(1,model="linear", mean.linear = 0, prec.linear = 100)', collapse = " + "),
    #                                      '+',
    #                                      paste0("Beta2_",covariates_to_include,'(1,model="linear", mean.linear = 0, prec.linear = 10000)', collapse = " + ")))
    # 
    sd_linear <- 0.2
    prec_linear = c(1/sd_linear^2,1/(sd_linear/2)^2)
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
    
    start <- Sys.time()
    fit_training_PC <- NULL
    while(is.null(fit_training_PC)){
      
      fit_training_PC <- try(bru(components = model_components,
                             
                             like(family = "nbinomial",
                                  formula = model_formula_PC,
                                  data = training_data_PC),
                             options = list(bru_verbose = 4)))
      if ("try-error" %in% class(fit_training_PC)) fit_training_PC <- NULL
    }
    
    end <- Sys.time() 
    print(end-start) # 1 min
    
    # ------------------------------------------------
    # Formula for making predictions
    # ------------------------------------------------
    
    pred_formula_PC = as.formula(paste0(' ~
                  Intercept_PC +
                  TSS +
                  log_QPAD_offset +
                  range_effect * distance_from_range +
                  spde_coarse +',
                                        paste0("Beta1_",covariates_to_include,'*',covariates_to_include, collapse = " + "),
                                        '+',
                                        paste0("Beta2_",covariates_to_include,'*',covariates_to_include, collapse = " + ")))
    
    # ------------------------------------------------
    # Predict
    # ------------------------------------------------
    
    start2 <- Sys.time()
    pred <- NULL
    pred <- generate(fit_training_PC,
                     validation_data,
                     formula =  pred_formula_PC,
                     n.samples = 500)
    end2 <- Sys.time()
    print(end2-start2) # 20 sec
    
    pred <- exp(pred)
    pred_INLA_PC <- apply(pred,1,function(x) median(x))
    RMSE_INLA_PC <- sqrt(mean((validation_data$count - pred_INLA_PC)^2))
    
    # Parameters defining the negative binomial distribution
    size <- fit_training_PC$summary.hyperpar$'0.5quant'[1]
    prob <- size/(size + pred_INLA_PC)
    
    # Probability of zero
    prob_zero_PC <- dnbinom(0,mu=pred_INLA_PC ,size=size)
    prob_nonzero_PC <- 1-prob_zero_PC
    AUC_INLA_PC <- auc(response = validation_data$presence,predictor = prob_nonzero_PC)
    
    lppd_INLA_PC <- sum(dnbinom(validation_data$count,mu=pred_INLA_PC ,size=size,log = TRUE))
    
    RMSE_INLA_PC
    AUC_INLA_PC
    
    # # ------------------------------------------------
    # # Fit to point counts + SC
    # # ------------------------------------------------
    # 
    # 
    # RMSE_INLA_PC_SC <- NA
    # AUC_INLA_PC_SC <- NA
    # ndet_SC <- sum(training_data_SC$count>0)
    # 
    # # Only fit integrated model if there are at least 10 additional detections
    # if (ndet_SC>10){
    #   model_components = as.formula(paste0('~
    #         Intercept_PC(1)+
    #         Intercept_SC(1)+
    #         range_effect(1,model="linear", mean.linear = -0.046, prec.linear = 10000)+
    #         TSS(main = Hours_Since_Sunrise,model = TSS_spde) +
    #         SC_duration(main = Survey_Duration_Minutes,model = SC_duration_spde) +
    #         spde_coarse(main = coordinates, model = matern_coarse) +',
    #                                        paste0("Beta1_",covariates_to_include,'(1,model="linear", mean.linear = 0, prec.linear = 100)', collapse = " + "),
    #                                        '+',
    #                                        paste0("Beta2_",covariates_to_include,'(1,model="linear", mean.linear = 0, prec.linear = 10000)', collapse = " + ")))
    #   
    #   model_formula_SC = as.formula(paste0('presence ~ log(1/exp(-exp(
    #               Intercept_SC +
    #               spde_coarse +
    #               TSS +
    #               range_effect * distance_from_range +
    #               SC_duration +',
    #                                        paste0("Beta1_",covariates_to_include,'*',covariates_to_include, collapse = " + "),
    #                                        '+',
    #                                        paste0("Beta2_",covariates_to_include,'*',covariates_to_include, collapse = " + "),"))-1)"))
    #   
    #   start <- Sys.time()
    #   
    #   fit_training_PC_SC <- NULL
    #   while(is.null(fit_training_PC_SC)){
    #     
    #     fit_training_PC_SC <- try(bru(components = model_components,
    #                                   
    #                                   like(family = "nbinomial",
    #                                        formula = model_formula_PC,
    #                                        data = training_data_PC),
    #                                   
    #                                   like(family = "binomial",
    #                                        formula = model_formula_SC,
    #                                        data = training_data_SC),
    #                                   
    #                                   options = list(bru_verbose = 4)))
    #     if ("try-error" %in% class(fit_training_PC_SC)) fit_training_PC_SC <- NULL
    #   }
    #   
    #   end <- Sys.time()
    #   print(end-start)
    #   
    #   start2 <- Sys.time()
    #   pred <- NULL
    #   pred <- generate(fit_training_PC_SC,
    #                    validation_data,
    #                    formula =  pred_formula_PC,
    #                    n.samples = 500)
    #   end2 <- Sys.time()
    #   print(end2-start2)
    #   
    #   pred <- exp(pred)
    #   pred_INLA_PC_SC <- apply(pred,1,function(x) median(x))
    #   RMSE_INLA_PC_SC <- sqrt(mean((validation_data$count - pred_INLA_PC_SC)^2))
    #   
    #   # Parameters defining the negative binomial distribution
    #   size_SC <- fit_training_PC_SC$summary.hyperpar$'0.5quant'[1]
    #   prob_SC <- size_SC/(size_SC + pred_INLA_PC_SC)
    #   
    #   # Probability of zero
    #   prob_zero_PC_SC <- dnbinom(0,mu=pred_INLA_PC_SC ,size=size_SC)
    #   prob_nonzero_PC_SC <- 1-prob_zero_PC_SC
    #   AUC_INLA_PC_SC <- auc(response = validation_data$presence,predictor = prob_nonzero_PC_SC)
    #   
    #   RMSE_INLA_PC_SC
    #   AUC_INLA_PC_SC
    # }
    # 
    # # ------------------------------------------------
    # # Fit to point counts + LT
    # # ------------------------------------------------
    # 
    # RMSE_INLA_PC_LT <- NA
    # AUC_INLA_PC_LT <- NA
    # ndet_LT <- sum(training_data_LT$count>0)
    # if (ndet_LT>10){
    #   
    #   # Only fit integrated model if there are at least 10 additional detections
    #   
    #   model_components = as.formula(paste0('~
    #         Intercept_PC(1)+
    #         Intercept_LT(1)+
    #         range_effect(1,model="linear", mean.linear = -0.046, prec.linear = 10000)+
    #         TSS(main = Hours_Since_Sunrise,model = TSS_spde) +
    #         LT_speed(main = Travel_Speed_mps,model = LT_speed_spde) +
    #         LT_distance(main = Travel_Distance_Metres,model = LT_distance_spde) +
    #         spde_coarse(main = coordinates, model = matern_coarse) +',
    #                                        paste0("Beta1_",covariates_to_include,'(1,model="linear", mean.linear = 0, prec.linear = 100)', collapse = " + "),
    #                                        '+',
    #                                        paste0("Beta2_",covariates_to_include,'(1,model="linear", mean.linear = 0, prec.linear = 10000)', collapse = " + ")))
    #   
    #   model_formula_LT = as.formula(paste0('presence ~ log(1/exp(-exp(
    #               Intercept_LT +
    #               spde_coarse +
    #               TSS +
    #               range_effect * distance_from_range +
    #               LT_distance +
    #               LT_speed +',
    #                                        paste0("Beta1_",covariates_to_include,'*',covariates_to_include, collapse = " + "),
    #                                        '+',
    #                                        paste0("Beta2_",covariates_to_include,'*',covariates_to_include, collapse = " + "),"))-1)"))
    #   
    #   start <- Sys.time()
    #   fit_training_PC_LT <- NULL
    #   while(is.null(fit_training_PC_LT)){
    #     
    #     fit_training_PC_LT <- try(bru(components = model_components,
    #                                   
    #                                   like(family = "nbinomial",
    #                                        formula = model_formula_PC,
    #                                        data = training_data_PC),
    #                                   
    #                                   like(family = "binomial",
    #                                        formula = model_formula_LT,
    #                                        data = training_data_LT),
    #                                   
    #                                   options = list(bru_verbose = 4)))
    #     if ("try-error" %in% class(fit_training_PC_LT)) fit_training_PC_LT <- NULL
    #   }
    #   end <- Sys.time()
    #   print(end-start)
    #   
    #   start2 <- Sys.time()
    #   pred <- NULL
    #   pred <- generate(fit_training_PC_LT,
    #                    validation_data,
    #                    formula =  pred_formula_PC,
    #                    n.samples = 500)
    #   end2 <- Sys.time()
    #   print(end2-start2)
    #   
    #   pred <- exp(pred)
    #   pred_INLA_PC_LT <- apply(pred,1,function(x) median(x))
    #   RMSE_INLA_PC_LT <- sqrt(mean((validation_data$count - pred_INLA_PC_LT)^2))
    #   
    #   # Parameters defining the negative binomial distribution
    #   size_LT <- fit_training_PC_LT$summary.hyperpar$'0.5quant'[1]
    #   prob_LT <- size_LT/(size_LT + pred_INLA_PC_LT)
    #   
    #   # Probability of zero
    #   prob_zero_PC_LT <- dnbinom(0,mu=pred_INLA_PC_LT ,size=size_LT)
    #   prob_nonzero_PC_LT <- 1-prob_zero_PC_LT
    #   AUC_INLA_PC_LT <- auc(response = validation_data$presence,predictor = prob_nonzero_PC_LT)
    #   
    #   RMSE_INLA_PC_LT
    #   AUC_INLA_PC_LT
    # }
    
    # *********************************************************************
    # FIT WITH BRT
    # *********************************************************************

    start3 <- Sys.time()
    BRT_covariates <- as.data.frame(training_data_PC) %>%
      dplyr::select(elevation_1km:PC10) %>%
      colnames()
    BRT_covariates <- c(BRT_covariates,"Hours_Since_Sunrise","distance_from_range","in_range")

    training_data_PC_df <- as.data.frame(training_data_PC)
    fit_training_BRT <-  try(fit_brt(model_data =  training_data_PC_df,
                                response_column = which(colnames( training_data_PC_df)=="count"),
                                covariate_columns = which(colnames( training_data_PC_df)%in% BRT_covariates)))
    if ("try-error" %in% class(fit_training_BRT)){
      fit_training_BRT <- NA
      RMSE_BRT <- NA
      AUC_BRT <- NA
      lppd_BRT <- NA
    } else{
      
    # Join with grid (for plotting)
    validation_df <- as.data.frame(validation_data)
    pred_BRT  <- predict(fit_training_BRT,
                         validation_df,
                         n.trees = fit_training_BRT$gbm.call$best.trees)

    pred_BRT <- exp(pred_BRT + validation_data$log_QPAD_offset)
    RMSE_BRT <- sqrt(mean((validation_df$count - pred_BRT)^2))
    prob_zero_BRT <- exp(-pred_BRT)
    prob_nonzero_BRT <- 1-prob_zero_BRT
    AUC_BRT <- auc(response = validation_data$presence,predictor = prob_nonzero_BRT)
    lppd_BRT <- sum(dpois(validation_data$count,lambda=pred_BRT,log = TRUE))
    }
    
    end3 <- Sys.time()
    
    # *********************************************************************
    # Save results
    # *********************************************************************
    
    run_results <- data.frame(Species = sp_code,
                              Crossval_Fold = xval_fold,
                              
                              RMSE_BRT = RMSE_BRT,
                              RMSE_INLA_PC = RMSE_INLA_PC,
                              #RMSE_INLA_PC_SC = RMSE_INLA_PC_SC,
                              #RMSE_INLA_PC_LT = RMSE_INLA_PC_LT,
                              
                              AUC_BRT = AUC_BRT,
                              AUC_INLA_PC = AUC_INLA_PC,
                              #AUC_INLA_PC_SC = AUC_INLA_PC_SC,
                              #AUC_INLA_PC_LT = AUC_INLA_PC_LT,
                              
                              lppd_BRT = lppd_BRT,
                              lppd_INLA_PC = lppd_INLA_PC,
                              
                              prior_range = paste0(prior_range, collapse = "_"),
                              prior_sigma = paste0(prior_sigma, collapse = "_"),
                              sd_linear = sd_linear
    )
    
    
    if (file.exists(results_path)) results <- readRDS(results_path)
    
    results <- rbind(results, run_results)
    print(results)
    saveRDS(results,results_path)
    
    # *********************************************************************
    # Plot fit comparisons
    # *********************************************************************
    # 
    # if (file.exists(results_path)) results <- readRDS(results_path)
    # 
    # results_summary <- results %>%
    #   group_by(Species) %>%
    #   summarize(AUC_INLA = mean(AUC_INLA_PC,na.rm = TRUE),
    #             RMSE_INLA = mean(RMSE_INLA_PC,na.rm = TRUE)
    #             #AUC_BRT = mean(AUC_BRT,na.rm = TRUE),
    #             #RMSE_BRT = mean(RMSE_BRT,na.rm = TRUE)
    #   ) #%>%
    #  mutate(percent_change_RMSE = 100*(RMSE_INLA - RMSE_BRT)/RMSE_BRT)
    
    # # Does INLA improve AUC relative to BRT?
    # fitplot1 <- ggplot(results_summary)+
    #   
    #   
    #   geom_segment(aes(y = Species, yend = Species,
    #                    x = AUC_BRT, xend = AUC_INLA,
    #                    col = factor(AUC_BRT > AUC_INLA)),
    #                size = 3,
    #                arrow = arrow(length = unit(0.5, "cm")))+
    #   
    #   
    #   theme_bw()+
    #   xlab("Cross-validation AUC")+
    #   scale_color_manual(values = c("dodgerblue","orangered"), guide = "none")+
    #   ggtitle("Comparison of AUC\n\nINLA relative to BRT")
    # 
    # # Does INLA reduce RMSE relative to BRT?
    # rng <- max(abs(results_summary$percent_change_RMSE))
    # fitplot2 <- ggplot(results_summary)+
    #   
    #   
    #   geom_segment(aes(y = Species, yend = Species,
    #                    x = 0, xend = percent_change_RMSE,
    #                    col = factor(percent_change_RMSE>0)),
    #                size = 3,
    #                arrow = arrow(length = unit(0.5, "cm")))+
    #   coord_cartesian(xlim=c(-rng,rng))+
    #   
    #   theme_bw()+
    #   xlab("% Difference in RMSE")+
    #   scale_color_manual(values = c("dodgerblue","orangered"), guide = "none")+
    #   ggtitle("Comparison of RMSE\n\nINLA relative to BRT")
    # 
    # fitplot <- ggarrange(fitplot1,fitplot2,nrow=1)
    # print(fitplot)
  }
  
} # close species loop


# ----------------------------------------------------------------
# Default settings
# ----------------------------------------------------------------

results_path <- "../Output/Crossvalidation/Crossval_results_Default.rds"
if (file.exists(results_path)) results_Default <- readRDS(results_path)

results_summary_Default <- results_Default %>%
  group_by(Species) %>%
  summarize(AUC_INLA_Default = mean(AUC_INLA_PC,na.rm = TRUE),
            RMSE_INLA_Default = mean(RMSE_INLA_PC,na.rm = TRUE),
            lppd_INLA_Default = mean(lppd_INLA_PC,na.rm = TRUE))

# ----------------------------------------------------------------
# Increase SD of beta parameters
# ----------------------------------------------------------------

results_path <- "../Output/Crossvalidation/Crossval_results_IncreaseBeta.rds"
if (file.exists(results_path)) results_IncreaseBeta <- readRDS(results_path)

results_summary_IncreaseBeta <- results_IncreaseBeta %>%
  group_by(Species) %>%
  summarize(AUC_INLA_IncreaseBeta = mean(AUC_INLA_PC,na.rm = TRUE),
            RMSE_INLA_IncreaseBeta = mean(RMSE_INLA_PC,na.rm = TRUE),
            lppd_INLA_IncreaseBeta = mean(lppd_INLA_PC,na.rm = TRUE))


# ----------------------------------------------------------------
# Join results
# ----------------------------------------------------------------
results_summary <- full_join(results_summary_Default,
                             results_summary_IncreaseBeta) 
 
mean(results_summary$AUC_INLA_IncreaseBeta > results_summary$AUC_INLA_Default, na.rm = TRUE)
mean(results_summary$RMSE_INLA_IncreaseBeta < results_summary$RMSE_INLA_Default, na.rm = TRUE)
mean(results_summary$lppd_INLA_IncreaseBeta > results_summary$lppd_INLA_Default, na.rm = TRUE)

range(log(results_summary$AUC_INLA_IncreaseBeta/results_summary$AUC_INLA_Default), na.rm = TRUE)
range(log(results_summary$RMSE_INLA_IncreaseBeta/results_summary$RMSE_INLA_Default), na.rm = TRUE)
range(results_summary$lppd_INLA_IncreaseBeta-results_summary$lppd_INLA_Default, na.rm = TRUE)

hist(log(results_summary$RMSE_INLA_IncreaseBeta/results_summary$RMSE_INLA_Default), breaks = 20)
hist(log(results_summary$AUC_INLA_IncreaseBeta/results_summary$AUC_INLA_Default), breaks = 20)
hist(results_summary$lppd_INLA_IncreaseBeta-results_summary$lppd_INLA_Default, breaks = 20)
