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
library(viridis)

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
    } else{
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

all_surveys <- all_surveys %>%
  subset(Survey_Type %in% c("Point_Count","ARU_SPT","ARU_SPM"))
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

species_for_crossval<- data.frame(Number_of_Detections = colSums(full_count_matrix>0),
                                  Species_Code_BSC = names(colSums(full_count_matrix>0))) %>%
  subset(Number_of_Detections>500) %>%
  subset(Species_Code_BSC %in% analysis_data$species_to_model$Species_Code_BSC) %>%
  left_join(.,analysis_data$species_to_model %>% dplyr::select(-Number_of_Detections)) %>%
  sample_n(25)%>% 
  arrange(desc(Number_of_Detections))
species_for_crossval

results <- data.frame()
results_path <- "../Output/Crossvalidation/Crossval_results_spline_nbinomial.rds"

for (xval_fold in rev(1:n_folds)){
  for (sp_code in rev(species_for_crossval$Species_Code_BSC)){
    
    if (file.exists(results_path)){
      results <- readRDS(results_path)
      
      if (nrow(subset(results, Species == sp_code & Crossval_Fold == xval_fold))>0) next
    }
    
    print(sp_code)
    
    # Prepare data for this species
    sp_dat <- all_surveys %>% mutate(count = full_count_matrix[,sp_code]) %>%
      # Only select point counts
      subset(Survey_Type %in% c("Point_Count","ARU_SPT","ARU_SPM")) %>%
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
    
    training_data <- subset(sp_dat, Crossval_Fold != xval_fold)
    validation_data <- subset(sp_dat, Crossval_Fold == xval_fold)
    
    # *********************************************************************
    # FIT WITH INLA
    # *********************************************************************
    
    {
      {
        {
          covariates_to_include <- c("PC1","PC2","PC3","PC4","PC5","PC6","PC7")
          
          # ------------------------------------------------
          # Convert to spatial objects
          # ------------------------------------------------
          
          training_data_sp <- as(training_data,'Spatial')
          validation_data_sp <- as(validation_data,'Spatial')
          
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
          #plot(mesh_spatial)
          
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
                                         prior.range = c(1,0.5),
                                         prior.sigma = c(1,0.1)) # 10% chance sd is larger than 1
          
          # ------------------------------------------------
          # Create mesh to model effect of time since sunrise (TSS)
          # ------------------------------------------------
          
          covar_range <- range(as.data.frame(sp_dat)[,covariates_to_include])
          covar_range[1] <- min(-10,covar_range[1])
          covar_range[2] <- max(10,covar_range[2])
          
          covar_meshpoints <- c(seq(covar_range[1],-3,1),seq(-3,3,0.2),seq(3,covar_range[2],1)) %>% unique() %>% sort()
          covar_mesh1D = inla.mesh.1d(covar_meshpoints,boundary="free")
          covar_spde = inla.spde2.pcmatern(covar_mesh1D,
                                           prior.range = c(1,0.5),
                                           prior.sigma = c(1,0.1)) # 10% chance sd is larger than 1
          
          # ------------------------------------------------
          # Model formulas
          # ------------------------------------------------
          
          model_components = as.formula(paste0('~
            Intercept_PC(1)+
            range_effect(1,model="linear", mean.linear = -0.046, prec.linear = 10000)+
            TSS(main = Hours_Since_Sunrise,model = TSS_spde) +
            spde_coarse(main = coordinates, model = matern_coarse) +',
                                               paste0(covariates_to_include,'_effect(main = ',covariates_to_include,",model=covar_spde)", collapse = " + ")))
          
          model_formula_PC = as.formula(paste0('count ~
                  Intercept_PC +
                  log_QPAD_offset +
                  TSS +
                  range_effect * distance_from_range +
                  spde_coarse +',
                                               paste0(covariates_to_include,"_effect", collapse = " + ")))
          
          pred_formula_PC = as.formula(paste0(' ~
                  Intercept_PC +
                  log_QPAD_offset +
                  TSS + 
                  range_effect * distance_from_range +
                  spde_coarse +',
                                              paste0(covariates_to_include,"_effect", collapse = " + ")))
          
          # ------------------------------------------------
          # Fit with nbinomial error
          # ------------------------------------------------
          
          start <- Sys.time()
          fit_training_INLA <- bru(components = model_components,
                                   like(family = "nbinomial",
                                        formula = model_formula_PC,
                                        data = training_data_sp),
                                   options = list(bru_verbose = 4))
          
          pred <- NULL
          pred <- generate(fit_training_INLA,
                           validation_data_sp,
                           formula =  pred_formula_PC,
                           n.samples = 500)
          
          pred <- exp(pred)
          pred_INLA <- apply(pred,1,function(x) median(x))
          RMSE_INLA_nbinomial <- sqrt(mean((validation_data$count - pred_INLA)^2))
          end <- Sys.time()
          print(end-start)
          
        }
      }
    }
    
    # *********************************************************************
    # FIT WITH BRT
    # *********************************************************************
    
    BRT_covariates <- as.data.frame(all_surveys) %>%
      dplyr::select(elevation_1km:PC20) %>%
      colnames()
    
    BRT_covariates <- c(BRT_covariates,"Hours_Since_Sunrise")
    
    training_df <- as.data.frame(training_data)
    fit_training_BRT <- fit_brt(model_data = training_df,
                                response_column = which(colnames(training_df)=="count"),
                                covariate_columns = which(colnames(training_df)%in% BRT_covariates))
    
    # Join with grid (for plotting)
    validation_df <- as.data.frame(validation_data)
    pred_BRT  <- predict(fit_training_BRT, 
                         validation_df, 
                         n.trees = fit_training_BRT$gbm.call$best.trees)
    
    pred_BRT <- exp(pred_BRT + validation_data$log_QPAD_offset)
    RMSE_BRT <- sqrt(mean((validation_df$count - pred_BRT)^2))
    
    # *********************************************************************
    # Save results
    # *********************************************************************
    
    run_results <- data.frame(Species = sp_code,
                              Crossval_Fold = xval_fold,
                              RMSE_INLA_poisson = RMSE_INLA_poisson,
                              RMSE_INLA_nbinomial = RMSE_INLA_nbinomial,
                              RMSE_BRT = RMSE_BRT)
    
    
    if (file.exists(results_path)) results <- readRDS(results_path)
    
    results <- rbind(results, run_results)
    print(results)
    saveRDS(results,results_path)
 
  }
  
} # close species loop



# *********************************************************************
# Plot fit comparisons
# *********************************************************************

if (file.exists(results_path)) results <- readRDS(results_path)

results_summary <- results %>%
  group_by(Species) %>%
  summarize(RMSE_INLA = mean(RMSE_INLA_nbinomial),
            RMSE_BRT = mean(RMSE_BRT),

            percent_change_RMSE = 100*((mean(RMSE_INLA) - mean(RMSE_BRT))/mean(RMSE_BRT)))

fitplot <- ggplot(results_summary)+
  geom_bar(aes(y = Species, x = percent_change_RMSE, fill = factor(sign(percent_change_RMSE))), stat = "identity")+
  theme_bw()+
  xlab("% change")+
  scale_fill_manual(values = c("dodgerblue","orangered"), guide = "none")+
  coord_cartesian(xlim=c(-50,50))+
  ggtitle("Percent Change in Root Mean Squared Error\n\nINLA relative to BRT")
print(fitplot)
