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
library(napops)

# Timeout INLA after 4 minutes (if it has not fit by then, it has likely stalled)
inla.setOption(inla.timeout = 60*10)

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
  subset(n_squares >= 50 & n_detections >= 100) %>%
  left_join(species_to_model %>% dplyr::select(-n_squares,-n_detections))%>%
  subset(offset_exists)%>%
  sample_n(50) %>%
  arrange(n_squares,n_detections) %>%
  as.data.frame()

results <- data.frame()
results_path <- "../Output/Crossvalidation/Crossval_results_TSS2.rds"

for (sp_code in rev(species_for_crossval$Species_Code_BSC)){
  
  # Prepare data for this species
  sp_dat <- all_surveys %>% mutate(count = full_count_matrix[,sp_code]) %>%
    st_intersection(Crossval_Grid) %>%
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
  
  # --------------------------------
  # Generate QPAD offsets for each survey based on NAPOPS
  # --------------------------------
  
  removal_models <- napops::coef_removal(species = sp_code)[c(1,4,5),]
  best <- c(1,4,5)[which.min(removal_models$AIC)]
  
  EDR <- edr(species = sp_code,road = FALSE, forest = 0.5,model = 1)[3] %>% as.numeric()
  A_metres <- pi*EDR^2
  
  cue_rate_TSS <- cue_rate(species = sp_code,
                           model=best,
                           od = 153,
                           tssr=as.numeric(sp_dat$Hours_Since_Sunrise))$CR_est
  
  p <- 1-exp(-sp_dat$Survey_Duration_Minutes*cue_rate_TSS)
  sp_dat$log_QPAD_offset_TSS <- log(A_metres * p)
  
  # Assume intercept-only
  cue_rate_intercept <- cue_rate(species = sp_code,
                                 model=1,
                                 od = 153,
                                 tssr=0)$CR_est
  
  p <- 1-exp(-sp_dat$Survey_Duration_Minutes*cue_rate_intercept)
  sp_dat$log_QPAD_offset <- log(A_metres * p)
  
  # # ------------------------------------------------
  # # Standardize 'time since sunrise' covariate (TSS)
  # # ------------------------------------------------
  # sp_dat$TSS_unscaled <- as.data.frame(sp_dat$Hours_Since_Sunrise)[,1] %>% as.numeric()
  # TSS_mean <- mean(sp_dat$TSS_unscaled)
  # TSS_sd <- sd(sp_dat$TSS_unscaled)
  # sp_dat$Hours_Since_Sunrise <- as.data.frame(sp_dat$Hours_Since_Sunrise)[,1] %>% scale()
  # 
  for (xval_fold in 1:n_folds){
    
    print(sp_code)
    if (file.exists(results_path)){
      results <- readRDS(results_path)
      
      if (nrow(subset(results, Species == sp_code & Crossval_Fold == xval_fold))>0) next
    }
    
    # ------------------------------------------------
    # Separate data types
    # ------------------------------------------------
    sp_dat$presence <- as.numeric(sp_dat$count>0)
    
    validation_data <- subset(sp_dat, Crossval_Fold == xval_fold) %>% as('Spatial')
    
    training_data <- subset(sp_dat, Crossval_Fold != xval_fold)
    training_data_PC <- training_data %>% as('Spatial')
    
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
                                   prior.range = c(6,0.1),
                                   prior.sigma = c(1,0.1)) # 10% chance sd is larger than 1
    
    
    sd_linear <- 0.2
    prec_linear = c(1/sd_linear^2,1/(sd_linear/2)^2)
    
    # ********************************************************
    # ********************************************************
    # Fit using flexible "time since sunrise" effect
    # ********************************************************
    # ********************************************************
    
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
    fit_METHOD1 <- NULL
    while(is.null(fit_METHOD1)){
      
      fit_METHOD1 <- try(bru(components = model_components,
                             
                             like(family = "nbinomial",
                                  formula = model_formula_PC,
                                  data = training_data_PC),
                             options = list(control.compute = list(waic = TRUE),
                                            bru_verbose = 4)))
      if ("try-error" %in% class(fit_METHOD1)) fit_METHOD1 <- NULL
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
    
    pred_METHOD1 <- NULL
    pred_METHOD1 <- generate(fit_METHOD1,
                             validation_data,
                             formula =  pred_formula_PC,
                             n.samples = 1000)
    
    pred_METHOD1 <- exp(pred_METHOD1)
    pred_METHOD1 <- apply(pred_METHOD1,1,function(x) median(x))
    
    # Parameters defining the negative binomial distribution
    size <- fit_METHOD1$summary.hyperpar$'0.5quant'[1]
    prob <- size/(size + pred_METHOD1)
    
    # Probability of zero
    prob_zero_PC <- dnbinom(0,mu=pred_METHOD1 ,size=size)
    prob_nonzero_PC <- 1-prob_zero_PC
    AUC_METHOD1 <- auc(response = validation_data$presence,predictor = prob_nonzero_PC)
    lppd_METHOD1 <- sum(dnbinom(validation_data$count,mu=pred_METHOD1 ,size=size,log = TRUE))
    WAIC_METHOD1 <- fit_METHOD1$waic$waic
    
    # ********************************************************
    # ********************************************************
    # Fit using only NAPOPS correction
    # ********************************************************
    # ********************************************************
    
    model_components = as.formula(paste0('~
            Intercept_PC(1)+
            range_effect(1,model="linear", mean.linear = -0.046, prec.linear = 10000)+
            spde_coarse(main = coordinates, model = matern_coarse) +',
                                         paste0("Beta1_",covariates_to_include,'(1,model="linear", mean.linear = 0, prec.linear = ', prec_linear[1],')', collapse = " + "),
                                         '+',
                                         paste0("Beta2_",covariates_to_include,'(1,model="linear", mean.linear = 0, prec.linear = ', prec_linear[2],')', collapse = " + ")))
    
    
    model_formula_PC = as.formula(paste0('count ~
                  Intercept_PC +
                  log_QPAD_offset_TSS +
                  range_effect * distance_from_range +
                  spde_coarse +',
                                         paste0("Beta1_",covariates_to_include,'*',covariates_to_include, collapse = " + "),
                                         '+',
                                         paste0("Beta2_",covariates_to_include,'*',covariates_to_include, collapse = " + ")))
    
    start <- Sys.time()
    fit_METHOD2 <- NULL
    while(is.null(fit_METHOD2)){
      
      fit_METHOD2 <- try(bru(components = model_components,
                             
                             like(family = "nbinomial",
                                  formula = model_formula_PC,
                                  data = training_data_PC),
                             options = list(control.compute = list(waic = TRUE),
                                            bru_verbose = 4)))
      if ("try-error" %in% class(fit_METHOD2)) fit_METHOD2 <- NULL
    }
    
    end <- Sys.time() 
    print(end-start) # 1 min
    
    # ------------------------------------------------
    # Formula for making predictions
    # ------------------------------------------------
    
    pred_formula_PC = as.formula(paste0(' ~
                  Intercept_PC +
                  log_QPAD_offset_TSS +
                  range_effect * distance_from_range +
                  spde_coarse +',
                                        paste0("Beta1_",covariates_to_include,'*',covariates_to_include, collapse = " + "),
                                        '+',
                                        paste0("Beta2_",covariates_to_include,'*',covariates_to_include, collapse = " + ")))
    
    # ------------------------------------------------
    # Predict
    # ------------------------------------------------
    
    pred_METHOD2 <- NULL
    pred_METHOD2 <- generate(fit_METHOD2,
                             validation_data,
                             formula =  pred_formula_PC,
                             n.samples = 1000)
    
    pred_METHOD2 <- exp(pred_METHOD2)
    pred_METHOD2 <- apply(pred_METHOD2,1,function(x) median(x))
    
    # Parameters defining the negative binomial distribution
    size <- fit_METHOD2$summary.hyperpar$'0.5quant'[1]
    prob <- size/(size + pred_METHOD2)
    
    # Probability of zero
    prob_zero_PC <- dnbinom(0,mu=pred_METHOD2 ,size=size)
    prob_nonzero_PC <- 1-prob_zero_PC
    AUC_METHOD2 <- auc(response = validation_data$presence,predictor = prob_nonzero_PC)
    lppd_METHOD2 <- sum(dnbinom(validation_data$count,mu=pred_METHOD2 ,size=size,log = TRUE))
    WAIC_METHOD2 <- fit_METHOD2$waic$waic
    
    lppd_METHOD1 - lppd_METHOD2
    
    # ********************************************************
    # ********************************************************
    # Plot estimated relationships with Time Since Sunrise
    # ********************************************************
    # ********************************************************
    # 
    # pred_TSSR <- data.frame(Hours_Since_Sunrise = seq(min(sp_dat$Hours_Since_Sunrise) %>% as.numeric(),
    #                                                   max(sp_dat$Hours_Since_Sunrise) %>% as.numeric(),
    #                                                   length.out = 100))
    # 
    # # Offset under intercept-only model
    # cue_rate_intercept <- cue_rate(species = sp_code,
    #                                model=1,
    #                                od = 153,
    #                                tssr=0)$CR_est
    # pred_TSSR$log_QPAD_offset <- log(A_metres * (1-exp(-5*cue_rate_intercept)))[1]
    # 
    # # Offset under best NAPOPS model
    # cue_rate_TSS <- cue_rate(species = sp_code,
    #                            model=best,
    #                            od = 153,
    #                            tssr=as.numeric(pred_TSSR$Hours_Since_Sunrise))$CR_est
    # 
    # pred_TSSR$cue_rate <- cue_rate_TSS
    # pred_TSSR$log_QPAD_offset_TSSR <- log(A_metres * (1-exp(-5*cue_rate_TSS)))
    # 
    # # predictions based on method 1 (flexible)
    # pred1 <- generate(fit_METHOD1,
    #                   pred_TSSR,
    #                   formula =  as.formula(' ~ Intercept_PC + TSS +log_QPAD_offset'),
    #                   n.samples = 1000)
    # 
    # # predictions based on method 2
    # pred2 <- generate(fit_METHOD2,
    #                   pred_TSSR,
    #                   formula =  as.formula(' ~ Intercept_PC + log_QPAD_offset_TSSR'),
    #                   n.samples = 1000)
    # 
    # pred1 <- pred1
    # pred2 <- pred2
    # 
    # # Expected count under method 1
    # pred_TSSR$Pred_METHOD1_q025 <- apply(pred1,1,function(x)quantile(x,0.05))
    # pred_TSSR$Pred_METHOD1_q500 <- apply(pred1,1,function(x)quantile(x,0.5))
    # pred_TSSR$Pred_METHOD1_q975 <- apply(pred1,1,function(x)quantile(x,0.95))
    # 
    # # Expected count under method 2
    # pred_TSSR$Pred_METHOD2_q025 <- apply(pred2,1,function(x)quantile(x,0.05))
    # pred_TSSR$Pred_METHOD2_q500 <- apply(pred2,1,function(x)quantile(x,0.5))
    # pred_TSSR$Pred_METHOD2_q975 <- apply(pred2,1,function(x)quantile(x,0.95))
    # 
    # ggplot(data = pred_TSSR)+
    #   geom_ribbon(aes(x = Hours_Since_Sunrise, ymin = exp(Pred_METHOD1_q025),ymax = exp(Pred_METHOD1_q975)), alpha = 0.5, fill = "dodgerblue")+
    #   geom_line(aes(x = Hours_Since_Sunrise, y = exp(Pred_METHOD1_q500)), col = "dodgerblue")+
    # 
    #   geom_ribbon(aes(x = Hours_Since_Sunrise, ymin = exp(Pred_METHOD2_q025),ymax = exp(Pred_METHOD2_q975)), alpha = 0.5, fill = "orangered")+
    #   geom_line(aes(x = Hours_Since_Sunrise, y = exp(Pred_METHOD2_q500)), col = "orangered")+
    # 
    #   theme_bw()+
    #   ylab("Predicted count")+
    #   xlab("Hours since sunrise")+
    #   scale_y_continuous(label = scales::comma, trans = "log10")+
    #   ggtitle(sp_code)
    # 
    # *********************************************************************
    # Save results
    # *********************************************************************
    
    run_results <- data.frame(Species = sp_code,
                              Crossval_Fold = xval_fold,
                              
                              best_removal_model = c("Null","TSS","TSS2")[which.min(removal_models$AIC)],
                              
                              WAIC_METHOD1 = WAIC_METHOD1,
                              WAIC_METHOD2 = WAIC_METHOD2,
                              
                              AUC_METHOD1 = AUC_METHOD1,
                              AUC_METHOD2 = AUC_METHOD2,
                              lppd_METHOD1 = lppd_METHOD1,
                              lppd_METHOD2 = lppd_METHOD2,
                              
                              
                              prior_range = paste0(prior_range, collapse = "_"),
                              prior_sigma = paste0(prior_sigma, collapse = "_"),
                              sd_linear = sd_linear
    )
    
    
    if (file.exists(results_path)) results <- readRDS(results_path)
    
    results <- rbind(results, run_results)
    print(results)
    saveRDS(results,results_path)
    
  }
  
} # close species loop

if (file.exists(results_path)) results <- readRDS(results_path)

results_summary <- results %>%
  group_by(Species,best_removal_model)%>%
  summarize(lppd_1 = mean(lppd_METHOD1),
            lppd_2 = mean(lppd_METHOD2),
            delta_lppd = mean(lppd_METHOD1 - lppd_METHOD2),
            
            AUC_1 = mean(AUC_METHOD1),
            AUC_2 = mean(AUC_METHOD2),
            delta_AUC = mean(AUC_METHOD1 - AUC_METHOD2))

# Does INLA increase lppd relative to BRT?
fitplot <- ggplot(results_summary)+
  
  geom_segment(aes(y = Species, yend = Species,
                   x = 0, xend = delta_lppd,
                   col = factor(delta_lppd<0)),
               size = 2,
               arrow = arrow(length = unit(0.1, "cm")))+
  
  
  theme_bw()+
  xlab("Change in Deviance\n\nPositive values indicate flexible model is better")+
  scale_color_manual(values = c("dodgerblue","orangered"), guide = "none")+
  ggtitle("Cross-validation Deviance\n\n'Flexible TSSR model' vs. 'Top model' from NAPOPS")
fitplot

# # Does INLA increase AUC relative to BRT?
# fitplot2 <- ggplot(results_summary)+
#   
#   geom_segment(aes(y = Species, yend = Species,
#                    x = 0, xend = delta_AUC,
#                    col = factor(delta_AUC<0)),
#                size = 2,
#                arrow = arrow(length = unit(0.1, "cm")))+
#   
#   
#   theme_bw()+
#   xlab("Change in AUC\n\nPositive values indicate flexible model is better")+
#   scale_color_manual(values = c("dodgerblue","orangered"), guide = "none")+
#   ggtitle("Cross-validation AUC\n\n'Flexible TSSR model' vs. 'Top model' from NAPOPS")
# fitplot2
