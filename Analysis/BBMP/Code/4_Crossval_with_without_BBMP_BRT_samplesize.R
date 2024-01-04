# ------------------------------------------------------------------
# This script evaluates the relative value of BBMP and non-BBMP data,
# while holding sample size constant
# ------------------------------------------------------------------

# ------------------------------------------------
# Load/install packages and set graphical themes / working directory
# ------------------------------------------------

my_packs = c('tidyverse',
             'magrittr',
             'sf',
             'terra',
             'ggspatial',
             'naturecounts',
             'suntools',
             "readxl",
             "viridis",
             "ebirdst",
             "terra",
             "colorspace",
             "gdata",
             "pROC",
             "exactextractr",
             "napops",
             "dismo")

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
# ggplot theme
# ------------------------------------------------

CustomTheme <- theme_set(theme_bw())
CustomTheme <- theme_update(legend.key = element_rect(colour = NA), 
                            legend.key.height = unit(1.2, "line"),
                            panel.grid.major = element_line(colour = 'gray95'),
                            #panel.grid.minor = element_line(colour = 'transparent'),
                            panel.border = element_rect(linetype = "solid",
                                                        colour = "black",
                                                        linewidth = 1, fill = NA),
                            axis.line = element_line(colour = "black"),
                            strip.text = element_text(size = 12, colour = "black"),
                            strip.background = element_rect(colour = "black",
                                                            fill = "lightblue2",
                                                            linetype = "solid"),
                            axis.title.y = element_text(margin = margin(0,10,0,0)),
                            axis.title.x = element_text(margin = margin(10,0,0,0)),
                            panel.background = element_rect(fill = "white"))

# ----------------------------------------------------------------
# Function to fit BRT
# ----------------------------------------------------------------

fit_brt <- function(model_data,
                    response_column = NA,
                    covariate_columns = NA){
  mod_brt <- NULL
  
  ntrees <- 50
  tcomplexity <- 5
  lrate <- 0.1
  m <- 0
  
  while(is.null(mod_brt)){
    
    m <- m + 1
    if(m < 11){
      ntrees <- 50
      lrate <- 0.1
    } else if(m < 21){
      lrate <- 0.05
    } else if(m < 31){
      ntrees <- 25
      lrate <- 0.05
    } else if(m < 41){
      ntrees <- 25
      lrate <- 0.025
    } else if(m < 51){
      ntrees <- 10
      lrate <- 0.025
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


# # -----------------------------------------------
# # Useful functions
# # -----------------------------------------------
# 
# `%!in%` <- Negate(`%in%`)
# 
# AEA_proj <- "+proj=aea +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-100 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs "
# 
# BCR_CAN <- st_read("../../../Data/Spatial/National/BCR/BCR_Terrestrial_master.shp")  %>%
#   subset(COUNTRY %in% "CANADA") %>%
#   st_make_valid() %>%
#   dplyr::select(BCR, PROVINCE_S) %>%
#   st_transform(crs = AEA_proj)
# 
# BCR_NA <- st_read("../../../Data/Spatial/National/BCR/BCR_Terrestrial_master.shp")  %>%
#   subset(COUNTRY %in% c("CANADA","USA")) %>%
#   st_make_valid() %>%
#   dplyr::select(BCR, PROVINCE_S) %>%
#   st_transform(crs = AEA_proj)
# 
# # ******************************************************************
# # LOAD BBMP DATA PACKAGE
# # ******************************************************************
# 
# # -----------------------------------------------
# # Data from WildTrax (prepared by previous script)
# # -----------------------------------------------
# 
# WT_dat <- readRDS(file = "../Data_Cleaned/WildTrax/WT_dat.rds")
# 
# BBMP_boundary <- WT_dat$BBMP_boundary %>%
#   st_transform(crs = AEA_proj) %>%
#   st_intersection(BCR_CAN %>% st_union())
# 
# BBMP_surveys <- WT_dat$WT_surveyinfo %>% st_transform(crs = AEA_proj) %>% mutate(Obs_Index = 1:nrow(.))
# BBMP_counts <- WT_dat$WT_matrix
# 
# # Survey locations
# in_BBMP <- BBMP_surveys %>%
#   st_transform(crs = AEA_proj) %>%
#   st_intersects(BBMP_boundary) %>%
#   as.numeric()
# 
# BBMP_surveys <- BBMP_surveys[!is.na(in_BBMP),]
# BBMP_counts <- BBMP_counts[BBMP_surveys$Obs_Index,]
# 
# # ******************************************************************
# # LOAD BAM DATA PACKAGE, CROP TO BBMP EXTENT
# # ******************************************************************
# 
# bam <- load("../../../Data/Bird_Data_Raw/BAM/04_NM5.0_data_stratify.Rdata")
# 
# bam_surveys <- visit %>%
#   mutate(Obs_Index = 1:nrow(.)) %>%
#   subset(year >= 2013) %>%
#   st_as_sf(coords = c("lon", "lat"),crs=4326, remove = FALSE) %>%
#   st_transform(AEA_proj)
# within_BBMP <- st_intersects(bam_surveys,BBMP_boundary) %>% as.numeric()
# bam_surveys <- bam_surveys[!is.na(within_BBMP),]
# 
# bam_counts <- bird[bam_surveys$Obs_Index,]
# bam_offsets <- offsets[bam_surveys$Obs_Index,]
# 
# # ******************************************************************
# # REMOVE BAM SURVEYS WITHIN 100 metres OF BBMP SURVEYS
# # ******************************************************************
# 
# # Identify surveys in the BAM dataset that are in the BBMP dataset
# # Criteria: within 100 m
# bam_surveys$Obs_Index <- 1:nrow(bam_surveys)
# BBMP_buff <- BBMP_surveys %>% st_buffer(100) %>% st_union()
# within_BBMP <- st_intersects(bam_surveys,BBMP_buff) %>% as.numeric()
# bam_surveys <- bam_surveys[is.na(within_BBMP),]
# bam_counts <- bam_counts[bam_surveys$Obs_Index,]
# 
# # ******************************************************************
# # PLOT DATA AVAILABILITY ACROSS THE BBMP SURVEY REGION
# # PLOT CWS data in a separate color from other surveys
# # ******************************************************************
# 
# cws_org <- c("CWS-NOR","CWS-ONT","CWS-PRA","ECCC")
# xlim <- range(as.data.frame(st_coordinates(BBMP_boundary))$X)
# ylim <- range(as.data.frame(st_coordinates(BBMP_boundary))$Y)
# 
# survey_map <- ggplot() +
#   geom_sf(data = BCR_NA, col = "white", fill = "white")+
#   geom_sf(data = BCR_CAN, col = "gray98", fill = "white")+
#   geom_sf(data = BBMP_boundary, col = "transparent", fill = "forestgreen", alpha = 0.1)+
#   geom_sf(data = subset(bam_surveys, organization %!in% cws_org), size = 0.1, col = "gray70")+
#   geom_sf(data = subset(bam_surveys, organization %in% cws_org), size = 0.1, col = "black")+
#   geom_sf(data = BBMP_surveys, size = 0.1, col = "black")+
# 
#   annotation_scale(style = "ticks")+
# 
#   theme(panel.background = element_rect(fill = desaturate("#bbfdfd", amount = 0.5),
#                                         colour = desaturate("#bbfdfd", amount = 0.5),
#                                         size = 0.5, linetype = "solid"),
#         legend.spacing.y = unit(0, "mm"),
#         panel.border = element_rect(colour = "black", fill=NA),
#         legend.background = element_blank(),
#         legend.box.background = element_rect(colour = "black"),
#         legend.position=c(0.98,0.02),
#         legend.justification = c("right","bottom"),
#         legend.title = element_text(size = 10),
#         legend.text = element_text(size = 6))+
#   coord_sf(xlim = xlim, ylim = ylim)+
#   ggtitle("")
# 
# survey_map
# 
# png(paste0("../Output/Survey_Maps/BBMP_surveys.png"), width=8, height=8, units="in", res=300, type="cairo")
# print(survey_map)
# dev.off()
# 
# # ******************************************************************
# # ******************************************************************
# # FIT MODELS TO DATA, WITH AND WITHOUT BBMP DATA
# # ******************************************************************
# # ******************************************************************
# 
# # ---------------------------------------------------
# # Spatial covariates
# # ---------------------------------------------------
# 
# # Path to spatial covariates
# covar_folder <- "../../../Data/Spatial/"
# 
# # Mean Annual Temperature
# MAT <- rast(paste0(covar_folder,"National/Bioclimate/Normal_1991_2020_MAT.tif")) #%>%
# #project(target_raster, align = TRUE, method = "bilinear")%>%
# #resample(y = target_raster,method="bilinear")
# 
# # Mean Annual Precipitation
# MAP <- rast(paste0(covar_folder,"National/Bioclimate/Normal_1991_2020_MAP.tif")) #%>%
# #project(target_raster, align = TRUE, method = "bilinear")%>%
# #resample(y = target_raster,method="bilinear")
# 
# # Land cover of Canada 2020
# lcc2020 <- rast(paste0(covar_folder,"National/LandCoverCanada2020/landcover-2020-classification.tif"))
# 
# raster_stack <- list(MAT,MAP) %>% rast()
# names(raster_stack) = c("MAT","MAP")
# 
# # ---------------------------------------------------
# # Only select relevant columns
# # ---------------------------------------------------
# 
# bam_dat <- bam_surveys %>%dplyr::select(organization,project,lat,lon,date,distance,duration)
# 
# bbmp_dat <- BBMP_surveys %>%
#   mutate(organization = "CWS",project="BBMP",lat = Latitude, lon = Longitude, date = Date_Time, distance = Inf, duration = Survey_Duration_Minutes) %>%
#   dplyr::select(organization,project,lat,lon,date,distance,duration)
# 
# alldat <- bind_rows(bam_dat,bbmp_dat)
# 
# # ---------------------------------------------------
# # Extract covariate values at each survey location
# # ---------------------------------------------------
# 
# # Extract continuous covariates at each survey
# alldat = terra::extract(raster_stack,vect(alldat) , bind = TRUE) %>%
#   st_as_sf() %>%
#   st_transform(AEA_proj)
# 
# alldat_1km <- alldat %>% st_buffer(1000)
# 
# # Proportion of each land cover class from lcc 2020
# prop_LCC_1km <- exact_extract(lcc2020,alldat_1km,"frac") %>% suppressWarnings()
# names(prop_LCC_1km) <- paste0(str_replace(names(prop_LCC_1km),"frac","LCC"),"_1km")
# prop_LCC_1km[setdiff(paste0("LCC_",seq(1,18),"_1km"),names(prop_LCC_1km))] <- 0
# prop_LCC_1km %<>% dplyr::select(sort(names(.)))
# 
# # Combine land cover class names
# prop_LCC_1km <- prop_LCC_1km %>%
#   mutate(
#     Needleleaf_forest_1km = LCC_1_1km + LCC_2_1km,
#     Mixed_forest_1km = LCC_5_1km + LCC_6_1km,
#     Grass_shrub_1km = LCC_8_1km + LCC_10_1km,
#     Crop_1km = LCC_15_1km,
#     Urban_1km = LCC_17_1km,
#     Wetland_1km = LCC_14_1km,
#     Water_1km = LCC_18_1km) %>%
#   dplyr::select(Needleleaf_forest_1km:Water_1km)
# 
# alldat <- bind_cols(alldat,prop_LCC_1km)
# 
# alldat$Obs_Index <- 1:nrow(alldat)
# 
# allcounts <- bind_rows(as.data.frame(bam_counts),as.data.frame(BBMP_counts))
# 
# alldat <- na.omit(alldat)
# allcounts <- allcounts[alldat$Obs_Index,]
# 
# # ---------------------------------------------------
# # Conduct PCA
# # ---------------------------------------------------
# covars_for_PCA <- c("MAT",
#                     "MAP",
#                     "Needleleaf_forest_1km",
#                     "Mixed_forest_1km",
#                     "Grass_shrub_1km",
#                     "Crop_1km",
#                     "Urban_1km",
#                     "Wetland_1km",
#                     'Water_1km')
# 
# dat_for_PCA <- alldat %>%
#   as.data.frame() %>%
#   dplyr::select(covars_for_PCA)
# 
# pca <- prcomp(dat_for_PCA, scale = TRUE)
# 
# alldat_PCA <- predict(pca, newdata = as.data.frame(alldat)[,names(pca$center)])
# 
# # Scale/center covariates
# for (covar in colnames(alldat_PCA)){
# 
#   covar_mean <- mean(alldat_PCA[,covar],na.rm = TRUE)
#   covar_sd <- sd(alldat_PCA[,covar],na.rm = TRUE)
# 
#   alldat_PCA[,covar] <- (as.data.frame(alldat_PCA)[,covar] - covar_mean)/covar_sd
# 
# }
# 
# alldat <- alldat %>% bind_cols(alldat_PCA)
# 
# # Skree plot
# factoextra::fviz_eig(pca)
# 
# # Biplot
# factoextra::fviz_pca_var(pca,
#              axes = c(1,2),
#              col.var = "contrib", # Color by contributions to the PC
#              gradient.cols = viridis(10),
#              repel = TRUE     # Avoid text overlapping
# )
# 
# # ---------------------------------------------------
# # Create cross-validation folds
# # ---------------------------------------------------
# 
# # Create 100 x 100 km spatial folds
# set.seed(999)
# n_folds <- 5
# Crossval_Grid <- st_make_grid(
#   st_buffer(BBMP_boundary,50000),
#   cellsize = units::set_units(100*100,km^2),
#   what = "polygons",
#   square = TRUE,
#   flat_topped = FALSE) %>%
#   st_as_sf() %>%
#   na.omit() %>%
#   mutate(id = sample(1:nrow(.),nrow(.))) %>%
#   mutate(Crossval_Fold = cut(id,breaks = seq(0,max(id)+1,length.out = n_folds+1)) %>% as.factor() %>% as.numeric())%>%
#   dplyr::select(Crossval_Fold,x)
# 
# alldat$Obs_Index <- 1:nrow(alldat)
# alldat <- alldat %>% st_intersection(Crossval_Grid) %>% arrange(Obs_Index)
# 
# # ---------------------------------------------------
# # Save relevant data so the stuff above doesn't have to be re-run
# # ---------------------------------------------------
# 
# analysis_data_package <- list(alldat = alldat,
#                               allcounts = allcounts,
#                               BBMP_counts = BBMP_counts,
#                               BBMP_boundary = BBMP_boundary,
#                               Crossval_Grid = Crossval_Grid,
#                               AEA_proj = AEA_proj,
#                               pca = pca)
# 
# saveRDS(analysis_data_package,file = "C:/Users/IlesD/OneDrive - EC-EC/Iles/Projects/Landbirds/Landbird-Distribution-Modeling-ECCC/Analysis/BBMP/Output/Crossvalidation/crossval_data_prepared.RDS")

# # *************************************************************
# # *************************************************************
# # Fit models
# # *************************************************************
# # *************************************************************

analysis_data_package <- readRDS("C:/Users/IlesD/OneDrive - EC-EC/Iles/Projects/Landbirds/Landbird-Distribution-Modeling-ECCC/Analysis/BBMP/Output/Crossvalidation/crossval_data_prepared.RDS")
attach(analysis_data_package)

`%!in%` <- Negate(`%in%`)

# ---------------------------------------------------
# Create mesh for INLA
# ---------------------------------------------------

napops_species <- list_species() %>% rename(Species_Code_NAPOPS = Species,
                                            Common_Name_NAPOPS = Common_Name,
                                            Scientific_Name_NAPOPS = Scientific_Name)

# ---------------------------------------------------
# Run through species and fit models
# ---------------------------------------------------

# Timeout INLA after 10 minutes (if it has not fit by then, it has likely stalled)
inla.setOption(inla.timeout = 60*20)

# Relative abundances of each species
species_relabund <- allcounts 
species_relabund[species_relabund>0] <- 1
species_relabund <- apply(species_relabund,2,sum) %>% sort(decreasing = TRUE)

results <- data.frame()
results_path <- "../Output/Crossvalidation/Crossval_results_with_without_BBMP_BRT_samplesize.rds"

for (xval_fold in (1:5)){
  for (sp_code in (names(species_relabund)[1:20])){
    
    if (file.exists(results_path)){
      results <- readRDS(results_path)
      if (nrow(results[which(results$sp_code == sp_code & results$Crossval_Fold == xval_fold),])>0) next
    }
    
    print(paste0(sp_code," fold ", xval_fold))
    if (nrow(results)>0){
      if (nrow(subset(results,Species == sp_code & fold == xval_fold))>0) next
    }
    
    # Prepare data for this species
    sp_dat <- alldat %>%
      mutate(count = allcounts[,sp_code]) %>%
      mutate(presence = as.numeric(count>0)) %>%
      subset(duration >=1 & duration <= 10)
    
    # ----------------------------------------------------
    # Generate QPAD offsets for each survey (assumes unlimited distance point counts)
    # ----------------------------------------------------
    sp_napops <- subset(napops_species,Species_Code_NAPOPS == sp_code)
    sp_dat$log_offset <- 0
    if (nrow(sp_napops)>0){
      if (sp_napops$Removal == 1 & sp_napops$Distance == 1){
        
        cue_rate <- cue_rate(species = sp_code,od = 153, tssr = 0, model = 1)[3] %>% as.numeric()
        EDR <- edr(species = sp_code,road = FALSE, forest = 0.5,model = 1)[3] %>% as.numeric()
        
        # Calculate A and p, which jointly determine offset
        A_metres <- c(pi*EDR^2)
        p <- 1-exp(-sp_dat$duration*cue_rate)
        
        sp_dat$log_offset <- log(A_metres * p)
        
      }
    }
    
    # ----------------------------------------------------
    # Separate Data types
    # ----------------------------------------------------
    
    
    sp_dat$Date <- lubridate::ymd_hms(sp_dat$date)
    sp_dat <- na.omit(sp_dat)
    sp_dat <- subset(sp_dat,
                     yday(Date) >= yday(ymd("2022-06-01")) &
                     yday(Date) <= yday(ymd("2022-07-15")) &
                     hour(Date) >= 5 &
                     hour(Date) <= 12)
    
    validation_data <- subset(sp_dat, Crossval_Fold == xval_fold & (project == "BBMP" | organization %in% c("CWS-NOR","CWS-ONT","CWS-PRA","ECCC"))) %>% as.data.frame()
    training_data <- subset(sp_dat, Crossval_Fold != xval_fold) %>% as.data.frame()
    rm(sp_dat)
    
    # ----------------------------------------------------
    # Create two datasets:
    # 1) 100% of data is from non-BBMP
    # 2) 50% of data is from non-BBMP, and 50% is from BBMP
    # Ensure sample size is the same
    # ----------------------------------------------------
    training_data$Obs_Index <- 1:nrow(training_data)
    training_data_bam <-  subset(training_data, project != "BBMP" & organization %!in% c("CWS-NOR","CWS-ONT","CWS-PRA","ECCC"))
    sample_size_bam <- nrow(training_data_bam)
    
    training_data_bbmp <-  subset(training_data, project == "BBMP" | organization %in% c("CWS-NOR","CWS-ONT","CWS-PRA","ECCC"))
    
    # create dataset with 50% bbmp and 50% non-bbmp data
    training_data_5050 <- sample_n(training_data_bam,size=round(sample_size_bam/2)) %>%
      bind_rows(sample_n(training_data_bbmp, size = round(sample_size_bam/2)))
    
    covariates_to_include <- training_data %>% dplyr::select(duration:PC9) %>% colnames()
    
    # par(mfrow=c(1,2))
    # plot(x = training_data$lon, y = training_data$lat, pch = 19, cex = 0.1)
    # points(x = training_data_bam$lon, y = training_data_bam$lat, pch = 19, cex = 0.1, col = "red")
    # points(x = validation_data$lon, y = validation_data$lat, pch = 19, cex = 0.2, col = "green")
    # 
    # plot(x = training_data$lon, y = training_data$lat, pch = 19, cex = 0.1)
    # points(x = training_data_bbmp$lon, y = training_data_bbmp$lat, pch = 19, cex = 0.1, col = "blue")
    # points(x = validation_data$lon, y = validation_data$lat, pch = 19, cex = 0.2, col = "green")
    # 
    # ----------------------------------------------------------------------------
    # Fit model to only non-bbmp data, assess cross-validation accuracy
    # ----------------------------------------------------------------------------
    start <- Sys.time()
    
    fit_bam <-  fit_brt(model_data = training_data_bam,
                        response_column = which(colnames(training_data_bam)=="count"),
                        covariate_columns = which(colnames(training_data_bam)%in% covariates_to_include))
    
    pred_bam  <- predict(fit_bam, 
                         validation_data, 
                         n.trees = fit_bam$gbm.call$best.trees)
    
    # ----------------------------------------------------------------------------
    # Fit model to all data (including bbmp)
    # ----------------------------------------------------------------------------
    
    fit_bbmp <-  fit_brt(model_data = training_data_5050,
                         response_column = which(colnames(training_data_5050)=="count"),
                         covariate_columns = which(colnames(training_data_5050)%in% covariates_to_include))
    
    pred_bbmp  <- predict(fit_bbmp, 
                          validation_data, 
                          n.trees = fit_bbmp$gbm.call$best.trees)
    
    end <- Sys.time()
    
    # ----------------------------------------------------------------------------
    # Evaluate cross-validation accuracy
    # ----------------------------------------------------------------------------
    
    pred_bam <- exp(pred_bam + validation_data$log_offset)
    pred_bbmp <- exp(pred_bbmp + validation_data$log_offset)
    
    prob_nonzero_bam <- 1-exp(-pred_bam)
    prob_nonzero_bbmp <- 1-exp(-pred_bbmp)
    
    AUC_bam <- auc(response = validation_data$presence,predictor = prob_nonzero_bam)
    AUC_bbmp <- auc(response = validation_data$presence,predictor = prob_nonzero_bbmp)
    
    lppd_bam <- sum(dpois(validation_data$count,lambda=pred_bam,log = TRUE))
    lppd_bbmp <- sum(dpois(validation_data$count,lambda=pred_bbmp,log = TRUE))
    
    
    sp_results <- data.frame(Species = sp_code,
                             fold = xval_fold,
                             mean_val = mean(validation_data$count),
                             
                             # Ability to predict presences and absences
                             AUC_bam = mean(AUC_bam),
                             AUC_bbmp = mean(AUC_bbmp),
                             
                             # Ability to predict counts
                             lppd_bam = mean(lppd_bam),
                             lppd_bbmp = mean(lppd_bbmp),
                             
                             # Actual predictions of mean counts in validation data
                             # could be used for precision, accuracy, and coverage assessments
                             pred_bam_mean = mean(pred_bam),
                             pred_bbmp_mean = mean(pred_bbmp))
    
    if (file.exists(results_path)) results <- readRDS(results_path)
    
    results <- rbind(results, sp_results)
    
    saveRDS(results,results_path)
    
    
    # ----------------------------------------------------------------------------
    # Plot results
    # ----------------------------------------------------------------------------
    if (file.exists(results_path)) results <- readRDS(results_path)
    
    results_summary <- results %>% 
      group_by(Species) %>% 
      summarize_all(mean) %>%
      arrange(mean_val)
    
    results_summary <- bind_rows(data.frame(Species = ""),results_summary)
    
    results_summary$Species <- factor(results_summary$Species, levels = results_summary$Species)
    
    label_df <- data.frame(Species = "", x = c(0.6,0.75,0.9), text = c("Poor","Adequate","Excellent"))
    label_df$text <- factor(label_df$text, levels = label_df$text)
    
    colors <- RColorBrewer::brewer.pal(3,name = "RdYlBu")
   
    comparison_plot2 <- ggplot(data = results_summary)+
      
      geom_rect(aes(xmin=0, xmax=0.7, ymin=-Inf, ymax=Inf), fill = colors[1], alpha = 0.05)+
      geom_rect(aes(xmin=0.7, xmax=0.8, ymin=-Inf, ymax=Inf), fill = colors[2], alpha = 0.05)+
      geom_rect(aes(xmin=0.8, xmax=1, ymin=-Inf, ymax=Inf), fill = colors[3], alpha = 0.05)+
      geom_text(data = label_df, aes(x = x, y = Species, label = text), fontface = "bold", alpha = 0.5)+
      
      scale_fill_manual(values = colors)+ 
      geom_point(data = results_summary,aes(y = Species, x = AUC_bam),size = 2,col = "black")+
      
      geom_segment(data = results_summary,aes(y = Species, yend = Species, x = AUC_bam, xend = AUC_bbmp),
                   size = 2,
                   arrow = arrow(length = unit(0.2, "cm")),
                   col = "black")+
      
      coord_cartesian(xlim = c(0.5,1))+
      ggtitle("Predictive performance after including BBMP data\n\n(Measured at continental scale)")+
      xlab("AUC score\n\n(cross-validation)")
    
    print(comparison_plot2)
    
  }
}


if (file.exists(results_path)) results <- readRDS(results_path)

results_summary <- results %>% 
  group_by(Species) %>% 
  summarize_all(mean) %>%
  arrange(mean_val)

results_summary <- bind_rows(data.frame(Species = ""),results_summary)

results_summary$Species <- factor(results_summary$Species, levels = results_summary$Species)

label_df <- data.frame(Species = "", x = c(0.6,0.75,0.9), text = c("Poor","Adequate","Excellent"))
label_df$text <- factor(label_df$text, levels = label_df$text)

colors <- RColorBrewer::brewer.pal(3,name = "RdYlBu")

comparison_plot2 <- ggplot(data = results_summary)+
  
  geom_rect(aes(xmin=0, xmax=0.7, ymin=-Inf, ymax=Inf), fill = colors[1], alpha = 0.05)+
  geom_rect(aes(xmin=0.7, xmax=0.8, ymin=-Inf, ymax=Inf), fill = colors[2], alpha = 0.05)+
  geom_rect(aes(xmin=0.8, xmax=1, ymin=-Inf, ymax=Inf), fill = colors[3], alpha = 0.05)+
  geom_text(data = label_df, aes(x = x, y = Species, label = text), fontface = "bold", alpha = 0.5)+
  
  scale_fill_manual(values = colors)+ 
  geom_point(aes(y = Species, x = AUC_bam),size = 2,col = "black")+
  
  geom_segment(aes(y = Species, yend = Species, x = AUC_bam, xend = AUC_bbmp),
               size = 2,
               arrow = arrow(length = unit(0.2, "cm")),
               col = "black")+
  
  coord_cartesian(xlim = c(0.5,1))+
  ggtitle("Predictive performance with BBMP data")+
  xlab("AUC score\n\n(cross-validation)")

print(comparison_plot2)
