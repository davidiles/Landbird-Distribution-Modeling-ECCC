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
             "napops")

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
# # cws_org <- c("CWS-NOR","CWS-ONT","CWS-PRA","ECCC")
# # xlim <- range(as.data.frame(st_coordinates(BBMP_boundary))$X)
# # ylim <- range(as.data.frame(st_coordinates(BBMP_boundary))$Y)
# # 
# # survey_map <- ggplot() +
# #   geom_sf(data = BCR_NA, col = "white", fill = "white")+
# #   geom_sf(data = BCR_CAN, col = "gray98", fill = "white")+
# #   geom_sf(data = BBMP_boundary, col = "transparent", fill = "forestgreen", alpha = 0.1)+
# #   geom_sf(data = subset(bam_surveys, organization %!in% cws_org), size = 0.1, col = "gray70")+
# #   geom_sf(data = subset(bam_surveys, organization %in% cws_org), size = 0.1, col = "black")+
# #   geom_sf(data = BBMP_surveys, size = 0.1, col = "black")+
# # 
# #   annotation_scale(style = "ticks")+
# # 
# #   theme(panel.background = element_rect(fill = desaturate("#bbfdfd", amount = 0.5),
# #                                         colour = desaturate("#bbfdfd", amount = 0.5),
# #                                         size = 0.5, linetype = "solid"),
# #         legend.spacing.y = unit(0, "mm"),
# #         panel.border = element_rect(colour = "black", fill=NA),
# #         legend.background = element_blank(),
# #         legend.box.background = element_rect(colour = "black"),
# #         legend.position=c(0.98,0.02),
# #         legend.justification = c("right","bottom"),
# #         legend.title = element_text(size = 10),
# #         legend.text = element_text(size = 6))+
# #   coord_sf(xlim = xlim, ylim = ylim)+
# #   ggtitle("")
# # 
# # survey_map
# # 
# # png(paste0("../Output/Survey_Maps/BBMP_surveys.png"), width=8, height=8, units="in", res=300, type="cairo")
# # print(survey_map)
# # dev.off()
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
# 
# # Mean Annual Precipitation
# MAP <- rast(paste0(covar_folder,"National/Bioclimate/Normal_1991_2020_MAP.tif")) #%>%
# 
# # Land cover of Canada 2020
# lcc2020 <- rast(paste0(covar_folder,"National/LandCoverCanada2020/landcover-2020-classification.tif"))
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
# # Extract covariate values within 1 km of each survey location
# # ---------------------------------------------------
# 
# alldat_1km <- alldat %>% st_buffer(1000)
# 
# # Extract continuous covariates at each survey location
# MAP_1km <- exact_extract(MAP,alldat_1km,"mean") %>% suppressWarnings()
# MAT_1km <- exact_extract(MAT,alldat_1km,"mean") %>% suppressWarnings()
# 
# alldat$MAP_1km <- MAP_1km
# alldat$MAT_1km <- MAT_1km
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
# # Extract covariate values within 25 km of each survey location
# # ---------------------------------------------------
# 
# alldat_25km <- alldat %>% st_buffer(25000)
# 
# # Extract continuous covariates at each survey location
# MAP_25km <- exact_extract(MAP,alldat_25km,"mean") %>% suppressWarnings()
# MAT_25km <- exact_extract(MAT,alldat_25km,"mean") %>% suppressWarnings()
# 
# alldat$MAP_25km <- MAP_25km
# alldat$MAT_25km <- MAT_25km
# 
# # Proportion of each land cover class from lcc 2020
# prop_LCC_25km <- exact_extract(lcc2020,alldat_25km,"frac") %>% suppressWarnings()
# names(prop_LCC_25km) <- paste0(str_replace(names(prop_LCC_25km),"frac","LCC"),"_25km")
# prop_LCC_25km[setdiff(paste0("LCC_",seq(1,18),"_25km"),names(prop_LCC_25km))] <- 0
# prop_LCC_25km %<>% dplyr::select(sort(names(.)))
# 
# # Combine land cover class names
# prop_LCC_25km <- prop_LCC_25km %>%
#   mutate(
#     Needleleaf_forest_25km = LCC_1_25km + LCC_2_25km,
#     Mixed_forest_25km = LCC_5_25km + LCC_6_25km,
#     Grass_shrub_25km = LCC_8_25km + LCC_10_25km,
#     Crop_25km = LCC_15_25km,
#     Urban_25km = LCC_17_25km,
#     Wetland_25km = LCC_14_25km,
#     Water_25km = LCC_18_25km) %>%
#   dplyr::select(Needleleaf_forest_25km:Water_25km)
# 
# alldat <- bind_cols(alldat,prop_LCC_25km)
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
# covars_for_PCA <- c("MAT_1km",
#                     "MAP_1km",
#                     "Needleleaf_forest_1km",
#                     "Mixed_forest_1km",
#                     "Grass_shrub_1km",
#                     "Crop_1km",
#                     "Urban_1km",
#                     "Wetland_1km",
#                     'Water_1km',
#                     "MAT_25km",
#                     "MAP_25km",
#                     "Needleleaf_forest_25km",
#                     "Mixed_forest_25km",
#                     "Grass_shrub_25km",
#                     "Crop_25km",
#                     "Urban_25km",
#                     "Wetland_25km",
#                     'Water_25km')
# 
# dat_for_PCA <- alldat %>%
#   as.data.frame() %>%
#   dplyr::select(covars_for_PCA)
# 
# pca <- prcomp(dat_for_PCA, scale = TRUE)
# 
# # Skree plot
# factoextra::fviz_eig(pca)
# 
# # Biplot
# factoextra::fviz_pca_var(pca,
#                          axes = c(2,3),
#                          col.var = "contrib", # Color by contributions to the PC
#                          gradient.cols = viridis(10),
#                          repel = TRUE     # Avoid text overlapping
# )
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
# species with qpad offsets
# ---------------------------------------------------

napops_species <- list_species() %>% rename(Species_Code_NAPOPS = Species,
                                            Common_Name_NAPOPS = Common_Name,
                                            Scientific_Name_NAPOPS = Scientific_Name)

# ---------------------------------------------------
# Create mesh for INLA
# ---------------------------------------------------

require(INLA)
require(inlabru)

proj <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
alldat <- alldat %>% st_transform(proj)
BBMP_boundary <- st_transform(BBMP_boundary, proj)

pt.bond <- inla.nonconvex.hull(coordinates(as(alldat,'Spatial')), 4, 2, resolution = c(64,21),crs=proj)
mesh <- inla.mesh.2d(loc=coordinates(as(alldat,'Spatial')),
                     boundary=pt.bond, 
                     max.edge=c(1,3), 
                     cut=0.5,
                     off=c(1e-5,4), crs=proj)
plot(mesh)
mesh_locs <- mesh$loc[,c(1,2)] %>% as.data.frame()
dim(mesh_locs)

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
results_path <- "../Output/Crossvalidation/Crossval_results_with_without_BBMP_INLA_fixedsamplesize.rds"
for (sp_code in (names(species_relabund)[1:100])){
  for (xval_fold in (1:5)){
    
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
    sp_dat <-  sp_dat %>% st_as_sf(coords = c("lon", "lat"),crs=4326, remove = FALSE)
    
    validation_data <- subset(sp_dat, Crossval_Fold == xval_fold)
    training_data <- subset(sp_dat, Crossval_Fold != xval_fold)
    
    # ----------------------------------------------------
    # Create two datasets:
    # 1) 100% of data is from non-BBMP
    # 2) 50% of data is from non-BBMP, and 50% is from BBMP
    # Ensure sample size is the same
    # ----------------------------------------------------
    
    training_data$Obs_Index <- 1:nrow(training_data)
    training_data_bam <-  subset(training_data, project != "BBMP" & organization %!in% c("CWS-NOR","CWS-ONT","CWS-PRA","ECCC"))
    training_data_bam$data_type <- "bam"
    sample_size_bam <- nrow(training_data_bam)
    
    training_data_bbmp <-  subset(training_data, project == "BBMP" | organization %in% c("CWS-NOR","CWS-ONT","CWS-PRA","ECCC"))
    training_data_bbmp$data_type <- "bbmp"
    
    # create dataset with 50% bbmp and 50% non-bbmp data
    training_data_5050 <- sample_n(training_data_bam,size=round(sample_size_bam/2)) %>%
      bind_rows(sample_n(training_data_bbmp, size = round(sample_size_bam/2)))
    covariates_to_include <- training_data %>% dplyr::select(duration:PC9) %>% colnames()
    
    
    validation_data$data_type <- "bam"
    validation_data$data_type[validation_data$project == "BBMP" | validation_data$organization %in% c("CWS-NOR","CWS-ONT","CWS-PRA","ECCC")] <- "bbmp"
    validation_data <- validation_data %>% as('Spatial')
    training_data_bam <- training_data_bam %>% as('Spatial')
    training_data_5050 <- training_data_5050 %>% as('Spatial')
    
    # ----------------------------------------------------
    # Prepare formulas and priors
    # ----------------------------------------------------
    
    # Controls the 'residual spatial field'.  This can be adjusted to create smoother surfaces.
    prior_range <- c(5,0.1)
    prior_sigma <- c(1,0.5)
    
    matern_coarse <- inla.spde2.pcmatern(mesh,
                                         prior.range = prior_range,
                                         prior.sigma = prior_sigma,
                                         constr = TRUE
    )
    
    # Mesh for survey duration effect
    duration_meshpoints <- seq(min(sp_dat$duration)-0.1,max(sp_dat$duration)+0.1,length.out = 11)
    duration_mesh1D = inla.mesh.1d(duration_meshpoints,boundary="free")
    duration_spde = inla.spde2.pcmatern(duration_mesh1D,
                                        prior.range = c(5,0.1),
                                        prior.sigma = c(1,0.5))
    
    
    covariates_to_include <- paste0("PC",1:11)
    
    # How much shrinkage should be applied to covariate effects?
    sd_linear <- 1   # Change to smaller value (e.g., 0.1), if you want to heavily shrink covariate effects and potentially create smoother surfaces
    prec_linear <-  c(1/sd_linear^2,1/(sd_linear/2)^2)
    
    model_components = as.formula(paste0('~
            Intercept_BBMP(1)+
            Intercept_BAM(1)+
            spde_coarse(main = coordinates, model = matern_coarse) +',
                                         paste0("Beta1_",covariates_to_include,'(1,model="linear", mean.linear = 0, prec.linear = ', prec_linear[1],')', collapse = " + ")))
    
    model_formula_BBMP = as.formula(paste0('count ~
                  Intercept_BBMP +
                  log_offset +
                  spde_coarse +',
                                         paste0("Beta1_",covariates_to_include,'*',covariates_to_include, collapse = " + ")))
    
    model_formula_BAM = as.formula(paste0('count ~
                  Intercept_BAM +
                  log_offset +
                  spde_coarse +',
                                         paste0("Beta1_",covariates_to_include,'*',covariates_to_include, collapse = " + ")))

    # ----------------------------------------------------------------------------
    # Fit model to only non-bbmp data, assess cross-validation accuracy
    # ----------------------------------------------------------------------------
    start <- Sys.time()
    fit_bam <- NULL
    while(is.null(fit_bam)){
      
      fit_model <- function(){
        tryCatch(expr = {bru(components = model_components,
                             
                             like(family = "nbinomial",
                                  formula = model_formula_BAM,
                                  data = training_data_bam),
                             options = list(control.compute = list(waic = FALSE, cpo = FALSE),
                                            bru_verbose = 4))},
                 error = function(e){NULL})
      }
      fit_bam <- fit_model()
      
      if ("try-error" %in% class(fit_bam)) fit_bam <- NULL
    }
    end <- Sys.time()
    runtime_bam <- difftime( end,start, units="mins") %>% round(2)
    print(paste0(sp_code," - ",runtime_bam," min to fit model")) 
    
    # ----------------------------------------------------------------------------
    # Fit model to all data (including bbmp)
    # ----------------------------------------------------------------------------
    
    start <- Sys.time()
    fit_bbmp <- NULL
    while(is.null(fit_bbmp)){
      
      fit_model <- function(){
        tryCatch(expr = {bru(components = model_components,
                             
                             like(family = "nbinomial",
                                  formula = model_formula_BAM,
                                  data = subset(training_data_5050,data_type == "bam")),
                             
                             like(family = "nbinomial",
                                  formula = model_formula_BBMP,
                                  data = subset(training_data_5050,data_type == "bbmp")),
                             
                             options = list(control.compute = list(waic = FALSE, cpo = FALSE),
                                            bru_verbose = 4))},
                 error = function(e){NULL})
      }
      fit_bbmp <- fit_model()
      
      if ("try-error" %in% class(fit_bbmp)) fit_bbmp <- NULL
    }
    end <- Sys.time()
    runtime_bbmp <- difftime( end,start, units="mins") %>% round(2)
    print(paste0(sp_code," - ",runtime_bbmp," min to fit model")) 
    
    # ------------------------------------------------
    # Predictions for model fit to non-bbmp data
    # ------------------------------------------------
    pred_formula = as.formula(paste0(' ~Intercept_BAM +log_offset +spde_coarse +',paste0("Beta1_",covariates_to_include,'*',covariates_to_include, collapse = " + ")))
  
    size <- fit_bam$summary.hyperpar$'0.5quant'[1]
    pred_bam <- NULL
    pred_bam <- generate(fit_bam,
                         validation_data,
                         formula =  pred_formula,
                         n.samples = 250)
    
    pred_bam <- exp(pred_bam)
    mean_bam <- apply(pred_bam,2,mean)
    
    pred_bam <- apply(pred_bam,1,median)
    prob_zero <- dnbinom(0,mu=pred_bam ,size=size)
    prob_nonzero_bam <- 1-prob_zero
    AUC_bam <- auc(response = validation_data$presence,predictor = prob_nonzero_bam)
    lppd_bam <- sum(dnbinom(validation_data$count,mu=pred_bam,size=size,log = TRUE))
    
    # ------------------------------------------------
    # Predictions for model fit to non-bbmp + bbmp data
    # ------------------------------------------------
    
    # Predictions for bam data
    pred_formula_1 = as.formula(paste0(' ~Intercept_BAM +log_offset +spde_coarse +',paste0("Beta1_",covariates_to_include,'*',covariates_to_include, collapse = " + ")))
    size_1 <- fit_bbmp$summary.hyperpar$'0.5quant'[1]
    pred_bbmp_1 <- NULL
    pred_bbmp_1 <- generate(fit_bbmp,
                         subset(validation_data, data_type == "bam"),
                         formula =  pred_formula_1,
                         n.samples = 250)
    
    pred_bbmp_1 <- exp(pred_bbmp_1)
    pred_bbmp_1 <- apply(pred_bbmp_1,1,median)
    prob_zero_1 <- dnbinom(0,mu=pred_bbmp_1 ,size=size_1)
    prob_nonzero_bbmp_1 <- 1-prob_zero_1
    
    # predictions for bbmp data
    pred_formula_2 = as.formula(paste0(' ~Intercept_BBMP +log_offset +spde_coarse +',paste0("Beta1_",covariates_to_include,'*',covariates_to_include, collapse = " + ")))
    size_2 <- fit_bbmp$summary.hyperpar$'0.5quant'[2]
    pred_bbmp_2 <- NULL
    pred_bbmp_2 <- generate(fit_bbmp,
                            subset(validation_data, data_type == "bbmp"),
                            formula =  pred_formula_2,
                            n.samples = 250)
    
    pred_bbmp_2 <- exp(pred_bbmp_2)
    pred_bbmp_2 <- apply(pred_bbmp_2,1,median)
    prob_zero_2 <- dnbinom(0,mu=pred_bbmp_2 ,size=size_2)
    prob_nonzero_bbmp_2 <- 1-prob_zero_2
    
    AUC_bbmp <- auc(response = c(subset(validation_data, data_type == "bam")$presence,subset(validation_data, data_type == "bbmp")$presence),predictor = c(prob_nonzero_bbmp_1,prob_nonzero_bbmp_2))
    lppd_bbmp <- sum(dnbinom(subset(validation_data, data_type == "bam")$count,mu=pred_bbmp_1,size=size_1,log = TRUE)) + sum(dnbinom(subset(validation_data, data_type == "bbmp")$count,mu=pred_bbmp_2,size=size_2,log = TRUE))
    
    sp_results <- data.frame(Species = sp_code,
                             fold = xval_fold,
                             mean_val = mean(validation_data$count),
                             
                             # Ability to predict presences and absences
                             AUC_bam = AUC_bam,
                             AUC_bbmp = AUC_bbmp,
                             
                             # Ability to predict counts
                             lppd_bam = lppd_bam,
                             lppd_bbmp = lppd_bbmp)
                             
    
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
      
      geom_rect(aes(xmin=0, xmax=0.7, ymin=-Inf, ymax=Inf), fill = colors[1], alpha = 0.3)+
      geom_rect(aes(xmin=0.7, xmax=0.8, ymin=-Inf, ymax=Inf), fill = colors[2], alpha = 0.3)+
      geom_rect(aes(xmin=0.8, xmax=1, ymin=-Inf, ymax=Inf), fill = colors[3], alpha = 0.3)+
      geom_text(data = label_df, aes(x = x, y = Species, label = text), fontface = "bold", alpha = 0.5)+
      
      scale_fill_manual(values = colors)+ 
      geom_point(data = results_summary,aes(y = Species, x = AUC_bam),size = 2,col = "black")+
      
      geom_segment(data = results_summary,aes(y = Species, yend = Species, x = AUC_bam, xend = AUC_bbmp),
                   size = 2,
                   arrow = arrow(length = unit(0.2, "cm")),
                   col = "black")+
      
      coord_cartesian(xlim = c(0.5,1))+
      ggtitle("Predictive performance after including BBMP data\n\n(Measured at national scale)")+
      xlab("AUC score\n\n(cross-validation)")
    
    print(comparison_plot2)
    
  }
}
