# ************************************************
# BAYESIAN ANALYSIS / SPECIES DISTRIBUTION MODELS FOR SAONATCHEWAN BREEDING BIRD ATLAS
# 
# 1) A 'processed data package' for analysis is prepared by previous script
# ************************************************

rm(list=ls())

# ------------------------------------------------
# Set working directory
# ------------------------------------------------

# dirname <- thisPath()
dirname <- "C:/Users/IlesD/OneDrive - EC-EC/Iles/Projects/Landbirds/Landbird-Distribution-Modeling-ECCC/Analysis/Ontario_BBA/Code"
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

analysis_data <- readRDS("../Data_Cleaned/analysis_data_package.rds")
attach(analysis_data)

# Which surveys are in the final dataset?
survey_summary <- all_surveys %>%
  as.data.frame() %>%
  group_by(Survey_Type,Project_Name) %>%
  summarize(n = n()) %>%
  pivot_wider(names_from = Survey_Type, values_from = n) %>%
  mutate(Total = rowSums(across(where(is.numeric)), na.rm = TRUE)) %>%
  relocate(Project_Name,Total,Point_Count,ARU_SPT) %>%
  arrange(desc(Total))

AEA_proj <- "+proj=aea +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-87 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs "

ONBoundary <- st_read("../../../Data/Spatial/National/BCR/BCR_Terrestrial_master.shp")  %>%
  subset(PROVINCE_S == "ONTARIO") %>%
  st_make_valid() %>%
  st_union() %>%
  st_transform(AEA_proj)

# ONGrid <- ONGrid %>% st_intersection(ONBoundary)
ONGrid <- ONGrid %>% rename(geometry = x)

# Raster with target properties (needed for plotting)
target_raster <- rast("../../../Data/Spatial/National/AnnualMeanTemperature/wc2.1_30s_bio_1.tif") %>% 
  crop(st_transform(ONBoundary,crs(.))) %>% 
  project(AEA_proj, res = 1000) %>%
  mask(vect(ONBoundary))

# ------------------------------------------
# Select Point Counts / ARUs to use
# Decision to remove "SPM" aru transcriptions, as they are considered less reliable
# ------------------------------------------

all_surveys$Hours_Since_Sunrise <- as.numeric(all_surveys$Hours_Since_Sunrise)
all_surveys$Obs_Index <- 1:nrow(all_surveys)

PC_to_use <- subset(all_surveys,
                    Survey_Type %in% c("Point_Count","ARU_SPT") &
                      
                      Survey_Duration_Minutes >= 1 &
                      Survey_Duration_Minutes <= 10 &
                      
                      Hours_Since_Sunrise >= -2 &
                      Hours_Since_Sunrise <= 4 &
                      
                      yday(Date_Time) >= yday(ymd("2022-05-15")) &
                      yday(Date_Time) <= yday(ymd("2022-07-15"))
)

dim(PC_to_use) # 124291

# ------------------------------------------
# Calculate mean covariate values within each Atlas square (covariates taken from ONGrid)
# Note: these will be merged with new stationary checklist dataframe ("ONSquareChecklist") within the species loop
# ------------------------------------------

# This vector should contain names of each covariate needed for analysis
covars <- paste0("PC",1:9)
cov_rast <- stars::st_rasterize(ONGrid %>% dplyr::select(covars, geometry)) %>% rast()
#plot(cov_rast)

# ------------------------------------------
# Subset
# ------------------------------------------

surveys_to_use <- c(PC_to_use$Obs_Index) # SC_to_use$Obs_Index, # LT_to_use$Obs_Index
all_surveys <- subset(all_surveys, Obs_Index %in% surveys_to_use) %>% arrange(Obs_Index)
full_count_matrix <- full_count_matrix[all_surveys$Obs_Index,]

# ******************************************************************
# Loop through species, fit models, generate maps
# ******************************************************************
#species_summary <- species_summary %>% subset(`PC/ARU` >=20)

species_list <- subset(species_to_model,
                       english_name %in% 
                         c("Chipping Sparrow")) %>%
  unique()

#species_list <- species_to_model[1:20,]

for (i in 1:nrow(species_list)){
  
  # ----------------------------------------------------
  # Extract counts/data for this species
  # ----------------------------------------------------
  
  sp_code = as.character(species_list$Species_Code_BSC[i])
  species_name = species_list$english_name[i]
  print(species_name)
  
  # Prepare data for this species
  sp_dat <- all_surveys %>% 
    mutate(count = full_count_matrix[,sp_code],
           presence = as.numeric(full_count_matrix[,sp_code]>0))
  
  # ----------------------------------------------------
  # Extract ebird range for this species (if it exists); prepared by previous script
  # ----------------------------------------------------
  range <- NA
  if (species_name %in% names(species_ranges)){
    
    range <- species_ranges[[species_name]] %>% st_transform(st_crs(sp_dat))
    
    # Identify distance of each survey to the edge of species range (in km)
    sp_dat$distance_from_range <- ((st_distance(sp_dat, range)) %>% as.numeric())/1000
  } else{
    sp_dat$distance_from_range <- 0
  }
  
  # ----------------------------------------------------
  # Generate QPAD offsets for each survey (assumes unlimited distance point counts)
  # ----------------------------------------------------
  
  # sp_dat$log_QPAD_offset
  
  # ----------------------------------------------------
  # Separate Data types
  # ----------------------------------------------------
  
  # Point count and ARU data treated as the same thing in this analysis
  PC_dat <- subset(sp_dat,Survey_Type %in% c("Point_Count","ARU_SPT")) %>% as('Spatial')
  
  # ----------------------------------------------------
  # Create a spatial mesh, which is used to fit the residual spatial field
  # NOTE: MAKING A COARSER MESH (LARGER VALUES OF MAX_EDGE) WILL SPEED UP COMPUTATIONS AND REDUCE RAM REQUIREMENTS, 
  # BUT WITH POTENTIAL LOSS OF ACCURACY/PRECISION
  # DEFAULT SHOULD CREATE APPROX 5000
  # ----------------------------------------------------
  
  # make a two extension hulls and mesh for spatial model
  hull <- fm_extensions(
    ONBoundary,
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
  
  # Controls the 'residual spatial field'.  This can be adjusted to create smoother surfaces.
  matern_abund <- inla.spde2.pcmatern(mesh_spatial,
                                       prior.range = c(300000, 0.1), # 10% chance range is smaller than 300000
                                       prior.sigma = c(0.5,0.1),    # 10% chance sd is larger than 0.5,
                                       constr = TRUE
  )
  
  # For every pixel on landscape, extract distance (in km) from eBird range limit
  ONGrid_species <- ONGrid %>%
    mutate(distance_from_range = (st_centroid(.) %>% 
                                    st_distance( . , range) %>% 
                                    as.numeric())/1000)
  
  covariates_to_include <- paste0("PC",1:2) 
  
  # ***************************************************************
  # OBBA - 2
  # ***************************************************************
  
  # How much shrinkage should be applied to covariate effects?
  sd_linear <- 0.1  # Change to smaller value (e.g., 0.1), if you want to heavily shrink covariate effects and potentially create smoother surfaces
  prec_linear <-  c(1/sd_linear^2,1/(sd_linear/2)^2)
  
  model_components = as.formula(paste0('~
            Intercept_OBBA2(1)+
            range_effect_OBBA2(1,model="linear", mean.linear = -0.046, prec.linear = 10000)+
            spde_OBBA2(main = coordinates, model = matern_abund) +',
                                       paste0("Beta1_",covariates_to_include,'(1,model="linear", mean.linear = 0, prec.linear = ', prec_linear[1],')', collapse = " + "),
                                       '+',
                                       paste0("Beta2_",covariates_to_include,'(1,model="linear", mean.linear = 0, prec.linear = ', prec_linear[2],')', collapse = " + ")))
  
  model_formula_OBBA2 = as.formula(paste0('count ~
                  Intercept_OBBA2 +
                  range_effect_OBBA2 * distance_from_range +
                  spde_OBBA2 +',
                                       paste0("Beta1_",covariates_to_include,'*',covariates_to_include, collapse = " + "),
                                       '+',
                                       paste0("Beta2_",covariates_to_include,'*',covariates_to_include,"^2", collapse = " + ")))
  
  start <- Sys.time()
  fit_INLA <- NULL
  while(is.null(fit_INLA)){
    
    fit_model <- function(){
      tryCatch(expr = {bru(components = model_components,
                           
                           like(family = "poisson",
                                formula = model_formula_OBBA2,
                                data = subset(PC_dat, year(Date_Time)%in% c(2001:2005))),
                           
                           options = list(bru_initial = list(lsig = 2,
                                                             Intercept_OBBA2 = -5),
                                          control.compute = list(waic = FALSE, cpo = FALSE),
                                          bru_verbose = 4))},
               error = function(e){NULL})
    }
    fit_INLA <- fit_model()
    
    if ("try-error" %in% class(fit_INLA)) fit_INLA <- NULL
  }
  
  end <- Sys.time()
  runtime_INLA <- difftime( end,start, units="mins") %>% round(2)
  print(paste0(species_name," - ",runtime_INLA," min to fit model"))  # ~6 min
  
  pred_formula_OBBA2 = as.formula(paste0(' ~
                  Intercept_OBBA2 +
                  log(15)+
                  range_effect_OBBA2 * distance_from_range +
                  spde_OBBA2 +',
                                      paste0("Beta1_",covariates_to_include,'*',covariates_to_include, collapse = " + "),
                                      '+',
                                      paste0("Beta2_",covariates_to_include,'*',covariates_to_include,"^2", collapse = " + ")))
  
  # Note that predictions are initially on log scale
  start2 <- Sys.time()
  
  pred_OBBA2 <- NULL
  pred_OBBA2 <- generate(fit_INLA,
                   as(ONGrid_species,'Spatial'),
                   formula =  pred_formula_OBBA2,
                   n.samples = 500,
                   seed = 123)
  
  pred_OBBA2 <- exp(pred_OBBA2)
  OBBA2_prediction_quantiles = apply(pred_OBBA2,1,function(x) quantile(x,c(0.05,0.5,0.95),na.rm = TRUE))
  ONGrid_species$OBBA2_pred_q05 <- OBBA2_prediction_quantiles[1,]
  ONGrid_species$OBBA2_pred_q50 <- OBBA2_prediction_quantiles[2,]
  ONGrid_species$OBBA2_pred_q95 <- OBBA2_prediction_quantiles[3,]
  ONGrid_species$OBBA2_pred_CI_width_90 <- OBBA2_prediction_quantiles[3,] - OBBA2_prediction_quantiles[1,]
  
  # ***************************************************************
  # OBBA - 3
  # ***************************************************************
  
  model_components = as.formula(paste0('~
            Intercept_OBBA3(1)+
            range_effect_OBBA3(1,model="linear", mean.linear = -0.046, prec.linear = 10000)+
            spde_OBBA3(main = coordinates, model = matern_abund) +',
                                       paste0("Beta1_",covariates_to_include,'(1,model="linear", mean.linear = 0, prec.linear = ', prec_linear[1],')', collapse = " + "),
                                       '+',
                                       paste0("Beta2_",covariates_to_include,'(1,model="linear", mean.linear = 0, prec.linear = ', prec_linear[2],')', collapse = " + ")))
  
  model_formula_OBBA3 = as.formula(paste0('count ~
                  Intercept_OBBA3 +
                  range_effect_OBBA3 * distance_from_range +
                  spde_OBBA3 +',
                                          paste0("Beta1_",covariates_to_include,'*',covariates_to_include, collapse = " + "),
                                          '+',
                                          paste0("Beta2_",covariates_to_include,'*',covariates_to_include,"^2", collapse = " + ")))
  
  start <- Sys.time()
  fit_INLA <- NULL
  while(is.null(fit_INLA)){
    
    fit_model <- function(){
      tryCatch(expr = {bru(components = model_components,
                           
                           like(family = "poisson",
                                formula = model_formula_OBBA3,
                                data = subset(PC_dat, year(Date_Time)%in% c(2021:2025))),
                           
                           options = list(bru_initial = list(lsig = 2,
                                                             Intercept_OBBA3 = -5),
                                          control.compute = list(waic = FALSE, cpo = FALSE),
                                          bru_verbose = 4))},
               error = function(e){NULL})
    }
    fit_INLA <- fit_model()
    
    if ("try-error" %in% class(fit_INLA)) fit_INLA <- NULL
  }
  
  end <- Sys.time()
  runtime_INLA <- difftime( end,start, units="mins") %>% round(2)
  print(paste0(species_name," - ",runtime_INLA," min to fit model"))  # ~6 min
  
  pred_formula_OBBA3 = as.formula(paste0(' ~
                  Intercept_OBBA3 +
                  log(15)+
                  range_effect_OBBA3 * distance_from_range +
                  spde_OBBA3 +',
                                         paste0("Beta1_",covariates_to_include,'*',covariates_to_include, collapse = " + "),
                                         '+',
                                         paste0("Beta2_",covariates_to_include,'*',covariates_to_include,"^2", collapse = " + ")))
  
  # Note that predictions are initially on log scale
  start2 <- Sys.time()
  
  pred_OBBA3 <- NULL
  pred_OBBA3 <- generate(fit_INLA,
                         as(ONGrid_species,'Spatial'),
                         formula =  pred_formula_OBBA3,
                         n.samples = 500,
                         seed = 123)
  
  pred_OBBA3 <- exp(pred_OBBA3)
  OBBA3_prediction_quantiles = apply(pred_OBBA3,1,function(x) quantile(x,c(0.05,0.5,0.95),na.rm = TRUE))
  ONGrid_species$OBBA3_pred_q05 <- OBBA3_prediction_quantiles[1,]
  ONGrid_species$OBBA3_pred_q50 <- OBBA3_prediction_quantiles[2,]
  ONGrid_species$OBBA3_pred_q95 <- OBBA3_prediction_quantiles[3,]
  ONGrid_species$OBBA3_pred_CI_width_90 <- OBBA3_prediction_quantiles[3,] - OBBA3_prediction_quantiles[1,]
  
  
  # ***************************************************************
  # Change between atlases
  # ***************************************************************
  
  pred_change <- pred_OBBA3 - pred_OBBA2
  change_prediction_quantiles = apply(pred_change,1,function(x) quantile(x,c(0.05,0.5,0.95),na.rm = TRUE))
  ONGrid_species$change_pred_q05 <- change_prediction_quantiles[1,]
  ONGrid_species$change_pred_q50 <- change_prediction_quantiles[2,]
  ONGrid_species$change_pred_q95 <- change_prediction_quantiles[3,]
  ONGrid_species$change_pred_CI_width_90 <- change_prediction_quantiles[3,] - change_prediction_quantiles[1,]
  
  end2 <- Sys.time() 
  runtime_pred <- difftime( end2,start2, units="mins") %>% round(2)
  print(paste0(species_name," - ",runtime_pred," min to generate predictions")) # 12 min
  
  # -----------------------------------------------------------------------------
  # Figures
  # -----------------------------------------------------------------------------
  
  # Convert to CRS of target raster
  ONGrid_species <- ONGrid_species %>% st_transform(st_crs(target_raster))
  
  colscale_relabund <- c("#FEFEFE", "#FBF7E2", "#FCF8D0", "#EEF7C2", "#CEF2B0", "#94E5A0", "#51C987", "#18A065", "#008C59", "#007F53", "#006344")
  colpal_relabund <- colorRampPalette(colscale_relabund)
  
  lower_bound <- 0.15
  upper_bound <- quantile(ONGrid_species$OBBA2_pred_q50,0.95,na.rm = TRUE) %>% signif(2)
  
  if (lower_bound >= (upper_bound/5)) lower_bound <- (upper_bound/5) %>% signif(2)
  
  sp_cut <- cut.fn(df = ONGrid_species,
                   target_raster = target_raster,
                   column_name = "OBBA2_pred_q50",
                   lower_bound = lower_bound,
                   upper_bound = upper_bound)
  
  raster_q50 <- sp_cut$raster 
  
  OBBA2_plot_q50 <- ggplot() +
    
    geom_stars(data = raster_q50) +
    scale_fill_manual(name = "Per 15 point counts",
                      values = colpal_relabund(length(levels(raster_q50$levs))), drop=FALSE,na.translate=FALSE)+
    
    geom_sf(data = ONBoundary,colour="black",fill=NA,lwd=0.3,show.legend = F) +
    
  coord_sf(clip = "off",xlim = range(as.data.frame(st_coordinates(ONBoundary))$X)) +
     theme(panel.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
     theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
     theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
     theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    annotate(geom="text",x=700000,y=1650000, label= paste0(species_name),lineheight = .85,hjust = 0,size=6,fontface =2) +
    annotate(geom="text",x=700000,y=1550000, label= "OBBA 2",lineheight = .85,hjust = 0,size=4,fontface =2) +
    annotate(geom="text",x=700000,y=1480000, label= paste0("Prepared on ",Sys.Date()),size=2,lineheight = .75,hjust = 0,color="gray75")
     
    print(OBBA2_plot_q50)
  
  png(paste0("../Output/Prediction_Maps/Relative_Abundance/",species_name,"_OBBA2_INDEP_q50.png"), width=6.5, height=8, units="in", res=300, type="cairo")
  print(OBBA2_plot_q50)
  dev.off()
  
  # ****************************************************************************
  # RELATIVE ABUNDANCE IN OBBA3
  # ****************************************************************************
  
  sp_cut <- cut.fn(df = ONGrid_species,
                   target_raster = target_raster,
                   column_name = "OBBA3_pred_q50",
                   lower_bound = lower_bound,
                   upper_bound = upper_bound)
  
  raster_q50 <- sp_cut$raster 
  
  OBBA3_plot_q50 <- ggplot() +
    
    geom_stars(data = raster_q50) +
    scale_fill_manual(name = "Per 15 point counts",
                      values = colpal_relabund(length(levels(raster_q50$levs))), drop=FALSE,na.translate=FALSE)+
    
    geom_sf(data = ONBoundary,colour="black",fill=NA,lwd=0.3,show.legend = F) +
    
  coord_sf(clip = "off",xlim = range(as.data.frame(st_coordinates(ONBoundary))$X)) +
    theme(panel.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    
    annotate(geom="text",x=700000,y=1650000, label= paste0(species_name),lineheight = .85,hjust = 0,size=6,fontface =2) +
    annotate(geom="text",x=700000,y=1550000, label= "OBBA 3",lineheight = .85,hjust = 0,size=4,fontface =2) +
    annotate(geom="text",x=700000,y=1480000, label= paste0("Prepared on ",Sys.Date()),size=2,lineheight = .75,hjust = 0,color="gray75")
  
  print(OBBA3_plot_q50)
  
  png(paste0("../Output/Prediction_Maps/Relative_Abundance/",species_name,"_OBBA3_INDEP_q50.png"), width=6.5, height=8, units="in", res=300, type="cairo")
  print(OBBA3_plot_q50)
  dev.off()
  
  # ****************************************************************************
  # CHANGE BETWEEN ATLASES
  # ****************************************************************************
  
  colscale_change <-c("orangered","white","dodgerblue")
  colpal_change <- colorRampPalette(colscale_change)
  
  cut_levs <- seq(-upper_bound/2,upper_bound/2,length.out = 6) %>% round(2)
  cut_levs <- c(-Inf,cut_levs,Inf)
  
  tmp <- ONGrid_species
  tmp$levs <- cut(as.data.frame(tmp)[,"change_pred_q50"], cut_levs)
  raster_q50 <- stars::st_rasterize(tmp %>% dplyr::select(levs, geometry))
  
  change_plot_q50 <- ggplot() +
    
    geom_stars(data = raster_q50) +
    scale_fill_manual(name = "Change",
                      values = colpal_change(length(levels(raster_q50$levs))), drop=FALSE,na.translate=FALSE)+
    
    #geom_sf(data = ONWater,colour=NA,fill="#59F3F3",show.legend = F)+
    
    geom_sf(data = ONBoundary,colour="black",fill=NA,lwd=0.3,show.legend = F) +
    
    #geom_sf(data = range,colour="darkred",fill=NA,show.legend = F) +
    
    # Surveyed squares
    #geom_sf(data = subset(ONSquares_centroids, !is.na(PC_detected) | !is.na(CL_detected)), col = "gray70", pch = 19, size = 0.1)+
    
    # Point count detections
    #geom_sf(data = subset(ONSquares_centroids, !is.na(PC_detected) & PC_detected > 0), col = "black", pch = 19, size = 0.2)+
    
    # Checklist detections
    #geom_sf(data = subset(ONSquares_species, !is.na(CL_detected) & CL_detected > 0), col = "black", size = 0.2, fill = "transparent")+
  
  coord_sf(clip = "off",xlim = range(as.data.frame(st_coordinates(ONBoundary))$X)) +
    theme(panel.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    # theme(legend.margin=margin(0,0,0,0),legend.box.margin=margin(5,10,5,-20),legend.title.align=0.5,
    #       legend.title = element_markdown(lineheight=.9,hjust = "left"))
    annotate(geom="text",x=700000,y=1650000, label= paste0(species_name),lineheight = .85,hjust = 0,size=6,fontface =2) +
    annotate(geom="text",x=700000,y=1550000, label= "OBBA 3",lineheight = .85,hjust = 0,size=4,fontface =2) +
    annotate(geom="text",x=700000,y=1480000, label= paste0("Prepared on ",Sys.Date()),size=2,lineheight = .75,hjust = 0,color="gray75")
  
  print(change_plot_q50)
  
  png(paste0("../Output/Prediction_Maps/Relative_Abundance/",species_name,"_change_INDEP_q50.png"), width=6.5, height=8, units="in", res=300, type="cairo")
  print(change_plot_q50)
  dev.off()
  
  
  # Estimate of total change
  total_change <- pred_change %>% apply(2,sum,na.rm = TRUE)
  hist(total_change, breaks = seq(-1000000,1000000,length.out = 25))
  
} # close species loop


