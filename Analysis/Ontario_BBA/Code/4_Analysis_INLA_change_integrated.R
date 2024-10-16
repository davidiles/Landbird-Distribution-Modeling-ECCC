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
  
  start <- Sys.time()
  
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
  
  # Controls the 'residual spatial field'.  This can be adjusted to create smoother surfaces.
  matern_change <- inla.spde2.pcmatern(mesh_spatial,
                                      prior.range = c(1000000, 0.1), # 10% chance range is smaller than 300000
                                      prior.sigma = c(0.1,0.1),    # 10% chance sd is larger than 0.1,
                                      constr = TRUE
  )
  
  # ----------------------------------------------------
  # Create mesh to model effect of time since sunrise (TSS)
  # Hours since sunrise is fit as a GAM-type effect.  
  # Recommend leaving these settings at these defaults
  # ----------------------------------------------------
  # 
  # TSS_range <- range(sp_dat$Hours_Since_Sunrise)
  # TSS_meshpoints <- seq(TSS_range[1]-0.1,TSS_range[2]+0.1,length.out = 11)
  # TSS_mesh1D = inla.mesh.1d(TSS_meshpoints,boundary="free")
  # TSS_spde = inla.spde2.pcmatern(TSS_mesh1D,
  #                                prior.range = c(6,0.1), # 10% range is smaller than 6
  #                                prior.sigma = c(1,0.1)) # 10% chance sd is larger than 1
  # 
  # ***************************************************************
  # Model formulas
  # ***************************************************************
  
  # Names of covariates to include in model. Ideally, covariates should be uncorrelated with each other.
  # Each covariate will include a quadratic effect to allow for 'intermediate optimum' effects
  covariates_to_include <- paste0("PC",1:2) 
  
  # How much shrinkage should be applied to covariate effects?
  sd_linear <- 0.1  # Change to smaller value (e.g., 0.1), if you want to heavily shrink covariate effects and potentially create smoother surfaces
  prec_linear <-  c(1/sd_linear^2,1/(sd_linear/2)^2)
  
  #TSS(main = Hours_Since_Sunrise,model = TSS_spde) +
  model_components = as.formula(paste0('~
            Intercept_OBBA2(1)+
            Intercept_OBBA3(1,model="linear", mean.linear = 0, prec.linear = 100)+
            range_effect_OBBA2(1,model="linear", mean.linear = -0.046, prec.linear = 10000)+
            range_effect_OBBA3(1,model="linear", mean.linear = 0, prec.linear = 10000)+
            
            spde_OBBA2(main = coordinates, model = matern_abund) +
            spde_change(main = coordinates, model = matern_change) +',
                                       paste0("Beta1_",covariates_to_include,'(1,model="linear", mean.linear = 0, prec.linear = ', prec_linear[1],')', collapse = " + "),
                                       '+',
                                       paste0("Beta2_",covariates_to_include,'(1,model="linear", mean.linear = 0, prec.linear = ', prec_linear[2],')', collapse = " + "),
                                       '+',
                                       paste0("Beta1_change_",covariates_to_include,'(1,model="linear", mean.linear = 0, prec.linear = 2500)', collapse = " + "),
                                       '+',
                                       paste0("Beta2_change_",covariates_to_include,'(1,model="linear", mean.linear = 0, prec.linear = 2500)', collapse = " + ")))
  
  # TSS +
  model_formula_OBBA2 = as.formula(paste0('count ~
                  Intercept_OBBA2 +
                  range_effect_OBBA2 * distance_from_range +
                  spde_OBBA2 +',
                                       paste0("Beta1_",covariates_to_include,'*',covariates_to_include, collapse = " + "),
                                       '+',
                                       paste0("Beta2_",covariates_to_include,'*',covariates_to_include,"^2", collapse = " + ")))
  
  # TSS +
  model_formula_OBBA3 = as.formula(paste0('count ~
                  Intercept_OBBA2 +
                  Intercept_OBBA3 +
                  (range_effect_OBBA2 + range_effect_OBBA3) * distance_from_range +
                  spde_OBBA2 + spde_change + ',
                                          paste0("Beta1_",covariates_to_include,'*',covariates_to_include, collapse = " + "),
                                          '+',
                                          paste0("Beta2_",covariates_to_include,'*',covariates_to_include,"^2", collapse = " + "),
                                          '+',
                                          paste0("Beta1_change_",covariates_to_include,'*',covariates_to_include, collapse = " + "),
                                          '+',
                                          paste0("Beta2_change_",covariates_to_include,'*',covariates_to_include,"^2", collapse = " + ")))
  
  # ----------------------------------------------------------------------------
  # Fit model to point counts and ARUs with negative binomial error
  # ----------------------------------------------------------------------------
  
  fit_INLA <- NULL
  while(is.null(fit_INLA)){
    
    fit_model <- function(){
      tryCatch(expr = {bru(components = model_components,
                           
                           like(family = "poisson",
                                formula = model_formula_OBBA2,
                                data = subset(PC_dat, year(Date_Time)%in% c(2001:2005))),
                           
                           like(family = "poisson",
                                formula = model_formula_OBBA3,
                                data = subset(PC_dat, year(Date_Time)%in% c(2021:2025))),
                           
                           
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
  
  # ****************************************************************************
  # ****************************************************************************
  # GENERATE MAPS
  # ****************************************************************************
  # ****************************************************************************
  
  # For every pixel on landscape, extract distance (in km) from eBird range limit
  ONGrid_species <- ONGrid %>%
    mutate(distance_from_range = (st_centroid(.) %>% 
                                    st_distance( . , range) %>% 
                                    as.numeric())/1000)
  
  # ----------------------------------------------------
  # QPAD offsets associated with a 5-minute unlimited distance survey
  # ----------------------------------------------------
  # 
  # species_offsets <- subset(species_to_model, Species_Code_BSC == sp_code)
  # log_offset_5min <- 0
  # if (species_offsets$offset_exists == TRUE) log_offset_5min <- species_offsets$log_offset_5min
  
  
  # ----------------------------------------------------
  # Generate predictions on ONGrid_species raster
  # Per 15 point counts; hence the log(15) in the prediction formula
  # ----------------------------------------------------
  
  pred_formula_OBBA2 = as.formula(paste0(' ~
                  Intercept_OBBA2 +
                  log(15)+
                  range_effect_OBBA2 * distance_from_range +
                  spde_OBBA2 +',
                                      paste0("Beta1_",covariates_to_include,'*',covariates_to_include, collapse = " + "),
                                      '+',
                                      paste0("Beta2_",covariates_to_include,'*',covariates_to_include,"^2", collapse = " + ")))
  
  pred_formula_OBBA3 = as.formula(paste0(' ~
                  Intercept_OBBA2 +
                  Intercept_OBBA3 +
                  log(15)+
                  (range_effect_OBBA2 + range_effect_OBBA3) * distance_from_range +
                  spde_OBBA2 + spde_change +',
                                         paste0("Beta1_",covariates_to_include,'*',covariates_to_include, collapse = " + "),
                                         '+',
                                         paste0("Beta2_",covariates_to_include,'*',covariates_to_include,"^2", collapse = " + "),
                                         '+',
                                         paste0("Beta1_change_",covariates_to_include,'*',covariates_to_include, collapse = " + "),
                                         '+',
                                         paste0("Beta2_change_",covariates_to_include,'*',covariates_to_include,"^2", collapse = " + ")))
  
  
  pred_formula_change = as.formula(paste0(' ~
                  exp(Intercept_OBBA2 +
                  Intercept_OBBA3 +
                  log(15)+
                  (range_effect_OBBA2 + range_effect_OBBA3) * distance_from_range +
                  spde_OBBA2 + spde_change +',
                                          paste0("Beta1_",covariates_to_include,'*',covariates_to_include, collapse = " + "),
                                          '+',
                                          paste0("Beta2_",covariates_to_include,'*',covariates_to_include,"^2", collapse = " + "),
                                          '+',
                                          paste0("Beta1_change_",covariates_to_include,'*',covariates_to_include, collapse = " + "),
                                          '+',
                                          paste0("Beta2_change_",covariates_to_include,'*',covariates_to_include,"^2", collapse = " + "),
                                         ') - 
                                         
                  exp(Intercept_OBBA2 +
                  log_offset_5min +
                  log(15)+
                  range_effect_OBBA2 * distance_from_range +
                  spde_OBBA2 +',
                                      paste0("Beta1_",covariates_to_include,'*',covariates_to_include, collapse = " + "),
                                      '+',
                                      paste0("Beta2_",covariates_to_include,'*',covariates_to_include,"^2", collapse = " + "),
                                      ')'))
  
  # Note that predictions are initially on log scale
  start2 <- Sys.time()
  
  # --------
  # OBBA-2
  # --------
  
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
  
  # --------
  # OBBA-3
  # --------
  
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

  # --------
  # Change
  # --------
  
  pred_change <- NULL
  pred_change <- generate(fit_INLA,
                         as(ONGrid_species,'Spatial'),
                         formula =  pred_formula_change,
                         n.samples = 500,
                         seed = 123)
  
  change_prediction_quantiles = apply(pred_change,1,function(x) quantile(x,c(0.05,0.5,0.95),na.rm = TRUE))
  ONGrid_species$change_pred_q05 <- change_prediction_quantiles[1,]
  ONGrid_species$change_pred_q50 <- change_prediction_quantiles[2,]
  ONGrid_species$change_pred_q95 <- change_prediction_quantiles[3,]
  ONGrid_species$change_pred_CI_width_90 <- change_prediction_quantiles[3,] - change_prediction_quantiles[1,]
  
  end2 <- Sys.time() 
  runtime_pred <- difftime( end2,start2, units="mins") %>% round(2)
  print(paste0(species_name," - ",runtime_pred," min to generate predictions")) # 12 min
  
  
  # ------------------------------------------------
  # Summarize ONSquares where species was detected (for plotting)
  # ------------------------------------------------
  # 
  # PC_detected <- sp_dat %>%
  #   subset(Survey_Type %in% c("Point_Count","ARU_SPT")) %>% 
  #   as.data.frame() %>%
  #   group_by(sq_id) %>%
  #   summarize(PC_detected = as.numeric(sum(count)>0),
  #             PC_mean_count = mean(count) %>% round(2))
  # 
  # CL_detected <- sp_dat %>%
  #   subset(Survey_Type %in% c("Linear transect")) %>%
  #   as.data.frame() %>%
  #   group_by(sq_id) %>%
  #   summarize(CL_detected = as.numeric(sum(count)>0),
  #             CL_mean_count = mean(count))
  # 
  # ONSquares_species <- ONSquares %>%
  #   relocate(geometry,.after = last_col()) %>%
  #   left_join(PC_detected) #%>% 
  #   #left_join(CL_detected)
  # 
  # ONSquares_centroids <- st_centroid(ONSquares_species)
  # 
  # ------------------------------------------------
  # Label for figure
  # ------------------------------------------------
  
  #species_name = ON_spcd$CommonName[which(ON_spcd$spcd == sp_code)]
  #species_label = ON_spcd$Label[which(ON_spcd$spcd == sp_code)]
  
  # ------------------------------------------------
  # sf object for ebird range limit (optional - not needed for plotting, but potentially helpful)
  # ------------------------------------------------
  # 
  # if (!is.na(range)) range <- range  %>% 
  #   st_transform(st_crs(ONBoundary)) %>%
  #   st_intersection(ONBoundary)
  
  # ****************************************************************************
  # RELATIVE ABUNDANCE IN OBBA2
  # ****************************************************************************
  
  # Convert to CRS of target raster
  ONGrid_species <- ONGrid_species %>% st_transform(st_crs(target_raster))
  
  colscale_relabund <-c("#FEFEFE", "#FBF7E2", "#FCF8D0", "#EEF7C2", "#CEF2B0", "#94E5A0", "#51C987", "#18A065", "#008C59", "#007F53", "#006344")
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
    annotate(geom="text",x=700000,y=1550000, label= "OBBA 2",lineheight = .85,hjust = 0,size=4,fontface =2) +
    annotate(geom="text",x=700000,y=1480000, label= paste0("Prepared on ",Sys.Date()),size=2,lineheight = .75,hjust = 0,color="gray75")
    
  #   
    print(OBBA2_plot_q50)
  
  png(paste0("../Output/Prediction_Maps/Relative_Abundance/",species_name,"_OBBA2_q50.png"), width=6.5, height=8, units="in", res=300, type="cairo")
  print(OBBA2_plot_q50)
  dev.off()
  
  # ****************************************************************************
  # RELATIVE ABUNDANCE IN OBBA3
  # ****************************************************************************
  
  colscale_relabund <-c("#FEFEFE", "#FBF7E2", "#FCF8D0", "#EEF7C2", "#CEF2B0", "#94E5A0", "#51C987", "#18A065", "#008C59", "#007F53", "#006344")
  colpal_relabund <- colorRampPalette(colscale_relabund)
  
  lower_bound <- 0.15
  upper_bound <- quantile(ONGrid_species$OBBA2_pred_q50,0.95,na.rm = TRUE) %>% signif(2)
  
  if (lower_bound >= (upper_bound/5)) lower_bound <- (upper_bound/5) %>% signif(2)
  
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
  
  #   
  print(OBBA3_plot_q50)
  
  png(paste0("../Output/Prediction_Maps/Relative_Abundance/",species_name,"_OBBA3_q50.png"), width=6.5, height=8, units="in", res=300, type="cairo")
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
  
  #   
  print(change_plot_q50)
  
  png(paste0("../Output/Prediction_Maps/Relative_Abundance/",species_name,"_change_q50.png"), width=6.5, height=8, units="in", res=300, type="cairo")
  print(change_plot_q50)
  dev.off()
  
  
  # Estimate of total change
  total_change <- pred_change %>% apply(2,sum,na.rm = TRUE)
  hist(total_change, breaks = seq(-1000000,1000000,length.out = 25))
  
  # # ****************************************************************************
  # # FIGURE 2: UNCERTAINTY EXPRESSED AS WIDTH OF 90% CREDIBLE INTERVAL
  # # ****************************************************************************
  # 
  # colscale_uncertainty <- c("#FEFEFE", "#FFF4B3", "#F5D271", "#F2B647", "#EC8E00", "#CA302A")
  # colpal_uncertainty <- colorRampPalette(colscale_uncertainty)
  # 
  # lower_bound <- 0.15
  # upper_bound <- 5#quantile(ONGrid_species$pred_CI_width_90,0.99,na.rm = TRUE) %>% signif(2)
  # if (lower_bound >= (upper_bound/5)) lower_bound <- (upper_bound/5) %>% signif(2)
  # 
  # raster_CI_width_90 <- cut.fn(df = ONGrid_species,
  #                              target_raster = target_raster,
  #                              column_name = "pred_CI_width_90",
  #                              lower_bound = lower_bound,
  #                              upper_bound = upper_bound)$raster
  # 
  # plot_CI_width_90 <- ggplot() +
  #   
  #   geom_stars(data = raster_CI_width_90) +
  #   scale_fill_manual(name = "<span style='font-size:13pt'>Uncertainty</span><br><span style='font-size:7pt'>Per 15 point counts</span><br><span style='font-size:7pt'>Width of 90% CI</span>",
  #                     values = colpal_uncertainty(length(levels(raster_CI_width_90$levs))), drop=FALSE,na.translate=FALSE)+
  #   
  #   geom_sf(data = ONWater,colour=NA,fill="#59F3F3",show.legend = F)+
  #   
  #   geom_sf(data = ONBoundary,colour="black",fill=NA,lwd=0.3,show.legend = F) +
  #   
  #   #geom_sf(data = range,colour="darkred",fill=NA,show.legend = F) +
  #   
  #   # Surveyed squares
  #   #geom_sf(data = subset(ONSquares_centroids, !is.na(PC_detected) | !is.na(CL_detected)), col = "gray70", pch = 19, size = 0.1)+
  #   
  #   # Point count detections
  #   #geom_sf(data = subset(ONSquares_centroids, !is.na(PC_detected) & PC_detected > 0), col = "black", pch = 19, size = 0.2)+
  #   
  #   # Checklist detections
  #   #geom_sf(data = subset(ONSquares_centroids, !is.na(CL_detected) & CL_detected > 0), col = "black", pch = 4, size = 0.2)+
  # 
  # coord_sf(clip = "off",xlim = range(as.data.frame(st_coordinates(ONBoundary))$X))+
  #   theme(panel.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  #   theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
  #   theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
  #   theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
  #   theme(legend.margin=margin(0,0,0,0),legend.box.margin=margin(5,10,5,-20),legend.title.align=0.5,
  #         legend.title = element_markdown(lineheight=.9,hjust = "left"))+
  #   annotate(geom="text",x=346000,y=1850000, label= paste0(species_label),lineheight = .85,hjust = 0,size=6,fontface =2) +
  #   annotate(geom="text",x=346000,y=1400000, label= paste0("Prepared on ",Sys.Date()),size=3,lineheight = .75,hjust = 0,color="#3b3b3b")
  # print(plot_CI_width_90)
  # 
  # png(paste0("../Output/Prediction_Maps/Uncertainty/",sp_code,"_CI_width_90.png"), width=6.5, height=8, units="in", res=300, type="cairo")
  # print(plot_CI_width_90)
  # dev.off()
  # 
  # # ****************************************************************************
  # # FIGURE 3: PROBABILITY OF OBSERVING SPECIES AFTER CONDUCTING 15 POINT COUNTS
  # # ****************************************************************************
  # 
  # colscale_pObs <- c("#FEFEFE",RColorBrewer::brewer.pal(5,"BuGn")[2:5])
  # colpal_pObs <- colorRampPalette(colscale_pObs)
  # 
  # cut_levs <- c(-0.1,0.05,0.125,0.25,0.5,1)
  # cut_levs_labs <- c("0 to 0.05", 
  #                    "0.05 to 0.125", 
  #                    "0.125 to 0.25", 
  #                    "0.25 to 0.50", 
  #                    "0.50 to 1")
  # 
  # ONGrid_species$pObs_levs <- cut(as.data.frame(ONGrid_species)[,"pObs_5min"], 
  #                                 cut_levs,labels=cut_levs_labs)
  # 
  # raster_pObs = stars::st_rasterize(ONGrid_species %>% dplyr::select(pObs_levs, geometry))
  # 
  # plot_pObs <- ggplot() +
  #   
  #   geom_stars(data = raster_pObs) +
  #   scale_fill_manual(name = "<span style='font-size:13pt'>Prob. of Observation</span><br><span style='font-size:7pt'>Per 15 point counts </span><br><span style='font-size:7pt'>(Posterior Median)</span>",
  #                     values = colpal_pObs(length(levels(raster_pObs$pObs_levs))), 
  #                     drop=FALSE,na.translate=FALSE)+
  #   
  #   geom_sf(data = ONWater,colour=NA,fill="#59F3F3",show.legend = F)+
  #   
  #   geom_sf(data = ONBoundary,colour="black",fill=NA,lwd=0.3,show.legend = F) +
  #   
  #   #geom_sf(data = range,colour="darkred",fill=NA,show.legend = F) +
  #   
  #   # Surveyed squares
  #   #geom_sf(data = subset(ONSquares_centroids, !is.na(PC_detected) | !is.na(CL_detected)), col = "gray70", pch = 19, size = 0.1)+
  #   
  #   # Point count detections
  #   #geom_sf(data = subset(ONSquares_centroids, !is.na(PC_detected) & PC_detected > 0), col = "black", pch = 19, size = 0.2)+
  #   
  #   # Checklist detections
  #   #geom_sf(data = subset(ONSquares_centroids, !is.na(CL_detected) & CL_detected > 0), col = "black", pch = 4, size = 0.2)+
  # 
  # coord_sf(clip = "off",xlim = range(as.data.frame(st_coordinates(ONBoundary))$X))+
  #   theme(panel.background = element_blank(),
  #         panel.grid.major = element_blank(), 
  #         panel.grid.minor = element_blank())+
  #   theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
  #   theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
  #   theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
  #   theme(legend.margin=margin(0,0,0,0),
  #         legend.box.margin=margin(5,10,5,-20),
  #         legend.title.align=0.5,
  #         legend.title = element_markdown(lineheight=.9,hjust = "left"))+
  #   annotate(geom="text",x=346000,y=1850000, label= paste0(species_label),lineheight = .85,hjust = 0,size=6,fontface =2) +
  #   annotate(geom="text",x=346000,y=1400000, label= paste0("Prepared on ",Sys.Date()),size=3,lineheight = .75,hjust = 0,color="#3b3b3b")
  # 
  # png(paste0("../Output/Prediction_Maps/PObs/",sp_code,"_PObs.png"), width=6.5, height=8, units="in", res=300, type="cairo")
  # print(plot_pObs)
  # dev.off()
  # 
  # # ****************************************************************************
  # # BONUS FIGURE: PObs based on stationary count checklists
  # # ****************************************************************************
  # 
  # # Mean effort per square (could instead use 1,10,100 or whatever other number is relevant)
  # ONGrid_species$total_effort <- mean(SC_dat$total_effort)
  # 
  # pred_formula_SC = as.formula(paste0(' ~
  #                 Intercept_SC +
  #                 SC_effort +
  #                 range_effect * distance_from_range +
  #                 spde_coarse +',
  #                                     paste0("Beta1_",covariates_to_include,'*',covariates_to_include, collapse = " + "),
  #                                     '+',
  #                                     paste0("Beta2_",covariates_to_include,'*',covariates_to_include, collapse = " + ")))
  # 
  # # Note that predictions are initially on log scale
  # pred_SC <- NULL
  # pred_SC <- generate(fit_INLA,
  #                     as(ONGrid_species,'Spatial'),
  #                     formula =  pred_formula_SC,
  #                     n.samples = 1000)
  # 
  # pred_SC <- exp(pred_SC)
  # 
  # prediction_quantiles_SC = apply(pred_SC,1,function(x) quantile(x,c(0.05,0.5,0.95),na.rm = TRUE))
  # 
  # # Median, upper,and lower 90% credible intervals on PObs
  # ONGrid_species$pObs_SC <- 1-dpois(0,lambda=prediction_quantiles_SC[2,])
  # ONGrid_species$pObs_SC_q05 <- 1-dpois(0,lambda=prediction_quantiles_SC[1,])
  # ONGrid_species$pObs_SC_q95 <- 1-dpois(0,lambda=prediction_quantiles_SC[3,])
  # 
  # ONGrid_species$pObs_SC_levs <- cut(as.data.frame(ONGrid_species)[,"pObs_SC"], 
  #                                    cut_levs,labels=cut_levs_labs)
  # 
  # raster_pObs_SC = stars::st_rasterize(ONGrid_species %>% dplyr::select(pObs_SC_levs, geometry))
  # 
  # plot_pObs_SC <- ggplot() +
  #   
  #   geom_stars(data = raster_pObs_SC) +
  #   scale_fill_manual(name = "<span style='font-size:13pt'>Prob. of Observation</span><br><span style='font-size:7pt'> Based on stationary checklists </span><br><span style='font-size:7pt'>(Posterior Median)</span>",
  #                     values = colpal_pObs(length(levels(raster_pObs_SC$pObs_SC_levs))), 
  #                     drop=FALSE,na.translate=FALSE)+
  #   
  #   geom_sf(data = ONWater,colour=NA,fill="#59F3F3",show.legend = F)+
  #   
  #   geom_sf(data = ONBoundary,colour="black",fill=NA,lwd=0.3,show.legend = F) +
  #   
  #   #geom_sf(data = range,colour="darkred",fill=NA,show.legend = F) +
  #   
  #   # Surveyed squares
  #   #geom_sf(data = subset(ONSquares_centroids, !is.na(PC_detected) | !is.na(CL_detected)), col = "gray70", pch = 19, size = 0.1)+
  #   
  #   # Point count detections
  #   #geom_sf(data = subset(ONSquares_centroids, !is.na(PC_detected) & PC_detected > 0), col = "black", pch = 19, size = 0.2)+
  #   
  #   # Checklist detections
  #   #geom_sf(data = subset(ONSquares_centroids, !is.na(CL_detected) & CL_detected > 0), col = "black", pch = 4, size = 0.2)+
  # 
  # coord_sf(clip = "off",xlim = range(as.data.frame(st_coordinates(ONBoundary))$X))+
  #   theme(panel.background = element_blank(),
  #         panel.grid.major = element_blank(), 
  #         panel.grid.minor = element_blank())+
  #   theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
  #   theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
  #   theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
  #   theme(legend.margin=margin(0,0,0,0),
  #         legend.box.margin=margin(5,10,5,-20),
  #         legend.title.align=0.5,
  #         legend.title = element_markdown(lineheight=.9,hjust = "left"))+
  #   annotate(geom="text",x=346000,y=1850000, label= paste0(species_label),lineheight = .85,hjust = 0,size=6,fontface =2) +
  #   annotate(geom="text",x=346000,y=1400000, label= paste0("Prepared on ",Sys.Date()),size=3,lineheight = .75,hjust = 0,color="#3b3b3b")
  # 
  # png(paste0("../Output/Prediction_Maps/PObs/",sp_code,"_PObs_SC.png"), width=6.5, height=8, units="in", res=300, type="cairo")
  # print(plot_pObs_SC)
  # dev.off()
  
} # close species loop


