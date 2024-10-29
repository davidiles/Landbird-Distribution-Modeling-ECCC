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
require(ggspatial)

# Timeout INLA after 10 minutes (if it has not fit by then, it has likely stalled)
inla.setOption(inla.timeout = 60*10)

# ------------------------------------------------
# Function to rasterize a series of spatial predictions (needed for plotting)
# ------------------------------------------------

cut.fn <- function(df = NA, 
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

AEA_proj <- "+proj=aea +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-87 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs "

ONBoundary <- st_read("../../../Data/Spatial/National/BCR/BCR_Terrestrial_master.shp")  %>%
  subset(PROVINCE_S == "ONTARIO") %>%
  st_make_valid() %>%
  st_union() %>%
  st_transform(AEA_proj) %>%
  st_buffer(5) %>%
  st_buffer(-5) 

# ONGrid <- ONGrid %>% st_intersection(ONBoundary)
ONGrid <- ONGrid %>% rename(geometry = x)

# Raster with target properties (needed for plotting)
target_raster <- rast("../../../Data/Spatial/National/AnnualMeanTemperature/wc2.1_30s_bio_1.tif") %>% 
  crop(st_transform(ONBoundary,crs(.))) %>% 
  project(AEA_proj, res = 1000) %>%
  mask(vect(ONBoundary))

# ----------------------------------------------------
# Create hexagonal grid across study area
# ----------------------------------------------------

hexgrid <- st_make_grid(ONBoundary, cellsize = 50, square=FALSE, what = "polygons") %>%
  st_as_sf() %>% mutate(hexid = 1:nrow(.)) %>% dplyr::rename(geometry = x) %>%
  st_intersection(ONBoundary)

hex_centroids <- st_centroid(hexgrid)

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
                      yday(Date_Time) <= yday(ymd("2022-07-15")) &
                      Project_Name %in% c("OBBA2",
                                          "OBBA3",
                                          "Ontario Breeding Bird Atlas Digital Point Counts 2023",
                                          "Ontario Breeding Bird Atlas Digital Point Counts 2022",
                                          "Ontario Breeding Bird Atlas Digital Point Counts 2021")
)

dim(PC_to_use) # 101206

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

set.seed(123)
species_list <- subset(species_to_model, n_detections > 500) %>%
  sample_n(50)

species_list <- rbind(subset(species_to_model,english_name %in% 
                               c("Canada Jay","Lesser Yellowlegs","Boreal Chickadee","Black-capped Chickadee","Rusty Blackbird","Tennessee Warbler")),
                      species_list) %>% unique()

# Prepare QPAD offsets for survey effort
require(napops)
napops_species <- list_species() %>% rename(Species_Code_NAPOPS = Species,Common_Name_NAPOPS = Common_Name,Scientific_Name_NAPOPS = Scientific_Name)

species_list <- left_join(species_list,napops_species[,c("Species_Code_NAPOPS","Common_Name_NAPOPS","Removal","Distance")], by = c("english_name" = "Common_Name_NAPOPS")) %>%
  mutate(EDR = NA, cue_rate = NA)

for (i in 1:nrow(species_list)){
  
  # ----------------------------------------------------
  # Extract counts/data for this species
  # ----------------------------------------------------
  
  sp_code = as.character(species_list$Species_Code_BSC[i])
  species_name = species_list$english_name[i]
  print(species_name)
  if (file.exists(paste0("../Output/Prediction_Maps/Relative_Abundance/",species_name,"_change_q50.png"))) next
  
  # Prepare data for this species
  sp_dat <- all_surveys %>% 
    mutate(count = full_count_matrix[,sp_code],
           presence = as.numeric(full_count_matrix[,sp_code]>0)) %>%
    st_transform(AEA_proj)
  
  # ----------------------------------------------------
  # Extract ebird range for this species (if it exists); prepared by previous script
  # ----------------------------------------------------
  
  range <- NA
  if (species_name %in% names(species_ranges)){
    
    range <- species_ranges[[species_name]] %>% st_transform(st_crs(sp_dat))
    
    # Identify distance of each survey to the edge of species range (in km)
    sp_dat$distance_from_range <- ((st_distance(sp_dat, range)) %>% as.numeric())
  } else{
    sp_dat$distance_from_range <- 0
  }
  
  # ----------------------------------------------------
  # Generate QPAD offsets for each survey (assumes unlimited distance point counts)
  # ----------------------------------------------------
  
  if (!is.na(species_list$Removal[i]) & !is.na(species_list$Distance[i]) &
      species_list$Removal[i] == 1 & species_list$Distance[i] == 1){
    species_list$cue_rate[i] <- cue_rate(species = species_list$Species_Code_NAPOPS[i],od = 153, tssr = 0, model = 1)[3] %>% as.numeric()
    species_list$EDR[i] <- edr(species = species_list$Species_Code_NAPOPS[i],road = FALSE, forest = 0.5,model = 1)[3] %>% as.numeric()
    
    # Calculate A and p, which jointly determine offset
    A_metres <- c(pi*species_list$EDR[i]^2)
    p <- 1-exp(-5*species_list$cue_rate[i])
    
    sp_dat$log_offset <- log(A_metres * (1-exp(-sp_dat$Survey_Duration_Minutes*species_list$cue_rate[i])))
    log_offset_5min <- log(A_metres * p)
  } else {
    sp_dat$log_offset <- 0
    log_offset_5min <- 0
  }
  
  # ----------------------------------------------------
  # Generate plot of observed data in each hexagon
  # ----------------------------------------------------
  
  # Intersect with survey information
  species_hex <- sp_dat %>% st_intersection(hexgrid)
  species_hex$Atlas <- "OBBA2"
  species_hex$Atlas[year(species_hex$Date_Time) >= 2020] <- "OBBA3"
  
  # Summarize total effort and mean count in each hexagon
  species_hex_centroid <- species_hex %>%
    as.data.frame() %>%
    group_by(hexid,Atlas) %>%
    summarize(count = sum(count),
              effort = sum(Survey_Duration_Minutes)) 
  species_hex_centroid <- full_join(hex_centroids,species_hex_centroid) %>% na.omit()
  
  col_lim <- c(0,max(species_hex_centroid$count,na.rm = TRUE))
  
  species_data_plot <- ggplot() +
    geom_sf(data=ONBoundary,colour="gray70", fill = "gray80")+
    geom_sf(data=species_hex_centroid, aes(col = count, size = effort))+
    #geom_sf(data=subset(species_hex_centroid,count>0), col = "black", fill = "transparent",size = 0.1)+
    
    annotation_scale(style = "ticks",
                     text_face = "bold")+
    theme_bw()+
    
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = 'transparent', colour = 'black'))+
    scale_color_gradientn(colors = c("black",viridis(10)),na.value=NA)+
    ggtitle(paste0(species_name," - Observed count per km^2"))+
    facet_grid(Atlas~.)
  species_data_plot
  
  # ----------------------------------------------------
  # Create a spatial mesh, which is used to fit the residual spatial field
  # NOTE: MAKING A COARSER MESH (LARGER VALUES OF MAX_EDGE) WILL SPEED UP COMPUTATIONS AND REDUCE RAM REQUIREMENTS, 
  # BUT WITH POTENTIAL LOSS OF ACCURACY/PRECISION
  # DEFAULT SHOULD CREATE APPROX 5000
  # ----------------------------------------------------
  
  # make a two extension hulls and mesh for spatial model
  hull <- fm_extensions(
    ONBoundary,
    convex = c(50, 200),
    concave = c(350, 500)
  )
  
  # Setting max.edge = c(30000,100000) requires longer to fit (7 min), but may increase model accuracy
  # Setting max.edge = c(50000,100000) allows the model to fit much more quickly (3.5 min), though possibly with reduced accuracy
  
  mesh_spatial <- fm_mesh_2d_inla(
    #loc = st_as_sfc(sp_dat),
    boundary = hull, 
    max.edge = c(50, 200), #c(50000, 200000), # km inside and outside  
    cutoff = 10, 
    crs = fm_crs(species_hex)
  )
  
  mesh_locs <- mesh_spatial$loc[,c(1,2)] %>% as.data.frame()
  dim(mesh_locs)
  
  # # Controls the 'residual spatial field'.  This can be adjusted to create smoother surfaces.
  # matern_abund <- inla.spde2.pcmatern(mesh_spatial,
  #                                     prior.range = c(300, 0.1), # 500, NA # 10% chance range is smaller than 500000
  #                                     prior.sigma = c(0.05,0.1),   # 0.1, 0.1 # 10% chance sd is larger than 1,
  #                                     constr = TRUE
  # )
  # 
  # # Controls the 'residual spatial field'.  This can be adjusted to create smoother surfaces.
  # matern_change <- inla.spde2.pcmatern(mesh_spatial,
  #                                      prior.range = c(600, 0.1), # 500, NA # 10% chance range is smaller than 500000
  #                                      prior.sigma = c(0.05,0.1),   # 0.1, 0.1# 10% chance sd is larger than 1,
  #                                      constr = TRUE
  # )
  
  # Controls the 'residual spatial field'.  This can be adjusted to create smoother surfaces.
  matern_abund <- inla.spde2.pcmatern(mesh_spatial,
                                      prior.range = c(250, 0.01), # 500, NA # 10% chance range is smaller than 500000
                                      prior.sigma = c(0.2,0.1),   # 0.1, 0.1 # 10% chance sd is larger than 1,
                                      constr = TRUE
  )
  
  # Controls the 'residual spatial field'.  This can be adjusted to create smoother surfaces.
  matern_change <- inla.spde2.pcmatern(mesh_spatial,
                                       prior.range = c(500,NA), # 500, NA # 10% chance range is smaller than 500000
                                       prior.sigma = c(0.2,0.1),   # 0.1, 0.1# 10% chance sd is larger than 1,
                                       constr = TRUE
  )
  
  # ----------------------------------------------------
  # iid random effect for 20-km blocks
  # ----------------------------------------------------
  pc_prec <- list(prior = "pcprec", param = c(0.1, 0.1))
  species_hex$hexid
  species_hex$hex_atlas <- as.numeric(as.factor(paste0(species_hex$hexid,"-",species_hex$Atlas)))
  
  # ----------------------------------------------------
  # Create mesh to model effect of time since sunrise (TSS)
  # Hours since sunrise is fit as a GAM-type effect.  
  # Recommend leaving these settings at these defaults
  # ----------------------------------------------------
  
  TSS_range <- range(species_hex$Hours_Since_Sunrise)
  TSS_meshpoints <- seq(TSS_range[1]-2,TSS_range[2]+2,length.out = 11)
  TSS_mesh1D = inla.mesh.1d(TSS_meshpoints,boundary="free")
  TSS_spde = inla.spde2.pcmatern(TSS_mesh1D,
                                 prior.range = c(6,0.1), # 10% range is smaller than 4
                                 prior.sigma = c(1,0.1),
                                 constr = TRUE) # 10% chance sd is larger than 1
  
  # ***************************************************************
  # Model formulas
  # ***************************************************************
  
  # Names of covariates to include in model. Ideally, covariates should be uncorrelated with each other.
  # Each covariate will include a quadratic effect to allow for 'intermediate optimum' effects
  covariates_to_include <- paste0("PC",1:5) 
  
  # How much shrinkage should be applied to covariate effects?
  sd_linear <- 0.2  # Change to smaller value (e.g., 0.1), if you want to heavily shrink covariate effects and potentially create smoother surfaces
  prec_linear <-  c(1/sd_linear^2,1/(sd_linear/2)^2)
  
  #TSS(main = Hours_Since_Sunrise,model = TSS_spde) +
  model_components = as.formula(paste0('~
            Intercept_OBBA2(1,model="linear",mean.linear = -10)+
            Intercept_OBBA3(1,model="linear", mean.linear = 0, prec.linear = 100)+
            range_effect_OBBA2(1,model="linear", mean.linear = -0.046, prec.linear = 10000)+
            kappa(hex_atlas, model = "iid", constr = TRUE, hyper = list(prec = pc_prec)) +
            spde_OBBA2(main = geometry, model = matern_abund) +
            spde_change(main = geometry, model = matern_change) +',
                                       paste0("Beta1_",covariates_to_include,'(1,model="linear", mean.linear = 0, prec.linear = ', prec_linear[1],')', collapse = " + "),
                                       '+',
                                       paste0("Beta2_",covariates_to_include,'(1,model="linear", mean.linear = 0, prec.linear = ', prec_linear[2],')', collapse = " + ")))
  
  #TSS +
  model_formula_OBBA2 = as.formula(paste0('count ~
                  Intercept_OBBA2 +
                  log_offset + 
                  kappa + 
                  range_effect_OBBA2 * distance_from_range +
                  spde_OBBA2 +',
                                          paste0("Beta1_",covariates_to_include,'*',covariates_to_include, collapse = " + "),
                                          '+',
                                          paste0("Beta2_",covariates_to_include,'*',covariates_to_include,"^2", collapse = " + ")))
  
  #TSS +
  model_formula_OBBA3 = as.formula(paste0('count ~
                  Intercept_OBBA2 +
                  Intercept_OBBA3 +
                  log_offset + 
                  kappa + 
                  range_effect_OBBA2 * distance_from_range +
                  spde_OBBA2 + spde_change + ',
                                          paste0("Beta1_",covariates_to_include,'*',covariates_to_include, collapse = " + "),
                                          '+',
                                          paste0("Beta2_",covariates_to_include,'*',covariates_to_include,"^2", collapse = " + ")))
  
  # ----------------------------------------------------------------------------
  # Fit model to point counts and ARUs with negative binomial error
  # ----------------------------------------------------------------------------
  start <- Sys.time()
  
  fit_INLA <- NULL
  while(is.null(fit_INLA)){
    
    fit_model <- function(){
      tryCatch(expr = {bru(components = model_components,
                           
                           like(family = "poisson",
                                formula = model_formula_OBBA2,
                                data = subset(species_hex, year(Date_Time)%in% c(2001:2005))),
                           
                           like(family = "poisson",
                                formula = model_formula_OBBA3,
                                data = subset(species_hex, year(Date_Time)%in% c(2021:2025))),
                           
                           
                           options = list(inla.mode = "experimental",
                                          bru_initial = list(lsig = 2,
                                                             Intercept_OBBA2 = -10,
                                                             Intercept_OBBA3 = 0,
                                                             Beta1_PC1 = 0,
                                                             Beta1_PC2 = 0,
                                                             Beta1_PC3 = 0,
                                                             Beta1_PC4 = 0,
                                                             Beta1_PC5 = 0,
                                                             
                                                             Beta2_PC1 = 0,
                                                             Beta2_PC2 = 0,
                                                             Beta2_PC3 = 0,
                                                             Beta2_PC4 = 0,
                                                             Beta2_PC5 = 0),
                                          control.compute = list(waic = FALSE, cpo = FALSE),
                                          bru_verbose = 4))},
               error = function(e){NULL})
    }
    fit_INLA <- fit_model()
    
    if ("try-error" %in% class(fit_INLA)) fit_INLA <- NULL
  }
  
  end <- Sys.time()
  runtime_INLA <- difftime( end,start, units="mins") %>% round(2)
  print(paste0(species_name," - ",runtime_INLA," min to fit model"))  # ~26 min
  
  # ****************************************************************************
  # ****************************************************************************
  # GENERATE MAPS
  # ****************************************************************************
  # ****************************************************************************
  
  # For every pixel on landscape, extract distance (in km) from eBird range limit
  ONGrid_species <- ONGrid %>%
    st_centroid() %>%
    st_transform(st_crs(species_hex)) %>%
    mutate(distance_from_range = (st_centroid(.) %>% 
                                    st_distance( . , range) %>% 
                                    as.numeric()))
  
  # ----------------------------------------------------
  # Generate predictions on ONGrid_species raster
  # ----------------------------------------------------
  
  pred_formula_OBBA2 = as.formula(paste0(' ~
                  Intercept_OBBA2 +
                  log_offset_5min +
                  range_effect_OBBA2 * distance_from_range +
                  spde_OBBA2 +',
                                         paste0("Beta1_",covariates_to_include,'*',covariates_to_include, collapse = " + "),
                                         '+',
                                         paste0("Beta2_",covariates_to_include,'*',covariates_to_include,"^2", collapse = " + ")))
  
  pred_formula_OBBA3 = as.formula(paste0(' ~
                  Intercept_OBBA2 +
                  Intercept_OBBA3 +
                  log_offset_5min +
                  range_effect_OBBA2 * distance_from_range +
                  spde_OBBA2 + spde_change +',
                                         paste0("Beta1_",covariates_to_include,'*',covariates_to_include, collapse = " + "),
                                         '+',
                                         paste0("Beta2_",covariates_to_include,'*',covariates_to_include,"^2", collapse = " + ")))
  
  pred_formula_change = as.formula(paste0(' ~
                  exp(Intercept_OBBA2 +
                  Intercept_OBBA3 +
                  log_offset_5min +
                  range_effect_OBBA2 * distance_from_range +
                  spde_OBBA2 + spde_change +',
                                          paste0("Beta1_",covariates_to_include,'*',covariates_to_include, collapse = " + "),
                                          '+',
                                          paste0("Beta2_",covariates_to_include,'*',covariates_to_include,"^2", collapse = " + "),
                                          ') - 
                                         
                  exp(Intercept_OBBA2 +
                  log_offset_5min +
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
                         ONGrid_species,
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
                         ONGrid_species,
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
  
  pred_change <- pred_OBBA3 - pred_OBBA2
  
  change_prediction_quantiles = apply(pred_change,1,function(x) quantile(x,c(0.05,0.5,0.95),na.rm = TRUE))
  ONGrid_species$change_pred_q05 <- change_prediction_quantiles[1,]
  ONGrid_species$change_pred_q50 <- change_prediction_quantiles[2,]
  ONGrid_species$change_pred_q95 <- change_prediction_quantiles[3,]
  ONGrid_species$change_pred_CI_width_90 <- change_prediction_quantiles[3,] - change_prediction_quantiles[1,]
  
  end2 <- Sys.time() 
  runtime_pred <- difftime( end2,start2, units="mins") %>% round(2)
  print(paste0(species_name," - ",runtime_pred," min to generate predictions")) # 14 min
  
  # ------------------------------------------------
  # sf object for ebird range limit (optional - not needed for plotting, but potentially helpful)
  # ------------------------------------------------
  
  if (!is.na(range)) range <- range  %>%
    st_transform(st_crs(ONBoundary)) %>%
    st_intersection(ONBoundary)
  
  # ****************************************************************************
  # RELATIVE ABUNDANCE IN OBBA2
  # ****************************************************************************
  
  # Convert to CRS of target raster
  tmp <- ONGrid_species %>% st_transform(st_crs(target_raster)) %>% st_centroid(.)
  
  colscale_relabund <-c("#FEFEFE", "#FBF7E2", "#FCF8D0", "#EEF7C2", "#CEF2B0", "#94E5A0", "#51C987", "#18A065", "#008C59", "#007F53", "#006344")
  colpal_relabund <- colorRampPalette(colscale_relabund)
  
  upper_bound <- quantile(ONGrid_species$OBBA2_pred_q50,0.95,na.rm = TRUE) %>% signif(2)
  lower_bound <- upper_bound/100
  
  tmp$OBBA2_pred_q50[tmp$OBBA2_pred_q50 > upper_bound] <-  upper_bound
  tmp$OBBA2_pred_q50[tmp$OBBA2_pred_q50 < lower_bound] <-  0
  
  OBBA2_plot_q50 <- ggplot() +
    
    geom_sf(data = tmp, aes(col = OBBA2_pred_q50), size = 0.1) +
    scale_color_gradientn(name = "Per point count",
                          colors = colpal_relabund(10),
                          trans = "log10",
                          na.value = "white")+
    
    
    geom_sf(data = ONBoundary,colour="black",fill=NA,lwd=0.3,show.legend = F) +
    
    geom_sf(data = range, linetype = 2, col = "gray50", fill = "transparent", size = 2)+
    
    coord_sf(clip = "off",xlim = range(as.data.frame(st_coordinates(ONBoundary))$X)) +
    theme(panel.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    # theme(legend.margin=margin(0,0,0,0),legend.box.margin=margin(5,10,5,-20),legend.title.align=0.5,
    #       legend.title = element_markdown(lineheight=.9,hjust = "left"))
    annotate(geom="text",x=400,y=1650, label= paste0(species_name),lineheight = .85,hjust = 0,size=6,fontface =2) +
    annotate(geom="text",x=400,y=1550, label= "OBBA 2",lineheight = .85,hjust = 0,size=4,fontface =2) +
    annotate(geom="text",x=400,y=1480, label= paste0("Prepared on ",Sys.Date()),size=2,lineheight = .75,hjust = 0,color="gray75")
  
  #print(OBBA2_plot_q50)
  
  png(paste0("../Output/Prediction_Maps/Relative_Abundance/",species_name,"_OBBA2_q50.png"), width=6.5, height=8, units="in", res=300, type="cairo")
  print(OBBA2_plot_q50)
  dev.off()
  
  # ****************************************************************************
  # RELATIVE ABUNDANCE IN OBBA3
  # ****************************************************************************
  
  tmp$OBBA3_pred_q50[tmp$OBBA3_pred_q50 > upper_bound] <-  upper_bound
  tmp$OBBA3_pred_q50[tmp$OBBA3_pred_q50 < lower_bound] <-  0
  
  OBBA3_plot_q50 <- ggplot() +
    
    geom_sf(data = tmp, aes(col = OBBA3_pred_q50), size = 0.1) +
    scale_color_gradientn(name = "Per point count",
                          colors = colpal_relabund(10),
                          trans = "log10",
                          na.value = "white")+
    
    
    geom_sf(data = ONBoundary,colour="black",fill=NA,lwd=0.3,show.legend = F) +
    geom_sf(data = range, linetype = 2, col = "gray50", fill = "transparent", size = 2)+
    
    coord_sf(clip = "off",xlim = range(as.data.frame(st_coordinates(ONBoundary))$X)) +
    theme(panel.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    # theme(legend.margin=margin(0,0,0,0),legend.box.margin=margin(5,10,5,-20),legend.title.align=0.5,
    #       legend.title = element_markdown(lineheight=.9,hjust = "left"))
    annotate(geom="text",x=400,y=1650, label= paste0(species_name),lineheight = .85,hjust = 0,size=6,fontface =2) +
    annotate(geom="text",x=400,y=1550, label= "OBBA 3",lineheight = .85,hjust = 0,size=4,fontface =2) +
    annotate(geom="text",x=400,y=1480, label= paste0("Prepared on ",Sys.Date()),size=2,lineheight = .75,hjust = 0,color="gray75")
  
  #   
  # print(OBBA3_plot_q50)
  
  png(paste0("../Output/Prediction_Maps/Relative_Abundance/",species_name,"_OBBA3_q50.png"), width=6.5, height=8, units="in", res=300, type="cairo")
  print(OBBA3_plot_q50)
  dev.off()
  
  # ****************************************************************************
  # CHANGE BETWEEN ATLASES
  # ****************************************************************************
  convert_numeric_to_char <- function(numeric_vector) {
    char_vector <- ifelse(numeric_vector > 0, paste0("+", numeric_vector), as.character(numeric_vector))
    return(char_vector)
  }
  
  colscale_change <-c("darkred","orangered","white","dodgerblue","darkblue")
  colpal_change <- colorRampPalette(colscale_change)
  
  tmp$change_pred_q50[tmp$change_pred_q50 > upper_bound] <-  upper_bound
  tmp$change_pred_q50[tmp$change_pred_q50 < -upper_bound] <-  -upper_bound
  
  percent_change <- 100*(apply(pred_OBBA3,2,sum,na.rm = TRUE) - apply(pred_OBBA2,2,sum,na.rm = TRUE))/apply(pred_OBBA2,2,sum,na.rm = TRUE)
  percent_change_quantiles <- quantile(percent_change,c(0.05,0.5,0.95)) %>% round() %>% convert_numeric_to_char()
  
  change_label <- paste0("Overall change = ",percent_change_quantiles[2],"% (",percent_change_quantiles[1]," to ",percent_change_quantiles[3],")")
  lim <- max(abs(tmp$change_pred_q50),na.rm = TRUE)
  change_plot_q50 <- ggplot() +
    
    geom_sf(data = tmp, aes(col = change_pred_q50), size = 0.1) +
    scale_color_gradientn(name = "Absolute\nchange\n",
                          colors = colpal_change(5),
                          na.value = "black",
                          limits = c(-lim,lim))+
    
    
    geom_sf(data = ONBoundary,colour="black",fill=NA,lwd=0.3,show.legend = F) +
    
    coord_sf(clip = "off",xlim = range(as.data.frame(st_coordinates(ONBoundary))$X)) +
    theme(panel.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    # theme(legend.margin=margin(0,0,0,0),legend.box.margin=margin(5,10,5,-20),legend.title.align=0.5,
    #       legend.title = element_markdown(lineheight=.9,hjust = "left"))
    annotate(geom="text",x=400,y=1650, label= paste0(species_name),lineheight = .85,hjust = 0,size=6,fontface =2) +
    annotate(geom="text",x=400,y=1550, label= change_label,lineheight = .85,hjust = 0,size=4,fontface =2) +
    annotate(geom="text",x=400,y=1480, label= paste0("Prepared on ",Sys.Date()),size=2,lineheight = .75,hjust = 0,color="gray75")+
    annotation_scale()
  # print(change_plot_q50)
  
  png(paste0("../Output/Prediction_Maps/Relative_Abundance/",species_name,"_change_q50.png"), width=6.5, height=8, units="in", res=300, type="cairo")
  print(change_plot_q50)
  dev.off()
  
  print(summary(fit_INLA))
  
} # close species loop