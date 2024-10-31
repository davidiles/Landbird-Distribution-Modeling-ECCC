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
require(ebirdst)

# Timeout INLA after 10 minutes (if it has not fit by then, it has likely stalled)
inla.setOption(inla.timeout = 60*10)

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
  st_transform(AEA_proj) %>%
  st_buffer(5) %>%
  st_buffer(-5) %>%
  st_union() %>%
  st_cast("POLYGON") %>%
  st_as_sf() %>%
  mutate(Area_km2 = as.numeric(st_area(.))) %>%
  subset(Area_km2 >= 100)

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

hexgrid <- st_make_grid(ONBoundary, cellsize = 25, square=FALSE, what = "polygons") %>%
  st_as_sf() %>% mutate(hexid = 1:nrow(.)) %>% dplyr::rename(geometry = x) %>%
  st_intersection(ONBoundary)

hex_centroids <- st_centroid(hexgrid)

dim(hexgrid)

# ------------------------------------------
# Select Point Counts / ARUs to use
# Decision to remove "SPM" aru transcriptions, as they are considered less reliable
# ------------------------------------------

all_surveys$Hours_Since_Sunrise <- as.numeric(all_surveys$Hours_Since_Sunrise)
all_surveys$Obs_Index <- 1:nrow(all_surveys)

all_surveys <- subset(all_surveys,
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
                                          "Ontario Breeding Bird Atlas Digital Point Counts 2021")) %>%
                      arrange(Obs_Index)

dim(all_surveys) # 101206

# ------------------------------------------
# Subset
# ------------------------------------------

full_count_matrix <- full_count_matrix[all_surveys$Obs_Index,]

# ******************************************************************
# Loop through species, fit models, generate maps
# ******************************************************************

# species_summary <- species_summary %>% subset(`PC/ARU` >=20)

set.seed(123)
species_list <- subset(species_to_model, n_detections > 500) %>%
  sample_n(50)

species_list <- rbind(subset(species_to_model,english_name %in% 
                               #c("Least Bittern","Short-billed Dowitcher"))) %>% unique()
                               c("Canada Jay","Blackpoll Warbler","Lesser Yellowlegs","Boreal Chickadee","Black-capped Chickadee","Rusty Blackbird","Tennessee Warbler")),
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
  # Extract/process ebird range for this species (if it exists)
  # ----------------------------------------------------
  
  check <- get_species(species_name) # does the species have a range estimated from eBird?
  
  if (length(check)==0 | (length(check)>0 & is.na(check))){
    
    print(paste0(species_name, " not available in eBird"))
    sp_dat$distance_from_range <- 0
    range <- NA
    
  }else{
    ebirdst_download_status(species_name, 
                            download_abundance = FALSE,
                            download_ranges = TRUE,
                            pattern = "_smooth_27km_") 
    
    path <- get_species_path(species_name)
    
    range <- load_ranges(species_name, resolution = "27km",smoothed = TRUE)
    
    range <- range %>% subset(season %in% c("resident","breeding")) %>% 
      st_transform(.,crs = st_crs(ONBoundary)) %>% 
      st_union()
    
    sp_dat$distance_from_range <- ((st_distance(sp_dat, range)) %>% as.numeric()) # should be in km
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
              effort = sum(Survey_Duration_Minutes)) %>%
    mutate(count_per_effort = count/effort)
  
  species_hex_centroid <- full_join(hex_centroids,species_hex_centroid) %>% na.omit()
  
  col_lim <- c(0,max(species_hex_centroid$count,na.rm = TRUE))
  
  species_data_plot <- ggplot() +
    geom_sf(data=ONBoundary,colour="gray90", fill = "white")+
    geom_sf(data=subset(species_hex_centroid,count==0), aes(size=effort), shape = 1, col = "gray80")+
    geom_sf(data=subset(species_hex_centroid,count>0), aes(col = count_per_effort, size = effort))+
    geom_sf(data = range %>% st_intersection(ONBoundary), col = "red", fill = "transparent")+
    annotation_scale(style = "ticks",
                     text_face = "bold")+
    theme_bw()+
    
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = 'transparent', colour = 'black'))+
    scale_color_gradientn(colors = c("black",viridis(10)),na.value=NA, trans = "log10")+
    scale_size_continuous(name = "Total Survey Effort (min)", range = c(1,4))+
    ggtitle(paste0(species_name," - Observed count per min"))+
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
                                      prior.range = c(500, 0.01), # 500, NA # 10% chance range is smaller than 500000
                                      prior.sigma = c(0.2,0.1),   # 0.1, 0.1 # 10% chance sd is larger than 1,
                                      constr = TRUE
  )
  
  # Controls the 'residual spatial field'.  This can be adjusted to create smoother surfaces.
  matern_change <- inla.spde2.pcmatern(mesh_spatial,
                                       prior.range = c(1000,0.01), # 500, NA # 10% chance range is smaller than 500000
                                       prior.sigma = c(0.2,0.1),   # 0.1, 0.1# 10% chance sd is larger than 1,
                                       constr = TRUE
  )
  
  # ----------------------------------------------------
  # iid random effect for 20-km blocks
  # ----------------------------------------------------
  pc_prec <- list(prior = "pcprec", param = c(0.1, 0.1))
  species_hex$hexid <- as.numeric(as.factor(species_hex$hexid))
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
  covariates_to_include <- paste0("PC",1:2) 
  
  # How much shrinkage should be applied to covariate effects?
  sd_linear <- 0.2  # Change to smaller value (e.g., 0.1), if you want to heavily shrink covariate effects and potentially create smoother surfaces
  prec_linear <-  c(1/sd_linear^2,1/(sd_linear/2)^2)
  
  #TSS(main = Hours_Since_Sunrise,model = TSS_spde) +
  model_components = as.formula(paste0('~
            Intercept_OBBA2(1,model="linear",mean.linear = -10)+
            Intercept_OBBA3(1,model="linear", mean.linear = 0, prec.linear = 100)+
            range_effect_OBBA2(1,model="linear", mean.linear = -0.046, prec.linear = 10000)+
            eta(hexid, model = "iid", constr = TRUE, hyper = list(prec = pc_prec)) +
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
                  eta +
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
                  eta +
                  kappa + 
                  range_effect_OBBA2 * distance_from_range +
                  spde_OBBA2 + spde_change + ',
                                          paste0("Beta1_",covariates_to_include,'*',covariates_to_include, collapse = " + "),
                                          '+',
                                          paste0("Beta2_",covariates_to_include,'*',covariates_to_include,"^2", collapse = " + ")))
  
  # ----------------------------------------------------------------------------
  # Fit model to point counts and ARUs with negative binomial error
  # ----------------------------------------------------------------------------
  
  error_type <- "poisson"
  if (sum(species_hex$count>4)>25 & sum(sp_dat$count>5)>10) error_type <- "nbinomial"
  
  start <- Sys.time()
  
  fit_INLA <- NULL
  while(is.null(fit_INLA)){
    
    fit_model <- function(){
      tryCatch(expr = {bru(components = model_components,
                           
                           like(family = error_type,
                                formula = model_formula_OBBA2,
                                data = subset(species_hex, year(Date_Time)%in% c(2001:2005))),
                           
                           like(family = error_type,
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
  print(paste0(species_name," - ",runtime_INLA," min to fit model"))  # ~8 min
  print(summary(fit_INLA))
  
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
  print(paste0(species_name," - ",runtime_pred," min to generate predictions")) # 7 min
  
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
  
  upper_bound <- quantile(ONGrid_species$OBBA2_pred_q50,0.99,na.rm = TRUE) %>% signif(2)
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
  
  tmp$trend_signif <- tmp$change_pred_q95 < 0 | tmp$change_pred_q05 > 0
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
    #geom_sf(data = subset(tmp, trend_signif == TRUE & change_pred_q50 < 0) %>% st_buffer(5) %>% st_union(), col = "red", fill = "red", alpha = 0.2)+
    #geom_sf(data = subset(tmp, trend_signif == TRUE & change_pred_q50 > 0) %>% st_buffer(5) %>% st_union(), col = "blue", fill = "blue", alpha = 0.2)+
    
    #geom_sf(data = subset(hex_for_trend, signif == TRUE & change_pred_q50 > 0) %>% st_union(), col = "blue", fill = "blue", alpha = 0.2)+
    
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
  
  
  # Combine all figures into a single one
  png(file = paste0("../Output/Prediction_Maps/Relative_Abundance/",species_name,"_full.png"), 
      units = "in", 
      width = 10, 
      height = 12, res = 600)
  grid.arrange(species_data_plot, OBBA2_plot_q50,OBBA3_plot_q50,change_plot_q50,
               ncol = 3,nrow=4, 
               layout_matrix = cbind(c(1,1,1,1), 
                                     c(2,2,3,3),
                                     c(NA,4,4,NA)))
  dev.off()
  
  print(summary(fit_INLA))
  
} # close species loop