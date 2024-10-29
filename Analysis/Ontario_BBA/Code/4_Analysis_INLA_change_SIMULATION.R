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
require(ebirdst)

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

# Add 50 km buffer so that covariates can extend slightly outside province if possible
ONBoundary_buffer <- ONBoundary %>% st_buffer(50000)


# ONGrid <- ONGrid %>% st_intersection(ONBoundary)
ONGrid <- ONGrid %>% rename(geometry = x)
ONGrid_centroids <- st_centroid(ONGrid)

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
                         c("Canada Jay")) %>%
  unique()

#species_list <- unique(species_to_model)[1:20,]

# Prepare QPAD offsets for survey effort
require(napops)
napops_species <- list_species() %>% rename(Species_Code_NAPOPS = Species,Common_Name_NAPOPS = Common_Name,Scientific_Name_NAPOPS = Scientific_Name)

species_list <- left_join(species_list,napops_species[,c("Species_Code_NAPOPS","Common_Name_NAPOPS","Removal","Distance")], 
                          by = c("english_name" = "Common_Name_NAPOPS")) %>%
  mutate(EDR = NA, cue_rate = NA)

for (i in 1:nrow(species_list)){
  
  start <- Sys.time()
  
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
  
  if (species_list$Removal[i] == 1 & species_list$Distance[i] == 1){
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
  
  # **********************************************************
  # **********************************************************
  # 
  # Simulate data
  # 
  # **********************************************************
  # **********************************************************
  
  # ----------------------------------------------------
  # Simulate change across the landscape
  # ----------------------------------------------------
  
  mean_change <- -0.3
  
  inla.seed = round(runif(1,1,1000))
  spde_sim = inla.spde2.pcmatern(mesh_spatial, prior.range = c(.5, .5), prior.sigma = c(.5, .5))
  
  sim_range <- 200000 # correlation range of spatial field
  sim_sd <- 0.5  # SD of spatial field
  
  # Example from https://haakonbakkagit.github.io/btopic122.html
  Qu = inla.spde.precision(spde_sim, theta=c(log(sim_range), log(sim_sd)))
  u = inla.qsample(n=1, Q=Qu, seed = inla.seed)
  u = u[ ,1]
  
  # Change measured at every grid location
  A = inla.spde.make.A(mesh=mesh_spatial, loc=ONGrid_centroids)
  ONGrid_centroids$log_change = drop(A %*% u) + mean_change
  
  # ----------------------------------------------------
  # Extract eBird raster as "truth" in OBBA-2
  # ----------------------------------------------------
  
  ebirdst_download_status(species_name,pattern = "_mean_9km_")
  
  ebirdSDM <- load_raster(species_name, product = "abundance", 
                          period = "seasonal", metric = "mean", 
                          resolution = "9km")
  
  if ("resident" %in% names(ebirdSDM)){
    ebirdSDM <- ebirdSDM$resident
  } else{
    ebirdSDM <- ebirdSDM$breeding
  }
  
  values(ebirdSDM)[is.na(values(ebirdSDM))] <- 0
  
  # ----------------------------------------------------
  # Expected counts during each atlas cycle
  # ----------------------------------------------------
  
  ONGrid_centroids$OBBA2_truth <- extract(ebirdSDM,vect(ONGrid_centroids %>% st_transform(crs(ebirdSDM))))[,2]
  ONGrid_centroids$OBBA3_truth <- ONGrid_centroids$OBBA2_truth * exp(ONGrid_centroids$log_change)
  # 
  # # Plot surfaces
  # lim <- max(as.data.frame(ONGrid_centroids)[,c("OBBA2_truth","OBBA3_truth")],na.rm = TRUE)
  # ggplot(sample_n(ONGrid_centroids,30000))+geom_sf(data = ONBoundary_buffer) + geom_sf(aes(col = OBBA2_truth), size = 2) + scale_color_gradientn(colors = viridis(10),lim=c(0,lim))
  # ggplot(sample_n(ONGrid_centroids,30000))+geom_sf(data = ONBoundary_buffer) + geom_sf(aes(col = OBBA3_truth), size = 2) + scale_color_gradientn(colors = viridis(10),lim=c(0,lim))
  # 
  # # lim <- max(abs(as.data.frame(ONGrid_centroids)[,c("log_change")]),na.rm = TRUE)
  # # ggplot(sample_n(ONGrid_centroids,30000))+geom_sf(data = ONBoundary_buffer) + geom_sf(aes(col = log_change), size = 2) + scale_color_gradientn(colors = inferno(10),lim=c(-lim,lim))
  # 
  # 
  N2 <- sum(ONGrid_centroids$OBBA2_truth,na.rm = TRUE)
  N3 <- sum(ONGrid_centroids$OBBA3_truth,na.rm = TRUE)
  N3/N2
  
  # ----------------------------------------------------
  # Simulate counts and change at survey locations
  # ----------------------------------------------------
  
  # Change measured at every survey location
  A = inla.spde.make.A(mesh=mesh_spatial, loc=sp_dat)
  sp_dat$log_change = drop(A %*% u) + mean_change
  
  sp_dat$OBBA2_truth <- extract(ebirdSDM,vect(sp_dat %>% st_transform(crs(ebirdSDM))))[,2]
  sp_dat$OBBA3_truth <- sp_dat$OBBA2_truth * exp(sp_dat$log_change)
  
  # # Plot surfaces
  # lim <- max(abs(as.data.frame(ONGrid_centroids)[,c("log_change")]),na.rm = TRUE)
  # ggplot()+
  #   geom_sf(data = ONBoundary_buffer) + 
  #   geom_sf(data = sample_n(ONGrid_centroids,30000), aes(col = log_change), size = 3) + 
  #   #geom_sf(data = sp_dat, aes(col = log_change), size = 2) + 
  #   scale_color_gradientn(colors = inferno(10),lim=c(-lim,lim))

  sp_dat$count_OBBA2 <- rpois(nrow(sp_dat),sp_dat$OBBA2_truth)
  sp_dat$count_OBBA3 <- rpois(nrow(sp_dat),sp_dat$OBBA3_truth)
  
  sp_dat$count <- sp_dat$count_OBBA2
  sp_dat$count[year(sp_dat$Date_Time) >= 2020] <- sp_dat$count_OBBA3[year(sp_dat$Date_Time) >= 2020]
  
  # **********************************************************
  # **********************************************************
  # 
  # Conduct analysis
  # 
  # **********************************************************
  # **********************************************************
  
  # ----------------------------------------------------
  # Prepare analysis
  # ----------------------------------------------------
  
  
  # Controls the 'residual spatial field'.  This can be adjusted to create smoother surfaces.
  matern_abund <- inla.spde2.pcmatern(mesh_spatial,
                                      prior.range = c(500000, 0.1), # 10% chance range is smaller than 300000
                                      prior.sigma = c(0.5,0.1),    # 10% chance sd is larger than 0.5,
                                      constr = TRUE
  )
  
  # Controls the 'residual spatial field'.  This can be adjusted to create smoother surfaces.
  matern_change <- inla.spde2.pcmatern(mesh_spatial,
                                       prior.range = c(500000, 0.1), # 10% chance range is smaller than 300000
                                       prior.sigma = c(0.1,0.1),    # 10% chance sd is larger than 0.1,
                                       constr = TRUE
  )
  
  # ----------------------------------------------------
  # Create mesh to model effect of time since sunrise (TSS)
  # Hours since sunrise is fit as a GAM-type effect.  
  # Recommend leaving these settings at these defaults
  # ----------------------------------------------------
   
  TSS_range <- range(sp_dat$Hours_Since_Sunrise)
  TSS_meshpoints <- seq(TSS_range[1]-2,TSS_range[2]+2,length.out = 11)
  TSS_mesh1D = inla.mesh.1d(TSS_meshpoints,boundary="free")
  TSS_spde = inla.spde2.pcmatern(TSS_mesh1D,
                                 prior.range = c(6,0.1), # 10% range is smaller than 4
                                 prior.sigma = c(1,0.1)) # 10% chance sd is larger than 1
  
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
            spde_change(main = coordinates, model = matern_change) +
            ',
                                       paste0("Beta1_",covariates_to_include,'(1,model="linear", mean.linear = 0, prec.linear = ', prec_linear[1],')', collapse = " + "),
                                       '+',
                                       paste0("Beta2_",covariates_to_include,'(1,model="linear", mean.linear = 0, prec.linear = ', prec_linear[2],')', collapse = " + "),
                                       '+',
                                       paste0("Beta1_change_",covariates_to_include,'(1,model="linear", mean.linear = 0, prec.linear = 2500)', collapse = " + "),
                                       '+',
                                       paste0("Beta2_change_",covariates_to_include,'(1,model="linear", mean.linear = 0, prec.linear = 2500)', collapse = " + ")))
  
  #TSS +
  model_formula_OBBA2 = as.formula(paste0('count ~
                  Intercept_OBBA2 +
                  range_effect_OBBA2 * distance_from_range +
                  spde_OBBA2 +',
                                          paste0("Beta1_",covariates_to_include,'*',covariates_to_include, collapse = " + "),
                                          '+',
                                          paste0("Beta2_",covariates_to_include,'*',covariates_to_include,"^2", collapse = " + ")))
  
  #TSS +
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
  
  start <- Sys.time()
  fit_INLA <- NULL
  while(is.null(fit_INLA)){
    
    fit_model <- function(){
      tryCatch(expr = {bru(components = model_components,
                           
                           like(family = "poisson",
                                formula = model_formula_OBBA2,
                                data = as(subset(sp_dat, year(Date_Time)%in% c(2001:2005)),'Spatial')),
                           
                            like(family = "poisson",
                                 formula = model_formula_OBBA3,
                                 data = as(subset(sp_dat, year(Date_Time)%in% c(2021:2025)),'Spatial')),
                           
                           options = list(bru_initial = list(lsig = 2,
                                                             Intercept_OBBA2 = -10),
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
  ONGrid <- ONGrid %>%
    mutate(distance_from_range = (st_centroid(.) %>% 
                                    st_distance( . , range) %>% 
                                    as.numeric())/1000)
  
  # ----------------------------------------------------
  # Generate predictions on ONGrid raster
  # ----------------------------------------------------
  
  pred_formula_OBBA2 = as.formula(paste0(' ~
                  Intercept_OBBA2 +
                  range_effect_OBBA2 * distance_from_range +
                  spde_OBBA2 +',
                                         paste0("Beta1_",covariates_to_include,'*',covariates_to_include, collapse = " + "),
                                         '+',
                                         paste0("Beta2_",covariates_to_include,'*',covariates_to_include,"^2", collapse = " + ")))
  
  pred_formula_OBBA3 = as.formula(paste0(' ~
                  Intercept_OBBA2 +
                  Intercept_OBBA3 +
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
                         as(ONGrid,'Spatial'),
                         formula =  pred_formula_OBBA2,
                         n.samples = 500,
                         seed = 123)
  
  pred_OBBA2 <- exp(pred_OBBA2)
  OBBA2_prediction_quantiles = apply(pred_OBBA2,1,function(x) quantile(x,c(0.05,0.5,0.95),na.rm = TRUE))
  ONGrid$OBBA2_pred_q05 <- OBBA2_prediction_quantiles[1,]
  ONGrid$OBBA2_pred_q50 <- OBBA2_prediction_quantiles[2,]
  ONGrid$OBBA2_pred_q95 <- OBBA2_prediction_quantiles[3,]
  ONGrid$OBBA2_pred_CI_width_90 <- OBBA2_prediction_quantiles[3,] - OBBA2_prediction_quantiles[1,]
  
  # --------
  # OBBA-3
  # --------

  pred_OBBA3 <- NULL
  pred_OBBA3 <- generate(fit_INLA,
                         as(ONGrid,'Spatial'),
                         formula =  pred_formula_OBBA3,
                         n.samples = 500,
                         seed = 123)

  pred_OBBA3 <- exp(pred_OBBA3)
  OBBA3_prediction_quantiles = apply(pred_OBBA3,1,function(x) quantile(x,c(0.05,0.5,0.95),na.rm = TRUE))
  ONGrid$OBBA3_pred_q05 <- OBBA3_prediction_quantiles[1,]
  ONGrid$OBBA3_pred_q50 <- OBBA3_prediction_quantiles[2,]
  ONGrid$OBBA3_pred_q95 <- OBBA3_prediction_quantiles[3,]
  ONGrid$OBBA3_pred_CI_width_90 <- OBBA3_prediction_quantiles[3,] - OBBA3_prediction_quantiles[1,]

  # --------
  # Change
  # --------

  pred_change <- NULL
  pred_change <- generate(fit_INLA,
                          as(ONGrid,'Spatial'),
                          formula =  pred_formula_change,
                          n.samples = 500,
                          seed = 123)

  change_prediction_quantiles = apply(pred_change,1,function(x) quantile(x,c(0.05,0.5,0.95),na.rm = TRUE))
  ONGrid$change_pred_q05 <- change_prediction_quantiles[1,]
  ONGrid$change_pred_q50 <- change_prediction_quantiles[2,]
  ONGrid$change_pred_q95 <- change_prediction_quantiles[3,]
  ONGrid$change_pred_CI_width_90 <- change_prediction_quantiles[3,] - change_prediction_quantiles[1,]

  end2 <- Sys.time()
  runtime_pred <- difftime( end2,start2, units="mins") %>% round(2)
  print(paste0(species_name," - ",runtime_pred," min to generate predictions")) # 12 min

  # ****************************************************************************
  # True versus predicted relative abundance surface
  # ****************************************************************************
  
  ONGrid_centroids$OBBA2_pred <- ONGrid$OBBA2_pred_q50
 
  # Plot abundance surface
  lim <- max(abs(as.data.frame(ONGrid_centroids)[,c("OBBA2_pred","OBBA2_truth")]),na.rm = TRUE)
  
  ggplot()+
     geom_sf(data = ONBoundary_buffer) + 
     geom_sf(data = sample_n(ONGrid_centroids,30000), aes(col = OBBA2_truth), size = 3) + 
     #geom_sf(data = sp_dat, aes(col = log_change), size = 2) + 
     scale_color_gradientn(colors = viridis(10),lim=c(0,lim))
  
  ggplot()+
    geom_sf(data = ONBoundary_buffer) + 
    geom_sf(data = sample_n(ONGrid_centroids,30000), aes(col = OBBA2_pred), size = 3) + 
    #geom_sf(data = sp_dat, aes(col = log_change), size = 2) + 
    scale_color_gradientn(colors = viridis(10),lim=c(0,lim))
  
  ONGrid_centroids$delta_pred <- ONGrid$change_pred_q50
  ONGrid_centroids$delta_truth <- ONGrid_centroids$OBBA3_truth - ONGrid_centroids$OBBA2_truth
  
  # Plot change surface
  lim <- max(abs(as.data.frame(ONGrid_centroids)[,c("delta_pred","delta_truth")]),na.rm = TRUE)
  
  ggplot()+
    geom_sf(data = ONBoundary_buffer) + 
    geom_sf(data = sample_n(ONGrid_centroids,30000), aes(col = delta_truth), size = 3) + 
    #geom_sf(data = sp_dat, aes(col = log_change), size = 2) + 
    scale_color_gradientn(colors = c("orangered","white","dodgerblue"),lim=c(-lim,lim))
  
  ggplot()+
    geom_sf(data = ONBoundary_buffer) + 
    geom_sf(data = sample_n(ONGrid_centroids,30000), aes(col = delta_pred), size = 3) + 
    #geom_sf(data = sp_dat, aes(col = log_change), size = 2) + 
    scale_color_gradientn(colors = c("orangered","white","dodgerblue"),lim=c(-lim,lim))
  
  
  
  # Estimate of total change
  total_change <- pred_change %>% apply(2,sum,na.rm = TRUE) # -71076510.5    264621.4
  xlim <- max(abs(total_change)) %>% ceiling()
  hist(total_change, breaks = seq(-xlim,xlim,length.out = 100))
  abline(v = (N3-N2),col="blue", lwd = 2)
  
  
} # close species loop


