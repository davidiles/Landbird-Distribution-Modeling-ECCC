# ************************************************
# BAYESIAN ANALYSIS / SPECIES DISTRIBUTION MODELS FOR ONTARIO BREEDING BIRD ATLAS
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
require(gridExtra)

# Timeout INLA after 20 minutes (if it has not fit by then, it has likely stalled)
inla.setOption(inla.timeout = 60*20)

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

ON_BCR <- st_read("../../../Data/Spatial/National/BCR/BCR_Terrestrial_master.shp")  %>%
  subset(PROVINCE_S == "ONTARIO") %>%
  st_make_valid() %>%
  group_by(BCR,BCRNAME) %>%
  summarize() %>%
  st_transform(AEA_proj)

# Smooth out the Ontario boundary a bit (mostly for visualization purposes)
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

# Identify the BCR in which each point is located (use st_intersects rather than st_intersection to save time)
# Note: ONGrid is a series of points 1 km apart at which predictions will be generated for mapping/summarizing

ONGrid <- ONGrid %>% 
  rename(geometry = x) %>% 
  st_centroid() %>%
  st_transform(AEA_proj) %>%
  mutate(BCR = NA)

for (bcr in unique(ON_BCR$BCR)){
  
  bcr_poly <- subset(ON_BCR, BCR == bcr) %>% st_transform(st_crs(ONGrid))
  
  tmp <- st_intersects(ONGrid,bcr_poly) %>% as.numeric()
  
  ONGrid$BCR[!is.na(tmp)] <- bcr
}

# # Raster with target properties (needed for plotting)
# target_raster <- rast("../../../Data/Spatial/National/AnnualMeanTemperature/wc2.1_30s_bio_1.tif") %>% 
#   crop(st_transform(ONBoundary,crs(.))) %>% 
#   project(AEA_proj, res = 1000) %>%
#   mask(vect(ONBoundary))

# ----------------------------------------------------
# Load shapefile containing atlas squares (10km x 10km)
# ----------------------------------------------------

squares <- st_read("../../../Data/Spatial/National/AtlasSquares/NationalSquares_FINAL.shp")  %>%
  subset(prov == "ON") %>%
  st_transform(AEA_proj)

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

set.seed(123)


# Prepare QPAD offsets for survey effort
require(napops)
napops_species <- list_species() %>% rename(Species_Code_NAPOPS = Species,Common_Name_NAPOPS = Common_Name,Scientific_Name_NAPOPS = Scientific_Name)

# ******************************************************************
# Loop through species, fit models, generate maps, store summary of total change overall and per BCR
# ******************************************************************

# Empty table to populate with summaries of change per BCR, and overall across the province
change_estimate_table <- data.frame()
if (file.exists(paste0("../Output/Tables_Summaries/OBBA3_change_estimate_table.csv"))) change_estimate_table <- read.csv(paste0("../Output/Tables_Summaries/OBBA3_change_estimate_table.csv"), header = TRUE)

# Empty table to populate with summaries of covariate effects
covariate_effect_table <- data.frame()
if (file.exists(paste0("../Output/Tables_Summaries/covariate_effect_table.csv"))) covariate_effect_table <- read.csv(paste0("../Output/Tables_Summaries/covariate_effect_table.csv"), header = TRUE)

species_to_model <- left_join(species_to_model,
                              napops_species[,c("Species_Code_NAPOPS","Common_Name_NAPOPS","Removal","Distance")], 
                              by = c("english_name" = "Common_Name_NAPOPS")) %>%
  mutate(EDR = NA, cue_rate = NA)

species_to_model <- subset(species_to_model, n_detections >= 600 & !is.na(Removal) & !is.na(Distance)) %>% unique()

species_list <- unique(species_to_model$english_name)

# Shuffle
species_list <- sample(species_list,length(species_list))

species_list <- c("Connecticut Warbler", "Canada Jay", "Boreal Chickadee", "Dark-eyed Junco")
for (spec in species_list){
  
  i <- which(species_to_model$english_name == spec)
  
  # ----------------------------------------------------
  # Extract counts/data for this species
  # ----------------------------------------------------
  
  sp_code = as.character(species_to_model$Species_Code_BSC[i])
  species_name = species_to_model$english_name[i]
  print(species_name)
  if (file.exists(paste0("../Output/Prediction_Maps/Relative_Abundance/",species_name,".png"))) next
  
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
    # ebirdst_download_status(species_name, 
    #                         download_abundance = FALSE,
    #                         download_ranges = TRUE,
    #                         pattern = "_smooth_27km_") 
    # 
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
  
  if (!is.na(species_to_model$Removal[i]) & !is.na(species_to_model$Distance[i]) &
      species_to_model$Removal[i] == 1 & species_to_model$Distance[i] == 1){
    species_to_model$cue_rate[i] <- cue_rate(species = species_to_model$Species_Code_NAPOPS[i],od = 153, tssr = 0, model = 1)[3] %>% as.numeric()
    species_to_model$EDR[i] <- edr(species = species_to_model$Species_Code_NAPOPS[i],road = FALSE, forest = 0.5,model = 1)[3] %>% as.numeric()
    
    # Calculate A and p, which jointly determine offset
    A_metres <- c(pi*species_to_model$EDR[i]^2)
    p <- 1-exp(-5*species_to_model$cue_rate[i])
    
    sp_dat$log_offset <- log(A_metres * (1-exp(-sp_dat$Survey_Duration_Minutes*species_to_model$cue_rate[i])))
    log_offset_5min <- log(A_metres * p)
  } else {
    sp_dat$log_offset <- 0
    log_offset_5min <- 0
  }
  
  # ----------------------------------------------------
  # Generate plot of observed data in each hexagon
  # ----------------------------------------------------
  
  # Intersect with survey information
  species_square <- sp_dat %>% st_intersection(squares)
  species_square$Atlas <- "OBBA2"
  species_square$Atlas[year(species_square$Date_Time) >= 2020] <- "OBBA3"
  
  # Summarize total effort and mean count in each squareagon
  species_square_centroid <- species_square %>%
    as.data.frame() %>%
    group_by(square_id,Atlas) %>%
    summarize(count = sum(count),
              effort = sum(Survey_Duration_Minutes)) %>%
    mutate(count_per_effort = count/effort)
  
  species_square_centroid <- full_join(st_centroid(squares) %>% select(square_id),species_square_centroid) %>% na.omit()
  
  col_lim <- c(min(species_square_centroid$count_per_effort[species_square_centroid$count_per_effort>0],na.rm = TRUE),max(species_square_centroid$count_per_effort,na.rm = TRUE))
  colscale_relabund <-c("#FEFEFE", "#FBF7E2", "#FCF8D0", "#EEF7C2", "#CEF2B0", "#94E5A0", "#51C987", "#18A065", "#008C59", "#007F53", "#006344")
  
  species_data_OBBA2 <- ggplot() +
    geom_sf(data=ONBoundary,colour="gray90", fill = "white")+
    geom_sf(data=subset(species_square_centroid,count==0 & Atlas == "OBBA2"), aes(size=effort), shape = 1, col = "gray80")+
    geom_sf(data=subset(species_square_centroid,count>0  & Atlas == "OBBA2"), aes(col = count_per_effort, size = effort))+
    geom_sf(data = range %>% st_intersection(ONBoundary), col = "gray50", fill = "transparent", linetype = 2)+
    geom_sf(data = ONBoundary,colour="black",fill=NA,lwd=0.3,show.legend = F) +
    
    coord_sf(clip = "off",xlim = range(as.data.frame(st_coordinates(ONBoundary))$X)) +
    theme(panel.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    
    annotate(geom="text",x=400,y=1650, label= "Observed Data OBBA2",lineheight = .85,hjust = 0,size=6,fontface =2) +
    #annotate(geom="text",x=400,y=1550, label= "OBBA 2 - Observed Data",lineheight = .85,hjust = 0,size=4,fontface =2) +
    annotate(geom="text",x=400,y=1480, label= paste0("Prepared on ",Sys.Date()),size=2,lineheight = .75,hjust = 0,color="gray75")+
    
    scale_color_gradientn(colors = colscale_relabund,na.value=NA, trans = "log10", name = "Count/Min", limits = col_lim)+
    scale_size_continuous(name = "Total Survey Effort (min)", range = c(1,4), guide = "none")
  
  #species_data_OBBA2
  
  species_data_OBBA3 <- ggplot() +
    geom_sf(data=ONBoundary,colour="gray90", fill = "white")+
    geom_sf(data=subset(species_square_centroid,count==0 & Atlas == "OBBA3"), aes(size=effort), shape = 1, col = "gray80")+
    geom_sf(data=subset(species_square_centroid,count>0  & Atlas == "OBBA3"), aes(col = count_per_effort, size = effort))+
    geom_sf(data = range %>% st_intersection(ONBoundary), col = "gray50", fill = "transparent", linetype = 2)+
    geom_sf(data = ONBoundary,colour="black",fill=NA,lwd=0.3,show.legend = F) +
    
    coord_sf(clip = "off",xlim = range(as.data.frame(st_coordinates(ONBoundary))$X)) +
    theme(panel.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    
    annotate(geom="text",x=400,y=1650, label= "Observed Data OBBA3",lineheight = .85,hjust = 0,size=6,fontface =2) +
    #annotate(geom="text",x=400,y=1550, label= "OBBA 3 - Observed Data",lineheight = .85,hjust = 0,size=4,fontface =2) +
    annotate(geom="text",x=400,y=1480, label= paste0("Prepared on ",Sys.Date()),size=2,lineheight = .75,hjust = 0,color="gray75")+
    
    scale_color_gradientn(colors = colscale_relabund,na.value=NA, trans = "log10", name = "Count/Min", limits = col_lim)+
    scale_size_continuous(name = "Total Survey Effort (min)", range = c(1,4), guide = "none")
  
  # species_data_OBBA3
  
  # ----------------------------------------------------
  # Only fit model if species was detected in at least 50 squares in both atlases
  # ----------------------------------------------------
  
  n_square_OBBA2 <- subset(species_square_centroid,Atlas=="OBBA2" & count > 0)
  n_square_OBBA3 <- subset(species_square_centroid,Atlas=="OBBA3" & count > 0)
  
  if (nrow(n_square_OBBA2) < 20 | nrow(n_square_OBBA3) < 20 ) next
  
  # ----------------------------------------------------
  # Create a spatial mesh, which is used to fit the residual spatial field
  # ----------------------------------------------------
  
  # make a two extension hulls and mesh for spatial model
  hull <- fm_extensions(
    ONBoundary,
    convex = c(50, 200),
    concave = c(350, 500)
  )
  
  mesh_spatial <- fm_mesh_2d_inla(
    boundary = hull, 
    max.edge = c(50, 200), # km inside and outside  
    cutoff = 10, 
    crs = fm_crs(species_square)
  )
  
  mesh_locs <- mesh_spatial$loc[,c(1,2)] %>% as.data.frame()
  dim(mesh_locs)
  
  # -------------------
  # Note: some hints in here about SPDE range / sd parameters
  #     https://tutorials.inbo.be/tutorials/r_inla/spatial.pdf
  # -------------------
  
  # Controls the 'residual spatial field'.  This can be adjusted to create smoother surfaces.
  matern_abund <- inla.spde2.pcmatern(mesh_spatial,
                                      prior.range = c(500,0.5),#c(500, 0.01), # 500, NA # 1% chance range is smaller than 500000
                                      prior.sigma = c(0.1,0.1),   # 10% chance sd is larger than 0.1,
                                      constr = TRUE
  )
  
  # Controls the 'residual spatial field'.  This can be adjusted to create smoother surfaces.
  matern_change <- inla.spde2.pcmatern(mesh_spatial,
                                       prior.range = c(500,0.5),#c(500,0.01), # 500, NA # 1% chance range is smaller than 500000
                                       prior.sigma = c(0.1,0.1),   # 10% chance sd is larger than 0.1,
                                       constr = TRUE
  )
  
  # ----------------------------------------------------
  # iid random effect for atlas squares
  # ----------------------------------------------------
  
  pc_prec <- list(prior = "pcprec", param = c(0.1, 0.1))
  species_square$square_id <- as.numeric(as.factor(species_square$square_id))
  species_square$square_atlas <- as.numeric(as.factor(paste0(species_square$square_id,"-",species_square$Atlas)))
  
  # ----------------------------------------------------
  # Create mesh to model effect of time since sunrise (HSS)
  # Hours since sunrise is fit as a GAM-type effect.  
  # Recommend leaving these settings at these defaults
  # ----------------------------------------------------
  
  HSS_range <- range(species_square$Hours_Since_Sunrise)
  HSS_meshpoints <- seq(HSS_range[1]-1,HSS_range[2]+1,length.out = 21)
  HSS_mesh1D = inla.mesh.1d(HSS_meshpoints,boundary="free")
  HSS_spde = inla.spde2.pcmatern(HSS_mesh1D,
                                 prior.range = c(3,0.1), # 10% chance range is smaller than 3
                                 prior.sigma = c(1,0.1), # 10% chance sd is larger than 1
                                 constr = TRUE) 
  
  # ***************************************************************
  # Model formulas
  # ***************************************************************
  
  # Names of covariates to include in model. Ideally, covariates should be uncorrelated with each other.
  # Each covariate will include a quadratic effect to allow for 'intermediate optimum' effects
  covariates_to_include <- paste0("PC",1:3) 
  
  # How much shrinkage should be applied to covariate effects?
  sd_linear <- 0.1  # Change to smaller value (e.g., 0.1), if you want to heavily shrink covariate effects and potentially create smoother surfaces
  prec_linear <-  c(1/sd_linear^2,1/(sd_linear/2)^2)
  
  
  # ----------------------------------------------------------------------------
  # Fit full model
  # ----------------------------------------------------------------------------
  
  model_components = as.formula(paste0('~
            Intercept_OBBA2(1,model="linear",mean.linear = -10)+
            Intercept_OBBA3(1,model="linear", mean.linear = -10)+
            range_effect_OBBA2(1,model="linear", mean.linear = -0.046, prec.linear = 10000)+
            HSS(main = Hours_Since_Sunrise, model = HSS_spde)+
            kappa(square_atlas, model = "iid", constr = TRUE, hyper = list(prec = pc_prec)) +
            spde_OBBA2(main = geometry, model = matern_abund) +
            spde_change(main = geometry, model = matern_change) +',
                                       paste0("Beta1_",covariates_to_include,'(1,model="linear", mean.linear = 0, prec.linear = ', prec_linear[1],')', collapse = " + "),
                                       '+',
                                       paste0("Beta2_",covariates_to_include,'(1,model="linear", mean.linear = 0, prec.linear = ', prec_linear[2],')', collapse = " + ")))
  
  
  model_formula_OBBA2 = as.formula(paste0('count ~
                  Intercept_OBBA2 +
                  log_offset + 
                  kappa + 
                  range_effect_OBBA2 * distance_from_range +
                  HSS +
                  spde_OBBA2 +',
                                          paste0("Beta1_",covariates_to_include,'*',covariates_to_include, collapse = " + "),
                                          '+',
                                          paste0("Beta2_",covariates_to_include,'*',covariates_to_include,"^2", collapse = " + ")))
  
  model_formula_OBBA3 = as.formula(paste0('count ~
                  Intercept_OBBA3 +
                  log_offset + 
                  kappa + 
                  range_effect_OBBA2 * distance_from_range +
                  HSS +
                  spde_OBBA2 + spde_change + ',
                                          paste0("Beta1_",covariates_to_include,'*',covariates_to_include, collapse = " + "),
                                          '+',
                                          paste0("Beta2_",covariates_to_include,'*',covariates_to_include,"^2", collapse = " + ")))
  
  
  fit_model <- function(error_type = NULL,
                        bru_init = list(Intercept_OBBA2 = -10,
                                        Intercept_OBBA3 = -10)){
    tryCatch(expr = {bru(components = model_components,
                         
                         like(family = error_type,
                              formula = model_formula_OBBA2,
                              data = subset(species_square, year(Date_Time)%in% c(2001:2005))),
                         
                         like(family = error_type,
                              formula = model_formula_OBBA3,
                              data = subset(species_square, year(Date_Time)%in% c(2021:2025))),
                         
                         options = list(inla.mode = "experimental",
                                        bru_initial = bru_init,
                                        control.compute = list(waic = FALSE, cpo = FALSE),
                                        bru_verbose = 4))},
             error = function(e){NULL})
  }
  
  start <- Sys.time()
  fit_INLA <- NULL
  for (tries in 2){
    
    if (tries == 1)  etype <- "nbinomial"
    if (tries == 2)  etype <- "poisson"
    
    fit_INLA <- fit_model(error_type = etype)
    
    if (!is.null(fit_INLA) & ("try-error" %!in% class(fit_INLA))) break
    
  } # tries
  
  if ("try-error" %in% class(fit_INLA) | is.null(fit_INLA)) next
  
  end <- Sys.time()
  runtime_INLA <- difftime( end,start, units="mins") %>% round(2)
  print(paste0(species_name," - ",runtime_INLA," min to fit model"))  # ~8 min
  print(summary(fit_INLA))
  
  #save(fit_INLA, file = paste0("../Output/Fitted_Models/",species_name,".RData"))
  

  # ****************************************************************************
  # ****************************************************************************
  # GENERATE AND STORE PREDICTIONS
  # ****************************************************************************
  # ****************************************************************************
  
  # ----------------------------------------------------------------------------
  # Predictions of covariate effects (holding all other covariates at zero)
  # ----------------------------------------------------------------------------
  
  covar_values <- as.data.frame(ONGrid_species)[,covariates_to_include]
  cov_df <- data.frame()
  for (covar in covariates_to_include){
    tmp <- data.frame(i = 1:1000)
    tmp[,covar] = seq(min(covar_values[,covar],na.rm = TRUE),max(covar_values[,covar],na.rm = TRUE),length.out = 1000)
    cov_df <- bind_rows(cov_df,tmp)
  }
  
  # add hours since sunrise covariate
  tmp <- data.frame(i = 1:1000)
  tmp[,"Hours_Since_Sunrise"] = seq(min(species_square$Hours_Since_Sunrise,na.rm = TRUE),max(species_square$Hours_Since_Sunrise,na.rm = TRUE),length.out = 1000)
  cov_df <- bind_rows(cov_df,tmp)
  
  cov_df <- cov_df[,-1]
  cov_df[is.na(cov_df)] <- 0
  
  cov_pred <- generate(fit_INLA,
                   cov_df,
                   formula =  
                     as.formula(paste0("~ Intercept_OBBA3 +log_offset_5min + HSS +" ,
                                       paste0("Beta1_",covariates_to_include,'*',covariates_to_include, collapse = " + "),
                                       '+',
                                       paste0("Beta2_",covariates_to_include,'*',covariates_to_include,"^2", collapse = " + "))) ,
                   n.samples = 500,
                   seed = 123) %>% exp()
  
  cov_df$pred_q50 <- apply(cov_pred,1, function(x)quantile(x,0.5))
  cov_df$pred_q05 <- apply(cov_pred,1, function(x)quantile(x,0.05))
  cov_df$pred_q95 <- apply(cov_pred,1, function(x)quantile(x,0.95))
  
  if (file.exists(paste0("../Output/Tables_Summaries/covariate_effect_table.csv"))) covariate_effect_table <- read.csv(paste0("../Output/Tables_Summaries/covariate_effect_table.csv"), header = TRUE)
  
  cov_df <- cov_df %>% 
    mutate(species_name = species_name, sp_code = sp_code, error_type = etype) %>%
    dplyr::relocate(species_name,sp_code,error_type)
  
  covariate_effect_table <- rbind(covariate_effect_table,cov_df)
  write.csv(covariate_effect_table,paste0("../Output/Tables_Summaries/covariate_effect_table.csv"), row.names = FALSE)
  
  ggplot(subset(cov_df, Hours_Since_Sunrise != 0), aes(x = Hours_Since_Sunrise, y = pred_q50, ymin = pred_q05, ymax = pred_q95))+
    geom_ribbon(alpha = 0.1)+
    geom_line()+
    theme_bw()
  
  # ----------------------------------------------------------------------------
  # Spatial predictions
  # ----------------------------------------------------------------------------
  
  # Colour scale for relative abundance maps
  colscale_relabund <-c("#FEFEFE", "#FBF7E2", "#FCF8D0", "#EEF7C2", "#CEF2B0", "#94E5A0", "#51C987", "#18A065", "#008C59", "#007F53", "#006344")
  colpal_relabund <- colorRampPalette(colscale_relabund)
  
  # For every pixel on landscape, extract distance (in km) from eBird range limit
  ONGrid_species <- ONGrid %>%
    st_transform(st_crs(species_square)) %>%
    mutate(distance_from_range = (st_distance( . , range) %>% as.numeric()),
           Hours_Since_Sunrise = 0)
  
  # Convert to CRS of target raster
  #ONGrid_species <- ONGrid_species %>% st_transform(st_crs(target_raster)) %>% st_centroid(.)
  
  start2 <- Sys.time()
  
  # Generate predictions for every pixel
  pred <- generate(fit_INLA,
                   ONGrid_species,
                   formula =  
                     as.formula(paste0("~ data.frame(OBBA2 = Intercept_OBBA2 +log_offset_5min + range_effect_OBBA2 * distance_from_range + spde_OBBA2 +" ,
                                       paste0("Beta1_",covariates_to_include,'*',covariates_to_include, collapse = " + "),
                                       '+',
                                       paste0("Beta2_",covariates_to_include,'*',covariates_to_include,"^2", collapse = " + "),
                                       ",OBBA3 = Intercept_OBBA3 + log_offset_5min + range_effect_OBBA2 * distance_from_range + spde_OBBA2 + spde_change +" ,
                                       paste0("Beta1_",covariates_to_include,'*',covariates_to_include, collapse = " + "),
                                       '+',
                                       paste0("Beta2_",covariates_to_include,'*',covariates_to_include,"^2", collapse = " + "),
                                       ")" )) ,
                   n.samples = 500,
                   seed = 123)
  
  # Separate into predictions for OBBA2 and OBBA3
  pred_OBBA2 = sapply(pred, function(x) exp(x$OBBA2))
  pred_OBBA3 = sapply(pred, function(x) exp(x$OBBA3))
  
  # Calculate absolute change between atlases (per pixel)
  pred_change = pred_OBBA3-pred_OBBA2
  
  # Estimate of overall change across Ontario
  percent_change_overall = 100*(apply(pred_OBBA3,2,sum,na.rm = TRUE)/apply(pred_OBBA2,2,sum,na.rm = TRUE)-1)
  
  end2 <- Sys.time() 
  runtime_pred <- difftime( end2,start2, units="mins") %>% round(2)
  print(paste0(species_name," - ",runtime_pred," min to generate predictions")) # 7 min
  
  # Median of posterior for each pixel
  ONGrid_species$OBBA2_pred_q50 <- apply(pred_OBBA2,1,median)
  ONGrid_species$OBBA3_pred_q50 <- apply(pred_OBBA3,1,median)
  
  # ------------------------------------------------
  # clipped ebird range limit
  # ------------------------------------------------
  
  if (!is.na(range)) range_clipped <- range  %>%
    st_transform(st_crs(ONBoundary)) %>%
    st_intersection(ONBoundary)
  
  # ------------------------------------------------
  # PLOT RELATIVE ABUNDANCE IN OBBA2
  # ------------------------------------------------
  
  # Bounds for plotting
  upper_bound <- quantile(ONGrid_species$OBBA3_pred_q50,0.99,na.rm = TRUE) %>% signif(2)
  lower_bound <- upper_bound/100
  
  ONGrid_species$OBBA2_pred_q50[ONGrid_species$OBBA2_pred_q50 > upper_bound] <-  upper_bound
  ONGrid_species$OBBA2_pred_q50[ONGrid_species$OBBA2_pred_q50 < lower_bound] <-  0
  
  OBBA2_plot_q50 <- ggplot() +
    
    geom_sf(data = ONGrid_species, aes(col = OBBA2_pred_q50), size = 0.1) +
    scale_color_gradientn(name = "Per point count",
                          colors = colpal_relabund(10),
                          trans = "log10",
                          na.value = "white")+
    
    geom_sf(data = ONBoundary,colour="black",fill=NA,lwd=0.3,show.legend = F) +
    
    geom_sf(data = range_clipped, linetype = 2, col = "gray50", fill = "transparent", size = 2)+
    
    coord_sf(clip = "off",xlim = range(as.data.frame(st_coordinates(ONBoundary))$X)) +
    theme(panel.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    
    annotate(geom="text",x=400,y=1650, label="Predicted OBBA2",lineheight = .85,hjust = 0,size=6,fontface =2) +
    annotate(geom="text",x=400,y=1480, label= paste0("Prepared on ",Sys.Date()),size=2,lineheight = .75,hjust = 0,color="gray75")
  
  # ------------------------------------------------
  # PLOT RELATIVE ABUNDANCE IN OBBA3
  # ------------------------------------------------
  
  # Median of posterior for each pixel
  ONGrid_species$OBBA3_pred_q50[ONGrid_species$OBBA3_pred_q50 > upper_bound] <-  upper_bound
  ONGrid_species$OBBA3_pred_q50[ONGrid_species$OBBA3_pred_q50 < lower_bound] <-  0
  
  OBBA3_plot_q50 <- ggplot() +
    
    geom_sf(data = ONGrid_species, aes(col = OBBA3_pred_q50), size = 0.1) +
    scale_color_gradientn(name = "Per point count",
                          colors = colpal_relabund(10),
                          trans = "log10",
                          na.value = "white")+
    
    geom_sf(data = ONBoundary,colour="black",fill=NA,lwd=0.3,show.legend = F) +
    
    geom_sf(data = range_clipped, linetype = 2, col = "gray50", fill = "transparent", size = 2)+
    
    coord_sf(clip = "off",xlim = range(as.data.frame(st_coordinates(ONBoundary))$X)) +
    theme(panel.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    
    annotate(geom="text",x=400,y=1650, label="Predicted OBBA3",lineheight = .85,hjust = 0,size=6,fontface =2) +
    annotate(geom="text",x=400,y=1480, label= paste0("Prepared on ",Sys.Date()),size=2,lineheight = .75,hjust = 0,color="gray75")
  
  # ------------------------------------------------
  # PLOT CHANGE BETWEEN ATLASES
  # ------------------------------------------------
  
  colscale_change <-c("darkred","orangered","white","dodgerblue","darkblue")
  colpal_change <- colorRampPalette(colscale_change)
  
  ONGrid_species$change_pred_q50 <- apply(pred_change,1,median)
  ONGrid_species$change_pred_q50[ONGrid_species$change_pred_q50 > upper_bound] <-  upper_bound
  ONGrid_species$change_pred_q50[ONGrid_species$change_pred_q50 < -upper_bound] <-  -upper_bound
  
  # Pixels where trend is "significantly different from zero"
  #ONGrid_species$trend_signif <- ONGrid_species$change_pred_q95 < 0 | ONGrid_species$change_pred2_q05 > 0
  
  # Label for total overall change in Ontario
  convert_numeric_to_char <- function(numeric_vector) {
    char_vector <- ifelse(numeric_vector > 0, paste0("+", numeric_vector), as.character(numeric_vector))
    return(char_vector)
  }
  percent_change_quantiles <- quantile(percent_change_overall,c(0.05,0.5,0.95)) %>% round() %>% convert_numeric_to_char()
  change_label <- paste0("Overall change = ",percent_change_quantiles[2],"% (",percent_change_quantiles[1]," to ",percent_change_quantiles[3],")")
  
  
  #lim <- max(abs(c(ONGrid_species$change_pred_q50,upper_bound)),na.rm = TRUE)
  change_plot_q50 <- ggplot() +
    
    geom_sf(data = ONGrid_species, aes(col = change_pred_q50), size = 0.1) +
    scale_color_gradientn(name = "Absolute\nchange\n",
                          colors = colpal_change(5),
                          na.value = "black",
                          limits = c(-upper_bound,upper_bound))+
    
    
    geom_sf(data = ONBoundary,colour="black",fill=NA,lwd=0.3,show.legend = F) +

    coord_sf(clip = "off",xlim = range(as.data.frame(st_coordinates(ONBoundary))$X)) +
    theme(panel.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    annotate(geom="text",x=400,y=1650, label= "Change between atlases",lineheight = .85,hjust = 0,size=6,fontface =2) +
    annotate(geom="text",x=400,y=1550, label= change_label,lineheight = .85,hjust = 0,size=4,fontface =2) +
    annotate(geom="text",x=400,y=1480, label= paste0("Prepared on ",Sys.Date()),size=2,lineheight = .75,hjust = 0,color="gray75")+
    annotation_scale()
  print(change_plot_q50)
  
  # ------------------------------------------------
  # Combine all figures into a single one
  # ------------------------------------------------
  
  png(file = paste0("../Output/Prediction_Maps/Relative_Abundance/",species_name,".png"),
      units = "in", 
      width = 30, 
      height = 14, res = 600)
  
  grid.arrange(species_data_OBBA2,species_data_OBBA3, OBBA2_plot_q50,OBBA3_plot_q50,change_plot_q50,
               ncol = 3,nrow=4, 
               layout_matrix = rbind(c(1,3,NA),
                                     c(1,3,5),
                                     c(2,4,5),
                                     c(2,4,NA)),
               top = species_name)
  
  dev.off()
  
  print(summary(fit_INLA))
  
  # ------------------------------------------------
  # Store estimates of change for comparisons to BBS and eBird
  # ------------------------------------------------
  
  if (file.exists(paste0("../Output/Tables_Summaries/OBBA3_change_estimate_table.csv"))) change_estimate_table <- read.csv(paste0("../Output/Tables_Summaries/OBBA3_change_estimate_table.csv"), header = TRUE)
  
  # Overall change across Ontario
  species_change <- data.frame(species_name = species_name,
                               sp_code = sp_code,
                               error_type = etype,
                               region = "ON",
                               mean_abund = (sum(ONGrid_species$OBBA2_pred_q50,na.rm=TRUE) + sum(ONGrid_species$OBBA3_pred_q50,na.rm=TRUE))/2,
                               change_q50 = quantile(percent_change_overall,0.5),
                               change_q05 = quantile(percent_change_overall,0.05),
                               change_q95 = quantile(percent_change_overall,0.95),
                               prob_decrease_0_percent = mean(percent_change_overall <= 0),
                               prob_decrease_30_percent = mean(percent_change_overall <= -30),
                               prob_decrease_50_percent = mean(percent_change_overall <= -50))
  
  # Change within each BCR
  for (bcr in unique(ON_BCR$BCR)){
    
    percent_change_bcr <- 100*(apply(pred_OBBA3[which(ONGrid$BCR == bcr),],2,sum,na.rm = TRUE) - apply(pred_OBBA2[which(ONGrid$BCR == bcr),],2,sum,na.rm = TRUE))/apply(pred_OBBA2[which(ONGrid$BCR == bcr),],2,sum,na.rm = TRUE)
    
    species_change <- rbind(species_change,data.frame(species_name = species_name,
                                                      sp_code = sp_code,
                                                      error_type = etype,
                                                      region = paste0("CA-ON-",bcr),
                                                      mean_abund = (sum(ONGrid_species$OBBA2_pred_q50[which(ONGrid$BCR == bcr)],na.rm=TRUE) + sum(ONGrid_species$OBBA3_pred_q50[which(ONGrid$BCR == bcr)],na.rm=TRUE))/2,
                                                      change_q50 = quantile(percent_change_bcr,0.5),
                                                      change_q05 = quantile(percent_change_bcr,0.05),
                                                      change_q95 = quantile(percent_change_bcr,0.95),
                                                      prob_decrease_0_percent = mean(percent_change_bcr <= 0),
                                                      prob_decrease_30_percent = mean(percent_change_bcr <= -30),
                                                      prob_decrease_50_percent = mean(percent_change_bcr <= -50)))
  }
  
  change_estimate_table <- rbind(change_estimate_table,species_change)
  write.csv(change_estimate_table,paste0("../Output/Tables_Summaries/OBBA3_change_estimate_table.csv"), row.names = FALSE)
  
} # close species loop
