# ************************************************
# BAYESIAN ANALYSIS / SPECIES DISTRIBUTION MODELS FOR ONTARIO BREEDING BIRD ATLAS
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
require(scales)

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
  
  tgt <- st_as_stars(target_raster)
  tmp = stars::st_rasterize(df %>% dplyr::select(levs, geometry),
                            nx = dim(tgt)[1],ny = dim(tgt)[2])
  
  return(list(raster = tmp,cut_levs = cut_levs))
}

# ------------------------------------------------
# Function to fit boosted regression tree models
# ------------------------------------------------

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

analysis_data <- readRDS("../Data_Cleaned/analysis_data_package.rds")
attach(analysis_data)

AEA_proj <- "+proj=aea +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-87 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs "

ONBoundary <- st_read("../../../Data/Spatial/National/BCR/BCR_Terrestrial_master.shp")  %>%
  subset(PROVINCE_S == "ONTARIO") %>%
  st_make_valid() %>%
  st_union() %>%
  st_transform(st_crs(AEA_proj))

ONBoundary_buffer <- ONBoundary %>% st_buffer(20000)

# Raster with target properties
target_raster <- rast("../../../Data/Spatial/National/AnnualMeanTemperature/wc2.1_30s_bio_1.tif") %>% 
  crop(st_transform(ONBoundary_buffer,crs(.))) %>% 
  project(AEA_proj, res = 500) %>%
  mask(vect(ONBoundary))

# ONWater <- read_sf("../../../Data/Spatial/ONatchewan/ONWater/ONWaterClip.shp") %>% 
#   st_transform(st_crs(target_raster))
# ONWater$area <- st_area(ONWater)
# ONWater$area <- as.numeric(ONWater$area)
# ONWater <- ONWater[ONWater$area>2.580e+06 ,]
#

# Atlas Squares # 10 km x 10 km grid
ONSquares <- st_make_grid(
  ONBoundary,
  cellsize = units::set_units(10*10,km^2),
  what = "polygons",
  square = TRUE,
  flat_topped = FALSE)%>%
  st_as_sf() %>%
  st_intersection(ONBoundary) %>%
  na.omit() %>%
  mutate(sq_id = 1:nrow(.)) %>%
  rename(geometry = x)

BCR_PROV <- st_read("../../../Data/Spatial/National/BCR/BCR_Terrestrial_master.shp")  %>%
  subset(PROVINCE_S == "ONTARIO") %>%
  st_make_valid() %>%
  group_by(BCR,PROVINCE_S)%>%
  summarize(geometry = st_union(geometry)) %>%
  st_transform(st_crs(AEA_proj))

ONGrid <- ONGrid %>%
  st_intersection(BCR_PROV) %>%
  dplyr::rename(geometry = x)

ONGrid_centroids <- st_centroid(ONGrid)

# ------------------------------------------------
# Simulate a spatial change process
# ------------------------------------------------

# Grid on which to simulate change (assume change will be a coarse process)
change_grid <- st_make_grid(
  ONBoundary,
  cellsize = units::set_units(15*15,km^2),
  what = "polygons",
  square = TRUE,
  flat_topped = FALSE)%>%
  st_as_sf() %>%
  rename(geometry = x)


coords <- change_grid %>% st_centroid() %>% st_coordinates()
change_field_sim <- geoR::grf(grid = coords,
                              cov.pars = c(0.5,1000000))

mean_change <- -0.5
sd_change <- 1
change_sim <- scale(change_field_sim$data)*sd_change + mean_change

change_grid$change <- change_sim

ONGrid <- st_intersection(ONGrid %>% st_centroid(),change_grid) %>%
  as.data.frame() %>%
  group_by(point_id)%>%
  summarize(change = mean(change)) %>%
  dplyr::select(point_id,change)%>%
  left_join(ONGrid,.)

change_rast <- stars::st_rasterize(ONGrid %>% dplyr::select(change, geometry),
                                   nx = dim(target_raster)[1],ny = dim(target_raster)[2])

ggplot()+
  geom_stars(data = change_rast)+
  scale_fill_gradient2(low=muted("red"),mid="white",high=muted("blue"),midpoint=0, guide = "none")+
  theme_bw()

# ------------------------------------------------
# Choose a focal species, fit a species distribution model to generate a plausible
# model of 'current' distribution, and then simulate data to re-fit models and evaluate
# if the change surface can be recovered
# ------------------------------------------------

species_to_model <- species_to_model %>%
  arrange(desc(n_squares),desc(n_detections)) %>%
  subset(n_detections <= 1000)

sp_code = 'CAWA'

print(sp_code)

# **************************************************************
# **************************************************************
# PART 1: EMPIRICAL ANALYSIS 
# **************************************************************
# **************************************************************

# Prepare data for this species
sp_dat <- all_surveys %>% 
  mutate(count = full_count_matrix[,sp_code]) %>%
  
  # Only select point counts
  subset(Survey_Type %in% c("Point_Count","ARU_SPT","ARU_SPM")) %>%
  st_transform(AEA_proj)

# n_det_sq <- subset(sp_dat,count>0) %>%
#   as.data.frame() %>%
# 
#   summarize(n_sq = length(unique(sq_id)),
#             n_det = sum(count>0))
# 
# if (n_det_sq$n_sq < 30 | n_det_sq$n_det < 60) next

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
# Generate QPAD offsets for each survey (assumes unlimited distance point counts)
# --------------------------------

species_offsets <- subset(species_to_model, Species_Code_BSC == sp_code)

if (species_offsets$offset_exists == FALSE) sp_dat$log_QPAD_offset <- 0

if (species_offsets$offset_exists == TRUE){
  A_metres <- pi*species_offsets$EDR^2
  p <- 1-exp(-sp_dat$Survey_Duration_Minutes*species_offsets$cue_rate)
  sp_dat$log_QPAD_offset <- log(A_metres * p)
}


covariates_to_include <- paste0("PC",1:8)

# ------------------------------------------------
# Create a spatial mesh, which is used to fit the residual spatial field
# ------------------------------------------------

# make a two extension hulls and mesh for spatial model
hull <- fm_extensions(
  ONBoundary,
  convex = c(50000, 200000),
  concave = c(350000, 500000)
)
mesh_spatial <- fm_mesh_2d_inla(
  boundary = hull, 
  max.edge = c(70000, 100000), # km inside and outside
  cutoff = 30000, 
  crs = fm_crs(sp_dat)
) # cutoff is min edge
mesh_locs <- mesh_spatial$loc[,c(1,2)] %>% as.data.frame()
dim(mesh_locs)
#plot(mesh_spatial)
#lines(BCR_PROV)


# # Note: mesh developed using tutorial at: https://rpubs.com/jafet089/886687
# max.edge = 50000 #diff(range(st_coordinates(sp_dat)[,1]))/30
# bound.outer = diff(range(st_coordinates(sp_dat)[,1]))/3
# cutoff = max.edge/5
# 
# mesh_spatial <- inla.mesh.2d(loc = st_coordinates(sp_dat),
#                              cutoff = max.edge/2,
#                              max.edge = c(1,2)*max.edge,
#                              offset=c(max.edge, bound.outer))
# mesh_locs <- mesh_spatial$loc[,c(1,2)] %>% as.data.frame()
# dim(mesh_locs)
# plot(mesh_spatial)
# lines(BCR_PROV)

prior_range <- c(300000,0.1) # 10% chance range is smaller than 300000
prior_sigma <- c(0.5,0.1) # 10% chance sd is larger than 0.5
matern_coarse <- inla.spde2.pcmatern(mesh_spatial,
                                     prior.range = prior_range, 
                                     prior.sigma = prior_sigma 
)

# ------------------------------------------------
# Create mesh to model effect of time since sunrise (TSS)
# ------------------------------------------------
sp_dat$Hours_Since_Sunrise <- as.numeric(sp_dat$Hours_Since_Sunrise)
TSS_range <- range(sp_dat$Hours_Since_Sunrise)
TSS_meshpoints <- seq(TSS_range[1]-0.1,TSS_range[2]+0.1,length.out = 11)
TSS_mesh1D = inla.mesh.1d(TSS_meshpoints,boundary="free")
TSS_spde = inla.spde2.pcmatern(TSS_mesh1D,
                               prior.range = c(6,0.1),
                               prior.sigma = c(1,0.1)) # 10% chance sd is larger than 1

# ------------------------------------------------
# Model formulas
# ------------------------------------------------

sd_linear <- 1
prec_linear <-  c(1/sd_linear^2,1/(sd_linear/2)^2)
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


# ------------------------------------------------
# Fit with nbinomial error
# ------------------------------------------------

PC_sp <- sp_dat %>% as('Spatial')
start <- Sys.time()
fit_INLA <- NULL
while(is.null(fit_INLA)){
  fit_INLA <- bru(components = model_components,
                  
                  like(family = "nbinomial",
                       formula = model_formula_PC,
                       data = PC_sp),
                  
                  options = list(
                    control.compute = list(waic = FALSE, cpo = FALSE),
                    bru_verbose = 4))
  if ("try-error" %in% class(fit_INLA)) fit_INLA <- NULL
}

end <- Sys.time()
runtime_INLA <- difftime( end,start, units="mins") %>% round(2)
print(paste0(sp_code," - ",runtime_INLA," min to fit model"))

# For every pixel on landscape, extract distance from eBird range limit
distance_from_range = (st_centroid(ONGrid) %>% 
                         st_distance( . , range) %>% 
                         as.numeric())/1000

ONGrid_species <- ONGrid %>% mutate(distance_from_range = distance_from_range)

# --------------------------------
# QPAD offsets associated with a 5-minute unlimited distance survey
# --------------------------------

species_offsets <- subset(species_to_model, Species_Code_BSC == sp_code)

log_offset_5min <- 0
if (species_offsets$offset_exists == TRUE) log_offset_5min <- species_offsets$log_offset_5min

# ------------------------------------------------
# Generate predictions on ONGrid_species raster (1 km x 1 km pixels)
# ------------------------------------------------

pred_formula_PC = as.formula(paste0(' ~
                  Intercept_PC +
                  log_offset_5min +
                  range_effect * distance_from_range +
                  spde_coarse +',
                                    paste0("Beta1_",covariates_to_include,'*',covariates_to_include, collapse = " + "),
                                    '+',
                                    paste0("Beta2_",covariates_to_include,'*',covariates_to_include, collapse = " + ")))

# Predictions are on log scale, and do not include variance components
pred <- NULL
pred <- generate(fit_INLA,
                 as(ONGrid_species,'Spatial'),
                 formula =  pred_formula_PC,
                 n.samples = 1)

pred <- exp(pred)

ONGrid_species$mean_cycle2 <- pred
ONGrid_species$mean_cycle1 <- exp(log(ONGrid_species$mean_cycle2) - ONGrid_species$change)

# ****************************************************************************
# ****************************************************************************
# Simulate data during each of two periods, at each existing survey location
# ****************************************************************************
# ****************************************************************************

colscale_relabund <-c("#FEFEFE", "#FBF7E2", "#FCF8D0", "#EEF7C2", "#CEF2B0", "#94E5A0", "#51C987", "#18A065", "#008C59", "#007F53", "#006344")
colpal_relabund <- colorRampPalette(colscale_relabund)

lower_bound <- 0.01
upper_bound <- quantile(c(ONGrid_species$mean_cycle2,ONGrid_species$mean_cycle1),0.99,na.rm = TRUE) %>% signif(2)

if (lower_bound >= (upper_bound/5)) lower_bound <- (upper_bound/5) %>% signif(2)

# -----------------------
# Simulated abundance in cycle 1
# -----------------------

cut1 <- cut.fn(df = ONGrid_species,
               target_raster = target_raster,
               column_name = "mean_cycle1",
               lower_bound = lower_bound,
               upper_bound = upper_bound)

raster_cycle1 <- cut1$raster

# Median of posterior
plot_cycle1 <- ggplot() +
  
  geom_stars(data = raster_cycle1) +
  scale_fill_manual(name = "<span style='font-size:13pt'>Relative abundance</span><br><span style='font-size:7pt'>Per 5 minute point count</span><br><span style='font-size:7pt'>Simulated</span>",
                    values = colpal_relabund(length(levels(raster_cycle1$levs))), drop=FALSE,na.translate=FALSE)+
  
  geom_sf(data = ONBoundary,colour="black",fill=NA,lwd=0.3,show.legend = F) +
  coord_sf(clip = "off",xlim = range(as.data.frame(st_coordinates(ONBoundary))$X))+
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
  theme(legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(5,10,5,-20),
        legend.title.align=0.5,
        legend.title = element_markdown(lineheight=.9,hjust = "left"))+
  theme(legend.key = element_rect(fill = "transparent", colour = "transparent"))

# -----------------------
# Simulated abundance in cycle 2
# -----------------------

cut2 <- cut.fn(df = ONGrid_species,
               target_raster = target_raster,
               column_name = "mean_cycle2",
               lower_bound = lower_bound,
               upper_bound = upper_bound)

raster_cycle2 <- cut2$raster

# Median of posterior
plot_cycle2 <- ggplot() +
  
  geom_stars(data = raster_cycle2) +
  scale_fill_manual(name = "<span style='font-size:13pt'>Relative abundance</span><br><span style='font-size:7pt'>Per 5 minute point count</span><br><span style='font-size:7pt'>Simulated</span>",
                    values = colpal_relabund(length(levels(raster_cycle2$levs))), drop=FALSE,na.translate=FALSE)+
  
  geom_sf(data = ONBoundary,colour="black",fill=NA,lwd=0.3,show.legend = F) +
  coord_sf(clip = "off",xlim = range(as.data.frame(st_coordinates(ONBoundary))$X))+
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
  theme(legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(5,10,5,-20),
        legend.title.align=0.5,
        legend.title = element_markdown(lineheight=.9,hjust = "left"))+
  theme(legend.key = element_rect(fill = "transparent", colour = "transparent"))

