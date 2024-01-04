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
BBMP_boundary <- st_transform(BBMP_boundary, proj)

max_lat <- st_coordinates(BBMP_boundary)[,2] %>% max()
min_lat <- st_coordinates(BBMP_boundary)[,2] %>% min()
mid_lat <- mean(c(min_lat,max_lat))

max_lon <- st_coordinates(BBMP_boundary)[,1] %>% max()
min_lon <- st_coordinates(BBMP_boundary)[,1] %>% min()
mid_lon <- mean(c(min_lon,max_lon))

proj <- paste(
  "+proj=aea +lat_0=",mid_lat,"+lon_0=",mid_lon,"+lat_1=",min_lat,"+lat_2=",max_lat,"+x_0=0 +y_0=0 +datum=NAD83",
  "+units=km +no_defs"
)

BBMP_boundary <- st_transform(BBMP_boundary,proj)
alldat <- alldat %>% st_transform(proj)

# pt.bond <- inla.nonconvex.hull(st_coordinates(BBMP_boundary)[,c(1,2)],crs=proj, convex = -0.07)
# mesh <- inla.mesh.2d(boundary=pt.bond, 
#                      max.edge=c(70,200),
#                      offset = c(10,100),
#                      cut=50,
#                      crs=proj)

hull <- fm_extensions(BBMP_boundary, crs = proj,
                      convex = -0.02,
                      concave = -20)

# Setting max.edge = c(30000,100000) requires longer to fit (7 min), but may increase model accuracy
# Setting max.edge = c(50000,100000) allows the model to fit much more quickly (3.5 min), though possibly with reduced accuracy
mesh <- fm_mesh_2d_inla(
  boundary = hull, 
  max.edge = c(50, 200),
  cutoff = 50, 
  crs = fm_crs(alldat)
)

plot(mesh)
lines(as(BBMP_boundary,'Spatial'))
mesh_locs <- mesh$loc[,c(1,2)] %>% as.data.frame()
dim(mesh_locs)


# ---------------------------------------------------
# Grid on which predictions will be made
# ---------------------------------------------------

# Path to spatial covariates
covar_folder <- "../../../Data/Spatial/"

# Mean Annual Temperature
MAT <- rast(paste0(covar_folder,"National/Bioclimate/Normal_1991_2020_MAT.tif"))
# Mean Annual Precipitation
MAP <- rast(paste0(covar_folder,"National/Bioclimate/Normal_1991_2020_MAP.tif"))

bnd <- st_transform(BBMP_boundary,crs(MAT)) %>% as('Spatial')
  
MAT <- MAT %>%
  crop(bnd)%>%
  mask(vect(bnd)) %>% 
  project(proj,method = "bilinear", res = c(5,5))

MAP <- MAP %>%
  crop(bnd)%>%
  mask(vect(bnd)) %>% 
  project(proj,method = "bilinear", res = c(5,5))

pred_grid <- rast(list(MAT,MAP)) %>% terra::extract(.,y = 1:((dim(MAT)[1]*dim(MAT)[2])),xy=TRUE)
colnames(pred_grid) <- c("x","y","MAT_1km","MAP_1km")
  
pred_grid <- pred_grid %>% st_as_sf(coords = c("x", "y"),crs=proj, remove = FALSE)

# ---------------------------------------------------
# Standardize covariates
# ---------------------------------------------------
covariates_to_include <- c("MAT_1km","MAP_1km")
alldat <- bind_cols(alldat,st_coordinates(alldat))  %>% as.data.frame(alldat)
pred_grid <- bind_cols(pred_grid,st_coordinates(pred_grid)) %>%as.data.frame()


for (cov in covariates_to_include){
  cov_mu <- mean(alldat[,cov])
  cov_sd <- sd(alldat[,cov])
  
  alldat[,cov] <- (alldat[,cov] - cov_mu)/cov_sd
  pred_grid[,cov] <- (pred_grid[,cov] - cov_mu)/cov_sd
  
}
pred_grid <- pred_grid %>% st_as_sf(coords = c("X", "Y"),crs=proj, remove = FALSE)
alldat <- alldat %>% st_as_sf(coords = c("X", "Y"),crs=proj, remove = FALSE)

# ---------------------------------------------------
# Fit models and generate predictions
# ---------------------------------------------------

# Timeout INLA after 10 minutes (if it has not fit by then, it has likely stalled)
inla.setOption(inla.timeout = 60*20)

# Relative abundances of each species
species_relabund <- allcounts 
species_relabund[species_relabund>0] <- 1
species_relabund <- apply(species_relabund,2,sum) %>% sort(decreasing = TRUE)

sp_code = "ATSP"

# Prepare data for this species
sp_dat <- alldat %>%
  mutate(count = allcounts[,sp_code]) %>%
  mutate(presence = as.numeric(count>0)) %>%
  subset(duration >=1 & duration <= 10)

# # ----------------------------------------------------
# # Generate QPAD offsets for each survey (assumes unlimited distance point counts)
# # ----------------------------------------------------
# sp_napops <- subset(napops_species,Species_Code_NAPOPS == sp_code)
# sp_dat$log_offset <- 0
# if (nrow(sp_napops)>0){
#   if (sp_napops$Removal == 1 & sp_napops$Distance == 1){
#     
#     cue_rate <- cue_rate(species = sp_code,od = 153, tssr = 0, model = 1)[3] %>% as.numeric()
#     EDR <- edr(species = sp_code,road = FALSE, forest = 0.5,model = 1)[3] %>% as.numeric()
#     
#     # Calculate A and p, which jointly determine offset
#     A_metres <- c(pi*EDR^2)
#     p <- 1-exp(-sp_dat$duration*cue_rate)
#     
#     sp_dat$log_offset <- log(A_metres * p)
#     
#   }
# }

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

training_data_bam <-  subset(sp_dat, project != "BBMP" & organization %!in% c("CWS-NOR","CWS-ONT","CWS-PRA","ECCC"))
training_data_bam$data_type <- "bam"

training_data_bbmp <-  subset(sp_dat, project == "BBMP" | organization %in% c("CWS-NOR","CWS-ONT","CWS-PRA","ECCC"))
training_data_bbmp$data_type <- "bbmp"

training_data_bam <- training_data_bam %>% as('Spatial')
training_data_bbmp <- training_data_bbmp %>% as('Spatial')

# ----------------------------------------------------
# Prepare formulas and priors
# ----------------------------------------------------

# Controls the 'residual spatial field'.  This can be adjusted to create smoother surfaces.
prior_range <- c(1000,0.01)
prior_sigma <- c(1,0.01)

matern_coarse <- inla.spde2.pcmatern(mesh,
                                     prior.range = prior_range,
                                     prior.sigma = prior_sigma
)

# Mesh for survey duration effect
duration_meshpoints <- seq(min(sp_dat$duration)-0.1,max(sp_dat$duration)+0.1,length.out = 11)
duration_mesh1D = inla.mesh.1d(duration_meshpoints,boundary="free")
duration_spde = inla.spde2.pcmatern(duration_mesh1D,
                                    prior.range = c(5,0.1),
                                    prior.sigma = c(1,0.5))

# How much shrinkage should be applied to covariate effects?
sd_linear <- 1   # Change to smaller value (e.g., 0.1), if you want to heavily shrink covariate effects and potentially create smoother surfaces
prec_linear <-  c(1/sd_linear^2,1/(sd_linear/2)^2)

model_components = as.formula(paste0('~
            Intercept_BBMP(1)+
            Intercept_BAM(1)+
            spde_coarse(main = coordinates, model = matern_coarse) +', paste0("Beta1_",covariates_to_include,'(1,model="linear", mean.linear = 0, prec.linear = ', prec_linear[1],')', collapse = " + ")))

model_formula_BBMP = as.formula(paste0('count ~
                  Intercept_BBMP +
                  spde_coarse +',paste0("Beta1_",covariates_to_include,'*',covariates_to_include, collapse = " + ")))

model_formula_BAM = as.formula(paste0('count ~
                  Intercept_BAM +
                  spde_coarse +',paste0("Beta1_",covariates_to_include,'*',covariates_to_include, collapse = " + ")))

# ----------------------------------------------------------------------------
# Fit model to only non-bbmp data, assess cross-validation accuracy
# ----------------------------------------------------------------------------
start <- Sys.time()
fit_bam <- NULL
while(is.null(fit_bam)){
  
  fit_model <- function(){
    tryCatch(expr = {bru(components = model_components,
                         
                         like(family = "poisson",
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
                         
                         like(family = "poisson",
                              formula = model_formula_BBMP,
                              data = training_data_bbmp),
                         
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
# Prepare for predictions on grid
# ------------------------------------------------

pred_grid <- na.omit(pred_grid)
#pred_grid$log_offset_bam <- mean(training_data_bam$log_offset)
#pred_grid$log_offset_bbmp <- mean(training_data_bbmp$log_offset)
pred_grid <- pred_grid %>% as('Spatial')

# ------------------------------------------------
# Predictions for model fit to non-bbmp data
# ------------------------------------------------


pred_formula = as.formula(paste0(' ~Intercept_BAM +spde_coarse +',paste0("Beta1_",covariates_to_include,'*',covariates_to_include, collapse = " + ")))
size <- fit_bam$summary.hyperpar$'0.5quant'[1]
pred_bam <- NULL
pred_bam <- generate(fit_bam,
                     pred_grid,
                     formula =  pred_formula,
                     n.samples = 1000)

pred_bam <- exp(pred_bam)
pred_bam <- apply(pred_bam,1,median)
pred_grid$pred_bam <- pred_bam

# ------------------------------------------------
# Predictions for model fit to bbmp data
# ------------------------------------------------

pred_formula = as.formula(paste0(' ~Intercept_BBMP +spde_coarse +',paste0("Beta1_",covariates_to_include,'*',covariates_to_include, collapse = " + ")))
size <- fit_bbmp$summary.hyperpar$'0.5quant'[1]
pred_bbmp <- NULL
pred_bbmp <- generate(fit_bbmp,
                     pred_grid,
                     formula =  pred_formula,
                     n.samples = 1000)

pred_bbmp <- exp(pred_bbmp)
pred_bbmp <- apply(pred_bbmp,1,median)
pred_grid$pred_bbmp <- pred_bbmp

# ------------------------------------------------
# Maps
# ------------------------------------------------

mesh_locs <- mesh$loc %>% as.data.frame()

ggplot(as.data.frame(pred_grid))+
  geom_raster(aes(x = x, y = y, fill = pred_bam))+
  scale_fill_gradientn(colors = viridis(10), trans = "log10")

pred_grid$pred_bbmp[pred_grid$pred_bbmp>20] <- 20
ggplot(as.data.frame(pred_grid))+
  geom_raster(aes(x = x, y = y, fill = pred_bbmp))+
  #geom_point(data = mesh_locs, aes(x = V1, y = V2), size = 0.1,col = "white") +
  geom_point(data = as.data.frame(training_data_bbmp), aes(x = X, y = Y, size = count)) +
  geom_point(data = as.data.frame(training_data_bam), aes(x = X, y = Y, size = count)) +
  
  scale_fill_gradientn(colors = viridis(10))
