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

# Timeout INLA after 10 minutes (if it has not fit by then, it has likely stalled)
inla.setOption(inla.timeout = 60*10)

# ------------------------------------------------
# LOAD DATA FROM NATURECOUNTS / WILDTRAX
# ------------------------------------------------

analysis_data <- readRDS("../Data_Cleaned/analysis_data_package.rds")
attach(analysis_data)

AEA_proj <- "+proj=aea +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-87 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs "

ONBoundary <- st_read("../../../Data/Spatial/National/BCR/BCR_Terrestrial_master.shp")  %>%
  subset(PROVINCE_S == "ONTARIO") %>%
  st_make_valid() %>%
  st_union() %>%
  st_transform(st_crs(AEA_proj))

# Raster with target properties
target_raster <- rast("../../../Data/Spatial/National/AnnualMeanTemperature/wc2.1_30s_bio_1.tif") %>% 
  crop(st_transform(ONBoundary,crs(.))) %>% 
  project(AEA_proj, res = 500) %>%
  mask(vect(ONBoundary))

BCR_PROV <- st_read("../../../Data/Spatial/National/BCR/BCR_Terrestrial_master.shp")  %>%
  subset(PROVINCE_S == "ONTARIO") %>%
  st_make_valid() %>%
  group_by(BCR,PROVINCE_S)%>%
  summarize(geometry = st_union(geometry)) %>%
  st_transform(st_crs(AEA_proj))

# ----------------------------------------------------------------
# Create hexagon layer across province (for spatial random effects)
# ----------------------------------------------------------------

hex25 <- st_make_grid(
  ONBoundary,
  cellsize = units::set_units(25*25,km^2),
  what = "polygons",
  square = FALSE,
  flat_topped = TRUE) %>%
  st_as_sf()
hex25$hexID <- 1:nrow(hex25)

# Hexagon membership of each survey
all_surveys <- all_surveys %>%
  st_intersection(hex25)

# ----------------------------------------------------------------
# Summarize number of point counts within each hexagon
# ----------------------------------------------------------------


nPC <- all_surveys %>%
  as.data.frame() %>%
  group_by(hexID)%>%
  summarize(nPC = n())

hex25 <- left_join(hex25, nPC)

ggplot(hex25)+geom_sf(aes(fill = nPC))

# ------------------------------------------------
# Plot uncertainty in prediction (coefficient of variation)
# ------------------------------------------------

colscale_uncertainty <- c("#FEFEFE", "#FFF4B3", "#F5D271", "#F2B647", "#EC8E00", "#CA302A")
colpal_uncertainty <- colorRampPalette(colscale_uncertainty)

cut_levs <- c(-0.1,5,25,50,1000)
cut_levs_labs <- c("1 to 5", 
                   "6 to 25", 
                   "26 to 50", 
                   ">50")

hex25$levs <- cut(as.data.frame(hex25)[,"nPC"], 
                              cut_levs,labels=cut_levs_labs)

raster_CV <- stars::st_rasterize(ONGrid_species %>% dplyr::select(CV_levs, x),
                                 nx = dim(raster_q50)[1],ny = dim(raster_q50)[2])


plot1 <- ggplot() +
  
  geom_stars(data = raster_CV) +
  scale_fill_manual(name = "<span style='font-size:13pt'>Relative Uncertainty</span><br><span style='font-size:7pt'>Per 5-minute point count</span><br><span style='font-size:7pt'>Coefficient of Variation</span>",
                    values = colpal_uncertainty(length(levels(raster_CV$CV_levs))), drop=FALSE,na.translate=FALSE)+
  
  # Surveyed squares
  geom_sf(data = subset(ONSquares_centroids, !is.na(PC_detected)), col = "gray40", size=0.4,stroke = 0, shape = 16)+
  # Point count detections
  geom_sf(data = subset(ONSquares_centroids, !is.na(PC_detected) & PC_detected > 0), col = "black",size=0.5,stroke = 0, shape = 16)+
  
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
  theme(legend.key = element_rect(fill = "transparent", colour = "transparent"))+
  
  annotate(geom="text",x=1040000,y=1500000, label= paste0(species_label),lineheight = .85,hjust = 0,size=6,fontface =2) +
  annotate(geom="text",x=1050000,y=680000, label= paste0("Prepared on ",Sys.Date()),size=2,lineheight = .75,hjust = 0,color="gray75")+
  
  guides(fill = guide_legend(order = 1), 
         size = guide_legend(order = 2))

png(paste0("../Output/Prediction_Maps/Uncertainty/",sp_code,"_CV.png"), width=8, height=6.5, units="in", res=1000, type="cairo")
print(plot_CV)
dev.off()

# ------------------------------------------------
# Loop through species, fit models, generate maps
# ------------------------------------------------

population_sums <- data.frame()

species_to_model <- species_to_model %>%
  arrange(desc(n_squares),desc(n_detections))

for (sp_code in species_to_model$Species_Code_BSC){
  
  print(sp_code)
  
  map_file <- paste0("../Output/Prediction_Maps/Relative_Abundance/",sp_code,"_q50.png")
  
  # Skip this species if already run
  #if (file.exists(map_file)) next
  
  # Prepare data for this species
  sp_dat <- all_surveys %>% 
    mutate(count = full_count_matrix[,sp_code]) %>%
    
    # Only select point counts
    subset(Survey_Type %in% c("Point_Count","ARU_SPT","ARU_SPM"))
  
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
  
  # # ------------------------------------------------
  # # Prepare iid random effects
  # # ------------------------------------------------
  # 
  # # Define 'square ID' covariate which will be treated as a random effect
  # sp_dat$sq_idx <- factor(sp_dat$sq_id) %>% as.numeric()
  # 
  # # Define 'square-day' covariate which will be treated as a random effect
  # sp_dat$square_day <- paste0(sp_dat$sq_id,"-",yday(sp_dat$Date_Time)) %>% factor() %>% as.numeric()
  # 
  # # Define unique survey location / year covariate
  # coords <- st_coordinates(sp_dat) %>% round() %>% as.data.frame()
  # X <- coords$X
  # Y <- coords$Y
  # Year <- year(sp_dat$Date_Time)
  # sp_dat$loc_year <- paste0(Year,"-",X,"-",Y) %>% factor() %>% as.numeric()
  # 
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
  
  # ****************************************************************************
  # ****************************************************************************
  # FIT MODEL WITH INLA
  # ****************************************************************************
  # ****************************************************************************
  
  covariates_to_include <- paste0("PC",1:8)
  
  # ------------------------------------------------
  # Create a spatial mesh, which is used to fit the residual spatial field
  # ------------------------------------------------
  
  # Note: mesh developed using tutorial at: https://rpubs.com/jafet089/886687
  max.edge = diff(range(st_coordinates(sp_dat)[,1]))/10
  bound.outer = diff(range(st_coordinates(sp_dat)[,1]))/3
  cutoff = max.edge/5
  bound.outer = diff(range(st_coordinates(sp_dat)[,1]))/3
  
  mesh_spatial <- inla.mesh.2d(loc = st_coordinates(sp_dat),
                               cutoff = max.edge/2,
                               max.edge = c(1,2)*max.edge,
                               offset=c(max.edge, bound.outer))
  mesh_locs <- mesh_spatial$loc[,c(1,2)] %>% as.data.frame()
  
  prior_range <- c(300000,0.1) # 10% chance range is smaller than 300000
  prior_sigma <- c(0.5,0.1) # 10% chance sd is larger than 0.5
  matern_coarse <- inla.spde2.pcmatern(mesh_spatial,
                                       prior.range = prior_range, 
                                       prior.sigma = prior_sigma 
  )
  
  # ------------------------------------------------
  # Create mesh to model effect of time since sunrise (TSS)
  # ------------------------------------------------
  
  TSS_range <- range(sp_dat$Hours_Since_Sunrise)
  TSS_meshpoints <- seq(TSS_range[1]-0.1,TSS_range[2]+0.1,length.out = 11)
  TSS_mesh1D = inla.mesh.1d(TSS_meshpoints,boundary="free")
  TSS_spde = inla.spde2.pcmatern(TSS_mesh1D,
                                 prior.range = c(6,0.1),
                                 prior.sigma = c(1,0.1)) # 10% chance sd is larger than 1
  
  # ------------------------------------------------
  # Model formulas
  # ------------------------------------------------
  
  sd_linear <- 0.2
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
  
  # ****************************************************************************
  # ****************************************************************************
  # GENERATE MAPS
  # ****************************************************************************
  # ****************************************************************************
  
  # For every pixel on landscape, extract distance from eBird range limit
  ONGrid_species <- ONGrid %>%
    mutate(distance_from_range = (st_centroid(.) %>% 
                                    st_distance( . , range) %>% 
                                    as.numeric())/1000)
  
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
  start2 <- Sys.time()
  pred <- NULL
  pred <- generate(fit_INLA,
                   as(ONGrid_species,'Spatial'),
                   formula =  pred_formula_PC,
                   n.samples = 1000)
  
  pred <- exp(pred)
  
  # Median and upper/lower credible intervals (90% CRI)
  prediction_quantiles = apply(pred,1,function(x) quantile(x,c(0.05,0.5,0.95),na.rm = TRUE))
  ONGrid_species$pred_q05 <- prediction_quantiles[1,]
  ONGrid_species$pred_q50 <- prediction_quantiles[2,]
  ONGrid_species$pred_q95 <- prediction_quantiles[3,]
  ONGrid_species$pred_CI_width_90 <- prediction_quantiles[3,] - prediction_quantiles[1,]
  ONGrid_species$CV <- apply(pred,1,function(x) sd(x,na.rm = TRUE)/mean(x,na.rm = TRUE))
  
  # Probability of observing species in 5-minute point count
  size <- fit_INLA$summary.hyperpar$'0.5quant'[1] # parameter of negative binomial
  
  # Probability of detecting species in a 5-minute point count
  prob_zero_PC <- dnbinom(0,mu=prediction_quantiles[2,],size=size)
  ONGrid_species$pObs_5min <- 1-prob_zero_PC
  
  # Convert to CRS of target raster
  ONGrid_species <- ONGrid_species %>% st_transform(st_crs(target_raster))
  
  end2 <- Sys.time() 
  runtime_pred <- difftime( end2,start2, units="mins") %>% round(2)
  print(paste0(sp_code," - ",runtime_pred," min to generate predictions"))
  
  # ------------------------------------------------
  # Summarize ONSquares where species was detected
  # ------------------------------------------------
  
  PC_detected <- sp_dat %>%
    subset(Survey_Type %in% c("Point_Count","ARU_SPT","ARU_SPM")) %>% 
    st_intersection(ONSquares)%>%
    as.data.frame() %>%
    group_by(sq_id) %>%
    summarize(PC_detected = as.numeric(sum(count)>0),
              PC_mean_count = mean(count) %>% round(2))
  
  # CL_detected <-sp_dat %>%
  #   subset(Survey_Type %in% c("Breeding Bird Atlas","Linear transect")) %>%
  #   as.data.frame() %>%
  #   group_by(sq_id) %>%
  #   summarize(CL_detected = as.numeric(sum(count)>0),
  #             CL_mean_count = mean(count))
  
  ONSquares_species <- ONSquares %>%
    relocate(geometry,.after = last_col()) %>%
    left_join(PC_detected) #%>% left_join(CL_detected)
  
  ONSquares_centroids <- st_centroid(ONSquares_species)
  
  # ------------------------------------------------
  # Label for figure and ebird range limit
  # ------------------------------------------------
  
  species_name = ON_spcd$CommonName[which(ON_spcd$spcd == sp_code)]
  species_label = ON_spcd$Label[which(ON_spcd$spcd == sp_code)]
  
  # ------------------------------------------------
  # sf object for ebird range limit (optional - not needed for plotting)
  # ------------------------------------------------
  
  range <- NA
  if (sp_code %in% names(species_ranges)){
    range <- species_ranges[[sp_code]]  %>% 
      st_transform(st_crs(ONBoundary)) %>%
      st_intersection(ONBoundary)
  }
  
  # ------------------------------------------------
  # Plot median prediction
  # ------------------------------------------------
  
  colscale_relabund <-c("#FEFEFE", "#FBF7E2", "#FCF8D0", "#EEF7C2", "#CEF2B0", "#94E5A0", "#51C987", "#18A065", "#008C59", "#007F53", "#006344")
  colpal_relabund <- colorRampPalette(colscale_relabund)

  lower_bound <- 0.01
  upper_bound <- quantile(ONGrid_species$pred_q50,0.99,na.rm = TRUE) %>% signif(2)
  if (lower_bound >= (upper_bound/5)) lower_bound <- (upper_bound/5) %>% signif(2)
  
  sp_cut <- cut.fn(df = ONGrid_species,
                   target_raster = target_raster,
                   column_name = "pred_q50",
                   lower_bound = lower_bound,
                   upper_bound = upper_bound)
  
  raster_q50 <- sp_cut$raster
  
  # Legend for size
  size_breaks <- c(0,mean(ONSquares_species$PC_mean_count,na.rm = TRUE),quantile(ONSquares_species$PC_mean_count,0.99,na.rm = TRUE),max(ONSquares_species$PC_mean_count,na.rm = TRUE)) %>% round(1) %>% as.numeric() %>% sort()
  
  # Median of posterior
  plot_q50 <- ggplot() +
    
    geom_stars(data = raster_q50) +
    scale_fill_manual(name = "<span style='font-size:13pt'>Relative Abundance</span><br><span style='font-size:7pt'>Per 5-minute point count</span><br><span style='font-size:7pt'>(Posterior Median)</span>",
                      values = colpal_relabund(length(levels(raster_q50$levs))), drop=FALSE,na.translate=FALSE)+
    
    
    # Surveyed squares
    geom_sf(data = subset(ONSquares_centroids, !is.na(PC_detected)), col = "gray70", size=0.4,stroke = 0, shape = 16)+
    
    # Point count detections
    geom_sf(data = subset(ONSquares_centroids, !is.na(PC_detected) & PC_detected > 0), col = "black",size=0.5,stroke = 0, shape = 16)+
    
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
    theme(legend.key = element_rect(fill = "transparent", colour = "transparent"))+
    
   annotate(geom="text",x=1050000,y=1500000, label= paste0(species_label),lineheight = .85,hjust = 0,size=6,fontface =2) +
    annotate(geom="text",x=1050000,y=680000, label= paste0("Prepared on ",Sys.Date()),size=2,lineheight = .75,hjust = 0,color="gray75")+
    
    guides(fill = guide_legend(order = 1), 
           size = guide_legend(order = 2))
  
  png(map_file, width=8, height=6.5, units="in", res=1000, type="cairo")
  print(plot_q50)
  dev.off()
  
  # ------------------------------------------------
  # Plot uncertainty in prediction (width of 90% CRI)
  # ------------------------------------------------
  
  colscale_uncertainty <- c("#FEFEFE", "#FFF4B3", "#F5D271", "#F2B647", "#EC8E00", "#CA302A")
  colpal_uncertainty <- colorRampPalette(colscale_uncertainty)
  
  lower_bound <- 0.01
  upper_bound <- quantile(ONGrid_species$pred_CI_width_90,0.99,na.rm = TRUE) %>% signif(2)
  if (lower_bound >= (upper_bound/5)) lower_bound <- (upper_bound/5) %>% signif(2)
  
  raster_CI_width_90 <- cut.fn(df = ONGrid_species,
                               target_raster = target_raster,
                               column_name = "pred_CI_width_90",
                               lower_bound = lower_bound,
                               upper_bound = upper_bound)$raster
  
  plot_CI_width_90 <- ggplot() +
    
    geom_stars(data = raster_CI_width_90) +
    scale_fill_manual(name = "<span style='font-size:13pt'>Relative Uncertainty</span><br><span style='font-size:7pt'>Per 5-minute point count</span><br><span style='font-size:7pt'>Width of 90% CI</span>",
                      values = colpal_uncertainty(length(levels(raster_CI_width_90$levs))), drop=FALSE,na.translate=FALSE)+
    
    # Surveyed squares
    geom_sf(data = subset(ONSquares_centroids, !is.na(PC_detected)), col = "gray40", size=0.4,stroke = 0, shape = 16)+
    # Point count detections
    geom_sf(data = subset(ONSquares_centroids, !is.na(PC_detected) & PC_detected > 0), col = "black",size=0.5,stroke = 0, shape = 16)+
    
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
    theme(legend.key = element_rect(fill = "transparent", colour = "transparent"))+
    
    annotate(geom="text",x=1040000,y=1500000, label= paste0(species_label),lineheight = .85,hjust = 0,size=6,fontface =2) +
    annotate(geom="text",x=1050000,y=680000, label= paste0("Prepared on ",Sys.Date()),size=2,lineheight = .75,hjust = 0,color="gray75")+
    
    guides(fill = guide_legend(order = 1), 
           size = guide_legend(order = 2))
  
  png(paste0("../Output/Prediction_Maps/Uncertainty/",sp_code,"_CI_width_90.png"), width=8, height=6.5, units="in", res=1000, type="cairo")
  print(plot_CI_width_90)
  dev.off()
  
  # ------------------------------------------------
  # Plot uncertainty in prediction (coefficient of variation)
  # ------------------------------------------------
  
  colscale_uncertainty <- c("#FEFEFE", "#FFF4B3", "#F5D271", "#F2B647", "#EC8E00", "#CA302A")
  colpal_uncertainty <- colorRampPalette(colscale_uncertainty)
  
  cut_levs <- c(-0.1,0.1,0.25,0.5,1,2,2000)
  cut_levs_labs <- c("0 to 0.10", 
                     "0.10 to 0.25", 
                     "0.25 to 0.50", 
                     "0.50 to 1", 
                     "1 to 2", 
                     "> 2")
  
  ONGrid_species$CV_levs <- cut(as.data.frame(ONGrid_species)[,"CV"], 
                                  cut_levs,labels=cut_levs_labs)
  raster_CV <- stars::st_rasterize(ONGrid_species %>% dplyr::select(CV_levs, x),
                                   nx = dim(raster_q50)[1],ny = dim(raster_q50)[2])
  
  
  plot_CV <- ggplot() +
    
    geom_stars(data = raster_CV) +
    scale_fill_manual(name = "<span style='font-size:13pt'>Relative Uncertainty</span><br><span style='font-size:7pt'>Per 5-minute point count</span><br><span style='font-size:7pt'>Coefficient of Variation</span>",
                      values = colpal_uncertainty(length(levels(raster_CV$CV_levs))), drop=FALSE,na.translate=FALSE)+
    
    # Surveyed squares
    geom_sf(data = subset(ONSquares_centroids, !is.na(PC_detected)), col = "gray40", size=0.4,stroke = 0, shape = 16)+
    # Point count detections
    geom_sf(data = subset(ONSquares_centroids, !is.na(PC_detected) & PC_detected > 0), col = "black",size=0.5,stroke = 0, shape = 16)+
    
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
    theme(legend.key = element_rect(fill = "transparent", colour = "transparent"))+
    
    annotate(geom="text",x=1040000,y=1500000, label= paste0(species_label),lineheight = .85,hjust = 0,size=6,fontface =2) +
    annotate(geom="text",x=1050000,y=680000, label= paste0("Prepared on ",Sys.Date()),size=2,lineheight = .75,hjust = 0,color="gray75")+
    
    guides(fill = guide_legend(order = 1), 
           size = guide_legend(order = 2))
  
  png(paste0("../Output/Prediction_Maps/Uncertainty/",sp_code,"_CV.png"), width=8, height=6.5, units="in", res=1000, type="cairo")
  print(plot_CV)
  dev.off()
  
  # ------------------------------------------------
  # Plot probability of observing species in a 5-minute point count
  # ------------------------------------------------
  
  colscale_pObs <- c("#FEFEFE",RColorBrewer::brewer.pal(5,"BuGn")[2:5])
  colpal_pObs <- colorRampPalette(colscale_pObs)
  
  cut_levs <- c(-0.1,0.01,0.05,0.125,0.5,1)
  cut_levs_labs <- c("0 to 0.01", 
                     "0.01 to 0.05", 
                     "0.05 to 0.125", 
                     "0.125 to 0.50", 
                     "0.50 to 1")
  
  ONGrid_species$pObs_levs <- cut(as.data.frame(ONGrid_species)[,"pObs_5min"], 
                                    cut_levs,labels=cut_levs_labs)
  
  raster_pObs = stars::st_rasterize(ONGrid_species %>% dplyr::select(pObs_levs, x),
                                    nx = dim(raster_q50)[1],ny = dim(raster_q50)[2])
  
  # Median of posterior
  plot_pObs <- ggplot() +
    
    geom_stars(data = raster_pObs) +
    scale_fill_manual(name = "<span style='font-size:13pt'>Prob. of Observation</span><br><span style='font-size:7pt'>Per 5-minute point count</span><br><span style='font-size:7pt'>(Posterior Median)</span>",
                      values = colpal_pObs(length(levels(raster_pObs$pObs_levs))), 
                      drop=FALSE,na.translate=FALSE)+
    
    # Surveyed squares
    geom_sf(data = subset(ONSquares_centroids, !is.na(PC_detected)), col = "gray40", size=0.4,stroke = 0, shape = 16)+
    # Point count detections
    geom_sf(data = subset(ONSquares_centroids, !is.na(PC_detected) & PC_detected > 0), col = "black",size=0.5,stroke = 0, shape = 16)+
    
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
    theme(legend.key = element_rect(fill = "transparent", colour = "transparent"))+
    
    annotate(geom="text",x=1040000,y=1500000, label= paste0(species_label),lineheight = .85,hjust = 0,size=6,fontface =2) +
    annotate(geom="text",x=1050000,y=680000, label= paste0("Prepared on ",Sys.Date()),size=2,lineheight = .75,hjust = 0,color="gray75")+
    
    guides(fill = guide_legend(order = 1), 
           size = guide_legend(order = 2))
  
  print(plot_pObs)
  
  png(paste0("../Output/Prediction_Maps/PObs/",sp_code,"_PObs.png"), width=8, height=6.5, units="in", res=1000, type="cairo")
  print(plot_pObs)
  dev.off()
  
  # # ------------------------------------------------
  # # Fit BRT for this species, for comparison
  # # ------------------------------------------------
  # 
  # PC_df <- as.data.frame(PC_sp)
  # PC_df$in_range <- (!(PC_df$distance_from_range>0)) %>% as.numeric()
  # 
  # BRT_covariates <- as.data.frame(PC_df) %>%
  #   dplyr::select(elevation_1km:PC10) %>%
  #   colnames()
  # BRT_covariates <- c(BRT_covariates,"Hours_Since_Sunrise","distance_from_range","in_range")
  # 
  # fit_BRT <- fit_brt(model_data =  PC_df,
  #                    response_column = which(colnames(PC_df)=="count"),
  #                    covariate_columns = which(colnames(PC_df)%in% BRT_covariates))
  # 
  # # Predictions across ON
  # ONGrid_df <- ONGrid %>%
  #   mutate(distance_from_range = (st_centroid(.) %>%
  #                                   st_distance( . , range) %>%
  #                                   as.numeric())/1000) %>%
  #   as.data.frame()
  # ONGrid_df$in_range <- (!(ONGrid_df$distance_from_range>0)) %>% as.numeric()
  # ONGrid_df$Hours_Since_Sunrise <- 0
  # pred_BRT  <- predict(fit_BRT,
  #                      ONGrid_df,
  #                      n.trees = fit_BRT$gbm.call$best.trees)
  # 
  # ONGrid_species$pred_BRT <- exp(pred_BRT + log_offset_5min)
  # 
  # lower_bound <- 0.01
  # upper_bound <- quantile(ONGrid_species$pred_q50,0.99,na.rm = TRUE) %>% signif(2)
  # if (lower_bound >= (upper_bound/5)) lower_bound <- (upper_bound/5) %>% signif(2)
  # 
  # raster_BRT <- cut.fn(df = ONGrid_species,
  #                      target_raster = target_raster,
  #                      column_name = "pred_BRT",
  #                      lower_bound = lower_bound,
  #                      upper_bound = upper_bound)$raster
  # 
  # # Median of posterior
  # plot_BRT <- ggplot() +
  # 
  #   geom_stars(data = raster_BRT) +
  #   scale_fill_manual(name = "<span style='font-size:13pt'>Relative Abundance</span><br><span style='font-size:7pt'>Per 5-minute point count</span><br><span style='font-size:7pt'>BRT model</span>",
  #                     values = colpal_relabund(length(levels(raster_BRT$levs))), drop=FALSE,na.translate=FALSE)+
  # 
  #   geom_sf(data = ONWater,colour=NA,fill="#59F3F3",show.legend = F)+
  # 
  #   geom_sf(data = ONBoundary,colour="black",fill=NA,lwd=0.3,show.legend = F) +
  # 
  #   #geom_sf(data = range,colour="darkred",fill=NA,show.legend = F) +
  # 
  #   # Surveyed squares
  #   geom_sf(data = subset(ONSquares_centroids, !is.na(PC_detected) | !is.na(CL_detected)), col = "gray70", pch = 19, size = 0.1)+
  # 
  #   # Point count detections
  #   geom_sf(data = subset(ONSquares_centroids, !is.na(PC_detected) & PC_detected > 0), col = "black", pch = 19, size = 0.2)+
  # 
  #   # Checklist detections
  #   #geom_sf(data = subset(ONSquares_centroids, !is.na(CL_detected) & CL_detected > 0), col = "black", pch = 4, size = 0.2)+
  # 
  #   coord_sf(clip = "off",xlim = range(as.data.frame(st_coordinates(ONBoundary))$X))+
  #   theme(panel.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  #   theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
  #   theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
  #   theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
  #   theme(legend.margin=margin(0,0,0,0),legend.box.margin=margin(5,10,5,-20),legend.title.align=0.5,
  #         legend.title = element_markdown(lineheight=.9,hjust = "left"))+
  #   annotate(geom="text",x=346000,y=1850000, label= paste0(species_label),lineheight = .85,hjust = 0,size=6,fontface =2) +
  #   annotate(geom="text",x=346000,y=1400000, label= paste0("Prepared on ",Sys.Date()),size=3,lineheight = .75,hjust = 0,color="#3b3b3b")
  # 
  # png(paste0("../Output/Prediction_Maps/Relative_Abundance/",sp_code,"_BRT.png"), width=6.5, height=8, units="in", res=300, type="cairo")
  # print(plot_BRT)
  # dev.off()
  # 
  # # ----------------------------------------------------------------
  # # Save estimates of relative population size
  # # ----------------------------------------------------------------
  # sum_INLA <- apply(pred,2,sum) %>% quantile(c(0.025,0.5,0.975)) %>% as.numeric()
  # sum_BRT <- sum(ONGrid_species$pred_BRT)
  # 
  # file_path <- "../Output/Population_Sums.csv"
  # if (file.exists(file_path)) population_sums <- read.csv(file_path)
  # population_sums <- rbind(population_sums,
  #                          data.frame(Species = sp_code,
  #                                     sum_BRT = sum_BRT,
  #                                     sum_INLA_q025 = sum_INLA[1],
  #                                     sum_INLA_q500 = sum_INLA[2],
  #                                     sum_INLA_q975 = sum_INLA[3]))
  # write.csv(population_sums,file = file_path, row.names = FALSE)
  # 
  # popsum_plot <- ggplot()+
  #   geom_errorbarh(data = population_sums, aes(y = Species, xmin = sum_INLA_q025, xmax = sum_INLA_q975), height = 0)+
  #   geom_point(data = population_sums, aes(y = Species, x = sum_INLA_q500))+
  #   geom_point(data = population_sums, aes(y = Species, x = sum_BRT), col = "blue", pch = 4)+
  #   scale_x_continuous(trans = "log10")+
  #   theme_bw()
  # print(popsum_plot)
  
} # close species loop


