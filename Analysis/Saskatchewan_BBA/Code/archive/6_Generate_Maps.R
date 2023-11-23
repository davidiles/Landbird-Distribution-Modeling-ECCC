# ************************************************
# GENERATE PREDICTIONS AND MAPS OF EXPECTED COUNTS
#  - requires model to have been fit by previous script
# ************************************************

# ------------------------------------------------
# Detach any loaded packages
# ------------------------------------------------

lapply(names(sessionInfo()$loadedOnly), require, character.only = TRUE)
invisible(lapply(paste0('package:', names(sessionInfo()$otherPkgs)), detach, character.only=TRUE, unload=TRUE, force=TRUE))

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
# Load packages
# ------------------------------------------------

require(tidyverse)
require(sf)

# Still uses some 'raster' and 'sp' functionality - needs to be fixed
require(terra)
#require(raster)
#require(fasterize)

# For plotting
require(ggtext)

require(INLA)     # INLA_22.05.07  
require(inlabru)  # inlabru_2.7.0

# ------------------------------------------------
# Load data
# ------------------------------------------------

# Raster with target properties
SaskRaster <- rast("../../../Data/Spatial/Saskatchewan/SaskGrid/SaskRaster.tif")
#target_raster <- raster(SaskRaster)

SaskBoundary <- SK_BCR <- st_read("../../../Data/Spatial/National/BCR/BCR_Terrestrial_master.shp")  %>%
  subset(PROVINCE_S == "SASKATCHEWAN") %>%
  st_make_valid() %>%
  st_union() %>%
  st_transform(st_crs(SaskRaster))

SaskWater <- read_sf("../../../Data/Spatial/Saskatchewan/SaskWater/SaskWaterClip.shp") %>% 
  st_transform(st_crs(SaskRaster))
SaskWater$area <- st_area(SaskWater)
SaskWater$area <- as.numeric(SaskWater$area)
SaskWater <- SaskWater[SaskWater$area>2.580e+06 ,]

# Atlas Squares
SaskSquares <- st_read("../../../Data/Spatial/Saskatchewan/SaskSquares/SaskSquares.shp") %>%
  st_transform(st_crs(all_surveys)) %>%
  dplyr::select(SQUARE_ID) %>%
  rename(sq_id = SQUARE_ID)

# Water_centroids <- SaskGrid %>%
#   dplyr::select(pointid) %>%
#   st_centroid() %>%
#   st_intersection(Water_centroids,SaskWater)
# saveRDS(Water_centroids,"!Data/!Processed/Water_centroids.rds")
# Water_centroids <- readRDS("!Data/!Processed/Water_centroids.rds")

# ------------------------------------------------
# Function to generate predictions for 1 km x 1 km grid cells (SaskGrid)
# ------------------------------------------------

pred_fn <- function(sp_code, sp_dat = sp_dat, path, nsamples){
  
  SaskGrid_species <- SaskGrid
  
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
  
  # For every pixel on landscape, extract distance from eBird range limit
  SaskGrid_species <- SaskGrid_species %>%
    mutate(distance_from_range = (st_centroid(SaskGrid_species) %>% 
                                    st_distance( . , range) %>% 
                                    as.numeric())/1000)
  
  # --------------------------------
  # QPAD offsets associated with a 5-minute unlimited distance survey
  # --------------------------------
  
  species_offsets <- subset(species_to_model, Species_Code_BSC == sp_code)
  
  log_offset_5min <- 0
  if (species_offsets$offset_exists == TRUE) log_offset_5min <- species_offsets$log_offset_5min
  
  # ------------------------------------------------
  # Generate predictions on SaskGrid_species raster (1 km x 1 km pixels)
  # ------------------------------------------------
  
  pred_formula_PC = as.formula(paste0(' ~

                  Intercept_PC +
                  spde_coarse +
                  range_effect * distance_from_range + 
                   ',
                                      paste0("Beta1_",covariates_to_include,'*',covariates_to_include, collapse = " + ")))
  
  # Predictions are on log scale, and do not include variance components
  pred <- generate(fit,
                   as(SaskGrid_species,'Spatial'),
                   formula =  pred_formula_PC,
                   n.samples = nsamples)
  
  # Convert to count scale and add log-normal variance component corrections
  Var_sq_idx <- 1/summary(fit)$inla$hyperpar["Precision for kappa_squareID","0.5quant"]
  
  pred <- exp(pred + log_offset_5min + 0.5*Var_sq_idx)
  
  # Median and upper/lower credible intervals
  SaskGrid_species$pred_q50 <- apply(pred,1,function(x) quantile(x,0.5,na.rm = TRUE))
  SaskGrid_species$pred_q05 <- apply(pred,1,function(x) quantile(x,0.05,na.rm = TRUE))
  SaskGrid_species$pred_q95 <- apply(pred,1,function(x) quantile(x,0.95,na.rm = TRUE))
  
  SaskGrid_species <- SaskGrid_species %>% 
    dplyr::select(pred_q50,pred_q05,pred_q95,distance_from_range)
  
  return(list(pred = pred, SaskGrid_species = SaskGrid_species))
}

# ------------------------------------------------
# Colour scales for plotting
# ------------------------------------------------

colscale_density <- c("#FEFEFE", "#FFF4B3", "#F5D271", "#F2B647", "#EC8E00", "#CA302A")
colpal_density <- colorRampPalette(colscale_density)

colscale_relabund <-c("#FEFEFE", "#FBF7E2", "#FCF8D0", "#EEF7C2", "#CEF2B0", "#94E5A0", "#51C987", "#18A065", "#008C59", "#007F53", "#006344")
colpal_relabund <- colorRampPalette(colscale_relabund)

# ------------------------------------------------
# Function to rasterize a series of spatial predictions (needed for plotting)
# ------------------------------------------------

cut.fn <- function(df, 
                   target_raster, 
                   region_boundary,
                   column_name, 
                   lower_bound = NA, 
                   upper_bound = NA){
  
  max_val <- upper_bound
  max_val <- ifelse(is.na(max_val), 0, max_val)
  
  max_lev <- ifelse(max_val > 1.6, 4,
                    ifelse(max_val > 0.8, 4, 3))
  
  cut_levs <- signif(max_val/(2^((max_lev-1):0)), 2)
  cut_levs <- unique(cut_levs)
  cut_levs <- ifelse(is.na(cut_levs), 0, cut_levs)
  
  if (lower_bound %in% cut_levs) cut_levs <- cut_levs[-which(cut_levs == lower_bound)]
  if (lower_bound > min(cut_levs)) cut_levs = cut_levs[-which(cut_levs < lower_bound)]
  
  max_lev <- length(cut_levs) ## do this because sometimes there are fewer levels
  
  cut_levs_labs <- c(paste0("0-",lower_bound),
                     paste(lower_bound, cut_levs[1], sep="-"),
                     paste(cut_levs[-max_lev], cut_levs[-1], sep="-"),
                     paste(cut_levs[max_lev], "+"))
  
  cut_levs <- c(-1, lower_bound, cut_levs, 1000) %>% unique()
  
  # ** MULTIPLY BY 15 FOR COMPARISON TO BSC
  df$pc_cut <- cut(as.data.frame(df)[,column_name], cut_levs, labels=cut_levs_labs, ordered=TRUE)
  
  sp_abund_raster <- terra:rasterize(df,SaskRaster,field = "pc_cut")
  sp_abund_raster <- terra::mask(sp_abund_raster,as(SaskBoundary,'Spatial'))
  
  # Convert to dataframe for plotting with ggplot
  sp_abund_raster_data <- as.data.frame(raster::rasterToPoints(sp_abund_raster))
  sp_abund_raster_data$layer[sp_abund_raster_data$layer==1] <- cut_levs_labs[1]
  sp_abund_raster_data$layer[sp_abund_raster_data$layer==2] <- cut_levs_labs[2]
  sp_abund_raster_data$layer[sp_abund_raster_data$layer==3] <- cut_levs_labs[3]
  sp_abund_raster_data$layer[sp_abund_raster_data$layer==4] <- cut_levs_labs[4]
  sp_abund_raster_data$layer[sp_abund_raster_data$layer==5] <- cut_levs_labs[5]
  sp_abund_raster_data$layer[sp_abund_raster_data$layer==6] <- cut_levs_labs[6]
  
  sp_abund_raster_data$layer <- factor(sp_abund_raster_data$layer,
                                       levels= c(cut_levs_labs[1],cut_levs_labs[2],cut_levs_labs[3],
                                                 cut_levs_labs[4],cut_levs_labs[5],cut_levs_labs[6]),
                                       ordered = T)
  return(sp_abund_raster_data)
}

# ------------------------------------------------
# Loop through species and generate predictions
# ------------------------------------------------

for (sp_code in species_to_model$Species_Code_BSC){
  
  # ------------------------------------------------
  # Skip species if model was not fit
  # ------------------------------------------------
  model_path <- paste0("../Output/Fitted_Models/Model_",sp_code,".rds")
  if (!file.exists(model_path)) next
  
  # ------------------------------------------------
  # Skip species if map is already generated
  # ------------------------------------------------
  
  #if (file.exists(paste0("../Output/Prediction_Maps/q50/",sp_code,"_q50.png"))) next
  
  # ------------------------------------------------
  # Extract fitted model for this species (and associated data that went into model)
  # ------------------------------------------------
  
  species_fit <- readRDS(paste0("../Output/Fitted_Models/Model_",sp_code,".rds"))
  
  sp_dat <- species_fit$sp_dat
  fit <- species_fit$fit_PConly
  
  # ------------------------------------------------
  # Generate predictions
  # ------------------------------------------------
  
  print(sp_code)
  
  start <- Sys.time()
  pred_sp <- pred_fn(sp_code = sp_code, 
                     sp_dat = sp_dat,
                     nsamples = 1)
  end <- Sys.time() # 4 min for 250 samples
  end-start
  
  # ------------------------------------------------
  # Create color scale limits
  # ------------------------------------------------
  
  lower_bound <- 0.01
  upper_bound <- quantile(pred_sp$SaskGrid_species$pred_q50,0.99,na.rm = TRUE) %>% signif(2)
  if (lower_bound >= (upper_bound/5)){
    upper_bound <- quantile(pred_sp$SaskGrid_species$pred_q50,0.99,na.rm = TRUE) %>% signif(2)
    lower_bound <- (upper_bound/5) %>% signif(2)
  }
  
  # ------------------------------------------------
  # Generate prediction rasters (use same scale limits for all 3 to visualize differences)
  # ------------------------------------------------
  
  pred_sp$SaskGrid_species <- pred_sp$SaskGrid_species %>% st_transform(st_crs(SaskRaster))
  raster_q50 <- cut.fn(df = pred_sp$SaskGrid_species, 
                       target_raster = SaskRaster, 
                       region_boundary = SaskBoundary,
                       column_name = "pred_q50", 
                       lower_bound = lower_bound, 
                       upper_bound = upper_bound)
  
  raster_q05 <- cut.fn(df = pred_sp$SaskGrid_species, 
                       target_raster = target_raster, 
                       region_boundary = SaskBoundary,
                       column_name = "pred_q05", 
                       lower_bound = lower_bound, 
                       upper_bound = upper_bound)
  
  raster_q95 <- cut.fn(df = pred_sp$SaskGrid_species, 
                       target_raster = target_raster, 
                       region_boundary = SaskBoundary,
                       column_name = "pred_q95", 
                       lower_bound = lower_bound, 
                       upper_bound = upper_bound)
  
  # ------------------------------------------------
  # Summarize SaskSquares where species was detected
  # ------------------------------------------------
  
  PC_detected <- sp_dat %>%
    subset(Survey_Type %in% c("Point_Count","ARU_SPT","ARU_SPM")) %>% 
    as.data.frame() %>%
    group_by(sq_id) %>%
    summarize(PC_detected = as.numeric(sum(count)>0))
  
  CL_detected <-sp_dat %>%
    subset(Survey_Type %in% c("Breeding Bird Atlas","Linear transect")) %>%
    as.data.frame() %>%
    group_by(sq_id) %>%
    summarize(CL_detected = as.numeric(sum(count)>0))
  
  SaskSquares_species <- SaskSquares %>%
    relocate(geometry,.after = last_col()) %>%
    left_join(PC_detected) %>% 
    left_join(CL_detected)
  
  SaskSquares_centroids <- st_centroid(SaskSquares_species)
  
  # ------------------------------------------------
  # Generate maps
  # ------------------------------------------------
  
  species_name = Sask_spcd$CommonName[which(Sask_spcd$spcd == sp_code)]
  species_label = Sask_spcd$Label[which(Sask_spcd$spcd == sp_code)]
  
  range <- NA
  if (sp_code %in% names(species_ranges)){
    range <- species_ranges[[sp_code]]  %>% 
      st_transform(st_crs(SaskBoundary)) %>%
      st_intersection(SaskBoundary)
  }
  
  # Median of posterior
  pred_q50 <- ggplot() +
    
    geom_raster(data = raster_q50 , aes(x = x, y = y, fill = layer)) +
    scale_fill_manual(name = "<span style='font-size:13pt'>Relative Abundance</span><br><span style='font-size:7pt'>Per 5-minute point count</span><br><span style='font-size:7pt'>(Posterior Median)</span>",
                      values = colpal_relabund(length(levels(raster_q50$layer))), drop=FALSE)+
    
    geom_sf(data = SaskWater,colour=NA,fill="#59F3F3",show.legend = F)+
    
    geom_sf(data = SaskBoundary,colour="black",fill=NA,lwd=0.3,show.legend = F) +
    geom_sf(data = range,colour="darkred",fill=NA,show.legend = F) +
    
    # Surveyed squares
    geom_sf(data = subset(SaskSquares_centroids, !is.na(PC_detected) | !is.na(CL_detected)), col = "gray70", pch = 19, size = 0.1)+
    
    # Point count detections
    geom_sf(data = subset(SaskSquares_centroids, !is.na(PC_detected) & PC_detected > 0), col = "black", pch = 19, size = 0.2)+
    
    # Checklist detections
    #geom_sf(data = subset(SaskSquares_centroids, !is.na(CL_detected) & CL_detected > 0), col = "black", pch = 4, size = 0.2)+
    
    coord_sf(clip = "off",xlim = range(raster_q50$x))+
    theme(panel.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    annotate(geom="text",x=346000,y=960000, label= paste0(species_label),lineheight = .85,hjust = 0,size=6,fontface =2) +
    annotate(geom="text",x=346000,y=550000, label= paste0("Prepared on ",Sys.Date()),size=3,lineheight = .75,hjust = 0,color="#3b3b3b")+
    theme(legend.margin=margin(0,0,0,0),legend.box.margin=margin(5,10,5,-20),legend.title.align=0.5,
          legend.title = element_markdown(lineheight=.9,hjust = "left"))
  
  # 5th percentile of posterior
  pred_q05 <- ggplot() +
    
    geom_raster(data = raster_q05 , aes(x = x, y = y, fill = layer)) +
    scale_fill_manual(name = "<span style='font-size:13pt'>Relative Abundance</span><br><span style='font-size:7pt'>Per 5-minute point count</span><br><span style='font-size:7pt'>(5% CRI)</span>",
                      values = colpal_relabund(length(levels(raster_q50$layer))), drop=FALSE)+
    
    geom_sf(data = SaskWater,colour=NA,fill="#59F3F3",show.legend = F)+
    
    geom_sf(data = SaskBoundary,colour="black",fill=NA,lwd=0.3,show.legend = F) +
    geom_sf(data = range,colour="darkred",fill=NA,show.legend = F) +
    
    # Surveyed squares
    geom_sf(data = subset(SaskSquares_centroids, !is.na(PC_detected) | !is.na(CL_detected)), col = "gray70", pch = 19, size = 0.1)+
    
    # Point count detections
    geom_sf(data = subset(SaskSquares_centroids, !is.na(PC_detected) & PC_detected > 0), col = "black", pch = 19, size = 0.2)+
    
    # Checklist detections
    #geom_sf(data = subset(SaskSquares_centroids, !is.na(CL_detected) & CL_detected > 0), col = "black", pch = 4, size = 0.2)+
    
    coord_sf(clip = "off",xlim = range(raster_q50$x))+
    theme(panel.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    annotate(geom="text",x=346000,y=960000, label= paste0(species_label),lineheight = .85,hjust = 0,size=6,fontface =2) +
    annotate(geom="text",x=346000,y=550000, label= paste0("Prepared on ",Sys.Date()),size=3,lineheight = .75,hjust = 0,color="#3b3b3b")+
    theme(legend.margin=margin(0,0,0,0),legend.box.margin=margin(5,10,5,-20),legend.title.align=0.5,
          legend.title = element_markdown(lineheight=.9,hjust = "left"))
  
  # 95% CI
  pred_q95 <- ggplot() +
    
    geom_raster(data = raster_q95 , aes(x = x, y = y, fill = layer)) +
    scale_fill_manual(name = "<span style='font-size:13pt'>Relative Abundance</span><br><span style='font-size:7pt'>Per 5-minute point count</span><br><span style='font-size:7pt'>(95% CRI)</span>",
                      values = colpal_relabund(length(levels(raster_q50$layer))), drop=FALSE)+
    
    geom_sf(data = SaskWater,colour=NA,fill="#59F3F3",show.legend = F)+
    
    geom_sf(data = SaskBoundary,colour="black",fill=NA,lwd=0.3,show.legend = F) +
    geom_sf(data = range,colour="darkred",fill=NA,show.legend = F) +
    
    # Surveyed squares
    geom_sf(data = subset(SaskSquares_centroids, !is.na(PC_detected) | !is.na(CL_detected)), col = "gray70", pch = 19, size = 0.1)+
    
    # Point count detections
    geom_sf(data = subset(SaskSquares_centroids, !is.na(PC_detected) & PC_detected > 0), col = "black", pch = 19, size = 0.2)+
    
    # Checklist detections
    #geom_sf(data = subset(SaskSquares_centroids, !is.na(CL_detected) & CL_detected > 0), col = "black", pch = 4, size = 0.2)+
    
    coord_sf(clip = "off",xlim = range(raster_q50$x))+
    theme(panel.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    annotate(geom="text",x=346000,y=960000, label= paste0(species_label),lineheight = .85,hjust = 0,size=6,fontface =2) +
    annotate(geom="text",x=346000,y=550000, label= paste0("Prepared on ",Sys.Date()),size=3,lineheight = .75,hjust = 0,color="#3b3b3b")+
    theme(legend.margin=margin(0,0,0,0),legend.box.margin=margin(5,10,5,-20),legend.title.align=0.5,
          legend.title = element_markdown(lineheight=.9,hjust = "left"))
  
  
  png(paste0("../Output/Prediction_Maps/q50/",sp_code,"_q50.png"), width=6.5, height=8, units="in", res=300, type="cairo")
  print(pred_q50)
  dev.off()
  
  png(paste0("../Output/Prediction_Maps/q05/",sp_code,"_q05.png"), width=6.5, height=8, units="in", res=300, type="cairo")
  print(pred_q05)
  dev.off()
  
  png(paste0("../Output/Prediction_Maps/q95/",sp_code,"_q95.png"), width=6.5, height=8, units="in", res=300, type="cairo")
  print(pred_q95)
  dev.off()
  
}
