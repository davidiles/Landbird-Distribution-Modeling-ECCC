# ************************************************
# BAYESIAN ANALYSIS / SPECIES DISTRIBUTION MODELS FOR SASKATCHEWAN BREEDING BIRD ATLAS
# 
# 1) A 'processed data package' for analysis is prepared by previous script
# ************************************************

rm(list=ls())

# ------------------------------------------------
# Set working directory
# ------------------------------------------------

dirname <- "C:/Users/IlesD/OneDrive - EC-EC/Iles/Projects/Landbirds/Landbird-Distribution-Modeling-ECCC/Analysis/Saskatchewan_BBA/Code"
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
library(pROC)

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

# Fix survey type labels
all_surveys$Survey_Type[all_surveys$Survey_Type == "IN_PERSON"] <- "Point_Count"
all_surveys$Survey_Type[all_surveys$Survey_Type %in% c("ARU_SM2","ARU_BAR_LT","ARU_SM_UNKN",
                                                       "ARU_SM4","ZOOM_H2N","ARU_IRIVER_E",
                                                       "ARU_MARANTZ","ARU_IRIVER","ARU_UNKNOWN")] <- "ARU_SPT"

# ------------------------------------------------------------------------
# Remove stationary count checklists that overlap spatially (within 500 m) and temporally (within 2 hours) with point counts
# ------------------------------------------------------------------------

SC_dat <- subset(all_surveys, Survey_Type == "Breeding Bird Atlas")
PC_dat <- subset(all_surveys, Survey_Type %in% c("Point_Count","ARU_SPM","ARU_SPT"))

PC_buff <- st_buffer(PC_dat,500) %>% st_union()
SC_dat$to_evaluate <- FALSE
SC_dat$to_evaluate[which(as.matrix(st_intersects(SC_dat,PC_buff)))] <- TRUE
SC_dat <- subset(SC_dat, to_evaluate)

SC_buff <- st_buffer(SC_dat,500) %>% st_union()
PC_dat$to_evaluate <- FALSE
PC_dat$to_evaluate[which(as.matrix(st_intersects(PC_dat,SC_buff)))] <- TRUE
PC_dat <- subset(PC_dat, to_evaluate)

distance_matrix <- st_distance(SC_dat,PC_dat)

# Flags for which observations to remove
SC_dat$to_remove <- FALSE

# Check for and flag overlap
for (i in 1:nrow(SC_dat)){
  
  obs_to_evaluate <- SC_dat[i,]
  dists <- distance_matrix[i,]
  PC_within_500m <- PC_dat[which(dists <= units::set_units(500,"m")),]
  time_diff <- abs(difftime(obs_to_evaluate$Date_Time, PC_within_500m$Date_Time, units='mins')) %>% as.numeric()
  
  # If the checklist was conducted within 12 hours of the point count, remove it
  if (min(time_diff)<=(60*12)) SC_dat$to_remove[i] <- TRUE
}

# SC to remove (approx 76% of spatially overlapping checklists)
mean(SC_dat$to_remove)
SC_to_remove <- subset(SC_dat,to_remove)$Obs_Index
SC_removed <- all_surveys[SC_to_remove,]





# Example of overlapping data
SC_ex <- SC_removed[2:10,]

PC_ex <- PC_dat %>% st_intersection(.,st_buffer(SC_ex,500))


ggplot()+
  geom_sf(data = st_buffer(SC_ex,1000), col = "transparent", fill = "transparent")+
  geom_sf(data = SC_ex, size = 5)+
  geom_sf(data = PC_ex, col = "red")

SC_summary <- SC_ex %>% 
  as.data.frame() %>% 
  dplyr::select(Survey_Type,survey_ID,Obs_Index,Latitude,Longitude,Date_Time,Survey_Duration_Minutes) %>%
  arrange(Latitude)

PC_summary <- PC_ex %>% 
  as.data.frame() %>% 
  dplyr::select(Survey_Type,survey_ID,Obs_Index,
                Latitude,
                Longitude,
                Date_Time,
                Survey_Duration_Minutes) %>%
  arrange(Latitude)

rownames(SC_summary) <- 1:nrow(SC_summary)
rownames(PC_summary) <- 1:nrow(PC_summary)

SC_summary
PC_summary


# Plot/compare species detections
plot_list <- list()

for (i in 1:nrow(SC_summary)){
  
  sp_SC <- full_count_matrix[SC_summary$Obs_Index[i],]
  sp_PC <- full_count_matrix[PC_summary$Obs_Index[i],]
  
  species_detected <- rbind(sp_SC,sp_PC) %>% colSums()
  species_detected <- species_detected[-which(species_detected == 0)] %>% names()
  
  sp_SC <- sp_SC[species_detected]
  sp_PC <- sp_PC[species_detected]
  
  data_comparison <- data.frame(Species = names(sp_SC),
                                n_SC = sp_SC,
                                n_PC = sp_PC)
  
  survey_plot = ggplot(data = data_comparison)+
    
    geom_segment(aes(y = Species, yend = Species, x = 0, xend = n_PC, col = "Point Count"),
                 size = 4,
                 arrow = arrow(length = unit(0.2, "cm")))+
    
    geom_segment(aes(y = Species, yend = Species, x = 0, xend = n_SC, col = "Checklist"),
                 size = 2,
                 arrow = arrow(length = unit(0.2, "cm")))+
    xlab("Count")+
    scale_color_manual(values=c("black","dodgerblue"), name = "Survey Type", guide = "none")+
    theme_bw()+
    theme(plot.title = element_text(size = 4))+
    ggtitle(paste0("Point count ID = ",PC_summary$survey_ID[i],"\nChecklist ID = ",SC_summary$survey_ID[i],"\n\nLat = ",PC_summary$Latitude[i],"\nLon = ",PC_summary$Longitude[i],"\nDate = ",PC_summary$Date_Time[i]))
    plot_list[[i]] = survey_plot
}

library(ggpubr)

ggarrange(plotlist = plot_list)


# Species lists
sp_SC <- full_count_matrix[SC_summary$Obs_Index,]
sp_PC <- full_count_matrix[PC_summary$Obs_Index,]

tmp <- rbind(sp_SC,sp_PC) %>% colSums()
species_detected <- tmp[-which(tmp == 0)] %>% names()

sp_SC <- sp_SC[,species_detected]
sp_PC <- sp_PC[,species_detected]

sp_SC
sp_PC
