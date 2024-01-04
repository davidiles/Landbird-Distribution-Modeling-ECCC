# ------------------------------------------------
# Load/install packages and set graphical themes / working directory
# ------------------------------------------------

# my_packs = c(
#   'tidyverse',
#              'magrittr',
#              'sf',
#              'terra',
#              'ggspatial',
#              'naturecounts',
#              'suntools',
#              "readxl",
#              "viridis",
#              "ebirdst",
#              "terra",
#              "colorspace",
#              "gdata",
#              "pROC",
#              "exactextractr",
#              "napops"
#   )
# 
# if (any(!my_packs %in% installed.packages()[, 'Package'])) {install.packages(my_packs[which(!my_packs %in% installed.packages()[, 'Package'])],dependencies = TRUE)}
# lapply(my_packs, require, character.only = TRUE)

library(tidyverse)
library(magrittr)
library(sf)
library(terra)
library(exactextractr)
library(ggspatial)
library(colorspace)

library(viridis)

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


# -----------------------------------------------
# Useful functions
# -----------------------------------------------

`%!in%` <- Negate(`%in%`)

AEA_proj <- "+proj=aea +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-100 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs "

BCR_CAN <- st_read("../../../Data/Spatial/National/BCR/BCR_Terrestrial_master.shp")  %>%
  subset(COUNTRY %in% "CANADA") %>%
  st_make_valid() %>%
  dplyr::select(BCR, PROVINCE_S) %>%
  st_transform(crs = AEA_proj)

BCR_NA <- st_read("../../../Data/Spatial/National/BCR/BCR_Terrestrial_master.shp")  %>%
  subset(COUNTRY %in% c("CANADA","USA")) %>%
  st_make_valid() %>%
  dplyr::select(BCR, PROVINCE_S) %>%
  st_transform(crs = AEA_proj)

# ******************************************************************
# LOAD BBMP DATA PACKAGE
# ******************************************************************

# -----------------------------------------------
# Data from WildTrax (prepared by previous script)
# -----------------------------------------------

WT_dat <- readRDS(file = "../Data_Cleaned/WildTrax/WT_dat.rds")

BBMP_boundary <- WT_dat$BBMP_boundary %>%
  st_transform(crs = AEA_proj) %>%
  st_intersection(BCR_CAN %>% st_union())

BBMP_surveys <- WT_dat$WT_surveyinfo %>% st_transform(crs = AEA_proj) %>% mutate(Obs_Index = 1:nrow(.))
BBMP_counts <- WT_dat$WT_matrix

# Survey locations
in_BBMP <- BBMP_surveys %>%
  st_transform(crs = AEA_proj) %>%
  st_intersects(BBMP_boundary) %>%
  as.numeric()

BBMP_surveys <- BBMP_surveys[!is.na(in_BBMP),]
BBMP_counts <- BBMP_counts[BBMP_surveys$Obs_Index,]

# ******************************************************************
# LOAD BAM DATA PACKAGE, CROP TO BBMP EXTENT
# ******************************************************************

bam <- load("../../../Data/Bird_Data_Raw/BAM/04_NM5.0_data_stratify.Rdata")

bam_surveys <- visit %>%
  mutate(Obs_Index = 1:nrow(.)) %>%
  subset(year >= 2013) %>%
  st_as_sf(coords = c("lon", "lat"),crs=4326, remove = FALSE) %>%
  st_transform(AEA_proj)
within_BBMP <- st_intersects(bam_surveys,BBMP_boundary) %>% as.numeric()
bam_surveys <- bam_surveys[!is.na(within_BBMP),]

bam_counts <- bird[bam_surveys$Obs_Index,]
bam_offsets <- offsets[bam_surveys$Obs_Index,]

# ******************************************************************
# REMOVE BAM SURVEYS WITHIN 100 metres OF BBMP SURVEYS
# ******************************************************************

# Identify surveys in the BAM dataset that are in the BBMP dataset
# Criteria: within 100 m
bam_surveys$Obs_Index <- 1:nrow(bam_surveys)
BBMP_buff <- BBMP_surveys %>% st_buffer(100) %>% st_union()
within_BBMP <- st_intersects(bam_surveys,BBMP_buff) %>% as.numeric()
bam_surveys <- bam_surveys[is.na(within_BBMP),]
bam_counts <- bam_counts[bam_surveys$Obs_Index,]

# ******************************************************************
# PLOT DATA AVAILABILITY ACROSS THE BBMP SURVEY REGION
# PLOT CWS data in a separate color from other surveys
# ******************************************************************

cws_org <- c("CWS-NOR","CWS-ONT","CWS-PRA","ECCC")
xlim <- range(as.data.frame(st_coordinates(BBMP_boundary))$X)
ylim <- range(as.data.frame(st_coordinates(BBMP_boundary))$Y)

survey_map <- ggplot() +
  geom_sf(data = BCR_NA, col = "white", fill = "white")+
  geom_sf(data = BCR_CAN, col = "gray98", fill = "white")+
  geom_sf(data = BBMP_boundary, col = "transparent", fill = "forestgreen", alpha = 0.1)+
  geom_sf(data = subset(bam_surveys, organization %!in% cws_org), size = 0.1, col = "gray70")+
  geom_sf(data = subset(bam_surveys, organization %in% cws_org), size = 0.1, col = "black")+
  geom_sf(data = BBMP_surveys, size = 0.1, col = "black")+

  annotation_scale(style = "ticks")+

  theme(panel.background = element_rect(fill = desaturate("#bbfdfd", amount = 0.5),
                                        colour = desaturate("#bbfdfd", amount = 0.5),
                                        size = 0.5, linetype = "solid"),
        legend.spacing.y = unit(0, "mm"),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        legend.position=c(0.98,0.02),
        legend.justification = c("right","bottom"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 6))+
  coord_sf(xlim = xlim, ylim = ylim)+
  ggtitle("")

survey_map

png(paste0("../Output/Survey_Maps/BBMP_surveys.png"), width=8, height=8, units="in", res=300, type="cairo")
print(survey_map)
dev.off()

# ******************************************************************
# Process Spatial covariates
# ******************************************************************

# Path to spatial covariates
covar_folder <- "../../../Data/Spatial/"

# Mean Annual Temperature
MAT <- rast(paste0(covar_folder,"National/Bioclimate/Normal_1991_2020_MAT.tif"))

# Mean Annual Precipitation
MAP <- rast(paste0(covar_folder,"National/Bioclimate/Normal_1991_2020_MAP.tif")) 

# Land cover of Canada 2020
lcc2020 <- rast(paste0(covar_folder,"National/LandCoverCanada2020/landcover-2020-classification.tif"))

raster_stack <- list(MAT,MAP) %>% rast()
names(raster_stack) = c("MAT","MAP")

# ---------------------------------------------------
# Only select relevant columns
# ---------------------------------------------------

bam_dat <- bam_surveys %>%dplyr::select(organization,project,lat,lon,date,distance,duration)

bbmp_dat <- BBMP_surveys %>%
  mutate(organization = "CWS",project="BBMP",lat = Latitude, lon = Longitude, date = Date_Time, distance = Inf, duration = Survey_Duration_Minutes) %>%
  dplyr::select(organization,project,lat,lon,date,distance,duration)

alldat <- bind_rows(bam_dat,bbmp_dat)

# ---------------------------------------------------
# Extract covariate values at each survey location
# ---------------------------------------------------

# Extract continuous covariates at each survey
alldat = terra::extract(raster_stack,vect(alldat) , bind = TRUE) %>%
  st_as_sf() %>%
  st_transform(AEA_proj)

alldat_1km <- alldat %>% st_buffer(1000)

# Proportion of each land cover class from lcc 2020
prop_LCC_1km <- exact_extract(lcc2020,alldat_1km,"frac") %>% suppressWarnings()
names(prop_LCC_1km) <- paste0(str_replace(names(prop_LCC_1km),"frac","LCC"),"_1km")
prop_LCC_1km[setdiff(paste0("LCC_",seq(1,18),"_1km"),names(prop_LCC_1km))] <- 0
prop_LCC_1km %<>% dplyr::select(sort(names(.)))

# Combine land cover class names
prop_LCC_1km <- prop_LCC_1km %>%
  mutate(
    Needleleaf_forest_1km = LCC_1_1km + LCC_2_1km,
    Mixed_forest_1km = LCC_5_1km + LCC_6_1km,
    Grass_shrub_1km = LCC_8_1km + LCC_10_1km,
    Crop_1km = LCC_15_1km,
    Urban_1km = LCC_17_1km,
    Wetland_1km = LCC_14_1km,
    Water_1km = LCC_18_1km) %>%
  dplyr::select(Needleleaf_forest_1km:Water_1km)

alldat <- bind_cols(alldat,prop_LCC_1km)

alldat$Obs_Index <- 1:nrow(alldat)

allcounts <- bind_rows(as.data.frame(bam_counts),as.data.frame(BBMP_counts))

alldat <- na.omit(alldat)
allcounts <- allcounts[alldat$Obs_Index,]

# ---------------------------------------------------
# Conduct PCA
# ---------------------------------------------------
covars_for_PCA <- c("MAT",
                    "MAP",
                    "Needleleaf_forest_1km",
                    "Mixed_forest_1km",
                    "Grass_shrub_1km",
                    "Crop_1km",
                    "Urban_1km",
                    "Wetland_1km",
                    'Water_1km')

dat_for_PCA <- alldat %>%
  as.data.frame() %>%
  dplyr::select(covars_for_PCA)

pca <- prcomp(dat_for_PCA, scale = TRUE)

alldat_PCA <- predict(pca, newdata = as.data.frame(alldat)[,names(pca$center)])

# Scale/center covariates
for (covar in colnames(alldat_PCA)){

  covar_mean <- mean(alldat_PCA[,covar],na.rm = TRUE)
  covar_sd <- sd(alldat_PCA[,covar],na.rm = TRUE)

  alldat_PCA[,covar] <- (as.data.frame(alldat_PCA)[,covar] - covar_mean)/covar_sd

}

alldat <- alldat %>% bind_cols(alldat_PCA)

# Skree plot
factoextra::fviz_eig(pca)

# Biplot
factoextra::fviz_pca_var(pca,
             axes = c(1,2),
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = viridis(10),
             repel = TRUE     # Avoid text overlapping
)

# ---------------------------------------------------
# Create cross-validation folds
# ---------------------------------------------------

# Create 100 x 100 km spatial folds
set.seed(999)
n_folds <- 5
Crossval_Grid <- st_make_grid(
  st_buffer(BBMP_boundary,50000),
  cellsize = units::set_units(50*50,km^2),
  what = "polygons",
  square = TRUE,
  flat_topped = FALSE) %>%
  st_as_sf() %>%
  na.omit() %>%
  mutate(id = sample(1:nrow(.),nrow(.))) %>%
  mutate(Crossval_Fold = cut(id,breaks = seq(0,max(id)+1,length.out = n_folds+1)) %>% as.factor() %>% as.numeric())%>%
  dplyr::select(Crossval_Fold,x)

alldat$Obs_Index <- 1:nrow(alldat)
alldat <- alldat %>% st_intersection(Crossval_Grid) %>% arrange(Obs_Index)

# ---------------------------------------------------
# Save relevant data so the stuff above doesn't have to be re-run
# ---------------------------------------------------

analysis_data_package <- list(alldat = alldat,
                              allcounts = allcounts,
                              BBMP_counts = BBMP_counts,
                              BBMP_boundary = BBMP_boundary,
                              Crossval_Grid = Crossval_Grid,
                              AEA_proj = AEA_proj,
                              pca = pca)

saveRDS(analysis_data_package,file = "../Data_Cleaned/analysis_data_package.RDS")

