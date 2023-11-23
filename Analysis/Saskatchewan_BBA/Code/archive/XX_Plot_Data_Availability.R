# ------------------------------------------------
# Load/install packages and set graphical themes / working directory
# ------------------------------------------------

my_packs = c('tidyverse',
             'sf',
             'terra',
             'ggspatial',
             'naturecounts',
             'wildRtrax',
             'viridis',
             'cowplot')

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

# -----------------------------------------------
# WildTrax credentials
# -----------------------------------------------

Sys.setenv(WT_USERNAME = 'diles', WT_PASSWORD = 'WildTrax15!')
wt_auth() # Need your wildtrax username

# -----------------------------------------------
# Useful functions
# -----------------------------------------------

`%!in%` <- Negate(`%in%`)

# -----------------------------------------------
# Load study area boundary
# -----------------------------------------------

AEA_proj <- "+proj=aea +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-106 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs "

SaskBoundary <- st_read("../../../Data/Spatial/National/BCR/BCR_Terrestrial_master.shp")  %>%
  subset(PROVINCE_S == "SASKATCHEWAN") %>%
  st_make_valid() %>%
  dplyr::select(BCR, PROVINCE_S) %>%
  st_union() %>%
  st_transform(st_crs(AEA_proj))

SaskBCR <- st_read("../../../Data/Spatial/National/BCR/BCR_Terrestrial_master.shp")  %>%
  subset(PROVINCE_S == "SASKATCHEWAN") %>%
  st_make_valid() %>%
  dplyr::select(BCR, PROVINCE_S) %>%
  st_transform(st_crs(AEA_proj))

# -----------------------------------------------
# Load Boreal Bird Monitoring Program sampling frame
# -----------------------------------------------

BBMP_boundary <- st_read("../../../Data/Spatial/National/BBMP_hexagons/Hexagons_wo_Attributes.shp") %>%
  st_union() %>%
  st_transform(st_crs(AEA_proj)) %>%
  st_intersection(SaskBoundary)

# -----------------------------------------------
# Create hexagon layer to summarize data availability across the study area
# -----------------------------------------------

hexgrid <- st_make_grid(
  SaskBoundary,
  cellsize = units::set_units(10*10,km^2),
  what = "polygons",
  square = FALSE,
  flat_topped = FALSE) %>%
  st_intersection(SaskBoundary) %>%
  st_make_valid() %>%
  st_as_sf()
hexgrid$hexID = 1:nrow(hexgrid)

hexgrid_centroids <- st_centroid(hexgrid)

# Hexagons in boreal
BBMP_hexagons <- st_intersection(hexgrid_centroids,BBMP_boundary)
hexgrid$Boreal <- "No"
hexgrid$Boreal[which(hexgrid$hexID %in% BBMP_hexagons$hexID)] <- "Yes"

# -----------------------------------------------
# Load BAM data package
# -----------------------------------------------

load("../../../Data/Bird_Data_Raw/BAM/BAM_data_package_August2019.RData")

BAM <- full_join(PKEYcombo,SScombo) %>%
  na.omit() %>%
  unique() %>%
  st_as_sf(coords = c("X","Y"), crs = LCC) %>%
  st_transform(crs = st_crs(AEA_proj))

rm(offcombo,LCC,PKEYcombo,PCcombo,SScombo,TAX)

# Remove BAM data more than 15 years old
n_BAM_old <- sum(BAM$YEAR < (2023-15))
BAM <- BAM %>% subset(YEAR >= (2023-15))

BAM$id <- paste0(BAM$PKEY,BAM$YEAR,BAM$ARU)

# Trim BAM dataset to boundary of sampling frame
BAM <- BAM %>% st_intersection(SaskBoundary)

# -----------------------------------------------
# Summarize BAM data within each square
# -----------------------------------------------

BAM_hexagons <- st_intersection(BAM, hexgrid) %>%
  as.data.frame() %>%
  group_by(hexID) %>%
  summarize(n_BAM = n())

hexgrid <- hexgrid %>% left_join(BAM_hexagons)

# -----------------------------------------------
# Summarize atlas data within each square
# -----------------------------------------------

SK_dat <- readRDS(file = "../Data_Cleaned/analysis_data.rds")
Checklist_data <- readRDS("../Data_Cleaned/NatureCounts/Checklist_data.rds")

atlas_PC_ARU <- SK_dat$PC_ARU_surveyinfo %>%
  mutate(Year = lubridate::year(Date_Time)) %>%
  subset(Year >= (2023-15)) %>%
  bind_rows(Checklist_data$DO_surveyinfo %>% st_transform(st_crs(SK_dat$PC_ARU_surveyinfo)))

atlas_PC_ARU <- st_intersection(atlas_PC_ARU, hexgrid) %>%
  as.data.frame() %>%
  group_by(hexID) %>%
  summarize(n_PC_ARU = n())

hexgrid <- hexgrid %>% left_join(atlas_PC_ARU)

# -----------------------------------------------
# Plot
# -----------------------------------------------

Map1 <- ggplot()+
  
  geom_sf(data = SaskBoundary, fill = "gray95", col = "black") +
  geom_sf(data = BBMP_boundary, col = "transparent", fill = "darkgreen", alpha = 0.1) +
  
  geom_sf(data = subset(hexgrid, !is.na(n_BAM)), col = "transparent", fill = "gray50")+
geom_sf(data = SaskBoundary, fill = "transparent", col = "black") +
  ggtitle("Prior data (within 15 years)")

Map2 <- ggplot()+
  
  geom_sf(data = SaskBoundary, fill = "gray95", col = "black") +
  geom_sf(data = BBMP_boundary, col = "transparent", fill = "darkgreen", alpha = 0.1) +
  
  geom_sf(data = subset(hexgrid, !is.na(n_PC_ARU) & Boreal == "Yes"), col = "transparent", fill = "darkgreen")+
  geom_sf(data = subset(hexgrid, !is.na(n_PC_ARU) & Boreal == "No"), col = "transparent", fill = "gray50")+
  
  geom_sf(data = SaskBoundary, fill = "transparent", col = "black") +
  ggtitle("New data (Atlas + BBMP)")

# Combine maps
plot_grid(Map1,Map2,nrow=1,align="hv")

