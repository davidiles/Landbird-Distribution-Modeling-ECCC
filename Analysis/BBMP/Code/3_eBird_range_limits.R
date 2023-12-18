# ---------------------------------------------------------------------------
# ***************************************************************************
# Download NatureCounts datasets within a study area boundary
#  - merge with WildTrax data, and remove duplicated records
# ***************************************************************************
# ---------------------------------------------------------------------------

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
             "terra")

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
# Useful functions
# -----------------------------------------------

`%!in%` <- Negate(`%in%`)

# -----------------------------------------------
# Data from WildTrax (prepared by previous script)
# -----------------------------------------------
AEA_proj <- "+proj=aea +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-100 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs "

WT_dat <- readRDS(file = "../Data_Cleaned/WildTrax/WT_dat.rds")
WT_surveyinfo <- WT_dat$WT_surveyinfo %>% st_transform(crs = AEA_proj)
WT_matrix <- WT_dat$WT_matrix
BBMP_boundary <- WT_dat$BBMP_boundary %>% st_transform(crs = AEA_proj)

# -----------------------------------------------
# Load list of species prepared by script 0_species_list.R
# -----------------------------------------------

all_species <- readRDS("../Data_Cleaned/all_species.RDS")
BSC_species <- readRDS("../Data_Cleaned/BSC_species.RDS")

# ------------------------------------------------
# Map of survey locations
# ------------------------------------------------

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

#Trim BBMP boundary to land border
BBMP_boundary <- BBMP_boundary %>% st_transform(crs = AEA_proj) %>%
  st_intersection(BCR_CAN %>% st_union())

# Survey locations
in_BBMP <- WT_dat$WT_surveyinfo %>%
  st_transform(crs = AEA_proj) %>%
  st_intersects(BBMP_boundary) %>%
  as.numeric()

survey_dat <- WT_dat$WT_surveyinfo %>%
  st_transform(crs = AEA_proj) %>%
  mutate(Obs_Index = 1:nrow(.))
survey_dat <- survey_dat[!is.na(in_BBMP),]
WT_matrix <- WT_matrix[survey_dat$Obs_Index,]

# ******************************************************************
# Evaluate detections outside the ebird range boundaries
# ******************************************************************

# Need to set access key (only once) before downloading ranges
# ebirdst::set_ebirdst_access_key()

# ------------------------------------------------
# Ensure path is correctly set
# ------------------------------------------------

usethis::edit_r_environ()

# This should read:
# EBIRDST_KEY='ntqm1ha68fov'
# EBIRDST_DATA_DIR='C:/Users/IlesD/OneDrive - EC-EC/Iles/Projects/Landbirds/Landbird-Distribution-Modeling-ECCC/Data/Spatial/eBird/'
# 

# # Download species ranges
# for (sp_code in colnames(WT_matrix)){
#   
#   species_name = BSC_species$english_name[which(BSC_species$BSC_spcd == sp_code)]
#   
#   # Check if species is available
#   check <- get_species(species_name)
#   if (length(check)==0 | is.na(check)){
#     print(paste0(sp_code, " not available in eBird"))
#     next
#   }
#   
#   ebirdst_download_status(species = species_name, download_abundance = FALSE,download_ranges = TRUE, force = FALSE)
# }

survey_dat$Obs_Index <- 1:nrow(survey_dat)
species_results <- data.frame()

for (sp_code in colnames(WT_matrix)){
  
  print(sp_code)
  
  # ------------------------------------------------
  # Download and trim ebird range for this species
  # ------------------------------------------------
  
  species_name = BSC_species$english_name[which(BSC_species$BSC_spcd == sp_code)]
  
  # Check if species is available
  check <- get_species(species_name)
  if (length(check)==0 | is.na(check)){
    print(paste0(sp_code, " not available in eBird"))
    next
  }
  
  range <- load_ranges(species = species_name, resolution = "27km") %>%
    subset(season %in% c("resident","breeding"))
  
  date_range <- yday(ymd(c(range$start_date, range$end_date))) %>% range()
  
  # Extract data within the appropriate (non-migratory) date range
  sp_dat <- survey_dat %>% mutate(count = WT_matrix[,sp_code]) %>%
    subset(count > 0)%>% 
    subset(
      yday(Date_Time) >= date_range[1] &
        yday(Date_Time) <= date_range[2])
  
  # Extract detections within and outside of date range
  range <- range %>%  st_transform(.,crs = st_crs(BBMP_boundary)) %>%
    st_union() %>%
    st_crop(BBMP_boundary)
  
  in_range <- st_intersects(sp_dat,range) %>% as.numeric()
  in_range[is.na(in_range)] <- 0
  sp_dat$in_range <- in_range
  
  # Save number of detections out of range
  n_det <- nrow(sp_dat)
  n_out <- sum(sp_dat$in_range == 0)
  species_results <- rbind(species_results, data.frame(sp_code = sp_code,
                                                       n_det = n_det,
                                                       n_out = n_out))
  # xlim <- range(as.data.frame(st_coordinates(BBMP_boundary))$X)
  # ylim <- range(as.data.frame(st_coordinates(BBMP_boundary))$Y)
  # 
  # # Plot detections
  #   sp_dat$in_range <- factor(sp_dat$in_range,levels = c(0,1), labels = c("No","Yes"))
  #   sp_map <- ggplot() +
  #     geom_sf(data = BCR_NA, col = "white", fill = "white")+
  #     geom_sf(data = BCR_CAN, col = "gray97", fill = "white")+
  #     geom_sf(data = range, col = "transparent", fill = "#ee2400", alpha = 0.1)+
  #     geom_sf(data = survey_dat, col = "gray80", size = 0.05)+
  #     geom_sf(data = sp_dat, aes(col = in_range, size = in_range))+
  #     scale_color_manual(values = c("black","orangered"), name = "Out of Range", labels = c("Yes","No"),drop = FALSE)+
  #     scale_size_manual(values = c(0.5,0.1), name = "Out of Range", labels = c("Yes","No"),drop = FALSE)+
  #     annotation_scale(style = "ticks")+
  #     
  #     theme(panel.background = element_rect(fill = "#e3fefe",
  #                                           colour = "#e3fefe",
  #                                           size = 0.5, linetype = "solid"),
  #           legend.spacing.y = unit(0, "mm"),
  #           panel.border = element_rect(colour = "black", fill=NA),
  #           #aspect.ratio = 0.2, axis.text = element_text(colour = 1, size = 12),
  #           legend.background = element_blank(),
  #           legend.box.background = element_rect(colour = "black"),
  #           legend.position=c(0.98,0.98),
  #           legend.justification = c("right","top"),
  #           legend.title = element_text(size = 10),
  #           legend.text = element_text(size = 6))+
  #     coord_sf(xlim = xlim, ylim = ylim)+
  #     ggtitle(sp_code)
  #   
  #   png(paste0("../Output/Range_Maps/",sp_code,"_range.png"), width=8, height=6.5, units="in", res=300, type="cairo")
  #   print(sp_map)
  #   dev.off()
} # close species loop



species_results$prop_out <- species_results$n_out/species_results$n_det
species_results <- species_results %>%
  arrange(n_out)
species_results$sp_code <- factor(species_results$sp_code, levels = species_results$sp_code)

barplot <- ggplot(data = subset(species_results, n_out >= 50))+
  geom_bar(aes(x = n_out, y = sp_code, fill = n_out),stat = "identity")+
  scale_fill_gradientn(colors = rev(viridis(10)), guide = "none")+
  theme_bw()+
  xlab("Number of detections")+
  ylab("Species")+
  ggtitle("Number of detections outside eBird breeding range")
print(barplot)


# 
barplot <- ggplot(data = subset(species_results, n_out >= 50))+
  geom_bar(aes(x = n_out, y = sp_code, fill = prop_out),stat = "identity")+
  scale_fill_gradientn(colors = rev(magma(10)), limits = c(0,1), name = "Proportion")+
  theme_bw()+
  xlab("Number of detections")+
  ylab("Species")+
  ggtitle("Number of detections outside eBird breeding range\n\nShading = proportion of total detetions outside range")
print(barplot)


