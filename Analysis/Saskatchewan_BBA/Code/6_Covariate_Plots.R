# ************************************************
# BAYESIAN ANALYSIS / SPECIES DISTRIBUTION MODELS FOR SASKATCHEWAN BREEDING BIRD ATLAS
# 
# 1) A 'processed data package' for analysis is prepared by previous script
# ************************************************

rm(list=ls())

# ------------------------------------------------
# Set working directory
# ------------------------------------------------
# 
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
dirname <- "C:/Users/IlesD/OneDrive - EC-EC/Iles/Backup/2023-11-09/SDM_ECCC/Analysis/Saskatchewan/Code"
setwd(dirname)

`%!in%` <- Negate(`%in%`)

# ------------------------------------------------
# Load packages
# ------------------------------------------------

require(tidyverse)
require(sf)
require(viridis)
require(ggpubr)
require(stars)

# ------------------------------------------------
# LOAD DATA
# ------------------------------------------------

analysis_data <- readRDS("../Data_Cleaned/analysis_data_package.rds")
attach(analysis_data)

AEA_proj <- "+proj=aea +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-106 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs "
SaskBoundary <- st_read("../../../Data/Spatial/National/BCR/BCR_Terrestrial_master.shp")  %>%
  subset(PROVINCE_S == "SASKATCHEWAN") %>%
  st_make_valid() %>%
  st_union() %>%
  st_transform(st_crs(AEA_proj))

SaskRast <- SaskGrid %>% 
  dplyr::select(elevation_1km:PC10) %>%
  stars::st_rasterize()

covar_to_plot <- names(SaskRast)

covar_plotlist <- list()

for (covar in covar_to_plot){
  cplot <- ggplot() + geom_stars(data = SaskRast, aes(fill = !!sym(covar)))+
    scale_fill_gradientn(colours = viridis(10), name = covar)+
    geom_sf(data = SaskBoundary,colour="black",fill=NA,lwd=0.3,show.legend = F)+
    ggtitle(covar)+
    theme_bw()+
    xlab("")+ylab("")
  
  covar_plotlist[[covar]] <- cplot
}

covar_plots <- ggarrange(plotlist = covar_plotlist,nrow=2,ncol=10)

png("../Output/Covariate_Maps/Covariate_Maps.png", width=50, height=10, units="in", res=300, type="cairo")
print(covar_plots)
dev.off()