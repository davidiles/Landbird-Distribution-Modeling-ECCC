# ---------------------------------------------------------------------------
# ***************************************************************************
# Download WildTrax datasets within a study area boundary
# ***************************************************************************
# ---------------------------------------------------------------------------

# ------------------------------------------------
# Load/install packages and set graphical themes / working directory
# ------------------------------------------------

# remotes::install_github("ABbiodiversity/wildRtrax@development")


my_packs = c('tidyverse',
             'openxlsx',
             'RColorBrewer',
             'viridis',
             'ggrepel',
             'scales',
             'wildrtrax',
             'lubridate',
             'sf'
)

if (any(!my_packs %in% installed.packages()[, 'Package'])) {install.packages(my_packs[which(!my_packs %in% installed.packages()[, 'Package'])],dependencies = TRUE)}
lapply(my_packs, require, character.only = TRUE)

rm(list=ls())

# -----------------------------------------------
# Useful functions
# -----------------------------------------------

`%!in%` <- Negate(`%in%`)

download_one_project <- function(project_id, sensor_type = "ARU") {
  start <- Sys.time()
  
  out <- tryCatch({
    
    aru_dat_list <- wildrtrax::wt_download_report(
      project_id = project_id,
      sensor_id  = sensor_type,
      report     = c("location", "recording"),
      max_seconds = 180
    )
    
    loc <- aru_dat_list[[1]] %>%
      dplyr::select(location_id, latitude, longitude)
    
    rec <- aru_dat_list[[2]] %>%
      dplyr::select(location_id, recording_id, task_method, task_duration)
    
    dat <- dplyr::full_join(loc, rec, by = "location_id") %>%
      dplyr::mutate(project_id = project_id)
    
    elapsed <- Sys.time() - start
    
    message(
      sprintf(
        "SUCCESS: project_id = %s | time = %.2f sec",
        project_id,
        as.numeric(elapsed, units = "secs")
      )
    )
    
    list(
      success = TRUE,
      project_id = project_id,
      data = dat,
      error = NA_character_,
      elapsed_sec = as.numeric(elapsed, units = "secs")
    )
    
  }, error = function(e) {
    elapsed <- Sys.time() - start
    
    message(
      sprintf(
        "FAILED: project_id = %s | time = %.2f sec | error = %s",
        project_id,
        as.numeric(elapsed, units = "secs"),
        conditionMessage(e)
      )
    )
    
    list(
      success = FALSE,
      project_id = project_id,
      data = NULL,
      error = conditionMessage(e),
      elapsed_sec = as.numeric(elapsed, units = "secs")
    )
  })
  
  out
}

# -----------------------------------------------
# Load list of species prepared by script 0_species_list.R
# -----------------------------------------------

all_species <- readRDS("Data_Cleaned/all_species.RDS")

# -----------------------------------------------
# Polygon delineating study area boundary
# -----------------------------------------------

# BCR <- st_read("../../Data/Spatial/National/BCR/BCR_Terrestrial_master.shp")  %>%
#   subset(COUNTRY == "CANADA") %>%
#   st_make_valid() %>%
#   dplyr::select(BCR, PROVINCE_S)
# 
# BBMP_boundary <- st_read("../../Data/Spatial/National/BBMP_hexagons/Hexagons_w_Attributes.shp") %>%
#   st_union() %>%
#   st_buffer(1000) %>%
#   st_buffer(-1000)
# 
# ggplot()+
#   geom_sf(data = BCR)+
#   geom_sf(data = BBMP_boundary, fill = "darkgreen", alpha = 0.5)

# -----------------------------------------------
# WildTrax credentials
# -----------------------------------------------

wt_auth() # Need your wildtrax username

# ---------------------------------------------------------
# Download projects with 'ARU' data
# ---------------------------------------------------------
aru_path <- "Data_Cleaned/WildTrax/wt_aru.csv"
projs <- wildrtrax::wt_get_projects(sensor = "ARU")

aru_fulldat <- data.frame()
if (file.exists(aru_path)) aru_fulldat <- read.csv(aru_path)
for (i in 1:nrow(projs)){
  print(i)
  if (projs$project_id[i] %in% aru_fulldat$project_id) next
  if (file.exists(aru_path)) aru_fulldat <- read.csv(aru_path)
  dat <- download_one_project(project_id = projs$project_id[i], sensor_type = projs$project_sensor[i])
  aru_fulldat <- bind_rows(aru_fulldat, dat$data)
  write.csv(aru_fulldat, file = "Data_Cleaned/WildTrax/wt_aru.csv", row.names = FALSE)
  
}

