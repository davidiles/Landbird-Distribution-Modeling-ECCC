# ============================================================================
# DOWNLOAD / FORMAT RASTERS; EXTRACTION DONE IN NEXT SCRIPT
# ============================================================================

# ---- Load/install packages
my_packs <- c('rgee',
              'tidyverse', 
              'sf',
              'terra')

if (any(!my_packs %in% installed.packages()[, 'Package'])){
  install.packages(my_packs[which(!my_packs %in% installed.packages()[, 'Package'])], dependencies = TRUE)
}
lapply(my_packs, require, character.only = TRUE)

rm(list = ls())

`%!in%` <- Negate(`%in%`)

# Initialize Earth Engine
ee_Initialize()

# ============================================================================
# Load data package prepared by script 1
# ============================================================================

analysis_data <- readRDS(file = "data_clean/analysis_data.rds")
all_surveys <- analysis_data$NC_surveyinfo
full_count_matrix <- analysis_data$full_count_matrix
all_species <- analysis_data$all_species
rm(analysis_data)

AEA_proj <- st_crs(all_surveys)

# ============================================================================
# Prepare polygon for data extraction (Ontario boundary + 100 km)
# ============================================================================

ON_boundary <- st_read("data/spatial/BCR/BCR_Terrestrial_master.shp") %>%
  subset(PROVINCE_S == "ONTARIO") %>%
  st_make_valid() %>%
  dplyr::select(BCR, PROVINCE_S) %>%
  st_transform(AEA_proj)%>%
  st_union() %>%
  st_set_precision(1)

# Add 25 km buffer so that rasters extend outside study area
ON_boundary_25km <- st_buffer(ON_boundary, 25000)

# Reproject boundary to crs 4326
boundary_ee <- ON_boundary_25km %>%
  st_transform(4326) %>%
  sf_as_ee()

# ****************************************************************************
# ****************************************************************************
# STEP 1: DOWNLOAD RASTERS FROM GOOGLE EARTH ENGINE
# ****************************************************************************
# ****************************************************************************

# ============================================================================
# Download MODIS MCD12Q1 - https://developers.google.com/earth-engine/datasets/catalog/MODIS_061_MCD12Q1
# ============================================================================

# Load MODIS MCD12Q1 v061 land cover
modis_lc <- ee$ImageCollection("MODIS/061/MCD12Q1")

# Select year of interest (e.g., 2001)
for (year in c(2001:2005, 2020:2024)){
  filename <- paste0("data/spatial/MODIS/modis_lc_", year, ".tif")
  if (!file.exists(filename)){
    
    # Get MODIS image
    modis_img <- modis_lc$
      filter(ee$Filter$calendarRange(year, year, "year"))$
      first()$
      select("LC_Type1")$
      clip(boundary_ee) %>%
      
      # Reproject in Earth Engine before download
      ee$Image$reproject(
        crs = "EPSG:3978",
        scale = 500
      )
    
    # Bring raster directly into R (stars object)
    modis_r <- ee_as_rast(
      image = modis_img,
      region = boundary_ee, 
      scale = 500,
      crs = "EPSG:3978",       
      via = "drive"          
    )
    
    # Save locally
    terra::writeRaster(modis_r, paste0("data/spatial/MODIS/modis_lc_", year, ".tif"), overwrite = TRUE)
  }
  
}

# # ============================================================================
# # Download Global Forest Change - https://developers.google.com/earth-engine/datasets/catalog/UMD_hansen_global_forest_change_2024_v1_12
# # - note: not implemented for this analysis
# # ============================================================================
# 
# # Province must be split into tiles to download
# tile_size_km <- 100  # size of tiles in km
# export_folder <- "rgee_hansen_tiles" # Google Drive folder for exports
# 
# # Convert tile size to degrees roughly (1 deg ~ 111 km)
# tile_size_deg <- tile_size_km / 111
# 
# # Split study area into tiles
# ON_boundary_25km <- ON_boundary_25km %>% st_transform(4326) %>% st_make_valid()
# 
# tiles <- st_make_grid(ON_boundary_25km,
#                       cellsize = tile_size_deg,
#                       what = "polygons") %>%
#   st_make_valid() %>%
#   st_as_sf()
# 
# tiles <- st_intersection(tiles, ON_boundary_25km)  # keep only intersection with boundary
# 
# # Convert tiles to Earth Engine objects
# tiles_ee <- lapply(seq_len(nrow(tiles)), function(i) sf_as_ee(tiles[i,]))
# 
# # Load Hansen dataset
# hansen <- ee$Image("UMD/hansen/global_forest_change_2024_v1_12")
# 
# # ---- Loop over tiles and bands
# # Select bands to download
# bands <- c("treecover2000", "loss", "gain", "lossyear")
# 
# for(i in seq_along(tiles_ee)){
#   
#   cat("Processing tile", i, "of", length(tiles_ee), "\n")
#   tile_ee <- tiles_ee[[i]]
#   tile_geom <- tile_ee$geometry()
#   
#   for(b in bands){
#     
#     # Clip band to tile
#     img <- hansen$select(b)$clip(tile_ee)
#     
#     # Local filename
#     filename <- paste0("data/spatial/Hansen_ForestChange/tiled_download/", b, "_tile_", i, ".tif")
#     
#     if(!file.exists(filename)){
#       
#       # Export raster via Drive
#       tile_r <- ee_as_rast(
#         image = img,
#         region = tile_geom,
#         scale = 30,
#         via = "drive",
#         maxPixels = 1e13
#       )
#       
#       # Save locally
#       terra::writeRaster(tile_r, filename, overwrite = TRUE)
#     }
#   }
# }
# 
# # Optional: Merge tiles per band
# for(b in bands){
#   files <- list.files("data/spatial/Hansen", pattern = paste0("^", b, "_tile_.*\\.tif$"), full.names = TRUE)
#   r_list <- lapply(files, terra::rast)
#   r_merged <- do.call(terra::merge, r_list)
#   
#   # Save merged raster
#   terra::writeRaster(r_merged, paste0("data/spatial/Hansen/", b, "_Ontario.tif"), overwrite = TRUE)
# }

# ============================================================================
# Download GHSL (global human settlement layer) Degree of Urbanization (SMOD-POP) - https://developers.google.com/earth-engine/datasets/catalog/JRC_GHSL_P2023A_GHS_SMOD_V2-0?hl=en
# ============================================================================

# Load product
ghsl_collection <- ee$ImageCollection("JRC/GHSL/P2023A/GHS_SMOD_V2-0")

# Check what years are available in the collection
available_years <- c()
collection_size <- ghsl_collection$size()$getInfo()
image_list <- ghsl_collection$toList(collection_size)

for(i in 0:(collection_size-1)) {
  image <- ee$Image(image_list$get(i))
  time_start <- image$get('system:time_start')
  date <- ee$Date(time_start)
  year <- date$get('year')$getInfo()
  available_years <- c(available_years, year)
}

print("Available years in GHSL collection:")
print(sort(available_years))

# Years of interest
for (year in c(2000, 2020)) {
  filename <- paste0("data/spatial/GHSL/ghsl_urban_", year, ".tif")
  
  if (!file.exists(filename)) {
    ghsl_img <- ghsl_collection$
      filterDate(paste0(year, '-01-01'), paste0(year, '-12-31'))$
      first()$
      select("smod_code")$        # "smod" band = Settlement classification
      clip(boundary_ee)  %>%
      
      # Reproject in Earth Engine before download
      ee$Image$reproject(
        crs = "EPSG:3978",
        scale = 500
      )      
    
    # Export raster via Google Drive, then pull into R
    ghsl_r <- ee_as_rast(
      image = ghsl_img,
      region = boundary_ee,    # FeatureCollection/geometry of Ontario
      scale = 1000,            # dataset native resolution = 1 km
      crs = "EPSG:3978",  
      via = "drive"
    )
    
    # Save locally
    terra::writeRaster(
      ghsl_r,
      filename,
      overwrite = TRUE
    )
    
    print(paste("Saved:", filename))
  }
}


# ****************************************************************************
# ****************************************************************************
# STEP 2: PREPARE COVARIATE LAYERS FOR ANALYSIS
# ****************************************************************************
# ****************************************************************************


# ============================================================================
# ---- Land cover (raster) from https://developers.google.com/earth-engine/datasets/catalog/MODIS_061_MCD12Q1
# ============================================================================

# Create a new synthetic layer for each atlas period that contains:
# 1) Proportion of years that each pixel was in each class
# 2) Modal composite: the most common land cover class in each pixel

rasters_OBBA2 <- rast(c("data/spatial/MODIS/modis_lc_2001.tif", 
                        "data/spatial/MODIS/modis_lc_2002.tif", 
                        "data/spatial/MODIS/modis_lc_2003.tif", 
                        "data/spatial/MODIS/modis_lc_2004.tif", 
                        "data/spatial/MODIS/modis_lc_2005.tif"))

rasters_OBBA3 <- rast(c("data/spatial/MODIS/modis_lc_2020.tif", 
                        "data/spatial/MODIS/modis_lc_2021.tif", 
                        "data/spatial/MODIS/modis_lc_2022.tif", 
                        "data/spatial/MODIS/modis_lc_2023.tif", 
                        "data/spatial/MODIS/modis_lc_2024.tif"))

# Function to compute the modal (most frequent) value across layers
modal_fun <- function(x, ...) {
  ux <- na.omit(x)
  if (length(ux) == 0) return(NA)
  tab <- table(ux)
  as.numeric(names(tab)[which.max(tab)])
}

## Function that, for each pixel (vector across years),
## returns the fraction of years classified in each class.
fractional_cover <- function(rasters, classes) {
  n_years <- nlyr(rasters)
  out_list <- list()
  for (k in classes) {
    # Each pixel: proportion of years == k
    frac <- sum(rasters == k, na.rm=TRUE) / n_years
    names(frac) <- paste0("frac_", k)
    out_list[[as.character(k)]] <- frac
  }
  rast(out_list)
}

#Use IGBP classes
classes <- c(1,4,5,8,9,11,12,13,14,17)

# Fraction of years that each pixel was classified as each class during each atlas period
frac_OBBA2 <- fractional_cover(rasters_OBBA2, classes)
frac_OBBA3 <- fractional_cover(rasters_OBBA3, classes)

## Compute modal class composites for each period
mode_OBBA2 <- modal(rasters_OBBA2, ties = "first", na.rm = TRUE)
mode_OBBA3 <- modal(rasters_OBBA3, ties = "first", na.rm = TRUE)

# Save outputs
writeRaster(frac_OBBA2, "data_clean/spatial/MCD12Q1_OBBA2_frac.tif", overwrite=TRUE)
writeRaster(frac_OBBA3, "data_clean/spatial/MCD12Q1_OBBA3_frac.tif", overwrite=TRUE)
writeRaster(mode_OBBA2, "data_clean/spatial/MCD12Q1_OBBA2_mode.tif", overwrite=TRUE)
writeRaster(mode_OBBA3, "data_clean/spatial/MCD12Q1_OBBA3_mode.tif", overwrite=TRUE)

table(values(mode_OBBA2))
table(values(mode_OBBA3))

# Notes: no class 2
# very few: deciduous 3: needleleaf larch forests
#                     6: closed shrublands
#                     15: permanent snow and ice 
#                     16: barren 
# those classes will be completely omitted


# ============================================================================
# ---- Climate (WorldClim 2.1)
# ============================================================================

# Function to compute mean climate raster over years and months
mean_climate <- function(variable = "prec",        # e.g., "prec", "tmax"
                         years = 2001:2005,        # vector of years
                         months = 5:7,             # vector of months (5=May,6=Jun,7=Jul)
                         base_dir = "data/spatial/WorldClim", 
                         boundary) {               # sf boundary object
  
  # figure out folder name (assumes years are consecutive in blocks)
  folder <- file.path(base_dir, variable)
  
  year_start <- min(years)
  year_end   <- max(years)
  # get CRS of first raster for reprojection
  test_file <- file.path(folder, 
                         sprintf("wc2.1_cruts4.09_2.5m_%s_%d-%02d.tif", 
                                 variable, years[1], months[1]))
  test_rast <- rast(test_file)
  
  # transform boundary to raster CRS
  boundary_trans <- st_transform(boundary, st_crs(test_rast))
  boundary_vect <- vect(boundary_trans)
  
  # container for rasters
  all_rasters <- list()
  
  for (yr in years) {
    for (mo in months) {
      f <- file.path(folder, 
                     sprintf("wc2.1_cruts4.09_2.5m_%s_%d-%02d.tif", 
                             variable, yr, mo))
      if (!file.exists(f)) stop("Missing file: ", f)
      
      r <- rast(f)
      r_crop <- crop(r, boundary_vect)
      all_rasters <- c(all_rasters, r_crop)
    }
  }
  
  # stack and compute mean
  all_stack <- rast(all_rasters)
  r_mean <- mean(all_stack, na.rm = TRUE)
  
  return(r_mean)
}

# Precipitation (5-year mean up to start of atlas)
mean_prec_OBBA2 <- mean_climate(variable = "prec", 
                                years = 2000:2004, 
                                months = 3:7,
                                base_dir = "data/spatial/WorldClim",
                                boundary = ON_boundary_25km)

mean_prec_OBBA3 <- mean_climate(variable = "prec", 
                                years = 2020:2024, 
                                months = 3:7,
                                base_dir = "data/spatial/WorldClim",
                                boundary = ON_boundary_25km)

zlim <- range(c(values(mean_prec_OBBA2), values(mean_prec_OBBA3)), na.rm = TRUE)
prec_stack <- c(mean_prec_OBBA2, mean_prec_OBBA3)
names(prec_stack) <- c("OBBA2", "OBBA3")
plot(prec_stack, range = zlim)  # range=TRUE forces same scale across layers

# Max Temperature (5-year mean up to start of atlas)
mean_tmax_OBBA2 <- mean_climate(variable = "tmax", 
                                years = 2000:2004, 
                                months = 3:7,
                                base_dir = "data/spatial/WorldClim",
                                boundary = ON_boundary_25km)

mean_tmax_OBBA3 <- mean_climate(variable = "tmax", 
                                years = 2020:2024, 
                                months = 3:7,
                                base_dir = "data/spatial/WorldClim",
                                boundary = ON_boundary_25km)

zlim <- range(c(values(mean_tmax_OBBA2), values(mean_tmax_OBBA3)), na.rm = TRUE)
tmax_stack <- c(mean_tmax_OBBA2, mean_tmax_OBBA3)
names(tmax_stack) <- c("OBBA2", "OBBA3")
plot(tmax_stack, range = zlim)  # range=TRUE forces same scale across layers

# Save rasters as GeoTIFFs
terra::writeRaster(mean_prec_OBBA2,  "data_clean/spatial/prec_OBBA2.tif",  overwrite = TRUE)
terra::writeRaster(mean_prec_OBBA3,  "data_clean/spatial/prec_OBBA3.tif",  overwrite = TRUE)
terra::writeRaster(mean_tmax_OBBA2,  "data_clean/spatial/tmax_OBBA2.tif",  overwrite = TRUE)
terra::writeRaster(mean_tmax_OBBA3,  "data_clean/spatial/tmax_OBBA3.tif",  overwrite = TRUE)

# ============================================================================
# ---- Urbanization (raster from global human settlement layer)
# ============================================================================

# Classes are:
# 0: no data

# Mostly uninhabited / low density
# 10: mostly uninhabited
# 11: rural grid cell (low density)

# 12: suburban / peri-urban (intermediate density) - fewer than 50 people / km2

# 13: urban center / city (high density urban core) - more than 300 people / km2
# 21: low density urban
# 22: intermediate density urban
# 23: high density urban
# 30: urban center / core

GHSL_2000 <- terra::rast("data/spatial/GHSL/ghsl_urban_2000.tif")
GHSL_2020 <- terra::rast("data/spatial/GHSL/ghsl_urban_2020.tif")

# Reclassify function
reclass_ghsl <- function(x) {
  x <- ifel(x %in% c(10,11), 1, x)   # mostly uninhabited
  x <- ifel(x == 12, 2, x)           # sparsely populated
  x <- ifel(x %in% c(13,21,22,23,30), 3, x)  # densely populated
  return(x)
}

GHSL_2000_reclass <- reclass_ghsl(GHSL_2000)
GHSL_2020_reclass <- reclass_ghsl(GHSL_2020)

# Save rasters as GeoTIFFs
terra::writeRaster(GHSL_2000_reclass,  "data_clean/spatial/Urban_2000.tif",  overwrite = TRUE)
terra::writeRaster(GHSL_2020_reclass,  "data_clean/spatial/Urban_2020.tif",  overwrite = TRUE)

# ============================================================================
# ---- Waterbodies
# ============================================================================

water <- st_read("data/spatial/Ontario_Hydro_Network_(OHN)_-_Waterbody/Ontario_Hydro_Network_(OHN)_-_Waterbody.shp")

# Reclassify into:
# 1 - Open water (Lake, Kettle Lake, Pond, Reservoir, Beaver Pond)
# 2 - River (River, Canal)

water <- water %>% 
  mutate(water_class = case_when(
    WATERBODY_ %in% c("Lake", "Kettle Lake", "Pond", "Reservoir", "Beaver Pond") ~ "Open Water",
    WATERBODY_ %in% c("Canal", "River") ~ "River",
    TRUE ~ NA_character_   # everything else set to NA
  )
  )

# Reproject water
water <- st_transform(water, st_crs(AEA_proj))

# Rasterize water polygons
res_m <- 30  # Set raster resolution

# Create empty raster covering Ontario
r <- rast(ext(water), res = res_m, crs = st_crs(water)$proj4string)

# Rasterize water polygons for each class
water_open <- water %>% filter(water_class == "Open Water")
water_river <- water %>% filter(water_class == "River")

# Simplify polygons using tolerance of 10 m = 0.01 km - roughly 1/3 raster resolution
tolerance <- 10
water_open_simple <- st_simplify(water_open, dTolerance = tolerance)
water_river_simple <- st_simplify(water_river, dTolerance = tolerance)

# Rasterize polygons for extracting % open water and % river within defined distance of each feature
r_open  <- rasterize(vect(water_open_simple), r, field = 1, background = 0)
r_river <- rasterize(vect(water_river_simple), r, field = 1, background = 0)

# Save rasters as GeoTIFFs
terra::writeRaster(r_open,  "data_clean/spatial/water_open.tif",  overwrite = TRUE)
terra::writeRaster(r_river, "data_clean/spatial/water_river.tif", overwrite = TRUE)

# ---- Save a simplified water layer that can be used for plotting and removing large waterbodies from predictions
water <- st_read("data/spatial/Ontario_Hydro_Network_(OHN)_-_Waterbody/Ontario_Hydro_Network_(OHN)_-_Waterbody.shp")
water <- water %>% 
  st_make_valid() %>%                               # fix invalid geometries
  st_cast("MULTIPOLYGON") %>%                       # ensure consistent geometry type
  st_collection_extract("POLYGON") %>%              # drop garbage non-polygon bits
  st_transform(AEA_proj) %>%                             # reproject to a metric CRS (metres)
  mutate(area_m2 = as.numeric(st_area(geometry))) %>%
  subset(area_m2 >= 1000000) # Only include waterbodies larger than 1000m x 1000m

st_write(water, "data_clean/spatial/water_filtered.shp", delete_layer = TRUE)

# ============================================================================
# ---- Road Network (polygons)
# ============================================================================

# Crop to study area (ON_boundary_25km must be vect, not sf)
ON_boundary_25km_v <- vect(ON_boundary_25km)
study_ext <- ext(ON_boundary_25km_v)

# 1. Read the roads directly into terra
road_2025 <- vect("data/spatial/RoadNetwork/2025/lrnf000r25a_e.shp")
road_2005 <- vect("data/spatial/RoadNetwork/2005/grnf035r05a_e.shp")

# 2. Reproject to study CRS (AEA_proj) and select roads in study area
road_2025 <- project(road_2025, st_crs(AEA_proj)$wkt)
road_2005 <- project(road_2005, st_crs(AEA_proj)$wkt)

road_2025 <- road_2025[intersect(road_2025, study_ext), ]
road_2005 <- road_2005[intersect(road_2005, study_ext), ]

# 3. Buffer 100 m
road_2025_buf <- buffer(road_2025, width = 100)
road_2005_buf <- buffer(road_2005, width = 100)

# Dissolve overlapping polygons into one
road_2025_buf <- aggregate(road_2025_buf, dissolve = TRUE)
road_2005_buf <- aggregate(road_2005_buf, dissolve = TRUE)

# 4. Make a raster template at 30 m resolution
r_template <- rast(ext(ON_boundary_25km_v), 
                   res = 30, 
                   crs = crs(ON_boundary_25km_v))

# Mask template to study area
values(r_template) <- 0
r_template <- mask(r_template, ON_boundary_25km_v)

# 5. Rasterize buffer polygons (1 = roadside, 0 = not roadside)
roadside_2025 <- rasterize(road_2025_buf, r_template, field = 1, background = 0)
roadside_2005 <- rasterize(road_2005_buf, r_template, field = 1, background = 0)

# 6. Save rasters
writeRaster(roadside_2025, "data_clean/spatial/roadside_2025_100m.tif", overwrite = TRUE)
writeRaster(roadside_2005, "data_clean/spatial/roadside_2005_100m.tif", overwrite = TRUE)

# ============================================================================
# ---- Insect damage, from https://geohub.lio.gov.on.ca/documents/308957e79fce445fa74997c4cf501437/about
# ============================================================================

# Create 100m raster across the study area
r_template <- rast(ext(ON_boundary_25km_v), 
                   res = 100, 
                   crs = crs(ON_boundary_25km_v))

# Mask template to study area
values(r_template) <- 0
r_template <- mask(r_template, ON_boundary_25km_v)

ON_boundary_25km_v <- vect(ON_boundary_25km)
study_ext <- ext(ON_boundary_25km_v)

insect <- st_read("data/spatial/Insect_Damage/FOREST_INSECT_DAMAGE_EVENT.shp") %>%
  st_transform(st_crs(AEA_proj)) %>%
  subset(INSECT %in% c("Gypsy Moth","Spruce Budworm","Forest Tent Caterpillar","Jack Pine Budworm"))
  
insect_Atlas2 <- subset(insect, EVENT_YEAR %in% c(2000:2004)) %>% 
  subset(RANKING != "Light")
# Create 100m raster across the study area
r_template <- rast(ext(ON_boundary_25km_v), 
                   res = 100, 
                   crs = crs(ON_boundary_25km_v))

# Mask template to study area
values(r_template) <- 0
r_template <- mask(r_template, ON_boundary_25km_v)

ON_boundary_25km_v <- vect(ON_boundary_25km)
study_ext <- ext(ON_boundary_25km_v)

insect <- st_read("data/spatial/Insect_Damage/FOREST_INSECT_DAMAGE_EVENT.shp") %>%
  st_transform(st_crs(AEA_proj)) %>%
  subset(INSECT %in% c("Gypsy Moth","Spruce Budworm","Forest Tent Caterpillar","Jack Pine Budworm"))
  
insect_Atlas2 <- subset(insect, EVENT_YEAR %in% c(2000:2004)) %>% 
  subset(RANKING != "Light") %>%
  st_simplify(dTolerance = 100)
insect_Atlas3 <- subset(insect, EVENT_YEAR %in% c(2020:2024)) %>% 
  subset(RANKING != "Light") %>%
  st_simplify(dTolerance = 100)

# Function to union polygons and rasterize by insect type
rasterize_insects <- function(insect_sf, r_template) {
  
  # Union polygons by insect type
  insect_union <- insect_sf %>%
    group_by(INSECT) %>%
    summarise(geometry = st_union(geometry), .groups = "drop")
  
  #  Rasterize each insect
  r_stack <- rast(r_template)  # copy template
  
  for (i in seq_len(nrow(insect_union))) {
    insect_name <- insect_union$INSECT[i]
    poly <- insect_union[i, ]
    
    # rasterize: 1 where insect occurs, 0 elsewhere
    r_tmp <- rasterize(vect(poly), r_template, field = 1, background = 0)
    names(r_tmp) <- insect_name
    
    if (i == 1) {
      r_stack <- r_tmp
    } else {
      r_stack <- c(r_stack, r_tmp)
    }
  }
  
  return(r_stack)
}

# Atlas 2000-2004
r_Atlas2 <- rasterize_insects(insect_Atlas2, r_template)

# Atlas 2020-2024
r_Atlas3 <- rasterize_insects(insect_Atlas3, r_template)
insect_Atlas3 <- subset(insect, EVENT_YEAR %in% c(2020:2024)) %>% subset(RANKING != "Light")

# Function to union polygons and rasterize by insect type
rasterize_insects <- function(insect_sf, r_template) {
  
  # Union polygons by insect type
  insect_union <- insect_sf %>%
    group_by(INSECT) %>%
    summarise(geometry = st_union(geometry), .groups = "drop")
  
  #  Rasterize each insect
  r_stack <- rast(r_template)  # copy template
  
  for (i in seq_len(nrow(insect_union))) {
    insect_name <- insect_union$INSECT[i]
    poly <- insect_union[i, ]
    
    # rasterize: 1 where insect occurs, 0 elsewhere
    r_tmp <- rasterize(vect(poly), r_template, field = 1, background = 0)
    names(r_tmp) <- insect_name
    
    if (i == 1) {
      r_stack <- r_tmp
    } else {
      r_stack <- c(r_stack, r_tmp)
    }
  }
  
  return(r_stack)
}

# Atlas 2000-2004
r_Atlas2 <- rasterize_insects(insect_Atlas2, r_template)

# Atlas 2020-2024
r_Atlas3 <- rasterize_insects(insect_Atlas3, r_template)

# Save rasters
writeRaster(r_Atlas2, "data_clean/spatial/insect_OBBA2.tif", overwrite = TRUE)
writeRaster(r_Atlas3, "data_clean/spatial/insect_OBBA3.tif", overwrite = TRUE)
