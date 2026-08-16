# ============================================================
# 04_process_vectors_as_rasters.R
#
# Strategy
#   - Rasterize all vector layers at 30 m resolution (binary 0/1).
#   - Percent cover per 1-km pixel is later computed as the mean of the 30 m
#     cells (exactextractr::exact_extract(..., "mean") in 05).
#   - Simplify geometries before buffering to cut vertex counts.
#   - Erase overlaps (e.g. small rivers minus large rivers) at the RASTER level
#     with terra::mask(), avoiding vector dissolve/difference complexity.
#   - Outputs are DEFLATE-compressed INT1U GeoTiffs.
#
# Outputs (data_clean/spatial/)
#   roads_2005_buf_125m.tif, roads_2025_buf_125m.tif, coastline_buf_500m.tif,
#   waterbodies.tif, rivers_large_buf_250m.tif, rivers_small_buf_125m.tif,
#   great_lakes_buf_500m.tif, lakes_small_buf_250m.tif, lakes_large_buf_250m.tif
# ============================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(terra)
  library(here)
  library(rnaturalearth)
})

source(here::here("R", "00_config_paths.R"))
source(file.path(paths$functions, "covariate_processing_utils.R"))

out_spatial_dir <- file.path(paths$data_clean, "spatial")
dir.create(out_spatial_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Config
# ------------------------------------------------------------

study_area_path  <- file.path(paths$data_clean, "spatial", "study_area.rds")

roads_2005_path  <- file.path(paths$data, "Spatial", "RoadNetwork", "2005", "grnf035r05a_e.shp")
roads_2025_path  <- file.path(paths$data, "Spatial", "RoadNetwork", "2025", "lrnf000r25a_e.shp")
waterbody_path   <- file.path(paths$data, "Spatial", "Ontario_Hydro_Network_(OHN)_-_Waterbody",
                              "Ontario_Hydro_Network_(OHN)_-_Waterbody.shp")
watercourse_path <- file.path(paths$data, "Spatial", "Ontario_Hydro_Network_(OHN)_-_Watercourse",
                              "Ontario_Hydro_Network_(OHN)_-_Watercourse.shp")

# Simplification tolerance before buffering (CRS units = km). 10 m removes
# redundant vertices without visible change at 1-km pixel resolution.
SIMPLIFY_KM <- 0.01

# Raster resolution (CRS units = km). 30 m.
RES_KM <- 0.03

# GeoTiff write options (applied to every output raster).
GTIFF_OPTS <- c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=6")

# ------------------------------------------------------------
# Study area and 30 m template raster
# ------------------------------------------------------------

study_area   <- readRDS(study_area_path)
boundary_buf <- study_area$boundary_buffer_5km
study_crs    <- sf::st_crs(boundary_buf)

# Study-area boundary as lines (used to identify coastline segments).
boundary_line <- study_area$boundary %>%
  st_make_valid() %>%
  st_union() %>%
  st_boundary()

# Single template reused for every rasterize() call, so all outputs share extent,
# resolution and CRS. `template_30m` / `SIMPLIFY_KM` are read by buf_rasterize_save().
template_30m <- terra::rast(
  terra::vect(boundary_buf), resolution = RES_KM, crs = study_crs$wkt, vals = 0L
)

message(sprintf("Template raster: %d rows x %d cols  (%.0f m resolution)",
                nrow(template_30m), ncol(template_30m), RES_KM * 1000))

# ------------------------------------------------------------
# Roads (buffer = 125 m)
# ------------------------------------------------------------

message("\n-- Roads --------------------------------------------------")

roads_2005 <- sf::st_read(roads_2005_path, quiet = TRUE) %>% sf::st_transform(study_crs)
buf_rasterize_save(roads_2005, buffer_dist_km = 0.125,
                   out_path = file.path(out_spatial_dir, "roads_2005_buf_125m.tif"),
                   label = "roads_2005")
rm(roads_2005); gc()

roads_2025 <- sf::st_read(roads_2025_path, quiet = TRUE) %>% sf::st_transform(study_crs)
buf_rasterize_save(roads_2025, buffer_dist_km = 0.125,
                   out_path = file.path(out_spatial_dir, "roads_2025_buf_125m.tif"),
                   label = "roads_2025")
rm(roads_2025); gc()

# ------------------------------------------------------------
# Hudson Bay coastline (buffer = 500 m)
# ------------------------------------------------------------

message("\n-- Coastline ----------------------------------------------")

# 1. Natural Earth coastline, clipped to the buffered study area and to >= 50 N.
coastline_ne <- rnaturalearth::ne_coastline(scale = 10) %>%
  st_transform(study_crs) %>%
  st_intersection(boundary_buf) %>%
  st_transform(4326) %>%
  st_crop(
    xmin = st_bbox(.)[["xmin"]], ymin = 50,
    xmax = st_bbox(.)[["xmax"]], ymax = st_bbox(.)[["ymax"]]
  ) %>%
  st_transform(study_crs) %>%
  st_union() %>%
  st_make_valid()

# 2. Buffer used only to pick out study-boundary segments near the NE coastline.
coastline_search_buf <- coastline_ne %>%
  st_buffer(dist = 5) %>%
  st_union() %>%
  st_make_valid()

# 3. Study-area boundary segments within 5 km of the NE Hudson Bay coastline.
coastline_boundary <- boundary_line %>%
  st_intersection(coastline_search_buf) %>%
  st_collection_extract("LINESTRING") %>%
  st_union() %>%
  st_make_valid()

# Trim a spurious "coastline" stretch at the SE corner of James Bay
# (~50.9589 N, -79.5222 W). Here the Ontario-Quebec land border runs within
# 5 km of the NE coastline and gets picked up as coast, but it is not shoreline.
trim_zone <- sf::st_sfc(sf::st_point(c(-79.5222, 50.9589)), crs = 4326) |>
  sf::st_transform(sf::st_crs(coastline_boundary)) |>
  sf::st_buffer(dist = 15)          # <- radius in KILOMETRES (study CRS is +units=km)

coastline_boundary <- coastline_boundary |>
  sf::st_difference(trim_zone) |>
  sf::st_collection_extract("LINESTRING") |>
  sf::st_make_valid()

# 4. Final 500-m coastal covariate from the boundary-derived coastline.
buf_rasterize_save(
  coastline_boundary, buffer_dist_km = 0.5,
  out_path = file.path(out_spatial_dir, "coastline_buf_500m.tif"),
  label = "coastline"
)

# ------------------------------------------------------------
# Load hydro source layers (once; filtered below)
# ------------------------------------------------------------

message("\n-- Loading hydro layers -----------------------------------")

watercourse <- sf::st_read(watercourse_path, quiet = TRUE) %>%
  sf::st_transform(study_crs)

waterbody <- sf::st_read(waterbody_path, quiet = TRUE) %>%
  sf::st_transform(study_crs) %>%
  dplyr::filter(VERIFICATI == "Verified")

# ------------------------------------------------------------
# All waterbodies (so open water can be removed from predictions)
# ------------------------------------------------------------

buf_rasterize_save(
  waterbody, buffer_dist_km = 0,
  out_path = file.path(out_spatial_dir, "waterbodies.tif"),
  label = "waterbodies"
)

# ------------------------------------------------------------
# Rivers
#   large rivers : waterbody polygons (River/Canal), buffered 250 m
#   small rivers : watercourse lines (Stream),       buffered 125 m,
#                  then large-river cells zeroed out
# ------------------------------------------------------------

message("\n-- Rivers -------------------------------------------------")

large_rivers <- waterbody %>%
  dplyr::filter(WATERBODY_ %in% c("River", "Canal"), PERMANENCY == "Permanent") %>%
  sf::st_make_valid()

large_rivers_r <- buf_rasterize_save(
  large_rivers, buffer_dist_km = 0.25,
  out_path = file.path(out_spatial_dir, "rivers_large_buf_250m.tif"),
  label = "large_rivers"
)

small_rivers <- watercourse %>%
  dplyr::filter(WATERCOURS == "Stream", PERMANENCY == "Permanent",
                FLOW_CLASS != "Flow Gap") %>%
  sf::st_make_valid()

# Rasterize small rivers, then erase cells already covered by large rivers.
small_rivers_r <- buf_rasterize_save(
  small_rivers, buffer_dist_km = 0.125,
  out_path = NULL, label = "small_rivers"        # saved after masking
) %>%
  terra::mask(large_rivers_r, maskvalues = 1L, updatevalue = 0L)

terra::writeRaster(small_rivers_r,
                   file.path(out_spatial_dir, "rivers_small_buf_125m.tif"),
                   datatype = "INT1U", overwrite = TRUE, gdal = GTIFF_OPTS)
message("[small_rivers] Saved -> rivers_small_buf_125m.tif")

rm(large_rivers, large_rivers_r, small_rivers, small_rivers_r); gc()

# ------------------------------------------------------------
# Lakes -- base dataset (shared across all lake layers)
# ------------------------------------------------------------

message("\n-- Lakes base dataset -------------------------------------")

lake_types <- c("Lake", "Kettle Lake", "Pond", "Reservoir", "Beaver Pond")

lakes_base <- waterbody %>%
  dplyr::filter(WATERBODY_ %in% lake_types, PERMANENCY == "Permanent") %>%
  sf::st_make_valid() %>%
  sf::st_filter(study_area$boundary_buffer_5km) %>%
  dplyr::mutate(area_km2 = as.numeric(sf::st_area(geometry))) %>%
  dplyr::arrange(dplyr::desc(area_km2))

# Remove ocean/Hudson Bay polygons that abut the coastline. The exclusion zone is
# the offshore side of a 10-km buffer of the Natural Earth Hudson Bay coastline
# (`coastline_ne`), unioned with a 100-m buffer of the coastline itself.
# (Previously this referenced an undefined `coastline`; `coastline_ne` is the
# unioned NE coastline built above.)
coast_buf_10km_outside <- sf::st_difference(
  sf::st_buffer(coastline_ne, dist = 10), study_area$boundary
)
coast_buf_100m  <- sf::st_buffer(coastline_ne, dist = 0.1)
coast_exclusion <- sf::st_union(
  c(sf::st_geometry(coast_buf_10km_outside), sf::st_geometry(coast_buf_100m))
) %>%
  sf::st_make_valid()

keep       <- sf::st_disjoint(lakes_base, coast_exclusion, sparse = FALSE)[, 1]
lakes_base <- lakes_base[keep, ]

rm(coast_buf_10km_outside, coast_buf_100m, coast_exclusion); gc()

# ------------------------------------------------------------
# Great Lakes (buffer = 500 m)
# ------------------------------------------------------------

message("\n-- Great Lakes --------------------------------------------")

GL_names <- c(
  "Lake Superior", "Lake Superior (lac Sup\u00e9rieur)",
  "Lake Michigan",
  "Lake Huron", "Lake Huron (lac Huron)", "Georgian Bay (baie Georgienne)", "Georgian Bay",
  "Lake Erie", "Lake Erie (lac \u00c9ri\u00e9)",
  "Lake Ontario", "Lake Ontario (lac Ontario)"
)

great_lakes_ne <- rnaturalearth::ne_download(
  scale = 10, type = "lakes", category = "physical", returnclass = "sf"
) %>%
  sf::st_transform(study_crs) %>%
  dplyr::filter(name %in% c("Lake Superior", "Lake Michigan", "Lake Huron",
                            "Lake Erie", "Lake Ontario")) %>%
  sf::st_intersection(boundary_buf) %>%
  sf::st_make_valid()

waterbody_GL <- lakes_base %>%
  dplyr::filter(OFFICIAL_N %in% GL_names)

great_lakes_combined <- bind_rows(great_lakes_ne, waterbody_GL)

# Study-boundary segments within 5 km of the combined Great Lakes shoreline.
great_lakes_search_buf <- great_lakes_combined %>%
  st_buffer(dist = 5) %>%
  st_union() %>%
  st_make_valid()

great_lakes_boundary <- boundary_line %>%
  st_intersection(great_lakes_search_buf) %>%
  st_collection_extract("LINESTRING") %>%
  st_union() %>%
  st_make_valid()

buf_rasterize_save(
  great_lakes_boundary, buffer_dist_km = 0.5,
  out_path = file.path(out_spatial_dir, "great_lakes_buf_500m.tif"),
  label = "great_lakes"
)

rm(great_lakes_ne, waterbody_GL, great_lakes_combined,
   great_lakes_search_buf, great_lakes_boundary); gc()

# ------------------------------------------------------------
# Non-GL lakes: small (< 1 km2) and large (>= 1 km2), both buffered 250 m
# ------------------------------------------------------------

message("\n-- Non-GL lakes -------------------------------------------")

large_lake_area_threshold_km2 <- 1

lakes_nonGL <- lakes_base %>% dplyr::filter(!(OFFICIAL_N %in% GL_names))

small_lakes <- lakes_nonGL %>% dplyr::filter(area_km2 < large_lake_area_threshold_km2)
buf_rasterize_save(small_lakes, buffer_dist_km = 0.25,
                   out_path = file.path(out_spatial_dir, "lakes_small_buf_250m.tif"),
                   label = "small_lakes")

large_lakes <- lakes_nonGL %>% dplyr::filter(area_km2 >= large_lake_area_threshold_km2)
buf_rasterize_save(large_lakes, buffer_dist_km = 0.25,
                   out_path = file.path(out_spatial_dir, "lakes_large_buf_250m.tif"),
                   label = "large_lakes")

rm(lakes_base, lakes_nonGL, small_lakes, large_lakes); gc()

message("\n04_process_vectors_as_rasters.R complete.")
message("Outputs: binary INT1U GeoTiffs at 30 m resolution.")
