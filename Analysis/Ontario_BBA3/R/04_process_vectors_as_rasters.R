# ============================================================
# 04_process_vector_layers.R
# ============================================================
# Strategy:
#   - Rasterize all layers at 30 m resolution (binary 0/1)
#   - Percent cover per 1-km pixel = mean of 30 m cells
#     (use exactextractr::exact_extract(..., fun = "mean"))
#   - Simplify geometries before buffering to reduce vertex counts
#   - Erasure (e.g. small rivers minus large rivers) is done at
#     the raster level with terra::mask() — avoids all vector
#     dissolve/difference complexity
#   - Outputs are GeoTiff with DEFLATE compression (~INT1U = 1 bit)
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
source(file.path(paths$functions, "spatial_utils.R"))
source(file.path(paths$functions, "covariate_processing_utils.R"))

out_spatial_dir <- file.path(paths$data_clean, "spatial")
dir.create(out_spatial_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Config
# ------------------------------------------------------------

study_area_path  <- file.path(paths$data_clean, "spatial", "study_area.rds")
grid_path        <- file.path(out_spatial_dir, "prediction_grid_1000_m.rds")

roads_2005_path  <- file.path(paths$data, "Spatial", "RoadNetwork", "2005", "grnf035r05a_e.shp")
roads_2025_path  <- file.path(paths$data, "Spatial", "RoadNetwork", "2025", "lrnf000r25a_e.shp")
waterbody_path   <- file.path(paths$data, "Spatial", "Ontario_Hydro_Network_(OHN)_-_Waterbody",
                              "Ontario_Hydro_Network_(OHN)_-_Waterbody.shp")
watercourse_path <- file.path(paths$data, "Spatial", "Ontario_Hydro_Network_(OHN)_-_Watercourse",
                              "Ontario_Hydro_Network_(OHN)_-_Watercourse.shp")

# Simplification tolerance before buffering (CRS units = km)
SIMPLIFY_KM  <- 0.01    # 10 m — removes redundant vertices without
# visible change at 1-km pixel resolution

# Raster resolution (CRS units = km)
RES_KM       <- 0.03   # 30 m

# GeoTiff write options (applied to every output raster)
GTIFF_OPTS   <- c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=6")

# ------------------------------------------------------------
# Study area & 30 m template raster
# ------------------------------------------------------------

study_area   <- readRDS(study_area_path)
boundary_buf <- study_area$boundary_buffer_5km
study_crs    <- sf::st_crs(boundary_buf)

# Single template reused for every rasterize() call.
# All output rasters will share the same extent, resolution, and CRS.
template_30m <- terra::rast(
  terra::vect(boundary_buf),
  resolution = RES_KM,
  crs        = study_crs$wkt,
  vals       = 0L
)

message(sprintf("Template raster: %d rows x %d cols  (%.0f m resolution)",
                nrow(template_30m), ncol(template_30m), RES_KM * 1000))

# ------------------------------------------------------------
# Helper: simplify → buffer → rasterize → save GeoTiff
# ------------------------------------------------------------
# sf_obj        : sf object (lines or polygons)
# buffer_dist_km: buffer radius in km (matches CRS units)
# out_path      : full .tif path, or NULL to skip writing
# label         : used in progress messages
# Returns        : SpatRaster (INT1U, 0/1)

buf_rasterize_save <- function(sf_obj, buffer_dist_km, out_path, label,
                               template     = template_30m,
                               simplify_tol = SIMPLIFY_KM) {
  
  message(sprintf("[%s] %d features | simplify %.0f m | buffer %.0f m | rasterizing...",
                  label, nrow(sf_obj), simplify_tol * 1000, buffer_dist_km * 1000))
  
  v <- sf_obj %>%
    sf::st_simplify(preserveTopology = TRUE, dTolerance = simplify_tol) %>%
    sf::st_make_valid() %>%
    sf::st_buffer(dist = buffer_dist_km) %>%
    sf::st_make_valid() %>%
    terra::vect()
  
  r <- terra::rasterize(v, template, field = 1L,
                        background = 0L, touches = TRUE)
  
  if (!is.null(out_path)) {
    terra::writeRaster(r, out_path, datatype = "INT1U",
                       overwrite = TRUE, gdal = GTIFF_OPTS)
    message(sprintf("[%s] Saved -> %s", label, basename(out_path)))
  }
  
  invisible(r)
}

# ------------------------------------------------------------
# Roads  (buffer = 125 m)
# ------------------------------------------------------------

message("\n── Roads ───────────────────────────────────────────────")

roads_2005 <- sf::st_read(roads_2005_path, quiet = TRUE) %>%
  sf::st_transform(study_crs)
buf_rasterize_save(roads_2005, buffer_dist_km = 0.125,
                   out_path = file.path(out_spatial_dir, "roads_2005_buf_125m.tif"),
                   label    = "roads_2005")
rm(roads_2005); gc()

roads_2025 <- sf::st_read(roads_2025_path, quiet = TRUE) %>%
  sf::st_transform(study_crs)
buf_rasterize_save(roads_2025, buffer_dist_km = 0.125,
                   out_path = file.path(out_spatial_dir, "roads_2025_buf_125m.tif"),
                   label    = "roads_2025")
rm(roads_2025); gc()

# ------------------------------------------------------------
# Hudson Bay coastline  (buffer = 500 m)
# ------------------------------------------------------------

message("\n── Coastline ───────────────────────────────────────────")

coastline <- rnaturalearth::ne_coastline(scale = 10) %>%
  sf::st_transform(study_crs) %>%
  sf::st_intersection(boundary_buf)

buf_rasterize_save(coastline, buffer_dist_km = 0.5,
                   out_path = file.path(out_spatial_dir, "coastline_buf_500m.tif"),
                   label    = "coastline")

# ------------------------------------------------------------
# Load hydro source layers (once; filtered below)
# ------------------------------------------------------------

message("\n── Loading hydro layers ────────────────────────────────")

watercourse <- sf::st_read(watercourse_path, quiet = TRUE) %>%
  sf::st_transform(study_crs)

waterbody <- sf::st_read(waterbody_path, quiet = TRUE) %>%
  sf::st_transform(study_crs) %>%
  dplyr::filter(VERIFICATI == "Verified")

# ------------------------------------------------------------
# Rasterize all waterbodies so open water can be removed from predictions 
# ------------------------------------------------------------
waterbody_raster <- buf_rasterize_save(
  waterbody, buffer_dist_km = 0,
  out_path = file.path(out_spatial_dir, "waterbodies.tif"),
  label    = "waterbodies"
)

# ------------------------------------------------------------
# Rivers
#   large rivers : waterbody polygons, buffered 250 m
#   small rivers : watercourse lines, buffered 125 m,
#                  then large-river cells zeroed out
# ------------------------------------------------------------

message("\n── Rivers ──────────────────────────────────────────────")

large_rivers <- waterbody %>%
  dplyr::filter(WATERBODY_ %in% c("River", "Canal"), PERMANENCY == "Permanent") %>%
  sf::st_make_valid()

large_rivers_r <- buf_rasterize_save(
  large_rivers, buffer_dist_km = 0.25,
  out_path = file.path(out_spatial_dir, "rivers_large_buf_250m.tif"),
  label    = "large_rivers"
)

small_rivers <- watercourse %>%
  dplyr::filter(WATERCOURS == "Stream", PERMANENCY == "Permanent",
                FLOW_CLASS != "Flow Gap") %>%
  sf::st_make_valid()

# Rasterize small rivers, then erase cells already covered by large rivers.
# terra::mask() with maskvalues = 1 sets those cells to updatevalue = 0.
small_rivers_r <- buf_rasterize_save(
  small_rivers, buffer_dist_km = 0.125,
  out_path = NULL,          # save after masking
  label    = "small_rivers"
) %>%
  terra::mask(large_rivers_r, maskvalues = 1L, updatevalue = 0L)

terra::writeRaster(small_rivers_r,
                   file.path(out_spatial_dir, "rivers_small_buf_125m.tif"),
                   datatype = "INT1U", overwrite = TRUE, gdal = GTIFF_OPTS)
message("[small_rivers] Saved -> rivers_small_buf_125m.tif")

rm(large_rivers, large_rivers_r, small_rivers, small_rivers_r); gc()

# ------------------------------------------------------------
# Lakes — base dataset (shared across all lake layers)
# ------------------------------------------------------------

message("\n── Lakes base dataset ──────────────────────────────────")

lake_types <- c("Lake", "Kettle Lake", "Pond", "Reservoir", "Beaver Pond")

lakes_base <- waterbody %>%
  dplyr::filter(WATERBODY_ %in% lake_types, PERMANENCY == "Permanent") %>%
  sf::st_make_valid() %>%
  sf::st_filter(study_area$boundary_buffer_5km) %>%
  dplyr::mutate(area_km2 = as.numeric(sf::st_area(geometry))) %>%
  dplyr::arrange(dplyr::desc(area_km2))

# Remove ocean/Hudson Bay polygons that abut the coastline
coast_buf_10km         <- sf::st_buffer(coastline, dist = 10)
coast_buf_10km_outside <- sf::st_difference(coast_buf_10km, study_area$boundary)
coast_buf_100m         <- sf::st_buffer(coastline, dist = 0.1)
coast_exclusion        <- dplyr::bind_rows(coast_buf_10km_outside, coast_buf_100m) %>%
  sf::st_union() %>%
  sf::st_make_valid()

keep       <- sf::st_disjoint(lakes_base, coast_exclusion, sparse = FALSE)[, 1]
lakes_base <- lakes_base[keep, ]

rm(coast_buf_10km, coast_buf_10km_outside, coast_buf_100m, coast_exclusion); gc()

# ------------------------------------------------------------
# Great Lakes  (buffer = 500 m)
# ------------------------------------------------------------

message("\n── Great Lakes ─────────────────────────────────────────")

GL_names <- c(
  "Lake Superior", "Lake Superior (lac Supérieur)",
  "Lake Michigan",
  "Lake Huron", "Lake Huron (lac Huron)", "Georgian Bay (baie Georgienne)", "Georgian Bay",
  "Lake Erie", "Lake Erie (lac Érié)",
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

# Rasterize the two Great Lakes sources separately, then combine
# with terra::merge() (union of two binary rasters).
great_lakes_ne_r <- buf_rasterize_save(
  great_lakes_ne, buffer_dist_km = 0.5,
  out_path = NULL,
  label    = "great_lakes_ne"
)

great_lakes_wb_r <- buf_rasterize_save(
  waterbody_GL, buffer_dist_km = 0.5,
  out_path = NULL,
  label    = "great_lakes_waterbody"
)

# Union: cell = 1 if either source has a 1
great_lakes_r <- terra::lapp(
  c(great_lakes_ne_r, great_lakes_wb_r),
  fun = function(a, b) pmax(a, b)
)

terra::writeRaster(great_lakes_r,
                   file.path(out_spatial_dir, "great_lakes_buf_500m.tif"),
                   datatype = "INT1U", overwrite = TRUE, gdal = GTIFF_OPTS)
message("[great_lakes] Saved -> great_lakes_buf_500m.tif")

rm(great_lakes_ne, waterbody_GL, great_lakes_ne_r, great_lakes_wb_r, great_lakes_r); gc()

# ------------------------------------------------------------
# Non-GL lakes: small (< 1 km²) and large (>= 1 km²)
# Both buffered 250 m — no erasure needed between them
# ------------------------------------------------------------

message("\n── Non-GL lakes ────────────────────────────────────────")

large_lake_area_threshold_km2 <- 1

lakes_nonGL <- lakes_base %>%
  dplyr::filter(!(OFFICIAL_N %in% GL_names))

small_lakes <- lakes_nonGL %>% dplyr::filter(area_km2 <  large_lake_area_threshold_km2)
buf_rasterize_save(small_lakes, buffer_dist_km = 0.25,
                   out_path = file.path(out_spatial_dir, "lakes_small_buf_250m.tif"),
                   label    = "small_lakes")

large_lakes <- lakes_nonGL %>% dplyr::filter(area_km2 >= large_lake_area_threshold_km2)
buf_rasterize_save(large_lakes, buffer_dist_km = 0.25,
                   out_path = file.path(out_spatial_dir, "lakes_large_buf_250m.tif"),
                   label    = "large_lakes")

rm(lakes_base, lakes_nonGL, small_lakes, large_lakes); gc()

# ------------------------------------------------------------

message("\n04_process_vector_layers.R complete.")
message("Outputs: binary INT1U GeoTiffs at 30 m resolution.")
message("Percent cover per 1-km pixel: exactextractr::exact_extract(..., fun = 'mean')")