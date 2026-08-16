# ============================================================
# 01_define_study_area.R
#
# Purpose
#   Build the Ontario study boundary (Ontario + Akimiski Island), a 5-km buffer
#   of it, and simplified Bird Conservation Region (BCR) polygons clipped to the
#   boundary. Saves them as a single study_area object for the rest of the
#   pipeline.
#
# Output
#   data_clean/spatial/study_area.rds
#     list(crs, boundary, boundary_buffer_5km, bcr)
# ============================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(here)
  library(smoothr)
})

source(here::here("R", "00_config_paths.R"))
source(file.path(paths$functions, "spatial_utils.R"))

crs_aea_km <- get_aea_km_crs()

# Simplification tolerances, in kilometres.
boundary_tolerance_km <- 0.1
bcr_tolerance_km      <- 0.1

# ------------------------------------------------------------
# Load full-resolution boundary data
# ------------------------------------------------------------

bcr_path <- file.path(
  paths$data, "Spatial", "National", "BCR_2026", "bcr2026_statprov.shp"
)

bcr_raw <- st_read(bcr_path, quiet = TRUE) %>%
  st_transform(crs_aea_km) %>%
  st_make_valid()

# ------------------------------------------------------------
# Extract Ontario at full resolution
# ------------------------------------------------------------

ontario_raw <- bcr_raw %>%
  filter(NAME_En == "Ontario") %>%
  summarise(geometry = st_union(geometry), .groups = "drop") %>%
  st_make_valid()

# ------------------------------------------------------------
# Extract Akimiski Island (a Nunavut polygon in James Bay) at full resolution
# ------------------------------------------------------------

nu_raw <- bcr_raw %>%
  filter(NAME_En == "Nunavut")

akimiski_pt <- st_as_sf(
  data.frame(lon = -81.363111, lat = 53.024846),
  coords = c("lon", "lat"),
  crs = 4326
) %>%
  st_transform(crs_aea_km)

nu_parts <- nu_raw %>%
  st_collection_extract("POLYGON") %>%
  st_cast("POLYGON", warn = FALSE) %>%
  mutate(part_id = row_number())

akimiski <- nu_parts[lengths(st_intersects(nu_parts, akimiski_pt)) > 0, ]

if (nrow(akimiski) != 1) {
  stop(
    "Expected exactly one Nunavut polygon to intersect the Akimiski point; found ",
    nrow(akimiski), "."
  )
}

# Keep only the geometry so Nunavut attributes are not carried forward.
akimiski <- akimiski %>% select(geometry)

# ------------------------------------------------------------
# Combine Ontario and Akimiski
# ------------------------------------------------------------

boundary_raw <- st_sf(
  geometry = st_union(c(st_geometry(ontario_raw), st_geometry(akimiski))),
  crs = crs_aea_km
) %>%
  st_make_valid()

# ------------------------------------------------------------
# Clean and simplify the final boundary
# ------------------------------------------------------------

boundary <- boundary_raw %>%
  fill_holes(threshold = units::set_units(10, km^2)) %>%
  drop_crumbs(threshold = units::set_units(1, km^2)) %>%
  st_simplify(dTolerance = boundary_tolerance_km, preserveTopology = TRUE) %>%
  st_make_valid() %>%
  # Re-union in case cleaning/validation returned multiple features.
  summarise(geometry = st_union(geometry), .groups = "drop") %>%
  st_make_valid()

# ------------------------------------------------------------
# 5-km study buffer
# ------------------------------------------------------------

boundary_buffer_5km <- boundary %>%
  st_buffer(dist = 5, nQuadSegs = 8) %>%
  st_make_valid() %>%
  st_simplify(dTolerance = boundary_tolerance_km, preserveTopology = TRUE) %>%
  summarise(geometry = st_union(geometry), .groups = "drop") %>%
  st_make_valid()

# ------------------------------------------------------------
# Simplified BCR polygons within the study boundary
# ------------------------------------------------------------

# Restrict the input BCR layer to features that may intersect the study area
# first, so st_intersection() is not run on irrelevant geometries.
bcr_candidates <- bcr_raw[lengths(st_intersects(bcr_raw, boundary)) > 0, ]

ontario_bcr <- st_intersection(bcr_candidates, boundary) %>%
  st_collection_extract("POLYGON") %>%
  st_make_valid()

# Replace BCR 23 attributes with BCR 13 attributes (23 is folded into 13).
bcr13_attrs <- ontario_bcr %>%
  st_drop_geometry() %>%
  filter(bcr == 13) %>%
  slice(1) %>%
  select(-any_of("fid"))

if (nrow(bcr13_attrs) != 1) {
  stop("Could not identify exactly one set of attributes for BCR 13.")
}

cols_to_overwrite <- names(bcr13_attrs)

ontario_bcr <- ontario_bcr %>%
  mutate(
    across(
      all_of(cols_to_overwrite),
      ~ if_else(bcr == 23, bcr13_attrs[[cur_column()]], .)
    )
  )

# Dissolve to one polygon per BCR, then simplify.
bcr <- ontario_bcr %>%
  group_by(bcr, bcr_label, bcr_label_, bcr_name_e, bcr_name_f, bcr_name_1) %>%
  summarise(geometry = st_union(geometry), .groups = "drop") %>%
  st_make_valid() %>%
  st_simplify(dTolerance = bcr_tolerance_km, preserveTopology = TRUE) %>%
  st_make_valid()

# ------------------------------------------------------------
# Sanity checks
# ------------------------------------------------------------

if (st_is_longlat(boundary)) {
  stop("Study boundary is in longitude/latitude; expected a projected CRS.")
}
if (st_crs(boundary) != crs_aea_km) {
  stop("Study boundary does not use the expected CRS.")
}
if (st_crs(boundary_buffer_5km) != crs_aea_km) {
  stop("Buffered boundary does not use the expected CRS.")
}
if (st_crs(bcr) != crs_aea_km) {
  stop("BCR polygons do not use the expected CRS.")
}

stopifnot(all(st_is_valid(boundary)))
stopifnot(all(st_is_valid(boundary_buffer_5km)))
stopifnot(all(st_is_valid(bcr)))

# Confirm Akimiski survived cleaning and simplification.
if (!any(st_intersects(akimiski_pt, boundary, sparse = FALSE))) {
  stop("Akimiski Island was lost while cleaning or simplifying the boundary.")
}

# ------------------------------------------------------------
# Geometry-complexity diagnostics
# ------------------------------------------------------------

count_vertices <- function(x) nrow(st_coordinates(st_geometry(x)))

complexity_summary <- tibble(
  object = c("Raw boundary", "Simplified boundary",
             "5-km boundary buffer", "Final BCR polygons"),
  n_features = c(nrow(boundary_raw), nrow(boundary),
                 nrow(boundary_buffer_5km), nrow(bcr)),
  n_vertices = c(count_vertices(boundary_raw), count_vertices(boundary),
                 count_vertices(boundary_buffer_5km), count_vertices(bcr))
)
print(complexity_summary)

# Boundary area should be nearly identical before/after simplification.
area_comparison <- tibble(
  object   = c("Raw boundary", "Simplified boundary"),
  area_km2 = c(as.numeric(st_area(boundary_raw)), as.numeric(st_area(boundary)))
) %>%
  mutate(pct_change = 100 * (area_km2 / first(area_km2) - 1))
print(area_comparison)

# ------------------------------------------------------------
# Save study-area object
# ------------------------------------------------------------

study_area <- list(
  crs                 = crs_aea_km,
  boundary            = boundary,
  boundary_buffer_5km = boundary_buffer_5km,
  bcr                 = bcr
)

out_dir <- file.path(paths$data_clean, "spatial")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

saveRDS(study_area, file = file.path(out_dir, "study_area.rds"))

message("01_define_study_area.R complete.")
