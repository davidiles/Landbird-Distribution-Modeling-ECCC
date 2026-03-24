# ============================================================
# 10_make_species_figures.R
#
# Purpose:
#   Generate PNG figures per species using the per-species
#   prediction summary objects created in script 08.
#
# Inputs:
#   data_clean/birds/data_ready_for_analysis.rds
#   data_clean/model_output/predictions/<sp>.rds  (many)
#   + external shapefiles for BCR, water, (optional) atlas squares
#
# Outputs:
#   data_clean/model_output/figures/<sp>_*.png
# ============================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(units)
  library(igraph)
  library(here)
})

source(here::here("R", "00_config_paths.R"))

# helper-functions
inla_utils_path <- file.path(paths$functions, "inla_model_utils.R")
figure_utils_path <- file.path(paths$functions, "figure_utils.R")

source(inla_utils_path)
source(figure_utils_path)

# ------------------------------------------------------------
# Config
# ------------------------------------------------------------

model_type <- "PC_ARU"

# File paths
in_data  <- file.path(paths$data_clean, "birds", "data_ready_for_analysis.rds")
pred_dir <- file.path(paths$model_output, paste0("predictions_", model_type))
fig_dir  <- file.path(paths$model_output, paste0("figures_", model_type))
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

in_water <- file.path(paths$data_clean, "spatial", "water_filtered.shp")
in_atlas_squares <- file.path(paths$data, "Spatial", "National", "AtlasSquares", "NationalSquares_FINAL.shp")

# Plot export settings
dpi <- 1000
width_in <- 10
height_in <- 8
ggsave_type <- "cairo"

# Raster resolution in map units (meters in EPSG:3978)
plot_res <- 1001

# set TRUE to plot square-level detection summaries (0/1 binary)
make_atlas_square_overlay <- TRUE

# ------------------------------------------------------------
# Load base data
# ------------------------------------------------------------

stopifnot(file.exists(in_data))
dat <- readRDS(in_data)

# raw data
all_surveys <- dat$all_surveys
counts      <- dat$counts

# prediction grid / study area
grid2 <- dat$grid_OBBA2
grid3 <- dat$grid_OBBA3
study_boundary <- dat$study_boundary %>% sf::st_as_sf()

hex_grid <- dat$hex_grid

stopifnot(inherits(grid2, "sf"), inherits(grid3, "sf"))

# Ensure consistent CRS
crs_use <- st_crs(grid2)
study_boundary <- st_transform(study_boundary, crs_use)
grid3 <- st_transform(grid3, crs_use)
grid2 <- st_transform(grid2, crs_use)

# BCR outlines
bcr_sf <- dat$bcr_sf

# Water
water_sf <- NULL
if (file.exists(in_water)) {
  water_sf <- st_read(in_water, quiet = TRUE) %>%
    st_make_valid() %>%
    st_transform(crs_use)
}

# Optional atlas squares overlay (centroids only) — we will filter to *surveyed* squares per atlas period
atlas_sq_centroids_all <- NULL
if (make_atlas_square_overlay && file.exists(in_atlas_squares)) {
  atlas_sq <- st_read(in_atlas_squares, quiet = TRUE) %>%
    st_transform(crs_use)
  
  # If you have a prov column in that file (as in your old script)
  if ("prov" %in% names(atlas_sq)) atlas_sq <- atlas_sq %>% filter(prov == "ON")
  
  atlas_sq_centroids_all <- st_centroid(atlas_sq)
}

# ------------------------------------------------------------
# Find prediction files
# ------------------------------------------------------------

pred_files <- list.files(pred_dir, pattern = "\\.rds$", full.names = TRUE)
stopifnot(length(pred_files) > 0)

message("Prediction files found: ", length(pred_files))

# ------------------------------------------------------------
# Loop species
# ------------------------------------------------------------

for (i in seq_len(length(pred_files))) {
  
  preds <- readRDS(pred_files[i])
  
  # ----------------------------------------------------------
  # Basic identifiers
  # ----------------------------------------------------------
  sp_english <- preds$sp_name
  sp_file <- sp_english |>
    stringr::str_to_lower() |>
    stringr::str_replace_all("[^a-z0-9]+", "_") |>
    stringr::str_replace_all("^_|_$", "")
  
  sp_code <- dat$species_to_model |>
    dplyr::filter(english_name == sp_english) |>
    dplyr::pull(species_id)
  
  fig_path_relabund_q50_OBBA2 <- file.path(fig_dir, paste0(sp_file, "_relabund_q50_OBBA2.png"))
  fig_path_relabund_q50_OBBA3 <- file.path(fig_dir, paste0(sp_file, "_relabund_q50_OBBA3.png"))
  fig_path_abs_change_q50     <- file.path(fig_dir, paste0(sp_file, "_abs_change_q50.png"))
  
  if (file.exists(fig_path_abs_change_q50)) {
    message("Skipping ", sp_english, ": change map already exists for this species")
    next
  }
  
  # ----------------------------------------------------------
  # Safety checks
  # ----------------------------------------------------------
  required_parts <- c("OBBA2", "OBBA3", "abs_change")
  if (!all(required_parts %in% names(preds)) ||
      any(vapply(preds[required_parts], is.null, logical(1)))) {
    message("Skipping (missing OBBA2/OBBA3/abs_change summaries): ", basename(pf))
    next
  }
  
  message("Mapping: ", sp_english)
  
  # ----------------------------------------------------------
  # Shared relative-abundance/change plotting scale
  # ----------------------------------------------------------
  zmax <- compute_shared_zmax(preds)
  rel_bounds <- list(lower = 0, upper = zmax)
  
  # ----------------------------------------------------------
  # Atlas square overlays
  # ----------------------------------------------------------
  atlas_sq_pts <- build_atlas_square_overlays(
    sp_square_summary = preds$sp_square_summary,
    atlas_sq_centroids_all = atlas_sq_centroids_all
  )
  
  # ----------------------------------------------------------
  # Atlas 2 / Atlas 3 relative abundance maps
  # ----------------------------------------------------------
  maps2 <- make_relabund_maps(
    species_name = sp_english,
    grid_sf = grid2,
    pred_summary = preds$OBBA2,
    prefix = "OBBA2",
    study_boundary = study_boundary,
    bcr_sf = bcr_sf,
    water_sf = water_sf,
    atlas_squares_centroids = atlas_sq_pts$OBBA2,
    title = "Atlas 2",
    subtitle = "Relative abundance",
    res = plot_res,
    bounds = rel_bounds
  )
  
  ggsave(
    filename = fig_path_relabund_q50_OBBA2,
    plot = maps2$q50_plot, width = width_in, height = height_in, units = "in",
    dpi = dpi, type = ggsave_type, limitsize = FALSE
  )
  
  maps3 <- make_relabund_maps(
    species_name = sp_english,
    grid_sf = grid3,
    pred_summary = preds$OBBA3,
    prefix = "OBBA3",
    study_boundary = study_boundary,
    bcr_sf = bcr_sf,
    water_sf = water_sf,
    atlas_squares_centroids = atlas_sq_pts$OBBA3,
    title = "Atlas 3",
    subtitle = "Relative abundance",
    res = plot_res,
    bounds = rel_bounds
  )
  
  ggsave(
    filename = fig_path_relabund_q50_OBBA3,
    plot = maps3$q50_plot, width = width_in, height = height_in, units = "in",
    dpi = dpi, type = ggsave_type, limitsize = FALSE
  )
  
  # ----------------------------------------------------------
  # Meaningful change polygons
  # ----------------------------------------------------------
  
  signif_change_polys <- build_meaningful_change_polys(
    eta_draws_per_hex = preds$eta_draws_per_hex,
    hexagon_sf = hex_grid,
    study_boundary = study_boundary,
    min_area_km2 = 10000,  # Area of contiguous polygons with "significant" change to be plotted
    smoothing_bandwidth_m = 75000,
    param = "abs_change",
    threshold = zmax / 20,
    prob_level = 0.90,
    direction = c("two_sided", "increase", "decrease"),
    ci_probs = c(0.05, 0.95),
    include_summary = TRUE,
    drop_holes = TRUE
  )
  #ggplot()+geom_sf(data = study_boundary)+geom_sf(data = signif_change_polys, fill = "black")
  
  chg <- make_abs_change_maps(
    species_name = sp_english,
    grid_sf = grid3,
    abs_change_summary = preds$abs_change,
    study_boundary = study_boundary,
    bcr_sf = bcr_sf,
    water_sf = water_sf,
    signif_poly = signif_change_polys,
    res = plot_res,
    max_abs = zmax,
    title = "Absolute change",
    subtitle = "OBBA2 to OBBA3"
  )
  
  ggsave(
    filename = fig_path_abs_change_q50,
    plot = chg$chg_plot, width = width_in, height = height_in, units = "in",
    dpi = dpi, type = ggsave_type, limitsize = FALSE
  )
  
  
  hex_change_sf <- summarize_abs_change_by_hex(
    eta_draws_per_hex = preds$eta_draws_per_hex,
    hexagon_sf = hex_grid,
    threshold = zmax / 20,
    ci_probs = c(0.05, 0.95)
  ) %>%
    st_centroid()
  
  ggplot(hex_change_sf) +
    geom_sf(data = study_boundary, fill = NA, colour = "black")+
    geom_sf(aes(col = mean_abs_change))
  
}

message("10_make_species_figures.R complete: ", fig_dir)