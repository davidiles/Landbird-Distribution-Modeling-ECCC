# ============================================================
# 09_make_species_figures.R
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

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(ggplot2)
})

# ---- utils
source("R/functions/figure_utils.R")

# ------------------------------------------------------------
# Config
# ------------------------------------------------------------

in_data <- "data_clean/birds/data_ready_for_analysis.rds"
pred_dir <- "data_clean/model_output/predictions"
fig_dir  <- "data_clean/model_output/figures"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# Shapefile inputs
in_bcr <- "../../Data/Spatial/BCR/BCR_Terrestrial_master.shp"
in_water <- "data_clean/spatial/water_filtered.shp"
in_atlas_squares <- "../../Data/Spatial/National/AtlasSquares/NationalSquares_FINAL.shp"  # optional

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

stopifnot(inherits(grid2, "sf"), inherits(grid3, "sf"))

# Ensure consistent CRS
crs_use <- st_crs(grid2)
study_boundary <- st_transform(study_boundary, crs_use)
grid3 <- st_transform(grid3, crs_use)
grid2 <- st_transform(grid2, crs_use)

# BCR outlines
bcr_sf <- st_read(in_bcr, quiet = TRUE) %>%
  st_make_valid() %>%
  st_transform(crs_use) %>%
  dplyr::filter(PROVINCE_S %in% c("ONTARIO", "ON", "Ontario") | is.na(PROVINCE_S)) %>%  # defensive
  dplyr::select(BCR, BCRNAME, PROVINCE_S) %>%
  group_by(BCR, BCRNAME) %>%
  summarise(geometry = st_union(geometry), .groups = "drop")

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

for (pf in pred_files) {
  
  preds <- readRDS(pf)
  sp_english <- preds$sp_english
  sp_file <- sp_english %>%
    str_to_lower() %>%
    str_replace_all("[^a-z0-9]+", "_") %>%
    str_replace_all("^_|_$", "")
  
  sp_code <- subset(dat$species_to_model, english_name == sp_english)$species_id
  
  sp_dat <- all_surveys %>%
    mutate(
      count = counts[[sp_code]],
      days_since_june15 = DayOfYear - 166,
      BCR_factor = as.numeric(factor(BCR)),
      Atlas3 = ifelse(Atlas == "OBBA2", 0, 1)
    )
  
  message("Mapping: ", sp_english)
  
  # ----------------------------------------------------------
  # Shared colour-scale limit across OBBA2 + OBBA3 + change map
  # zmax = max(99th percentile of q50 surface in OBBA2, OBBA3)
  # ----------------------------------------------------------
  zmax2 <- as.numeric(stats::quantile(preds$OBBA2$OBBA2_q50, 0.99, na.rm = TRUE))
  zmax3 <- as.numeric(stats::quantile(preds$OBBA3$OBBA3_q50, 0.99, na.rm = TRUE))
  zmax <- max(zmax2, zmax3, na.rm = TRUE)
  if (!is.finite(zmax) || zmax <= 0) zmax <- 1
  rel_bounds <- list(lower = 0, upper = zmax)
  
  # ----------------------------------------------------------
  # Optional atlas square overlays (surveyed squares only;
  # black = detected, gray = not detected)
  # ----------------------------------------------------------
  atlas_sq_pts2 <- NULL
  atlas_sq_pts3 <- NULL
  if (make_atlas_square_overlay && !is.null(atlas_sq_centroids_all)) {
    if (!("square_id" %in% names(sp_dat))) {
      stop("sp_dat/all_surveys must contain a 'square_id' column to build atlas-square overlays.")
    }
    
    sq2 <- sp_dat %>%
      sf::st_drop_geometry() %>%
      filter(Atlas == "OBBA2") %>%
      group_by(square_id) %>%
      summarise(detected = any(count > 0, na.rm = TRUE), .groups = "drop")
    
    sq3 <- sp_dat %>%
      sf::st_drop_geometry() %>%
      filter(Atlas == "OBBA3") %>%
      group_by(square_id) %>%
      summarise(detected = any(count > 0, na.rm = TRUE), .groups = "drop")
    
    # Attribute join by square_id (no spatial join)
    atlas_sq_pts2 <- atlas_sq_centroids_all %>%
      left_join(sq2, by = "square_id") %>%
      filter(!is.na(detected))
    
    atlas_sq_pts3 <- atlas_sq_centroids_all %>%
      left_join(sq3, by = "square_id") %>%
      filter(!is.na(detected))
  }
  
  # Safety checks
  if (is.null(preds$OBBA2) || is.null(preds$OBBA3) || is.null(preds$abs_change)) {
    message("  Skipping (missing OBBA2/OBBA3/abs_change summaries): ", basename(pf))
    next
  }
  
  # --- Relabund maps
  maps2 <- make_relabund_maps(
    species_name = sp_english,
    grid_sf = grid2,
    pred_summary = preds$OBBA2,
    prefix = "OBBA2",
    study_boundary = study_boundary,
    bcr_sf = bcr_sf,
    water_sf = water_sf,
    atlas_squares_centroids = atlas_sq_pts2,
    title = "Atlas 2",
    subtitle = "Relative Abundance",
    res = plot_res,
    bounds = rel_bounds
  )
  
  maps3 <- make_relabund_maps(
    species_name = sp_english,
    grid_sf = grid3,
    pred_summary = preds$OBBA3,
    prefix = "OBBA3",
    study_boundary = study_boundary,
    bcr_sf = bcr_sf,
    water_sf = water_sf,
    atlas_squares_centroids = atlas_sq_pts3,
    title = "Atlas 3",
    subtitle = "Relative Abundance",
    res = plot_res,
    bounds = rel_bounds
  )
  
  # --- Change maps (absolute change + CI width)
  chg <- make_abs_change_maps(
    species_name = sp_english,
    grid_sf = grid3,
    abs_change_summary = preds$abs_change,
    study_boundary = study_boundary,
    bcr_sf = bcr_sf,
    water_sf = water_sf,
    res = plot_res,
    max_abs = zmax,
    title = "Absolute change",
    subtitle = "OBBA2 to OBBA3"
  )
  
  # --- Save
  
  ggsave(
    filename = file.path(fig_dir, paste0(sp_file, "_relabund_q50_OBBA2.png")),
    plot = maps2$q50_plot, width = width_in, height = height_in, units = "in",
    dpi = dpi, type = ggsave_type, limitsize = FALSE
  )
  
  # ggsave(
  #   filename = file.path(fig_dir, paste0(sp_file, "_relabund_CV_OBBA2.png")),
  #   plot = maps2$cv_plot, width = width_in, height = height_in, units = "in",
  #   dpi = dpi, type = ggsave_type, limitsize = FALSE
  # )
  
  ggsave(
    filename = file.path(fig_dir, paste0(sp_file, "_relabund_q50_OBBA3.png")),
    plot = maps3$q50_plot, width = width_in, height = height_in, units = "in",
    dpi = dpi, type = ggsave_type, limitsize = FALSE
  )
  
  # ggsave(
  #   filename = file.path(fig_dir, paste0(sp_file, "_relabund_CV_OBBA3.png")),
  #   plot = maps3$cv_plot, width = width_in, height = height_in, units = "in",
  #   dpi = dpi, type = ggsave_type, limitsize = FALSE
  # )
  
  ggsave(
    filename = file.path(fig_dir, paste0(sp_file, "_abs_change_q50.png")),
    plot = chg$chg_plot, width = width_in, height = height_in, units = "in",
    dpi = dpi, type = ggsave_type, limitsize = FALSE
  )
  
  # ggsave(
  #   filename = file.path(fig_dir, paste0(sp_file, "_abs_change_ciw.png")),
  #   plot = chg$ciw_plot, width = width_in, height = height_in, units = "in",
  #   dpi = dpi, type = ggsave_type, limitsize = FALSE
  # )
  
}

message("09_make_species_figures.R complete: ", fig_dir)
