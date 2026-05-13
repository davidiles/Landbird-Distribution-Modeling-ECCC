# ============================================================
# 08_generate_model_products.R
#
# Purpose
#   Convert fitted model prediction RDS files into map-ready rasters and
#   multi-page species assessment PDFs.
#
# Main outputs
#   - model_output/rasters_<model_name>/<species>_a2.tif
#   - model_output/rasters_<model_name>/<species>_a3.tif
#   - model_output/rasters_<model_name>/<species>_chg.tif
#   - model_output/figures_<model_name>/<species>.pdf
#
# Notes
#   - This script preserves the prediction objects produced upstream.
#   - Relative abundance means are water-corrected; CI widths are currently
#     calculated from the uncorrected summaries, matching the original script.
#   - CRS choices are inherited from the saved analysis-ready data and are not
#     modified beyond aligning objects to the prediction-grid CRS.
# ============================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(ggplot2)
  library(terra)
  library(patchwork)
  library(magick)
  library(viridis)
})

source(here::here("R", "00_config_paths.R"))
source(file.path(paths$functions, "model_product_utils.R"))

# ------------------------------------------------------------
# Configuration and paths
# ------------------------------------------------------------

model_name <- "m2"
plot_res <- 1001  # raster resolution in map units; meters in EPSG:3978

in_data  <- file.path(paths$data_clean, "birds", "data_ready_for_analysis.rds")
pred_dir <- file.path(paths$model_output, paste0("predictions_", model_name))
rast_dir <- file.path(paths$model_output, paste0("rasters_", model_name))
fig_dir  <- file.path(paths$model_output, paste0("figures_", model_name))
png_dir  <- file.path(fig_dir, "PNGs")
out_dir  <- paths$model_output

model_summaries_path <- file.path(
  out_dir,
  paste0("summaries_", model_name),
  "model_summaries.rds"
)

for (dir_i in c(rast_dir, fig_dir, png_dir)) {
  dir.create(dir_i, recursive = TRUE, showWarnings = FALSE)
}

stopifnot(file.exists(in_data))
stopifnot(file.exists(model_summaries_path))

# ------------------------------------------------------------
# Load analysis-ready data and model summaries
# ------------------------------------------------------------

dat <- readRDS(in_data)
model_summaries <- readRDS(model_summaries_path)

grid2 <- dat$grid_OBBA2 %>% na.omit()
grid3 <- dat$grid_OBBA3 %>% na.omit()
study_boundary <- dat$study_boundary %>% sf::st_as_sf()
hex_grid <- dat$hex_grid_25km

# Align all map layers to the Atlas 2 prediction-grid CRS.
crs_use <- sf::st_crs(grid2)
study_boundary <- sf::st_transform(study_boundary, crs_use)
grid2 <- sf::st_transform(grid2, crs_use)
grid3 <- sf::st_transform(grid3, crs_use)

# ------------------------------------------------------------
# Load or create water polygons used as map overlays
# ------------------------------------------------------------

in_water <- file.path(paths$data_clean, "spatial", "water_to_plot.rds")

if (file.exists(in_water)) {
  water_to_plot <- readRDS(in_water)
} else {
  in_water_raw <- file.path(
    paths$data,
    "Spatial",
    "Ontario_Hydro_Network_(OHN)_-_Waterbody",
    "Ontario_Hydro_Network_(OHN)_-_Waterbody.shp"
  )
  
  water_to_plot <- sf::st_read(in_water_raw, quiet = TRUE) %>%
    sf::st_make_valid() %>%
    sf::st_transform(crs_use) %>%
    dplyr::mutate(area_m2 = as.numeric(sf::st_area(geometry))) %>%
    dplyr::filter(VERIFICATI == "Verified", area_m2 >= (5000 * 5000)) %>%
    sf::st_intersection(study_boundary)
  
  saveRDS(water_to_plot, in_water)
}

# ------------------------------------------------------------
# Generate raster and figure products for each species
# ------------------------------------------------------------

pred_files <- list.files(pred_dir, pattern = "\\.rds$", full.names = TRUE)
stopifnot(length(pred_files) > 0)

message("Prediction files found: ", length(pred_files))

for (i in seq_along(pred_files)) {
  
  preds <- readRDS(pred_files[i])
  
  sp_english <- preds$sp_name
  sp_file <- sp_english |>
    stringr::str_to_lower() |>
    stringr::str_replace_all("[^a-z0-9]+", "_") |>
    stringr::str_replace_all("^_|_$", "")
  
  sp_code <- dat$species_to_model |>
    dplyr::filter(english_name == sp_english) |>
    dplyr::pull(species_id)
  
  if (length(sp_code) == 0) {
    warning("No species_id found for ", sp_english, "; using NA in raster metadata.")
    sp_code <- NA_character_
  }
  
  sp_model_summary <- model_summaries[[sp_english]]
  if (is.null(sp_model_summary)) {
    warning("No model summary found for ", sp_english, ".")
  }
  
  message("\nGenerating products for: ", sp_english)
  
  # Species-specific output paths.
  rast_path_a2  <- file.path(rast_dir, paste0(sp_file, "_a2.tif"))
  rast_path_a3  <- file.path(rast_dir, paste0(sp_file, "_a3.tif"))
  rast_path_chg <- file.path(rast_dir, paste0(sp_file, "_chg.tif"))
  dat_path <- file.path(
    out_dir,
    paste0("data_used_", model_name),
    paste0(sp_file, "_1km.rds")
  )
  
  stopifnot(file.exists(dat_path))
  sp_dat <- readRDS(dat_path)
  
  # ----------------------------------------------------------
  # Raster products: Atlas 2, Atlas 3, and absolute change
  # ----------------------------------------------------------
  
  a2 <- grid2 %>%
    dplyr::mutate(
      mu_mean = preds$OBBA2_Corrected_for_Water$OBBA2_mean,
      CI_95_width = preds$OBBA2_Corrected_for_Water$OBBA2_upper - preds$OBBA2_Corrected_for_Water$OBBA2_lower
    )
  
  r2 <- rasterize_sf(
    grid_sf = a2,
    field = c("mu_mean", "CI_95_width"),
    res = plot_res,
    metadata = c(
      species_name = sp_english,
      species_id = sp_code,
      species_filename = sp_file,
      model_name = model_name,
      units = "Expected count (Atlas 2)"
    )
  )
  
  terra::writeRaster(r2, filename = rast_path_a2, overwrite = TRUE)
  
  a3 <- grid3 %>%
    dplyr::mutate(
      mu_mean = preds$OBBA3_Corrected_for_Water$OBBA3_mean,
      CI_95_width = preds$OBBA3_Corrected_for_Water$OBBA3_upper - preds$OBBA3_Corrected_for_Water$OBBA3_lower
    )
  
  r3 <- rasterize_sf(
    grid_sf = a3,
    field = c("mu_mean", "CI_95_width"),
    res = plot_res,
    metadata = c(
      species_name = sp_english,
      species_id = sp_code,
      species_filename = sp_file,
      model_name = model_name,
      units = "Expected count (Atlas 3)"
    )
  )
  
  terra::writeRaster(r3, filename = rast_path_a3, overwrite = TRUE)
  
  chg_sf <- grid3 %>%
    dplyr::mutate(
      chg_mean = preds$abs_change_Corrected_for_Water$abs_change_mean,
      CI_95_width = preds$abs_change$abs_change_upper - preds$abs_change$abs_change_lower
    )
  
  rchg <- rasterize_sf(
    grid_sf = chg_sf,
    field = c("chg_mean", "CI_95_width"),
    res = plot_res,
    metadata = c(
      species_name = sp_english,
      species_id = sp_code,
      species_filename = sp_file,
      model_name = model_name,
      units = "Change in expected counts (Atlas 2 to Atlas 3)"
    )
  )
  
  terra::writeRaster(rchg, filename = rast_path_chg, overwrite = TRUE)
  
  # ----------------------------------------------------------
  # Relative abundance maps
  # ----------------------------------------------------------
  
  # Assume the species is effectively absent if model predicts expected count less than 1 per 150 point counts
  relabund_absent_limit <- 1 / 150
  
  rasters_relabund_prepared <- prepare_relative_abundance_rasters(
    Atlas2 = r2,
    Atlas3 = r3,
    rast_absent_limit = relabund_absent_limit,
    rast_max_quantile = 0.99
  )
  
  
  colpal_relabund <- grDevices::colorRampPalette(c(
    "#FBF7E2", "#FCF8D0", "#EEF7C2", "#CEF2B0",
    "#94E5A0", "#51C987", "#18A065", "#008C59",
    "#007F53", "#006344"
  ))(100)
  
  
  Relabund_Map_Atlas2 <- make_map(
    species_name = sp_english,
    subtitle = "Relative Abundance - Atlas 2",
    legend_title = "Expected count\nper 5-min\npoint count",
    rast = rasters_relabund_prepared$rasters$Atlas2,
    region = study_boundary,
    water = water_to_plot,
    colpal = colpal_relabund,
    water_fill = "white",
    transform = "sqrt",
    zlim = rasters_relabund_prepared$zlim,
    zbreaks = rasters_relabund_prepared$zbreaks
  )
  
  Relabund_Map_Atlas3 <- make_map(
    species_name = sp_english,
    subtitle = "Relative Abundance - Atlas 3",
    legend_title = "Expected count\nper 5-min\npoint count",
    rast = rasters_relabund_prepared$rasters$Atlas3,
    region = study_boundary,
    water = water_to_plot,
    colpal = colpal_relabund,
    water_fill = "white",
    transform = "sqrt",
    zlim = rasters_relabund_prepared$zlim,
    zbreaks = rasters_relabund_prepared$zbreaks
  )
  
  # ----------------------------------------------------------
  # Hex-level change map
  # ----------------------------------------------------------
  
  prov_change <- summarize_polygon_hex_draw_change(
    hex_draws = preds$hex_draws_Corrected_for_Water,
    hex_grid = hex_grid,
    polygon = study_boundary,
    ci_level = 0.90
  ) |>
    dplyr::mutate(
      pct_change_median = 100 * prop_change_median,
      pct_change_qlow = 100 * prop_change_qlow,
      pct_change_qhigh = 100 * prop_change_qhigh
    )
  
  hex_change_sf <- summarize_hex_draw_change(
    hex_grid = hex_grid,
    hex_draws = preds$hex_draws_Corrected_for_Water,
    ci_level = 0.90
  ) |>
    classify_min_supported_change()
  
  Hex_Change_Map <- make_hex_abs_change_map(
    species_name = sp_english,
    hex_change_sf = hex_change_sf,
    region = study_boundary,
    prov_change = prov_change,
    ci_level = 0.90,
    zlim = rasters_relabund_prepared$zlim,
    water = water_to_plot,
    water_fill = "gray97"
  )
  
  # ----------------------------------------------------------
  # PObs Maps
  # ----------------------------------------------------------
  
  pobs2 <- 1-exp(-r2)
  pobs3 <- 1-exp(-r2)
  
  pobs_absent_limit <- 1/250
  
  rasters_pobs_prepared <- prepare_relative_abundance_rasters(
    Atlas2 = pobs2,
    Atlas3 = pobs3,
    rast_absent_limit = pobs_absent_limit,
    rast_max_quantile = 0.99
  )
  
  colpal_pobs <- viridis::mako(100, direction = -1) 
  
  pobs_Map_Atlas2 <- make_map(
    species_name = sp_english,
    subtitle = "Probability of Observation - Atlas 2",
    legend_title = "P(Obs)\nper 5-min\npoint count",
    rast = rasters_pobs_prepared$rasters$Atlas2,
    region = study_boundary,
    water = water_to_plot,
    colpal = colpal_pobs,
    water_fill = "white",
    transform = "identity",
    zlim = c(0,1),
    zbreaks = rasters_pobs_prepared$zbreaks
  )
  
  pobs_Map_Atlas3 <- make_map(
    species_name = sp_english,
    subtitle = "Probability of Observation - Atlas 3",
    legend_title = "P(Obs)\nper 5-min\npoint count",
    rast = rasters_pobs_prepared$rasters$Atlas3,
    region = study_boundary,
    water = water_to_plot,
    colpal = colpal_pobs,
    water_fill = "white",
    transform = "identity",
    zlim = c(0,1),
    zbreaks = rasters_pobs_prepared$zbreaks
  )
  
  # ----------------------------------------------------------
  # Model-assessment panels for each atlas period
  # ----------------------------------------------------------
  
  r2_clamped <- r2$mu_mean
  r2_clamped[r2_clamped < relabund_absent_limit] <- 0
  
  A2_dat_to_plot <- sp_dat$sp_dat %>%
    dplyr::filter(Survey_Type %in% c("Point_Count", "ARU"), Atlas == "OBBA2") %>%
    dplyr::mutate(count_per_effort = count)
  
  A2_assessment <- assess_region(
    region = study_boundary,
    sp_dat = A2_dat_to_plot,
    rast = r2_clamped,
    n_hexagons = 2000,
    water = NULL,
    rast_max_q = 0.99,
    transform = "identity",
    title = NULL,
    model_source = "Atlas 2 Predictions",
    data_source = "Atlas 2 Point Counts + ARUs"
  )
  
  r3_clamped <- r3$mu_mean
  r3_clamped[r3_clamped < relabund_absent_limit] <- 0
  
  A3_dat_to_plot <- sp_dat$sp_dat %>%
    dplyr::filter(Survey_Type %in% c("Point_Count", "ARU"), Atlas == "OBBA3") %>%
    dplyr::mutate(count_per_effort = count)
  
  A3_assessment <- assess_region(
    region = study_boundary,
    sp_dat = A3_dat_to_plot,
    rast = r3_clamped,
    n_hexagons = 2000,
    water = NULL,
    rast_max_q = 0.99,
    transform = "identity",
    title = NULL,
    model_source = "Atlas 3 Predictions",
    data_source = "Atlas 3 Point Counts + ARUs"
  )
  
  # ----------------------------------------------------------
  # Save one multi-page PDF per species
  # ----------------------------------------------------------
  
  page1_png <- file.path(png_dir, paste0(sp_file, "_page1.png"))
  page2_png <- file.path(png_dir, paste0(sp_file, "_page2.png"))
  page3_png <- file.path(png_dir, paste0(sp_file, "_page3.png"))
  page4_png <- file.path(png_dir, paste0(sp_file, "_page4.png"))
  pdf_path  <- file.path(fig_dir, paste0(sp_file, ".pdf"))
  
  # --- Page 1: relative abundance and change maps
  page1 <-
    Relabund_Map_Atlas2 +
    Relabund_Map_Atlas3 +
    Hex_Change_Map +
    patchwork::plot_layout(ncol = 3, widths = c(1, 1, 1))
  
  ragg::agg_png(
    filename = page1_png,
    width = 24,
    height = 8,
    units = "in",
    res = 300,
    background = "white"
  )
  print(page1)
  grDevices::dev.off()
  
  # --- Page 2: PObs maps
  page2 <-
    pobs_Map_Atlas2 +
    pobs_Map_Atlas3 +
    patchwork::plot_layout(ncol = 2, widths = c(1, 1))
  
  ragg::agg_png(
    filename = page2_png,
    width = 24,
    height = 8,
    units = "in",
    res = 300,
    background = "white"
  )
  print(page2)
  grDevices::dev.off()
  
  # --- Page 3: Model assessment figures for Atlas 2
  ragg::agg_png(
    filename = page3_png,
    width = 24,
    height = 8,
    units = "in",
    res = 300,
    background = "white"
  )
  print(A2_assessment$plot_combined)
  grDevices::dev.off()
  
  # --- Page 4: Model assessment figures for Atlas 3
  ragg::agg_png(
    filename = page4_png,
    width = 24,
    height = 8,
    units = "in",
    res = 300,
    background = "white"
  )
  print(A3_assessment$plot_combined)
  grDevices::dev.off()
  
  imgs <- magick::image_read(c(page1_png, page2_png, page3_png, page4_png))
  magick::image_write(imgs, path = pdf_path, format = "pdf")
}

message("08_generate_model_products.R complete")
