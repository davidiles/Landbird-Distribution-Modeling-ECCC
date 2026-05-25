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
source(file.path(paths$functions, "model_product_utils_updated.R"))

# ------------------------------------------------------------
# Configuration and paths
# ------------------------------------------------------------

model_name <- "PC_ARU_only"
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
    sf::st_intersection(study_boundary)|>
    sf::st_collection_extract("POLYGON")
  
  saveRDS(water_to_plot, in_water)
}

# ------------------------------------------------------------
# Load atlas regions for summarizing relative abundance and change estimates
# ------------------------------------------------------------

atlas_regions <- st_read(file.path(paths$data, "Spatial", "Ontario_Atlas_biol_regions", "Atlas_biol_regions.shp")) %>%
  st_transform(st_crs(study_boundary)) %>%
  arrange(Biol_Regio)

# ------------------------------------------------------------
# Load estimates of change from paired analysis
# ------------------------------------------------------------

paired_analysis_path <- file.path(
  out_dir,
  paste0("summaries_", model_name),
  "paired_summaries.rds"
)

if (file.exists(paired_analysis_path)){
  paired_summaries <- readRDS(paired_analysis_path)
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
  
  # Estimates of population change based on paired surveys
  regional_estimates_paired <- paired_summaries[[sp_english]]
  if (is.null(regional_estimates_paired)){
    message("Missing paired survey change analysis for ", sp_english,"; skipping")
    next
  }
  
  pdf_path  <- file.path(fig_dir, paste0(sp_file, ".pdf"))
  # if (file.exists(pdf_path)){
  #   message("Outputs already prepared for ", sp_english,"; skipping")
  #   next
  # }
  
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
      mu_q50 = preds$OBBA2_Corrected_for_Water$OBBA2_q50,
      CI_95_width = preds$OBBA2_Corrected_for_Water$OBBA2_upper - preds$OBBA2_Corrected_for_Water$OBBA2_lower
    )
  
  r2 <- rasterize_sf(
    grid_sf = a2,
    field = c("mu_q50", "CI_95_width"),
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
      mu_q50 = preds$OBBA3_Corrected_for_Water$OBBA3_q50,
      CI_95_width = preds$OBBA3_Corrected_for_Water$OBBA3_upper - preds$OBBA3_Corrected_for_Water$OBBA3_lower
    )
  
  r3 <- rasterize_sf(
    grid_sf = a3,
    field = c("mu_q50", "CI_95_width"),
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
      chg_q50 = preds$abs_change_Corrected_for_Water$abs_change_q50,
      CI_95_width = preds$abs_change$abs_change_upper - preds$abs_change$abs_change_lower
    )
  
  rchg <- rasterize_sf(
    grid_sf = chg_sf,
    field = c("chg_q50", "CI_95_width"),
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
  
  # Assume the species is effectively absent if model predicts expected count less than 1 per 200 point counts
  relabund_absent_limit <- 1 / 200
  
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
    water = NULL,
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
    water = NULL,
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
    region_boundaries = atlas_regions,
    prov_change = prov_change,
    ci_level = 0.90,
    zlim = rasters_relabund_prepared$zlim,
    water = NULL,
    water_fill = "gray97"
  )
  
  # ----------------------------------------------------------
  # Estimates of change within each atlas region
  # ----------------------------------------------------------
  
  # ---- Full spatial model
  regional_estimates_FullModel <- data.frame()
  
  for (j in 1:nrow(atlas_regions)){
    
    a_reg <- atlas_regions[j,]
    
    reg_est <- summarize_polygon_hex_draw_change(
      hex_draws = preds$hex_draws_Corrected_for_Water,
      hex_grid = hex_grid,
      polygon = a_reg,
      ci_level = 0.90
    ) |>
      dplyr::mutate(
        pct_change_median = 100 * prop_change_median,
        pct_change_qlow = 100 * prop_change_qlow,
        pct_change_qhigh = 100 * prop_change_qhigh
      ) %>%
      mutate(Region_Number = a_reg$Biol_Regio,
             Region_Name = a_reg$ECOZONE_NA) %>%
      relocate(Region_Name,Region_Number)
    
    regional_estimates_FullModel <- bind_rows(regional_estimates_FullModel,reg_est)
  }
  
  # ---- Paired analysis (only using surveys within shared 50 m footprint)
  paired_data_summary <- paired_summaries[[sp_english]]$shared_data_summary
  paired_change_summary <- paired_summaries[[sp_english]]$shared_change_summary
  paired_data <- paired_summaries[[sp_english]]$shared_data
  
  change_comparison_plot <- make_change_comparison_plot(
    regional_estimates_FullModel = regional_estimates_FullModel,
    paired_change_summary = paired_change_summary,
    title = "Regional population change estimates",
    subtitle = "Comparison of full atlas model vs repeated surveys only"
  )
  
  # Plot locations of shared data
  paired_locations_map <-
    ggplot2::ggplot() +
    ggplot2::geom_sf(data = study_boundary, colour = "black") +
    ggplot2::geom_sf(data = paired_data, size = 0.1, col = "#1FAA59") +
    ggplot2::geom_sf(
      data = atlas_regions %>% sf::st_intersection(study_boundary),
      fill = "transparent",
      colour = "black",
      linewidth = 0.3
    ) +
    ggplot2::ggtitle("Repeated survey locations") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(size = 13),
      plot.margin = ggplot2::margin(5, 5, 5, 5)
    )
  
  paired_data_summary <- paired_data %>%
    as.data.frame() %>%
    group_by(Biol_Region,Atlas) %>%
    summarize(mean_count = mean(count),
              PObs = mean(count>0),
              n_surveys = n(),
              n_squares = length(unique(square_id))) %>%
    left_join(as.data.frame(atlas_regions)[,c("Biol_Regio","ECOZONE_NA")], by = c("Biol_Region" = "Biol_Regio"))
  
  paired_summary_table <- make_paired_summary_table(paired_data_summary)
  
  
  # ----------------------------------------------------------
  # Provincial model assessment panels for each atlas period
  # ----------------------------------------------------------
  
  n_hexagons <- 2000
  
  r2_clamped <- r2$mu_q50
  r2_clamped[r2_clamped < relabund_absent_limit] <- 0
  
  r3_clamped <- r3$mu_q50
  r3_clamped[r3_clamped < relabund_absent_limit] <- 0
  
  A2_dat_to_plot <- sp_dat$sp_dat %>%
    dplyr::filter(Survey_Type %in% c("Point_Count", "ARU"), Atlas == "OBBA2") %>%
    dplyr::mutate(count_per_effort = count)
  
  A3_dat_to_plot <- sp_dat$sp_dat %>%
    dplyr::filter(Survey_Type %in% c("Point_Count", "ARU"), Atlas == "OBBA3") %>%
    dplyr::mutate(count_per_effort = count)
  
  # Use one shared assessment hex grid for both atlases
  assessment_hex_grid <- make_hex_grid(
    study_area = study_boundary,
    n_hexagons = n_hexagons
  )
  
  # First create hex summaries only, so shared display limits are based on
  # both atlas periods together.
  A2_hex_summary <- summarize_hex(
    dat = A2_dat_to_plot,
    hex_grid = assessment_hex_grid,
    rast = r2_clamped,
    rast_max_q = 0.99
  )
  
  A3_hex_summary <- summarize_hex(
    dat = A3_dat_to_plot,
    hex_grid = assessment_hex_grid,
    rast = r3_clamped,
    rast_max_q = 0.99
  )
  
  assessment_scale_limits <- compute_assessment_scale_limits(
    hex_summaries = list(A2_hex_summary, A3_hex_summary),
    rasts = list(r2_clamped, r3_clamped),
    rast_max_q = 0.99,
    count_max_q = 0.99,
    max_surveys_q = 0.80
  )
  
  A2_assessment <- assess_region(
    region = study_boundary,
    region_boundaries = atlas_regions,
    sp_dat = A2_dat_to_plot,
    rast = r2_clamped,
    hex_grid = assessment_hex_grid,
    water = NULL,
    rast_max = assessment_scale_limits$rast_max,
    max_count_per_effort = assessment_scale_limits$max_count_per_effort,
    max_surveys = assessment_scale_limits$max_surveys,
    transform = "sqrt",
    title = paste0(sp_english," - Atlas 2"),
    model_source = NULL,
    data_source = NULL
  )
  
  A3_assessment <- assess_region(
    region = study_boundary,
    region_boundaries = atlas_regions,
    sp_dat = A3_dat_to_plot,
    rast = r3_clamped,
    hex_grid = assessment_hex_grid,
    water = NULL,
    rast_max = assessment_scale_limits$rast_max,
    max_count_per_effort = assessment_scale_limits$max_count_per_effort,
    max_surveys = assessment_scale_limits$max_surveys,
    transform = "sqrt",
    title = paste0(sp_english," - Atlas 3"),
    model_source = NULL,
    data_source = NULL
  )
  
  # ----------------------------------------------------------
  # Southern Ontario assessment panels for each atlas period
  # ----------------------------------------------------------
  
  southern_ontario <- subset(atlas_regions, Biol_Regio %in% c(1,2,3)) %>%
    st_union() %>%
    st_intersection(study_boundary) %>% 
    st_transform(st_crs(A2_dat_to_plot)) %>%
    st_as_sf()
  
  southern_ontario_vect <- terra::vect(southern_ontario %>% st_transform(crs(r2_clamped)))
  
  r2_southern <- r2_clamped |>
    terra::crop(southern_ontario_vect) |>
    terra::mask(southern_ontario_vect)
  
  r3_southern <- r3_clamped |>
    terra::crop(southern_ontario_vect) |>
    terra::mask(southern_ontario_vect)
  
  A2_southern <- A2_dat_to_plot |>
    sf::st_filter(southern_ontario, .predicate = sf::st_intersects)
  
  A3_southern <- A3_dat_to_plot |>
    sf::st_filter(southern_ontario, .predicate = sf::st_intersects)
  
  # Use one shared assessment hex grid for both atlases
  southern_assessment_hex_grid <- make_hex_grid(
    study_area = southern_ontario,
    n_hexagons = n_hexagons
  )
  
  # First create hex summaries only, so shared display limits are based on
  # both atlas periods together.
  A2_southern_hex_summary <- summarize_hex(
    dat = A2_southern,
    hex_grid = southern_assessment_hex_grid,
    rast = r2_southern,
    rast_max_q = 0.99
  )
  
  A3_southern_hex_summary <- summarize_hex(
    dat = A3_southern,
    hex_grid = southern_assessment_hex_grid,
    rast = r3_southern,
    rast_max_q = 0.99
  )
  
  southern_assessment_scale_limits <- compute_assessment_scale_limits(
    hex_summaries = list(A2_southern_hex_summary, A3_southern_hex_summary),
    rasts = list(r2_southern, r3_southern),
    rast_max_q = 0.99,
    count_max_q = 0.99,
    max_surveys_q = 0.80
  )
  
  A2_southern_assessment <- assess_region(
    region = southern_ontario,
    region_boundaries = atlas_regions,
    sp_dat = A2_southern,
    rast = r2_southern,
    hex_grid = southern_assessment_hex_grid,
    water = NULL,
    rast_max = southern_assessment_scale_limits$rast_max,
    max_count_per_effort = southern_assessment_scale_limits$max_count_per_effort,
    max_surveys = southern_assessment_scale_limits$max_surveys,
    transform = "identity",
    title = paste0(sp_english," - Atlas 2 - Southern Ontario"),
    model_source = NULL,
    data_source = NULL
  )
  
  A3_southern_assessment <- assess_region(
    region = southern_ontario,
    region_boundaries = atlas_regions,
    sp_dat = A3_southern,
    rast = r3_southern,
    hex_grid = southern_assessment_hex_grid,
    water = NULL,
    rast_max = southern_assessment_scale_limits$rast_max,
    max_count_per_effort = southern_assessment_scale_limits$max_count_per_effort,
    max_surveys = southern_assessment_scale_limits$max_surveys,
    transform = "identity",
    title = paste0(sp_english," - Atlas 3 - Southern Ontario"),
    model_source = NULL,
    data_source = NULL
  )
  
  # ----------------------------------------------------------
  # Boreal assessment panels for each atlas period
  # ----------------------------------------------------------
  
  boreal_ontario <- subset(atlas_regions, Biol_Regio %in% c(4)) %>%
    st_union() %>%
    st_intersection(study_boundary) %>% 
    st_transform(st_crs(A2_dat_to_plot)) %>%
    st_as_sf()
  
  boreal_ontario_vect <- terra::vect(boreal_ontario %>%
                                       st_transform(crs(r2_clamped)))
  
  r2_boreal <- r2_clamped |>
    terra::crop(boreal_ontario_vect) |>
    terra::mask(boreal_ontario_vect)
  
  r3_boreal <- r3_clamped |>
    terra::crop(boreal_ontario_vect) |>
    terra::mask(boreal_ontario_vect)
  
  
  A2_boreal <- A2_dat_to_plot |>
    sf::st_filter(boreal_ontario, .predicate = sf::st_intersects)
  
  A3_boreal <- A3_dat_to_plot |>
    sf::st_filter(boreal_ontario, .predicate = sf::st_intersects)
  
  # Use one shared assessment hex grid for both atlases
  boreal_assessment_hex_grid <- make_hex_grid(
    study_area = boreal_ontario,
    n_hexagons = n_hexagons
  )
  
  # First create hex summaries only, so shared display limits are based on
  # both atlas periods together.
  A2_boreal_hex_summary <- summarize_hex(
    dat = A2_boreal,
    hex_grid = boreal_assessment_hex_grid,
    rast = r2_boreal,
    rast_max_q = 0.99
  )
  
  A3_boreal_hex_summary <- summarize_hex(
    dat = A3_boreal,
    hex_grid = boreal_assessment_hex_grid,
    rast = r3_boreal,
    rast_max_q = 0.99
  )
  
  boreal_assessment_scale_limits <- compute_assessment_scale_limits(
    hex_summaries = list(A2_boreal_hex_summary, A3_boreal_hex_summary),
    rasts = list(r2_boreal, r3_boreal),
    rast_max_q = 0.99,
    count_max_q = 0.99,
    max_surveys_q = 0.80
  )
  
  A2_boreal_assessment <- assess_region(
    region = boreal_ontario,
    region_boundaries = atlas_regions,
    sp_dat = A2_boreal,
    rast = r2_boreal,
    hex_grid = boreal_assessment_hex_grid,
    water = NULL,
    rast_max = boreal_assessment_scale_limits$rast_max,
    max_count_per_effort = boreal_assessment_scale_limits$max_count_per_effort,
    max_surveys = boreal_assessment_scale_limits$max_surveys,
    transform = "identity",
    title = paste0(sp_english," - Atlas 2 - Boreal Ontario"),
    model_source = NULL,
    data_source = NULL
  )
  
  A3_boreal_assessment <- assess_region(
    region = boreal_ontario,
    region_boundaries = atlas_regions,
    sp_dat = A3_boreal,
    rast = r3_boreal,
    hex_grid = boreal_assessment_hex_grid,
    water = NULL,
    rast_max = boreal_assessment_scale_limits$rast_max,
    max_count_per_effort = boreal_assessment_scale_limits$max_count_per_effort,
    max_surveys = boreal_assessment_scale_limits$max_surveys,
    transform = "identity",
    title = paste0(sp_english," - Atlas 3 - Boreal Ontario"),
    model_source = NULL,
    data_source = NULL
  )
  
  # ----------------------------------------------------------
  # HB assessment panels for each atlas period
  # ----------------------------------------------------------
  
  HB_ontario <- subset(atlas_regions, Biol_Regio %in% c(5)) %>%
    st_union() %>%
    st_intersection(study_boundary) %>% 
    st_transform(st_crs(A2_dat_to_plot)) %>%
    st_as_sf()
  
  HB_ontario_vect <- terra::vect(HB_ontario %>%
                                   st_transform(crs(r2_clamped)))
  
  r2_HB <- r2_clamped |>
    terra::crop(HB_ontario_vect) |>
    terra::mask(HB_ontario_vect)
  
  r3_HB <- r3_clamped |>
    terra::crop(HB_ontario_vect) |>
    terra::mask(HB_ontario_vect)
  
  
  A2_HB <- A2_dat_to_plot |>
    sf::st_filter(HB_ontario, .predicate = sf::st_intersects)
  
  A3_HB <- A3_dat_to_plot |>
    sf::st_filter(HB_ontario, .predicate = sf::st_intersects)
  
  # Use one shared assessment hex grid for both atlases
  HB_assessment_hex_grid <- make_hex_grid(
    study_area = HB_ontario,
    n_hexagons = n_hexagons
  )
  
  # First create hex summaries only, so shared display limits are based on
  # both atlas periods together.
  A2_HB_hex_summary <- summarize_hex(
    dat = A2_HB,
    hex_grid = HB_assessment_hex_grid,
    rast = r2_HB,
    rast_max_q = 0.99
  )
  
  A3_HB_hex_summary <- summarize_hex(
    dat = A3_HB,
    hex_grid = HB_assessment_hex_grid,
    rast = r3_HB,
    rast_max_q = 0.99
  )
  
  HB_assessment_scale_limits <- compute_assessment_scale_limits(
    hex_summaries = list(A2_HB_hex_summary, A3_HB_hex_summary),
    rasts = list(r2_HB, r3_HB),
    rast_max_q = 0.99,
    count_max_q = 0.99,
    max_surveys_q = 0.80
  )
  
  A2_HB_assessment <- assess_region(
    region = HB_ontario,
    region_boundaries = atlas_regions,
    sp_dat = A2_HB,
    rast = r2_HB,
    hex_grid = HB_assessment_hex_grid,
    water = NULL,
    rast_max = HB_assessment_scale_limits$rast_max,
    max_count_per_effort = HB_assessment_scale_limits$max_count_per_effort,
    max_surveys = HB_assessment_scale_limits$max_surveys,
    transform = "identity",
    title = paste0(sp_english," - Atlas 2 - Hudson Bay"),
    model_source = NULL,
    data_source = NULL
  )
  
  A3_HB_assessment <- assess_region(
    region = HB_ontario,
    region_boundaries = atlas_regions,
    sp_dat = A3_HB,
    rast = r3_HB,
    hex_grid = HB_assessment_hex_grid,
    water = NULL,
    rast_max = HB_assessment_scale_limits$rast_max,
    max_count_per_effort = HB_assessment_scale_limits$max_count_per_effort,
    max_surveys = HB_assessment_scale_limits$max_surveys,
    transform = "identity",
    title = paste0(sp_english," - Atlas 3 - Hudson Bay"),
    model_source = NULL,
    data_source = NULL
  )
  
  
  # ----------------------------------------------------------
  # PObs Maps
  # ----------------------------------------------------------
  
  pobs2 <- 1-exp(-r2)
  pobs3 <- 1-exp(-r3)
  
  pobs_absent_limit <- relabund_absent_limit
  
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
    water = NULL,
    colpal = colpal_pobs,
    water_fill = "white",
    transform = "identity",
    zlim = c(0,1),
    zbreaks = seq(0,1,length.out = 5)
  )
  
  pobs_Map_Atlas3 <- make_map(
    species_name = sp_english,
    subtitle = "Probability of Observation - Atlas 3",
    legend_title = "P(Obs)\nper 5-min\npoint count",
    rast = rasters_pobs_prepared$rasters$Atlas3,
    region = study_boundary,
    water = NULL,
    colpal = colpal_pobs,
    water_fill = "white",
    transform = "identity",
    zlim = c(0,1),
    zbreaks = seq(0,1,length.out = 5)
  )
  
  
  
  
  # ----------------------------------------------------------
  # Save one multi-page PDF per species
  # ----------------------------------------------------------
  
  # --- Page 1: Atlas 3 relative abundance maps
  page1_png <- file.path(png_dir, paste0(sp_file, "_page1.png"))
  
  ragg::agg_png(
    filename = page1_png,
    width = 20,
    height = 10,
    units = "in",
    res = 300,
    background = "white"
  )
  print(A3_assessment$plot_combined)
  grDevices::dev.off()
  
  # --- Page 2: Atlas 2 relative abundance maps
  page2_png <- file.path(png_dir, paste0(sp_file, "_page2.png"))
  
  ragg::agg_png(
    filename = page2_png,
    width = 20,
    height = 10,
    units = "in",
    res = 300,
    background = "white"
  )
  print(A2_assessment$plot_combined)
  grDevices::dev.off()
  
  # --- Page 3: Population change between atlases
  page3_png <- file.path(png_dir, paste0(sp_file, "_page3.png"))
  
  bottom_right_panel <-
    paired_locations_map |
    paired_summary_table +
    patchwork::plot_layout(widths = c(1,1))
  
  right_panel <-
    change_comparison_plot /
    bottom_right_panel +
    patchwork::plot_layout(heights = c(0.5,1))
  
  population_change_page <-
    Hex_Change_Map |
    right_panel +
    patchwork::plot_layout(widths = c(1.5,1))
  
  ragg::agg_png(
    filename = page3_png,
    width = 20,
    height = 10,
    units = "in",
    res = 300,
    background = "white"
  )
  print(population_change_page)
  grDevices::dev.off()
  
  
  # --- Page 4: Southern Ontario - atlas 3
  page4_png <- file.path(png_dir, paste0(sp_file, "_page4.png"))
  
  ragg::agg_png(
    filename = page4_png,
    width = 20,
    height = 10,
    units = "in",
    res = 300,
    background = "white"
  )
  print(A3_southern_assessment$plot_combined)
  grDevices::dev.off()
  
  # --- Page 5: Southern Ontario - atlas 2
  page5_png <- file.path(png_dir, paste0(sp_file, "_page5.png"))
  
  ragg::agg_png(
    filename = page5_png,
    width = 20,
    height = 10,
    units = "in",
    res = 300,
    background = "white"
  )
  print(A2_southern_assessment$plot_combined)
  grDevices::dev.off()
  
  # --- Page 6: Boreal Ontario - atlas 3
  page6_png <- file.path(png_dir, paste0(sp_file, "_page6.png"))
  
  ragg::agg_png(
    filename = page6_png,
    width = 20,
    height = 10,
    units = "in",
    res = 300,
    background = "white"
  )
  print(A3_boreal_assessment$plot_combined)
  grDevices::dev.off()
  
  # --- Page 7: Boreal Ontario - atlas 2
  page7_png <- file.path(png_dir, paste0(sp_file, "_page7.png"))
  
  ragg::agg_png(
    filename = page7_png,
    width = 20,
    height = 10,
    units = "in",
    res = 300,
    background = "white"
  )
  print(A2_boreal_assessment$plot_combined)
  grDevices::dev.off()
  
  
  
  
  
  
  # --- Page 8: Hudson Bay - atlas 3
  page8_png <- file.path(png_dir, paste0(sp_file, "_page8.png"))
  
  ragg::agg_png(
    filename = page8_png,
    width = 20,
    height = 10,
    units = "in",
    res = 300,
    background = "white"
  )
  print(A3_HB_assessment$plot_combined)
  grDevices::dev.off()
  
  # --- Page 9: Hudson Bay - atlas 2
  page9_png <- file.path(png_dir, paste0(sp_file, "_page9.png"))
  
  ragg::agg_png(
    filename = page9_png,
    width = 20,
    height = 10,
    units = "in",
    res = 300,
    background = "white"
  )
  print(A2_HB_assessment$plot_combined)
  grDevices::dev.off()
  
  
  # --- Page 10: PObs Maps
  page10_png <- file.path(png_dir, paste0(sp_file, "_page10.png"))
  
  pobs_page <-
    pobs_Map_Atlas3 +
    pobs_Map_Atlas2 +
    patchwork::plot_layout(ncol = 2, widths = c(1, 1))
  
  ragg::agg_png(
    filename = page10_png,
    width = 20,
    height = 10,
    units = "in",
    res = 300,
    background = "white"
  )
  print(pobs_page)
  grDevices::dev.off()
  
  
  imgs <- magick::image_read(c(page1_png, 
                               page2_png, 
                               page3_png, 
                               page4_png,
                               page5_png,
                               page6_png,
                               page7_png, 
                               page8_png,
                               page9_png,
                               page10_png))
  magick::image_write(imgs, path = pdf_path, format = "pdf")
  
}

message("08_generate_model_products.R complete")
