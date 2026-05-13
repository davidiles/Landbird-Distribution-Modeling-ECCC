# ============================================================
# 8_generate_model_products.R
# ============================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(ggplot2)
  library(terra)
  library(patchwork)
  library(magick)
})

source(here::here("R", "00_config_paths.R"))

# helper-functions
source(file.path(paths$functions, "model_product_utils.R"))

# ------------------------------------------------------------
# Config
# ------------------------------------------------------------

model_name <- "m2"

# File paths
in_data  <- file.path(paths$data_clean, "birds", "data_ready_for_analysis.rds")
pred_dir <- file.path(paths$model_output, paste0("predictions_", model_name))
rast_dir <- file.path(paths$model_output, paste0("rasters_", model_name))
fig_dir  <- file.path(paths$model_output, paste0("figures_", model_name))

dir.create(rast_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- paths$model_output

model_summaries_path  <- file.path(out_dir, paste0("summaries_", model_name), "model_summaries.rds")
model_summaries <- readRDS(model_summaries_path)

# Raster resolution in map units (meters in EPSG:3978)
plot_res <- 1001

# ------------------------------------------------------------
# Load base data
# ------------------------------------------------------------

stopifnot(file.exists(in_data))
dat <- readRDS(in_data)

# raw data
all_surveys <- dat$all_surveys
counts      <- dat$counts

# prediction grid / study area
grid2 <- dat$grid_OBBA2 %>% na.omit()
grid3 <- dat$grid_OBBA3 %>% na.omit()
study_boundary <- dat$study_boundary %>% sf::st_as_sf()
bcr_sf <- dat$bcr_sf

# Ensure consistent CRS
crs_use <- st_crs(grid2)
study_boundary <- st_transform(study_boundary, crs_use)
grid3 <- st_transform(grid3, crs_use)
grid2 <- st_transform(grid2, crs_use)

# Hexagon grid in which to summarize change estimates (25 km hexagons)
hex_grid <- dat$hex_grid_25km

# ------------------------------------------------------------
# Create a "water_to_plot" object
# Only plot waterbodies exceeding 25 km^2 in size on maps
# ------------------------------------------------------------

in_water <- file.path(paths$data_clean, "spatial", "water_to_plot.rds")
if (file.exists(in_water)){
  water_to_plot <- readRDS(in_water)
} else{
  in_water_raw <- file.path(paths$data, "Spatial", "Ontario_Hydro_Network_(OHN)_-_Waterbody","Ontario_Hydro_Network_(OHN)_-_Waterbody.shp")
  water_to_plot <- st_read(in_water, quiet = TRUE) %>%
    st_make_valid() %>%
    st_transform(crs_use) %>%
    dplyr::mutate(area_m2 = as.numeric(sf::st_area(geometry))) 
  
  water_to_plot <- water_to_plot %>%
    filter(VERIFICATI == "Verified" & area_m2 >= (5000*5000)) %>%
    st_intersection(study_boundary)
  
  saveRDS(water_to_plot,in_water)
}

# ------------------------------------------------------------
# Loop through species
# ------------------------------------------------------------

pred_files <- list.files(pred_dir, pattern = "\\.rds$", full.names = TRUE)
stopifnot(length(pred_files) > 0)

message("Prediction files found: ", length(pred_files))

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
  
  # Paths to store outputs for this species
  rast_path_a2 <- file.path(rast_dir, paste0(sp_file, "_a2.tif"))
  rast_path_a3 <- file.path(rast_dir, paste0(sp_file, "_a3.tif"))
  rast_path_chg     <- file.path(rast_dir, paste0(sp_file, "_chg.tif"))
  fig_path <- file.path(fig_dir,paste0(sp_file,".png"))
  dat_path   <- file.path(out_dir, paste0("data_used_", model_name), paste0(sp_file, "_1km.rds"))
  
  # ----------------------------------------------------------
  # Load model fitted effects to aid interpretation
  # ----------------------------------------------------------
  
  sp_model_summary <- model_summaries[[sp_english]]
  
  # ----------------------------------------------------------
  # Load data that was used for this species
  # ----------------------------------------------------------
  
  sp_dat <- readRDS(dat_path)
  
  # ------------------------------------------------------------
  # Create rasters
  # ------------------------------------------------------------
  message("Generating rasters for: ", sp_english)
  
  # ---- Atlas 2 relative abundance
  a2 <- grid2 %>% 
    mutate(mu_mean = preds$OBBA2_Corrected_for_Water$OBBA2_mean,
           CI_95_width = preds$OBBA2$OBBA2_upper - preds$OBBA2$OBBA2_lower)
  
  vars_to_rasterize <- c("mu_mean","CI_95_width")
  r2 = rasterize_sf(a2,vars_to_rasterize,res = plot_res,
                    metadata = c(
                      species_name = sp_english,
                      species_id   = sp_code,
                      species_filename = sp_file,
                      model_name        = model_name,
                      units        = "Expected count (Atlas 2)"
                    ))
  
  writeRaster(r2,
              filename = rast_path_a2,
              overwrite = TRUE)
  
  
  # ---- Atlas 3 relative abundance
  a3 <- grid3 %>% 
    mutate(mu_mean = preds$OBBA3_Corrected_for_Water$OBBA3_mean,
           CI_95_width = preds$OBBA3$OBBA3_upper - preds$OBBA3$OBBA3_lower)
  
  vars_to_rasterize <- c("mu_mean","CI_95_width")
  r3 = rasterize_sf(a3,vars_to_rasterize,res = plot_res,
                    metadata = c(
                      species_name = sp_english,
                      species_id   = sp_code,
                      species_filename = sp_file,
                      model_name        = model_name,
                      units        = "Expected count (Atlas 3)"
                    ))
  
  writeRaster(r3,
              filename = rast_path_a3,
              overwrite = TRUE)
  
  # ---- Change from Atlas 2 to Atlas 3
  chg_sf <- grid3 %>% 
    mutate(chg_mean = preds$abs_change_Corrected_for_Water$abs_change_mean,
           CI_95_width = preds$abs_change$abs_change_upper - preds$abs_change$abs_change_lower)
  
  vars_to_rasterize <- c("chg_mean","CI_95_width")
  rchg = rasterize_sf(chg_sf,vars_to_rasterize,res = plot_res,
                      metadata = c(
                        species_name = sp_english,
                        species_id   = sp_code,
                        species_filename = sp_file,
                        model_name        = model_name,
                        units        = "Change in expected counts (Atlas 2 to Atlas 3)"
                      ))
  
  writeRaster(rchg,
              filename = rast_path_chg,
              overwrite = TRUE)
  
  # ------------------------------------------------------------
  # Create maps of relative abundance, change, and PObs
  # ------------------------------------------------------------
  
  # Clamps rasters to desired upper and lower limits
  # By default: considered "absent" if rarer than 1 detection every 250 point counts
  #             - this could be adjusted separately for each species
  rast_absent_limit <- 1/250
  rasters_prepared <- prepare_relative_abundance_rasters(
    Atlas2 = r2,
    Atlas3 = r3,
    rast_absent_limit = rast_absent_limit,
    rast_max_quantile = 0.99
  )
  
  Relabund_Map_Atlas2 <- make_relabund_map(
    species_name = sp_english,
    atlas_label = "Atlas 2",
    rast = rasters_prepared$rasters$Atlas2,
    zlim = rasters_prepared$zlim,
    zbreaks = rasters_prepared$zbreaks,
    region = study_boundary,
    water = water_to_plot,
    water_fill = "white",
    transform = "identity"
  )
  
  Relabund_Map_Atlas3 <- make_relabund_map(
    species_name = sp_english,
    atlas_label = "Atlas 3",
    rast = rasters_prepared$rasters$Atlas3,
    zlim = rasters_prepared$zlim,
    zbreaks = rasters_prepared$zbreaks,
    region = study_boundary,
    water = water_to_plot,
    water_fill = "white",
    transform = "identity"
  )


  # ------------------------------------------------------------
  # Create change map
  # ------------------------------------------------------------
  
  prov_change <- summarize_polygon_hex_draw_change(
    hex_draws = preds$hex_draws_Corrected_for_Water,
    hex_grid = hex_grid,
    polygon = study_boundary,
    ci_level = 0.90
  ) |>
    dplyr::mutate(
      pct_change_median = 100 * prop_change_median,
      pct_change_qlow   = 100 * prop_change_qlow,
      pct_change_qhigh  = 100 * prop_change_qhigh
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
    zlim = rasters_prepared$zlim,
    water = water_to_plot,
    water_fill = "gray97"
  )
  
  # ------------------------------------------------------------
  # Model Assessment Map Atlas 2
  # ------------------------------------------------------------
  
  r2_clamped <- r2$mu_mean
  r2_clamped[r2_clamped<rast_absent_limit] <- 0
  
  A2_dat_to_plot <- sp_dat$sp_dat %>%
    filter(Survey_Type %in% c("Point_Count","ARU")) %>%
    subset(Atlas == "OBBA2") %>%
    mutate(count_per_effort = count)
  
  A2_assessment <- assess_region(
    region     = study_boundary,
    sp_dat     = A2_dat_to_plot,
    rast       = r2_clamped,
    n_hexagons = 2000,
    water      = NULL,  # Do not plot water at this scale
    rast_max_q = 0.99,
    transform  = "identity",
    title      = ,
    model_source = "",
    data_source  = ""
  )
  
  # ------------------------------------------------------------
  # Model Assessment Map Atlas 3
  # ------------------------------------------------------------
  
  r3_clamped <- r3$mu_mean
  r3_clamped[r3_clamped<rast_absent_limit] <- 0
  
  A3_dat_to_plot <- sp_dat$sp_dat %>%
    filter(Survey_Type %in% c("Point_Count","ARU")) %>%
    subset(Atlas == "OBBA3") %>%
    mutate(count_per_effort = count)
  
  A3_assessment <- assess_region(
    region     = study_boundary,
    sp_dat     = A3_dat_to_plot,
    rast       = r3_clamped,
    n_hexagons = 2000,
    water      = NULL,  # Do not plot water at this scale
    rast_max_q = 0.99,
    transform  = "identity",
    title      = ,
    model_source = "",
    data_source  = ""
  )
  
  # ------------------------------------------------------------
  # Save outputs as a multi-page PDF for each species
  # ------------------------------------------------------------
  
  page1_png <- file.path(fig_dir, "PNGs",paste0(sp_file, "_page1.png"))
  page2_png <- file.path(fig_dir, "PNGs",paste0(sp_file, "_page2.png"))
  page3_png <- file.path(fig_dir, "PNGs",paste0(sp_file, "_page3.png"))
  
  pdf_path  <- file.path(fig_dir, paste0(sp_file, ".pdf"))
  
  # ---- Relative Abundance and Change Maps
  page1 <-
    Relabund_Map_Atlas2 +
    Relabund_Map_Atlas3 +
    Hex_Change_Map +
    patchwork::plot_layout(
      ncol = 3,
      widths = c(1, 1, 1)
    )
  
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
  
  # ---- Atlas 2 assessment figures
  page2 <- A2_assessment
  
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
  
  # ---- Atlas 3 assessment figures
  page3 <- A3_assessment
  
  ragg::agg_png(
    filename = page3_png,
    width = 24,
    height = 8,
    units = "in",
    res = 300,
    background = "white"
  )
  print(page3)
  grDevices::dev.off()
  
  # ---- Combine into multi-page PDF
  imgs <- magick::image_read(c(page1_png, page2_png, page3_png))
  
  magick::image_write(imgs, path = pdf_path, format = "pdf")
  
}