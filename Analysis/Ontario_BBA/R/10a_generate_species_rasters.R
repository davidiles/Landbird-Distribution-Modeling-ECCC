# ============================================================
# 10_make_species_rasters.R
#
# Purpose:
#   Generate geotiff rasters per species using the per-species
#   prediction summary objects created in script 08.
#
# Inputs:
#   data_clean/birds/data_ready_for_analysis.rds
#   data_clean/model_output/predictions/<sp>.rds
#
# Outputs:
#   rasters
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
  library(RColorBrewer)
  library(viridis)
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

model_type <- "PC_ARU_CL"

# File paths
in_data  <- file.path(paths$data_clean, "birds", "data_ready_for_analysis.rds")
pred_dir <- file.path(paths$model_output, paste0("predictions_", model_type))

rast_dir  <- file.path(paths$model_output, paste0("rasters_", model_type))
dir.create(rast_dir, recursive = TRUE, showWarnings = FALSE)

fig_dir  <- file.path(paths$model_output, paste0("figures_", model_type),"simplified_maps")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

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
grid2 <- dat$grid_OBBA2 %>% na.omit()
grid3 <- dat$grid_OBBA3 %>% na.omit()
study_boundary <- dat$study_boundary %>% sf::st_as_sf()

stopifnot(inherits(grid2, "sf"), inherits(grid3, "sf"))

# Ensure consistent CRS
crs_use <- st_crs(grid2)
study_boundary <- st_transform(study_boundary, crs_use)
grid3 <- st_transform(grid3, crs_use)
grid2 <- st_transform(grid2, crs_use)

# ------------------------------------------------------------
# Find prediction files
# ------------------------------------------------------------

pred_files <- list.files(pred_dir, pattern = "\\.rds$", full.names = TRUE)
stopifnot(length(pred_files) > 0)

message("Prediction files found: ", length(pred_files))

# ------------------------------------------------------------
# Create species rasters
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
  
  rast_path_a2 <- file.path(rast_dir, paste0(sp_file, "_a2.tif"))
  rast_path_a3 <- file.path(rast_dir, paste0(sp_file, "_a3.tif"))
  rast_path_chg     <- file.path(rast_dir, paste0(sp_file, "_chg.tif"))
  fig_path <- file.path(fig_dir,paste0(sp_file,".png"))
  
  
   # if (file.exists(fig_path)) {
   #   message("Skipping ", sp_english, ": rasters/figures already exist for this species")
   #   next
   # }
  
  # ----------------------------------------------------------
  # Safety checks
  # ----------------------------------------------------------
  required_parts <- c("OBBA2", "OBBA3", "abs_change")
  if (!all(required_parts %in% names(preds)) ||
      any(vapply(preds[required_parts], is.null, logical(1)))) {
    message("Skipping (missing OBBA2/OBBA3/abs_change summaries): ", basename(pf))
    next
  }
  
  message("Generating rasters for: ", sp_english)
  
  # ----------------------------------------------------------
  # Atlas 2 relative abundance
  # ----------------------------------------------------------
  
  a2 <- grid2 %>% 
    mutate(mu_q50 = preds$OBBA2$OBBA2_q50,
           CI_95_width = preds$OBBA2$OBBA2_upper - preds$OBBA2$OBBA2_lower)
  
  vars_to_rasterize <- c("mu_q50","CI_95_width")
  r2 = rasterize_sf(a2,vars_to_rasterize,res = plot_res,
                    metadata = c(
                      species_name = sp_english,
                      species_id   = sp_code,
                      species_filename = sp_file,
                      model_type        = model_type,
                      units        = "Expected count (Atlas 2)"
                    ))
  
  writeRaster(r2,
              filename = rast_path_a2,
              overwrite = TRUE)
  
  # ----------------------------------------------------------
  # Atlas 3 relative abundance
  # ----------------------------------------------------------
  
  # Prepare raster file
  a3 <- grid3 %>% 
    mutate(mu_q50 = preds$OBBA3$OBBA3_q50,
           CI_95_width = preds$OBBA3$OBBA3_upper - preds$OBBA3$OBBA3_lower)
  
  vars_to_rasterize <- c("mu_q50","CI_95_width")
  r3 = rasterize_sf(a3,vars_to_rasterize,res = plot_res,
                    metadata = c(
                      species_name = sp_english,
                      species_id   = sp_code,
                      species_filename = sp_file,
                      model_type        = model_type,
                      units        = "Expected count (Atlas 3)"
                    ))
  
  writeRaster(r3,
              filename = rast_path_a3,
              overwrite = TRUE)
  
  # ----------------------------------------------------------
  # Change between atlases
  # ----------------------------------------------------------
  
  chg_sf <- grid3 %>% 
    mutate(chg_q50 = preds$abs_change$abs_change_q50,
           CI_95_width = preds$abs_change$abs_change_upper - preds$abs_change$abs_change_lower)
  
  vars_to_rasterize <- c("chg_q50","CI_95_width")
  rchg = rasterize_sf(chg_sf,vars_to_rasterize,res = plot_res,
                      metadata = c(
                        species_name = sp_english,
                        species_id   = sp_code,
                        species_filename = sp_file,
                        model_type        = model_type,
                        units        = "Change in expected counts (Atlas 2 to Atlas 3)"
                      ))
  
  writeRaster(rchg,
              filename = rast_path_chg,
              overwrite = TRUE)
  
}

# ------------------------------------------------------------
# Generate figures
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
  
  rast_path_a2 <- file.path(rast_dir, paste0(sp_file, "_a2.tif"))
  rast_path_a3 <- file.path(rast_dir, paste0(sp_file, "_a3.tif"))
  rast_path_chg     <- file.path(rast_dir, paste0(sp_file, "_chg.tif"))
  fig_path <- file.path(fig_dir,paste0(sp_file,".png"))
  
  r2 <- rast(rast_path_a2)
  r3 <- rast(rast_path_a3)
  rchg <- rast(rast_path_chg)
  
  
  # if (file.exists(fig_path)) {
  #   message("Skipping ", sp_english, ": figures already exist for this species")
  #   next
  # }
  
  message("Generating figures for ",sp_english)
  
  # ----------------------------------------------------------
  # Generate simple raster plots
  # ----------------------------------------------------------
  
  zmax <- c(values(r2$mu_q50), values(r3$mu_q50)) %>%
    quantile(0.99, na.rm = TRUE)
  
  png(fig_path, width = 12, height = 4, units = "in", res = 600)
  par(mfrow = c(1, 3))
  par(mar = c(0, 0, 0, 0))
  
  # -------------------------
  # abundance maps
  # -------------------------
  
  # cap values at zmax for plotting
  r2_plot <- clamp(r2$mu_q50, lower = 0, upper = zmax)
  r3_plot <- clamp(r3$mu_q50, lower = 0, upper = zmax)
  
  # continuous-looking palette
  colramp_abund <- colorRampPalette(c("gray98", RColorBrewer::brewer.pal(7, "Greens")))(200)
  
  plot(
    r2_plot,
    col = colramp_abund,
    colNA = "gray95",
    main = paste0(sp_english, " - Atlas 2")
  )
  plot(vect(study_boundary), add = TRUE)
  
  plot(
    r3_plot,
    col = colramp_abund,
    colNA = "gray95",
    main = paste0(sp_english, " - Atlas 3")
  )
  plot(vect(study_boundary), add = TRUE)
  
  # -------------------------
  # change map
  # -------------------------
  
  # cap change between -zmax and +zmax
  rchg_plot <- clamp(rchg$chg_q50, lower = -zmax, upper = zmax)
  
  cols_chg_base <- (RColorBrewer::brewer.pal(11, "RdBu"))
  cols_chg_base[6] <- "white"
  colramp_chg <- colorRampPalette(cols_chg_base)(201)
  
  plot(
    rchg_plot,
    col = colramp_chg,
    colNA = "gray95",
    main = paste0(sp_english, " - Change"),
    range = c(-zmax,zmax)
  )
  plot(vect(study_boundary), add = TRUE)
  
  dev.off()
  
}