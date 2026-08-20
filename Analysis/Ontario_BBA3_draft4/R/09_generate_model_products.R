# ============================================================
# 09_generate_model_products.R
#
# Purpose
#   Convert fitted-model prediction RDS files into map-ready rasters and
#   multi-page species assessment PDFs.
#
# Running order and dependencies
#   Runs AFTER 08_export_regional_change_estimates.R. The province-wide and
#   per-BCR posterior change numbers shown on page 3 are READ from the table 08
#   saves (regional_change_estimates_<model_name>.rds) rather than recomputed
#   here. Recomputing them meant six polygon x hex st_intersection() calls PER
#   SPECIES (province + five BCRs), which dominated this script's runtime; the
#   numbers are identical either way, so 09 now just looks them up. The per-hex
#   change map (summarize_hex_draw_change) is still computed here -- it is cheap
#   and species-specific.
#
#   Also requires 07b's paired_summaries.rds for the page-3 comparison panel
#   (repeated-survey change, locations map, and summary table use the raw shared
#   survey points, which live only in that file).
#
# Main outputs (unchanged)
#   - model_output/rasters_<model_name>/<species>_{a2,a3,chg}.tif
#   - model_output/figures_<model_name>/<species>.pdf   (14 pages, same order)
#
# Notes
#   - Both the medians and the CI widths of the mapped surfaces are taken from
#     the water-corrected summaries.
#   - CRS choices are inherited from the saved analysis-ready data; objects are
#     aligned to the prediction-grid CRS only.
# ============================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(terra)
  library(patchwork)
  library(magick)
  library(viridis)
})

source(here::here("R", "00_config_paths.R"))
source(file.path(paths$functions, "model_product_utils.R"))
source(file.path(paths$functions, "inla_model_utils.R"))   # shared helpers

# ============================================================
# CONFIG
# ============================================================

model_name <- "PC_ARU_CL_nosite"

cfg <- list(
  # Rasterization
  plot_res = 1001,            # raster resolution in map units (m in EPSG:3978)

  # Assessment hex grids / scaling
  n_hexagons    = 2000,
  rast_max_q    = 0.99,
  count_max_q   = 0.99,
  max_surveys_q = 0.80,

  # Colour-scale transform for the relative-abundance (raster) panel. One of
  # "identity", "sqrt", "log", etc. "sqrt" compresses a long upper tail and
  # reveals low-abundance structure.
  relabund_transform = "sqrt",

  # Change summaries. ci_level MUST match the value 08 used to build the change
  # table (checked below), so the page-3 intervals match the exported table.
  ci_level   = 0.90,
  min_sq_det = 1,            # hide BCR-level change unless detected in >= this many squares

  # Relative-abundance colouring uses a CUMULATIVE-POPULATION floor: pixels are
  # whited out from the low-density tail until the coloured area retains this
  # fraction of the predicted population (see prepare_relative_abundance_rasters).
  relabund_coverage = 0.99,   # colour the pixels holding >= 99% of the birds

  # Biological absence floor (1 expected bird per 1000 point counts). Under the
  # cumulative method this is passed as a LOWER GUARD (min_absent_limit): white
  # always covers at least sub-detection densities, even if the cumulative floor
  # falls below it. The probability-of-observation panel still uses it directly
  # (that call keeps the default fixed-floor method).
  relabund_absent_limit = 1/1000, # 1/1000

  # PDF page rendering
  page_width  = 20,
  page_height = 10,
  page_res    = 300
)

# Bird Conservation Regions used for the sub-region panels.
# Codes come from bcr_regions$BCR (character). BCR 76 and 77 were formerly merged
# into a single "76_77" region; they are now separate, so each gets its own pair
# of assessment pages. make_region_geom() stops if a code listed here is absent
# from bcr_regions, catching a data file written by an older version of 06.
# List order sets the order of the region pages in the output PDF.
region_defs <- list(
  province = list(bcr = NULL,          suffix = ""),   # whole study area
  bcr12    = list(bcr = c("12"),       suffix = "Temperate Mixed"),
  bcr13    = list(bcr = c("13"),       suffix = "Lower Great Lakes"),
  bcr76    = list(bcr = c("76"),       suffix = "Boreal Shield West"),
  bcr77    = list(bcr = c("77"),       suffix = "Boreal Shield East"),
  hb       = list(bcr = c("74"),       suffix = "Hudson Plains")
)

# Display order (south -> north) for the regional change-comparison plot and the
# repeated-survey summary table on page 3. Defined once so the two panels cannot
# disagree, and checked against bcr_regions$BCR_Label below. These are BCR_Label
# values, not the shorter page-title suffixes used in region_defs above.
region_label_order <- c(
  "Lower Great Lakes / St. Lawrence Plain",
  "Temperate Mixed",
  "Boreal Shield West",
  "Boreal Shield East",
  "Hudson Plains"
)

# ------------------------------------------------------------
# Paths and directories
# ------------------------------------------------------------

in_data  <- file.path(paths$data_clean, "birds", "data_ready_for_analysis.rds")
pred_dir <- file.path(paths$model_output, paste0("predictions_", model_name))
rast_dir <- file.path(paths$model_output, paste0("rasters_", model_name))
fig_dir  <- file.path(paths$model_output, paste0("figures_", model_name))
png_dir  <- file.path(fig_dir, "PNGs")
data_used_dir <- file.path(paths$model_output, paste0("data_used_", model_name))
out_dir  <- paths$model_output

model_summaries_path <- file.path(
  out_dir, paste0("summaries_", model_name), "model_summaries.rds"
)
paired_analysis_path <- file.path(out_dir, "paired_summaries", "paired_summaries.rds")

# Regional change table produced by script 08 (read, not recomputed).
change_table_path <- file.path(
  out_dir, paste0("change_estimates_", model_name),
  paste0("regional_change_estimates_", model_name, ".rds")
)

for (dir_i in c(rast_dir, fig_dir, png_dir)) {
  dir.create(dir_i, recursive = TRUE, showWarnings = FALSE)
}

stopifnot(file.exists(in_data))
stopifnot(file.exists(model_summaries_path))

if (!file.exists(change_table_path)) {
  stop("Regional change table not found at:\n  ", change_table_path,
       "\nRun 08_export_regional_change_estimates.R first (same model_name).",
       call. = FALSE)
}

# ============================================================
# Load analysis-ready data, model summaries, and the change table
# ============================================================

dat             <- readRDS(in_data)
model_summaries <- readRDS(model_summaries_path)
regional_change_all <- readRDS(change_table_path)

# The exported table carries the ci_level 08 used. If it differs from cfg here,
# the page-3 intervals would silently disagree with the numbers plotted, so stop.
if ("ci_level" %in% names(regional_change_all)) {
  tbl_ci <- unique(regional_change_all$ci_level)
  if (length(tbl_ci) == 1 && !isTRUE(all.equal(tbl_ci, cfg$ci_level))) {
    stop("cfg$ci_level (", cfg$ci_level, ") does not match the change table's ci_level (",
         tbl_ci, "). Re-run 08 or align cfg$ci_level.", call. = FALSE)
  }
}

grid2          <- dat$grid_OBBA2 %>% na.omit()
grid3          <- dat$grid_OBBA3 %>% na.omit()
study_boundary <- dat$study_boundary
hex_grid       <- dat$hex_grid_25km

# Align all map layers to the Atlas 2 prediction-grid CRS.
crs_use        <- sf::st_crs(grid2)
study_boundary <- sf::st_transform(study_boundary, crs_use)
grid2          <- sf::st_transform(grid2, crs_use)
grid3          <- sf::st_transform(grid3, crs_use)

# ------------------------------------------------------------
# BCRs (panel boundaries; also used to mask/annotate change)
# ------------------------------------------------------------

bcr_regions <- dat$bcr_sf %>%
  sf::st_transform(crs_use) %>%
  sf::st_make_valid() %>%
  dplyr::arrange(BCR)

stopifnot(all(c("BCR", "BCR_Label") %in% names(bcr_regions)))

bcr_regions$BCR       <- as.character(bcr_regions$BCR)
bcr_regions$BCR_Label <- as.character(bcr_regions$BCR_Label)

# Guard against a data_ready_for_analysis.rds written before BCR 76_77 was split.
if (!setequal(region_label_order, bcr_regions$BCR_Label)) {
  stop(
    "region_label_order does not match the BCR labels in dat$bcr_sf.\n",
    "  region_label_order: ", paste(region_label_order, collapse = ", "), "\n",
    "  bcr_sf BCR_Label:   ", paste(sort(unique(bcr_regions$BCR_Label)), collapse = ", "), "\n",
    "Re-run 06_filter_and_finalize_surveys.R if bcr_sf still carries the merged ",
    "76_77 region.",
    call. = FALSE
  )
}

# Clipped copy reused by the repeated-survey locations map (built once).
bcr_regions_clipped <- suppressWarnings(
  bcr_regions %>% sf::st_intersection(study_boundary)
)

# Shared survey table. Per-species data_used files store only survey_id + count;
# geometry and attributes live here once and are rejoined by survey_id.
all_surveys <- dat$all_surveys %>% sf::st_transform(crs_use)

stopifnot("BCR" %in% names(all_surveys))
all_surveys$BCR <- as.character(all_surveys$BCR)

if (any(is.na(all_surveys$BCR))) {
  message(sum(is.na(all_surveys$BCR)),
          " surveys have no BCR assigned; they will be dropped from regional ",
          "change masking.")
}

# Pixel -> Biol_Region lookups for safe-date trimming (straight off 06's grids).
pixel_region_OBBA2 <- grid_region_lookup(dat$grid_OBBA2)
pixel_region_OBBA3 <- grid_region_lookup(dat$grid_OBBA3)

# ------------------------------------------------------------
# Water polygons used as map overlays (load or build once)
# ------------------------------------------------------------

in_water <- file.path(paths$data_clean, "spatial", "water_to_plot.rds")

if (file.exists(in_water)) {
  water_to_plot <- readRDS(in_water)
} else {
  in_water_raw <- file.path(
    paths$data, "Spatial",
    "Ontario_Hydro_Network_(OHN)_-_Waterbody",
    "Ontario_Hydro_Network_(OHN)_-_Waterbody.shp"
  )

  water_to_plot <- sf::st_read(in_water_raw, quiet = TRUE) %>%
    sf::st_make_valid() %>%
    sf::st_transform(crs_use) %>%
    dplyr::mutate(area_m2 = as.numeric(sf::st_area(geometry))) %>%
    dplyr::filter(VERIFICATI == "Verified", area_m2 >= (5000 * 5000)) %>%
    sf::st_intersection(study_boundary) |>
    sf::st_collection_extract("POLYGON")

  saveRDS(water_to_plot, in_water)
}

# ============================================================
# Precompute species-independent objects (built ONCE)
# ============================================================

# Relative-abundance colour palette (deuteranomaly-friendly green ramp).
colpal_relabund <- grDevices::colorRampPalette(c(
  "#FBF7E2", "#FCF8D0", "#EEF7C2",
  "#CEF2B0",
  "#94E5A0", "#51C987", "#18A065", "#008C59",
  "#007F53", "#006344"
))(100)

# Probability-of-observation palette.
colpal_pobs <- c("white", viridis::mako(100, direction = -1))

# Region geometries + their assessment hex grids. These depend only on the study
# area and n_hexagons, so they are constant across all species.
region_specs <- lapply(names(region_defs), function(nm) {
  d    <- region_defs[[nm]]
  geom <- make_region_geom(d$bcr, bcr_regions, study_boundary)
  list(
    geom   = geom,
    suffix = d$suffix,
    hex    = make_hex_grid(study_area = geom, n_hexagons = cfg$n_hexagons),
    bounds = suppressWarnings(sf::st_intersection(bcr_regions, geom))
  )
})
names(region_specs) <- names(region_defs)

# ============================================================
# Per-species products
#
# Iterate EVERY species that has a saved honeycomb record (data_used). Species
# with a fitted model AND a paired-survey change analysis get the full multi-page
# assessment PDF; all others get a honeycomb-only PDF.
# ============================================================

dat_used_files <- list.files(data_used_dir, pattern = "\\.rds$", full.names = TRUE)
stopifnot(length(dat_used_files) > 0)
message("Species honeycomb records found: ", length(dat_used_files))

# Paired-survey change estimates (raw shared points live here; the page-3 panel
# needs them, so they are read directly rather than from 08's export).
if (file.exists(paired_analysis_path)) {
  paired_summaries <- readRDS(paired_analysis_path)
} else {
  message("Paired summaries not found at: ", paired_analysis_path,
          "; species without paired results will receive honeycomb-only PDFs.")
  paired_summaries <- list()
}

for (i in seq_along(dat_used_files)) {

  dat_used <- readRDS(dat_used_files[i])

  # sp_file is taken from the filename so it always matches the names script 07
  # used for the prediction and raster files.
  sp_file    <- basename(dat_used_files[i]) |> stringr::str_remove("_1km\\.rds$")
  sp_english <- dat_used$sp_english

  pdf_path  <- file.path(fig_dir, paste0(sp_file, ".pdf"))
  pred_path <- file.path(pred_dir, paste0(sp_file, "_1km.rds"))

  # ----------------------------------------------------------
  # Decide what, if anything, needs to be generated
  #   1. PDF already exists            -> skip.
  #   2. Modelable but no predictions  -> skip (don't make a honeycomb-only PDF).
  #   3. Not modelable                 -> honeycomb-only PDF.
  #   4. Modelable + predictions exist -> full package (needs paired summary).
  # ----------------------------------------------------------

  if (file.exists(pdf_path)) {
    message("Skipping ", sp_english, "; PDF already exists at: ", pdf_path)
    next
  }

  if (!"is_modelable" %in% names(dat_used)) {
    dat_used$is_modelable <- file.exists(pred_path)
  }

  is_modelable <- isTRUE(dat_used$is_modelable)
  has_model    <- file.exists(pred_path)
  paired       <- paired_summaries[[sp_english]]
  has_paired   <- !is.null(paired)

  if (is_modelable && !has_model) {
    message("Skipping ", sp_english,
            "; species is modelable, but model predictions were not found.")
    next
  }

  if (is_modelable && has_model && !has_paired) {
    message("Skipping ", sp_english,
            "; species is modelable and predictions exist, but paired-survey ",
            "summaries were not found. Run 07b_fit_models_to_shared_footprint.R first.")
    next
  }

  # Reconstruct this species' survey sf (geometry + attributes + count) by
  # rejoining the lightweight count table to the shared all_surveys table.
  sp_surveys <- load_species_surveys(dat_used, all_surveys)

  # ----------------------------------------------------------
  # Honeycomb-only branch for species that were not modelable
  # ----------------------------------------------------------

  if (!is_modelable) {
    message("Honeycomb-only maps for ", sp_english, " (species not modelable).")

    honey_prep <- prepare_honeycomb_surveys(sp_surveys)
    A2_obs     <- honey_prep$data %>% dplyr::filter(Atlas == "OBBA2")
    A3_obs     <- honey_prep$data %>% dplyr::filter(Atlas == "OBBA3")

    honey_size_label <- honeycomb_size_label(honey_prep$per_minute)
    page_subtitle    <- if (honey_prep$per_minute) {
      "Observed survey effort and mean counts per survey-minute (no model fitted)"
    } else {
      "Observed survey effort and mean counts (no model fitted)"
    }

    honey_pages <- lapply(region_specs, function(spec) {
      h <- build_period_honeycombs(
        region_geom      = spec$geom,
        region_hex       = spec$hex,
        region_bounds    = spec$bounds,
        A2_dat           = A2_obs,
        A3_dat           = A3_obs,
        sp_english       = sp_english,
        title_suffix     = spec$suffix,
        max_surveys_q    = cfg$max_surveys_q,
        count_max_q      = cfg$count_max_q,
        count_size_label = honey_size_label
      )

      page_title <- if (nzchar(spec$suffix)) paste0(sp_english, " - ", spec$suffix) else sp_english

      (h$A3 + margin_theme) + (h$A2 + margin_theme) +
        patchwork::plot_layout(ncol = 2) +
        patchwork::plot_annotation(
          title    = page_title,
          subtitle = page_subtitle,
          theme = ggplot2::theme(
            plot.title    = ggplot2::element_text(size = 18, face = "bold"),
            plot.subtitle = ggplot2::element_text(size = 12, colour = "grey30")
          )
        )
    })

    png_paths <- file.path(
      png_dir, paste0(sp_file, "_honey_page", seq_along(honey_pages), ".png")
    )
    for (k in seq_along(honey_pages)) {
      save_page(honey_pages[[k]], png_paths[k],
                width = cfg$page_width, height = cfg$page_height, res = cfg$page_res)
    }

    imgs <- magick::image_read(png_paths)
    magick::image_write(imgs, path = pdf_path, format = "pdf")
    next
  }

  # ----------------------------------------------------------
  # Full model products
  # ----------------------------------------------------------

  preds <- readRDS(pred_path)

  # Trim out predictions for areas with no safe dates
  # NOTE: this ensures a species cannot occur in a region where it has no safe dates
  preds <- trim_predictions_to_safe_dates(
    preds        = preds,
    lookup_obba2 = pixel_region_OBBA2,
    lookup_obba3 = pixel_region_OBBA3,
    hex_no_safe_threshold = 0.5, # zero a hex when >= this fraction of its pixels are no-safe
    trim_unclassified     = FALSE
  )
  
  # Change numbers for this species come from 08's saved table.
  rc_sp <- regional_change_all %>% dplyr::filter(sp_english == !!sp_english)
  if (nrow(rc_sp) == 0) {
    message("Skipping ", sp_english,
            "; no rows in the regional change table. Re-run ",
            "08_export_regional_change_estimates.R for model ", model_name, ".")
    next
  }

  sp_code <- dat$species_detection_summaries$species_id[
    which(dat$species_detection_summaries$sp_english == sp_english)
  ][1]

  if (length(sp_code) == 0) {
    warning("No species_id found for ", sp_english, "; using NA in raster metadata.")
    sp_code <- NA_character_
  }

  # Retained as a diagnostic only.
  sp_model_summary <- model_summaries[[sp_english]]
  if (is.null(sp_model_summary)) {
    warning("No model summary found for ", sp_english, ".")
  }

  message("\nGenerating products for: ", sp_english)

  rast_path_a2  <- file.path(rast_dir, paste0(sp_file, "_a2.tif"))
  rast_path_a3  <- file.path(rast_dir, paste0(sp_file, "_a3.tif"))
  rast_path_chg <- file.path(rast_dir, paste0(sp_file, "_chg.tif"))

  # ----------------------------------------------------------
  # Raster products: Atlas 2, Atlas 3, absolute change
  # ----------------------------------------------------------

  a2 <- grid2 %>%
    dplyr::mutate(
      mu_q50      = preds$OBBA2_Corrected_for_Water$OBBA2_q50,
      CI_95_width = preds$OBBA2_Corrected_for_Water$OBBA2_upper -
        preds$OBBA2_Corrected_for_Water$OBBA2_lower
    )

  r2 <- rasterize_sf(
    grid_sf = a2, field = c("mu_q50", "CI_95_width"), res = cfg$plot_res,
    metadata = c(species_name = sp_english, species_id = sp_code,
                 species_filename = sp_file, model_name = model_name,
                 units = "Expected count (Atlas 2)")
  )
  terra::writeRaster(r2, filename = rast_path_a2, overwrite = TRUE)

  a3 <- grid3 %>%
    dplyr::mutate(
      mu_q50      = preds$OBBA3_Corrected_for_Water$OBBA3_q50,
      CI_95_width = preds$OBBA3_Corrected_for_Water$OBBA3_upper -
        preds$OBBA3_Corrected_for_Water$OBBA3_lower
    )

  r3 <- rasterize_sf(
    grid_sf = a3, field = c("mu_q50", "CI_95_width"), res = cfg$plot_res,
    metadata = c(species_name = sp_english, species_id = sp_code,
                 species_filename = sp_file, model_name = model_name,
                 units = "Expected count (Atlas 3)")
  )
  terra::writeRaster(r3, filename = rast_path_a3, overwrite = TRUE)

  # Both the median and the interval width come from the water-corrected summaries.
  chg_sf <- grid3 %>%
    dplyr::mutate(
      chg_q50     = preds$abs_change_Corrected_for_Water$abs_change_q50,
      CI_95_width = preds$abs_change_Corrected_for_Water$abs_change_upper -
        preds$abs_change_Corrected_for_Water$abs_change_lower
    )

  rchg <- rasterize_sf(
    grid_sf = chg_sf, field = c("chg_q50", "CI_95_width"), res = cfg$plot_res,
    metadata = c(species_name = sp_english, species_id = sp_code,
                 species_filename = sp_file, model_name = model_name,
                 units = "Change in expected counts (Atlas 2 to Atlas 3)")
  )
  terra::writeRaster(rchg, filename = rast_path_chg, overwrite = TRUE)

  # ----------------------------------------------------------
  # Relative-abundance rasters (shared absence floor + scale)
  # ----------------------------------------------------------

  rasters_relabund_prepared <- prepare_relative_abundance_rasters(
    Atlas2            = r2[["mu_q50"]],
    Atlas3            = r3[["mu_q50"]],
    absent_method     = "cumulative",
    coverage          = cfg$relabund_coverage,
    coverage_basis    = "each",     # retain >= coverage in BOTH periods (never erase a period's signal)
    min_absent_limit  = cfg$relabund_absent_limit,  # biological floor as a lower guard
    rast_max_quantile = 0.99
  )

  r2_clamped <- rasters_relabund_prepared$rasters$Atlas2
  r3_clamped <- rasters_relabund_prepared$rasters$Atlas3

  # ----------------------------------------------------------
  # Hex-level change map (Page 3, left)
  # ----------------------------------------------------------

  # Province annotation comes from 08's table (no polygon x hex intersection here).
  prov_change <- rc_sp %>% dplyr::filter(region_type == "province")
  stopifnot(nrow(prov_change) == 1)

  # Per-hex change surface IS still computed here (cheap; species-specific).
  hex_change_sf <- summarize_hex_draw_change(
    hex_grid  = hex_grid,
    hex_draws = preds$hex_draws_Corrected_for_Water,
    ci_level  = cfg$ci_level
  ) |>
    classify_min_supported_change()

  Hex_Change_Map <- make_hex_abs_change_map(
    species_name      = sp_english,
    hex_change_sf     = hex_change_sf,
    region            = study_boundary,
    region_boundaries = bcr_regions,
    prov_change       = prov_change,
    ci_level          = cfg$ci_level,
    zlim              = rasters_relabund_prepared$zlim,
    water             = NULL,
    water_fill        = "gray97"
  )

  # ----------------------------------------------------------
  # Regional change estimates: full model vs paired surveys (Page 3)
  # ----------------------------------------------------------

  # ---- Full spatial model, per BCR (read from 08's table) ----
  regional_estimates_FullModel <- rc_sp %>%
    dplyr::filter(region_type == "BCR") %>%
    dplyr::transmute(
      Region_Name, Region_Number,
      pct_change_median, pct_change_qlow, pct_change_qhigh,
      direction
    ) %>%
    dplyr::relocate(Region_Name, Region_Number)

  # ---- Paired analysis (surveys within shared 100 m footprint) ----
  # NOTE: 07b must also summarize by BCR; require_bcr_key() fails loudly on legacy
  # Biol_Region-keyed files.
  paired_change_summary <- paired_summaries[[sp_english]]$shared_change_summary %>%
    require_bcr_key("shared_change_summary", sp_english) %>%
    dplyr::rename(
      pct_change_median = pct_change_q50,
      pct_change_qlow   = pct_change_q05,
      pct_change_qhigh  = pct_change_q95
    )

  paired_data <- paired_summaries[[sp_english]]$shared_data %>%
    require_bcr_key("shared_data", sp_english)

  # Hide change estimates where a species was detected in too few squares. The
  # detection counts are read from 08's table (identical to recomputing them from
  # sp_surveys), keyed by BCR.
  data_availability <- rc_sp %>%
    dplyr::filter(region_type == "BCR") %>%
    dplyr::transmute(
      BCR = Region_Number,
      n_sq_det_OBBA2, n_sq_det_OBBA3,
      sufficient_data = n_sq_det_OBBA2 >= cfg$min_sq_det |
        n_sq_det_OBBA3 >= cfg$min_sq_det
    )

  change_estimate_cols <- c("pct_change_median", "pct_change_qlow", "pct_change_qhigh")

  regional_estimates_FullModel_masked <- regional_estimates_FullModel %>%
    dplyr::left_join(data_availability, by = c("Region_Number" = "BCR")) %>%
    dplyr::mutate(sufficient_data = dplyr::coalesce(sufficient_data, FALSE)) %>%
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(change_estimate_cols),
        ~ dplyr::if_else(sufficient_data, .x, NA_real_)
      ),
      direction = dplyr::if_else(sufficient_data, direction, NA_character_)
    )

  paired_change_summary_masked <- paired_change_summary %>%
    dplyr::left_join(data_availability, by = c("BCR" = "BCR")) %>%
    dplyr::mutate(sufficient_data = dplyr::coalesce(sufficient_data, FALSE)) %>%
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(change_estimate_cols),
        ~ dplyr::if_else(sufficient_data, .x, NA_real_)
      )
    )

  change_comparison_plot <- make_change_comparison_plot(
    regional_estimates_FullModel = regional_estimates_FullModel_masked,
    paired_change_summary        = paired_change_summary_masked,
    region_order = region_label_order,
    title    = "Regional population change estimates",
    subtitle = "Comparison of full atlas model vs repeated surveys only"
  )

  paired_locations_map <-
    ggplot2::ggplot() +
    ggplot2::geom_sf(data = study_boundary, colour = "black") +
    ggplot2::geom_sf(data = paired_data, size = 0.1, col = "#1FAA59") +
    ggplot2::geom_sf(data = bcr_regions_clipped, fill = "transparent",
                     colour = "black", linewidth = 0.3) +
    ggplot2::ggtitle("Repeated survey locations") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(size = 13),
      plot.margin = ggplot2::margin(5, 5, 5, 5)
    )

  paired_data_summary <- paired_data %>%
    as.data.frame() %>%
    dplyr::group_by(BCR, Atlas) %>%
    dplyr::summarize(
      mean_count = mean(count),
      PObs       = mean(count > 0),
      n_surveys  = dplyr::n(),
      n_squares  = length(unique(square_id)),
      .groups    = "drop"
    ) %>%
    dplyr::left_join(
      as.data.frame(bcr_regions)[, c("BCR", "BCR_Label")],
      by = c("BCR" = "BCR")
    ) %>%
    dplyr::mutate(BCR_Label = factor(BCR_Label, levels = region_label_order))

  paired_summary_table <- make_paired_summary_table(paired_data_summary)

  # ----------------------------------------------------------
  # Survey data prepared once, then assessed per region
  # ----------------------------------------------------------

  honey_prep     <- prepare_honeycomb_surveys(sp_surveys)
  A2_dat_to_plot <- honey_prep$data %>% dplyr::filter(Atlas == "OBBA2")
  A3_dat_to_plot <- honey_prep$data %>% dplyr::filter(Atlas == "OBBA3")
  honey_size_label <- honeycomb_size_label(honey_prep$per_minute)

  # One unified call per region defined in region_defs.
  assessments <- lapply(region_specs, function(spec) {
    build_period_assessments(
      region_geom       = spec$geom,
      region_hex        = spec$hex,
      r2                = r2_clamped,
      r3                = r3_clamped,
      A2_dat            = A2_dat_to_plot,
      A3_dat            = A3_dat_to_plot,
      region_boundaries = bcr_regions,
      sp_english        = sp_english,
      title_suffix      = spec$suffix,
      colpal            = colpal_relabund,
      rast_max_q        = cfg$rast_max_q,
      count_max_q       = cfg$count_max_q,
      max_surveys_q     = cfg$max_surveys_q,
      transform         = cfg$relabund_transform,
      count_size_label  = honey_size_label
    )
  })
  names(assessments) <- names(region_specs)

  # ----------------------------------------------------------
  # Probability-of-observation maps (Page 14)
  # ----------------------------------------------------------

  pobs2 <- 1 - exp(-r2)
  pobs3 <- 1 - exp(-r3)

  rasters_pobs_prepared <- prepare_relative_abundance_rasters(
    Atlas2            = pobs2,
    Atlas3            = pobs3,
    rast_absent_limit = cfg$relabund_absent_limit,
    rast_max_quantile = 0.99
  )

  pobs_Map_Atlas2 <- make_map(
    species_name = sp_english, subtitle = "Probability of Observation - Atlas 2",
    legend_title = "P(Obs)\nper 5-min\npoint count",
    rast = rasters_pobs_prepared$rasters$Atlas2, region = study_boundary,
    water = NULL, colpal = colpal_pobs, water_fill = "white",
    transform = "identity", zlim = c(0, 1), zbreaks = seq(0, 1, length.out = 5)
  )

  pobs_Map_Atlas3 <- make_map(
    species_name = sp_english, subtitle = "Probability of Observation - Atlas 3",
    legend_title = "P(Obs)\nper 5-min\npoint count",
    rast = rasters_pobs_prepared$rasters$Atlas3, region = study_boundary,
    water = NULL, colpal = colpal_pobs, water_fill = "white",
    transform = "identity", zlim = c(0, 1), zbreaks = seq(0, 1, length.out = 5)
  )

  # ----------------------------------------------------------
  # Compose pages (BCR 76 and 77 get separate page pairs)
  # ----------------------------------------------------------

  # Page 3 layout. Intermediate variables are kept (rather than collapsing into
  # one expression) so the patchwork operator precedence is preserved.
  bottom_right_panel <-
    paired_locations_map |
    paired_summary_table +
    patchwork::plot_layout(widths = c(1, 1))

  right_panel <-
    change_comparison_plot /
    bottom_right_panel +
    patchwork::plot_layout(heights = c(0.5, 1))

  population_change_page <-
    Hex_Change_Map |
    right_panel +
    patchwork::plot_layout(widths = c(1.5, 1))

  pobs_page <-
    pobs_Map_Atlas3 + pobs_Map_Atlas2 +
    patchwork::plot_layout(ncol = 2, widths = c(1, 1))

  # Page plots in final PDF order.
  page_plots <- list(
    assessments$province$A3$plot_combined,   # page 1
    assessments$province$A2$plot_combined,   # page 2
    population_change_page,                  # page 3
    assessments$bcr12$A3$plot_combined,      # page 4   (BCR 12, Temperate Mixed)
    assessments$bcr12$A2$plot_combined,      # page 5
    assessments$bcr13$A3$plot_combined,      # page 6   (BCR 13, Lower Great Lakes)
    assessments$bcr13$A2$plot_combined,      # page 7
    assessments$bcr76$A3$plot_combined,      # page 8   (BCR 76, Boreal Shield West)
    assessments$bcr76$A2$plot_combined,      # page 9
    assessments$bcr77$A3$plot_combined,      # page 10  (BCR 77, Boreal Shield East)
    assessments$bcr77$A2$plot_combined,      # page 11
    assessments$hb$A3$plot_combined,         # page 12  (BCR 74, Hudson Plains)
    assessments$hb$A2$plot_combined,         # page 13
    pobs_page                                # page 14
  )

  png_paths <- file.path(
    png_dir, paste0(sp_file, "_page", seq_along(page_plots), ".png")
  )

  for (k in seq_along(page_plots)) {
    save_page(page_plots[[k]], png_paths[k],
              width = cfg$page_width, height = cfg$page_height, res = cfg$page_res)
  }

  imgs <- magick::image_read(png_paths)
  magick::image_write(imgs, path = pdf_path, format = "pdf")
}

message("09_generate_model_products.R complete")
