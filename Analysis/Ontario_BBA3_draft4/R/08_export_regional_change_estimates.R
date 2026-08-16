# ============================================================
# 08_export_regional_change_estimates.R
#
# Purpose
#   Compute and SAVE the province-wide and BCR-specific posterior population
#   change estimates for every modelled species. These tables are the
#   quantitative backbone of the workflow: script 09 READS the regional table
#   produced here to annotate its figures (instead of recomputing the same
#   geometry-heavy summaries), and script 10 compares the exported estimates
#   against independent trends (e.g. the Breeding Bird Survey, BBS).
#
#   Running order: this script comes BEFORE 09_generate_model_products.R. It only
#   needs 07's prediction files (and, optionally, 07b's paired summaries); it does
#   not depend on any figures.
#
# Speed
#   The per-region area weights used to aggregate hex draws are PURELY geometric
#   (they depend on hex_grid + region polygon, not on any species). They are
#   therefore computed ONCE per region here and reused across all species via the
#   `hex_weights` argument to summarize_polygon_hex_draw_change(), which skips the
#   per-species st_intersection() that used to dominate runtime.
#
# Main outputs
#   model_output/change_estimates_<model_name>/
#     regional_change_estimates_<model_name>.csv   (tidy, one row per species x region)
#     regional_change_estimates_<model_name>.rds   (same tibble)
#     regional_change_draws_<model_name>.rds        (optional; posterior draws)
#   model_output/change_estimates/
#     paired_change_estimates.csv / .rds            (optional; repeated-survey only,
#                                                    MODEL-INDEPENDENT -- see below)
#
# Notes
#   - "province" is the whole study boundary; the per-BCR rows loop over every
#     row of bcr_regions (BCR 12, 13, 74, 76, 77).
#   - Estimates are NOT masked by default. Script 09 hides a BCR's change on the
#     figure when the species was detected in too few squares (min_sq_det); for a
#     quantitative export it is more useful to keep the raw estimate and carry the
#     detection counts + a `sufficient_data` flag as columns, so you can filter
#     yourself. Set MASK_INSUFFICIENT <- TRUE to reproduce the figure masking.
#   - pct_change_* columns are 100 * prop_change_* (proportional change, %).
#   - Optional annualised trend columns are derived from the symmetric log-ratio
#     sym_change = log(mu3 / mu2), which annualises cleanly through its quantiles.
#     Produced only if BOTH atlas reference years are set below -- VERIFY those
#     years against your atlas definitions before using them.
#   - The paired ("shared-footprint") export from 07b is a STANDALONE analysis,
#     not a variant of this model, so it is written to a model-independent folder
#     (change_estimates/) with no model_name suffix. Script 10 loads it as the
#     plain source "paired".
# ============================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(here)
})

source(here::here("R", "00_config_paths.R"))
source(file.path(paths$functions, "model_product_utils.R"))
source(file.path(paths$functions, "inla_model_utils.R"))

# ============================================================
# CONFIG
# ============================================================

model_name <- "PC_ARU_nosite"

# Credible-interval level for the change summaries. 0.90 -> qlow = 5th
# percentile, qhigh = 95th percentile. Script 09 reads these intervals for its
# figures, so keep this consistent across a run.
ci_level <- 0.90

# Reproduce the figure-level masking of low-information BCR estimates?
#   FALSE (default): keep every raw estimate; add detection counts +
#                    `sufficient_data` so you can filter downstream.
#   TRUE           : NA-out pct/prop/abs/sym change for BCRs where the species
#                    was detected in < min_sq_det squares in BOTH atlases.
MASK_INSUFFICIENT <- FALSE
min_sq_det        <- 1        # squares-with-detection threshold (per BCR)

# Also export the repeated-survey-only ("paired" / shared-footprint) change
# estimates from 07b, when available. These are often the MOST directly
# comparable quantity to BBS, since both are trends from revisited stations.
INCLUDE_PAIRED <- TRUE

# Also save the raw per-region posterior draw vectors (prop_change / sym_change
# / abs_change) so you can propagate atlas uncertainty into a posterior BBS
# comparison rather than relying on point estimates + CIs. Off by default
# because the file is large (species x region x n_draws).
SAVE_DRAWS <- FALSE

# Reference years for the two atlas periods, used ONLY to derive optional
# annualised %/year trends (comparable to BBS annual trends). Set either to NA
# to skip annualisation. These are midpoints of the atlas windows -- CONFIRM
# they match your atlas definitions before trusting the annualised columns.
atlas_year_OBBA2 <- 2003
atlas_year_OBBA3 <- 2023

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------

in_data       <- file.path(paths$data_clean, "birds", "data_ready_for_analysis.rds")
pred_dir      <- file.path(paths$model_output, paste0("predictions_", model_name))
data_used_dir <- file.path(paths$model_output, paste0("data_used_", model_name))
out_dir       <- file.path(paths$model_output, paste0("change_estimates_", model_name))

# Paired estimates are model-independent, so they live in their own folder.
paired_out_dir       <- file.path(paths$model_output, "change_estimates")
paired_analysis_path <- file.path(paths$model_output, "paired_summaries",
                                  "paired_summaries.rds")

dir.create(out_dir,        recursive = TRUE, showWarnings = FALSE)
dir.create(paired_out_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(in_data))
stopifnot(dir.exists(pred_dir))

# ============================================================
# Load shared inputs (once)
# ============================================================

dat <- readRDS(in_data)

study_boundary <- dat$study_boundary
hex_grid       <- dat$hex_grid_25km          # native CRS; weights are computed in it
all_surveys    <- dat$all_surveys

# Align vector layers to the prediction-grid CRS.
crs_use        <- sf::st_crs(dat$grid_OBBA2)
study_boundary <- sf::st_transform(study_boundary, crs_use)

bcr_regions <- dat$bcr_sf %>%
  sf::st_transform(crs_use) %>%
  sf::st_make_valid() %>%
  dplyr::arrange(BCR)

stopifnot(all(c("BCR", "BCR_Label") %in% names(bcr_regions)))
bcr_regions$BCR       <- as.character(bcr_regions$BCR)
bcr_regions$BCR_Label <- as.character(bcr_regions$BCR_Label)

# For the sufficiency flag.
all_surveys <- sf::st_transform(all_surveys, crs_use)
stopifnot("BCR" %in% names(all_surveys))
all_surveys$BCR <- as.character(all_surveys$BCR)

# Ordering vector (south -> north) for tidy output only; not load-bearing.
region_order <- c("province", bcr_regions$BCR)

# ============================================================
# Precompute per-region hex area weights ONCE (species-independent)
# ============================================================
# Each entry is a tibble of hex_id + area_weight for one region. 

message("Precomputing hex area weights for province + ", nrow(bcr_regions), " BCRs ...")

region_weights <- c(
  list(province = compute_hex_polygon_weights(hex_grid, study_boundary)),
  stats::setNames(
    lapply(seq_len(nrow(bcr_regions)),
           function(j) compute_hex_polygon_weights(hex_grid, bcr_regions[j, ])),
    bcr_regions$BCR
  )
)

# ============================================================
# Helper: all region change estimates for one species
# ============================================================

# Returns a tibble with one row for the province plus one row per BCR, carrying
# every column from summarize_polygon_hex_draw_change() plus region identifiers.
# Uses the precomputed region_weights (no per-call geometry).
summarize_species_regions <- function(hex_draws, ci_level) {

  prov <- summarize_polygon_hex_draw_change(
    hex_draws   = hex_draws,
    hex_grid    = hex_grid,
    hex_weights = region_weights[["province"]],
    ci_level    = ci_level
  ) %>%
    dplyr::mutate(region_type   = "province",
                  Region_Number = "province",
                  Region_Name   = "Province-wide",
                  .before = 1)

  per_bcr <- purrr::map_dfr(seq_len(nrow(bcr_regions)), function(j) {
    a_reg <- bcr_regions[j, ]
    summarize_polygon_hex_draw_change(
      hex_draws   = hex_draws,
      hex_grid    = hex_grid,
      hex_weights = region_weights[[as.character(a_reg$BCR)]],
      ci_level    = ci_level
    ) %>%
      dplyr::mutate(region_type   = "BCR",
                    Region_Number = as.character(a_reg$BCR),
                    Region_Name   = a_reg$BCR_Label,
                    .before = 1)
  })

  dplyr::bind_rows(prov, per_bcr)
}

# Per-BCR (and province) squares-with-detection counts for the sufficiency flag.
species_detection_availability <- function(sp_surveys) {
  df <- sp_surveys %>%
    sf::st_drop_geometry() %>%
    dplyr::filter(count > 0, !is.na(BCR))

  by_bcr <- df %>%
    dplyr::group_by(Atlas, BCR) %>%
    dplyr::summarise(n_sq_det = dplyr::n_distinct(square_id), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = Atlas, values_from = n_sq_det,
                       names_prefix = "n_sq_det_", values_fill = 0) %>%
    dplyr::rename(Region_Number = BCR)

  prov <- df %>%
    dplyr::group_by(Atlas) %>%
    dplyr::summarise(n_sq_det = dplyr::n_distinct(square_id), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = Atlas, values_from = n_sq_det,
                       names_prefix = "n_sq_det_", values_fill = 0) %>%
    dplyr::mutate(Region_Number = "province")

  dplyr::bind_rows(prov, by_bcr)
}

# ============================================================
# Iterate species
#
# Iterate the data_used records (which carry sp_english and the exact sp_file
# naming), and process a species only if its matching prediction file exists and
# carries water-corrected hex draws.
# ============================================================

dat_used_files <- list.files(data_used_dir, pattern = "\\.rds$", full.names = TRUE)
stopifnot(length(dat_used_files) > 0)
message("Species honeycomb records found: ", length(dat_used_files))

change_estimate_cols <- c(
  "abs_change_mean", "abs_change_median", "abs_change_qlow", "abs_change_qhigh",
  "prop_change_mean", "prop_change_median", "prop_change_qlow", "prop_change_qhigh",
  "sym_change_mean", "sym_change_median", "sym_change_qlow", "sym_change_qhigh",
  "pct_change_mean", "pct_change_median", "pct_change_qlow", "pct_change_qhigh"
)

all_rows   <- vector("list", length(dat_used_files))
draws_list <- if (SAVE_DRAWS) vector("list", length(dat_used_files)) else NULL

for (i in seq_along(dat_used_files)) {

  dat_used   <- readRDS(dat_used_files[i])
  sp_file    <- basename(dat_used_files[i]) |> stringr::str_remove("_1km\\.rds$")
  sp_english <- dat_used$sp_english
  pred_path  <- file.path(pred_dir, paste0(sp_file, "_1km.rds"))

  if (!file.exists(pred_path)) {
    message("Skipping ", sp_english, "; no prediction file (", basename(pred_path), ").")
    next
  }

  preds     <- readRDS(pred_path)
  hex_draws <- preds$hex_draws_Corrected_for_Water

  if (is.null(hex_draws)) {
    message("Skipping ", sp_english,
            "; prediction file has no hex_draws_Corrected_for_Water ",
            "(predates the script 07 water correction). Rerun 07 for this species.")
    next
  }

  sp_code <- dat$species_detection_summaries$species_id[
    which(dat$species_detection_summaries$sp_english == sp_english)
  ][1]
  if (length(sp_code) == 0 || is.na(sp_code)) sp_code <- NA_character_

  message("Summarising change for: ", sp_english)

  # ---- Core: province + per-BCR posterior change ----
  est <- summarize_species_regions(hex_draws, ci_level = ci_level) %>%
    dplyr::mutate(
      pct_change_mean   = 100 * prop_change_mean,
      pct_change_median = 100 * prop_change_median,
      pct_change_qlow   = 100 * prop_change_qlow,
      pct_change_qhigh  = 100 * prop_change_qhigh
    )

  # ---- Detection counts + sufficiency flag ----
  sp_surveys <- load_species_surveys(dat_used, all_surveys)
  avail      <- species_detection_availability(sp_surveys)

  est <- est %>%
    dplyr::left_join(avail, by = "Region_Number") %>%
    dplyr::mutate(
      n_sq_det_OBBA2 = dplyr::coalesce(n_sq_det_OBBA2, 0L),
      n_sq_det_OBBA3 = dplyr::coalesce(n_sq_det_OBBA3, 0L),
      sufficient_data = n_sq_det_OBBA2 >= min_sq_det | n_sq_det_OBBA3 >= min_sq_det
    )

  # Optionally reproduce the figure masking (province is never masked).
  if (isTRUE(MASK_INSUFFICIENT)) {
    est <- est %>%
      dplyr::mutate(
        mask = region_type == "BCR" & !sufficient_data,
        dplyr::across(dplyr::all_of(change_estimate_cols),
                      ~ dplyr::if_else(mask, NA_real_, .x)),
        direction = dplyr::if_else(mask, NA_character_, direction)
      ) %>%
      dplyr::select(-mask)
  }

  # ---- Optional annualised %/year trend (from the symmetric log-ratio) ----
  if (!is.na(atlas_year_OBBA2) && !is.na(atlas_year_OBBA3)) {
    n_years <- atlas_year_OBBA3 - atlas_year_OBBA2
    stopifnot(n_years > 0)
    annualise <- function(sym) 100 * (exp(sym / n_years) - 1)
    est <- est %>%
      dplyr::mutate(
        interval_years          = n_years,
        annual_trend_pct_median = annualise(sym_change_median),
        annual_trend_pct_qlow   = annualise(sym_change_qlow),
        annual_trend_pct_qhigh  = annualise(sym_change_qhigh)
      )
  }

  # ---- Species / model identifiers, put first ----
  est <- est %>%
    dplyr::mutate(
      sp_english = sp_english,
      sp_file    = sp_file,
      species_id = sp_code,
      model_name = model_name,
      ci_level   = ci_level,
      .before = 1
    ) %>%
    dplyr::mutate(
      Region_Number = factor(Region_Number, levels = region_order)
    ) %>%
    dplyr::arrange(Region_Number) %>%
    dplyr::mutate(Region_Number = as.character(Region_Number))

  all_rows[[i]] <- est

  # ---- Optional posterior draws for full uncertainty propagation ----
  if (SAVE_DRAWS) {
    prov_d <- summarize_polygon_hex_draw_change(
      hex_draws = hex_draws, hex_grid = hex_grid,
      hex_weights = region_weights[["province"]],
      ci_level = ci_level, return_draws = TRUE
    )$draws
    bcr_d <- lapply(seq_len(nrow(bcr_regions)), function(j) {
      summarize_polygon_hex_draw_change(
        hex_draws = hex_draws, hex_grid = hex_grid,
        hex_weights = region_weights[[as.character(bcr_regions$BCR[j])]],
        ci_level = ci_level, return_draws = TRUE
      )$draws
    })
    names(bcr_d) <- bcr_regions$BCR
    draws_list[[i]] <- list(
      sp_english = sp_english, sp_file = sp_file, species_id = sp_code,
      province = prov_d, bcr = bcr_d
    )
  }
}

# ============================================================
# Assemble and save the full-model regional change table
# ============================================================

regional_change <- dplyr::bind_rows(all_rows)
if (nrow(regional_change) == 0) stop("No species produced change estimates.")

n_sp <- dplyr::n_distinct(regional_change$sp_english)
message("Regional change estimates assembled for ", n_sp, " species (",
        nrow(regional_change), " species x region rows).")

csv_path <- file.path(out_dir, paste0("regional_change_estimates_", model_name, ".csv"))
rds_path <- file.path(out_dir, paste0("regional_change_estimates_", model_name, ".rds"))

utils::write.csv(regional_change, csv_path, row.names = FALSE)
saveRDS(regional_change, rds_path)
message("Wrote: ", csv_path)
message("Wrote: ", rds_path)

if (SAVE_DRAWS) {
  draws_list <- Filter(Negate(is.null), draws_list)
  draws_path <- file.path(out_dir, paste0("regional_change_draws_", model_name, ".rds"))
  saveRDS(draws_list, draws_path)
  message("Wrote: ", draws_path)
}

# ============================================================
# Optional: repeated-survey (paired / shared-footprint) estimates
#
# These come straight from 07b's paired_summaries and are per-BCR. They are a
# separate, more directly BBS-comparable quantity that does NOT belong to any
# single 07 model, so they are written to the model-independent change_estimates/
# folder with no model_name suffix. Script 10 loads them as the source "paired".
# ============================================================

if (isTRUE(INCLUDE_PAIRED) && file.exists(paired_analysis_path)) {

  paired_summaries <- readRDS(paired_analysis_path)
  bcr_lookup <- as.data.frame(bcr_regions)[, c("BCR", "BCR_Label")]

  paired_rows <- purrr::imap_dfr(paired_summaries, function(ps, sp_english) {
    scs <- ps$shared_change_summary
    if (is.null(scs)) return(NULL)
    scs <- require_bcr_key(scs, "shared_change_summary", sp_english)
    scs %>%
      dplyr::left_join(bcr_lookup, by = c("BCR" = "BCR")) %>%
      dplyr::mutate(sp_english = sp_english, .before = 1) %>%
      dplyr::rename(Region_Number = BCR, Region_Name = BCR_Label)
  })

  if (nrow(paired_rows) > 0) {
    paired_csv <- file.path(paired_out_dir, "paired_change_estimates.csv")
    paired_rds <- file.path(paired_out_dir, "paired_change_estimates.rds")
    utils::write.csv(paired_rows, paired_csv, row.names = FALSE)
    saveRDS(paired_rows, paired_rds)
    message("Wrote: ", paired_csv)
    message("Wrote: ", paired_rds)
  } else {
    message("Paired summaries present but held no shared_change_summary rows; ",
            "skipped paired export.")
  }

} else if (isTRUE(INCLUDE_PAIRED)) {
  message("Paired summaries not found at ", paired_analysis_path,
          "; skipped paired export.")
}

message("08_export_regional_change_estimates.R complete")
