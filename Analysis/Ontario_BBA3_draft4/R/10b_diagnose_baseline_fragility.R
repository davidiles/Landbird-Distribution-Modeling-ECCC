# ============================================================
# 10b_diagnose_baseline_fragility.R
#
# Purpose
#   One-look screen for the "near-zero baseline" problem in the full-model
#   province/BCR change estimates: proportional change = mu3/mu2 - 1 is unstable
#   and upward-skewed when the OBBA2 baseline mu2 is small AND poorly identified,
#   which is exactly the situation in sparsely surveyed (northern) regions and in
#   the province total of a northern-range species.
#
#   For every species x region it reports five fragility signals:
#     (1) baseline spread   mu2_qhigh / mu2_qlow  -- how uncertain the denominator is
#     (2) mean-minus-median gap in proportional change -- fingerprint of a heavy
#         right tail from a small/uncertain denominator (you report the median, but
#         a large gap means the median itself sits on a skewed posterior)
#     (3) baseline share    region mu2_median / province mu2_median -- how much of
#         the province population the model places in this (possibly fragile) region
#     (4) detection cross-check -- a positive modelled baseline where OBBA2 had ~no
#         detections is baseline extrapolation, not data
#     (5) floor-sensitivity delta -- how much the change estimate moves when the
#         lowest-density hexes (bottom `floor_cum_pop` of OBBA2 population) are
#         dropped from the ratio-of-sums. THE decisive test. (optional; recomputes)
#
# Inputs
#   Part A (fast, no recompute):
#     model_output/change_estimates_<model_name>/regional_change_estimates_<model_name>.rds
#   Part B (floor sensitivity; set RECOMPUTE_FLOOR = TRUE):
#     model_output/predictions_<model_name>/<species>_1km.rds   (hex_draws_Corrected_for_Water)
#     data_ready_for_analysis.rds                               (hex_grid, bcr_sf, boundary)
#
# Outputs
#   model_output/diagnostics_<model_name>/
#     baseline_fragility_by_region.csv    (one row per species x region)
#     baseline_fragility_province.csv     (one row per species: province screen)
# ============================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(sf)
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

# ---- Fragility thresholds (used only to set the composite `fragile` flag) ----
# baseline 90% interval spanning > this many-fold  -> untrustworthy denominator
spread_thresh      <- 8
# |prop_change_mean - prop_change_median| in percentage POINTS above this
# -> heavy right tail even though you report the median
gap_thresh_pct     <- 15
# a BCR with fewer than this many OBBA2 detection-squares is "low coverage"
low_cov_det_thresh <- 5
# floor-sensitivity: province median change moving by more than this many
# percentage POINTS when low-density hexes are dropped -> baseline-driven
floor_delta_thresh <- 20

# ---- Floor-sensitivity (Part B) ----
# Drop the lowest-density hexes that together hold the bottom this-fraction of
# the region's OBBA2 population, then recompute the change. Mirrors the
# cumulative-population absence floor you already use elsewhere.
RECOMPUTE_FLOOR <- TRUE
floor_cum_pop   <- 0.05      # 5% of OBBA2 population held in the lowest-density tail
floor_regions   <- "province"   # "province" (fast) or "all" (province + every BCR)
floor_species   <- NULL      # NULL = all modelled species; or a character vector

# baseline treated as effectively absent (matches summarize_polygon_hex_draw_change)
baseline_min <- 1e-5

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------

in_data  <- file.path(paths$data_clean, "birds", "data_ready_for_analysis.rds")
est_rds  <- file.path(paths$model_output, paste0("change_estimates_", model_name),
                      paste0("regional_change_estimates_", model_name, ".rds"))
pred_dir <- file.path(paths$model_output, paste0("predictions_", model_name))
out_dir  <- file.path(paths$model_output, paste0("diagnostics_", model_name))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(est_rds))

# ============================================================
# PART A -- fast screen straight off the script-08 export
# ============================================================
# Every column used here is already produced by script 08:
#   mu2_median, mu2_qlow, mu2_qhigh, prop_change_mean, prop_change_median,
#   n_sq_det_OBBA2, n_sq_det_OBBA3, Region_Number, sp_english.

est <- readRDS(est_rds)

needed <- c("sp_english", "Region_Number", "region_type",
            "mu2_median", "mu2_qlow", "mu2_qhigh",
            "prop_change_mean", "prop_change_median",
            "n_sq_det_OBBA2", "n_sq_det_OBBA3")
miss <- setdiff(needed, names(est))
if (length(miss)) {
  stop("script-08 export is missing columns: ", paste(miss, collapse = ", "),
       "\nRe-run 08 with this model; these are standard outputs of ",
       "summarize_polygon_hex_draw_change().")
}

# Province baseline per species, for the baseline-share denominator.
prov_base <- est %>%
  dplyr::filter(Region_Number == "province") %>%
  dplyr::select(sp_english, prov_mu2_median = mu2_median)

screen <- est %>%
  dplyr::left_join(prov_base, by = "sp_english") %>%
  dplyr::mutate(
    # (1) how many-fold does the 90% baseline interval span?
    mu2_spread_90 = dplyr::if_else(mu2_qlow > 0, mu2_qhigh / mu2_qlow, Inf),

    # (2) heavy-right-tail fingerprint, in percentage POINTS
    mean_minus_median_pct = 100 * (prop_change_mean - prop_change_median),

    # (3) fraction of the province population the model places in this region
    baseline_share = dplyr::if_else(prov_mu2_median > 0,
                                    mu2_median / prov_mu2_median, NA_real_),

    # (4) positive modelled baseline where OBBA2 had ~no detections
    detection_artifact = region_type == "BCR" &
      mu2_median > baseline_min & n_sq_det_OBBA2 == 0,

    low_coverage = region_type == "BCR" & n_sq_det_OBBA2 < low_cov_det_thresh
  )

# ============================================================
# PART B -- floor sensitivity (optional; recomputes from prediction files)
# ============================================================

# Recompute a region's median proportional change after dropping the
# lowest-density hexes that together hold the bottom `cum_frac` of its OBBA2
# population. Returns the floored median (proportion) and the number kept.
region_change_floored <- function(hex_draws, hex_weights,
                                  cum_frac = 0.05, baseline_min = 1e-5) {

  hw <- hex_draws %>% dplyr::inner_join(hex_weights, by = "hex_id")
  if (nrow(hw) == 0) return(list(prop_median = NA_real_, n_kept = 0L, n_drop = 0L))

  mu2_mat <- do.call(rbind, hw$mu_OBBA2)          # hex x draws
  mu3_mat <- do.call(rbind, hw$mu_OBBA3)
  w       <- hw$area_weight

  # per-hex baseline density (posterior median) and its population contribution
  hex_base <- if (requireNamespace("matrixStats", quietly = TRUE)) {
    matrixStats::rowMedians(mu2_mat)
  } else {
    apply(mu2_mat, 1, stats::median)
  }
  hex_pop <- pmax(hex_base, 0) * w

  # order low -> high density; drop the lowest-density hexes holding the bottom
  # cum_frac of population
  ord      <- order(hex_base, decreasing = FALSE)
  cum_frac_pop <- cumsum(hex_pop[ord]) / sum(hex_pop)
  drop_ord <- ord[cum_frac_pop <= cum_frac]
  keep     <- rep(TRUE, nrow(mu2_mat)); keep[drop_ord] <- FALSE

  m2 <- colSums((mu2_mat * w)[keep, , drop = FALSE], na.rm = TRUE)
  m3 <- colSums((mu3_mat * w)[keep, , drop = FALSE], na.rm = TRUE)

  ratio <- m3 / m2
  ratio[m2 <= baseline_min] <- NA_real_

  list(prop_median = stats::median(ratio, na.rm = TRUE) - 1,
       n_kept = sum(keep), n_drop = sum(!keep))
}

floor_tbl <- NULL

if (isTRUE(RECOMPUTE_FLOOR)) {

  stopifnot(file.exists(in_data), dir.exists(pred_dir))
  dat        <- readRDS(in_data)
  hex_grid   <- dat$hex_grid_25km
  crs_use    <- sf::st_crs(dat$grid_OBBA2)
  boundary   <- sf::st_transform(dat$study_boundary, crs_use)
  bcr_regions <- dat$bcr_sf %>% sf::st_transform(crs_use) %>%
    sf::st_make_valid() %>% dplyr::arrange(BCR) %>%
    dplyr::mutate(BCR = as.character(BCR))

  # precompute weights once, exactly as script 08 does
  region_weights <- c(
    list(province = compute_hex_polygon_weights(hex_grid, boundary)),
    stats::setNames(
      lapply(seq_len(nrow(bcr_regions)),
             function(j) compute_hex_polygon_weights(hex_grid, bcr_regions[j, ])),
      bcr_regions$BCR
    )
  )

  regions_to_do <- if (identical(floor_regions, "all")) {
    c("province", bcr_regions$BCR)
  } else {
    "province"
  }

  pred_files <- list.files(pred_dir, pattern = "_1km\\.rds$", full.names = TRUE)

  floor_tbl <- purrr::map_dfr(pred_files, function(pf) {
    preds <- readRDS(pf)
    hd    <- preds$hex_draws_Corrected_for_Water
    sp    <- preds$sp_name
    if (is.null(hd) || is.null(sp)) return(NULL)
    if (!is.null(floor_species) && !(sp %in% floor_species)) return(NULL)

    purrr::map_dfr(regions_to_do, function(rk) {
      w <- region_weights[[rk]]
      if (is.null(w)) return(NULL)
      fl <- region_change_floored(hd, w, cum_frac = floor_cum_pop,
                                  baseline_min = baseline_min)
      tibble::tibble(
        sp_english    = sp,
        Region_Number = rk,
        floored_prop_change_median = fl$prop_median,
        n_hex_kept  = fl$n_kept,
        n_hex_dropped = fl$n_drop
      )
    })
  })
}

# ============================================================
# Assemble the by-region table
# ============================================================

by_region <- screen %>%
  dplyr::select(sp_english, Region_Number, region_type,
                mu2_median, mu2_spread_90,
                prop_change_median, mean_minus_median_pct,
                baseline_share, n_sq_det_OBBA2, n_sq_det_OBBA3,
                detection_artifact, low_coverage)

if (!is.null(floor_tbl)) {
  by_region <- by_region %>%
    dplyr::left_join(floor_tbl, by = c("sp_english", "Region_Number")) %>%
    dplyr::mutate(
      floor_delta_pct = 100 * (floored_prop_change_median - prop_change_median)
    )
} else {
  by_region <- by_region %>%
    dplyr::mutate(floored_prop_change_median = NA_real_,
                  floor_delta_pct = NA_real_,
                  n_hex_kept = NA_integer_, n_hex_dropped = NA_integer_)
}

by_region <- by_region %>%
  dplyr::mutate(
    # composite flag: any single strong signal marks the estimate fragile
    fragile =
      (is.finite(mu2_spread_90) & mu2_spread_90 > spread_thresh) |
      (abs(mean_minus_median_pct) > gap_thresh_pct) |
      detection_artifact |
      (!is.na(floor_delta_pct) & abs(floor_delta_pct) > floor_delta_thresh)
  ) %>%
  dplyr::arrange(sp_english, dplyr::desc(region_type == "province"), Region_Number)

# ============================================================
# Assemble the per-species province screen
# ============================================================
# The headline number is the province estimate. This table says, per species,
# how exposed that number is to fragile northern extrapolation.

# how much province population sits in low-coverage BCRs
north_share <- by_region %>%
  dplyr::filter(region_type == "BCR", low_coverage) %>%
  dplyr::group_by(sp_english) %>%
  dplyr::summarise(north_baseline_share = sum(baseline_share, na.rm = TRUE),
                   .groups = "drop")

province <- screen %>%
  dplyr::filter(Region_Number == "province") %>%
  dplyr::transmute(
    sp_english,
    prov_prop_change_median = prop_change_median,
    prov_mu2_spread_90      = mu2_spread_90,
    prov_mean_minus_median_pct = mean_minus_median_pct
  ) %>%
  dplyr::left_join(north_share, by = "sp_english") %>%
  dplyr::mutate(north_baseline_share = dplyr::coalesce(north_baseline_share, 0))

if (!is.null(floor_tbl)) {
  prov_floor <- floor_tbl %>%
    dplyr::filter(Region_Number == "province") %>%
    dplyr::select(sp_english, prov_floored = floored_prop_change_median)
  province <- province %>%
    dplyr::left_join(prov_floor, by = "sp_english") %>%
    dplyr::mutate(prov_floor_delta_pct = 100 * (prov_floored - prov_prop_change_median))
} else {
  province <- province %>%
    dplyr::mutate(prov_floored = NA_real_, prov_floor_delta_pct = NA_real_)
}

province <- province %>%
  dplyr::mutate(
    verdict = dplyr::case_when(
      (is.finite(prov_mu2_spread_90) & prov_mu2_spread_90 > spread_thresh) |
        (abs(prov_mean_minus_median_pct) > gap_thresh_pct) |
        (!is.na(prov_floor_delta_pct) & abs(prov_floor_delta_pct) > floor_delta_thresh) |
        north_baseline_share > 0.5                         ~ "fragile",
      north_baseline_share > 0.2 |
        (is.finite(prov_mu2_spread_90) & prov_mu2_spread_90 > spread_thresh / 2) ~ "check",
      TRUE                                                 ~ "robust"
    )
  ) %>%
  dplyr::arrange(factor(verdict, levels = c("fragile", "check", "robust")),
                 dplyr::desc(north_baseline_share))

# ============================================================
# Write
# ============================================================

# reg_csv  <- file.path(out_dir, "baseline_fragility_by_region.csv")
# prov_csv <- file.path(out_dir, "baseline_fragility_province.csv")
# 
# utils::write.csv(by_region, reg_csv, row.names = FALSE)
# utils::write.csv(province,  prov_csv, row.names = FALSE)

# message("Wrote: ", reg_csv)
# message("Wrote: ", prov_csv)

message("\nProvince verdict counts:")
print(province %>% dplyr::count(verdict))

message("\nMost baseline-fragile province estimates:")
print(utils::head(province %>%
  dplyr::select(sp_english, verdict, prov_prop_change_median,
                north_baseline_share, prov_mu2_spread_90,
                prov_mean_minus_median_pct, prov_floor_delta_pct), 15))

message("\n10b_diagnose_baseline_fragility.R complete.")
