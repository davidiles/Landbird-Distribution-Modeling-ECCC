rm(list = ls())

suppressPackageStartupMessages({
  library(sf)
  library(ggplot2)
  library(terra)
  library(patchwork)
  library(magick)
  library(viridis)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(patchwork)
})

source(here::here("R", "00_config_paths.R"))
source(file.path(paths$functions, "model_product_utils_updated.R"))

# ------------------------------------------------------------
# Configuration and paths
# ------------------------------------------------------------

model_name <- "PC_ARU_only"

in_data  <- file.path(paths$data_clean, "birds", "data_ready_for_analysis.rds")
pred_dir <- file.path(paths$model_output, paste0("predictions_", model_name))
out_dir  <- paths$model_output

dat <- readRDS(in_data)
hex_grid <- dat$hex_grid_25km

# ------------------------------------------------------------
# Load atlas regions for summarizing relative abundance and change estimates
# ------------------------------------------------------------

atlas_regions <- st_read(file.path(paths$data, "Spatial", "Ontario_Atlas_biol_regions", "Atlas_biol_regions.shp"))

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
# Conduct comparisons
# ------------------------------------------------------------

pred_files <- list.files(pred_dir, pattern = "\\.rds$", full.names = TRUE)
stopifnot(length(pred_files) > 0)

message("Prediction files found: ", length(pred_files))

full_model_data = data.frame()
full_model_estimates = data.frame()
repeated_model_data = data.frame()
repeated_model_estimates = data.frame()


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
  
  sp_code <- dat$species_to_model |>
    dplyr::filter(english_name == sp_english) |>
    dplyr::pull(species_id)
  
  message("\nGenerating estimates for: ", sp_english)
  
  # Species-specific output paths.
  dat_path <- file.path(
    out_dir,
    paste0("data_used_", model_name),
    paste0(sp_file, "_1km.rds")
  )
  
  stopifnot(file.exists(dat_path))
  sp_dat <- readRDS(dat_path)$sp_dat
  
  # ----------------------------------------------------------
  # Estimates of change from full spatial model
  # ----------------------------------------------------------
  
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
  
  # Summary of data going into the full model
  full_detection_summary <- sp_dat %>%
    as.data.frame() %>%
    filter(count > 0) %>%
    group_by(Biol_Region,Atlas) %>%
    summarize(n_sq_detected_full = length(unique(square_id)))
  
  full_data_summary <- sp_dat %>%
    as.data.frame() %>%
    group_by(Biol_Region,Atlas, square_id) %>%
    summarize(mean_count = mean(count)) %>%
    group_by(Biol_Region,Atlas) %>%
    summarize(mean_count = mean(mean_count),
              n_squares = length(unique(square_id))) %>%
    left_join(.,full_detection_summary) %>%
    left_join(as.data.frame(atlas_regions)[,c("Biol_Regio","ECOZONE_NA")], by = c("Biol_Region" = "Biol_Regio")) 
  
  # ----------------------------------------------------------
  # Estimates of change from repeated survey analysis
  # ----------------------------------------------------------
  
  paired_data_summary <- paired_summaries[[sp_english]]$shared_data_summary
  paired_change_summary <- paired_summaries[[sp_english]]$shared_change_summary
  paired_data <- paired_summaries[[sp_english]]$shared_data
  
  repeated_detection_summary <- paired_data %>%
    as.data.frame() %>%
    dplyr::filter(count > 0) %>%
    dplyr::group_by(Biol_Region, Atlas) %>%
    dplyr::summarise(
      n_sq_detected_repeated = dplyr::n_distinct(square_id),
      .groups = "drop"
    )
  
  repeated_data_summary <- paired_data %>%
    as.data.frame() %>%
    group_by(Biol_Region,Atlas, square_id) %>%
    summarize(mean_count = mean(count)) %>%
    group_by(Biol_Region,Atlas) %>%
    summarize(mean_count = mean(mean_count),
              n_squares = length(unique(square_id))) %>%
    left_join(.,repeated_detection_summary) %>%
    left_join(as.data.frame(atlas_regions)[,c("Biol_Regio","ECOZONE_NA")], by = c("Biol_Region" = "Biol_Regio")) 
  
  # ----------------------------------------------------------
  # Compile into tables
  # ----------------------------------------------------------
  
  full_model_data <- full_model_data %>%
    bind_rows(full_data_summary %>% mutate(sp_english = sp_english, model = "Full spatial model"))
  
  repeated_model_data <- repeated_model_data %>%
    bind_rows(repeated_data_summary %>% mutate(sp_english = sp_english, model = "Repeated surveys (50 m)"))
  
  full_model_estimates <- full_model_estimates %>%
    bind_rows(regional_estimates_FullModel %>% mutate(sp_english = sp_english, model = "Full spatial model"))
  
  repeated_model_estimates <- repeated_model_estimates %>%
    bind_rows(paired_change_summary %>% mutate(sp_english = sp_english, model = "Repeated surveys (50 m)"))
  
  
}

# ============================================================
# Diagnostic comparison: full spatial model vs repeated surveys
# ============================================================

# This section compares:
#   1. Raw observed change in the full atlas dataset
#   2. Raw observed change in the repeated-survey subset
#   3. Estimated change from the full spatial model
#   4. Estimated change from the repeated-survey model
#
# Interpretation goal:
#   Determine whether the full model is more positive because:
#     A) the underlying full dataset is already more positive than the repeated subset;
#     B) the full model adds positivity beyond the raw full data;
#     C) the discrepancy is strongest in sparse regions/species;
#     D) the discrepancy is mainly abundance-based or also present in PObs.

# ============================================================
# Helper: calculate raw Atlas 3 vs Atlas 2 change
# ============================================================

calc_raw_change <- function(data_summary,
                            value_cols = c("mean_count")) {
  
  extra_cols <- intersect(
    c("n_sq_detected_full", "n_sq_detected_repeated"),
    names(data_summary)
  )
  
  data_summary %>%
    dplyr::select(
      sp_english,
      Biol_Region,
      ECOZONE_NA,
      Atlas,
      dplyr::all_of(value_cols),
      n_squares,
      dplyr::all_of(extra_cols)
    ) %>%
    tidyr::pivot_wider(
      names_from = Atlas,
      values_from = c(
        dplyr::all_of(value_cols),
        n_squares,
        dplyr::all_of(extra_cols)
      )
    ) %>%
    dplyr::mutate(
      raw_pct_change_mean_count =
        100 * (mean_count_OBBA3 - mean_count_OBBA2) / mean_count_OBBA2,
      
      raw_log_change_mean_count =
        log(mean_count_OBBA3 / mean_count_OBBA2),
      
      delta_n_squares =
        n_squares_OBBA3 - n_squares_OBBA2,
      
      ratio_n_squares =
        n_squares_OBBA3 / n_squares_OBBA2
    )
}

# ============================================================
# Raw observed change in full vs repeated datasets
# ============================================================

raw_full_change <- full_model_data %>%
  dplyr::select(-model) %>%
  calc_raw_change() %>%
  dplyr::rename_with(
    ~ paste0("full_", .x),
    -c(sp_english, Biol_Region, ECOZONE_NA)
  ) %>%
  dplyr::rename(
    full_detected_squares_OBBA2 = full_n_sq_detected_full_OBBA2,
    full_detected_squares_OBBA3 = full_n_sq_detected_full_OBBA3
  )

raw_repeated_change <- repeated_model_data %>%
  dplyr::select(-model) %>%
  calc_raw_change() %>%
  dplyr::rename_with(
    ~ paste0("repeated_", .x),
    -c(sp_english, Biol_Region, ECOZONE_NA)
  ) %>%
  dplyr::rename(
    repeated_detected_squares_OBBA2 = repeated_n_sq_detected_repeated_OBBA2,
    repeated_detected_squares_OBBA3 = repeated_n_sq_detected_repeated_OBBA3
  )

# ============================================================
# Model-estimated change
# ============================================================

full_model_change <- full_model_estimates %>%
  dplyr::transmute(
    sp_english,
    Biol_Region = Region_Number,
    Region_Name,
    full_model_pct_change = pct_change_median,
    full_model_lwr = pct_change_qlow,
    full_model_upr = pct_change_qhigh,
    full_model_ci_width = pct_change_qhigh - pct_change_qlow
  )

repeated_model_change <- repeated_model_estimates %>%
  dplyr::transmute(
    sp_english,
    Biol_Region,
    repeated_model_pct_change = pct_change_q50,
    repeated_model_lwr = pct_change_q025,
    repeated_model_upr = pct_change_q975,
    repeated_model_ci_width = pct_change_q975 - pct_change_q025
  )

diagnostic_comparison <- full_model_change %>%
  dplyr::left_join(
    repeated_model_change,
    by = c("sp_english", "Biol_Region")
  ) %>%
  dplyr::left_join(
    raw_full_change,
    by = c("sp_english", "Biol_Region")
  ) %>%
  dplyr::left_join(
    raw_repeated_change,
    by = c("sp_english", "Biol_Region")
  ) %>%
  dplyr::mutate(
    model_difference =
      full_model_pct_change - repeated_model_pct_change,
    
    raw_data_difference =
      full_raw_pct_change_mean_count - repeated_raw_pct_change_mean_count,
    
    full_model_minus_raw_full =
      full_model_pct_change - full_raw_pct_change_mean_count,
    
    repeated_model_minus_raw_repeated =
      repeated_model_pct_change - repeated_raw_pct_change_mean_count,
    
    mean_full_abundance =
      rowMeans(
        cbind(full_mean_count_OBBA2, full_mean_count_OBBA3),
        na.rm = TRUE
      ),
    
    min_full_squares =
      pmin(full_n_squares_OBBA2, full_n_squares_OBBA3, na.rm = TRUE),
    
    min_repeated_squares =
      pmin(repeated_n_squares_OBBA2, repeated_n_squares_OBBA3, na.rm = TRUE),
    
    min_full_detected_squares =
      pmin(
        full_detected_squares_OBBA2,
        full_detected_squares_OBBA3,
        na.rm = TRUE
      ),
    
    min_repeated_detected_squares =
      pmin(
        repeated_detected_squares_OBBA2,
        repeated_detected_squares_OBBA3,
        na.rm = TRUE
      )
  )

# ============================================================
# Filter to species-region combinations with enough support
# ============================================================

min_squares_needed <- 5

diagnostic_comparison_filtered <- diagnostic_comparison %>%
  dplyr::filter(
    min_repeated_detected_squares >= min_squares_needed,
    
    is.finite(raw_data_difference),
    is.finite(model_difference)
  )

region_levels <- c(
  "Hudson Bay",
  "Boreal",
  "Shield",
  "South Central",
  "Carolinian"
)

diagnostic_comparison_filtered <- diagnostic_comparison_filtered %>%
  dplyr::mutate(
    Region_Name = factor(Region_Name, levels = region_levels)
  )

# ============================================================
# Diagnostic 1: Is the model discrepancy already present in raw data?
# ============================================================
#
# Interpretation:
#   - Points near the 1:1 line:
#       The full-model vs repeated-model discrepancy mostly reflects differences
#       between the full dataset and the repeated-survey subset.
#
#   - Points centered near x = 0 but y > 0:
#       The raw full and repeated datasets show similar change, but the full
#       spatial model is still more positive. This points toward model structure,
#       prediction, smoothing, or aggregation effects.
#
#   - Positive intercept in the linear model:
#       The full model is systematically more positive after accounting for
#       raw data differences.

p_raw_vs_model_discrepancy <- ggplot(
  diagnostic_comparison_filtered,
  aes(
    x = raw_data_difference,
    y = model_difference
  )
) +
  geom_hline(yintercept = 0, lty = 2, colour = "grey50") +
  geom_vline(xintercept = 0, lty = 2, colour = "grey50") +
  geom_abline(slope = 1, intercept = 0, colour = "blue") +
  geom_point(alpha = 0.7) +
  facet_wrap(~ Region_Name) +
  coord_cartesian(xlim = c(-500,500),
                  ylim = c(-500,500))+
  theme_bw() +
  labs(
    title = "Do model discrepancies mirror raw-data discrepancies?",
    x = "Raw observed change difference\n(full data - repeated subset)",
    y = "Model-estimated change difference\n(full model - repeated model)"
  )
p_raw_vs_model_discrepancy

quantify_correspondence <- list()
for (RN in levels(diagnostic_comparison_filtered$Region_Name)){
  quantify_correspondence[[RN]] <- lm(
    model_difference ~ raw_data_difference,
    data = diagnostic_comparison_filtered %>% subset(Region_Name == RN)) %>% 
      summary()
}


#### INTERPRETATION OF FIGURE
# x axis:
#   full_raw_pct_change_mean_count - repeated_raw_pct_change_mean_count
#
#   Difference in observed Atlas 3 vs Atlas 2 change BEFORE fitting either model.
#
#   Positive x-values imply that the full atlas dataset itself shows more
#   positive change than the repeated-survey subset.
#
#   This reflects differences in the underlying datasets, including:
#     - survey footprint
#     - habitat representation
#     - spatial coverage
#     - accessibility bias
#     - differences between persistent sites and the full landscape
#
# y axis:
#   full_model_pct_change - repeated_model_pct_change
#
#   Difference between the modeled estimates of population change.
#
#   Positive y-values imply that the full spatial model estimated more positive
#   change than the repeated-survey analysis.
#
# Interpretation:
#
#   Each point asks:
#
#     "If the full model is more positive than the repeated model,
#      is that because the full dataset was already more positive
#      than the repeated subset?"
#
#   Points near the 1:1 line:
#     The discrepancy between models closely mirrors differences already present
#     in the underlying datasets. This suggests the discrepancy is primarily
#     data-driven rather than model-driven.
#
#   Points above the 1:1 line:
#     The full model is more positive than expected based on raw-data
#     differences alone. This may reflect:
#       - spatial smoothing
#       - extrapolation into poorly sampled areas
#       - habitat-based prediction over the full landscape
#       - prediction standardization
#       - aggregation effects
#
#   Points below the 1:1 line:
#     The modeled discrepancy is smaller than the raw-data discrepancy,
#     suggesting shrinkage or partial regularization by the models.
#
#   Upper-right quadrant:
#     The full dataset is more positive AND the full model is more positive.
#     This supports the idea that differing survey footprints or estimands
#     contribute to the discrepancy.
#
#   Upper-left quadrant:
#     The repeated raw data are more positive, but the full model is still more
#     positive. These points are especially important because they suggest
#     model structure or prediction behavior contributes to the discrepancy.
#
# Statistical interpretation:
#
#   A strong positive relationship between x and y indicates that raw-data
#   differences contribute to the discrepancy between analyses.
#   A slope near 1 implies the model discrepancy is almost entirely due to raw data discrepancy
#   A slope near 0 implies the model discrepancy is unrelated to the raw data discrepancy
#
#   A positive intercept would suggest that the full model tends to produce
#   more positive estimates even after accounting for differences in the
#   underlying datasets.
#
#   Scatter around the 1:1 line indicates that multiple processes likely
#   contribute simultaneously, including:
#     - differing survey footprints
#     - nonrepresentative repeated sites
#     - spatial smoothing
#     - extrapolation
#     - species-specific responses

# ============================================================
# Diagnostic 2: Does the full model differ from its own raw data?
# ============================================================
#
# Interpretation:
#   - Points near 1:1:
#       Full spatial model broadly reflects the raw full atlas data.
#
#   - Points mostly above 1:1:
#       Full spatial model estimates are more positive than the raw full data.
#       This suggests prediction, spatial smoothing, standardization, or
#       aggregation effects.
#
#   - Strong deviations in specific regions:
#       Region-specific extrapolation or data-coverage effects.

p_full_model_vs_raw <- ggplot(
  diagnostic_comparison_filtered,
  aes(
    x = full_raw_pct_change_mean_count,
    y = full_model_pct_change
  )
) +
  geom_hline(yintercept = 0, lty = 2, colour = "grey70") +
  geom_vline(xintercept = 0, lty = 2, colour = "grey70") +
  geom_abline(slope = 1, intercept = 0, lty = 2) +
  geom_point(alpha = 0.7) +
  facet_wrap(~ Region_Name) +
  theme_bw() +
  labs(
    title = "Full model estimates vs raw full-data change",
    x = "Raw full-data change in mean count (%)",
    y = "Full spatial model estimated change (%)"
  )

p_full_model_vs_raw

# ============================================================
# Diagnostic 3: Does the repeated model differ from its own raw data?
# ============================================================
#
# Interpretation:
#   - If repeated model is close to raw repeated change, then it behaves like a
#     relatively direct estimator of repeated-site change.
#
#   - If repeated model is strongly shrunk toward 0, then disagreement may partly
#     reflect shrinkage/uncertainty in the repeated analysis.

p_repeated_model_vs_raw <- ggplot(
  diagnostic_comparison_filtered,
  aes(
    x = repeated_raw_pct_change_mean_count,
    y = repeated_model_pct_change
  )
) +
  geom_hline(yintercept = 0, lty = 2, colour = "grey70") +
  geom_vline(xintercept = 0, lty = 2, colour = "grey70") +
  geom_abline(slope = 1, intercept = 0, lty = 2) +
  geom_point(alpha = 0.7) +
  facet_wrap(~ Region_Name) +
  theme_bw() +
  labs(
    title = "Repeated model estimates vs raw repeated-data change",
    x = "Raw repeated-subset change in mean count (%)",
    y = "Repeated-survey model estimated change (%)"
  )

p_repeated_model_vs_raw

# ============================================================
# Diagnostic 4: Is the discrepancy strongest where repeated data are sparse?
# ============================================================
#
# Interpretation:
#   - More positive discrepancies at low sample size:
#       Consistent with spatial smoothing / extrapolation in the full model, or
#       instability in the repeated analysis.
#
#   - No relationship with sample size:
#       Discrepancy is less likely to be simply sparse-data driven.

p_discrepancy_vs_repeated_squares <- ggplot(
  diagnostic_comparison_filtered,
  aes(
    x = log(min_repeated_squares),
    y = model_difference
  )
) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "loess", se = TRUE) +
  theme_bw() +
  labs(
    title = "Model discrepancy vs repeated-survey spatial support",
    x = "log(min repeated squares across atlases)",
    y = "Model-estimated change difference\n(full model - repeated model)"
  )

p_discrepancy_vs_repeated_squares

p_discrepancy_vs_repeated_surveys <- ggplot(
  diagnostic_comparison_filtered,
  aes(
    x = log(min_repeated_surveys),
    y = model_difference
  )
) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "loess", se = TRUE) +
  theme_bw() +
  labs(
    title = "Model discrepancy vs repeated-survey sample size",
    x = "log(min repeated surveys across atlases)",
    y = "Model-estimated change difference\n(full model - repeated model)"
  )

p_discrepancy_vs_repeated_surveys

# ============================================================
# Diagnostic 5: Is the discrepancy region-specific?
# ============================================================
#
# Interpretation:
#   - Larger positive discrepancies in northern regions:
#       Supports spatial extrapolation / sparse-data smoothing hypothesis.
#
#   - Larger positive discrepancies in southern regions:
#       Could indicate site-selection or habitat-composition differences between
#       full and repeated datasets.
#
#   - Similar positive discrepancy everywhere:
#       Suggests a more general estimand or model-standardization difference.


p_region_discrepancy <- ggplot(
  diagnostic_comparison_filtered,
  aes(
    x = Region_Name,
    y = model_difference
  )
) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.4) +
  coord_flip() +
  theme_bw() +
  labs(
    title = "Distribution of model discrepancies by region",
    x = NULL,
    y = "Model-estimated change difference\n(full model - repeated model)"
  )

p_region_discrepancy

# ============================================================
# Diagnostic 6: Does the discrepancy depend on species abundance?
# ============================================================
#
# Interpretation:
#   - Larger positive discrepancies for low-abundance species:
#       Consistent with smoothing, sparse-data effects, or instability in raw
#       percent-change calculations.
#
#   - Similar discrepancy across abundance:
#       More consistent with a broad estimand difference.

p_discrepancy_vs_abundance <- ggplot(
  diagnostic_comparison_filtered,
  aes(
    x = log(mean_full_abundance),
    y = model_difference
  )
) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "loess", se = TRUE) +
  theme_bw() +
  labs(
    title = "Model discrepancy vs species abundance",
    x = "log(mean full-data abundance)",
    y = "Model-estimated change difference\n(full model - repeated model)"
  )

p_discrepancy_vs_abundance

# ============================================================
# Diagnostic 7: Is the raw-data discrepancy abundance-based or occupancy-based?
# ============================================================
#
# Interpretation:
#   - If raw mean-count discrepancies are positive but PObs discrepancies are not:
#       The disagreement may be driven by abundance among occupied sites or
#       high-count observations, not broad distributional change.
#
#   - If both mean-count and PObs discrepancies are positive:
#       The full dataset likely differs from the repeated subset in both
#       abundance and occurrence.

p_raw_count_vs_pobs_discrepancy <- ggplot(
  diagnostic_comparison_filtered,
  aes(
    x = raw_PObs_difference,
    y = raw_data_difference
  )
) +
  geom_hline(yintercept = 0, lty = 2, colour = "grey60") +
  geom_vline(xintercept = 0, lty = 2, colour = "grey60") +
  geom_point(alpha = 0.7) +
  facet_wrap(~ Region_Name) +
  theme_bw() +
  labs(
    title = "Raw-data discrepancy: mean count vs PObs",
    x = "Raw PObs change difference\n(full data - repeated subset)",
    y = "Raw mean-count change difference\n(full data - repeated subset)"
  )

p_raw_count_vs_pobs_discrepancy

# ============================================================
# Diagnostic 8: Does repeated-analysis uncertainty explain discrepancies?
# ============================================================
#
# Interpretation:
#   - Larger discrepancies when repeated CI is wide:
#       Some disagreement may simply reflect weak repeated-survey information.
#
#   - Discrepancy remains positive even with narrow repeated CIs:
#       Stronger evidence for a systematic difference.

p_discrepancy_vs_repeated_uncertainty <- ggplot(
  diagnostic_comparison_filtered,
  aes(
    x = repeated_model_ci_width,
    y = model_difference
  )
) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "loess", se = TRUE) +
  theme_bw() +
  labs(
    title = "Model discrepancy vs repeated-analysis uncertainty",
    x = "Repeated model 95% CI width",
    y = "Model-estimated change difference\n(full model - repeated model)"
  )

p_discrepancy_vs_repeated_uncertainty

# ============================================================
# Compact numerical summaries
# ============================================================

diagnostic_summary_overall <- diagnostic_comparison_filtered %>%
  dplyr::summarise(
    n_species_regions = dplyr::n(),
    
    median_model_difference = median(model_difference, na.rm = TRUE),
    mean_model_difference = mean(model_difference, na.rm = TRUE),
    
    median_raw_data_difference = median(raw_data_difference, na.rm = TRUE),
    mean_raw_data_difference = mean(raw_data_difference, na.rm = TRUE),
    
    median_full_model_minus_raw_full =
      median(full_model_minus_raw_full, na.rm = TRUE),
    
    median_repeated_model_minus_raw_repeated =
      median(repeated_model_minus_raw_repeated, na.rm = TRUE),
    
    prop_full_model_more_positive =
      mean(model_difference > 0, na.rm = TRUE),
    
    prop_raw_full_more_positive =
      mean(raw_data_difference > 0, na.rm = TRUE)
  )

diagnostic_summary_overall

diagnostic_summary_by_region <- diagnostic_comparison_filtered %>%
  dplyr::group_by(Region_Name) %>%
  dplyr::summarise(
    n_species = dplyr::n(),
    
    median_model_difference = median(model_difference, na.rm = TRUE),
    median_raw_data_difference = median(raw_data_difference, na.rm = TRUE),
    
    median_full_model_minus_raw_full =
      median(full_model_minus_raw_full, na.rm = TRUE),
    
    prop_full_model_more_positive =
      mean(model_difference > 0, na.rm = TRUE),
    
    prop_raw_full_more_positive =
      mean(raw_data_difference > 0, na.rm = TRUE),
    
    .groups = "drop"
  )

diagnostic_summary_by_region











# ============================================================
# PROCESS RESULTS FOR COMPARISON
# ============================================================

# ============================================================
# 1. Prepare full-model estimates
# ============================================================

full_estimates_clean <- full_model_estimates %>%
  dplyr::transmute(
    sp_english,
    Biol_Region = Region_Number,
    Region_Name,
    
    full_pct_change = pct_change_median,
    full_lwr = pct_change_qlow,
    full_upr = pct_change_qhigh
  )

# ============================================================
# 2. Prepare repeated-survey estimates
# ============================================================

repeated_estimates_clean <- repeated_model_estimates %>%
  dplyr::transmute(
    sp_english,
    Biol_Region,
    
    repeated_pct_change = pct_change_q50,
    repeated_lwr = pct_change_q025,
    repeated_upr = pct_change_q975
  )

# ============================================================
# 3. Combine estimate tables
#    (one row per species x region)
# ============================================================

combined_estimates <- full_estimates_clean %>%
  dplyr::left_join(
    repeated_estimates_clean,
    by = c("sp_english", "Biol_Region")
  )

# ============================================================
# 4. Extract filtering information from full_model_data
# ============================================================

filtering_info <- full_model_data %>%
  dplyr::group_by(
    sp_english,
    Biol_Region
  ) %>%
  dplyr::summarise(
    n_sq_detected_full = max(n_sq_detected_full, na.rm = TRUE),
    .groups = "drop"
  )

# ============================================================
# 5. Join filtering info and filter
# ============================================================

combined_estimates_filtered <- combined_estimates %>%
  dplyr::left_join(
    filtering_info,
    by = c("sp_english", "Biol_Region")
  ) %>%
  dplyr::filter(n_sq_detected_full >= 20)

# ============================================================
# 5. Prepare for plotting
# ============================================================

# percent change -> symmetric log-change scale
pct_to_log_change <- function(x) {
  log1p(x / 100)
}

# log-change scale -> percent labels
log_change_to_pct_label <- function(x) {
  scales::percent(exp(x) - 1, accuracy = 1)
}

region_levels <- c(
  "Hudson Bay",
  "Boreal",
  "Shield",
  "South Central",
  "Carolinian"
)

combined_estimates_filtered_plot <- combined_estimates_filtered %>%
  dplyr::mutate(
    Region_Name = factor(Region_Name, levels = region_levels),
    full_pct_change_plot = pct_to_log_change(full_pct_change),
    full_lwr_plot        = pct_to_log_change(full_lwr),
    full_upr_plot        = pct_to_log_change(full_upr),
    
    repeated_pct_change_plot = pct_to_log_change(repeated_pct_change),
    repeated_lwr_plot        = pct_to_log_change(repeated_lwr),
    repeated_upr_plot        = pct_to_log_change(repeated_upr)
  )


# ============================================================
# 6. Create plotting function
# ============================================================

make_region_plot <- function(region_df) {
  
  region_name <- unique(region_df$Region_Name)
  
  ggplot(
    region_df,
    aes(
      x = repeated_pct_change_plot,
      y = full_pct_change_plot
    )
  ) +
    
    # 1:1 agreement line
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      colour = "grey50"
    ) +
    
    # Horizontal CI for repeated analysis (x-axis)
    geom_errorbarh(
      aes(
        xmin = repeated_lwr_plot,
        xmax = repeated_upr_plot
      ),
      height = 0,
      alpha = 0.5
    ) +
    
    # Vertical CI for full model (y-axis)
    geom_errorbar(
      aes(
        ymin = full_lwr_plot,
        ymax = full_upr_plot
      ),
      width = 0,
      alpha = 0.5
    ) +
    
    # Point estimates
    geom_point(size = 2) +
    
    # Back-transform axis labels to percent change
    scale_x_continuous(
      labels = log_change_to_pct_label
    ) +
    
    scale_y_continuous(
      labels = log_change_to_pct_label
    ) +
    
    coord_equal() +
    
    labs(
      title = region_name,
      x = "Repeated surveys (% change)",
      y = "Full spatial model (% change)"
    ) +
    
    theme_bw()
}

# ============================================================
# Generate one plot per region
# ============================================================

region_plots <- purrr::map(
  region_levels,
  ~ make_region_plot(
    dplyr::filter(
      combined_estimates_filtered_plot,
      Region_Name == .x
    )
  )
)

combined_plot <- patchwork::wrap_plots(region_plots)

combined_plot


subset(combined_estimates_filtered_plot,
       Region_Name == "Carolinian" & full_pct_change >= 500)
