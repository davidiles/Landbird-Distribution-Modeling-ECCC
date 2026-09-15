# ============================================================
# 07b_fit_models_to_shared_footprint.R   (paired analysis)
#
# Purpose
#   Estimate OBBA2 -> OBBA3 change from the subset of surveys that fall inside
#   the footprint shared by both atlases, using point counts only. This is a
#   design-based counterpart to the full spatio-temporal model in script 07, and
#   carries the same detectability smooths as that model: a day-of-year (DOY)
#   smooth alongside the time-of-day (TOD) smooth, so the two analyses adjust for
#   within-season and within-day survey timing the same way.
#
#   This "paired" analysis is a standalone alternative to the 07 model, NOT a
#   variant of it: it is run separately and never tied to a specific 07 model
#   (PC_ARU / PC_ARU_CL). Its output is therefore stored model-independently and
#   is referred to simply as "paired" downstream (08 exports it; 10 compares it).
#
# Per-species data
#   Reads the same per-species records written by script 06 as script 07 does,
#   but narrows them to point counts (5-minute surveys), then to the shared
#   footprint. Script 07 does NOT need to have been run first.
#
# Output
#   model_output/paired_summaries/paired_summaries.rds
#     one entry per species: shared_data_summary, shared_change_summary,
#     shared_data
# ============================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(sf)
  library(stringr)
  library(purrr)
  library(INLA)
  library(inlabru)
  library(here)
})

# ============================================================
# 1. Paths, utilities, and global configuration
# ============================================================

source(here::here("R", "00_config_paths.R"))
source(file.path(paths$functions, "inla_model_utils.R"))

min_detections <- 10
min_squares    <- 10

# Radius used to define the footprint shared by both atlases.
shared_radius_km <- 0.100

# Negative-binomial switch thresholds (see choose_error_family()). Evaluated on
# the point counts this script actually fits, so it may differ from script 07.
nb_switch_count <- 20
nb_switch_min_n <- 20

# Timing smooth priors.
# The paired analysis mirrors the full model's two detectability smooths: a
# time-of-day smooth AND a day-of-year smooth, each with the range held fixed and
# only sigma estimated (see fit_inla_shared_footprint_change()).
fixed_TOD_range <- 9                 # hours; held fixed, no prior
prior_TOD_sigma <- c(2, 0.1)         # P(sigma > 2) = 0.1 -> prior median ~0.6

fixed_DOY_range <- 45                # days; held fixed, no prior
prior_DOY_sigma <- c(2, 0.1)         # P(sigma > 2) = 0.1 -> prior median ~0.6

# Random-effect priors.
square_pcprec <- c(log(2), 0.1)

# INLA approximation settings used inside fit_inla_shared_footprint_change().
int_strategy <- "ccd"
strategy     <- "laplace"

rerun_species <- FALSE

# ============================================================
# 2. Input/output locations
# ============================================================

in_file     <- file.path(paths$data_clean, "birds", "data_ready_for_analysis.rds")
species_dir <- file.path(paths$data_clean, "birds", "species_data")   # written by 06

if (!file.exists(in_file)) {
  stop("Cannot find input at: ", in_file,
       "\nHave you run 06_filter_and_finalize_surveys.R?")
}

out_dir    <- paths$model_output
paired_dir <- file.path(out_dir, "paired_summaries")

purrr::walk(c(out_dir, paired_dir), dir.create, recursive = TRUE, showWarnings = FALSE)

paired_summaries_path <- file.path(paired_dir, "paired_summaries.rds")
paired_summaries      <- load_or_empty_list(paired_summaries_path)

# ============================================================
# 3. Load finalized data
# ============================================================

dat <- readRDS(in_file)

all_surveys                 <- dat$all_surveys
study_boundary              <- dat$study_boundary %>% sf::st_as_sf()
species_detection_summaries <- dat$species_detection_summaries

# Landscape-level grids
grid_OBBA2 <- dat$grid_OBBA2
grid_OBBA3 <- dat$grid_OBBA3

colnames(grid_OBBA3)

if (!dir.exists(species_dir)) {
  stop("Per-species data directory not found: ", species_dir,
       "\nRerun 06_filter_and_finalize_surveys.R.")
}

# ============================================================
# 4. Select species to model
# ============================================================
# Same selection rule as script 07, so the two analyses cover the same species.

species_sel <- select_modelable_species(
  species_detection_summaries,
  min_detections = min_detections,
  min_squares    = min_squares
)

species_all <- species_sel$species_all
model_ids   <- species_sel$model_ids

message("Species to model: ", length(model_ids))

# ============================================================
# 5. Main species loop: analysis within the shared survey footprint
# ============================================================

for (i in seq_along(model_ids)) {

  # ---- 5.1 Species identifiers and output paths ----
  sp_code <- as.character(model_ids[i])
  sp_name <- species_all$sp_english[as.character(species_all$species_id) == sp_code]

  message("\n====================\n",
          i, "/", length(model_ids), ": ", sp_name,
          " (species_id = ", sp_code, ")\n====================")

  # Skip if this species has already been run and saved.
  if (!rerun_species &&
      !is.null(paired_summaries[[sp_name]]) &&
      length(paired_summaries[[sp_name]]) > 0) {
    message("Skipping ", sp_name, "; paired analysis already exists.")
    next
  }

  sp_path <- sp_data_path(species_dir, sp_name)

  if (!file.exists(sp_path)) {
    message("Skipping ", sp_name, "; no per-species data file from script 06.")
    next
  }

  # ---- 5.2 Load species data, restricted to 5-minute point counts ----
  # Safe-date filtering and days_midpoint come from the record written by 06.
  sp <- load_sp_dat(sp_path, all_surveys, survey_types = "Point_Count")

  sp_dat <- sp$sp_dat %>% filter(Survey_Duration_Minutes == 5)

  if (nrow(sp_dat) == 0) {
    message("Skipping ", sp_name, "; no 5-minute point counts remain.")
    next
  }

  # ---- 5.3 Restrict to the footprint shared by both atlases ----
  sp_dat_shared <- make_shared_footprint_dataset(
    dat          = sp_dat,
    atlas_col    = "Atlas",
    atlas_levels = c("OBBA2", "OBBA3"),
    buffer_km    = shared_radius_km
  )

  if (nrow(sp_dat_shared) == 0) {
    message("Skipping ", sp_name, "; no surveys in the shared footprint.")
    next
  }

  n_bcr_shared <- length(unique(sp_dat_shared$BCR))

  if (n_bcr_shared == 1) {
    message(sp_name, ": shared footprint spans a single BCR (",
            unique(sp_dat_shared$BCR),
            "); fitting a single-region change estimate.")
  }

  error_family <- choose_error_family(
    count = subset(sp_dat, Survey_Type %in% c("ARU","Point_Count"))$count,
    nb_switch_count = nb_switch_count,
    nb_switch_min_n = nb_switch_min_n,
    top1pct_share   = 0.10
  )

  if (is.na(error_family)) {
    message("Skipping ", sp_name, "; no positive counts in the shared footprint.")
    next
  }

  # ---- 5.4 Fit the paired change model ----
  mod_shared <- try(
    fit_inla_shared_footprint_change(
      sp_dat_shared = sp_dat_shared,
      region_col    = "BCR",
      family        = error_family,

      fixed_TOD_range = fixed_TOD_range,
      prior_TOD_sigma = prior_TOD_sigma,

      fixed_DOY_range = fixed_DOY_range,
      prior_DOY_sigma = prior_DOY_sigma,

      square_pcprec = square_pcprec,

      int_strategy = int_strategy,
      strategy     = strategy,

      waic = TRUE,
      cpo  = FALSE
    ),
    silent = TRUE
  )

  if (inherits(mod_shared, "try-error") || is.null(mod_shared)) {
    message("Paired model failed for ", sp_name, "; skipping this species.")
    print(mod_shared)
    next
  }

  print(summary(mod_shared))

  # ---- 5.5 Summarize and save ----
  shared_change_summary <- summarize_shared_footprint_change(mod_shared) %>%
    mutate(sp_name = sp_name, shared_radius_km = shared_radius_km) %>%
    relocate(sp_name, shared_radius_km)

  shared_data_summary <- sp_dat_shared %>%
    st_drop_geometry() %>%
    mutate(sp_name = sp_name, shared_radius_km = shared_radius_km) %>%
    group_by(sp_name, shared_radius_km, BCR, Atlas) %>%
    summarize(
      n_svy      = n(),
      mean_count = mean(count),
      PObs       = mean(count > 0),
      mean_DOY   = mean(DayOfYear),
      min_DOY    = min(DayOfYear),
      max_DOY    = max(DayOfYear),
      mean_TOD   = mean(Hours_After_Reference),
      min_TOD    = min(Hours_After_Reference),
      max_TOD    = max(Hours_After_Reference),
      .groups    = "drop"
    )

  paired_summaries[[sp_name]][["shared_data_summary"]]   <- shared_data_summary
  paired_summaries[[sp_name]][["shared_change_summary"]] <- shared_change_summary
  paired_summaries[[sp_name]][["shared_data"]] <- sp_dat_shared %>%
    dplyr::select(Date_Time, Survey_Type, count, Hours_After_Reference,
                  DayOfYear, Atlas, square_id, BCR)

  save_atomic(paired_summaries, paired_summaries_path)
}

message("\n07b_fit_models_to_shared_footprint.R complete.")
















# 
# 
# # ==============================================================================
# # OPTIONAL SUPPLEMENTARY ANALYSIS TO EVALUATE HABITAT REPRESENTATIVENESS OF PAIRED OBBA2 / OBBA3 POINT-COUNT LOCATIONS
# #
# # Questions:
# #
# #   1. Are paired locations representative of ALL point-count locations?
# #   2. Are paired locations representative of the LANDSCAPE within each BCR?
# #   3. Are all point-count locations themselves representative of the landscape?
# #
# # Comparisons are made separately by:
# #   - BCR
# #   - Atlas (OBBA2 / OBBA3)
# #   - habitat covariate
# #
# # Assumptions:
# #   - dat, all_surveys, grid_OBBA2, grid_OBBA3 are already loaded as in script 07b
# #   - make_shared_footprint_dataset() has been sourced from inla_model_utils.R
# #   - CRS units are kilometres
# # ==============================================================================
# 
# library(dplyr)
# library(tidyr)
# library(purrr)
# library(sf)
# library(ggplot2)
# 
# 
# # ==============================================================================
# # 1. Habitat variables to evaluate
# # ==============================================================================
# 
# # These are the underlying landscape variables from the prediction grids.
# 
# habitat_covars <- c(
#   "ForestNeedleleaf",
#   "ForestBroadleaf",
#   "ForestMixed",
#   "Cropland",
#   "Urban",
#   "On_Road",
#   "Lake_Lg",
#   "Lake_Sm",
#   "GreatLakes",
#   "HudsonBayCoast",
#   "River_Lg",
#   "River_Sm",
#   
#   # Modeled with separate north/south effects
#   "Grassland_South",
#   "Grassland_North",
#   "Shrubland_South",
#   "Shrubland_North",
#   "Wetland_South",
#   "Wetland_North"
# )
# 
# # Keep only variables that actually occur in all required objects.
# habitat_covars <- Reduce(
#   intersect,
#   list(
#     habitat_covars,
#     names(all_surveys),
#     names(grid_OBBA2),
#     names(grid_OBBA3)
#   )
# )
# 
# message(
#   "Evaluating ", length(habitat_covars),
#   " habitat variables:\n",
#   paste(habitat_covars, collapse = ", ")
# )
# 
# 
# # ==============================================================================
# # 2. Construct the full point count and ARU datasets
# # ==============================================================================
# 
# surveys_all <- all_surveys %>%
#   filter(
#     Atlas %in% c("OBBA2", "OBBA3"),
#     Survey_Type %in% c("Point_Count", "ARU"),
#     
#     # Retain 5-minute point counts, but retain ARUs regardless of duration
#     Survey_Type == "ARU" |
#       (Survey_Type == "Point_Count" & Survey_Duration_Minutes == 5)
#   )
# 
# pc_all <- all_surveys %>%
#   filter(
#     Survey_Type == "Point_Count",
#     Survey_Duration_Minutes == 5,
#     Atlas %in% c("OBBA2", "OBBA3")
#   )
# # ==============================================================================
# # 3. Identify the paired/shared-footprint surveys
# # ==============================================================================
# 
# # Same distance tolerance used in script 07b:
# # 0.100 km = 100 m
# 
# paired_radius_km <- 0.100
# 
# pc_paired <- make_shared_footprint_dataset(
#   dat          = pc_all,
#   atlas_col    = "Atlas",
#   atlas_levels = c("OBBA2", "OBBA3"),
#   buffer_km    = paired_radius_km
# )
# 
# message(
#   "\nAll 5-min point counts: ", nrow(pc_all),
#   "\nPaired/shared point counts: ", nrow(pc_paired)
# )
# 
# 
# # ==============================================================================
# # 4. Prepare landscape grids
# # ==============================================================================
# 
# landscape <- bind_rows(
#   grid_OBBA2 %>%
#     mutate(Atlas = "OBBA2"),
#   
#   grid_OBBA3 %>%
#     mutate(Atlas = "OBBA3")
# )
# 
# 
# # ==============================================================================
# # 5. OPTIONAL BUT RECOMMENDED:
# #    avoid giving repeatedly surveyed coordinates extra habitat weight
# # ==============================================================================
# 
# # Habitat representativeness is fundamentally about LOCATION rather than
# # number of repeat observations.
# #
# # If the same coordinate was surveyed repeatedly, counting every survey would
# # make that habitat appear more prevalent in the sampling design merely because
# # that point received more visits.
# #
# # We therefore create unique location x Atlas datasets.
# #
# # Coordinates are rounded very finely here only to identify exact/nearly exact
# # repeated survey locations. This is NOT the 100-m pairing criterion.
# 
# make_unique_locations <- function(x) {
#   
#   xy <- sf::st_coordinates(x)
#   
#   x %>%
#     mutate(
#       .x = round(xy[, 1], 3),   # 0.001 km = 1 m
#       .y = round(xy[, 2], 3)
#     ) %>%
#     group_by(Atlas, .x, .y) %>%
#     slice(1) %>%
#     ungroup() %>%
#     select(-.x, -.y)
# }
# 
# surveys_all_locations    <- make_unique_locations(surveys_all)
# pc_paired_locations <- make_unique_locations(pc_paired)
# 
# 
# # ==============================================================================
# # 6. Sample sizes by BCR and atlas
# # ==============================================================================
# 
# sampling_summary <- surveys_all_locations %>%
#   st_drop_geometry() %>%
#   count(BCR, Atlas, name = "n_all_locations") %>%
#   full_join(
#     pc_paired_locations %>%
#       st_drop_geometry() %>%
#       count(BCR, Atlas, name = "n_paired_locations"),
#     by = c("BCR", "Atlas")
#   ) %>%
#   left_join(
#     landscape %>%
#       st_drop_geometry() %>%
#       count(BCR, Atlas, name = "n_landscape_pixels"),
#     by = c("BCR", "Atlas")
#   ) %>%
#   mutate(
#     across(
#       c(n_all_locations, n_paired_locations, n_landscape_pixels),
#       ~replace_na(.x, 0)
#     ),
#     
#     prop_survey_locations_paired =
#       n_paired_locations / n_all_locations
#   ) %>%
#   arrange(BCR, Atlas)
# 
# print(sampling_summary)
# 
# 
# # ==============================================================================
# # 7. Function for standardized mean difference
# # ==============================================================================
# 
# # SMD = difference in means expressed in pooled SD units.
# #
# # Interpretation:
# #       0 = identical means
# #   +/-0.1 = small difference
# #   +/-0.2 = noticeable difference
# #   +/-0.5 = large difference
# #
# # Sign:
# #   positive = first dataset has HIGHER values than reference dataset
# #   negative = first dataset has LOWER values
# 
# calc_smd <- function(x, ref) {
#   
#   x   <- x[is.finite(x)]
#   ref <- ref[is.finite(ref)]
#   
#   if (length(x) < 2 || length(ref) < 2) {
#     return(NA_real_)
#   }
#   
#   sd_pool <- sqrt(
#     (var(x) + var(ref)) / 2
#   )
#   
#   if (!is.finite(sd_pool) || sd_pool == 0) {
#     return(NA_real_)
#   }
#   
#   (mean(x) - mean(ref)) / sd_pool
# }
# 
# 
# # ==============================================================================
# # 8. Convert the three datasets to long format
# # ==============================================================================
# 
# paired_long <- pc_paired_locations %>%
#   st_drop_geometry() %>%
#   select(BCR, Atlas, all_of(habitat_covars)) %>%
#   pivot_longer(
#     cols = all_of(habitat_covars),
#     names_to = "covariate",
#     values_to = "value"
#   )
# 
# all_long <- surveys_all_locations %>%
#   st_drop_geometry() %>%
#   select(BCR, Atlas, all_of(habitat_covars)) %>%
#   pivot_longer(
#     cols = all_of(habitat_covars),
#     names_to = "covariate",
#     values_to = "value"
#   )
# 
# landscape_long <- landscape %>%
#   st_drop_geometry() %>%
#   select(BCR, Atlas, all_of(habitat_covars)) %>%
#   pivot_longer(
#     cols = all_of(habitat_covars),
#     names_to = "covariate",
#     values_to = "value"
#   )
# 
# 
# # ==============================================================================
# # 9. Univariate comparison 1:
# #    PAIRED LOCATIONS vs ALL SURVEY LOCATIONS
# # ==============================================================================
# 
# paired_vs_all <- paired_long %>%
#   group_by(BCR, Atlas, covariate) %>%
#   summarise(
#     n_paired   = sum(is.finite(value)),
#     mean_paired = mean(value, na.rm = TRUE),
#     sd_paired   = sd(value, na.rm = TRUE),
#     .groups = "drop"
#   ) %>%
#   
#   left_join(
#     all_long %>%
#       group_by(BCR, Atlas, covariate) %>%
#       summarise(
#         n_all   = sum(is.finite(value)),
#         mean_all = mean(value, na.rm = TRUE),
#         sd_all   = sd(value, na.rm = TRUE),
#         .groups = "drop"
#       ),
#     by = c("BCR", "Atlas", "covariate")
#   ) %>%
#   
#   rowwise() %>%
#   mutate(
#     pooled_sd = sqrt((sd_paired^2 + sd_all^2) / 2),
#     
#     SMD_paired_vs_all =
#       ifelse(
#         is.finite(pooled_sd) & pooled_sd > 0,
#         (mean_paired - mean_all) / pooled_sd,
#         NA_real_
#       )
#   ) %>%
#   ungroup()
# 
# 
# # ==============================================================================
# # 10. Univariate comparison 2:
# #     PAIRED LOCATIONS vs LANDSCAPE
# # ==============================================================================
# 
# paired_vs_landscape <- paired_long %>%
#   group_by(BCR, Atlas, covariate) %>%
#   summarise(
#     n_paired    = sum(is.finite(value)),
#     mean_paired = mean(value, na.rm = TRUE),
#     sd_paired   = sd(value, na.rm = TRUE),
#     .groups = "drop"
#   ) %>%
#   
#   left_join(
#     landscape_long %>%
#       group_by(BCR, Atlas, covariate) %>%
#       summarise(
#         n_landscape    = sum(is.finite(value)),
#         mean_landscape = mean(value, na.rm = TRUE),
#         sd_landscape   = sd(value, na.rm = TRUE),
#         .groups = "drop"
#       ),
#     by = c("BCR", "Atlas", "covariate")
#   ) %>%
#   
#   rowwise() %>%
#   mutate(
#     pooled_sd =
#       sqrt((sd_paired^2 + sd_landscape^2) / 2),
#     
#     SMD_paired_vs_landscape =
#       ifelse(
#         is.finite(pooled_sd) & pooled_sd > 0,
#         (mean_paired - mean_landscape) / pooled_sd,
#         NA_real_
#       )
#   ) %>%
#   ungroup()
# 
# 
# # ==============================================================================
# # 11. Univariate comparison 3:
# #     ALL SURVEY LOCATIONS vs LANDSCAPE
# # ==============================================================================
# 
# all_vs_landscape <- all_long %>%
#   group_by(BCR, Atlas, covariate) %>%
#   summarise(
#     n_all    = sum(is.finite(value)),
#     mean_all = mean(value, na.rm = TRUE),
#     sd_all   = sd(value, na.rm = TRUE),
#     .groups = "drop"
#   ) %>%
#   
#   left_join(
#     landscape_long %>%
#       group_by(BCR, Atlas, covariate) %>%
#       summarise(
#         n_landscape    = sum(is.finite(value)),
#         mean_landscape = mean(value, na.rm = TRUE),
#         sd_landscape   = sd(value, na.rm = TRUE),
#         .groups = "drop"
#       ),
#     by = c("BCR", "Atlas", "covariate")
#   ) %>%
#   
#   rowwise() %>%
#   mutate(
#     pooled_sd =
#       sqrt((sd_all^2 + sd_landscape^2) / 2),
#     
#     SMD_all_vs_landscape =
#       ifelse(
#         is.finite(pooled_sd) & pooled_sd > 0,
#         (mean_all - mean_landscape) / pooled_sd,
#         NA_real_
#       )
#   ) %>%
#   ungroup()
# 
# 
# # ==============================================================================
# # 12. Put the three comparisons together
# # ==============================================================================
# 
# habitat_balance <- paired_vs_landscape %>%
#   select(
#     BCR, Atlas, covariate,
#     mean_paired,
#     mean_landscape,
#     SMD_paired_vs_landscape
#   ) %>%
#   
#   left_join(
#     paired_vs_all %>%
#       select(
#         BCR, Atlas, covariate,
#         mean_all,
#         SMD_paired_vs_all
#       ),
#     by = c("BCR", "Atlas", "covariate")
#   ) %>%
#   
#   left_join(
#     all_vs_landscape %>%
#       select(
#         BCR, Atlas, covariate,
#         SMD_all_vs_landscape
#       ),
#     by = c("BCR", "Atlas", "covariate")
#   ) %>%
#   
#   mutate(
#     abs_SMD_paired_vs_landscape =
#       abs(SMD_paired_vs_landscape),
#     
#     abs_SMD_paired_vs_all =
#       abs(SMD_paired_vs_all),
#     
#     abs_SMD_all_vs_landscape =
#       abs(SMD_all_vs_landscape)
#   ) %>%
#   
#   mutate(BCR = factor(BCR, levels = c("13","12","76","77","74"))) %>%
#   
#   arrange(
#     BCR,
#     Atlas,
#     desc(abs_SMD_paired_vs_landscape)
#   )
# 
# print(habitat_balance, n = 100)
# 
# 
# # ==============================================================================
# # 13. BCR-level summary
# # ==============================================================================
# 
# # This produces a compact measure of the overall severity of habitat imbalance.
# #
# # mean_abs_SMD = average imbalance across covariates
# # max_abs_SMD  = worst individual covariate
# # n_SMD_gt_0.2 = number showing appreciable imbalance
# # n_SMD_gt_0.5 = number showing strong imbalance
# 
# bcr_balance_summary <- habitat_balance %>%
#   group_by(BCR, Atlas) %>%
#   summarise(
#     
#     n_covariates = sum(!is.na(SMD_paired_vs_landscape)),
#     
#     mean_abs_SMD_paired_landscape =
#       mean(abs(SMD_paired_vs_landscape), na.rm = TRUE),
#     
#     max_abs_SMD_paired_landscape =
#       max(abs(SMD_paired_vs_landscape), na.rm = TRUE),
#     
#     n_SMD_gt_0.2 =
#       sum(abs(SMD_paired_vs_landscape) > 0.2, na.rm = TRUE),
#     
#     n_SMD_gt_0.5 =
#       sum(abs(SMD_paired_vs_landscape) > 0.5, na.rm = TRUE),
#     
#     mean_abs_SMD_paired_all =
#       mean(abs(SMD_paired_vs_all), na.rm = TRUE),
#     
#     mean_abs_SMD_all_landscape =
#       mean(abs(SMD_all_vs_landscape), na.rm = TRUE),
#     
#     .groups = "drop"
#   ) %>%
#   
#   left_join(
#     sampling_summary,
#     by = c("BCR", "Atlas")
#   )
# 
# print(bcr_balance_summary)
# 
# 
# # ==============================================================================
# # 14. Plot: paired sites vs landscape
# # ==============================================================================
# 
# # This is probably the single most useful diagnostic figure.
# #
# # Each dot is one habitat variable.
# # Vertical zero line = perfect balance.
# #
# # BCRs with dots clustered tightly around zero have environmentally
# # representative paired samples.
# #
# # Systematic displacement in one direction indicates habitat selection.
# 
# ggplot(
#   habitat_balance,
#   aes(
#     x = SMD_paired_vs_landscape,
#     y = reorder(covariate, SMD_paired_vs_landscape)
#   )
# ) +
#   geom_vline(
#     xintercept = 0,
#     linetype = 2
#   ) +
#   geom_vline(
#     xintercept = c(-0.2, 0.2),
#     linetype = 3
#   ) +
#   geom_point() +
#   facet_grid(
#     BCR ~ Atlas,
#     scales = "free_y",
#     space = "free_y"
#   ) +
#   labs(
#     x = "Standardized mean difference: paired sites minus landscape",
#     y = NULL,
#     title = "Habitat representativeness of paired OBBA survey locations",
#     subtitle = "0 = paired locations have the same mean habitat composition as the BCR landscape"
#   ) +
#   theme_bw()
# 
# 
# # ==============================================================================
# # 15. Plot: distinguish PAIRING bias from GENERAL survey bias
# # ==============================================================================
# 
# balance_arrow <- habitat_balance %>%
#   select(
#     BCR,
#     Atlas,
#     covariate,
#     SMD_paired_vs_landscape,
#     SMD_all_vs_landscape
#   ) %>%
#   mutate(
#     improvement =
#       abs(SMD_paired_vs_landscape) -
#       abs(SMD_all_vs_landscape)
#   )
# 
# 
# ggplot(
#   balance_arrow,
#   aes(y = covariate)
# ) +
# 
# # Zero = perfect agreement with landscape
# geom_vline(
#   xintercept = 0,
#   linetype = 2,
#   linewidth = 0.6
# ) +
#   
#   # +/- 0.2 SMD reference lines
#   geom_vline(
#     xintercept = c(-0.2, 0.2),
#     linetype = 3,
#     linewidth = 0.4
#   ) +
#   
# geom_segment(
#   aes(
#     x = SMD_paired_vs_landscape,
#     xend = SMD_all_vs_landscape,
#     yend = covariate
#   ),
#   linewidth = 0.7,
#   arrow = arrow(
#     type = "closed",
#     length = unit(0.10, "inches")
#   )
# ) +
# 
# geom_point(
#   aes(
#     x = SMD_paired_vs_landscape,
#     colour = "Paired point counts"
#   ),
#   size = 2.8
# ) +
#   
# geom_point(
#   aes(
#     x = SMD_all_vs_landscape,
#     colour = "Full PC + ARU network"
#   ),
#   size = 2.8
# ) +
#   
# scale_colour_manual(
#   values = c(
#     "Paired point counts" = "black",
#     "Full PC + ARU network" = "orangered"
#   )
# ) +
#   
# 
# facet_grid(
#   BCR ~ Atlas,
#   scales = "free_y",
#   space = "free_y"
# ) +
#   
# labs(
#   x = "Standardized mean difference from landscape",
#   y = NULL,
#   colour = "Sampling set",
#   title = "Habitat representativeness of paired and full survey networks",
#   subtitle = paste0(
#     "Arrows point from paired point counts (black) to the full survey network (red);\n",
#     "movement toward zero indicates improved representativeness\n",
#     "\n",
#     "Negative values on x axis mean covariate is under-represented in sample\n",
#     "Positive values on x axis mean covariate is over-represented in sample"
#   )
# ) +
#   
#   theme_bw()
# 
# 
# # ==============================================================================
# # 16. Identify the most strongly under/over-sampled habitats
# # ==============================================================================
# 
# worst_covariates <- habitat_balance %>%
#   group_by(BCR, Atlas) %>%
#   slice_max(
#     abs_SMD_paired_vs_landscape,
#     n = 5,
#     with_ties = FALSE
#   ) %>%
#   ungroup() %>%
#   select(
#     BCR,
#     Atlas,
#     covariate,
#     mean_paired,
#     mean_all,
#     mean_landscape,
#     SMD_paired_vs_landscape,
#     SMD_paired_vs_all,
#     SMD_all_vs_landscape
#   ) %>%
#   mutate(
#     region_name = case_when(
#       BCR == "12" ~ "Temperate Mixed",
#       BCR == "13" ~ "Lower Great Lakes",
#       BCR == "74" ~ "Hudson Plains",
#       BCR == "76" ~ "Boreal Shield West",
#       BCR == "77" ~ "Boreal Shield East",
#       TRUE ~ NA_character_
#     )
#   )
# 
# worst_covariates %>%
#   arrange(BCR, desc(abs(SMD_paired_vs_landscape))) %>%
#   as.data.frame()
# 
# 
# # Key takeaways:
# # In BCR 13, paired and all surveys are highly representative of landscape
# # In BCR 12, the entire survey network is highly road-biased
#       # the "full" survey network is less road-biased than the paired sampling
# 
# 
# # ==============================================================================
# # CHANGE IN HABITAT COMPOSITION BETWEEN OBBA2 AND OBBA3
# #
# # Question:
# # Did habitat at the paired survey locations change differently from habitat
# # across the BCR landscape as a whole?
# # ==============================================================================
# 
# 
# # ------------------------------------------------------------------------------
# # 1. Mean habitat values at paired sites in each atlas
# # ------------------------------------------------------------------------------
# 
# paired_atlas_means <- paired_long %>%
#   group_by(BCR, Atlas, covariate) %>%
#   summarise(
#     mean_paired = mean(value, na.rm = TRUE),
#     .groups = "drop"
#   ) %>%
#   pivot_wider(
#     names_from = Atlas,
#     values_from = mean_paired,
#     names_prefix = "paired_"
#   ) %>%
#   mutate(
#     # Positive = habitat increased at paired sites
#     # Negative = habitat decreased at paired sites
#     change_paired = paired_OBBA3 - paired_OBBA2
#   )
# 
# 
# # ------------------------------------------------------------------------------
# # 2. Mean habitat values across the landscape in each atlas
# # ------------------------------------------------------------------------------
# 
# landscape_atlas_means <- landscape_long %>%
#   group_by(BCR, Atlas, covariate) %>%
#   summarise(
#     mean_landscape = mean(value, na.rm = TRUE),
#     .groups = "drop"
#   ) %>%
#   pivot_wider(
#     names_from = Atlas,
#     values_from = mean_landscape,
#     names_prefix = "landscape_"
#   ) %>%
#   mutate(
#     # Positive = habitat increased across the BCR
#     # Negative = habitat decreased across the BCR
#     change_landscape = landscape_OBBA3 - landscape_OBBA2
#   )
# 
# 
# # ------------------------------------------------------------------------------
# # 3. Compare habitat change at paired sites with habitat change across landscape
# # ------------------------------------------------------------------------------
# 
# habitat_change_comparison <- paired_atlas_means %>%
#   left_join(
#     landscape_atlas_means,
#     by = c("BCR", "covariate")
#   ) %>%
#   mutate(
#     # Positive:
#     # habitat increased MORE (or declined LESS) at paired sites than landscape
#     #
#     # Negative:
#     # habitat increased LESS (or declined MORE) at paired sites than landscape
#     differential_change =
#       change_paired - change_landscape
#   ) %>%
#   arrange(
#     BCR,
#     desc(abs(differential_change))
#   )
# 
# 
# print(habitat_change_comparison, n = Inf)
