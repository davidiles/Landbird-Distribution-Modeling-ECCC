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

min_detections <- 50
min_squares    <- 10

# Radius used to define the footprint shared by both atlases.
shared_radius_km <- 0.100

# Negative-binomial switch thresholds (see choose_error_family()). Evaluated on
# the point counts this script actually fits, so it may differ from script 07.
nb_switch_count <- 50
nb_switch_min_n <- 50

# Timing smooth priors.
# The paired analysis mirrors the full model's two detectability smooths: a
# time-of-day smooth AND a day-of-year smooth, each with the range held fixed and
# only sigma estimated (see fit_inla_shared_footprint_change()).
fixed_TOD_range <- 8                 # hours; held fixed, no prior
prior_TOD_sigma <- c(2, 0.1)         # P(sigma > 2) = 0.1 -> prior median ~0.6

fixed_DOY_range <- 30                # days; held fixed, no prior
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

if (!dir.exists(species_dir)) {
  stop("Per-species data directory not found: ", species_dir,
       "\nRerun 06_filter_and_finalize_surveys.R.")
}

# ============================================================
# 4. Select species to model
# ============================================================
# Same selection rule as script 07, so the two analyses cover the same species.

species_sel <- select_modelable_species(
  species_detection_summaries %>% subset(sp_english %in% c("Purple Martin","Barn Swallow","Wild Turkey","Sandhill Crane","Palm Warbler","Olive-sided Flycatcher","Bald Eagle")),
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
    sp_dat_shared$count,
    nb_switch_count = nb_switch_count,
    nb_switch_min_n = nb_switch_min_n
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
