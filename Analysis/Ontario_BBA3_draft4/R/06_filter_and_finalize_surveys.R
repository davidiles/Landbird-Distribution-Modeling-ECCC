# ============================================================
# 06_filter_and_finalize_surveys.R
#
# Build the analysis-ready survey, count, species-summary, prediction-grid and
# per-species products used by the downstream atlas models (07 and 07b).
#
# Outputs
#   data_clean/birds/data_ready_for_analysis.rds
#       study_boundary, bcr_sf, all_surveys, counts, grid_OBBA2, grid_OBBA3,
#       species_detection_summaries, checklist_candidates, safe_dates_breeding,
#       hex_grid_25km
#
#   data_clean/birds/species_data/<species>_sp_dat.rds        (one per species)
#       Lightweight per-species record: survey_id + count for every survey inside
#       the species' safe-date window (ALL survey types; 07/07b narrow on load),
#       plus pred_doy and sp_safe_dates.
#
# Important: obs_idx links filtered surveys back to rows of count_matrix. Keep
# that alignment unchanged when editing this script.
# ============================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(here)
  library(readxl)
  library(ggplot2)
})

source(here::here("R", "00_config_paths.R"))
source(file.path(paths$functions, "survey_processing_utils.R"))
source(file.path(paths$functions, "covariate_processing_utils.R"))
source(file.path(paths$functions, "inla_model_utils.R"))  # sp_data_path(), make_sp_dat_record(), save_atomic()

in_file     <- file.path(paths$data_clean, "birds", "analysis_data_covariates.rds")
out_file    <- file.path(paths$data_clean, "birds", "data_ready_for_analysis.rds")
species_dir <- file.path(paths$data_clean, "birds", "species_data")

dir.create(file.path(paths$data_clean, "birds"), recursive = TRUE, showWarnings = FALSE)
dir.create(species_dir, recursive = TRUE, showWarnings = FALSE)

safe_date_path <- file.path(paths$data, "Bird_Data_Raw", "OBBA3",
                            "OnAtlasSafeDates_update_for_NatureCounts_2026-07-20.xlsx")
Biol_Regions_path <- file.path(paths$data, "Spatial", "Ontario_Atlas_biol_regions",
                               "Atlas_biol_regions.shp")

# ------------------------------------------------------------
# Analysis settings
# ------------------------------------------------------------

safe_levels  <- c(1, 2)   # safe-date levels to include (core breeding + shoulder)
min_safe_doy <- 1
max_safe_doy <- 365

make_diagnostic_plots <- TRUE   # FALSE to run non-interactively without plots
rebuild_species_data  <- TRUE   # overwrite existing per-species files?

# ------------------------------------------------------------
# Local helper: derived habitat covariates
# ------------------------------------------------------------
# Several covariates represent different vegetation communities, or are
# structurally distinct across BCRs / north-south. Derived covariates let those
# effects be modelled separately.
add_derived_covariates <- function(df,
                                   urban_var     = "LCC_17",
                                   grassland_var = "LCC_10",
                                   shrubland_var = "LCC_8") {

  if (!("BCR" %in% names(df))) {
    stop("BCR missing; required for derived covariates.")
  }

  get_or_zero <- function(nm) if (nm %in% names(df)) df[[nm]] else rep(0, nrow(df))

  df %>%
    mutate(
      Urban   = as.integer(get_or_zero(urban_var) > 0.5),
      On_Road = as.integer(get_or_zero("road") > 0 & Urban == 0),

      Grassland_North = if_else(BCR %in% c("76", "77", "74"),   get_or_zero(grassland_var), 0),
      Grassland_South = if_else(!(BCR %in% c("76", "77", "74")), get_or_zero(grassland_var), 0),

      Shrubland_South = if_else(BCR == "13", get_or_zero(shrubland_var), 0),
      Shrubland_North = if_else(BCR != "13", get_or_zero(shrubland_var), 0),

      River_Lg_North  = if_else(BCR == "74", get_or_zero("rivers_large"), 0),
      River_Lg_South  = if_else(BCR != "74", get_or_zero("rivers_large"), 0),

      River_Sm_North  = if_else(BCR == "74", get_or_zero("rivers_small"), 0),
      River_Sm_South  = if_else(BCR != "74", get_or_zero("rivers_small"), 0)
    ) %>%
    dplyr::rename(
      "GreatLakes"       = "great_lakes",
      "HudsonBayCoast"   = "coastline",
      "Lake_Lg"          = "lakes_large",
      "Lake_Sm"          = "lakes_small",
      "River_Lg"         = "rivers_large",
      "River_Sm"         = "rivers_small",
      "ForestNeedleleaf" = "LCC_1",
      "ForestBroadleaf"  = "LCC_5",
      "ForestMixed"      = "LCC_6",
      "Wetland"          = "LCC_14",
      "Cropland"         = "LCC_15"
    )
}

# ============================================================
# 1. Load study area and covariate-enriched survey objects
# ============================================================

study_area     <- readRDS(file.path(paths$data_clean, "spatial", "study_area.rds"))
bcr            <- study_area$bcr
study_boundary <- study_area$boundary

# Retain BCRs 12, 13, 74, 76 and 77 as separate regions. BCR 76 (Boreal Shield
# West) and 77 (Boreal Shield East) were previously merged into "76_77"; they are
# now kept apart, so every downstream product is reported for five regions. The
# group_by()/st_union() below only dissolves the multiple input features that
# belong to each individual BCR.
bcr_split <- bcr %>%
  filter(bcr %in% c(12, 13, 74, 76, 77)) %>%
  mutate(
    BCR = case_when(
      bcr == 12 ~ "12", bcr == 13 ~ "13", bcr == 74 ~ "74",
      bcr == 76 ~ "76", bcr == 77 ~ "77"
    ),
    BCR_Label = case_when(
      bcr == 12 ~ "Temperate Mixed",
      bcr == 13 ~ "Lower Great Lakes / St. Lawrence Plain",
      bcr == 74 ~ "Hudson Plains",
      bcr == 76 ~ "Boreal Shield West",
      bcr == 77 ~ "Boreal Shield East"
    )
  ) %>%
  group_by(BCR, BCR_Label) %>%
  summarise(geometry = st_union(geometry), .groups = "drop")

dat <- readRDS(in_file)

all_surveys  <- dat$all_surveys_with_covs
count_matrix <- dat$count_matrix
grid_OBBA2   <- dat$grid_OBBA2
grid_OBBA3   <- dat$grid_OBBA3

stopifnot(inherits(all_surveys, "sf"))
stopifnot(inherits(grid_OBBA2, "sf"))
stopifnot(inherits(grid_OBBA3, "sf"))

all_species_unique <- dat$all_species %>%
  distinct(species_id, .keep_all = TRUE) %>%
  transmute(species_id = as.character(species_id), english_name)

# ============================================================
# 2. Validate required fields and build a stable count table
# ============================================================

if (!("obs_idx" %in% names(all_surveys))) {
  all_surveys <- all_surveys %>% mutate(obs_idx = row_number(), .before = 1)
} else {
  all_surveys <- all_surveys %>% arrange(obs_idx)
}

if (!("DayOfYear" %in% names(all_surveys))) {
  all_surveys <- all_surveys %>% mutate(DayOfYear = yday(Date_Time))
}

# Align count_matrix to survey order BY KEY, then carry survey_id into the count
# table so every downstream subset joins on the key rather than trusting row
# position. (02 and 05 already enforce the ordering; this keeps 06 self-contained.)
stopifnot(
  setequal(rownames(count_matrix), all_surveys$survey_id),
  anyDuplicated(all_surveys$survey_id) == 0
)
count_matrix <- count_matrix[all_surveys$survey_id, , drop = FALSE]

counts <- as_counts_tbl(count_matrix)
counts$survey_id <- rownames(count_matrix)

required_cols <- c("survey_id", "Survey_Duration_Minutes", "DayOfYear",
                   "Date_Time", "Atlas", "Survey_Type")
missing_cols <- setdiff(required_cols, names(all_surveys))
if (length(missing_cols) > 0) {
  stop("Missing required survey columns: ", paste(missing_cols, collapse = ", "))
}

# survey_id is the natural key used to rejoin per-species counts in 07 / 07b.
if (anyDuplicated(all_surveys$survey_id) > 0) {
  stop("survey_id is not unique in all_surveys; the per-species join in 07 will fail.")
}

# ============================================================
# 3. Survey filter
# ============================================================

# --- Tunable constants (kept together so the recode and the stationary rule can
#     never drift apart and reopen a dead zone). ---
STATIONARY_TOL_M    <- 25    # <= this distance counts as stationary
MAX_TRANSECT_M      <- 500   # upper bound on transect length
MAX_SPEED_M_PER_MIN <- 50    # upper bound on travel speed
MIN_CHECKLIST_MIN   <- 3     # minimum duration for SC and LT
MIN_DURATION_MIN    <- 0.5   # general lower bound
MAX_DURATION_MIN    <- 30    # general upper bound
PC_DURATION_MIN     <- 5     # point counts must be exactly this
ARU_DURATIONS_MIN   <- c(1, 3, 5)
MAX_YEAR            <- 2025

CHECKLIST_TYPES  <- c("Breeding Bird Atlas", "Linear transect")
STATIONARY_TYPES <- c("ARU", "Breeding Bird Atlas", "Point_Count")
RETAINED_TYPES   <- c("Point_Count", "ARU", "Breeding Bird Atlas", "Linear transect")

surveys_f <- all_surveys %>%
  mutate(
    Survey_Type_raw = Survey_Type,   # provenance: what the provider called it

    # Negative distances are impossible: treat as unknown rather than letting
    # them slip past a `<= tolerance` test as if they were stationary.
    Distance_Traveled_m = if_else(
      !is.na(Distance_Traveled_m) & Distance_Traveled_m < 0,
      NA_real_, Distance_Traveled_m
    ),

    walking_speed_metres_per_min = Distance_Traveled_m / Survey_Duration_Minutes,

    # --- Standardize survey types before filtering ---
    # Within the checklist family, distance decides the protocol, not the
    # incoming label. Unknown distance keeps its raw label (a checklist is
    # assumed stationary; a transect with no distance is dropped downstream).
    Survey_Type = case_when(
      !Survey_Type %in% CHECKLIST_TYPES      ~ Survey_Type,
      is.na(Distance_Traveled_m)             ~ Survey_Type,
      Distance_Traveled_m > STATIONARY_TOL_M ~ "Linear transect",
      TRUE                                   ~ "Breeding Bird Atlas"
    )
  ) %>%
  filter(
    # Retained survey types (first, so the cascade reads correctly).
    Survey_Type %in% RETAINED_TYPES,

    # NA effort is missing-by-design for several providers, NOT a "Special"
    # protocol designation, so retain it.
    is.na(EffortMeasurement1) | EffortMeasurement1 != "Special",

    lubridate::year(Date_Time) <= MAX_YEAR,

    # NOTE: drops NA durations silently (~30k checklists), pending the
    # EffortUnits* investigation into whether duration lives elsewhere.
    between(Survey_Duration_Minutes, MIN_DURATION_MIN, MAX_DURATION_MIN),

    !is.na(rivers_large),

    # Protocol-specific duration rules.
    Survey_Type != "Point_Count" | Survey_Duration_Minutes == PC_DURATION_MIN,
    Survey_Type != "ARU"         | Survey_Duration_Minutes %in% ARU_DURATIONS_MIN,

    # Stationary types must not have travelled beyond tolerance.
    !Survey_Type %in% STATIONARY_TYPES |
      is.na(Distance_Traveled_m) |
      Distance_Traveled_m <= STATIONARY_TOL_M,

    # Linear-transect distance and speed restrictions.
    Survey_Type != "Linear transect" |
      (!is.na(Distance_Traveled_m) &
         Distance_Traveled_m > STATIONARY_TOL_M &
         Distance_Traveled_m <= MAX_TRANSECT_M &
         Survey_Duration_Minutes >= MIN_CHECKLIST_MIN &
         walking_speed_metres_per_min <= MAX_SPEED_M_PER_MIN),

    # Stationary checklist duration restriction.
    Survey_Type != "Breeding Bird Atlas" | Survey_Duration_Minutes >= MIN_CHECKLIST_MIN
  ) %>%
  arrange(obs_idx) %>%
  mutate(
    SC        = as.integer(Survey_Type == "Breeding Bird Atlas"),
    LT        = as.integer(Survey_Type == "Linear transect"),
    ARU       = as.integer(Survey_Type == "ARU"),
    Checklist = as.integer(Survey_Type %in% CHECKLIST_TYPES),

    # Keep the character key alongside the integer code. Integer codes are
    # assigned post-filter, so they shift whenever the filter changes; join on
    # the key, never a cached index.
    square_atlas_key = paste(square_id, Atlas, sep = "-"),
    square_atlas     = as.integer(factor(square_atlas_key))
  )

print(table(as.data.frame(surveys_f)[, c("Atlas", "Survey_Type")]))

# ------------------------------------------------------------
# Subset counts to the retained survey rows BY KEY (survey_id), never by row
# position. The join is 1:1 (survey_id is unique on both sides).
# ------------------------------------------------------------

counts_f <- tibble(obs_idx = surveys_f$obs_idx, survey_id = surveys_f$survey_id) %>%
  left_join(counts %>% select(-obs_idx), by = "survey_id")

stopifnot(
  nrow(counts_f) == nrow(surveys_f),
  identical(counts_f$survey_id, surveys_f$survey_id),
  !anyNA(counts_f$survey_id)
)

# ============================================================
# 4. Solar time, site IDs, and region covariates
# ============================================================

# Corrects time zones, adds sunrise/sunset times, and adds a time-of-day
# covariate (Hours_After_Reference, referenced to 3 h before sunrise).
surveys_f <- add_solar_time_covariates_reference(surveys_sf = surveys_f, reference_hours = -3)

# Site (site-year) identifiers for the site-level iid random effect.
surveys_f <- add_site_ids(surveys_f, tolerance_m = 100)
surveys_f$site_id <- surveys_f$site_year
summarise_site_structure(surveys_f)

# --- BCR ---
surveys_f$BCR <- assign_poly_id(surveys_f, bcr_split, id_col = "BCR", nearest_fallback = TRUE)

grid_bcr <- assign_poly_id(grid_OBBA2, bcr_split %>% st_transform(st_crs(grid_OBBA2)),
                           id_col = "BCR", nearest_fallback = TRUE)
grid_OBBA2$BCR <- grid_bcr
grid_OBBA3$BCR <- grid_bcr

# --- Atlas biological regions (used by safe-date summaries) ---
if (!file.exists(Biol_Regions_path)) {
  stop("Biological regions shapefile not found at: ", Biol_Regions_path)
}

Biol_Regions_ON <- st_read(Biol_Regions_path, quiet = TRUE) %>%
  st_make_valid() %>%
  dplyr::rename(Biol_Region = Biol_Regio)

surveys_f$Biol_Region <- assign_poly_id(surveys_f, Biol_Regions_ON,
                                        id_col = "Biol_Region", nearest_fallback = TRUE)

grid_Biol_Region <- assign_poly_id(grid_OBBA2,
                                   Biol_Regions_ON %>% st_transform(st_crs(grid_OBBA2)),
                                   id_col = "Biol_Region", nearest_fallback = TRUE)
grid_OBBA2$Biol_Region <- grid_Biol_Region
grid_OBBA3$Biol_Region <- grid_Biol_Region

# --- Survey-level model helper variables ---
surveys_f <- surveys_f %>%
  mutate(
    Atlas3        = if_else(Atlas == "OBBA3", 1L, 0L),
    Atlas3_c      = Atlas3 - 0.5,
    days_rescaled = DayOfYear - 166,
    BCR_idx       = as.integer(factor(BCR))
  )

# --- Derived habitat covariates ---
surveys_f  <- add_derived_covariates(surveys_f)
grid_OBBA2 <- add_derived_covariates(grid_OBBA2)
grid_OBBA3 <- add_derived_covariates(grid_OBBA3)

# ============================================================
# 5. Diagnostics
# ============================================================

surveys_f %>%
  st_drop_geometry() %>%
  group_by(Atlas, Survey_Type) %>%
  summarize(n = n(), .groups = "drop") %>%
  print()

surveys_f %>%
  st_drop_geometry() %>%
  group_by(Survey_Type) %>%
  summarise(
    across(Hours_After_Reference,
           list(q0     = ~ quantile(.x, 0,     na.rm = TRUE),
                q0.025 = ~ quantile(.x, 0.025, na.rm = TRUE),
                q0.5   = ~ quantile(.x, 0.5,   na.rm = TRUE),
                q0.975 = ~ quantile(.x, 0.975, na.rm = TRUE),
                q1     = ~ quantile(.x, 1,     na.rm = TRUE)),
           .names = "{.col}_{.fn}"),
    .groups = "drop"
  ) %>%
  print()

if (make_diagnostic_plots) {

  print(
    ggplot(surveys_f) +
      geom_histogram(aes(x = hour(Date_Time_Local)), bins = 50) +
      facet_grid(Survey_Type ~ Atlas) + theme_bw() + ggtitle("Local hour of day")
  )

  print(
    ggplot(surveys_f %>% filter(Survey_Duration_Minutes %in% c(1, 3, 5) &
                                  Survey_Type %in% c("Point_Count", "ARU"))) +
      geom_histogram(aes(x = Survey_Duration_Minutes)) +
      facet_grid(Atlas ~ Survey_Type) + theme_bw() +
      ggtitle("Survey duration (ARU & point counts only)")
  )

  print(
    ggplot(surveys_f) +
      geom_histogram(aes(x = Hours_After_Reference), bins = 24 * 2) +
      facet_grid(Survey_Type ~ .) +
      geom_vline(xintercept = 3, linetype = 2) + theme_bw()
  )

  print(
    ggplot() +
      geom_sf(data = study_boundary) +
      geom_sf(data = surveys_f, size = 0.1) +
      facet_grid(Atlas ~ Survey_Type) + theme_bw()
  )
}

# ============================================================
# 6. Species safe-date windows
# ============================================================

Biol_Regions_df <- Biol_Regions_ON %>%
  st_drop_geometry() %>%
  dplyr::select(Biol_Region, ECOZONE_NA)

safe_dates <- readxl::read_xlsx(safe_date_path) %>%
  dplyr::rename(
    sp_english  = english_name,
    Biol_Region = biol_region,
    start_doy   = `start_dt (julian)`,
    end_doy     = `end_dt (julian)`
  )

# Collapse selected levels to one window per species x region.
safe_date_windows <- safe_dates %>%
  dplyr::filter(level %in% safe_levels) %>%
  dplyr::group_by(sp_english, Biol_Region) %>%
  dplyr::summarise(
    start_doy = min(start_doy, na.rm = TRUE),
    end_doy   = max(end_doy,   na.rm = TRUE),
    .groups   = "drop"
  ) %>%
  dplyr::mutate(
    start_doy = pmax(start_doy, min_safe_doy),
    end_doy   = pmin(end_doy,   max_safe_doy)
  ) %>%
  dplyr::filter(start_doy < end_doy)

safe_dates_breeding <- safe_dates %>%
  dplyr::distinct(sp_english, Biol_Region) %>%
  dplyr::left_join(safe_date_windows, by = c("sp_english", "Biol_Region")) %>%
  dplyr::mutate(midpoint = (start_doy + end_doy) / 2) %>%
  dplyr::select(sp_english, Biol_Region, start_doy, end_doy, midpoint) %>%
  dplyr::full_join(Biol_Regions_df, by = "Biol_Region") %>%
  dplyr::left_join(all_species_unique, by = c("sp_english" = "english_name")) %>%
  dplyr::relocate(sp_english, species_id, Biol_Region, ECOZONE_NA, start_doy, end_doy) %>%
  arrange(sp_english)

# ============================================================
# 7. Per-species datasets and detection summaries
# ============================================================
#   (a) species_detection_summaries -- detections AND total surveys per
#       Atlas x Survey_Type, counted over surveys INSIDE the safe-date window
#       (i.e. exactly the surveys stored in <species>_sp_dat.rds, which 07/07b
#       load and model).
#   (b) <species>_sp_dat.rds -- survey_id + count for surveys INSIDE the
#       safe-date window, all survey types retained.

survey_df <- surveys_f %>% st_drop_geometry()
species_detection_summaries <- data.frame()
species_with_no_detections  <- character(0)
species_with_no_safe_dates  <- character(0)
species_with_no_surveys     <- character(0)

species_ids <- unique(safe_dates_breeding$species_id)
species_ids <- species_ids[!is.na(species_ids)]

for (sp_id in species_ids) {
  sp_safe_dates <- safe_dates_breeding %>% filter(species_id == sp_id)
  sp_name       <- sp_safe_dates$sp_english[1]
  sp_safe_dates <- na.omit(sp_safe_dates)

  if (nrow(sp_safe_dates) == 0) {
    message("!!!! ", sp_name, " (species_id = ", sp_id,
            ") has no safe dates listed; skipping")
    species_with_no_safe_dates <- c(species_with_no_safe_dates, sp_name)
    next
  }
  if (!(as.character(sp_id) %in% colnames(counts_f))) {
    message("!!!! ", sp_name, " (species_id = ", sp_id,
            ") was not observed in atlas dataset; skipping")
    species_with_no_detections <- c(species_with_no_detections, sp_name)
    next
  }
  message(sp_name)
  count_vec <- counts_f[[as.character(sp_id)]]

  # ---- (b) per-species record ----
  # Built (or reloaded) FIRST, because the detection summary now derives from it.
  sp_path <- sp_data_path(species_dir, sp_name)
  if (file.exists(sp_path) && !rebuild_species_data) {
    sp_record <- readRDS(sp_path)
  } else {
    sp_record <- make_sp_dat_record(
      surveys_f     = surveys_f,
      count_vec     = count_vec,
      sp_english    = sp_name,
      species_id    = sp_id,
      sp_safe_dates = sp_safe_dates
    )
    if (is.null(sp_record)) {
      message("     no surveys remain inside the safe-date window; no file written")
      species_with_no_surveys <- c(species_with_no_surveys, sp_name)
      next
    }
    save_atomic(sp_record, sp_path)
  }

  # ---- (a) detection summary (safe-date window only) ----
  sp_summary <- sp_record$survey_counts %>%
    dplyr::left_join(
      dplyr::select(survey_df, survey_id, Atlas, square_id, Survey_Type),
      by = "survey_id"
    ) %>%
    dplyr::group_by(Atlas, Survey_Type) %>%
    dplyr::summarize(
      n_det         = sum(count > 0),
      n_surveys     = dplyr::n(),
      n_sq          = dplyr::n_distinct(square_id[count > 0]),
      n_sq_surveyed = dplyr::n_distinct(square_id),
      .groups       = "drop"
    ) %>%
    dplyr::mutate(sp_english = sp_name, species_id = sp_id) %>%
    dplyr::relocate(sp_english, species_id)
  species_detection_summaries <- bind_rows(species_detection_summaries, sp_summary)
}

message(
  "\nPer-species data written to: ", species_dir,
  "\n  files created:                ", length(list.files(species_dir, pattern = "_sp_dat\\.rds$")),
  "\n  species with no safe dates:   ", length(species_with_no_safe_dates),
  "\n  species never detected:       ", length(species_with_no_detections),
  "\n  species with no safe surveys: ", length(species_with_no_surveys)
)
print(head(species_detection_summaries, 20))

# ------------------------------------------------------------
# Checklist triage: species disproportionately detected in checklists
# ------------------------------------------------------------

checklist_types  <- c("Linear transect", "Breeding Bird Atlas")
structured_types <- c("ARU", "Point_Count")
min_det <- 50   # per class; below this the ratio is mostly noise

checklist_triage <- species_detection_summaries %>%
  filter(Atlas == "OBBA3") %>%                        # checklist types are OBBA3-only
  mutate(class = case_when(
    Survey_Type %in% checklist_types  ~ "chk",
    Survey_Type %in% structured_types ~ "str"
  )) %>%
  filter(!is.na(class)) %>%
  group_by(sp_english, species_id, class) %>%
  summarize(det = sum(n_det), srv = sum(n_surveys), .groups = "drop") %>%
  pivot_wider(names_from = class, values_from = c(det, srv)) %>%
  mutate(across(c(det_chk, det_str), ~ replace_na(., 0L))) %>%
  mutate(
    rate_chk  = det_chk / srv_chk,
    rate_str  = det_str / srv_str,
    # +0.5 continuity correction so zero-detection classes stay finite.
    log_rr    = log((det_chk + 0.5) / srv_chk) - log((det_str + 0.5) / srv_str),
    se_log_rr = sqrt(1 / (det_chk + 0.5) + 1 / (det_str + 0.5)),
    rr        = exp(log_rr),
    rr_lo     = exp(log_rr - 1.96 * se_log_rr),
    rr_hi     = exp(log_rr + 1.96 * se_log_rr)
  )

med_log_rr <- median(checklist_triage$log_rr[checklist_triage$det_chk >= min_det &
                                               checklist_triage$det_str >= min_det], na.rm = TRUE)
checklist_triage <- checklist_triage %>% mutate(rel_log_rr = log_rr - med_log_rr)

checklist_candidates <- checklist_triage %>%
  filter(det_chk >= min_det) %>%
  mutate(
    over_rep = rr_lo > 1,        # checklist rate reliably exceeds structured
    thin_str = det_str < 250,    # structured surveys alone can't carry the fit
    priority = over_rep & thin_str
  ) %>%
  arrange(desc(rel_log_rr)) %>%
  select(sp_english, species_id, det_str, det_chk,
         rate_str, rate_chk, rr, rr_lo, rr_hi, rel_log_rr, priority)
print(head(as.data.frame(checklist_candidates), 50))

# ============================================================
# 8. Build 25-km hex grid and pixel-to-hex lookup
# ============================================================

set.seed(123)
hex_grid_25km <- make_hex_grid(study_boundary, width_km = 25)

grid_OBBA2 <- grid_OBBA2 %>% mutate(Atlas = "OBBA2")
grid_OBBA3 <- grid_OBBA3 %>% mutate(Atlas = "OBBA3")

hex_idx <- build_pixel_polygon_index(
  grid_sf     = grid_OBBA2,
  polygons_sf = hex_grid_25km,
  poly_id_col = "hex_id",
  join        = "within"
)

prediction_hex_lookup <- tibble(
  pixel_id = grid_OBBA2$pixel_id,
  hex_id   = hex_idx$pix_poly_id
) %>%
  filter(!is.na(hex_id))

grid_OBBA2 <- grid_OBBA2 %>% left_join(prediction_hex_lookup, by = "pixel_id") %>% na.omit()
grid_OBBA3 <- grid_OBBA3 %>% left_join(prediction_hex_lookup, by = "pixel_id") %>% na.omit()

# ============================================================
# 9. Save analysis-ready objects
# ============================================================

save_atomic(
  list(
    study_boundary = study_boundary,
    bcr_sf         = bcr_split,

    all_surveys = surveys_f,
    counts      = counts_f,
    grid_OBBA2  = grid_OBBA2,
    grid_OBBA3  = grid_OBBA3,

    species_detection_summaries = species_detection_summaries,
    checklist_candidates        = checklist_candidates,
    safe_dates_breeding         = safe_dates_breeding,
    hex_grid_25km               = hex_grid_25km,

    date_created = Sys.time()
  ),
  out_file
)

message("06_filter_and_finalize_surveys.R complete")
