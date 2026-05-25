# ============================================================
# Paired analysis
# ============================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(sf)
  library(ggplot2)
  library(INLA)
  library(inlabru)
  library(fmesher)
  library(here)
  library(mgcv)
})

# ============================================================
# 1. Paths, utilities, and global configuration
# ============================================================

source(here::here("R", "00_config_paths.R"))
source(file.path(paths$functions, "inla_model_utils_testing.R"))

model_name <- "PC_ARU_only"

# INLA approximation settings used inside fit_inla_multi_atlas().
int_strategy <- "ccd"
strategy <- "laplace"

# ============================================================
# 2. Input/output locations
# ============================================================

in_file <- file.path(paths$data_clean, "birds", "data_ready_for_analysis.rds")

if (!file.exists(in_file)) {
  stop(
    "Cannot find input at: ", in_file,
    "\nHave you run 06_filter_and_finalize_surveys.R?"
  )
}

out_dir <- paths$model_output

summary_dir <- file.path(out_dir, paste0("summaries_", model_name))
data_used_dir <- file.path(out_dir, paste0("data_used_", model_name))

purrr::walk(
  c(out_dir, summary_dir, data_used_dir),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
)

model_summaries_path <- file.path(summary_dir, "model_summaries.rds")
model_summaries <- load_or_empty_list(model_summaries_path)

# ============================================================
# 3. Load finalized data
# ============================================================

dat <- readRDS(in_file)

all_surveys <- dat$all_surveys
counts <- dat$counts
grid_OBBA2 <- dat$grid_OBBA2
grid_OBBA3 <- dat$grid_OBBA3
study_boundary <- dat$study_boundary %>% sf::st_as_sf()
species_to_model <- dat$species_to_model
safe_dates_breeding <- dat$safe_dates_breeding

# ----------------------------------------------------------
# Conduct analysis only in shared survey footprint
# ----------------------------------------------------------
species_run <- species_to_model

paired_summaries <- list()

if (file.exists(file.path(summary_dir,"paired_summaries.rds"))){
  paired_summaries <- readRDS(file.path(summary_dir,"paired_summaries.rds"))
}
for (i in seq_len(nrow(species_run))) {
  
  # ----------------------------------------------------------
  # 5.1 Species identifiers and output paths
  # ----------------------------------------------------------
  
  sp_name <- species_run$english_name[i]
  sp_code <- as.character(species_run$species_id[i])
  sp_file <- sp_filename(sp_name)
  
  if (!(sp_name %in% names(model_summaries))) next
  
  message(
    "\n====================\n",
    i, "/", nrow(species_run), ": ", sp_name,
    " (species_id = ", sp_code, ")\n",
    "===================="
  )
  
  dat_path <- file.path(data_used_dir, paste0(sp_file, "_1km.rds"))
  
  if (!file.exists(dat_path)) next
  
  sp_dat_used <- readRDS(dat_path)
  error_family <- sp_dat_used$error_family
  
  sp_dat <- sp_dat_used$sp_dat
  
  # 50 m buffer
  shared_radius_km = 0.05
  
  summary_name <- paste0("shared_radius_km_",shared_radius_km)
  
  # Skip if already run
  if (length(paired_summaries[[sp_name]])>0) next
  
  sp_dat_shared <- make_shared_footprint_dataset(
    dat = sp_dat %>% subset(Survey_Type %in% c("Point_Count","ARU")),
    atlas_col = "Atlas",
    atlas_levels = c("OBBA2", "OBBA3"),
    buffer_km = shared_radius_km
  ) 
  
  sp_dat_shared <- sp_dat_shared |>
    dplyr::mutate(
      region_factor = factor(Biol_Region)
    )
  
  mod_shared <- fit_inla_shared_footprint_change(
    
    sp_dat_shared = sp_dat_shared,
    
    region_col = "Biol_Region",
    
    family = error_family,
    
    # Timing smooth priors
    prior_HSS_range = c(5, 0.9),
    prior_HSS_sigma = c(3, 0.1),
    
    prior_DOY_range = c(7, 0.9),
    prior_DOY_sigma = c(3, 0.1),
    
    # Random-effect priors
    square_pcprec = c(1, 0.1),
    
    # INLA settings
    int_strategy = int_strategy,
    strategy = strategy,
    
    waic = TRUE,
    cpo = FALSE
  )
  
  summary(mod_shared)
  
  # Summarize model estimates
  shared_change_summary <- summarize_shared_footprint_change(mod_shared) %>%
    mutate(sp_name = sp_name,
           shared_radius_km = shared_radius_km
    ) %>%
    relocate(sp_name,shared_radius_km)
  
  # Summarize data
  shared_data_summary <- sp_dat_shared %>%
    as.data.frame() %>%
    mutate(sp_name = sp_name,
           shared_radius_km = shared_radius_km) %>%
    group_by(sp_name,shared_radius_km,Biol_Region,Atlas) %>%
    summarize(
      n_svy = n(),
      
      mean_count = mean(count),
      PObs = mean(count>0),
      
      mean_DOY = mean(DayOfYear),
      min_DOY = min(DayOfYear),
      max_DOY = max(DayOfYear),
      
      mean_HSS = mean(Hours_Since_Sunrise),
      min_HSS = min(Hours_Since_Sunrise),
      max_HSS = max(Hours_Since_Sunrise)
      
    )
  
  data_to_save <- sp_dat_shared %>%
    dplyr::select(Date_Time,Survey_Type,count,Hours_Since_Sunrise,DayOfYear,Atlas,square_id,Biol_Region)
  
  paired_summaries[[sp_name]][["shared_data_summary"]]   <- shared_data_summary
  paired_summaries[[sp_name]][["shared_change_summary"]] <- shared_change_summary
  paired_summaries[[sp_name]][["shared_data"]] <- data_to_save
  
  saveRDS(paired_summaries, file = file.path(summary_dir,"paired_summaries.rds"))

}
