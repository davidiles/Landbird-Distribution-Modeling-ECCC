# ============================================================
# 09_fit_models_and_predict_PC_ARU_chunked.R
#
# Purpose:
#   Fit INLA/inlabru models for selected species and generate
#   predictions on the OBBA2/OBBA3 grids in memory-safe chunks.
#
# Assumptions from script 07:
#   - safe_dates_bcr is saved in data_ready_for_analysis.rds
#   - grid_OBBA2 / grid_OBBA3 already contain:
#       * pixel_id
#       * Atlas labels
#       * hex_id (or a prediction_chunk_lookup is saved separately)
#   - hex_grid is saved
#   - prediction_chunk_lookup maps pixel_id -> hex_id
#
# Key implementation details:
#   - Predictions are generated chunk-by-chunk (here, hex-by-hex)
#   - The SAME prediction seed is used for every chunk so posterior
#     draw k is aligned across chunks
#   - This avoids constructing a full province-wide n_pixel x n_draw
#     object in memory
#
# Inputs:
#   data_clean/birds/data_ready_for_analysis.rds
#
# Outputs (written incrementally per species):
#   data_clean/model_output/models_<model_name>/<sp>_1km.rds
#   data_clean/model_output/predictions_<model_name>/<sp>_1km.rds
#   data_clean/model_output/summaries_<model_name>/model_summaries.rds
#   data_clean/model_output/summaries_<model_name>/change_summaries.rds
#   data_clean/model_output/summaries_<model_name>/HSS_DOY_summaries.rds
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
})

# ------------------------------------------------------------
# Centralized paths
# ------------------------------------------------------------

source(here::here("R", "00_config_paths.R"))
source(file.path(paths$functions, "inla_model_utils2.R"))

# ------------------------------------------------------------
# Config
# ------------------------------------------------------------

model_name <- "PC_ARU_CL_new"
rerun_models <- TRUE
rerun_predictions <- TRUE

n_prediction_draws <- 500
prediction_seed <- 123
chunk_id_col <- "chunk_id"

in_file <- file.path(paths$data_clean, "birds", "data_ready_for_analysis.rds")
if (!file.exists(in_file)) {
  stop("Cannot find input at: ", in_file,
       "\nHave you run 07_filter_and_finalize_surveys.R?")
}

out_dir <- paths$model_output
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, paste0("predictions_", model_name)), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, paste0("summaries_", model_name)), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, paste0("data_used_", model_name)), recursive = TRUE, showWarnings = FALSE)

model_summaries_path  <- file.path(out_dir, paste0("summaries_", model_name), "model_summaries.rds")
change_summaries_path <- file.path(out_dir, paste0("summaries_", model_name), "change_summaries.rds")
hss_doy_path          <- file.path(out_dir, paste0("summaries_", model_name), "HSS_DOY_summaries.rds")

model_summaries   <- load_or_empty_list(model_summaries_path)
change_summaries  <- load_or_empty_list(change_summaries_path)
hss_doy_summaries <- load_or_empty_list(hss_doy_path)

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

summarize_chunk_predictions <- function(mu2, mu3) {
  abs_change <- mu3 - mu2
  list(
    OBBA2 = summarize_posterior(mu2, CI_probs = c(0.05, 0.95), prefix = "OBBA2"),
    OBBA3 = summarize_posterior(mu3, CI_probs = c(0.05, 0.95), prefix = "OBBA3"),
    abs_change = summarize_posterior(abs_change, CI_probs = c(0.05, 0.95), prefix = "abs_change")
  )
}

bind_chunk_summaries <- function(chunk_summaries, atlas = c("OBBA2", "OBBA3", "abs_change")) {
  atlas <- match.arg(atlas)
  out <- purrr::map(chunk_summaries, atlas)
  out <- purrr::compact(out)
  if (length(out) == 0) return(NULL)
  dplyr::bind_rows(out)
}

make_pred_grid <- function(grid_obba2, grid_obba3) {
  bind_rows(
    grid_obba2 %>% mutate(Atlas3 = 0L),
    grid_obba3 %>% mutate(Atlas3 = 1L)
  ) %>%
    mutate(
      Hours_Since_Sunrise = 0,
      days_midpoint = 0,
      Atlas3_c = Atlas3 - 0.5,
      BCR_idx = as.integer(BCR)
    )
}

predict_one_chunk <- function(mod, pred_grid_chunk, pred_formula,
                              n.samples, seed,
                              on_water_col = "on_water") {
  
  preds <- predict_inla(
    mod = mod,
    grid = pred_grid_chunk,
    pred_formula = pred_formula,
    n.samples = n.samples,
    seed = seed
  )
  
  # # Force deterministic zero abundance on open water
  # if (on_water_col %in% names(pred_grid_chunk)) {
  #   water_idx <- which(!is.na(pred_grid_chunk[[on_water_col]]) & pred_grid_chunk[[on_water_col]])
  #   if (length(water_idx) > 0) {
  #     preds$eta[water_idx, ] <- -Inf
  #   }
  # }
  
  idx2 <- which(pred_grid_chunk$Atlas == "OBBA2")
  idx3 <- which(pred_grid_chunk$Atlas == "OBBA3")
  
  eta <- preds$eta
  mu2 <- exp(eta[idx2, , drop = FALSE])
  mu3 <- exp(eta[idx3, , drop = FALSE])
  
  if (nrow(mu2) != nrow(mu3)) {
    stop("Chunk does not contain matched OBBA2/OBBA3 rows. Check chunk construction.")
  }
  
  list(
    eta = eta,
    mu2 = mu2,
    mu3 = mu3,
    summary = summarize_chunk_predictions(mu2 = mu2, mu3 = mu3)
  )
}

# ------------------------------------------------------------
# Load data
# ------------------------------------------------------------

dat <- readRDS(in_file)

all_surveys    <- dat$all_surveys
counts         <- dat$counts
grid_OBBA2     <- dat$grid_OBBA2
grid_OBBA3     <- dat$grid_OBBA3
study_boundary <- dat$study_boundary %>% st_as_sf()
species_to_model <- dat$species_to_model
hex_grid       <- dat$hex_grid

# Safe dates expanded to BCR should already exist from script 07
if ("safe_dates_bcr" %in% names(dat)) {
  safe_dates_bcr <- dat$safe_dates_bcr
} else {
  safe_dates_bcr <- dat$safe_dates %>%
    mutate(
      BCR = case_when(
        ecoregion == "Mixedwood Plains" ~ "13",
        ecoregion == "Hudson Plains"    ~ "7",
        ecoregion == "Boreal Shield"    ~ "8, 12",
        TRUE ~ NA_character_
      )
    ) %>%
    separate_rows(BCR, sep = ",") %>%
    mutate(BCR = as.integer(trimws(BCR))) %>%
    select(sp_english, BCR, start_doy, end_doy) %>%
    arrange(sp_english, BCR) %>%
    mutate(midpoint = (start_doy + end_doy) / 2)
}

# Pixel -> chunk mapping
if ("prediction_chunk_lookup" %in% names(dat)) {
  pred_chunk_lookup <- dat$prediction_chunk_lookup
} else {
  idx <- build_pixel_polygon_index(
    grid_sf = grid_OBBA2,
    polygons_sf = hex_grid,
    poly_id_col = chunk_id_col,
    join = "within"
  )
  
  pred_chunk_lookup <- tibble(
    pixel_id = grid_OBBA2$pixel_id,
    !!chunk_id_col := idx$pix_poly_id
  ) %>%
    filter(!is.na(.data[[chunk_id_col]]))
}

stopifnot(nrow(all_surveys) == nrow(counts))
stopifnot("pixel_id" %in% names(grid_OBBA2), "pixel_id" %in% names(grid_OBBA3))
stopifnot(chunk_id_col %in% names(pred_chunk_lookup))

lookup_cols <- c("pixel_id", chunk_id_col)

# Only join chunk IDs if not already present
if (!(chunk_id_col %in% names(grid_OBBA2))) {
  grid_OBBA2 <- grid_OBBA2 %>% left_join(pred_chunk_lookup[, lookup_cols], by = "pixel_id")
}
if (!(chunk_id_col %in% names(grid_OBBA3))) {
  grid_OBBA3 <- grid_OBBA3 %>% left_join(pred_chunk_lookup[, lookup_cols], by = "pixel_id")
}

# Keep only pixels assigned to a chunk
grid_OBBA2 <- grid_OBBA2 %>% filter(!is.na(.data[[chunk_id_col]]))
grid_OBBA3 <- grid_OBBA3 %>% filter(!is.na(.data[[chunk_id_col]]))

chunk_ids <- sort(unique(grid_OBBA2[[chunk_id_col]]))

# ------------------------------------------------------------
# Main loop
# ------------------------------------------------------------

base_covars <- c(
  "ForestNeedleleaf",
  "ForestBroadleaf",
  "ForestMixed",
  "Wetland",
  "Cropland",
  
  "InsectBroadleaf",
  "InsectNeedleleaf",
  "Urban",
  
  "On_Road",
  "Grassland_BCR7_8",
  "Grassland_BCR12_13",
  "Shrubland_BCR13",
  "Shrubland_BCR7_8_12",
  
  "River_Sm",
  "River_Lg",
  "Lake_Lg",
  "Lake_Sm"
)

min_detections <- 250
min_squares <- 50

set.seed(123)
species_run <- species_to_model %>%
  filter(
    (detections_safe_OBBA2 >= min_detections |
       detections_safe_OBBA3 >= min_detections) &
      (n_squares_safe_OBBA2 >= min_squares |
         n_squares_safe_OBBA3 >= min_squares)
  ) %>%
  na.omit() %>%
  mutate(
    delta_dets_safe = log(detections_safe_OBBA3 / detections_safe_OBBA2),
    delta_squares_safe = log(n_squares_safe_OBBA3 / n_squares_safe_OBBA2)
  ) %>%
  arrange(abs(delta_squares_safe))

species_run <- species_run %>%
  subset(english_name %in% 
           c(#"Belted Kingfisher"
             #"Winter Wren",
             # "Spotted Sandpiper",           # BAM makes different predictions in HBL
             # "Wilson's Warbler",            # eBird, BAM, and Atlas make different predictions in HBL
             # "Yellow-rumped Warbler",       # eBird and BAM make very different predictions in HBL
             # "Wilson's Snipe",              # Atlas and eBird make different predictions in HBL
             # "Solitary Sandpiper",          # Good example of eBird and Atlas showing different hotspot locations; BAM and eBird also make very different predictions
             # "Savannah Sparrow",            # Check for correspondence with eBird in revised Atlas model
             # "Philadelphia Vireo",          # BAM predicts some extreme hotspots in Northern Ontario
             # "Hermit Thrush",               # Atlas shows HBL as hotspot, BAM shows almost none, eBird intermediate
             # "Dark-eyed Junco",             # Atlas shows HBL as hotspot
             #"Connecticut Warbler",
             #"Black-and-white Warbler",
             #"Bobolink",
             #"Palm Warbler",
             #"Sandhill Crane",
             # "Bank Swallow",
             # "Eastern Meadlowlark",
             # "Northern Cardinal",
             # "Olive-sided Flycatcher",
             # "Dark-eyed Junco",
             # "Osprey",
             "Rock Pigeon (Feral Pigeon)"
             # "Swainson's Thrush",
             # "American Crow",
             # "American Goldfinch",
             # "Bald Eagle",
             # "Boreal Chickadee",
             # "Common Yellowthroat"
             ))
print(species_run)
message("Species queued: ", nrow(species_run))

for (i in seq_len(nrow(species_run))) {
  
  # To choose a particular species manually, uncomment and edit:
  # i <- which(species_run$english_name == "Belted Kingfisher")
  
  sp_name <- species_run$english_name[i]
  sp_code <- as.character(species_run$species_id[i])
  sp_file <- sp_filename(sp_name)
  
  message(
    "\n====================\n",
    i, "/", nrow(species_run), ": ", sp_name,
    " (species_id = ", sp_code, ")\n",
    "===================="
  )
  
  model_path <- file.path(out_dir, paste0("models_", model_name), paste0(sp_file, "_1km.rds"))
  pred_path  <- file.path(out_dir, paste0("predictions_", model_name), paste0(sp_file, "_1km.rds"))
  dat_path   <- file.path(out_dir, paste0("data_used_", model_name), paste0(sp_file, "_1km.rds"))
  
  if (!(sp_code %in% names(counts))) {
    message("Skipping (species_id not found in counts columns): ", sp_code)
    next
  }
  
  sp_dat <- all_surveys %>%
    mutate(count = counts[[sp_code]])
  
  # ------------------------------------------------------------
  # Filter data to species-safe dates
  # ------------------------------------------------------------
  
  sp_safe_dates <- safe_dates_bcr %>%
    filter(sp_english == sp_name)
  
  if (nrow(sp_safe_dates) == 0) {
    pred_doy <- NA_real_
    warning(paste0("Species '", sp_name, "' has no BCR safe dates at all."), call. = FALSE)
  } else {
    fallback_start <- max(sp_safe_dates$start_doy, na.rm = TRUE)
    fallback_end   <- min(sp_safe_dates$end_doy,   na.rm = TRUE)
    
    if (fallback_start <= fallback_end) {
      all_bcrs <- sp_dat %>% distinct(BCR)
      safe_bcrs <- sp_safe_dates %>% distinct(BCR)
      missing_bcrs <- all_bcrs %>% anti_join(safe_bcrs, by = "BCR")
      
      if (nrow(missing_bcrs) > 0) {
        sp_safe_dates <- bind_rows(
          sp_safe_dates,
          missing_bcrs %>%
            mutate(
              sp_english = sp_name,
              start_doy = fallback_start,
              end_doy   = fallback_end,
              midpoint  = floor((fallback_start + fallback_end) / 2)
            ) %>%
            select(sp_english, BCR, start_doy, end_doy, midpoint)
        )
      }
      
      pred_doy <- floor((fallback_start + fallback_end) / 2)
    } else {
      pred_doy <- NA_real_
    }
  }
  
  sp_dat <- sp_dat %>%
    left_join(sp_safe_dates, by = "BCR") %>%
    filter(
      !is.na(start_doy),
      !is.na(end_doy),
      DayOfYear >= start_doy,
      DayOfYear <= end_doy
    ) %>%
    mutate(days_midpoint = DayOfYear - pred_doy)
  
  
  # ------------------------------------------------------------
  # Error family triage
  # ------------------------------------------------------------
  
  large_counts <- sp_dat %>% filter(count > 50)
  sp_dat <- sp_dat %>% filter(count <= 50)
  
  n_det <- sum(sp_dat$count > 0)
  y_pos <- sp_dat$count[sp_dat$count > 0]
  
  if (length(y_pos) == 0) {
    message("Skipping ", sp_name, " because there are no positive counts after filtering.")
    next
  }
  
  n_top <- max(1, ceiling(0.01 * n_det))
  prop_total_top1pct_nonzero <- sum(sort(y_pos, decreasing = TRUE)[seq_len(n_top)]) / sum(y_pos)
  
  error_family <- "poisson"
  if (prop_total_top1pct_nonzero >= 0.10 || sum(y_pos > 15) > 10) {
    error_family <- "nbinomial"
  }
  
  dat_for_review <- sp_dat %>%
    select(Date_Time, Survey_Type, count, Hours_Since_Sunrise, Atlas)
  
  saveRDS(
    list(
      sp_english = sp_name,
      sp_code = sp_code,
      sp_safe_dates = sp_safe_dates,
      sp_dat = dat_for_review
    ),
    dat_path
  )
  
  # ------------------------------------------------------------
  # Fit model
  # ------------------------------------------------------------
  
  if (file.exists(pred_path) && !rerun_predictions) {
    message("Predictions already exist for ", sp_name, "; skipping")
    next
  }
  
  covars_present <- intersect(base_covars, names(sp_dat))
  
  # sd covariate values
  cov_sd <- sp_dat %>%
    as.data.frame() %>%
    select(covars_present) %>%
    apply(2,function(x)length(unique(x)))
  covars_present <- names(cov_sd)[which(cov_sd>0)]
  cov_df_sp <- make_cov_df(covars_present, mean = 0, sd_linear = 0.5)
  
  priors_list <- list(
    prior_range_abund = c(250, 0.5),
    prior_sigma_abund = c(1, 0.1),
    prior_range_change = c(250, 0.5),
    prior_sigma_change = c(0.1, 0.1),
    prior_HSS_range = c(3, 0.9),
    prior_HSS_sigma = c(3, 0.1),
    prior_DOY_range_global = c(7, 0.9),
    prior_DOY_sigma_global = c(3, 0.1),
    kappa_pcprec_diff = c(log(2), 0.1)
  )
  
  # ------------------------------------------------------------
  # Prepare mesh
  # ------------------------------------------------------------
  mesh_max_edge = c(70, 150)
  mesh_cutoff   = 20
  mesh_convex   = c(100, 200)
  mesh_concave  = c(50, 200)
  
  hull <- fmesher::fm_extensions(
    study_boundary,
    convex  = mesh_convex,
    concave = mesh_concave
  )
  
  mesh_spatial <- fmesher::fm_mesh_2d_inla(
    loc      = sf::st_as_sfc(sp_dat),
    boundary = hull,
    max.edge = mesh_max_edge,
    cutoff   = mesh_cutoff,
    crs      = sf::st_crs(sp_dat)
  )
  
  plot(mesh_spatial)
  
  plot(
    as(study_boundary, "Spatial"),
    add = TRUE,
    border = "red",
    lwd = 2
  )
  
  dim(mesh_spatial$loc)
  
  int_strategy <- 'eb' # "auto" # 
  strategy <- "simplified.laplace" # "laplace" # 
  
  start_model <- Sys.time()
  
  mod <- NULL
  if (file.exists(model_path) && !rerun_models) {
    mod <- readRDS(model_path)
    message("Loaded existing model: ", model_path)
  } else {
    mod <- try(
      fit_inla_multi_atlas(
        sp_dat = sp_dat,
        study_boundary = study_boundary,
        covariates = cov_df_sp,
        mesh_spatial = mesh_spatial,
        
        prior_range_abund = priors_list$prior_range_abund,
        prior_sigma_abund = priors_list$prior_sigma_abund,
        prior_range_change = priors_list$prior_range_change,
        prior_sigma_change = priors_list$prior_sigma_change,
        prior_HSS_range = priors_list$prior_HSS_range,
        prior_HSS_sigma = priors_list$prior_HSS_sigma,
        prior_DOY_range_global = priors_list$prior_DOY_range_global,
        prior_DOY_sigma_global = priors_list$prior_DOY_sigma_global,
        kappa_pcprec_diff = priors_list$kappa_pcprec_diff,
        int_strategy = int_strategy,
        strategy = strategy,
        family = error_family
      ),
      silent = TRUE
    )
    
    if (inherits(mod, "try-error") || is.null(mod)) {
      message("Model failed for ", sp_name, "; continuing.")
      next
    }
    
    # Do not save model, as file sizes are prohibitive
    # save_atomic(mod, model_path)
  }
  
  end_model <- Sys.time()
  fit_minutes <- round(as.numeric(end_model - start_model, units = "mins"))

  
  model_summaries[[sp_name]] <- list(
    sp_name = sp_name,
    sp_code = sp_code,
    error_family = error_family,
    priors = priors_list,
    summary_fixed = mod$summary.fixed,
    summary_hyperpar = mod$summary.hyperpar
  )
  save_atomic(model_summaries, model_summaries_path)
  
  message(
    "\n====================\n",
    i, "/", nrow(species_run), ": ", sp_name,
    " (species_id = ", sp_code, "); ",fit_minutes," min to fit model\n",
    "===================="
  )
  
  print(summary(mod))
  
  # ------------------------------------------------------------
  # Chunked predictions
  # ------------------------------------------------------------
  
  if (file.exists(pred_path) && !rerun_predictions) {
    message("Loaded existing predictions: ", pred_path)
  } else {
    start_prediction <- Sys.time()
    message("Generating predictions in chunks using column: ", chunk_id_col)
    
    pred_formula <- make_pred_formula_multiatlas(cov_df_sp)
    
    pixel_summary_chunks <- vector("list", length(chunk_ids))
    names(pixel_summary_chunks) <- as.character(chunk_ids)
    
    # one row per hex, each containing posterior draw vectors
    hex_draw_chunks <- vector("list", length(chunk_ids))
    names(hex_draw_chunks) <- as.character(chunk_ids)
    
    for (chunk_k in seq_along(chunk_ids)) {
      chunk_id <- chunk_ids[[chunk_k]]
      message("  chunk ", chunk_k, "/", length(chunk_ids), " (", chunk_id, ")")
      
      g2_chunk <- grid_OBBA2 %>%
        filter(.data[[chunk_id_col]] == chunk_id) %>%
        arrange(pixel_id)
      
      g3_chunk <- grid_OBBA3 %>%
        filter(.data[[chunk_id_col]] == chunk_id) %>%
        arrange(pixel_id)
      
      if (nrow(g2_chunk) == 0 && nrow(g3_chunk) == 0) next
      
      if (nrow(g2_chunk) != nrow(g3_chunk)) {
        stop("Chunk ", chunk_id, " has mismatched numbers of OBBA2 and OBBA3 pixels.")
      }
      
      if (!identical(g2_chunk$pixel_id, g3_chunk$pixel_id)) {
        stop("Chunk ", chunk_id, " has non-matching pixel_id order between OBBA2 and OBBA3.")
      }
      
      if (!identical(g2_chunk$hex_id, g3_chunk$hex_id)) {
        stop("Chunk ", chunk_id, " has non-matching hex_id order between OBBA2 and OBBA3.")
      }
      
      pred_grid_chunk <- make_pred_grid(g2_chunk, g3_chunk)
      
      # Use the SAME seed for every chunk so draw j stays aligned across chunks
      pred_chunk <- predict_one_chunk(
        mod = mod,
        pred_grid_chunk = pred_grid_chunk,
        pred_formula = pred_formula,
        n.samples = n_prediction_draws,
        seed = prediction_seed
      )
      
      # --------------------------------------------------------
      # 1) Per-pixel posterior summaries
      # --------------------------------------------------------
      
      pixel_summary_chunks[[chunk_k]] <- list(
        OBBA2 = bind_cols(
          g2_chunk %>% st_drop_geometry() %>% select(pixel_id, hex_id, all_of(chunk_id_col)),
          pred_chunk$summary$OBBA2
        ),
        OBBA3 = bind_cols(
          g3_chunk %>% st_drop_geometry() %>% select(pixel_id, hex_id, all_of(chunk_id_col)),
          pred_chunk$summary$OBBA3
        ),
        abs_change = bind_cols(
          g2_chunk %>% st_drop_geometry() %>% select(pixel_id, hex_id, all_of(chunk_id_col)),
          pred_chunk$summary$abs_change
        )
      )
      
      # --------------------------------------------------------
      # 2) Aggregate posterior draws to the hex level
      #
      # Store one draw vector per hex, which is sufficient for
      # later aggregation to any region larger than a hexagon.
      # --------------------------------------------------------
      
      hex_ids_chunk <- g2_chunk$hex_id
      u_hex <- unique(hex_ids_chunk)
      
      hex_draw_chunks[[chunk_k]] <- bind_rows(
        lapply(u_hex, function(hx) {
          
          idx_hex <- which(hex_ids_chunk == hx)
          
          mu2_hex <- colMeans(pred_chunk$mu2[idx_hex, , drop = FALSE])
          mu3_hex <- colMeans(pred_chunk$mu3[idx_hex, , drop = FALSE])
          
          tibble(
            hex_id = hx,
            !!chunk_id_col := chunk_id,
            n_pixels = length(idx_hex),
            mu_OBBA2 = list(mu2_hex),
            mu_OBBA3 = list(mu3_hex),
            abs_change = list(mu3_hex - mu2_hex)
          )
        })
      )
      
      rm(pred_chunk, pred_grid_chunk, g2_chunk, g3_chunk)
      gc(verbose = FALSE)
    }
    
    preds_OBBA2_summary <- bind_chunk_summaries(pixel_summary_chunks, "OBBA2")
    preds_OBBA3_summary <- bind_chunk_summaries(pixel_summary_chunks, "OBBA3")
    preds_abs_change_summary <- bind_chunk_summaries(pixel_summary_chunks, "abs_change")
    hex_draws <- bind_rows(hex_draw_chunks)
    
    end_prediction <- Sys.time()
    pred_minutes <- round(as.numeric(end_prediction - start_prediction, units = "mins"))
    
    message()
    sp_square_summary <- sp_dat %>%
      as.data.frame() %>%
      group_by(Atlas, square_id) %>%
      summarise(
        n_surveys    = n(),
        total_count  = sum(count),
        n_detections = sum(count > 0),
        BCR = names(which.max(table(BCR))),
        .groups = "drop"
      )
    
    save_atomic(
      list(
        sp_name = sp_name,
        sp_code = sp_code,
        sp_safe_dates = sp_safe_dates,
        sp_square_summary = sp_square_summary,
        error_family = error_family,
        priors = priors_list,
        fit_minutes = fit_minutes,
        summary_fixed = mod$summary.fixed,
        summary_hyperpar = mod$summary.hyperpar,
        pred_minutes = pred_minutes,
        pred_doy = pred_doy,
        prediction_seed = prediction_seed,
        n_prediction_draws = n_prediction_draws,
        chunk_id_col = chunk_id_col,
        OBBA2 = preds_OBBA2_summary,
        OBBA3 = preds_OBBA3_summary,
        abs_change = preds_abs_change_summary,
        hex_draws = hex_draws
      ),
      pred_path
    )
    
    message(
      "\n====================\n",
      i, "/", nrow(species_run), ": ", sp_name,
      " (species_id = ", sp_code, "); ",pred_minutes," min to generate predictions\n",
      "===================="
    )
    
  }
}

message("\n09_fit_models_and_predict_PC_ARU_chunked.R complete.")