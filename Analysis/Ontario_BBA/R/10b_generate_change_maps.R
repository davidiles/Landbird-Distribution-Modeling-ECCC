# ============================================================
# 10_make_species_figures.R
#
# Purpose:
#   Generate PNG figures per species using the per-species
#   prediction summary objects created in script 08.
#
# Inputs:
#   data_clean/birds/data_ready_for_analysis.rds
#   data_clean/model_output/predictions/<sp>.rds  (many)
#   + external shapefiles for BCR, water, (optional) atlas squares
#
# Outputs:
#   data_clean/model_output/figures/<sp>_*.png
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
  library(scales)
  library(patchwork)
  library(RColorBrewer)
  library(colorspace)
  library(ggtext)
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

model_type <- "PC_ARU"

# File paths
in_data  <- file.path(paths$data_clean, "birds", "data_ready_for_analysis.rds")
pred_dir <- file.path(paths$model_output, paste0("predictions_", model_type))
fig_dir  <- file.path(paths$model_output, paste0("figures_", model_type))
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

in_water <- file.path(paths$data_clean, "spatial", "water_filtered.shp")
in_atlas_squares <- file.path(paths$data, "Spatial", "National", "AtlasSquares", "NationalSquares_FINAL.shp")

# Plot export settings
dpi <- 1000
width_in <- 10
height_in <- 8
ggsave_type <- "cairo"

# ------------------------------------------------------------
# Custom function for generating a change map
# ------------------------------------------------------------

summarize_hex_change <- function(hex_grid,
                                 eta_draws_per_hex,
                                 min_pix_number = 100,
                                 baseline_min = 1e-5,
                                 ci_level = 0.95) {
  
  # ---------------------------------------------------------------------------
  # PURPOSE
  # ---------------------------------------------------------------------------
  # This function summarizes posterior draws of abundance and change for each
  # hexagon, and joins those summaries back onto an sf hex grid.
  #
  # INPUTS
  #   hex_grid:
  #     An sf object containing one row per hexagon, with a column called hex_id.
  #
  #   eta_draws_per_hex:
  #     A list containing posterior draws for each hexagon. At minimum, it must
  #     contain:
  #       - Mu2: matrix of posterior draws for Atlas 2 abundance
  #       - Mu3: matrix of posterior draws for Atlas 3 abundance
  #       - meta$n_pixels: number of 1-km pixels contributing to each hexagon
  #
  #     Mu2 and Mu3 are expected to have:
  #       - one row per hexagon
  #       - one column per posterior draw
  #       - rownames equal to hex_id values
  #
  #   min_pix_number:
  #     Minimum number of pixels required for a hexagon to be retained.
  #     Hexagons with fewer pixels than this are assigned NA for all summaries.
  #
  #   baseline_min:
  #     Small positive number used to avoid unstable proportional or log change
  #     calculations when Mu2 is extremely close to zero.
  #
  #   ci_level:
  #     Width of the credible interval to calculate. For example:
  #       0.95 = 95% credible interval
  #       0.90 = 90% credible interval
  #
  # OUTPUT
  #   Returns the original sf object, with posterior summaries joined by hex_id.
  #
  # ---------------------------------------------------------------------------
  # BASIC INPUT CHECKS
  # ---------------------------------------------------------------------------
  
  # Check that ci_level is sensible
  if (!is.numeric(ci_level) || length(ci_level) != 1 ||
      ci_level <= 0 || ci_level >= 1) {
    stop("ci_level must be a single number strictly between 0 and 1.")
  }
  
  # Convert the requested interval width into lower and upper tail probabilities.
  # Example:
  #   ci_level = 0.95  -> alpha = 0.025
  #   lower quantile   = 0.025
  #   upper quantile   = 0.975
  alpha <- (1 - ci_level) / 2
  
  # ---------------------------------------------------------------------------
  # EXTRACT REQUIRED OBJECTS
  # ---------------------------------------------------------------------------
  
  Mu2  <- eta_draws_per_hex$Mu2
  Mu3  <- eta_draws_per_hex$Mu3
  npix <- eta_draws_per_hex$meta$n_pixels
  
  # Basic checks on posterior draw matrices
  if (is.null(Mu2) || is.null(Mu3)) {
    stop("eta_draws_per_hex must contain elements named 'Mu2' and 'Mu3'.")
  }
  
  if (!is.matrix(Mu2) || !is.matrix(Mu3)) {
    stop("Mu2 and Mu3 must both be matrices.")
  }
  
  if (!identical(dim(Mu2), dim(Mu3))) {
    stop("Mu2 and Mu3 must have identical dimensions.")
  }
  
  if (is.null(rownames(Mu2)) || is.null(rownames(Mu3))) {
    stop("Mu2 and Mu3 must both have rownames corresponding to hex_id.")
  }
  
  if (!identical(rownames(Mu2), rownames(Mu3))) {
    stop("Rownames of Mu2 and Mu3 must match exactly.")
  }
  
  # Extract hex IDs from the matrix rownames
  hex_ids <- as.integer(rownames(Mu2))
  
  # ---------------------------------------------------------------------------
  # ALIGN THE NUMBER OF PIXELS TO THE POSTERIOR MATRICES
  # ---------------------------------------------------------------------------
  # We want one pixel count per row of Mu2 / Mu3.
  #
  # If n_pixels is named by hex_id, align by name.
  # If it is unnamed, assume it is already in the same order as Mu2 / Mu3.
  
  if (is.null(npix)) {
    stop("eta_draws_per_hex$meta$n_pixels is missing.")
  }
  
  if (is.null(names(npix))) {
    if (length(npix) != nrow(Mu2)) {
      stop("n_pixels has no names and its length does not match nrow(Mu2).")
    }
    npix_aligned <- npix
  } else {
    npix_aligned <- npix[as.character(hex_ids)]
  }
  
  # Identify hexagons with enough underlying pixels to keep
  keep_hex <- !is.na(npix_aligned) & npix_aligned >= min_pix_number
  
  # ---------------------------------------------------------------------------
  # SUMMARIZE POSTERIOR ABUNDANCE FOR EACH ATLAS
  # ---------------------------------------------------------------------------
  # We calculate the posterior median abundance within each hexagon for Atlas 2
  # and Atlas 3.
  
  mu2_med <- apply(Mu2, 1, median, na.rm = TRUE)
  mu3_med <- apply(Mu3, 1, median, na.rm = TRUE)
  
  # Absolute change is simply the difference in posterior medians
  abs_change_med <- mu3_med - mu2_med
  
  # ---------------------------------------------------------------------------
  # CALCULATE RAW PROPORTIONAL CHANGE DRAW-BY-DRAW
  # ---------------------------------------------------------------------------
  # Proportional change is defined as:
  #
  #   (Mu3 - Mu2) / Mu2
  #
  # This gives values such as:
  #   0.25  = +25% change
  #  -0.40  = -40% change
  #
  # We calculate this for every posterior draw, then summarize the resulting
  # distribution for each hexagon.
  #
  # If Mu2 is extremely small, the ratio becomes unstable, so we set those
  # entries to NA.
  
  prop_change_draws <- (Mu3 - Mu2) / Mu2
  prop_change_draws[Mu2 <= baseline_min] <- NA_real_
  
  # Posterior summaries for proportional change
  prop_change_mean   <- rowMeans(prop_change_draws, na.rm = TRUE)
  prop_change_median <- apply(prop_change_draws, 1, median,   na.rm = TRUE)
  prop_change_qlow   <- apply(prop_change_draws, 1, quantile, probs = alpha,       na.rm = TRUE)
  prop_change_qhigh  <- apply(prop_change_draws, 1, quantile, probs = 1 - alpha,   na.rm = TRUE)
  
  # ---------------------------------------------------------------------------
  # CALCULATE SYMMETRIC LOG CHANGE DRAW-BY-DRAW
  # ---------------------------------------------------------------------------
  # Symmetric log change is defined as:
  #
  #   log(Mu3 / Mu2)
  #
  # This is often preferable for change mapping because it is symmetric around 0:
  #   no change   -> 0
  #   doubling    -> +log(2)
  #   halving     -> -log(2)
  #
  # As above, values are set to NA when Mu2 is too close to zero.
  
  sym_change_draws <- log(Mu3 / Mu2)
  sym_change_draws[Mu2 <= baseline_min] <- NA_real_
  
  # Posterior summaries for symmetric log change
  sym_change_mean   <- rowMeans(sym_change_draws, na.rm = TRUE)
  sym_change_median <- apply(sym_change_draws, 1, median,   na.rm = TRUE)
  sym_change_qlow   <- apply(sym_change_draws, 1, quantile, probs = alpha,       na.rm = TRUE)
  sym_change_qhigh  <- apply(sym_change_draws, 1, quantile, probs = 1 - alpha,   na.rm = TRUE)
  
  # ---------------------------------------------------------------------------
  # CALCULATE POSTERIOR PROBABILITY OF INCREASE OR DECREASE
  # ---------------------------------------------------------------------------
  # These probabilities describe how much posterior support there is for the
  # direction of change.
  #
  # For example:
  #   p_increase = 0.99 means 99% of posterior draws imply an increase
  #   p_decrease = 0.97 means 97% of posterior draws imply a decrease
  
  p_increase <- rowMeans(prop_change_draws > 0, na.rm = TRUE)
  p_decrease <- rowMeans(prop_change_draws < 0, na.rm = TRUE)
  
  # A simple direction label based on strong posterior support.
  # Here, "strong" is defined as > 0.975, but you could change that if desired.
  direction <- dplyr::case_when(
    p_increase > 0.975 ~ "increase",
    p_decrease > 0.975 ~ "decrease",
    TRUE ~ "uncertain"
  )
  
  # ---------------------------------------------------------------------------
  # APPLY THE MINIMUM-PIXEL FILTER
  # ---------------------------------------------------------------------------
  # Hexagons with too few underlying pixels are treated as unreliable and are
  # assigned NA for all summary quantities.
  
  mu2_med[!keep_hex]            <- NA_real_
  mu3_med[!keep_hex]            <- NA_real_
  abs_change_med[!keep_hex]     <- NA_real_
  
  prop_change_mean[!keep_hex]   <- NA_real_
  prop_change_median[!keep_hex] <- NA_real_
  prop_change_qlow[!keep_hex]   <- NA_real_
  prop_change_qhigh[!keep_hex]  <- NA_real_
  
  sym_change_mean[!keep_hex]    <- NA_real_
  sym_change_median[!keep_hex]  <- NA_real_
  sym_change_qlow[!keep_hex]    <- NA_real_
  sym_change_qhigh[!keep_hex]   <- NA_real_
  
  p_increase[!keep_hex]         <- NA_real_
  p_decrease[!keep_hex]         <- NA_real_
  direction[!keep_hex]          <- NA_character_
  
  # ---------------------------------------------------------------------------
  # BUILD A SUMMARY TABLE
  # ---------------------------------------------------------------------------
  # This is a plain data frame/tibble with one row per hexagon and the posterior
  # summaries we just calculated.
  
  hex_summary <- tibble::tibble(
    hex_id = hex_ids,
    n_pixels = npix_aligned,
    mu2_med = mu2_med,
    mu3_med = mu3_med,
    abs_change_med = abs_change_med,
    
    prop_change_mean = prop_change_mean,
    prop_change_median = prop_change_median,
    prop_change_qlow = prop_change_qlow,
    prop_change_qhigh = prop_change_qhigh,
    
    sym_change_mean = sym_change_mean,
    sym_change_median = sym_change_median,
    sym_change_qlow = sym_change_qlow,
    sym_change_qhigh = sym_change_qhigh,
    
    p_increase = p_increase,
    p_decrease = p_decrease,
    direction = direction
  )
  
  # ---------------------------------------------------------------------------
  # JOIN THE SUMMARIES BACK TO THE HEX GRID
  # ---------------------------------------------------------------------------
  # The result is an sf object that contains all the original geometry, plus the
  # new summary columns for mapping or further analysis.
  
  hex_grid_out <- hex_grid %>%
    dplyr::left_join(hex_summary, by = "hex_id")
  
  return(hex_grid_out)
}


classify_min_supported_change <- function(hex_sf) {
  
  hex_sf %>%
    dplyr::mutate(
      min_supported_change = dplyr::case_when(
        
        is.na(sym_change_qlow) | is.na(sym_change_qhigh) ~ NA_character_,
        
        # ---- Increases (lower CI bound above threshold) ----
        sym_change_qlow >= log(2.0)  ~ "> 2x",
        sym_change_qlow >= log(1.5)  ~ "1.5x to 2x",
        sym_change_qlow >= log(1.1)  ~ "1.1x to 1.5x",
        
        # ---- Declines (upper CI bound below threshold) ----
        sym_change_qhigh <= log(1/2)   ~ "< 0.5x",        # mirror of > 2x
        sym_change_qhigh <= log(1/1.5) ~ "0.5x to 0.67x", # mirror of 1.5x to 2x
        sym_change_qhigh <= log(1/1.1) ~ "0.67x to 0.9x", # mirror of 1.1x to 1.5x
        
        # ---- Neutral: CI overlaps the no-change zone ----
        TRUE ~ "0.9x to 1.1x"
        
      ),
      
      min_supported_change = factor(
        min_supported_change,
        levels = c(
          "< 0.5x",
          "0.5x to 0.67x",
          "0.67x to 0.9x",
          "0.9x to 1.1x",
          "1.1x to 1.5x",
          "1.5x to 2x",
          "> 2x"
        ),
        ordered = TRUE
      )
    )
}


summarize_provincial_change <- function(eta_draws_per_hex,
                                        min_pix_number = 100,
                                        baseline_min = 1e-5,
                                        ci_level = 0.95,
                                        return_draws = FALSE) {
  
  # ----------------------------------------------------------
  # PURPOSE
  # ----------------------------------------------------------
  # Summarize province-wide abundance and change by combining
  # posterior draws across hexagons, weighted by the number of
  # 1-km pixels contributing to each hexagon.
  #
  # ASSUMPTION:
  #   Mu2 and Mu3 are posterior expected abundance values on a
  #   per-pixel (or average-within-hex) basis, so province-wide
  #   totals are obtained as weighted sums across hexagons.
  #
  # INPUTS
  #   eta_draws_per_hex:
  #     List with:
  #       - Mu2: matrix [n_hex x n_draw]
  #       - Mu3: matrix [n_hex x n_draw]
  #       - meta$n_pixels: vector of pixel counts per hex
  #
  #   min_pix_number:
  #     Hexagons with fewer than this many pixels are excluded.
  #
  #   baseline_min:
  #     Small positive number used to avoid unstable proportional
  #     and log-change calculations when province-wide Mu2 is near 0.
  #
  #   ci_level:
  #     Credible interval width, e.g. 0.95 for 95% CI.
  #
  #   return_draws:
  #     If TRUE, also return the full province-wide posterior draws.
  #
  # OUTPUT
  #   A tibble with one row summarizing province-wide change.
  #   Optionally returns the underlying draw vectors as well.
  # ----------------------------------------------------------
  
  # ---- checks
  if (!is.numeric(ci_level) || length(ci_level) != 1 ||
      ci_level <= 0 || ci_level >= 1) {
    stop("ci_level must be a single number strictly between 0 and 1.")
  }
  
  alpha <- (1 - ci_level) / 2
  
  Mu2  <- eta_draws_per_hex$Mu2
  Mu3  <- eta_draws_per_hex$Mu3
  npix <- eta_draws_per_hex$meta$n_pixels
  
  if (is.null(Mu2) || is.null(Mu3)) {
    stop("eta_draws_per_hex must contain elements named 'Mu2' and 'Mu3'.")
  }
  
  if (!is.matrix(Mu2) || !is.matrix(Mu3)) {
    stop("Mu2 and Mu3 must both be matrices.")
  }
  
  if (!identical(dim(Mu2), dim(Mu3))) {
    stop("Mu2 and Mu3 must have identical dimensions.")
  }
  
  if (is.null(rownames(Mu2)) || is.null(rownames(Mu3))) {
    stop("Mu2 and Mu3 must both have rownames corresponding to hex_id.")
  }
  
  if (!identical(rownames(Mu2), rownames(Mu3))) {
    stop("Rownames of Mu2 and Mu3 must match exactly.")
  }
  
  if (is.null(npix)) {
    stop("eta_draws_per_hex$meta$n_pixels is missing.")
  }
  
  hex_ids <- rownames(Mu2)
  
  # ---- align pixel counts
  if (is.null(names(npix))) {
    if (length(npix) != nrow(Mu2)) {
      stop("n_pixels has no names and its length does not match nrow(Mu2).")
    }
    npix_aligned <- npix
  } else {
    npix_aligned <- npix[hex_ids]
  }
  
  # ---- keep only sufficiently large hexes
  keep_hex <- !is.na(npix_aligned) & npix_aligned >= min_pix_number
  
  if (!any(keep_hex)) {
    stop("No hexagons remain after applying min_pix_number.")
  }
  
  Mu2_keep <- Mu2[keep_hex, , drop = FALSE]
  Mu3_keep <- Mu3[keep_hex, , drop = FALSE]
  w_keep   <- npix_aligned[keep_hex]
  
  # ---- province-wide posterior totals for each draw
  # weighted sum across hexagons
  prov_mu2_draws <- colSums(Mu2_keep * w_keep, na.rm = TRUE)
  prov_mu3_draws <- colSums(Mu3_keep * w_keep, na.rm = TRUE)
  
  # ---- draw-wise change metrics
  abs_change_draws  <- prov_mu3_draws - prov_mu2_draws
  prop_change_draws <- (prov_mu3_draws - prov_mu2_draws) / prov_mu2_draws
  sym_change_draws  <- log(prov_mu3_draws / prov_mu2_draws)
  
  prop_change_draws[prov_mu2_draws <= baseline_min] <- NA_real_
  sym_change_draws[prov_mu2_draws <= baseline_min]  <- NA_real_
  
  # ---- summarize
  out <- tibble::tibble(
    n_hex_total = nrow(Mu2),
    n_hex_used = sum(keep_hex),
    n_pixels_total_used = sum(w_keep, na.rm = TRUE),
    
    mu2_mean = mean(prov_mu2_draws, na.rm = TRUE),
    mu2_median = median(prov_mu2_draws, na.rm = TRUE),
    mu2_qlow = unname(stats::quantile(prov_mu2_draws, probs = alpha, na.rm = TRUE)),
    mu2_qhigh = unname(stats::quantile(prov_mu2_draws, probs = 1 - alpha, na.rm = TRUE)),
    
    mu3_mean = mean(prov_mu3_draws, na.rm = TRUE),
    mu3_median = median(prov_mu3_draws, na.rm = TRUE),
    mu3_qlow = unname(stats::quantile(prov_mu3_draws, probs = alpha, na.rm = TRUE)),
    mu3_qhigh = unname(stats::quantile(prov_mu3_draws, probs = 1 - alpha, na.rm = TRUE)),
    
    abs_change_mean = mean(abs_change_draws, na.rm = TRUE),
    abs_change_median = median(abs_change_draws, na.rm = TRUE),
    abs_change_qlow = unname(stats::quantile(abs_change_draws, probs = alpha, na.rm = TRUE)),
    abs_change_qhigh = unname(stats::quantile(abs_change_draws, probs = 1 - alpha, na.rm = TRUE)),
    
    prop_change_mean = mean(prop_change_draws, na.rm = TRUE),
    prop_change_median = median(prop_change_draws, na.rm = TRUE),
    prop_change_qlow = unname(stats::quantile(prop_change_draws, probs = alpha, na.rm = TRUE)),
    prop_change_qhigh = unname(stats::quantile(prop_change_draws, probs = 1 - alpha, na.rm = TRUE)),
    
    sym_change_mean = mean(sym_change_draws, na.rm = TRUE),
    sym_change_median = median(sym_change_draws, na.rm = TRUE),
    sym_change_qlow = unname(stats::quantile(sym_change_draws, probs = alpha, na.rm = TRUE)),
    sym_change_qhigh = unname(stats::quantile(sym_change_draws, probs = 1 - alpha, na.rm = TRUE)),
    
    p_increase = mean(abs_change_draws > 0, na.rm = TRUE),
    p_decrease = mean(abs_change_draws < 0, na.rm = TRUE),
    
    direction = dplyr::case_when(
      mean(abs_change_draws > 0, na.rm = TRUE) > 0.975 ~ "increase",
      mean(abs_change_draws < 0, na.rm = TRUE) > 0.975 ~ "decrease",
      TRUE ~ "uncertain"
    )
  )
  
  if (return_draws) {
    return(list(
      summary = out,
      draws = list(
        mu2 = prov_mu2_draws,
        mu3 = prov_mu3_draws,
        abs_change = abs_change_draws,
        prop_change = prop_change_draws,
        sym_change = sym_change_draws
      )
    ))
  }
  
  return(out)
}

# ------------------------------------------------------------
# Load base data
# ------------------------------------------------------------

stopifnot(file.exists(in_data))
dat <- readRDS(in_data)

# raw data
all_surveys <- dat$all_surveys
counts      <- dat$counts

# prediction grid / study area
grid2 <- dat$grid_OBBA2
grid3 <- dat$grid_OBBA3
study_boundary <- dat$study_boundary %>% sf::st_as_sf()

hex_grid <- dat$hex_grid

stopifnot(inherits(grid2, "sf"), inherits(grid3, "sf"))

# Ensure consistent CRS
crs_use <- st_crs(grid2)
study_boundary <- st_transform(study_boundary, crs_use)
grid3 <- st_transform(grid3, crs_use)
grid2 <- st_transform(grid2, crs_use)

# BCR outlines
bcr_sf <- dat$bcr_sf

# # Water
# water_sf <- NULL
# if (file.exists(in_water)) {
#   water_sf <- st_read(in_water, quiet = TRUE) %>%
#     st_make_valid() %>%
#     st_transform(crs_use)
# }

# ------------------------------------------------------------
# Find prediction files
# ------------------------------------------------------------

pred_files <- list.files(pred_dir, pattern = "\\.rds$", full.names = TRUE)
stopifnot(length(pred_files) > 0)

message("Prediction files found: ", length(pred_files))

# ------------------------------------------------------------
# Loop species
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
  
  fig_path <- file.path(fig_dir, paste0(sp_file, "_chg.png"))
  
  # if (file.exists(fig_path)) {
  #   message("Skipping ", sp_english, ": change map already exists for this species")
  #   next
  # }
  
  # ----------------------------------------------------------
  # Safety checks
  # ----------------------------------------------------------
  
  message("Mapping: ", sp_english)
  
  # Estimate of overall provincial change
  prov_change <- summarize_provincial_change(
    eta_draws_per_hex = preds$eta_draws_per_hex,
    min_pix_number = 100,
    ci_level = 0.90
  )%>%
    dplyr::mutate(
      pct_change_median = 100 * prop_change_median,
      pct_change_qlow   = 100 * prop_change_qlow,
      pct_change_qhigh  = 100 * prop_change_qhigh
    )
  
  
  hex_change_sf <- summarize_hex_change(
    hex_grid = hex_grid,
    eta_draws_per_hex = preds$eta_draws_per_hex,
    ci_level = 0.90
  )
  
  hex_change_sf <- hex_change_sf %>%
    classify_min_supported_change()
  
  zmax <- max(abs(hex_change_sf$mu3_med),na.rm = TRUE)
  zmin <- min(c(0.01,(zmax/20)))
  
  hex_change_sf <- hex_change_sf %>%
    filter(mu2_med>zmin)
  
  
  # Define levels and colours in one place
  change_levels <- c(
    "< 0.5x", "0.5x to 0.67x", "0.67x to 0.9x",
    "0.9x to 1.1x",
    "1.1x to 1.5x", "1.5x to 2x", "> 2x"
  )
  change_labels <- c(
    "Large decrease\n(at least 50%)", "Moderate decrease\n(at least 33%)", "Small decrease\n(at least 10%)",
    "Little change\n(-10% to +10%)",
    "Small increase\n(at least 10%)", "Moderate increase\n(at least 50%)", "Large increase\n(at least 100%)"
  )
  
  tweak_colour <- function(col, lighten_amount = -0.15, saturate_amount = 0.2) {
    hsl <- as(hex2RGB(col), "HLS")
    hsl@coords[, "L"] <- pmin(1, pmax(0, hsl@coords[, "L"] + lighten_amount))
    hsl@coords[, "S"] <- pmin(1, pmax(0, hsl@coords[, "S"] + saturate_amount))
    hex(hsl)
  }
  
  change_colours <- RColorBrewer::brewer.pal(7, "RdBu")
  change_colours[4] <- "white"
  change_colours[-4] <- sapply(change_colours[-4], tweak_colour,
                               lighten_amount = -0.1,   # negative = darker
                               saturate_amount = 0.1)
  
  
  # Build the legend as its own ggplot
  legend_df <- tibble(
    level   = factor(change_levels, levels = rev(change_levels)),  # rev so top = strongest increase
    label   = factor(change_labels, levels = rev(change_labels)),  # keep in same order
    colour  = change_colours
  )
  
  legend_plot <- ggplot(legend_df, aes(x = 1, y = level)) +
    geom_tile(
      aes(fill = level),
      colour    = "grey40",
      width     = 0.6,
      height    = 0.8,
      linewidth = 0.3
    ) +
    geom_text(
      aes(
        x     = 1.55,
        label = label
      ),
      hjust = 0,
      size  = 3
    ) +
    scale_fill_manual(values  = setNames(change_colours, change_levels), guide = "none") +
    labs(title = "Minimum supported\npopulation change") +
    scale_x_continuous(limits = c(0.6, 5.5), expand = c(0, 0)) +  # wider for longer labels
    theme_void() +
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0.15, vjust = 1, margin = margin(b = 4))
    )
  
  # Main map — no legend
  
  # create title describing percent change
  fmt_pct <- function(x, digits = 1) {
    ifelse(
      is.na(x),
      NA_character_,
      sprintf("%+.1f", round(100 * x, digits))
    )
  }
  
  title_text <- paste0(
    "<span style='font-size:18pt;'><b>", sp_english, "</b></span><br><br>",
    "<span style='font-size:14pt;'>",
    "Overall change = ",
    fmt_pct(prov_change$prop_change_median), "% [",
    fmt_pct(prov_change$prop_change_qlow), "% to ",
    fmt_pct(prov_change$prop_change_qhigh), "%]",
    "</span>"
  )
  
  chg_plot <- ggplot() +
    geom_sf(data = study_boundary, fill = "gray90") +
    geom_sf(
      data = st_centroid(hex_change_sf),
      aes(colour = min_supported_change, size = mu3_med)
    ) +
    scale_colour_manual(
      values = setNames(change_colours, change_levels),
      guide  = "none"
    ) +
    scale_size_continuous(range = c(0, 2), guide = "none") +
    theme_bw()+
    ggtitle(title_text) +
    theme(
      plot.title = ggtext::element_markdown(lineheight = 1.1)
    )

  
  final_plot <- chg_plot + 
    inset_element(
      legend_plot,
      left   = 0.8,
      bottom = 0.5,
      right  = 1.0,
      top    = 0.99
    )
  final_plot
  
  ggsave(
    filename = fig_path,
    plot = final_plot, width = width_in, height = height_in, units = "in",
    dpi = dpi, type = ggsave_type, limitsize = FALSE
  )
  
}

message("10_make_species_figures.R complete: ", fig_dir)