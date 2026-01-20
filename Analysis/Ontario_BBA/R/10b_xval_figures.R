suppressPackageStartupMessages({
  library(dplyr)
  library(sf)
  library(ggplot2)
  library(stringr)
  library(readr)
  library(scales)
})

out_dir <- "data_clean/model_output/cv"

# Helper: list species block summary files that exist so far
list_completed_species <- function() {
  files <- list.files(file.path(out_dir, "block_summaries"), pattern = "\\.rds$", full.names = TRUE)
  tibble(
    path = files,
    sp_file = basename(files) %>% str_replace("\\.rds$", "")
  ) %>%
    arrange(sp_file)
}

# Helper: safe read (file may be mid-write briefly)
safe_read_rds <- function(path, tries = 3, wait_sec = 0.5) {
  for (k in seq_len(tries)) {
    x <- try(readRDS(path), silent = TRUE)
    if (!inherits(x, "try-error")) return(x)
    Sys.sleep(wait_sec)
  }
  stop("Could not read: ", path)
}

# Load one species' block summaries (all reps/folds completed so far)
load_block_summaries <- function(sp_file) {
  path <- file.path(out_dir, "block_summaries", paste0(sp_file, ".rds"))
  safe_read_rds(path)
}

# Filter to a specific rep/fold for “single experiment” views
filter_rep_fold <- function(block_summ, rep = 1, fold = 1) {
  block_summ %>% filter(rep == rep, fold == fold)
}

# Plot 1: calibration scatter for a rep/fold
plot_block_calibration <- function(block_summ_rf, log1p_axes = FALSE) {
  p <- ggplot(block_summ_rf, aes(x = mean_obs, y = mean_pred_q50)) +
    geom_abline(intercept = 0, slope = 1, linetype = 2) +
    geom_errorbar(aes(ymin = mean_pred_q025, ymax = mean_pred_q975), width = 0, alpha = 0.5) +
    geom_point(aes(size = n, color = bias), alpha = 0.9) +
    scale_color_gradient2(midpoint = 0, low = "blue", mid = "grey90", high = "red") +
    coord_equal() +
    theme_bw()
  
  if (log1p_axes) {
    p <- p + scale_x_continuous(trans = "log1p") + scale_y_continuous(trans = "log1p")
  }
  p
}

# Plot 2: bias map (symmetric scaling around 0)
plot_bias_map <- function(block_summ_rf) {
  bmax <- max(abs(block_summ_rf$bias), na.rm = TRUE)
  
  ggplot(block_summ_rf) +
    geom_sf(aes(fill = bias), color = NA) +
    scale_fill_gradient2(
      midpoint = 0,
      limits = c(-bmax, bmax),
      low = "blue", mid = "grey95", high = "red",
      oob = squish
    ) +
    theme_bw()
}

# Plot 3: “outside 95%” map
plot_outside_95_map <- function(block_summ_rf) {
  x <- block_summ_rf %>%
    mutate(outside_95 = (mean_obs < mean_pred_q025) | (mean_obs > mean_pred_q975))
  
  ggplot(x) +
    geom_sf(aes(fill = outside_95), color = NA) +
    scale_fill_manual(values = c(`FALSE` = "grey90", `TRUE` = "black")) +
    theme_bw()
}

# ---- Example interactive usage ----
# 1) See what species are available so far:
avail <- list_completed_species()
print(avail)

# 2) Choose one species file name from avail$sp_file
sp_file <- avail$sp_file[1]

block_summ <- load_block_summaries(sp_file)

# 3) Pick a rep/fold you want to inspect
rep_pick <- 1
fold_pick <- 1
block_summ_rf <- filter_rep_fold(block_summ, rep = rep_pick, fold = fold_pick)

# 4) Make plots
print(plot_block_calibration(block_summ_rf, log1p_axes = TRUE))
print(plot_bias_map(block_summ_rf))
print(plot_outside_95_map(block_summ_rf))
