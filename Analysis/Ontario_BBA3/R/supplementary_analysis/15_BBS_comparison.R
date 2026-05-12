# ============================================================
# Compare BCR-level atlas change estimates to BBS
# Scale: overall percent change across 2002-2022
# ============================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(here)
  library(tidyr)
  library(ggrepel)
  library(purrr)
})

source(here::here("R", "00_config_paths.R"))

# helper-functions
source(file.path(paths$functions, "inla_model_utils.R"))
source(file.path(paths$functions, "figure_utils.R"))

# ------------------------------------------------------------
# Config
# ------------------------------------------------------------

model_type <- "joint2"
atlas_year_gap <- 20
plot_regions <- c("CA-ON-13", "CA-ON-12", "CA-ON-8", "CA-ON-7")
eps <- 1e-9

in_data  <- file.path(paths$data_clean, "birds", "data_ready_for_analysis.rds")
pred_dir <- file.path(paths$model_output, paste0("predictions_", model_type))
fig_dir  <- file.path(paths$model_output, paste0("figures_", model_type))
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

pct_to_log <- function(pct) log1p(pct / 100)
log_to_pct <- function(x) 100 * (exp(x) - 1)

weighted_cor <- function(x, y, w) {
  mx <- weighted.mean(x, w, na.rm = TRUE)
  my <- weighted.mean(y, w, na.rm = TRUE)
  cov_xy <- sum(w * (x - mx) * (y - my), na.rm = TRUE) / sum(w, na.rm = TRUE)
  var_x  <- sum(w * (x - mx)^2,       na.rm = TRUE) / sum(w, na.rm = TRUE)
  var_y  <- sum(w * (y - my)^2,       na.rm = TRUE) / sum(w, na.rm = TRUE)
  cov_xy / sqrt(var_x * var_y)
}

summarize_draws <- function(mat, prefix) {
  tibble(
    !!paste0(prefix, "_q50")  := apply(mat, 1, median),
    !!paste0(prefix, "_q025") := apply(mat, 1, quantile, probs = 0.025),
    !!paste0(prefix, "_q975") := apply(mat, 1, quantile, probs = 0.975)
  )
}

align_hex_to_draws <- function(hex_sf, poly_ids) {
  out <- hex_sf %>%
    mutate(hex_id_chr = as.character(hex_id)) %>%
    slice(match(as.character(poly_ids), hex_id_chr)) %>%
    select(-hex_id_chr) %>%
    mutate(mu_row = row_number())
  
  stopifnot(nrow(out) == length(poly_ids))
  stopifnot(all(out$hex_id == poly_ids))
  out
}

build_bcr_weight_matrix <- function(sp_hex, bcr_sf) {
  bcr_sf_use <- st_transform(bcr_sf, st_crs(sp_hex)) %>%
    mutate(bcr_id = row_number())
  
  hex_bcr <- st_intersection(sp_hex, bcr_sf_use) %>%
    mutate(area_piece = st_area(geometry))
  
  hex_area <- sp_hex %>%
    transmute(mu_row, hex_area = st_area(geometry)) %>%
    st_drop_geometry()
  
  hex_bcr <- hex_bcr %>%
    left_join(hex_area, by = "mu_row") %>%
    mutate(weight = as.numeric(area_piece / hex_area))
  
  W <- matrix(0, nrow = nrow(bcr_sf_use), ncol = nrow(sp_hex))
  for (i in seq_len(nrow(hex_bcr))) {
    W[hex_bcr$bcr_id[i], hex_bcr$mu_row[i]] <-
      W[hex_bcr$bcr_id[i], hex_bcr$mu_row[i]] + hex_bcr$weight[i]
  }
  
  list(W = W, bcr_sf_use = bcr_sf_use)
}

summarize_bcr_detections <- function(sp_dat) {
  sp_dat %>%
    st_drop_geometry() %>%
    filter(count > 0) %>%
    group_by(Atlas, BCR) %>%
    summarise(n_sq = n_distinct(square_id), .groups = "drop") %>%
    pivot_wider(
      names_from = Atlas,
      values_from = n_sq,
      names_prefix = "n_sq_"
    ) %>%
    mutate(
      n_sq_OBBA2 = coalesce(n_sq_OBBA2, 0L),
      n_sq_OBBA3 = coalesce(n_sq_OBBA3, 0L)
    )
}

process_species_file <- function(pf, dat, all_surveys, counts, bcr_sf,hex_sf,
                                 atlas_year_gap = 20, eps = 1e-9) {
  preds <- readRDS(pf)
  sp_english <- preds$sp_name
  message(sp_english)
  
  sp_code <- dat$species_to_model %>%
    filter(english_name == sp_english) %>%
    pull(species_id)
  
  sp_dat <- all_surveys %>%
    mutate(count = counts[[sp_code]])
  
  eta <- preds$eta_draws_per_hex
  Mu2 <- eta$Mu2
  Mu3 <- eta$Mu3
  poly_ids <- eta$meta$poly_ids
  
  sp_hex <- align_hex_to_draws(hex_sf, poly_ids)
  stopifnot(nrow(sp_hex) == nrow(Mu2), nrow(sp_hex) == nrow(Mu3))
  
  w_obj <- build_bcr_weight_matrix(sp_hex, bcr_sf)
  W <- w_obj$W
  bcr_sf_use <- w_obj$bcr_sf_use
  
  Mu2_sum <- W %*% Mu2
  Mu3_sum <- W %*% Mu3
  
  prop_change <- (Mu3_sum - Mu2_sum) / pmax(Mu2_sum, eps)
  annual_change <- (Mu3_sum / pmax(Mu2_sum, eps))^(1 / atlas_year_gap) - 1
  
  change_summary <- tibble(
    english_name = sp_english,
    bcr_id = bcr_sf_use$bcr_id
  ) %>%
    bind_cols(summarize_draws(prop_change, "propchg")) %>%
    bind_cols(summarize_draws(annual_change, "amchg"))
  
  det_summary <- summarize_bcr_detections(sp_dat)
  
  bcr_sf_use %>%
    st_drop_geometry() %>%
    left_join(change_summary, by = "bcr_id") %>%
    left_join(det_summary, by = "BCR") %>%
    mutate(
      n_sq_OBBA2 = coalesce(n_sq_OBBA2, 0L),
      n_sq_OBBA3 = coalesce(n_sq_OBBA3, 0L),
      
      # Atlas overall percent change across the full interval
      atlas_median = 100 * propchg_q50,
      atlas_q025   = 100 * propchg_q025,
      atlas_q975   = 100 * propchg_q975
    ) %>%
    select(
      english_name, BCR, n_sq_OBBA2, n_sq_OBBA3,
      atlas_median, atlas_q025, atlas_q975,
      propchg_q50, propchg_q025, propchg_q975,
      amchg_q50, amchg_q025, amchg_q975
    )
}

# ------------------------------------------------------------
# Load base data
# ------------------------------------------------------------

dat <- readRDS(in_data)

all_surveys <- dat$all_surveys
counts      <- dat$counts
bcr_sf      <- st_transform(dat$bcr_sf, st_crs(dat$grid_OBBA2))
hex_sf      <- dat$hex_grid

pred_files <- list.files(pred_dir, pattern = "\\.rds$", full.names = TRUE)
stopifnot(length(pred_files) > 0)
message("Prediction files found: ", length(pred_files))

# ------------------------------------------------------------
# Build atlas summaries
# ------------------------------------------------------------

atlas_change_summary <- map_dfr(
  pred_files,
  process_species_file,
  dat = dat,
  all_surveys = all_surveys,
  counts = counts,
  bcr_sf = bcr_sf,
  hex_sf = hex_sf,
  atlas_year_gap = atlas_year_gap,
  eps = eps
) %>%
  mutate(region = paste0("CA-ON-", BCR))

saveRDS(atlas_change_summary,file = paste0("data_clean/model_output/summaries_",model_type,"/atlas_BCR_change_estimates.rds"))

# ------------------------------------------------------------
# Load and clean BBS summaries
# ------------------------------------------------------------

BBS_trend_summary <- read.csv(
  "../../Data/bbs_trend_estimates/Ontario_Atlas_trends_means_by_atlas_period_2002_2022.csv"
) %>%
  rename(
    BBS_q025 = percent_change_q_0.025,
    BBS_mean = percent_change,
    BBS_q975 = percent_change_q_0.975,
    n_BBS_routes_mean = mean_n_routes
  ) %>%
  select(
    english_name, region, n_routes, n_BBS_routes_mean,
    BBS_mean, BBS_q025, BBS_q975
  )

subset(BBS_trend_summary, english_name == "Bank Swallow")

# ------------------------------------------------------------
# Combine atlas and BBS
# ------------------------------------------------------------

comparison <- atlas_change_summary %>%
  left_join(BBS_trend_summary, by = c("english_name", "region")) %>%
  filter(
    !is.na(BBS_mean),
    !is.na(BBS_q025),
    !is.na(BBS_q975)
  ) %>%
  mutate(
    region = factor(region, levels = rev(plot_regions))
  ) %>%
  filter(english_name != c("Bank Swallow"))

# ------------------------------------------------------------
# Transform to log-multiplicative change scale
# ------------------------------------------------------------

comparison_log <- comparison %>%
  mutate(
    atlas_median_log = pct_to_log(atlas_median),
    atlas_q025_log   = pct_to_log(atlas_q025),
    atlas_q975_log   = pct_to_log(atlas_q975),
    atlas_CIW        = atlas_q975 - atlas_q025,
    
    BBS_mean_log = pct_to_log(BBS_mean),
    BBS_q025_log = pct_to_log(BBS_q025),
    BBS_q975_log = pct_to_log(BBS_q975),
    BBS_CIW      = BBS_q975 - BBS_q025
  )

# ------------------------------------------------------------
# Precision weights and filtered comparison set
# ------------------------------------------------------------

zval <- qnorm(0.975)

comparison_log <- comparison_log %>%
  mutate(
    atlas_sd  = (atlas_q975_log - atlas_q025_log) / (2 * zval),
    bbs_sd    = (BBS_q975_log   - BBS_q025_log)   / (2 * zval),
    diff_mean = atlas_median_log - BBS_mean_log,
    diff_sd   = sqrt(atlas_sd^2 + bbs_sd^2),
    weight    = 1 / (atlas_sd^2 + bbs_sd^2),
    diff_CIW  = BBS_CIW - atlas_CIW
  )

# ------------------------------------------------------------
# Filter species for further consideration
# ------------------------------------------------------------

# Must have at least 25 squares per atlas to be included in the figure
min_sq_per_atlas <- 25
comparison_plot <- comparison_log %>%
  filter(
    n_sq_OBBA2 >= min_sq_per_atlas,
    n_sq_OBBA3 >= min_sq_per_atlas
  )

# ------------------------------------------------------------
# Identify species for which counts are potentially too variable for Poisson
# ------------------------------------------------------------
nb_species <- subset(dat$species_to_model, error_family == "nbinomial")

comparison_plot$consider_nbinomial <- "No"
comparison_plot$consider_nbinomial[comparison_plot$english_name %in% nb_species$english_name] <- "Yes"

# ------------------------------------------------------------
# Summaries by BCR
# ------------------------------------------------------------

bcr_summaries <- comparison_plot %>%
  group_by(BCR, region) %>%
  summarise(
    n_species = n(),
    perc_higher = mean(atlas_median_log > BBS_mean_log, na.rm = TRUE) * 100,
    mean_diff = weighted.mean(diff_mean, weight, na.rm = TRUE),
    correlation = weighted_cor(atlas_median_log, BBS_mean_log, weight),
    mean_diff_CIW = mean(diff_CIW, na.rm = TRUE),
    median_diff_CIW = median(diff_CIW, na.rm = TRUE),
    perc_atlas_more_precise = mean(diff_CIW > 0, na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  mutate(cor_label = paste0("r = ", round(correlation, 2)))

# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------

lim <- range(
  comparison_plot[, c("BBS_q025_log", "BBS_q975_log", "atlas_q025_log", "atlas_q975_log")],
  na.rm = TRUE
)

fig3 <- ggplot(comparison_plot) +
  geom_abline(slope = 1, col = "gray80") +
  
  # geom_text_repel(
  #   aes(
  #     x = BBS_mean_log,
  #     y = atlas_median_log,
  #     label = english_name
  #   ),
  #   size = 2,
  #   max.overlaps = 50,
  #   col = "gray50"
  # ) +
  geom_text(
    data = bcr_summaries,
    aes(x = -Inf, y = Inf, label = cor_label),
    hjust = -0.25,
    vjust = 2,
    inherit.aes = FALSE
  ) +

  # geom_errorbar(
  #   aes(
  #     x = BBS_mean_log,
  #     ymin = atlas_q025_log,
  #     ymax = atlas_q975_log,
  #     col = consider_nbinomial
  #   ),
  #   width = 0
  # ) +
  # geom_errorbarh(
  #   aes(
  #     xmin = BBS_q025_log,
  #     xmax = BBS_q975_log,
  #     y = atlas_median_log,
  #     col = consider_nbinomial
  #   ),
  #   height = 0
  # ) +
  geom_point(
    aes(
      x = BBS_mean_log,
      y = atlas_median_log,
      col = consider_nbinomial
    ),
  ) +
  facet_grid(region ~ .) +
  scale_x_continuous(
    name = "BBS overall percent change from 2002 to 2022",
    labels = function(x) sprintf("%+d%%", round(log_to_pct(x)))
  ) +
  scale_y_continuous(
    name = "Atlas overall percent change from Atlas 2 to 3",
    labels = function(x) sprintf("%+d%%", round(log_to_pct(x)))
  ) +
  scale_alpha_continuous(
    range = c(0.2, 1),
    guide = "none",
    trans = "log"
  ) +
  scale_color_manual(
    values = c("black","red"), guide = "none"
  ) +
  
  coord_cartesian(xlim = lim, ylim = lim) +
  theme_bw()

fig3
