# ============================================================
# sampling_change_summaries_by_BCR.R
#
# Purpose:
#   Generate within-BCR summaries of how sampling changed between
#   OBBA2 and OBBA3:
#     (A) Effort/coverage summaries (survey counts, unique days, etc.)
#     (B) Method mix (ARU vs Point_Count) + on_road/on_river
#     (C) Temporal footprint (DayOfYear, Hours_Since_Sunrise)
#     (D) Covariate-space novelty metric (distance to nearest surveyed covariates)
#
# Inputs (in memory or read from disk):
#   - all_surveys : sf with survey points & attributes (Atlas, BCR, ARU, etc.)
#   - grid_OBBA2  : sf grid with pixel_id, BCR, covariates for Atlas 2
#   - grid_OBBA3  : sf grid with pixel_id, BCR, covariates for Atlas 3
#
# Outputs:
#   - out_dir/tables/*.csv and *.rds
#   - out_dir/figures/*.png
#   - out_dir/maps/*.rds  (covariate novelty by pixel; optional gpkg)
# ============================================================

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(RANN)
  library(readr)
  library(stringr)
  library(viridis)
  library(grid)
  library(gridBase)
})

source("R/functions/figure_utils.R")
source("R/functions/spatial_utils.R")

# ---------------------------
# Helpers
# ---------------------------

as_numeric_df <- function(df, cols) {
  df %>% mutate(across(all_of(cols), ~ suppressWarnings(as.numeric(.x))))
}

# helper: draw a ggplot into the current base graphics panel
draw_ggplot_in_base_panel <- function(p) {
  plot.new()
  vp <- gridBase::baseViewports()
  pushViewport(vp$figure)
  grid.draw(ggplotGrob(p))
  popViewport()
}

# ---------------------------
# User settings
# ---------------------------

in_data <- "data_clean/birds/data_ready_for_analysis.rds"
dat <- readRDS(in_data)

all_surveys <- dat$all_surveys
counts      <- dat$counts
grid_OBBA2  <- dat$grid_OBBA2
grid_OBBA3  <- dat$grid_OBBA3
study_boundary <- dat$study_boundary %>% st_as_sf()
species_to_model <- dat$species_to_model

# Temporal variables (names in all_surveys)
doy_var  <- "DayOfYear"
hsr_var  <- "Hours_Since_Sunrise"  # already in your data

# # Create a "Method" variable that is robust to how Survey_Type is recorded
# # Priority: ARU flag; else use Survey_Type
# all_surveys <- all_surveys %>%
#   mutate(
#     Method = case_when(
#       .data[[aru_var]] == 1 ~ "ARU",
#       TRUE ~ "Human"
#     ),
#     Method = factor(Method, levels = c("Human", "ARU"))
#   )

# Method variables
aru_var  <- "ARU"          # 1/0
stype_var <- "Survey_Type" 
atlas_var <- "Project_Name"       # "OBBA2", "OBBA3"
bcr_var   <- "BCR"

# Output dirs
out_dir <- "data_clean/sampling_diagnostics"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "maps"), recursive = TRUE, showWarnings = FALSE)

# Optional: cap for novelty maps
cap_quantile <- 0.99

# Figure sizes
fig_w <- 10
fig_h <- 6
fig_dpi <- 300

# ============================================================
# Map of survey coverage
# ============================================================

svy_OBBA2 <- subset(all_surveys, Project_Name == "OBBA2")
svy_OBBA3 <- subset(all_surveys, Project_Name == "OBBA3")

atlas2_plot <- ggplot() +
  geom_sf(data = study_boundary, fill = "gray90")+
  theme_map() +
  geom_sf(data = subset(all_surveys, Project_Name == "OBBA2"), size = 0.01, shape = 16)+
  ggspatial::annotation_scale(location = "br", width_hint = 0.3) +
  ggspatial::annotation_north_arrow(
    location = "tr",
    which_north = "true",
    pad_x = unit(0.2, "in"),
    pad_y = unit(0.2, "in"),
    style = ggspatial::north_arrow_fancy_orienteering()
  )

atlas3_plot <- ggplot() +
  geom_sf(data = study_boundary, fill = "gray90")+
  theme_map() +
  geom_sf(data = subset(all_surveys, Project_Name == "OBBA3" & Survey_Type == "ARU"), size = 0.01, col = "red", shape = 16)+
  geom_sf(data = subset(all_surveys, Project_Name == "OBBA3" & Survey_Type == "Point_Count"), size = 0.05, shape = 16)+
  ggspatial::annotation_scale(location = "br", width_hint = 0.3) +
  ggspatial::annotation_north_arrow(
    location = "tr",
    which_north = "true",
    pad_x = unit(0.2, "in"),
    pad_y = unit(0.2, "in"),
    style = ggspatial::north_arrow_fancy_orienteering()
  )

# print(atlas2_plot)
# print(atlas3_plot)

ggsave(file.path(out_dir, "figures/map_atlas2_surveys.png"),
       atlas2_plot, width = 10, height = 10, dpi = fig_dpi)

ggsave(file.path(out_dir, "figures/map_atlas3_surveys.png"),
       atlas3_plot, width = 10, height = 10, dpi = fig_dpi)


# ============================================================
# Summaries of covariate coverage within each BCR
# Highlight landscape conditions that were poorly sampled during surveys — 
# they fall outside the most extreme ~5% of the survey distribution.
# ============================================================

covar_names <- readxl::read_xlsx("data_clean/metadata/covariate_names.xlsx")

# Prepare raster template for plotting maps
r_template <- rast(vect(grid_OBBA2), res = 1001)

# Covariates to *potentially* include; per-species covariate variance screening can happen later
covars_to_evaluate <- c(
  "insect_broadleaf",
  "insect_needleleaf",
   "on_river",
   "on_road",
  "urban_2",
   "urban_3",
   "lc_1","lc_4","lc_5",
   "lc_8","lc_8",
   "lc_9","lc_10",
   "lc_11","lc_12",
   "lc_14","lc_17",
  "prec","tmax"
)

# Sample size to flag:
# If density is 0.025 birds per ha, how many surveys would be required for 1 detection within a
# 250m point count radius, assuming perfect detection: 
# n = 1/(D*A*qp) where qp = 1
# n_flag = 1/(0.025 * pi*0.25^2) # about 200
n_flag <- 200

# Loop through covariates
for (lc_to_eval in covars_to_evaluate[12:21]){
  print(lc_to_eval)
  cname <- subset(covar_names,covariate_code==lc_to_eval)$covariate_name
  
  if (lc_to_eval %in% c("prec","tmax","insect_broadleaf","insect_needleleaf")){
    grid_OBBA2$region <- "Province"
    svy_OBBA2$region <- "Province"
    
    grid_OBBA3$region <- "Province"
    svy_OBBA3$region <- "Province"
    
  } else{
    # Support will be evaluated separately in north and southern regions of study area for land cover covariates
    grid_OBBA2$region <- ifelse(grid_OBBA2$BCR %in% c(13,12),"South","North")
    svy_OBBA2$region <- ifelse(svy_OBBA2$BCR %in% c(13,12),"South","North")
    
    grid_OBBA3$region <- ifelse(grid_OBBA3$BCR %in% c(13,12),"South","North")
    svy_OBBA3$region <- ifelse(svy_OBBA3$BCR %in% c(13,12),"South","North")
  }
  
  
  hist_OBBA2 <- plot_covariate_hist_overlay(
    grid_OBBA2, svy_OBBA2, lc_to_eval,
    region_col = "region",
    bins = 30,                 # or set binwidth = 0.5, etc.
    survey_width_frac = 0.5,    # red bars half-width,
    landscape_fill = "black",
    survey_fill = "red",
    title = paste0("Atlas 2 survey (red) vs landscape (black)")
  )
  
  hist_OBBA3 <- plot_covariate_hist_overlay(
    grid_OBBA3, svy_OBBA3, lc_to_eval,
    region_col = "region",
    bins = 30,                 # or set binwidth = 0.5, etc.
    survey_width_frac = 0.5,    # red bars half-width,
    landscape_fill = "black",
    survey_fill = "red",
    title = paste0("Atlas 3 survey (red) vs landscape (black)")
  )
  
  # Prepare rasters of the raw covariate value in each atlas
  rast_OBBA2 <- rasterize(vect(grid_OBBA2), r_template, field = lc_to_eval, fun = "mean")
  rast_OBBA3 <- rasterize(vect(grid_OBBA3), r_template, field = lc_to_eval, fun = "mean")
  range_OBBA2 <- global(rast_OBBA2, fun = "range", na.rm = TRUE)
  range_OBBA3 <- global(rast_OBBA3, fun = "range", na.rm = TRUE)
  common_min <- min(range_OBBA2[1], range_OBBA3[1])
  common_max <- max(range_OBBA2[2], range_OBBA3[2])
  common_range <- c(common_min, common_max)
  
  # --- Evaluate support in Atlas 2
  support_vals <- covariate_support_by_bcr(
    grid_sf   = grid_OBBA2,
    survey_sf = svy_OBBA2,
    covariate = lc_to_eval,
    bcr_col   = "region"
  )
  
  grid_OBBA2$support <- support_vals$support
  grid_OBBA2$q <- support_vals$q
  grid_OBBA2$n_surveys <- support_vals$n_surveys
  grid_OBBA2$n_tail <- grid_OBBA2$n_surveys *grid_OBBA2$support / 2
  grid_OBBA2$flag_fewK <- grid_OBBA2$n_tail < n_flag
  
  # --- Evaluate support in Atlas 3
  
  support_vals <- covariate_support_by_bcr(
    grid_sf   = grid_OBBA3,
    survey_sf = svy_OBBA3,
    covariate = lc_to_eval,
    bcr_col   = "region"
  )
  grid_OBBA3$support <- support_vals$support
  grid_OBBA3$q <- support_vals$q
  grid_OBBA3$n_surveys <- support_vals$n_surveys
  grid_OBBA3$n_tail <- grid_OBBA3$n_surveys *grid_OBBA3$support / 2
  grid_OBBA3$flag_fewK <- grid_OBBA3$n_tail < n_flag
 
  # --- Rasterize support metrics
  r_support_OBBA2 <- rasterize(vect(grid_OBBA2), r_template, field = "n_tail", fun = "mean")
  r_support_OBBA3 <- rasterize(vect(grid_OBBA3), r_template, field = "n_tail", fun = "mean")
  r_flag_OBBA2 <- rasterize(vect(grid_OBBA2), r_template, field = "flag_fewK", fun = "max")
  r_flag_OBBA3 <- rasterize(vect(grid_OBBA3), r_template, field = "flag_fewK", fun = "max")
  
  # ---- Generate plots
  pdf(
    file = paste0(out_dir, "/maps/covariate_support_", lc_to_eval, ".pdf"),
    width = 11,
    height = 14
  )
  
  layout(
    matrix(c(
      1, 1,   # title spans both columns
      2, 3,   # raw maps
      4, 5,   # histograms
      6, 7    # support maps
    ), nrow = 4, byrow = TRUE),
    heights = c(0.3, 1.2, 0.8, 1.2)
  )
  
  # Title row
  par(mar = c(0,0,0,0))
  plot.new()
  text(
    0.5, 0.5,
    paste0(cname,"\n", 
           "covariate code = ",lc_to_eval),
    cex = 1.8,
    font = 2
  )
  
  # Raw covariate maps
  par(mar = c(3,3,3,5))
  plot(rast_OBBA2, col = viridis(10), range = common_range,
       main = paste0("Atlas 2: ", lc_to_eval), legend = FALSE)
  
  par(mar = c(3,3,3,5))
  plot(rast_OBBA3, col = viridis(10), range = common_range,
       main = paste0("Atlas 3: ", lc_to_eval), legend = TRUE)
  
  # covariate histograms
  par(mar = c(2,2,2,2))
  draw_ggplot_in_base_panel(hist_OBBA2)
  
  par(mar = c(2,2,2,2))
  draw_ggplot_in_base_panel(hist_OBBA3)
  
  # Support maps
  max_scale <- max(c(values(r_support_OBBA2),values(r_support_OBBA3)),na.rm = TRUE)
  par(mar = c(3,3,3,5))
  plot(r_support_OBBA2, col = gray.colors(100), range = c(0,max_scale),
       main = paste0("Atlas 2 support: ", lc_to_eval), legend = TRUE)
  
  # Must have at least one "flagged" cell to plot "flag" overlay
  r_flag_OBBA2[r_flag_OBBA2 == 0] <- NA
  n_notna2 <- terra::global(!is.na(r_flag_OBBA2), "sum", na.rm = TRUE)[1,1]
  if (is.finite(n_notna2) && n_notna2 > 0) {
    plot(r_flag_OBBA2, col = "blue", add = TRUE, legend = FALSE)
  }
  
  par(mar = c(3,3,3,5))
  plot(r_support_OBBA3, col = gray.colors(100),range=c(0,max_scale),
       main = paste0("Atlas 3 support: ", lc_to_eval), legend = TRUE)
  
  # Must have at least one "flagged" cell to plot overlay
  r_flag_OBBA3[r_flag_OBBA3 == 0] <- NA
  n_notna3 <- terra::global(!is.na(r_flag_OBBA3), "sum", na.rm = TRUE)[1,1]
  if (is.finite(n_notna3) && n_notna3 > 0) {
    plot(r_flag_OBBA3, col = "blue", add = TRUE, legend = FALSE)
  }
  
  dev.off()
  
}

# ============================================================
# A) Effort / coverage summaries (by BCR x Atlas)
# ============================================================

effort_tbl <- all_surveys %>%
  st_drop_geometry() %>%
  group_by(BCR = .data[[bcr_var]], Atlas = .data[[atlas_var]]) %>%
  summarise(
    n_surveys = n(),
    n_pixels_sampled = n_distinct(pixel_id),
    n_squares = n_distinct(square_id),
    n_unique_days = n_distinct(as.Date(Date_Time)),
    total_minutes = sum(Survey_Duration_Minutes, na.rm = TRUE),
    median_minutes = median(Survey_Duration_Minutes, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(BCR, Atlas)

effort_wide <- effort_tbl %>%
  pivot_wider(
    names_from = Atlas,
    values_from = c(n_surveys, n_pixels_sampled, n_squares, n_unique_days, total_minutes, median_minutes)
  ) %>%
  mutate(
    surveys_ratio = n_surveys_OBBA3 / pmax(n_surveys_OBBA2, 1),
    pixels_ratio  = n_pixels_sampled_OBBA3 / pmax(n_pixels_sampled_OBBA2, 1),
    minutes_ratio = total_minutes_OBBA3 / pmax(total_minutes_OBBA2, 1)
  )

saveRDS(effort_tbl,  file.path(out_dir, "tables/effort_by_BCR_atlas.rds"))
write_csv(effort_tbl, file.path(out_dir, "tables/effort_by_BCR_atlas.csv"))
saveRDS(effort_wide, file.path(out_dir, "tables/effort_by_BCR_change.rds"))
write_csv(effort_wide, file.path(out_dir, "tables/effort_by_BCR_change.csv"))

# Figure: effort bars
p_effort <- effort_tbl %>%
  ggplot(aes(x = factor(BCR), y = n_surveys, fill = Atlas)) +
  geom_col(position = "dodge") +
  labs(x = "BCR", y = "Number of surveys", title = "Survey volume by BCR and Atlas") +
  theme_minimal()

ggsave(file.path(out_dir, "figures/effort_n_surveys_by_BCR.png"),
       p_effort, width = fig_w, height = fig_h, dpi = fig_dpi)

# ============================================================
# B) Method mix by BCR x Atlas
# ============================================================

method_tbl <- all_surveys %>%
  sf::st_drop_geometry() %>%
  dplyr::group_by(BCR = .data[[bcr_var]], Atlas = .data[[atlas_var]], Method) %>%
  dplyr::summarise(
    n = dplyr::n(),
    pct_on_road  = mean(on_road == 1, na.rm = TRUE),
    pct_on_river = mean(on_river == 1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::group_by(BCR, Atlas) %>%
  dplyr::mutate(pct = n / sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(BCR, Atlas, Method)

saveRDS(method_tbl, file.path(out_dir, "tables/method_mix_by_BCR_atlas.rds"))
write_csv(method_tbl, file.path(out_dir, "tables/method_mix_by_BCR_atlas.csv"))

p_method <- method_tbl %>%
  ggplot(aes(x = factor(BCR), y = pct, fill = Method)) +
  geom_col(position = "stack") +
  facet_wrap(~Atlas, nrow = 1) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "BCR", y = "Proportion of surveys", title = "Method composition (Human vs ARU) by BCR") +
  theme_minimal()

ggsave(file.path(out_dir, "figures/method_mix_by_BCR.png"),
       p_method, width = fig_w, height = fig_h, dpi = fig_dpi)

# ============================================================
# C) Temporal footprint (DayOfYear & Hours_Since_Sunrise)
#    - overall by BCR x Atlas
#    - and stratified by Method
# ============================================================

temporal_tbl <- all_surveys %>%
  st_drop_geometry() %>%
  group_by(BCR = .data[[bcr_var]], Atlas = .data[[atlas_var]]) %>%
  summarise(
    n = n(),
    doy_p05 = quantile(.data[[doy_var]], 0.05, na.rm = TRUE),
    doy_p50 = quantile(.data[[doy_var]], 0.50, na.rm = TRUE),
    doy_p95 = quantile(.data[[doy_var]], 0.95, na.rm = TRUE),
    hsr_p05 = quantile(.data[[hsr_var]], 0.05, na.rm = TRUE),
    hsr_p50 = quantile(.data[[hsr_var]], 0.50, na.rm = TRUE),
    hsr_p95 = quantile(.data[[hsr_var]], 0.95, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(BCR, Atlas)

saveRDS(temporal_tbl, file.path(out_dir, "tables/temporal_summary_by_BCR_atlas.rds"))
write_csv(temporal_tbl, file.path(out_dir, "tables/temporal_summary_by_BCR_atlas.csv"))

# Plots: DOY density by BCR
p_doy <- all_surveys %>%
  st_drop_geometry() %>%
  ggplot(aes(x = .data[[doy_var]], group = Atlas, linetype = Atlas)) +
  geom_density(adjust = 1) +
  facet_wrap(~BCR, scales = "free_y") +
  labs(x = "Day of Year", y = "Density", title = "Day-of-year distribution by BCR (Atlas 2 vs 3)") +
  theme_minimal()

ggsave(file.path(out_dir, "figures/temporal_DOY_density_by_BCR.png"),
       p_doy, width = fig_w, height = fig_h, dpi = fig_dpi)

# Plots: Hours since sunrise density by BCR
p_hsr <- all_surveys %>%
  st_drop_geometry() %>%
  ggplot(aes(x = .data[[hsr_var]], group = Atlas, linetype = Atlas)) +
  geom_density(adjust = 1) +
  facet_wrap(~BCR, scales = "free_y") +
  labs(x = "Hours Since Sunrise", y = "Density",
       title = "Time-of-day distribution by BCR (Atlas 2 vs 3)") +
  theme_minimal()

ggsave(file.path(out_dir, "figures/temporal_HSR_density_by_BCR.png"),
       p_hsr, width = fig_w, height = fig_h, dpi = fig_dpi)

# Stratified by Method (Human vs ARU)
p_doy_method <- all_surveys %>%
  st_drop_geometry() %>%
  ggplot(aes(x = .data[[doy_var]], group = Atlas, linetype = Atlas)) +
  geom_density(adjust = 1) +
  facet_grid(Method ~ BCR, scales = "free_y") +
  labs(x = "Day of Year", y = "Density",
       title = "Day-of-year distributions by BCR, stratified by method") +
  theme_minimal()

ggsave(file.path(out_dir, "figures/temporal_DOY_density_by_BCR_byMethod.png"),
       p_doy_method, width = 14, height = 7, dpi = fig_dpi)

p_hsr_method <- all_surveys %>%
  st_drop_geometry() %>%
  ggplot(aes(x = .data[[hsr_var]], group = Atlas, linetype = Atlas)) +
  geom_density(adjust = 1) +
  facet_grid(Method ~ BCR, scales = "free_y") +
  labs(x = "Hours Since Sunrise", y = "Density",
       title = "Time-of-day distributions by BCR, stratified by method") +
  theme_minimal()

ggsave(file.path(out_dir, "figures/temporal_HSR_density_by_BCR_byMethod.png"),
       p_hsr_method, width = 14, height = 7, dpi = fig_dpi)

# 2D heatmap counts (TOD x DOY) per BCR x Atlas
# (coarse binning keeps figure reasonable)
bin_doy <- 5
bin_hsr <- 0.5

tod_doy_tbl <- all_surveys %>%
  st_drop_geometry() %>%
  mutate(
    doy_bin = floor(.data[[doy_var]] / bin_doy) * bin_doy,
    hsr_bin = floor(.data[[hsr_var]] / bin_hsr) * bin_hsr
  ) %>%
  group_by(BCR = .data[[bcr_var]], Atlas = .data[[atlas_var]], Method, doy_bin, hsr_bin) %>%
  summarise(n = n(), .groups = "drop")

saveRDS(tod_doy_tbl, file.path(out_dir, "tables/tod_doy_binned_counts.rds"))

p_tod_doy <- tod_doy_tbl %>%
  ggplot(aes(x = doy_bin, y = hsr_bin, fill = n)) +
  geom_raster() +
  facet_grid(Method ~ BCR + Atlas, scales = "free") +
  labs(x = "Day of Year (binned)", y = "Hours Since Sunrise (binned)",
       title = "Survey timing footprint (counts) by BCR × Atlas × Method") +
  scale_fill_gradientn(colours = viridis(10))+
  theme_minimal()

ggsave(file.path(out_dir, "figures/temporal_TODxDOY_heatmap.png"),
       p_tod_doy, width = 16, height = 7, dpi = fig_dpi)

# ============================================================
# D) Covariate-space novelty metric (by BCR)
# ============================================================

compute_novelty_for_one_bcr <- function(bcr_id) {
  s_bcr  <- all_surveys %>% filter(.data[[bcr_var]] == bcr_id)
  g2_bcr <- grid_OBBA2  %>% filter(.data[[bcr_var]] == bcr_id)
  g3_bcr <- grid_OBBA3  %>% filter(.data[[bcr_var]] == bcr_id)
  
  s2 <- s_bcr %>% filter(.data[[atlas_var]] == "OBBA2")
  s3 <- s_bcr %>% filter(.data[[atlas_var]] == "OBBA3")
  
  scaler <- compute_scaler_from_grid(g2_bcr, g3_bcr, base_covars)
  
  Xg2 <- scale_matrix(st_drop_geometry(g2_bcr), base_covars, scaler)
  Xg3 <- scale_matrix(st_drop_geometry(g3_bcr), base_covars, scaler)
  Xs2 <- scale_matrix(st_drop_geometry(s2), base_covars, scaler)
  Xs3 <- scale_matrix(st_drop_geometry(s3), base_covars, scaler)
  
  d2 <- nearest_distance(Xg2, Xs2)
  d3 <- nearest_distance(Xg3, Xs3)
  
  novelty_df <- tibble(
    pixel_id = g3_bcr$pixel_id,
    BCR = bcr_id
  ) %>%
    left_join(tibble(pixel_id = g2_bcr$pixel_id, d_OBBA2 = d2), by = "pixel_id") %>%
    mutate(d_OBBA3 = d3,
           delta_d = d_OBBA2 - d_OBBA3)
  
  novelty_sf <- g3_bcr %>%
    select(pixel_id, BCR) %>%
    left_join(novelty_df, by = c("pixel_id", "BCR"))
  
  # capped columns for mapping
  cap2 <- quantile(novelty_sf$d_OBBA2, cap_quantile, na.rm = TRUE, names = FALSE)
  cap3 <- quantile(novelty_sf$d_OBBA3, cap_quantile, na.rm = TRUE, names = FALSE)
  capd <- quantile(abs(novelty_sf$delta_d), cap_quantile, na.rm = TRUE, names = FALSE)
  
  novelty_sf <- novelty_sf %>%
    mutate(
      d_OBBA2_capped = pmin(d_OBBA2, cap2),
      d_OBBA3_capped = pmin(d_OBBA3, cap3),
      delta_d_capped = pmax(pmin(delta_d, capd), -capd)
    )
  
  sum_tbl <- tibble(BCR = bcr_id) %>%
    bind_cols(
      summarize_dist_vec(novelty_sf$d_OBBA2) %>% rename_with(~ paste0("d2_", .x)),
      summarize_dist_vec(novelty_sf$d_OBBA3) %>% rename_with(~ paste0("d3_", .x)),
      summarize_dist_vec(novelty_sf$delta_d) %>% rename_with(~ paste0("dd_", .x))
    ) %>%
    mutate(
      pct_pixels_improved = mean(novelty_sf$delta_d > 0, na.rm = TRUE),
      pct_pixels_worse    = mean(novelty_sf$delta_d < 0, na.rm = TRUE),
      n_surveys_obba2 = nrow(s2),
      n_surveys_obba3 = nrow(s3)
    )
  
  list(novelty_sf = novelty_sf, summary_tbl = sum_tbl)
}

nov_results <- purrr::map(bcrs, compute_novelty_for_one_bcr)
novelty_by_pixel <- purrr::map_dfr(nov_results, "novelty_sf") %>% st_as_sf()
novelty_summary_by_BCR <- purrr::map_dfr(nov_results, "summary_tbl")

saveRDS(novelty_by_pixel, file.path(out_dir, "maps/sampling_novelty_by_pixel.rds"))
saveRDS(novelty_summary_by_BCR, file.path(out_dir, "tables/sampling_novelty_summary_by_BCR.rds"))
write_csv(novelty_summary_by_BCR, file.path(out_dir, "tables/sampling_novelty_summary_by_BCR.csv"))

# Figure: novelty distributions by BCR
lims <- quantile(novelty_by_pixel$delta_d,
                 c(0.01, 0.99),
                 na.rm = TRUE)

p_dd <- ggplot(novelty_by_pixel %>% sf::st_drop_geometry(),
               aes(x = delta_d)) +
  geom_histogram(bins = 5000, color = "white") +
  coord_cartesian(xlim = lims) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    x = expression(Delta~"covariate support (OBBA2 – OBBA3)"),
    y = "Number of pixels",
    title = "Change in covariate-space support (central 98%)",
    subtitle = "Positive values indicate better covariate support in atlas 3\nTails truncated for interpretability; full distribution summarized in text"
  ) +
  theme_minimal()+
  facet_wrap(BCR~., scales = "free_y")
p_dd

ggsave(file.path(out_dir, "figures/novelty_delta_hist_by_BCR.png"),
       p_dd, width = 12, height = 6, dpi = fig_dpi)


# ---------------------------
# Univariate covariate coverage
# ---------------------------

compute_decile_coverage <- function(grid_vals, survey_vals, n_bins = 10) {
  # Remove NAs
  grid_vals   <- grid_vals[is.finite(grid_vals)]
  survey_vals <- survey_vals[is.finite(survey_vals)]
  
  if (length(grid_vals) < n_bins || length(survey_vals) == 0) {
    return(NULL)
  }
  
  # Atlas-specific decile breaks from GRID
  probs  <- seq(0, 1, length.out = n_bins + 1)
  breaks <- quantile(grid_vals, probs = probs, na.rm = TRUE, names = FALSE)
  
  # Guard against duplicate breaks
  breaks <- unique(breaks)
  if (length(breaks) <= 2) return(NULL)
  
  # Bin surveys
  bins <- cut(
    survey_vals,
    breaks = breaks,
    include.lowest = TRUE,
    labels = FALSE
  )
  
  p <- tabulate(bins, nbins = length(breaks) - 1)
  p <- p / sum(p)
  
  # Ideal distribution
  p0 <- rep(1 / length(p), length(p))
  
  tibble(
    max_dev = max(abs(p - p0)),
    rmse    = sqrt(mean((p - p0)^2)),
    entropy = {
      pe <- p[p > 0]
      H  <- -sum(pe * log(pe))
      H / log(length(p))  # normalized 0–1
    }
  )
}

univariate_coverage_tbl <- purrr::map_dfr(base_covars, function(v) {
  
  purrr::map_dfr(bcrs, function(b) {
    
    # --- OBBA2 ---
    g2 <- grid_OBBA2 %>%
      filter(BCR == b) %>%
      st_drop_geometry()
    
    s2 <- all_surveys %>%
      filter(BCR == b, Atlas == "OBBA2") %>%
      st_drop_geometry()
    
    res2 <- compute_decile_coverage(g2[[v]], s2[[v]])
    
    # --- OBBA3 ---
    g3 <- grid_OBBA3 %>%
      filter(BCR == b) %>%
      st_drop_geometry()
    
    s3 <- all_surveys %>%
      filter(BCR == b, Atlas == "OBBA3") %>%
      st_drop_geometry()
    
    res3 <- compute_decile_coverage(g3[[v]], s3[[v]])
    
    bind_rows(
      if (!is.null(res2)) res2 %>% mutate(BCR = b, Atlas = "OBBA2", covariate = v),
      if (!is.null(res3)) res3 %>% mutate(BCR = b, Atlas = "OBBA3", covariate = v)
    )
  })
})

# Order covariates by average "badness" to make patterns pop
cov_order <- univariate_coverage_tbl %>%
  group_by(covariate) %>%
  summarise(avg_rmse = mean(rmse, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(avg_rmse)) %>%
  pull(covariate)

plot_heatmap_metric <- function(metric = c("rmse", "max_dev", "entropy")) {
  metric <- match.arg(metric)
  
  df <- univariate_coverage_tbl %>%
    mutate(
      covariate = factor(covariate, levels = cov_order),
      BCR = factor(BCR),
      Atlas = factor(Atlas, levels = c("OBBA2", "OBBA3"))
    )
  
  ggplot(df, aes(x = BCR, y = covariate, fill = .data[[metric]])) +
    geom_tile(color = "white", linewidth = 0.2) +
    facet_wrap(~ Atlas, nrow = 1) +
    labs(
      x = "BCR",
      y = NULL,
      fill = metric,
      title = paste("Univariate sampling representativeness by covariate (", metric, ")", sep = "")
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_text(size = 9)
    )
}

p_rmse <- plot_heatmap_metric("rmse")
p_rmse

# ---------------------------
# Done
# ---------------------------

message("All sampling-change summaries written to: ", out_dir)
