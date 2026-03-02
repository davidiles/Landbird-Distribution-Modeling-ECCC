# ============================================================
# 10b_generate_crossvalidation_figures with cross-validation flags.R
#
# ============================================================

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(cowplot)
  library(pROC)
})

# ---- utils
source("R/functions/figure_utils.R")

# ------------------------------------------------------------
# Cross-validation overlay settings
# ------------------------------------------------------------
cv_dir <- "data_clean/model_output/cv"
include_cv_flags_default <- TRUE

aggregate_cv_blocks <- function(block_summ) {
  req <- c("sp_english","sp_code","block_size_km","block_id","Atlas",
           "rep","fold","n","mean_obs","bias","rmse","mae","geometry",
           "mean_pred_yrep_q50","mean_pred_yrep_q025","mean_pred_yrep_q975","covered_95",
           "p_any_detect","obs_any_detect")
  stopifnot(all(req %in% names(block_summ)))
  
  block_summ %>%
    dplyr::group_by(sp_english, sp_code, block_size_km, block_id, Atlas) %>%
    dplyr::summarize(
      n_eval = dplyr::n_distinct(interaction(rep, fold, drop = TRUE)),
      n_surveys_total = mean(n, na.rm = TRUE),
      
      mean_obs = mean(mean_obs, na.rm = TRUE),
      
      mean_pred_yrep_q025 = mean(mean_pred_yrep_q025, na.rm = TRUE),
      mean_pred_yrep_q50  = mean(mean_pred_yrep_q50,  na.rm = TRUE),
      mean_pred_yrep_q975 = mean(mean_pred_yrep_q975, na.rm = TRUE),
      mean_pred_yrep_mean = mean(mean_pred_mean, na.rm = TRUE),
      
      covered_95 = mean(covered_95, na.rm = TRUE),
      
      rmse = mean(rmse, na.rm = TRUE),
      mae  = mean(mae,  na.rm = TRUE),
      bias = mean(bias, na.rm = TRUE),
      
      # presence/absence channel
      p_any_detect   = mean(p_any_detect, na.rm = TRUE),         # model probability of ≥1 detection
      obs_any_detect = as.integer(mean(obs_any_detect, na.rm = TRUE) >= 0.5),  # should be constant; this is defensive
      
      geometry = dplyr::first(geometry),
      .groups = "drop"
    ) %>%
    sf::st_as_sf()
}


make_cv_flag_layers <- function(
    block_avg,
    atlas = c("OBBA2","OBBA3"),
    block_size_km = NULL,
    n_min = 25,
    fail_rate_outline = 0.5
) {
  atlas <- match.arg(atlas)
  
  x <- block_avg %>%
    dplyr::filter(Atlas == atlas) %>%
    { if (!is.null(block_size_km)) dplyr::filter(., .data$block_size_km == block_size_km) else . } %>%
    dplyr::filter(n_surveys_total >= n_min) %>%
    dplyr::mutate(
      fail_rate = 1 - covered_95,
      poor_coverage = fail_rate >= fail_rate_outline,
      
      # Direction of bias
      tri_dir = dplyr::case_when(
        bias >  0 ~ "up",
        bias < 0 ~ "down",
        TRUE ~ NA_character_
      )
    )
  
  all_polys  <- x
  poor_polys <- x %>% dplyr::filter(poor_coverage)
  
  # Points for triangles, regardless of coverage
  tri_pts <- sf::st_point_on_surface(poor_polys) %>%
    dplyr::filter(!is.na(tri_dir))
  
  list(
    all_polys  = all_polys,
    poor_polys = poor_polys,
    tri_pts    = tri_pts,
    flagged_tbl = x
  )
}

add_cv_flags_to_plot <- function(
    base_plot,
    all_polys,
    poor_polys,
    tri_pts,
    study_boundary,
    bias_positive_color = "blue",
    bias_negative_color = "red",
    
    grey_lwd = 0.20,
    black_lwd = 0.45,
    tri_size = 2.5,
    tri_stroke = 0.5,
    pad_frac = 0.10
) {
  
  bb <- sf::st_bbox(study_boundary)
  bb <- expand_bbox(bb, pad_frac)
  
  # Split triangles for stable shape mapping
  tri_up   <- tri_pts %>% dplyr::filter(tri_dir == "up")
  tri_down <- tri_pts %>% dplyr::filter(tri_dir == "down")
  
  p <- base_plot +
    
    # 1) Grey borders for all evaluated blocks (draw FIRST)
    ggplot2::geom_sf(
      data = all_polys,
      fill = "black",
      alpha = 0.2,
      colour = "gray75",
      linewidth = grey_lwd,
      inherit.aes = FALSE
    ) +
    
    # 2) Blue borders for blocks where predictions were too high
    ggplot2::geom_sf(
      data = poor_polys %>% subset(bias > 0),
      fill = NA,
      colour = bias_positive_color,
      linewidth = black_lwd,
      inherit.aes = FALSE
    ) +
    
    # 2) Red borders for blocks where predictions were too low
    ggplot2::geom_sf(
      data = poor_polys %>% subset(bias < 0),
      fill = NA,
      colour = bias_negative_color,
      linewidth = black_lwd,
      inherit.aes = FALSE
    ) +
    
    # 3) Triangles for positive bias
    ggplot2::geom_sf(
      data = tri_up,
      shape = 24,               # filled triangle up
      size = tri_size,
      stroke = tri_stroke,
      colour = "black",
      fill = bias_positive_color,
      inherit.aes = FALSE
    ) +
    
    # Triangles for negative bias
    ggplot2::geom_sf(
      data = tri_down,
      shape = 25,               # filled triangle down
      size = tri_size,
      stroke = tri_stroke,
      colour = "black",
      fill = bias_negative_color,
      inherit.aes = FALSE
    )
  
  p + ggplot2::coord_sf(
    crs = sf::st_crs(study_boundary),
    xlim = c(bb["xmin"], bb["xmax"]),
    ylim = c(bb["ymin"], bb["ymax"]),
    expand = FALSE
  )
}

plot_cv_calibration <- function(
    block_avg_atlas,
    atlas_label = NULL,
    log_scale = TRUE,
    point_alpha = 0.8,
    point_size_range = c(1.2, 5),
    n_min = NULL,
    # controls for the in-plot text
    show_metrics = TRUE,
    metrics_digits = 2,
    metrics_x = Inf,
    metrics_y = -Inf,
    metrics_hjust = 1.02,
    metrics_vjust = -0.2,
    metrics_label_size = 3.2
) {
  x <- block_avg_atlas
  
  if (!is.null(n_min)) {
    x <- x %>% dplyr::filter(n_surveys_total >= n_min)
  }
  
  # Safety: drop missing
  x <- x %>%
    dplyr::filter(
      is.finite(mean_obs),
      is.finite(mean_pred_yrep_q50),
      is.finite(mean_pred_yrep_q025),
      is.finite(mean_pred_yrep_q975),
      is.finite(n_surveys_total)
    )
  
  # Pick a subtitle label
  if (is.null(atlas_label) && "Atlas" %in% names(x)) {
    atlas_label <- unique(x$Atlas)
    atlas_label <- if (length(atlas_label) == 1) atlas_label else NULL
  }
  
  # ---- base plot (your original)
  p <- ggplot2::ggplot(x, ggplot2::aes(x = mean_obs, y = mean_pred_yrep_q50)) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = 2, linewidth = 0.4) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = mean_pred_yrep_q025, ymax = mean_pred_yrep_q975),
      width = 0,
      alpha = 0.35,
      linewidth = 0.25
    ) +
    ggplot2::geom_point(
      ggplot2::aes(size = n_surveys_total, colour = bias),
      alpha = point_alpha
    ) +
    ggplot2::scale_size_continuous(range = point_size_range, guide = "none") +
    ggplot2::scale_colour_gradient2(
      midpoint = 0,
      low = "red",
      mid = "grey85",
      high = "blue",
      name = "Bias\n(pred − obs)"
    ) +
    ggplot2::labs(
      title = "Obs vs predicted counts per xval block",
      subtitle = atlas_label,
      x = "Observed mean count (block)",
      y = "Predicted median count (block)"
    ) +
    ggplot2::coord_equal() +
    ggplot2::theme_bw()
  
  if (log_scale) {
    p <- p +
      ggplot2::scale_x_continuous(trans = "log1p") +
      ggplot2::scale_y_continuous(trans = "log1p")
  }
  
  # ---- metrics (computed on the exact plotted data)
  if (show_metrics) {
    
    n_blocks <- nrow(x)
    
    # --- Presence/absence metrics for AUC (block-level "any detection")
    auc <- NA_real_
    n_pos <- NA_integer_
    n_neg <- NA_integer_
    
    auc <- NA_real_
    auc_ci_txt <- NULL
    
    if (all(c("p_any_detect","obs_any_detect") %in% names(x))) {
      y_bin <- x$obs_any_detect
      s_bin <- x$p_any_detect
      
      # need both classes
      if (length(unique(y_bin[is.finite(y_bin) & is.finite(s_bin)])) == 2) {
        roc_obj <- tryCatch(
          pROC::roc(response = y_bin, predictor = s_bin, quiet = TRUE),
          error = function(e) NULL
        )
        if (!is.null(roc_obj)) {
          auc <- as.numeric(pROC::auc(roc_obj))
          
          ci <- tryCatch(pROC::ci.auc(roc_obj), error = function(e) NULL)
          ci_num <- as.numeric(ci)
          if (!is.null(ci)) {
            auc_ci_txt <- paste0(" (", formatC(ci_num[1], format="f", digits=2), "–",
                                 formatC(ci_num[3], format="f", digits=2), ")")
          }
        }
      }
    }
    
    # detection context
    n_detect_blocks <- sum(x$mean_obs > 0, na.rm = TRUE)
    p_detect <- if (n_blocks > 0) n_detect_blocks / n_blocks else NA_real_
    p_zero   <- if (n_blocks > 0) sum(x$mean_obs == 0, na.rm = TRUE) / n_blocks else NA_real_
    
    mean_obs <- mean(x$mean_obs, na.rm = TRUE)
    med_obs  <- stats::median(x$mean_obs, na.rm = TRUE)
    
    # errors on count scale 
    err  <- x$mean_pred_yrep_q50 - x$mean_obs
    mae  <- mean(abs(err), na.rm = TRUE)
    rmse <- sqrt(mean(err^2, na.rm = TRUE))
    bias_mean <- mean(err, na.rm = TRUE)
    
    # --- Conditional bias summaries (helpful for rare species)
    is_detect <- x$mean_obs > 0
    is_zero   <- x$mean_obs == 0
    
    # Mean bias for blocks where species was detected
    bias_detect <- if (any(is_detect, na.rm = TRUE)) mean(err[is_detect], na.rm = TRUE) else NA_real_
    
    # Mean bias for blocks where species was not detected
    bias_zero   <- if (any(is_zero,   na.rm = TRUE)) mean(err[is_zero],   na.rm = TRUE) else NA_real_
    
    # optional context (often very informative)
    mean_pred_in_detect <- if (any(is_detect, na.rm = TRUE)) mean(x$mean_pred_yrep_q50[is_detect], na.rm = TRUE) else NA_real_
    mean_obs_in_detect  <- if (any(is_detect, na.rm = TRUE)) mean(x$mean_obs[is_detect], na.rm = TRUE) else NA_real_
    
    mean_pred_in_zero <- if (any(is_zero, na.rm = TRUE)) mean(x$mean_pred_yrep_q50[is_zero], na.rm = TRUE) else NA_real_
    
    # PI coverage on count scale
    cover <- mean(
      x$mean_obs >= x$mean_pred_yrep_q025 & x$mean_obs <= x$mean_pred_yrep_q975,
      na.rm = TRUE
    )
    
    # R2 computed on the same scale as the plot axes (log1p if requested)
    y    <- x$mean_obs
    yhat <- x$mean_pred_yrep_q50
    w    <- x$n_surveys_total
    
    if (log_scale) {
      y    <- log1p(y)
      yhat <- log1p(yhat)
    }
    
    # observed variation (same scale as R²)
    y_var <- stats::var(y, na.rm = TRUE)
    y_sd  <- stats::sd(y, na.rm = TRUE)
    
    # R² (unweighted)
    ss_res <- sum((y - yhat)^2, na.rm = TRUE)
    ss_tot <- sum((y - mean(y, na.rm = TRUE))^2, na.rm = TRUE)
    r2 <- if (is.finite(ss_tot) && ss_tot > 0) 1 - ss_res / ss_tot else NA_real_
    
    # R² (weighted by n_surveys_total)
    ybar_w <- stats::weighted.mean(y, w = w, na.rm = TRUE)
    ss_res_w <- sum(w * (y - yhat)^2, na.rm = TRUE)
    ss_tot_w <- sum(w * (y - ybar_w)^2, na.rm = TRUE)
    r2_w <- if (is.finite(ss_tot_w) && ss_tot_w > 0) 1 - ss_res_w / ss_tot_w else NA_real_
    
    # calibration slope (weighted) on same scale as plot
    slope <- NA_real_
    if (all(is.finite(y), is.finite(yhat), is.finite(w)) && length(y) >= 2) {
      lm_fit <- stats::lm(y ~ yhat, weights = w)
      slope <- unname(stats::coef(lm_fit)[2])
    }
    
    # formatting helper
    fmt <- function(val, digits = metrics_digits) {
      if (!is.finite(val)) return("NA")
      formatC(val, format = "f", digits = digits)
    }
    pct <- function(val, digits = 1) {
      if (!is.finite(val)) return("NA")
      paste0(formatC(100 * val, format = "f", digits = digits), "%")
    }
    
    metrics_text <- paste0(
      "Blocks: ", n_blocks, "\n",
      "Detected blocks: ", n_detect_blocks, " (", pct(p_detect), ")\n",
      "Zero blocks: ", pct(p_zero), "\n",
      "Mean obs: ", fmt(mean_obs), " | Median: ", fmt(med_obs), "\n",
      if (log_scale) {
        paste0("SD[log1p(obs)]: ", fmt(y_sd), " | Var: ", fmt(y_var), "\n")
      } else {
        paste0("SD[obs]: ", fmt(y_sd), " | Var: ", fmt(y_var), "\n")
      },
      "R²: ", fmt(r2), " (wtd: ", fmt(r2_w), ")\n",
      "MAE: ", fmt(mae), " | RMSE: ", fmt(rmse), "\n",
      "Mean bias (pred−obs): ", fmt(bias_mean), "\n",
      "  Bias | detected blocks: ", fmt(bias_detect), "\n",
      "  Bias | zero blocks: ", fmt(bias_zero), "\n",
      "PI coverage: ", pct(cover), "\n",
      "Cal slope (obs~pred): ", fmt(slope), "\n",
      "AUC (any detect): ", fmt(auc), if (!is.null(auc_ci_txt)) auc_ci_txt else ""
    )
    
    # Put the label inside the visible panel (robust for log scales + cowplot)
    x_left <- 0  # safe with log1p
    y_top  <- max(c(x$mean_pred_yrep_q975, x$mean_obs, x$mean_pred_yrep_q50), na.rm = TRUE)
    
    p <- p +
      ggplot2::annotate(
        "label",
        x = x_left,
        y = y_top,
        label = metrics_text,
        hjust = 0,   # left
        vjust = 1,   # top
        size = metrics_label_size,
        fill = "white"
      )
  }
  
  p
}


# ------------------------------------------------------------
# Config
# ------------------------------------------------------------

in_data <- "data_clean/birds/data_ready_for_analysis.rds"
pred_dir <- "data_clean/model_output/predictions"
fig_dir  <- "data_clean/model_output/cv/figures/"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

in_bcr <- "../../Data/Spatial/BCR/BCR_Terrestrial_master.shp"
in_water <- "data_clean/spatial/water_filtered.shp"
in_atlas_squares <- "../../Data/Spatial/National/AtlasSquares/NationalSquares_FINAL.shp"

dpi <- 1000
width_in <- 10
height_in <- 8
ggsave_type <- "cairo"

plot_res <- 1.001
make_atlas_square_overlay <- TRUE

# ------------------------------------------------------------
# Load base data
# ------------------------------------------------------------
stopifnot(file.exists(in_data))
dat <- readRDS(in_data)

all_surveys <- dat$all_surveys
counts      <- dat$counts

grid2 <- dat$grid_OBBA2
grid3 <- dat$grid_OBBA3
study_boundary <- dat$study_boundary %>% sf::st_as_sf()

stopifnot(inherits(grid2, "sf"), inherits(grid3, "sf"))

crs_use <- st_crs(all_surveys)
study_boundary <- st_transform(study_boundary, crs_use)
grid3 <- st_transform(grid3, crs_use)
grid2 <- st_transform(grid2, crs_use)

bcr_sf <- st_read(in_bcr, quiet = TRUE) %>%
  st_make_valid() %>%
  st_transform(crs_use) %>%
  dplyr::filter(PROVINCE_S %in% c("ONTARIO", "ON", "Ontario") | is.na(PROVINCE_S)) %>%
  dplyr::select(BCR, BCRNAME, PROVINCE_S) %>%
  group_by(BCR, BCRNAME) %>%
  summarise(geometry = st_union(geometry), .groups = "drop")

water_sf <- NULL
if (file.exists(in_water)) {
  water_sf <- st_read(in_water, quiet = TRUE) %>%
    st_make_valid() %>%
    st_transform(crs_use)
}

atlas_sq_centroids_all <- NULL
if (make_atlas_square_overlay && file.exists(in_atlas_squares)) {
  atlas_sq <- st_read(in_atlas_squares, quiet = TRUE) %>%
    st_transform(crs_use)
  if ("prov" %in% names(atlas_sq)) atlas_sq <- atlas_sq %>% filter(prov == "ON")
  atlas_sq_centroids_all <- st_centroid(atlas_sq)
}

# ------------------------------------------------------------
# Find prediction files
# ------------------------------------------------------------
pred_files <- list.files(pred_dir, pattern = "\\.rds$", full.names = TRUE)
stopifnot(length(pred_files) > 0)
message("Prediction files found: ", length(pred_files))

# ------------------------------------------------------------
# Loop species
# ------------------------------------------------------------

species_to_plot <- c("White-throated Sparrow","Boreal Chickadee")

for (pf in pred_files) {
  
  preds <- readRDS(pf)
  sp_english <- preds$sp_english
  
  if ((sp_english %in% species_to_plot) == FALSE) next
  
  sp_file <- sp_english %>%
    str_to_lower() %>%
    str_replace_all("[^a-z0-9]+", "_") %>%
    str_replace_all("^_|_$", "")
  
  cv_path_survey_summary <- file.path(cv_dir, "predictions", paste0(sp_file, ".rds"))     # Predictions at individual survey level
  cv_path_block_summary <- file.path(cv_dir, "block_summaries", paste0(sp_file, ".rds"))  # Predictions summarized at block level
  if (!file.exists(cv_path_block_summary)) next
  
  message("Mapping: ", sp_english)
  
  # Shared colour-scale limit across OBBA2 + OBBA3
  zmax2 <- as.numeric(stats::quantile(preds$OBBA2$OBBA2_q50, 0.99, na.rm = TRUE))
  zmax3 <- as.numeric(stats::quantile(preds$OBBA3$OBBA3_q50, 0.99, na.rm = TRUE))
  zmax <- max(zmax2, zmax3, na.rm = TRUE)
  if (!is.finite(zmax) || zmax <= 0) zmax <- 1
  rel_bounds <- list(lower = 0, upper = zmax)
  
  if (is.null(preds$OBBA2) || is.null(preds$OBBA3)) {
    message("  Skipping (missing OBBA2/OBBA3 summaries): ", basename(pf))
    next
  }
  
  maps2 <- make_relabund_maps(
    species_name = sp_english,
    grid_sf = grid2,
    pred_summary = preds$OBBA2,
    prefix = "OBBA2",
    study_boundary = study_boundary,
    bcr_sf = bcr_sf,
    water_sf = NULL,
    atlas_squares_centroids = NULL,
    title = "Atlas 2",
    subtitle = "Relative Abundance",
    res = plot_res,
    bounds = rel_bounds
  )
  
  maps3 <- make_relabund_maps(
    species_name = sp_english,
    grid_sf = grid3,
    pred_summary = preds$OBBA3,
    prefix = "OBBA3",
    study_boundary = study_boundary,
    bcr_sf = bcr_sf,
    water_sf = NULL,
    atlas_squares_centroids = NULL,
    title = "Atlas 3",
    subtitle = "Relative Abundance",
    res = plot_res,
    bounds = rel_bounds
  )
  
  # ------------------------------------------------------------
  # Overlay cross-validation flags
  # ------------------------------------------------------------
  
  cv_block_summ <- readRDS(cv_path_block_summary)
  cv_block_avg  <- aggregate_cv_blocks(cv_block_summ)
  
  block_sizes <- sort(unique(cv_block_avg$block_size_km))
  block_size_use <- if (length(block_sizes) > 0) max(block_sizes) else NULL
  
  cv2 <- make_cv_flag_layers(
    cv_block_avg, atlas = "OBBA2",
    block_size_km = block_size_use,
    n_min = 10, fail_rate_outline = 0.5
  )
  
  cv3 <- make_cv_flag_layers(
    cv_block_avg, atlas = "OBBA3",
    block_size_km = block_size_use,
    n_min = 10, fail_rate_outline = 0.5
  )
  
  # CRS alignment
  cv2$all_polys  <- sf::st_transform(cv2$all_polys,  sf::st_crs(study_boundary))
  cv2$poor_polys <- sf::st_transform(cv2$poor_polys, sf::st_crs(study_boundary))
  cv2$tri_pts    <- sf::st_transform(cv2$tri_pts,    sf::st_crs(study_boundary))
  
  cv3$all_polys  <- sf::st_transform(cv3$all_polys,  sf::st_crs(study_boundary))
  cv3$poor_polys <- sf::st_transform(cv3$poor_polys, sf::st_crs(study_boundary))
  cv3$tri_pts    <- sf::st_transform(cv3$tri_pts,    sf::st_crs(study_boundary))
  
  maps2$xval_plot <- add_cv_flags_to_plot(
    base_plot = maps2$q50_plot,
    all_polys = cv2$all_polys,
    poor_polys = cv2$poor_polys,
    tri_pts = cv2$tri_pts,
    tri_size = 1.5,
    study_boundary = study_boundary
  )
  
  maps3$xval_plot <- add_cv_flags_to_plot(
    base_plot = maps3$q50_plot,
    all_polys = cv3$all_polys,
    poor_polys = cv3$poor_polys,
    tri_pts = cv3$tri_pts,
    tri_size = 1.5,
    study_boundary = study_boundary
  )
  
  # Calibration plots
  cv2_tbl <- cv_block_avg %>% dplyr::filter(Atlas == "OBBA2")
  cv3_tbl <- cv_block_avg %>% dplyr::filter(Atlas == "OBBA3")
  
  cal2 <- plot_cv_calibration(cv2_tbl, atlas_label = "Atlas 2", log_scale = TRUE, n_min = 10,metrics_digits = 3, metrics_label_size = 2)
  cal3 <- plot_cv_calibration(cv3_tbl, atlas_label = "Atlas 3", log_scale = TRUE, n_min = 10,metrics_digits = 3, metrics_label_size = 2)
  
  # Calibration scatterplots
  cal_plot <- plot_grid(cal2,cal3,nrow=1,align="hv")
  #print(cal_plot)
  
  pdf(file.path(fig_dir, paste0(sp_file, "_xval.pdf")), height = 8, width = 12)
  print(maps2$xval_plot)
  print(maps3$xval_plot)
  print(cal_plot)
  dev.off()
  
}

message("09_make_species_figures.R complete: ", fig_dir)
