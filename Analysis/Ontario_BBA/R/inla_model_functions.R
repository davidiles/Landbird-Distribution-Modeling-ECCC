# -----------------------------------
# functions for inla model workflow
# -----------------------------------

fit_inla_testing <- function(sp_dat,
                     study_boundary,
                     covariates,
                     timeout_min,
                     prior_range_abund = c(150,0.1),  
                     prior_sigma_abund = c(0.5,0.1),    
                     prior_range_change = c(500,0.1), 
                     prior_sigma_change = c(0.1,0.1)
) {
  
  # Timeout INLA after 10 minutes (if it has not fit by then, it has likely stalled)
  inla.setOption(inla.timeout = 60*timeout_min)
  
  # Create spatial mesh
  hull <- fm_extensions(
    study_boundary,
    convex = c(50, 150),
    concave = c(50, 150)
  )
  
  mesh_spatial <- fm_mesh_2d_inla(
    loc = st_as_sfc(sp_dat),
    boundary = hull,
    max.edge = c(100, 200), # triangle size near data
    cutoff = 50,           # minimum node separation
    crs = st_crs(sp_dat)
  )
  
  #dim(mesh_spatial$loc)
  #plot(mesh_spatial)
  
  
  # Controls residual spatial field for abundance
  matern_abund <- inla.spde2.pcmatern(mesh_spatial,
                                      prior.range = prior_range_abund,
                                      prior.sigma = prior_sigma_abund,
                                      constr = TRUE
  )
  
  # Controls residual spatial field for change over time
  matern_change <- inla.spde2.pcmatern(mesh_spatial,
                                       prior.range = prior_range_change,
                                       prior.sigma = prior_sigma_change,
                                       constr = TRUE
  )
  
  # Time-of-day effects
  sp_dat$Hours_Since_Sunrise <- as.numeric(sp_dat$Hours_Since_Sunrise)
  HSS_range <- range(sp_dat$Hours_Since_Sunrise)
  HSS_meshpoints <- seq(HSS_range[1] - 1, HSS_range[2] + 1, length.out = 21)
  HSS_mesh1D <- INLA::inla.mesh.1d(HSS_meshpoints, boundary = "free")
  HSS_spde <- INLA::inla.spde2.pcmatern(HSS_mesh1D,
                                        prior.range = c(2, 0.1),
                                        prior.sigma = c(2, 0.5),
                                        constr = TRUE
  )
  
  # Day-of-year effects
  DOY_range <- range(sp_dat$days_since_june15, na.rm = TRUE)
  DOY_meshpoints <- seq(DOY_range[1] - 5,
                        DOY_range[2] + 5,
                        length.out = 21)
  
  DOY_mesh1D <- INLA::inla.mesh.1d(
    DOY_meshpoints,
    boundary = "free"
  )
  
  DOY_spde_global <- INLA::inla.spde2.pcmatern(
    DOY_mesh1D,
    prior.range = c(30, 0.1),   # 10% chance range < 30 days
    prior.sigma = c(2, 0.5),
    constr = TRUE
  )
  
  DOY_spde_BCR <- INLA::inla.spde2.pcmatern(
    DOY_mesh1D,
    prior.range = c(30, 0.1),
    prior.sigma = c(0.3, 0.05),  # shrinkage towards global mean
    constr = TRUE
  )
  
  # iid random effect for atlas squares
  pc_prec <- list(prior = "pcprec", param = c(1, 0.1)) #10% chance SD is larger than 1
  
  
  # Model formulas (if using normal priors)
  covariates <- covariates %>%
    mutate(
      components = paste0(
        "Beta", beta, "_", covariate, '(1,model=\"', model,
        '\"', ", mean.linear = ", mean, ", prec.linear = ",
        prec, ")"
      ),
      formula = paste0("Beta", beta, "_", covariate, "*I(", covariate, "^", beta,")")
    )
  
  
  # DOY_BCR(main = days_since_june15, model = DOY_spde_BCR,group = BCR_factor, control.group = list(model = "exchangeable")) +
  # DOY_BCR +
  
  model_components <- as.formula(paste0(
    '~Intercept(1) +
   effect_Atlas3(1, model="linear", mean.linear = 0, prec.linear = 1) +
   effect_ARU(1, model="linear", mean.linear = 0, prec.linear = 100) +
   
   spde_abund(main = geometry, model = matern_abund) +
   spde_change(main = geometry, model = matern_change) +
   
   HSS(main = Hours_Since_Sunrise, model = HSS_spde) +
   
   DOY_global(main = days_since_june15, model = DOY_spde_global) +
    
   kappa(square_atlas, model = "iid", constr = TRUE, hyper = list(prec = pc_prec))',
    if (length(covariates$components) > 0)
      paste0(" + ", paste(covariates$components, collapse = " + "))
  ))
  
  covar_terms <- paste(covariates$formula, collapse = " + ")
  
  
  # 
  model_formula <- as.formula(paste0(
    "count ~
     Intercept +
     ARU * effect_ARU +
     HSS +
     DOY_global +
     kappa +
     spde_abund +
     Atlas3 * spde_change +
     Atlas3 * effect_Atlas3",
    if (nchar(covar_terms) > 0) paste0(" + ", covar_terms)
  ))
  
  # Specify penalize complexity prior for negative binomial overdispersion
  # (50% chance the size parameter is larger than 5, favouring little overdispersion)
  lambda <- calibrate_pc_lambda(target_prob = 0.5,threshold_theta = 5)
  
  # Fit model
  start <- Sys.time()
  fit_INLA <- NULL
  while (is.null(fit_INLA)) {
    
    fit_INLA <- inlabru::bru(
      
      components = model_components,
      
      # # For point counts
      # inlabru::like(
      #   family = "nbinomial",
      #   formula = model_formula,
      #   data = subset(sp_dat, ARU == 0),
      #   control.family = list(
      #     hyper = list(theta = list(prior = "pc.gamma",
      #                               param = c(lambda)))
      #   )),
      # 
      # # For ARUs
      # inlabru::like(
      #   family = "nbinomial",
      #   formula = model_formula,
      #   data = subset(sp_dat, ARU == 1),
      #   control.family = list(
      #     hyper = list(theta = list(prior = "pc.gamma",
      #                               param = c(lambda)))
      #   )),
      
      # For point counts
      inlabru::like(
        family = "poisson",
        formula = model_formula,
        data = subset(sp_dat, ARU == 0)),
      
      # For ARUs
      inlabru::like(
        family = "poisson",
        formula = model_formula,
        data = subset(sp_dat, ARU == 1)),
      
      options = list(
        inla.mode = "experimental",
        control.compute = list(waic = FALSE, cpo = FALSE),
        control.inla = list(
          int.strategy = "eb",
          strategy = "simplified.laplace"
        ),
        bru_verbose = 4)
      
    )
    
    if ("try-error" %in% class(fit_INLA)) fit_INLA <- NULL
  }
  
  end <- Sys.time()
  runtime_INLA <- difftime(end, start, units = "mins") %>% round(2)
  message(paste0(sp_code, " - ", runtime_INLA, " min to fit model"))
  
  return(fit_INLA)
}

#' Calibrate lambda for a PC prior on the negative binomial size parameter
#'
#' This function finds the \eqn{\lambda} parameter for a penalized complexity (PC) prior
#' such that the prior probability of the negative binomial size parameter
#' (\eqn{\theta}) exceeding a user-defined threshold matches a target probability.
#'
#' Internally, the function uses Monte Carlo simulation and root finding to match
#' the desired tail probability constraint.
#'
#' @param target_prob Numeric. The desired tail probability, i.e. the probability that
#'   \eqn{\theta > threshold}. For example, `0.9` means 90% prior probability that
#'   \eqn{\theta} is greater than the threshold.
#' @param threshold_theta Numeric. The threshold value on the natural \eqn{\theta} scale.
#'   For example, `5` means the constraint applies to \eqn{P(\theta > 5)}.
#' @param n_mc Integer. Number of Monte Carlo draws used to approximate the prior
#'   distribution. Default is `20000`.
#' @param lower Numeric. Lower bound for the root-finding interval on \eqn{\lambda}.
#'   Default is `1e-6`.
#' @param upper Numeric. Upper bound for the root-finding interval on \eqn{\lambda}.
#'   Default is `20`.
#' @param tol Numeric. Tolerance for the root-finding procedure. Default is `1e-3`.
#'
#' @return A single numeric value: the calibrated \eqn{\lambda} parameter for the PC prior.
#'
#' @details
#' The function works by:
#' 1. Sampling from the PC prior distribution using `inla.pc.rgamma()`.
#' 2. Computing the Monte Carlo estimate of \eqn{P(\theta > threshold)}.
#' 3. Using `uniroot()` to find the value of \eqn{\lambda} that matches the user-specified
#'    target probability.
#'
#' Note that depending on the INLA version, `inla.pc.rgamma()` may return values on the
#' log scale. If this is the case, uncomment the transformation `samp <- exp(samp)`
#' inside the helper function.
#'
#' @examples
#' \dontrun{
#' # Example: find lambda so that P(theta > 5) = 0.9
#' lambda <- calibrate_pc_lambda(target_prob = 0.9, threshold_theta = 5)
#' lambda
#' }
#'
#' @export
calibrate_pc_lambda <- function(target_prob,
                                threshold_theta,
                                n_mc = 20000,
                                lower = 1e-6,
                                upper = 20,
                                tol = 1e-3) {
  # Helper: given lambda, estimate tail probability by Monte Carlo
  pc_tail_prob <- function(lambda, n = n_mc, threshold = threshold_theta) {
    samp <- inla.pc.rgamma(n = n, lambda = lambda)
    mean(samp > threshold)
  }

  # Root function: difference between estimated prob and target
  f_for_root <- function(lambda) pc_tail_prob(lambda) - target_prob

  # Root-finding for lambda
  sol <- try(
    uniroot(f_for_root, lower = lower, upper = upper, tol = tol),
    silent = TRUE
  )
  if (inherits(sol, "try-error")) {
    stop("Could not find lambda in the specified interval. Try increasing `upper` or `n_mc`.")
  }
  sol$root
}


summarize_posterior <- function(mat, CI_probs = c(0.05, 0.95), prefix = "var") {
  stopifnot(is.matrix(mat))

  mean_vals   <- matrixStats::rowMeans2(mat, na.rm = TRUE)
  median_vals <- matrixStats::rowMedians(mat, na.rm = TRUE)
  sd_vals     <- matrixStats::rowSds(mat, na.rm = TRUE)
  cv_vals     <- sd_vals / median_vals

  lower_vals <- matrixStats::rowQuantiles(mat, probs = CI_probs[1], na.rm = TRUE)
  upper_vals <- matrixStats::rowQuantiles(mat, probs = CI_probs[2], na.rm = TRUE)

  out <- data.frame(
    setNames(list(mean_vals),   paste0(prefix, "_mean")),
    setNames(list(median_vals), paste0(prefix, "_q50")),
    setNames(list(sd_vals),     paste0(prefix, "_sd")),
    setNames(list(cv_vals),     paste0(prefix, "_cv_median")),
    setNames(list(lower_vals),  paste0(prefix, "_lower")),
    setNames(list(upper_vals),  paste0(prefix, "_upper"))
  )

  return(out)
}


predict_inla <- function(mod, grid, pred_formula) {

  start <- Sys.time()

  pred <- inlabru::generate(mod,
                            grid,
                            formula = pred_formula,
                            n.samples = 1000,
                            seed = 123
  )

  error_type <- mod$bru_info$lhoods[[1]]$family

  # Reformat pred such to a named list, with prediction matrices for each object (n_grid_cells x n_samples)
  pred_vars <- names(pred[[1]])
  pred <- lapply(pred_vars, function(v) sapply(pred, function(x) x[[v]]))
  names(pred) <- pred_vars

  end <- Sys.time()

  runtime_pred <- difftime(end, start, units = "mins") %>% round(2)
  message(paste0(sp_code, " - ", runtime_pred, " min to generate predictions"))
  return(pred)
}


map_relabund <- function(species_name,
                         species_filename,
                         obs_dat,
                         grid,
                         preds_summarized,
                         atlas_squares,
                         study_boundary,
                         BCR_shapefile = allBCR,
                         water_shapefile = water,
                         map_dir = "figures/species_maps/",
                         train_dat_filter = "TRUE",
                         prefix = "OBBA3",
                         plot_obs_data = TRUE,
                         title = "Relative Abundance",
                         subtitle = "Per 5-minute point count",
                         upper_bound = 1,
                         lower_bound = 0.01,
                         res = 1.1) {

  proj_use <- st_crs(obs_dat)

  # Helper to split species name if it's too long
  wrap_species_label <- function(label, max_length = 15) {
    if (nchar(label) <= max_length) return(label)

    words <- strsplit(label, " ")[[1]]
    if (length(words) == 1) return(label)  # Single word, don't split

    # Put everything except the last word on the first line
    paste0(paste(words[-length(words)], collapse = " "), "<br>", words[length(words)])
  }

  # Summarize atlas_squares where species was detected
  sp_detected <- obs_dat %>%
    sf::st_intersection(atlas_squares %>% st_transform(st_crs(obs_dat))) %>%
    as.data.frame() %>%
    group_by(square_id_) %>%
    summarize(
      sp_detected = as.numeric(sum(count) > 0),
      sp_mean_count = mean(count) %>% round(2)
    )

  atlas_squares_species <- atlas_squares %>%
    relocate(geometry, .after = last_col()) %>%
    left_join(sp_detected, by = join_by(square_id_)) # %>% left_join(CL_detected)

  atlas_squares_centroids <- sf::st_centroid(atlas_squares_species)

  # ---- Bind summarized predictions to grid
  q50_col <- names(preds_summarized)[endsWith(names(preds_summarized), "_q50")]
  grid$pred_q50 <- preds_summarized[[q50_col]]

  CV_col <- names(preds_summarized)[endsWith(names(preds_summarized), "_cv_median")]
  grid$pred_CV <- preds_summarized[[CV_col]]

  # ---- Create bounding box for study area
  bb <- st_bbox(study_boundary)

  # Expand by 10%
  xrange <- bb$xmax - bb$xmin
  yrange <- bb$ymax - bb$ymin

  xbuf <- xrange * 0.1
  ybuf <- yrange * 0.1

  bb_expanded <- bb
  bb_expanded$xmin <- bb$xmin - xbuf
  bb_expanded$xmax <- bb$xmax + xbuf
  bb_expanded$ymin <- bb$ymin - ybuf
  bb_expanded$ymax <- bb$ymax + ybuf

  # ---- Plot median predictions

  # Bounds for plotting and associated labels
  breaks <- 10^seq(log10(lower_bound),log10(upper_bound),length.out = 5) %>% signif(2)

  # Cap values at upper/lower bounds
  grid$pred_capped <- as.numeric(pmax(pmin(grid$pred_q50, upper_bound), lower_bound))

  # Convert sf to SpatVector
  v <- terra::vect(grid)

  # Create a raster template with desired resolution
  r_template <- terra::rast(v, res = res)

  # Rasterize pred_capped values using mean within each cell
  pred_rast <- terra::rasterize(v, r_template, field = "pred_capped", fun = mean)

  # Convert raster to stars for plotting
  pred_rast_stars <- stars::st_as_stars(pred_rast)

  # Set legend labels
  break_labels <- as.character(breaks)
  break_labels[1] <- paste0("<", break_labels[1])
  break_labels[length(break_labels)] <- paste0(">", break_labels[length(break_labels)])

  colscale_q50 <- c(
    "#FEFEFE", "#FBF7E2", "#FCF8D0", "#EEF7C2", "#CEF2B0",
    "#94E5A0", "#51C987", "#18A065", "#008C59", "#007F53", "#006344"
  )
  colpal_q50 <- colorRampPalette(colscale_q50)

  q50_plot <- ggplot() +
    geom_sf(data = BCR_shapefile, colour = "gray80", fill = NA, lwd = 0.3) +
    stars::geom_stars(data = pred_rast_stars) +
    
    scale_fill_gradientn(
      name = paste0(
        "<span style='font-size:20pt; font-weight:bold'>", wrap_species_label(species_name), "</span><br><br>",
        "<span style='font-size:14pt'>", title, "</span><br>",
        "<span style='font-size:7pt'>", subtitle, "</span><br>",
        "<span style='font-size:7pt'>Posterior Median</span>"
      ),
      colors = colpal_q50(10),
      na.value = "transparent",
      
      limits = c(lower_bound,round(upper_bound*1.05,3))
      
      # # --- If using log scale, uncomment below
      #trans = "log10",
      #breaks = breaks,
      #labels = break_labels,
      #limits = c(min(breaks) / 1.1, max(breaks) * 1.1)
      
    ) +

    geom_sf(data = water_shapefile, fill = "#EDF7FB", col = "transparent")+
    geom_sf(data = atlas_squares_centroids %>% subset(!is.na(sp_mean_count)), colour = "black", size = 0.75, shape = 1, stroke = 0.1, alpha = 0.5) +
    geom_sf(data = atlas_squares_centroids %>% subset(sp_mean_count>0), colour = "black", size = 0.75, shape = 4, stroke = 0.2, alpha = 0.5) +
    geom_sf(data = study_boundary, colour = "black", fill = NA, lwd = 0.5, show.legend = FALSE) +

    coord_sf(
      xlim = c(bb_expanded$xmin, bb_expanded$xmax),
      ylim = c(bb_expanded$ymin, bb_expanded$ymax),
      expand = FALSE
    ) +

    theme_void() +
    theme(
      panel.background = element_rect(fill = "#F5F5F5", color = NA),
      plot.margin = unit(c(0, 0, 0, 0), "pt"),
      legend.title = ggtext::element_markdown(lineheight = .9),
      legend.position = c(0.02,0.02),
      legend.justification = c(0,0),
      legend.background = element_rect(fill = "white", color = "black"),
      legend.margin = margin(5, 5, 5, 5),
      panel.border = element_rect(color = "black", fill = NA, size = 0.8)
    ) +
    ggspatial::annotation_scale(
      location = "br",      # bottom right
      width_hint = 0.3      # relative width of scale bar
    ) +
    ggspatial::annotation_north_arrow(
      location = "tr",      # top right
      which_north = "true",
      pad_x = unit(0.2, "in"),
      pad_y = unit(0.2, "in"),
      style = ggspatial::north_arrow_fancy_orienteering()
    )

  ggsave(
    filename = paste0(map_dir, "/", species_filename, "_relabund_q50_", prefix, ".png"),
    plot = q50_plot,
    width = 10,
    height = 8,
    units = "in",
    dpi = 1000,
    type = "cairo",
    limitsize = FALSE
  )


  # ---- Plot uncertainty in predictions (width of 90% CRI)

  colscale_uncertainty <- c("#FEFEFE", "#FFF4B3", "#F5D271", "#F2B647", "#EC8E00", "#CA302A")
  colpal_uncertainty <- colorRampPalette(colscale_uncertainty)

  # Cap values at upper/lower bounds
  grid$pred_capped <- as.numeric(pmax(pmin(grid$pred_CV, 1), 0))

  # Convert sf to SpatVector
  v <- terra::vect(grid)

  # Create a raster template with desired resolution
  r_template <- terra::rast(v, res = res)

  # Rasterize pred_capped values using mean within each cell
  pred_rast <- terra::rasterize(v, r_template, field = "pred_capped", fun = mean)

  # Convert raster to stars for plotting
  pred_rast_stars <- stars::st_as_stars(pred_rast)

  # Set legend labels/breaks
  breaks <- seq(0,1,length.out = 5) %>% signif(2)
  break_labels <- as.character(breaks)
  break_labels[length(break_labels)] <- paste0(">", break_labels[length(break_labels)])

  CV_plot <- ggplot() +
    geom_sf(data = BCR_shapefile, colour = "gray80", fill = NA, lwd = 0.3) +
    stars::geom_stars(data = pred_rast_stars) +
    scale_fill_gradientn(
      name = paste0(
        "<span style='font-size:20pt; font-weight:bold'>", wrap_species_label(species_name), "</span><br><br>",
        "<span style='font-size:14pt'>", title, "</span><br>",
        "<span style='font-size:7pt'>", subtitle, "</span><br>",
        "<span style='font-size:7pt'>Prediction CV</span>"
      ),
      colors = colpal_uncertainty(11),
      na.value = "transparent",
      breaks = breaks,
      labels = break_labels,
      limits = c(min(breaks) / 1.1, max(breaks) * 1.1)
    ) +
    geom_sf(data = water_shapefile, fill = "#EDF7FB", col = "transparent")+
    geom_sf(data = atlas_squares_centroids %>% subset(!is.na(sp_mean_count)), colour = "black", size = 0.75, shape = 1, stroke = 0.1, alpha = 0.5) +
    geom_sf(data = atlas_squares_centroids %>% subset(sp_mean_count>0), colour = "black", size = 0.75, shape = 4, stroke = 0.2, alpha = 0.5) +
    geom_sf(data = study_boundary, colour = "black", fill = NA, lwd = 0.5, show.legend = FALSE) +

    coord_sf(
      xlim = c(bb_expanded$xmin, bb_expanded$xmax),
      ylim = c(bb_expanded$ymin, bb_expanded$ymax),
      expand = FALSE
    ) +

    theme_void() +
    theme(
      panel.background = element_rect(fill = "#F5F5F5", color = NA),
      plot.margin = unit(c(0, 0, 0, 0), "cm"),
      legend.title = ggtext::element_markdown(lineheight = .9),
      legend.position = c(0.02,0.02),
      legend.justification = c(0,0),
      legend.background = element_rect(fill = "white", color = "black"),
      legend.margin = margin(5, 5, 5, 5),
      panel.border = element_rect(color = "black", fill = NA, size = 0.8)
    ) +
    ggspatial::annotation_scale(
      location = "br",      # bottom right
      width_hint = 0.3      # relative width of scale bar
    ) +
    ggspatial::annotation_north_arrow(
      location = "tr",      # top right
      which_north = "true",
      pad_x = unit(0.2, "in"),
      pad_y = unit(0.2, "in"),
      style = ggspatial::north_arrow_fancy_orienteering()
    )

  #png(paste0(map_dir,"/",species_name,"_relabund_CV_",prefix,".png"), width = 10, height = 8, units = "in", res = 1000, type = "cairo")
  #print(CV_plot)
  #dev.off()

  ggsave(
    filename = paste0(map_dir,"/",species_filename,"_relabund_CV_",prefix,".png"),
    plot = CV_plot,
    width = 10,
    height = 8,
    units = "in",
    dpi = 1000,
    type = "cairo",
    limitsize = FALSE
  )

}

map_change <- function(species_name,
                       species_filename,
                       grid,
                       preds_summarized,
                       study_boundary,
                       BCR_shapefile = allBCR,
                       water_shapefile = water,
                       map_dir = "figures/species_maps/",
                       upper_bound = -1,
                       lower_bound = 1,
                       res = 1.1,
                       change_type = "Percent") {

  proj_use <- st_crs(grid)

  # Helper to split species name if it's too long
  wrap_species_label <- function(label, max_length = 15) {
    if (nchar(label) <= max_length) return(label)

    words <- strsplit(label, " ")[[1]]
    if (length(words) == 1) return(label)  # Single word, don't split

    # Put everything except the last word on the first line
    paste0(paste(words[-length(words)], collapse = " "), "<br>", words[length(words)])
  }

  # ---- Bind summarized predictions to grid
  q50_col <- names(preds_summarized)[endsWith(names(preds_summarized), "_q50")]
  grid$pred_q50 <- preds_summarized[[q50_col]]

  CV_col <- names(preds_summarized)[endsWith(names(preds_summarized), "_cv_median")]
  grid$pred_CV <- preds_summarized[[CV_col]]

  lcl_col <- names(preds_summarized)[endsWith(names(preds_summarized), "_lower")]
  grid$pred_lcl <- preds_summarized[[lcl_col]]

  ucl_col <- names(preds_summarized)[endsWith(names(preds_summarized), "_upper")]
  grid$pred_ucl <- preds_summarized[[ucl_col]]

  grid$pred_ciw <- grid$pred_ucl - grid$pred_lcl

  # ---- Create bounding box for study area
  bb <- st_bbox(study_boundary)

  # Expand by 10%
  xrange <- bb$xmax - bb$xmin
  yrange <- bb$ymax - bb$ymin

  xbuf <- xrange * 0.1
  ybuf <- yrange * 0.1

  bb_expanded <- bb
  bb_expanded$xmin <- bb$xmin - xbuf
  bb_expanded$xmax <- bb$xmax + xbuf
  bb_expanded$ymin <- bb$ymin - ybuf
  bb_expanded$ymax <- bb$ymax + ybuf

  # ---- Plot median predictions

  grid$pred_capped <- as.numeric(pmax(pmin(grid$pred_q50, upper_bound), lower_bound))

  # Convert sf to SpatVector
  v <- terra::vect(grid)

  # Create a raster template with desired resolution
  r_template <- terra::rast(v, res = res)

  # Rasterize pred_capped values using mean within each cell
  pred_rast <- terra::rasterize(v, r_template, field = "pred_capped", fun = mean)

  # Convert raster to stars for plotting
  pred_rast_stars <- stars::st_as_stars(pred_rast)

  # Bounds for plotting and associated labels
  breaks <- seq(lower_bound,upper_bound,length.out = 7)
  breaks[4] <- 0

  # Set legend labels
  if (change_type == "Percent"){
    break_labels <- (100 * (exp(breaks) - 1)) %>% signif(2)
    break_labels <- paste0(break_labels,"%")
    prefix = "pct_change"
    title = "Percent change"
  } else{
    break_labels = breaks %>% signif(2)
    prefix = "abs_change"
    title = "Absolute change"
  }
  break_labels[breaks>0] <- paste0("+",break_labels[breaks>0])
  break_labels[1] <- paste0("< ", break_labels[1])
  break_labels[length(break_labels)] <- paste0("> ", break_labels[length(break_labels)])

  colscale_q50 <- RColorBrewer::brewer.pal(11,"RdBu")
  colpal_q50 <- colorRampPalette(colscale_q50)

  q50_plot <- ggplot() +
    geom_sf(data = BCR_shapefile, colour = "gray80", fill = NA, lwd = 0.3) +
    stars::geom_stars(data = pred_rast_stars) +
    scale_fill_gradientn(
      name = paste0(
        "<span style='font-size:20pt; font-weight:bold'>", wrap_species_label(species_name), "</span><br><br>",
        "<span style='font-size:14pt'>", title, "</span><br>",
        "<span style='font-size:7pt'>OBBA2 to OBBA3</span><br>",
        "<span style='font-size:7pt'>Posterior Median</span>"
      ),
      colors = colpal_q50(11),
      na.value = "transparent",
      breaks = breaks,
      labels = break_labels,
      limits = c(min(breaks) * 1.1, max(breaks) * 1.1)
    ) +

    geom_sf(data = water_shapefile, fill = "#F5F5F5", col = "transparent")+
    geom_sf(data = study_boundary, colour = "black", fill = NA, lwd = 0.5, show.legend = FALSE) +

    coord_sf(
      xlim = c(bb_expanded$xmin, bb_expanded$xmax),
      ylim = c(bb_expanded$ymin, bb_expanded$ymax),
      expand = FALSE
    ) +

    theme_void() +
    theme(
      panel.background = element_rect(fill = "#F5F5F5", color = NA),
      plot.margin = unit(c(0, 0, 0, 0), "cm"),
      legend.title = ggtext::element_markdown(lineheight = .9),
      legend.position = c(0.02,0.02),
      legend.justification = c(0,0),
      legend.background = element_rect(fill = "white", color = "black"),
      legend.margin = margin(5, 5, 5, 5),
      panel.border = element_rect(color = "black", fill = NA, size = 0.8)
    ) +
    ggspatial::annotation_scale(
      location = "br",      # bottom right
      width_hint = 0.3      # relative width of scale bar
    ) +
    ggspatial::annotation_north_arrow(
      location = "tr",      # top right
      which_north = "true",
      pad_x = unit(0.2, "in"),
      pad_y = unit(0.2, "in"),
      style = ggspatial::north_arrow_fancy_orienteering()
    )

  ggsave(
    filename = paste0(map_dir,"/",species_filename,"_",prefix,"_q50.png"),
    plot = q50_plot,
    width = 10,
    height = 8,
    units = "in",
    dpi = 1000,
    type = "cairo",
    limitsize = FALSE
  )

  # ---- Plot uncertainty in predictions (width of 90% CRI)

  colscale_uncertainty <- c("#FEFEFE", "#FFF4B3", "#F5D271", "#F2B647", "#EC8E00", "#CA302A")
  colpal_uncertainty <- colorRampPalette(colscale_uncertainty)

  # Cap values at upper/lower bounds
  grid$ci_width <- grid$pred
  grid$pred_capped <- as.numeric(pmax(pmin(grid$pred_ciw, upper_bound), 0))

  # Convert sf to SpatVector
  v <- terra::vect(grid)

  # Create a raster template with desired resolution
  r_template <- terra::rast(v, res = res)

  # Rasterize pred_capped values using mean within each cell
  pred_rast <- terra::rasterize(v, r_template, field = "pred_capped", fun = mean)

  # Convert raster to stars for plotting
  pred_rast_stars <- stars::st_as_stars(pred_rast)

  # Bounds for plotting and associated labels
  breaks <- seq(0,upper_bound,length.out = 7)

  # Set legend labels
  if (change_type == "Percent"){
    break_labels <- (100 * (exp(breaks) - 1)) %>% signif(2)
    break_labels <- paste0(break_labels,"%")
    prefix = "pct_change"
    title = "Percent change"
  } else{
    break_labels = breaks %>% signif(2)
    prefix = "abs_change"
    title = "Absolute Change"
  }
  break_labels[breaks>0] <- paste0("+",break_labels[breaks>0])
  break_labels[length(break_labels)] <- paste0("> ", break_labels[length(break_labels)])

  ciw_plot <- ggplot() +
    geom_sf(data = BCR_shapefile, colour = "gray80", fill = NA, lwd = 0.3) +
    stars::geom_stars(data = pred_rast_stars) +
    scale_fill_gradientn(
      name = paste0(
        "<span style='font-size:20pt; font-weight:bold'>", wrap_species_label(species_name), "</span><br><br>",
        "<span style='font-size:14pt'>", title, "</span><br>",
        "<span style='font-size:7pt'>OBBA2 to OBBA3</span><br>",
        "<span style='font-size:7pt'>90% CI width</span>"
      ),
      colors = colpal_uncertainty(11),
      na.value = "transparent",
      breaks = breaks,
      labels = break_labels,
      limits = c(min(breaks) * 1.1, max(breaks) * 1.1)
    ) +

    geom_sf(data = water_shapefile, fill = "#F5F5F5", col = "transparent")+
    geom_sf(data = study_boundary, colour = "black", fill = NA, lwd = 0.5, show.legend = FALSE) +

    coord_sf(
      xlim = c(bb_expanded$xmin, bb_expanded$xmax),
      ylim = c(bb_expanded$ymin, bb_expanded$ymax),
      expand = FALSE
    ) +

    theme_void() +
    theme(
      panel.background = element_rect(fill = "#F5F5F5", color = NA),
      plot.margin = unit(c(0, 0, 0, 0), "cm"),
      legend.title = ggtext::element_markdown(lineheight = .9),
      legend.position = c(0.02,0.02),
      legend.justification = c(0,0),
      legend.background = element_rect(fill = "white", color = "black"),
      legend.margin = margin(5, 5, 5, 5),
      panel.border = element_rect(color = "black", fill = NA, size = 0.8)
    ) +
    ggspatial::annotation_scale(
      location = "br",      # bottom right
      width_hint = 0.3      # relative width of scale bar
    ) +
    ggspatial::annotation_north_arrow(
      location = "tr",      # top right
      which_north = "true",
      pad_x = unit(0.2, "in"),
      pad_y = unit(0.2, "in"),
      style = ggspatial::north_arrow_fancy_orienteering()
    )

  ggsave(
    filename = paste0(map_dir,"/",species_filename,"_",prefix,"_ciw.png"),
    plot = ciw_plot,
    width = 10,
    height = 8,
    units = "in",
    dpi = 1000,
    type = "cairo",
    limitsize = FALSE
  )


}
