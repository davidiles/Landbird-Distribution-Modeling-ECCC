library(sf)
library(tidyverse)
require(INLA)     # INLA_22.05.07  
require(inlabru)  # inlabru_2.7.0
require(fmesher)
library(viridis)

rm(list=ls())


study_boundary <- data.frame(x = c(-1000,-1000,1000,1000),
                             y = c(-1000,1000,-1000,1000)) %>%  
  st_as_sf(coords = c("x", "y"), crs = 32610) %>% 
  st_bbox() %>% st_as_sfc() 

hull <- fm_extensions(study_boundary)

mesh <- fm_mesh_2d_inla(
  boundary = hull, 
  max.edge = c(100, 500), # km inside and outside  
  cutoff = 5,
  crs = 32610
)

# ------------------------------------
# Simulation to confirm that spatially biased repeated sampling does not introduce bias into estimates of density
# ------------------------------------

n_simulations <- 100

# Measuring the probability that east part of region has higher abundance than west
results <- data.frame(sim_number = 1:n_simulations,prob_difference = NA)


for (sim_rep in 1:n_simulations){
  
  # Select sample locations
  n = 500
  sample_locs <- data.frame(site = 1:n,x = runif(n,-1000,1000),y = runif(n,-1000,1000)) 
  
  # Simulate sampling (sites with x > 0 are surveyed 10 times, sites with x<0 are only surveyed once)
  dat <- expand.grid(site = 1:n, survey = 1:10)
  dat$x <- sample_locs$x[dat$site]
  dat$y<- sample_locs$y[dat$site]
  
  dat <- dat %>% subset((x<0 & survey == 1) | (x>0)) %>%  st_as_sf(coords = c("x", "y"), crs = 32610)
  
  # Simulate counts
  dat$count <- rpois(nrow(dat),1)
  
  
  #ggplot(data = dat %>% group_by(site) %>% summarize(n_surveys = n())) + geom_sf(aes(col = n_surveys))
  
  # ------------------------------------
  # Fit a spatial model to the simulated data
  # ------------------------------------
  
  mesh_locs <- mesh$loc[,c(1,2)] %>% as.data.frame()
  #dim(mesh_locs)
  #plot(mesh)
  
  ggplot()+geom_sf(data = study_boundary)+geom_sf(data = dat)
  
  # Controls the 'residual spatial field'.  This can be adjusted to create smoother surfaces.
  matern <- inla.spde2.pcmatern(mesh,
                                prior.range = c(50,0.01), # 500, NA # 1% chance range is smaller than 500000
                                prior.sigma = c(0.1,0.1),   # 10% chance sd is larger than 0.1,
                                constr = TRUE
  )
  
  comp <- ~ spde(main = geometry, model = matern)
  form <- count ~ .
  
  fit <- bru(components = comp,
             like(family = "poisson",
                  formula = form,
                  data = dat),
             
             options = list(inla.mode = "experimental",
                            bru_initial = list(),
                            control.compute = list(waic = FALSE, cpo = FALSE),
                            bru_verbose = 4))
  
  
  # ------------------------------------
  # Predictions
  # ------------------------------------
  
  pred_df <- expand.grid(x = seq(-1000,1000,length.out = 100),
                         y = seq(-1000,1000,length.out = 100)) %>%
    st_as_sf(coords = c("x", "y"), crs = st_crs(dat), remove = FALSE)
  
  pred <- generate(fit,pred_df,formula = ~Intercept + spde, n.samples = 100) %>% exp()
  
  pred_df$pred_q50 <- apply(pred,1,function(x)quantile(x,0.5))
  pred_df$pred_q05 <- apply(pred,1,function(x)quantile(x,0.05))
  pred_df$pred_q95 <- apply(pred,1,function(x)quantile(x,0.95))
  pred_df$ci_width <- pred_df$pred_q95 - pred_df$pred_q05
  
  # Plot mean predictions
  ggplot(data = pred_df, aes(col = pred_q50))+geom_sf(size=2, shape = 15) + scale_color_gradientn(colours = viridis(10))
  
  # Plot uncertainty
  ggplot(data = pred_df, aes(col = ci_width))+geom_sf(size=2, shape = 15) + scale_color_gradientn(colours = viridis(10))
  
  # Estimate probability that abundance in region where x>0 is greater than x<0
  rows <- which(pred_df$x > 0)
  
  diff_estimate <- c()
  
  for (sample in 1:dim(pred)[2]){
    diff_estimate[sample] <- sum(pred[rows,sample]) - sum(pred[-rows,sample])
  }
  
  results$prob_difference[sim_rep] <- mean(diff_estimate>0)
}

hist(results$prob_difference)

# No bias