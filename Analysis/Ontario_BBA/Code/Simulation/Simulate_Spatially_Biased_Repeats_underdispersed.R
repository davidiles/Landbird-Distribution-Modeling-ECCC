
# ------------------------------------
# In these simulations, assume each bird on the landscape does not move (will be repeated detected across multiple surveys)
# Are estimates of mean count biased?
# Are estimates of differences in mean count (between east and west) biased, given that eastern sites are surveyed 10 times while western sites are only surveyed once
# ------------------------------------


# Results: if there is extreme underdispersion at survey locations (e.g., binomial detections with p = 0.2) the model overestimates abundance by about 10%
         # if there is moderate underdispersion (p = 0.5), the model is nearly unbiased
         # spatial unbalance in repeat surveys does not affect results

library(sf)
library(tidyverse)
require(INLA)     # INLA_22.05.07  
require(inlabru)  # inlabru_2.7.0
require(fmesher)
library(viridis)

rm(list=ls())

# ------------------------------------
# Relationship between mean and variance for binomial 
# ------------------------------------
mean <- 1
p <- seq(0,1,length.out = 100)
n <- mean/p
var <- n*p*(1-p)

plot(var/mean~p)

# ------------------------------------
# Simulations
# ------------------------------------

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

n_simulations <- 100

# Measuring the probability that east part of region has higher abundance than west
results <- data.frame(sim_number = 1:n_simulations,mean_count = NA,prob_difference = NA)

for (sim_rep in 1:n_simulations){
  
  # Select sample locations
  n = 500
  sample_locs <- data.frame(site = 1:n,x = runif(n,-1000,1000),y = runif(n,-1000,1000)) 
  
  # Simulate sampling
  dat <- expand.grid(site = 1:n, survey = 1:10)
  dat$x <- sample_locs$x[dat$site]
  dat$y<- sample_locs$y[dat$site]
  
  # Sites with only a single survey (rather than 10)
  single_surveys <- subset(dat, x<0)$site %>% unique() # If all repeat surveys are in east
  #single_surveys <- unique(dat$site) %>% sample(.,round(length(.)/2)) # If repeat surveys are random
  
  dat <- dat %>% subset((site %in% single_surveys & survey == 1) | !(site %in% single_surveys)) %>%  st_as_sf(coords = c("x", "y"), crs = 32610, remove = FALSE)
  
  dat <- dat %>% arrange(site)
  
  # Simulate counts that vary among repeated visits, but are potentially underdispersed (as p gets large, counts are more underdispersed)
  p <- 0.9 # Controls underdispersion for methods 1 and 2. Higher values of p result in more underdispersion (more similar repeated counts at locations)
  
  #~~~ Method #1; pure binomial
  #mean <- 1
  #N <- round(mean/p)
  #dat$count <- rbinom(nrow(dat),N,p)
  
  #~~~ Method #2; N-mixture model
   N <- rpois(n,mean/p) # True super-population abundance at each site
   dat$count <- rbinom(nrow(dat),N[dat$site],p)
  
  #~~~ Method #3; Pure poisson
  #dat$count <- rpois(nrow(dat))

  data_summary <- dat %>% group_by(site) %>% summarize(n_surveys = n(), mean_count = mean(count))
  
  # Fix counts to be identical to the first visit, even on repeat visits
  ggplot(data =  data_summary) + geom_sf(aes(col = n_surveys))
  ggplot(data =  data_summary) + geom_sf(aes(col = mean_count))
  
  # ------------------------------------
  # Fit a spatial model to the simulated data
  # ------------------------------------
  
  mesh_locs <- mesh$loc[,c(1,2)] %>% as.data.frame()

  # ----------------------------------------------------
  # iid random effect for survey location
  # ----------------------------------------------------
  
  pc_prec <- list(prior = "pcprec", param = c(0.1, 0.1))
  dat$site_id <- as.numeric(as.factor(dat$site))
  
  ggplot()+geom_sf(data = study_boundary)+geom_sf(data = dat)
  
  # Controls the 'residual spatial field'.  This can be adjusted to create smoother surfaces.
  matern <- inla.spde2.pcmatern(mesh,
                                prior.range = c(500,0.01), # 500, NA # 1% chance range is smaller than 500000
                                prior.sigma = c(0.1,0.1),   # 10% chance sd is larger than 0.1,
                                constr = TRUE
  )
  
  comp <- ~ spde(main = geometry, model = matern) + kappa(site_id, model = "iid", constr = TRUE, hyper = list(prec = pc_prec))
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
  
  # Add uncertainty from kappa
  
  kappa_var <- 1/fit$summary.hyperpar$mean[3]
  
  pred_df <- expand.grid(x = seq(-1000,1000,length.out = 100),
                         y = seq(-1000,1000,length.out = 100)) %>%
    st_as_sf(coords = c("x", "y"), crs = st_crs(dat), remove = FALSE)
  
  pred <- (generate(fit,pred_df,formula = ~ Intercept + spde, n.samples = 100) + 0.5*kappa_var) %>% exp()
  #pred <- predict(fit,pred_df, n.samples = 100) 
  
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
  results$mean_count[sim_rep] <- apply(pred,2,mean) %>% quantile(0.5)
  
  hist(results$prob_difference, breaks = seq(0,1,length.out = 20))
  hist(results$mean_count, breaks = seq(0.5,1.5,length.out = 20))
}

mean(results$mean_count,na.rm = TRUE)
mean(results$prob_difference,na.rm = TRUE)

# No bias