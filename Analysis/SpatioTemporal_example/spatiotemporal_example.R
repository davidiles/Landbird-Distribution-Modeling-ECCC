# Example adapted from https://inlabru-org.github.io/inlabru/articles/2d_lgcp_spatiotemporal.html

rm(list=ls())

library(inlabru)
library(INLA)
library(ggplot2)

# ---------------------------------------------------------
# Load data
# ---------------------------------------------------------

data(mrsea, package = "inlabru")

ggplot() +
  gg(mrsea$mesh) +
  gg(mrsea$boundary) +
  gg(mrsea$samplers) +
  gg(mrsea$points, size = 0.5) +
  facet_wrap(~season) +
  ggtitle("MRSea observation seasons")

ips <- fm_int(
  domain = list(coordinates = mrsea$mesh, season = 1:4),
  samplers = mrsea$samplers
)

ggplot() +
  gg(mrsea$mesh) +
  gg(ips, aes(size = weight)) +
  scale_size_area(max_size = 1) +
  facet_wrap(~season)

# Specify matern
matern <- inla.spde2.pcmatern(mrsea$mesh,
                              prior.sigma = c(0.1, 0.01),
                              prior.range = c(10, 0.01)
)

# ---------------------------------------------------------
# Fit model separately for each season
# ---------------------------------------------------------

cmp <- coordinates ~ Intercept(1) +
  mySmooth(coordinates,model = matern)

fit1 <- lgcp(cmp,
            data = subset(mrsea$points, season == 1),
            samplers = mrsea$samplers,
            domain = list(coordinates = mrsea$mesh))

fit2 <- lgcp(cmp,
             data = subset(mrsea$points, season == 2),
             samplers = mrsea$samplers,
             domain = list(coordinates = mrsea$mesh))

fit3 <- lgcp(cmp,
             data = subset(mrsea$points, season == 3),
             samplers = mrsea$samplers,
             domain = list(coordinates = mrsea$mesh))

fit4 <- lgcp(cmp,
             data = subset(mrsea$points, season == 4),
             samplers = mrsea$samplers,
             domain = list(coordinates = mrsea$mesh))

# Predictions
ppxl <- fm_pixels(mrsea$mesh, mask = mrsea$boundary, format = "sp")

lambda1 <- predict(
  fit1,
  ppxl,
  ~ data.frame(lambda = exp(mySmooth + Intercept))
)

lambda2 <- predict(
  fit2,
  ppxl,
  ~ data.frame(lambda = exp(mySmooth + Intercept))
)

lambda3 <- predict(
  fit3,
  ppxl,
  ~ data.frame(lambda = exp(mySmooth + Intercept))
)

lambda4 <- predict(
  fit4,
  ppxl,
  ~ data.frame(lambda = exp(mySmooth + Intercept))
)

# ---------------------------------------------------------
# Fit spatio-temporal model
# ---------------------------------------------------------

cmp <- coordinates ~ Intercept(1) +
  smooth1(coordinates,model = matern)+
  smooth2(coordinates,model = matern)+
  smooth3(coordinates,model = matern)+
  smooth4(coordinates,model = matern)
  

start <- Sys.time()
fit_INLA <- NULL
while(is.null(fit_INLA)){
  
  fit_model <- function(){
    tryCatch(expr = {bru(components = cmp,
                         
                         like(family = "cp",
                              formula = as.formula(coordinates ~ Intercept + smooth1),
                              data = subset(mrsea$points, season == 1),
                              samplers = mrsea$samplers,
                              domain = list(coordinates = mrsea$mesh)),
                         
                         like(family = "cp",
                              formula = as.formula(coordinates ~ Intercept + smooth1 + smooth2),
                              data = subset(mrsea$points, season == 2),
                              samplers = mrsea$samplers,
                              domain = list(coordinates = mrsea$mesh)),
                         
                         like(family = "cp",
                              formula = as.formula(coordinates ~ Intercept + smooth1 + smooth2 + smooth3),
                              data = subset(mrsea$points, season == 3),
                              samplers = mrsea$samplers,
                              domain = list(coordinates = mrsea$mesh)),
                         
                         like(family = "cp",
                              formula = as.formula(coordinates ~ Intercept + smooth1 + smooth2 + smooth3 + smooth4),
                              data = subset(mrsea$points, season == 4),
                              samplers = mrsea$samplers,
                              domain = list(coordinates = mrsea$mesh)),
                         
                         options = list(control.compute = list(waic = FALSE, cpo = FALSE),
                                        bru_verbose = 4))},
             error = function(e){NULL})
  }
  fit_INLA <- fit_model()
  
  if ("try-error" %in% class(fit_INLA)) fit_INLA <- NULL
}

end <- Sys.time()
runtime_INLA <- difftime( end,start, units="mins") %>% round(2)

# -------------------------
# Generate predictions
# -------------------------
ppxl <- fm_pixels(mrsea$mesh, mask = mrsea$boundary, format = "sp")

# Season 1
pred1 <- generate(fit_INLA,ppxl,formula = as.formula(coordinates ~ Intercept + smooth1)) %>% exp()

# Season 2
pred2 <- generate(fit_INLA,ppxl,formula = as.formula(coordinates ~ Intercept + smooth1 + smooth2)) %>% exp()

# Season 3
pred3 <- generate(fit_INLA,ppxl,formula = as.formula(coordinates ~ Intercept + smooth1 + smooth2 + smooth3)) %>% exp()

# Season 4
pred4 <- generate(fit_INLA,ppxl,formula = as.formula(coordinates ~ Intercept + smooth1 + smooth2 + smooth3 + smooth4)) %>% exp()

# -------------------------
# Plot predictions
# -------------------------

mean1 <- apply(pred1,1,mean)
mean2 <- apply(pred2,1,mean)
mean3 <- apply(pred3,1,mean)
mean4 <- apply(pred4,1,mean)

ppxl$mean1 <- mean1
ppxl$mean2 <- mean2
ppxl$mean3 <- mean3
ppxl$mean4 <- mean4

pl1 <- ggplot() +
  gg(as(ppxl, "SpatialPixelsDataFrame"), aes(fill = mean1)) +
  coord_equal()
pl1


pl2 <- ggplot() +
  gg(as(ppxl, "SpatialPixelsDataFrame"), aes(fill = mean2)) +
  coord_equal()
pl2

pl3 <- ggplot() +
  gg(as(ppxl, "SpatialPixelsDataFrame"), aes(fill = mean3)) +
  coord_equal()
pl3

pl4 <- ggplot() +
  gg(as(ppxl, "SpatialPixelsDataFrame"), aes(fill = mean4)) +
  coord_equal()
pl4

# -------------------------
# Compare standard errors when fit to each season separately, versus using a spatio-temporal model
# -------------------------

par(mfrow=c(2,2))

ST_sd1 = apply(pred1,1,sd)

plot(ST_sd1 ~ lambda1$sd)
abline(a = 0, b = 1)

# Season 2
pred2 <- generate(fit_INLA,ppxl,
                  formula = as.formula(coordinates ~ Intercept + smooth1 + smooth2)) %>% exp()


ST_sd2 = apply(pred2,1,sd)

plot(ST_sd2 ~ lambda1$sd)
abline(a = 0, b = 1)


plot(ST_sd ~ lambda1$sd)
abline(a = 0, b = 1)

plot(ST_sd ~ lambda1$sd)
abline(a = 0, b = 1)

plot(ST_sd ~ lambda1$sd)
abline(a = 0, b = 1)

mean(ST_sd < lambda1$sd)
