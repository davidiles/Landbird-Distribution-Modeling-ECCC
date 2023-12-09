
results <- data.frame()

for (run in 1:100){
  
  set.seed(run)
  
  # Simulate data
  log_intercept = log(runif(1,0,5))
  log_slope = runif(1,-1,1)
  n = 500
  x = rnorm(n)
  mu = 4
  size = 0.1
  
  mu = exp(log_intercept + log_slope*x)
  
  y = rnbinom(n,mu = mu, size = size)
  
  dat = data.frame(y = y,x=x)
  
  # Fit model
  model_components = as.formula('~ Intercept(1)+ Beta1(1,model="linear", mean.linear = 0, prec.linear = 1)')
  model_formula = as.formula('y ~Intercept +Beta1*x')
  fit <- bru(components = model_components,
             like(family = "nbinomial",
                  formula = model_formula,
                  data = dat))
  summary(fit)
  
  
  # Validation data
  preddat <- data.frame(x = rnorm(500))
  
  # Predictions for each location in validation data
  pred <- generate(fit,
                   preddat,
                   formula =  as.formula(' ~Intercept +Beta1*x'),
                   n.samples = 1000)
  
  size_fit <- fit$summary.hyperpar$'0.5quant'[1]
  
  pred <- exp(pred)
  for (i in 1:ncol(pred)) pred[,i] <- rnbinom(nrow(pred),mu = pred[,i], size = size_fit)
  
  med <- apply(pred,1,function(x) median(x))
  mean <- apply(pred,1,function(x) mean(x))
  
  # True values in validation data
  preddat$mu = exp(log_intercept + log_slope*preddat$x)
  preddat$y = rnbinom(nrow(preddat),mu = preddat$mu, size = size)
  
  # Store results (predicted sum, versus true sum)
  results <- rbind(results, data.frame(run = run,
                                       sum_pred = mean(apply(pred,2,sum)),
                                       lci = quantile(apply(pred,2,sum),0.025),
                                       uci = quantile(apply(pred,2,sum),0.975),
                                       sum_true = sum(preddat$y)))
  
  
  plot <- ggplot(data = results, aes(x = sum_true, y = sum_pred, ymin = lci, ymax = uci))+
    geom_errorbar(width=0)+
    geom_point()+
    geom_abline(intercept = 0, slope = 1)+
    theme_bw()
  print(plot)
  
}


mean(results$sum_pred - results$sum_true)
cov = mean(results$lci <= results$sum_true & results$uci >= results$sum_true)
cov
