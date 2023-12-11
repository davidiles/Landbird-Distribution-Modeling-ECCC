library(tidyverse)
library(ggplot2)

setwd("C:/Users/IlesD/OneDrive - EC-EC/Iles/Projects/Landbirds/Landbird-Distribution-Modeling-ECCC/Analysis/Saskatchewan_BBA/Code")

results_PConly <- readRDS("../Output/Crossvalidation/Crossval_results_PConly.rds") %>%
  dplyr::select(-mean_count_val_pred)
results_integrated <- readRDS("../Output/Crossvalidation/Crossval_results_integrated_LT.rds") %>%
  dplyr::select(-mean_count_val_pred) %>% rename(lppd_integrated = lppd_integrated_LT)

results <- full_join(results_PConly,results_integrated) %>%
  na.omit() %>%
  group_by(sp_code) %>%
  summarize_all(mean)

ggplot(results, aes(y = sp_code, yend = sp_code, x = 0, xend = lppd_integrated - lppd_PConly))+
  geom_segment(aes(col = factor(lppd_PConly > lppd_integrated)),
               size = 2,
               arrow = arrow(length = unit(0.2, "cm")))+
  
  
  theme_bw()+
  xlab("Likelihood(integrated) - Likelihood(PConly)")+
  scale_color_manual(values = c("dodgerblue","orangered"), guide = "none")+
  ggtitle("Likelihood of integrated model versus PConly model")+
  coord_cartesian(xlim=c(-100,100))

