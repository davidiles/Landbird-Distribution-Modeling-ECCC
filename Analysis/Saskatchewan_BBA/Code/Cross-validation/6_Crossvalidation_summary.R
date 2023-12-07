library(tidyverse)
library(ggplot2)

setwd("C:/Users/IlesD/OneDrive - EC-EC/Iles/Projects/Landbirds/Landbird-Distribution-Modeling-ECCC/Analysis/Saskatchewan_BBA/Code")

results_PConly <- readRDS("../Output/Crossvalidation/Crossval_results_PConly.rds") %>%
  dplyr::select(-mean_count_val_pred)
results_integrated <- readRDS("../Output/Crossvalidation/Crossval_results_integrated.rds")%>%
  dplyr::select(-mean_count_val_pred)

results <- full_join(results_PConly,results_integrated) %>%
  na.omit() %>%
  group_by(sp_code) %>%
  summarize_all(mean)

ggplot(results, aes(y = sp_code, yend = sp_code, x = 0, xend = lppd_integrated - lppd_PConly))+
  geom_segment(aes(col = factor(lppd_PConly > lppd_integrated)),
               size = 2,
               arrow = arrow(length = unit(0.2, "cm")))+
  
  
  theme_bw()+
  xlab("Cross-validation likelihood")+
  scale_color_manual(values = c("dodgerblue","orangered"), guide = "none")+
  ggtitle("Difference in cross-validation likelihood")
