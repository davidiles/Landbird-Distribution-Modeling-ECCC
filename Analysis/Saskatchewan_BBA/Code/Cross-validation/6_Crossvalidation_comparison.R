library(tidyverse)
library(ggplot2)
rm(list=ls())
setwd("C:/Users/IlesD/OneDrive - EC-EC/Iles/Projects/Landbirds/Landbird-Distribution-Modeling-ECCC/Analysis/Saskatchewan_BBA/")

# Fitting to PC only (model "A")
results_A <- readRDS("Output/Crossvalidation/Crossval_results_integrated_PC_fixed.rds") %>%
  dplyr::rename(pred_A = mean_count_val_pred, lppd_A = lppd, AUC_A = AUC, RMSE_A = RMSE)

# Fitting to PC & LT  (model "B")
results_B <- readRDS("Output/Crossvalidation/Crossval_results_integrated_PC_LT_fixed.rds") %>%
  dplyr::rename(pred_B = mean_count_val_pred, lppd_B = lppd, AUC_B = AUC, RMSE_B = RMSE)

# Fitting to PC & LT & SC (model "C")
results_C <- readRDS("Output/Crossvalidation/Crossval_results_integrated_PC_SC_fixed.rds") %>%
  dplyr::rename(pred_C = mean_count_val_pred, lppd_C = lppd, AUC_C = AUC, RMSE_C = RMSE)

results <- full_join(results_A,results_B) %>%
  full_join(results_C) %>%
  na.omit() %>%
  group_by(sp_code) %>%
  summarize_all(mean)

# --------------------------------------------------------
# Plot comparisons based on cross-validation likelihood
# --------------------------------------------------------

results$best <- "PC only"
results$best[results$lppd_B > results$lppd_A & results$lppd_B > results$lppd_C] <- "PC + LT"
results$best[results$lppd_C > results$lppd_B & results$lppd_C > results$lppd_A] <- "PC + SC"

ggplot(data = results)+
  
  geom_vline(xintercept = 0, col = "gray", lwd = 2, alpha = 0.5)+
  geom_segment(aes(y = sp_code, yend = sp_code, x = 0, xend = lppd_C - lppd_A, col = "Effect of including SC"),
               size = 5,
               arrow = arrow(length = unit(0, "cm")))+

  geom_segment(aes(y = sp_code, yend = sp_code, x = 0, xend = lppd_B - lppd_A, col = "Effect of including LT"),
               size = 3,
               arrow = arrow(length = unit(0, "cm")))+
  
  # Summarize best model
  geom_text(data = subset(results, best == "PC only"), aes(x = -200, y = sp_code, label = best), col = "gray", hjust=0)+
  geom_text(data = subset(results, best == "PC + LT"), aes(x = -200, y = sp_code, label = best), col = "black", hjust=0)+
  geom_text(data = subset(results, best == "PC + SC"), aes(x = -200, y = sp_code, label = best), col = "dodgerblue", hjust=0)+
  
  scale_color_manual(values=c("black","dodgerblue"),name = "Which data are included?")+
  theme_bw()+
  xlab("Change in cross-validation likelihood\n\n(relative to PC only model)")+
  ylab("Species")+
  ggtitle("Do checklists improve cross-validation accuracy, compared to analysis of point counts only?\n\nPositive values indicate better performance\n\nBest model is written on left")+
  coord_cartesian(xlim=c(-200,200))

# --------------------------------------------------------
# Plot comparisons based on cross-validation AUC
# --------------------------------------------------------

results$best <- "PC only"
results$best[results$AUC_B > results$AUC_A & results$AUC_B > results$AUC_C] <- "PC + LT"
results$best[results$AUC_C > results$AUC_B & results$AUC_C > results$AUC_A] <- "PC + SC"

ggplot(data = results)+
  
  geom_vline(xintercept = 0, col = "gray", lwd = 2, alpha = 0.5)+
  geom_segment(aes(y = sp_code, yend = sp_code, x = 0, xend = AUC_C - AUC_A, col = "Effect of including LT + SC"),
               size = 5,
               arrow = arrow(length = unit(0, "cm")))+
  
  geom_segment(aes(y = sp_code, yend = sp_code, x = 0, xend = AUC_B - AUC_A, col = "Effect of including LT"),
               size = 3,
               arrow = arrow(length = unit(0, "cm")))+
  
  # Summarize best model
  geom_text(data = subset(results, best == "PC only"), aes(x = -0.05, y = sp_code, label = best), col = "gray", hjust=0)+
  geom_text(data = subset(results, best == "PC + LT"), aes(x = -0.05, y = sp_code, label = best), col = "black", hjust=0)+
  geom_text(data = subset(results, best == "PC + SC"), aes(x = -0.05, y = sp_code, label = best), col = "dodgerblue", hjust=0)+
  
  scale_color_manual(values=c("black","dodgerblue"),name = "Which data are included?")+
  theme_bw()+
  xlab("Change in cross-validation AUC\n\n(relative to PC only model)")+
  ylab("Species")+
  ggtitle("Do checklists improve cross-validation accuracy, compared to analysis of point counts only?\n\nPositive values indicate better performance\n\nBest model is noted in text")+
  coord_cartesian(xlim=c(-0.05,0.05))
