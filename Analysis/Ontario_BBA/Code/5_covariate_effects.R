# ************************************************
# Visualize covariate effects for each species
# ************************************************

# ------------------------------------------------
# Load packages
# ------------------------------------------------

require(tidyverse)
require(viridis)
require(gridExtra)

rm(list=ls())

# ------------------------------------------------
# Set working directory
# ------------------------------------------------

dirname <- "C:/Users/IlesD/OneDrive - EC-EC/Iles/Projects/Landbirds/Landbird-Distribution-Modeling-ECCC/Analysis/Ontario_BBA/Code"
setwd(dirname)

`%!in%` <- Negate(`%in%`)

# ------------------------------------------------
# Generate plots
# ------------------------------------------------

covariate_effect_table <- read.csv(paste0("../Output/Tables_Summaries/covariate_effect_table.csv"), header = TRUE)

covar <- "PC3"
tmp <- covariate_effect_table[which(covariate_effect_table[,covar] != 0),]
colnum <- which(colnames(tmp) == covar)
ggplot(data = tmp, aes(x = tmp[,colnum], y = pred_q50, ymin = pred_q05, ymax = pred_q95))+
  geom_ribbon(alpha = 0.1, col = "transparent")+
  geom_line()+
  xlab(covar)+
  #scale_y_continuous(trans="log10")+
  ylab("Predicted Count")+
  ggtitle("Effect of Hours-Since-Sunrise")+
  theme_bw()+
  theme(legend.position = "none")+
  facet_wrap(species_name~error_type, scales = "free_y")+
  geom_vline(xintercept = 0, linetype = 2)

