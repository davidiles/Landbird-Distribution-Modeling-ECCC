# ************************************************
# BAYESIAN ANALYSIS / SPECIES DISTRIBUTION MODELS FOR SASKATCHEWAN BREEDING BIRD ATLAS
# 
# 1) A 'processed data package' for analysis is prepared by previous script
# ************************************************

rm(list=ls())

# ------------------------------------------------
# Set working directory
# ------------------------------------------------

# stub <- function() {}
# thisPath <- function() {
#   cmdArgs <- commandArgs(trailingOnly = FALSE)
#   if (length(grep("^-f$", cmdArgs)) > 0) {
#     # R console option
#     normalizePath(dirname(cmdArgs[grep("^-f", cmdArgs) + 1]))[1]
#   } else if (length(grep("^--file=", cmdArgs)) > 0) {
#     # Rscript/R console option
#     scriptPath <- normalizePath(dirname(sub("^--file=", "", cmdArgs[grep("^--file=", cmdArgs)])))[1]
#   } else if (Sys.getenv("RSTUDIO") == "1") {
#     # RStudio
#     dirname(rstudioapi::getSourceEditorContext()$path)
#   } else if (is.null(attr(stub, "srcref")) == FALSE) {
#     # 'source'd via R console
#     dirname(normalizePath(attr(attr(stub, "srcref"), "srcfile")$filename))
#   } else {
#     stop("Cannot find file path")
#   }
# }
# 
# dirname <- thisPath()
dirname <- "C:/Users/IlesD/OneDrive - EC-EC/Iles/Backup/2023-11-09/SDM_ECCC/Analysis/Saskatchewan/Code"
setwd(dirname)

`%!in%` <- Negate(`%in%`)

# ------------------------------------------------
# Load packages
# ------------------------------------------------

require(tidyverse)
require(sf)
require(ggpubr)

# ----------------------------------------------------------------
# Default settings
# ----------------------------------------------------------------

results_path <- "../Output/Crossvalidation/Crossval_results_Default.rds"
if (file.exists(results_path)) results_Default <- readRDS(results_path)

results_Default <- results_Default %>%
  dplyr::rename(AUC_INLA_Default = AUC_INLA_PC,
            RMSE_INLA_Default = RMSE_INLA_PC,
            lppd_INLA_Default = lppd_INLA_PC)

results_summary_Default <- results_Default %>%
  group_by(Species) %>%
  summarize(AUC_INLA = mean(AUC_INLA_Default,na.rm = TRUE),
            RMSE_INLA = mean(RMSE_INLA_Default,na.rm = TRUE),
            lppd_INLA = mean(lppd_INLA_Default, na.rm = TRUE),
            AUC_BRT = mean(AUC_BRT,na.rm = TRUE),
            RMSE_BRT = mean(RMSE_BRT,na.rm = TRUE),
            lppd_BRT = mean(lppd_BRT, na.rm = TRUE)
  ) %>%
mutate(percent_change_RMSE = 100*(RMSE_INLA - RMSE_BRT)/RMSE_BRT)

mean(results_summary_Default$lppd_INLA > results_summary_Default$lppd_BRT)
mean(results_summary_Default$AUC_INLA> results_summary_Default$AUC_BRT)
mean(results_summary_Default$RMSE_INLA < results_summary_Default$RMSE_BRT)

#hist(results_summary_Default$lppd_INLA - results_summary_Default$lppd_BRT, breaks = 20)
#hist(results_summary_Default$AUC_INLA- results_summary_Default$AUC_BRT, breaks = 20)
#hist(log(results_summary_Default$RMSE_INLA/results_summary_Default$RMSE_BRT), breaks = 20)


# Does INLA improve AUC relative to BRT?
fitplot1 <- ggplot(results_summary_Default)+
  
  
  geom_segment(aes(y = Species, yend = Species,
                   x = AUC_BRT, xend = AUC_INLA,
                   col = factor(AUC_BRT > AUC_INLA)),
               size = 2,
               arrow = arrow(length = unit(0.1, "cm")))+
  
  
  theme_bw()+
  xlab("Change in AUC")+
  scale_color_manual(values = c("dodgerblue","orangered"), guide = "none")+
  ggtitle("Cross-validation AUC\n\nINLA relative to BRT")

# Does INLA increase lppd relative to BRT?
fitplot3 <- ggplot(results_summary_Default)+
  
  geom_segment(aes(y = Species, yend = Species,
                   x = lppd_BRT, xend = lppd_INLA,
                   col = factor(lppd_BRT > lppd_INLA)),
               size = 2,
               arrow = arrow(length = unit(0.1, "cm")))+
  
  
  theme_bw()+
  xlab("Change in Deviance")+
  scale_color_manual(values = c("dodgerblue","orangered"), guide = "none")+
  ggtitle("Cross-validation Deviance\n\nINLA relative to BRT")

fitplot <- ggarrange(fitplot1,fitplot3,nrow=1)
print(fitplot)

# ----------------------------------------------------------------
# Compare default to alternative
# ----------------------------------------------------------------

results_path <- "../Output/Crossvalidation/Crossval_results_Default.rds"
if (file.exists(results_path)) results_Default <- readRDS(results_path)

results_Default <- results_Default %>%
  dplyr::rename(AUC_INLA_Default = AUC_INLA_PC,
                RMSE_INLA_Default = RMSE_INLA_PC,
                lppd_INLA_Default = lppd_INLA_PC) %>%
  dplyr::select(Species,Crossval_Fold,lppd_INLA_Default)

results_path <- "../Output/Crossvalidation/Crossval_results_alteredPrior.rds"
if (file.exists(results_path)) results_Alternate <- readRDS(results_path)

results_Alternate <- results_Alternate %>%
  dplyr::rename(AUC_INLA_Alternate = AUC_INLA_PC,
                RMSE_INLA_Alternate = RMSE_INLA_PC,
                lppd_INLA_Alternate = lppd_INLA_PC)%>%
  dplyr::select(Species,Crossval_Fold,lppd_INLA_Alternate)

results_compare <- full_join(results_Default, results_Alternate) %>%
  na.omit() %>%
  group_by(Species) %>%
  summarize(delta_lppd = mean(lppd_INLA_Alternate - lppd_INLA_Default))

results_compare