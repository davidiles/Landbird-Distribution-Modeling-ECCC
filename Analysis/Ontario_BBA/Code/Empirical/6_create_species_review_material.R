# ===============================================================================
# BAYESIAN ANALYSIS / SPECIES DISTRIBUTION MODELS FOR ONTARIO BREEDING BIRD ATLAS
# ===============================================================================

library(tidyverse)
library(sf)
library(terra)
library(viridis)

rm(list=ls())

setwd("C:/Users/IlesD/OneDrive - EC-EC/Iles/Projects/Landbirds/Landbird-Distribution-Modeling-ECCC/")

# ------------------------------------------------------------------------------
# -------- Load data / results
# ------------------------------------------------------------------------------

# Load prepared data
dat <- readRDS(file = "Analysis/Ontario_BBA/Data_Cleaned/data_ready_for_analysis.rds")
all_surveys <- dat$all_surveys
species_to_model <- dat$species_to_model
counts <- dat$counts
study_boundary <- dat$study_boundary
species_to_model <- dat$species_to_model

# BCR boundaries
allBCR <- st_read("Data/Spatial/BCR/BCR_Terrestrial_master.shp") %>%
  st_transform(st_crs(all_surveys))

NorthAmerica = allBCR %>% 
  subset(BCR != 100) %>%
  subset(COUNTRY %in% c("USA","CANADA")) %>%
  group_by(PROVINCE_S) %>%
  summarise()

BCR <- st_read("Data/Spatial/BCR/BCR_Terrestrial_master.shp") %>%
  subset(PROVINCE_S == "ONTARIO") %>%
  st_transform(st_crs(study_boundary)) %>%
  st_make_valid() %>%
  st_buffer(0.01) %>%
  group_by(BCR,BCRNAME,PROVINCE_S) %>%
  summarise() %>%
  st_transform(st_crs(all_surveys))

# hexagon grid across study area
hex <- st_make_grid(study_boundary, cellsize = 40, square=TRUE) %>% 
  st_as_sf() %>%
  dplyr::rename(geometry = x)
hex <- hex[st_intersects(hex, study_boundary, sparse = FALSE), ]
hex$hexid <- 1:nrow(hex)

hex_pts <- hex %>% st_point_on_surface()

nearest_bcr_idx <- st_nearest_feature(hex_pts, BCR)

hex <- hex %>%mutate(BCR = BCR$BCR[nearest_bcr_idx])

table(hex$BCR, useNA = "always")

ggplot()+geom_sf(data = hex, aes(fill = factor(BCR)), col = "transparent")+theme_bw()

# ------------------------------------------------------------------------------
# PART 1: SUMMARY OF OVERALL ATLAS SAMPLING EFFORT
# ------------------------------------------------------------------------------

n_survey_per_sq <- dat$all_surveys %>%
  as.data.frame() %>%
  group_by(Atlas,square_atlas,BCR) %>%
  summarize(n = n())

n_survey_per_sq %>%
  group_by(BCR,Atlas) %>%
  summarize(ntot = sum(n),
            nsq = n(),
            max_n_per_sq = max(n),
            min_n_per_sq = min(n),
            mean_n_per_sq = mean(n),
            med_n_per_sq = median(n)) %>%
  as.data.frame()

# Very similar effort in southern BCRs

# Much higher spatial coverage (nsq) in Atlas 3 vs Atlas 2 (especially BCR 7)
# Higher total number of point counts in BCR 7 in atlas 3 vs atlas 2
# Lower total number of point counts in BCR 8 in atlas 3 vs atlas 2
# Much higher mean n per square in atlas 2 (nearly double)
# Much higher median n per square in atlas 2 (4-5 times higher)

# Atlas 2 placed more emphasis on high density of point counts per square
# Atlas 3 placed more emphasis on high spatial coverage at expense of nsurvey per square
#  - this reflects CWS-ONT emphasis on improving spatial coverage, learning new things about the north, etc.

# Maps of survey effort
hex_summary <- st_join(
  hex,
  all_surveys[, c("Atlas","DayOfYear","Hours_Since_Sunrise","Survey_Type")],
  join = st_intersects,
  left = TRUE
) %>%
  group_by(geometry, Atlas,BCR) %>%
  summarise(
    n_surveys = sum(!is.na(Atlas)),
    n_ARU = sum(Survey_Type == "ARU"),
    n_PC = sum(Survey_Type == "Point_Count"),
    mean_HSS = mean(Hours_Since_Sunrise),
    mean_dsj15 = mean(DayOfYear),
    .groups = "drop"
  ) %>%
  na.omit()

head(hex_summary)

focal_BCR <- c(7,8)
dat_for_plot <- hex_summary %>% subset(BCR %in% focal_BCR)
lim <- max(dat_for_plot$n_surveys)
ggplot()+
  
  geom_sf(data = study_boundary, fill = "transparent")+
  geom_sf(data = BCR %>% subset(BCR %in% focal_BCR) , fill = "transparent")+
  
  geom_sf(data = dat_for_plot, 
          aes(fill = n_surveys),
          col = "transparent") +
  #geom_sf(data = dat_for_plot %>% dplyr::select(-Atlas),
  #        fill = "transparent",
  #        col = "gray90")+
  scale_fill_gradientn(colours = c("gray80","black"),
                       #trans = "log10",
                       limits = c(1,lim))+
  facet_grid(.~Atlas)+
  theme_bw()

# ------------------------------------------------------------------------------
# -------- Species-level results summaries
# ------------------------------------------------------------------------------
model_summaries <- readRDS("Analysis/Ontario_BBA/Output/Tables_Summaries/model_summaries.rds")
#HSS_DOY_summaries <- readRDS("Analysis/Ontario_BBA/Output/Tables_Summaries/HSS_DOY_summaries.rds")

sp_list <- names(model_summaries)

sp_list = c("Lincoln's Sparrow","Olive-sided Flycatcher")

#for (i in 1:length(sp_list)){

#sp_english <- sp_list[i]

sp_english <- "Palm Warbler"

sp_filename <- sp_english %>%
  str_to_lower() %>%                       # lower case
  str_replace_all("[^a-z0-9]+", "_") %>%   # replace non-alphanumeric with _
  str_replace_all("^_|_$", "")             # trim leading/trailing underscores

print(sp_english)

sp_code <- subset(species_to_model, english_name == sp_english)$species_id[1] %>% as.character()

# ----------------------------------------------------------------------------
# -------- Maps of relative abundance
# ----------------------------------------------------------------------------

# Prepare species survey data
spdat <- all_surveys %>% 
  mutate(count = counts[[sp_code]],
         days_since_june15 = DayOfYear - 166,
         BCR_factor = as.numeric(factor(BCR))) 


# Summarize survey and species-level information
spdat_summary <- st_join(
  hex,
  spdat[, c("Atlas","days_since_june15","Hours_Since_Sunrise","Survey_Type","count")],
  join = st_intersects,
  left = TRUE
) %>%
  group_by(geometry, Atlas,BCR) %>%
  summarise(
    n_surveys = sum(!is.na(Atlas)),
    n_ARU = sum(Survey_Type == "ARU"),
    n_PC = sum(Survey_Type == "Point_Count"),
    mean_HSS = mean(Hours_Since_Sunrise),
    mean_dsj15 = mean(days_since_june15),
    mean_count = mean(count),
    .groups = "drop"
  ) %>%
  na.omit()


# Summary plots
focal_BCR <- c(7,8,12,13) #c(7,8)
dat_for_plot <- spdat_summary %>% 
  subset(BCR %in% focal_BCR) %>%
  subset(n_surveys >= 1)

cent <- st_centroid(dat_for_plot)

max_val <- quantile(subset(cent, mean_count > 0)$mean_count, 0.95) %>%
  round(2)
cent$mean_count_trunc <- cent$mean_count
cent$mean_count_trunc[cent$mean_count_trunc < 0.005 & cent$mean_count_trunc > 0] <- 0.005

cent$mean_count_trunc[cent$mean_count_trunc > max_val] <- max_val

breaks <-  c(0.005, 0.1, 0.25, 0.5, 1)
labels <- c("<0.005", "0.1", "0.25", "0.5", ">1")

if (max_val > max(breaks)){
  breaks = c(breaks,max_val)
  labels = c(labels,paste0(">",max_val))
} else{
  
  labels <- labels[-which(breaks>max_val)]
  breaks <- breaks[-which(breaks>max_val)]
  breaks <- c(breaks,max_val)
  labels <- c(labels,paste0(">",max_val))
}

colscale <- c("#FEFEFE", "#FBF7E2", "#FCF8D0", "#EEF7C2", "#CEF2B0",
              "#94E5A0", "#51C987", "#18A065", "#008C59", "#007F53", "#006344")

study_bbox <- st_bbox(study_boundary)
dat_plot <- ggplot() +
  geom_sf(data = NorthAmerica, 
          fill = "gray99", col = "gray90") +
  geom_sf(data = BCR, 
          fill = "gray99", col = "gray90") +
  geom_sf(data = study_boundary, 
          fill = "transparent", col = "gray75") +
  
  geom_sf(data = cent, aes(size = n_surveys), 
          col = "gray80",
          fill = "white",
          shape = 21) +
  geom_sf(data = cent %>% subset(mean_count_trunc > 0), 
          aes(fill = mean_count_trunc, 
              size = n_surveys),
          shape = 21,
          col = "gray80") +
  
  scale_fill_gradientn(
    colours = colscale[3:length(colscale)],
    name = "Mean\nCount",
    # Set breaks and labels for legend
    breaks = breaks,
    labels = labels
  ) +
  
  scale_size_area(name = "# surveys",
                  trans = "sqrt",
                  max_size = 3.5,
                  breaks = c(1, 50, 100, 250, 500, 750),
                  labels = c("1", "50", "100", "250", "500", "750")) +
  
  coord_sf(
    xlim = c(study_bbox["xmin"], study_bbox["xmax"]),
    ylim = c(study_bbox["ymin"], study_bbox["ymax"])
  ) +
  
  facet_grid(.~Atlas, labeller = label_both) +
  theme_bw() +
  theme(
    panel.background = element_rect(fill = "#F5FBFF", color = NA),
    panel.grid = element_blank()
  )+
  ggtitle(sp_english)
print(dat_plot)
# 12 x 6

# # ----------------------------------------------------------------------------
# # -------- DOY and HSS results
# # ----------------------------------------------------------------------------
# 
# pred_df <- HSS_DOY_summaries[[sp_english]]
# 
# # Simple method for calculating peak DOY
# peak_DOY <- pred_df %>%
#   filter(effect == "DOY" & days_since_june15 >= -15 & days_since_june15 <= 25) %>%
#   group_by(BCR_factor) %>%
#   slice_max(`50%`, n = 1, with_ties = FALSE) %>%
#   ungroup() %>%
#   transmute(
#     BCR_factor,
#     days_since_june15_peak = days_since_june15
#   )
# 
# # Plots
# HSS_plot <- ggplot(subset(pred_df, effect == "HSS"), aes(x = Hours_Since_Sunrise, y = `50%`, ymin = `2.5%`, ymax = `97.5%`))+
#   geom_ribbon(fill = "gray90") +
#   geom_line() +
#   geom_vline(xintercept = 0, col = "dodgerblue", linewidth = 2, alpha = 0.5)+
#   theme_bw()+
#   xlab("Hours Since Sunrise")+
#   ylab("Relative abundance")+
#   ggtitle(sp_english)
# 
# print(HSS_plot)
# 
# DOY_plot <- ggplot(subset(pred_df, effect == "DOY" & days_since_june15 >= -15 & days_since_june15 <= 25), aes(x = days_since_june15, y = `50%`, ymin = `2.5%`, ymax = `97.5%`))+
#   geom_ribbon(fill = "gray90") +
#   geom_line() +
#   geom_vline(xintercept = -7, col = "dodgerblue", linewidth = 2, alpha = 0.5)+
#   facet_wrap(BCR_factor~., scales = "free_y")+
#   theme_bw()+
#   xlab("Days since June 15")+
#   ylab("Relative abundance")+
#   ggtitle(sp_english)
# 
# print(DOY_plot)
# 
# # ----------------------------------------------------------------------------
# # -------- Species habitat effects
# # ----------------------------------------------------------------------------
# 
# mod <- model_summaries[[sp_english]]
# 
# fixef <- mod$summary.fixed[4:nrow(mod$summary.fixed),]
# fixef <- as.data.frame(fixef) %>%
#   mutate(effect = rownames(fixef))
# 
# ggplot(fixef,aes(y = effect,
#                  x = mean, 
#                  xmin = `0.025quant`,
#                  xmax = `0.975quant`))+
#   geom_vline(xintercept = 0, linetype = 2)+
#   geom_point()+
#   geom_errorbarh(height=0)+
#   theme_bw()

#}


# # ------------------------------------------------------------------------------
# # -------- Covariate maps
# # ------------------------------------------------------------------------------
# 
# create_raster <- function(grid, covar_name,res) {
#   
#   v <- grid[,covar_name] %>%
#     vect()
#   
#   r_template <- rast(v, res = res)
#   r <- rasterize(v, r_template, field = covar_name, fun = mean)
#   names(r) <- covar_name
#   
#   return(r)
# }
# 
# covars <- c("insect_broadleaf",
#             "insect_needleleaf",
#             "road",
#             "water_river",
#             "prec",
#             "tmax",
#             "urban_2",
#             "urban_3",
#             "lc_1",  # Evergreen Needleaf Forest
#             "lc_4",  # Deciduous Broadleaf Forest
#             "lc_5",  # Mixed Forest
#             "lc_8",  # Woody Savanna; Tree cover 30-60% (canopy >2m)
#             "lc_9",  # Savanna; Tree cover 10-30% (canopy >2m)
#             "lc_10", # Grassland
#             "lc_11", # Permanent Wetland
#             "lc_12", # Cropland
#             "lc_14",  # Cropland / Natural Vegetation Mosaic'
#             "lc_17"  # Water
# )
# 
# covar_names <- c("Broadleaf Defoliation",
#                           "Needleleaf Defoliation",
#                           "Road density",
#                           "River density",
#                           "Precipitation",
#                           "Max Temperature",
#                           "Rural Settlement (50-100 people/km^2)",
#                           "Urban Settlement (>500 people/km^2)",
#                           "Evergreen Needleleaf",
#                           "Deciduous Broadleaf",
#                           "Mixed Forest",
#                           "Woody Savanna (Tree cover 30-60%; canopy >2m)",
#                           "Savanna (Tree cover 10-30%; canopy >2m)",
#                           "Grassland",
#                           "Permanent Wetland",
#                           "Cropland",
#                           "Cropland / Natural Veg Mosaic",
#                           "Water")
# 
# covar_df <- data.frame(covar_code = covars,
#                        covar_name = covar_names)
# 
# grid_OBBA2 <- dat$grid_OBBA2
# grid_OBBA3 <- dat$grid_OBBA3
# 
# for (i in 1:nrow(covar_df)){
#   
#   cov <- covar_df$covar_code[i]
#   cov_name <- covar_df$covar_name[i]
#   
#   r_OBBA2 <- create_raster(grid = grid_OBBA2,
#                            covar_name = cov,
#                            res = 1.00001)
#   r_OBBA3 <- create_raster(grid = grid_OBBA3,
#                            covar_name = cov,
#                            res = 1.00001)
#   
#   r_delta <- r_OBBA3 - r_OBBA2
#   
#   
#   pdf(paste0("Analysis/Ontario_BBA/Output/Covariate_Maps/",cov,".pdf"),height = 6, width = 15)
#   par(mfrow = c(1, 3), oma = c(0, 0, 2, 0))
#   
#   # Covariates themselves
#   colpal <- colorRampPalette(c("white",rev(viridis(5)), "black"))
#   zlim <- range(c(values(r_OBBA2),values(r_OBBA3)),na.rm = TRUE)
#   zlim[1] <- zlim[1] - 0.01
#   zlim[2] <- zlim[2] + 0.01
#   
#   plot(r_OBBA2, col = colpal(100), range = zlim, maxcell = Inf,
#        main = "Atlas 2")
#   plot(vect(dat$study_boundary), add = TRUE, border = "gray80")
#   plot(r_OBBA3, col = colpal(100), range = zlim, maxcell = Inf,
#        main = "Atlas 3")
#   plot(vect(dat$study_boundary), add = TRUE, border = "gray80")
#   # Change in covariate between atlases
#   colpal <- colorRampPalette(c("red","white","blue"))
#   zlim <- max(abs(values(r_delta)),na.rm = TRUE)
#   plot(r_delta, col = colpal(100), range = c(-zlim,zlim), maxcell = Inf,
#        main = "Change (Atlas 3 - Atlas 2)")
#   plot(vect(dat$study_boundary), add = TRUE, border = "gray80")
#   par(mfrow=c(1,1))
#   # Overarching title for the species
#   mtext(cov_name, outer = TRUE, line = 0.5, cex = 1.5)
#   mtext(cov, outer = TRUE, line = -1, cex = 1.5)
#   dev.off()
#   
# }
