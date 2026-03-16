rm(list=ls())

library(BAMexploreR)
library(sf)
library(dplyr)
library(ggplot2)
library(terra)
library(ggrepel)

# ------------------------------------------------------------------------------
# Download BAM rasters (only needs to be done once)
# ------------------------------------------------------------------------------

bam_sp <- BAMexploreR::bam_spp_list(version = "v5")
bam_spp_tbl <- BAMexploreR::spp_tbl %>% filter(speciesCode %in% bam_sp)

# Download species prediction rasters
dest <- "../../Data/Spatial/BAMv5/"

years_to_download <- c("2000","2020")

for (sp in bam_sp){
  for (yr in years_to_download){
    
    # Skip if file already exists
    if (file.exists(paste0(dest,sp,"_mosaic_",yr,".tif"))) next
    
    # Download
    BAMexploreR::bam_get_layer(spList = sp,
                               version = "v5",
                               destfile = dest,
                               year = yr)
  }
}

# ------------------------------------------------------------------------------
# Define study region
# ------------------------------------------------------------------------------

# Load study region
bcr <- st_read("../../Data/Spatial/BCR/BCR_Terrestrial_master.shp", quiet = TRUE)

study_region <- bcr %>%
  filter(PROVINCE_S == "ONTARIO") %>%
  st_make_valid() %>%
  mutate(region = paste0("CA-ON-", BCR)) %>%
  group_by(region) %>%
  summarise(do_union = TRUE, .groups = "drop") %>%
  st_make_valid()

# ------------------------------------------------------------------------------
# Load BBS trend estimates (2002 - 2022)
# ------------------------------------------------------------------------------

BBS_trend_summary <- read.csv(
  "../../Data/bbs_trend_estimates/Ontario_Atlas_trends_means_by_atlas_period_2002_2022.csv"
) %>%
  rename(
    BBS_q025 = percent_change_q_0.025,
    BBS_mean = percent_change,
    BBS_q975 = percent_change_q_0.975,
    n_BBS_routes_mean = mean_n_routes
  ) %>%
  select(
    english_name, region, n_routes, n_BBS_routes_mean,
    BBS_mean, BBS_q025, BBS_q975
  ) %>%
  filter(region %in% unique(study_region$region))

# ------------------------------------------------------------------------------
# Calculate BAM-based change estimates
# (only for species shared with BBS to save time)
# ------------------------------------------------------------------------------

# Species shared with BBS
shared_spp<- intersect(BBS_trend_summary$english_name, bam_spp_tbl$commonName)

sp_to_run <- filter(bam_spp_tbl, commonName %in% shared_spp)

# ---- colour ramps
abund_cols <- colorRampPalette(c(
  "#FEFEFE", "#FBF7E2", "#FCF8D0", "#EEF7C2", "#CEF2B0",
  "#94E5A0", "#51C987", "#18A065", "#008C59", "#007F53", "#006344"
))(10)

chg_cols <- colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(101)


bam_summary_list <- list()
k <- 0

for (j in seq_len(nrow(study_region))) {
  
  region_sf <- study_region[j, ]
  region_id <- region_sf$region
  
  message("Processing region: ", region_id)
  
  for (i in seq_len(nrow(sp_to_run))) {
    
    sp <- sp_to_run$speciesCode[i]
    cn <- sp_to_run$commonName[i]
    
    bam_file_2000 <- paste0(dest, sp, "_mosaic_2000.tif")
    bam_file_2020 <- paste0(dest, sp, "_mosaic_2020.tif")
    
    if (!file.exists(bam_file_2000) || !file.exists(bam_file_2020)) next
    
    message("  ", i, "/", nrow(sp_to_run), ": ", cn, " (", sp, ")")
    
    # Load rasters
    bam_rast_2000 <- rast(bam_file_2000)[["mean"]]
    bam_rast_2020 <- rast(bam_file_2020)[["mean"]]
    
    # Reproject 2020 onto 2000 grid
    bam_rast_2020_reproj <- project(bam_rast_2020, bam_rast_2000, method = "bilinear")
    
    # Transform just this BCR polygon to raster CRS
    study_vect <- region_sf %>%
      st_transform(crs(bam_rast_2000)) %>%
      vect()
    
    # Crop/mask
    bam_rast_2000_reg <- crop(bam_rast_2000, study_vect) |> mask(study_vect)
    bam_rast_2020_reg <- crop(bam_rast_2020_reproj, study_vect) |> mask(study_vect)
    
    # Totals
    N2000 <- global(bam_rast_2000_reg, "sum", na.rm = TRUE)[1, 1]
    N2020 <- global(bam_rast_2020_reg, "sum", na.rm = TRUE)[1, 1]
    
    prop_change <- if (isTRUE(N2000 > 0)) (N2020 - N2000) / N2000 else NA_real_
    percent_change <- 100 * prop_change
    
    k <- k + 1
    bam_summary_list[[k]] <- data.frame(
      bam_commonName = cn,
      bam_speciesCode = sp,
      region = region_id,
      bam_N2000 = N2000,
      bam_N2020 = N2020,
      bam_prop_change = prop_change,
      bam_percent_change = percent_change,
      stringsAsFactors = FALSE
    )
  }
}

bam_summary <- bind_rows(bam_summary_list)

saveRDS(bam_summary,"output/bam_change_summary_2000_2020.rds")

# ------------------------------------------------------------------------------
# Code to generate figures
# ------------------------------------------------------------------------------

# # Use common abundance scale for both years within species
# abund_max <- max(
#   values(bam_rast_2000),
#   values(bam_rast_2020_reproj),
#   na.rm = TRUE
# )
# 
# # png(
# #   filename = file.path("../Figures/BAM_maps", paste0(sp, "_bam_change.png")),
# #   width = 1800,
# #   height = 700,
# #   res = 150
# # )
# 
# par(mfrow = c(1, 3))
# 
# plot(
#   bam_rast_2000,
#   col = abund_cols,
#   range = c(0, abund_max),
#   main = paste0(sp, " - 2000")
# )
# 
# plot(
#   bam_rast_2020_reproj,
#   col = abund_cols,
#   range = c(0, abund_max),
#   main = paste0(sp, " - 2020")
# )
# 
# # -------------------------
# # Plot change surface
# # -------------------------
# 
# max_abs <- quantile(abs(values(bam_chg)), 0.995, na.rm = TRUE)
# zlim <- c(-max_abs, max_abs)
# 
# plot(
#   bam_chg,
#   col = chg_cols,   # makes blue = positive, red = negative
#   range = zlim,
#   main = paste0(sp, " - Change (2020 - 2000)")
# )
# dev.off()
# 
# par(mfrow = c(1, 1))
# 