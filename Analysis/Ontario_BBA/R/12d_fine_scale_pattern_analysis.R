# ============================================================
# 09_make_species_figures.R
#
# Purpose:
#   Generate PNG figures per species using the per-species
#   prediction summary objects created in script 08.
#
# Inputs:
#   data_clean/birds/data_ready_for_analysis.rds
#   data_clean/model_output/predictions/<sp>.rds  (many)
#   + external shapefiles for BCR, water, (optional) atlas squares
#
# Outputs:
#   data_clean/model_output/figures/<sp>_*.png
# ============================================================

suppressPackageStartupMessages({
  library(sf)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(dplyr)
  library(lubridate)
  library(tidyr)
  library(mgcv)
  library(naturecounts)
})

# ------------------------------------------------------------
# Config
# ------------------------------------------------------------

in_data <- "data_clean/birds/data_ready_for_analysis.rds"

# Shapefile inputs
in_bcr <- "../../Data/Spatial/BCR/BCR_Terrestrial_master.shp"
in_water <- "data_clean/spatial/water_filtered.shp"
in_atlas_squares <- "../../Data/Spatial/National/AtlasSquares/NationalSquares_FINAL.shp"  # optional

# ------------------------------------------------------------
# Load base data
# ------------------------------------------------------------
NC_names <- naturecounts::search_species()

dat <- readRDS(in_data)

# raw data
all_surveys <- dat$all_surveys
counts      <- dat$counts

# prediction grid / study area
grid2 <- dat$grid_OBBA2
grid3 <- dat$grid_OBBA3

study_boundary <- dat$study_boundary %>% sf::st_as_sf()

# Ensure consistent CRS
crs_use <- st_crs(grid2)
study_boundary <- st_transform(study_boundary, crs_use)
grid3 <- st_transform(grid3, crs_use)
grid2 <- st_transform(grid2, crs_use)

# BCR outlines
bcr_sf <- st_read(in_bcr, quiet = TRUE) %>%
  st_make_valid() %>%
  st_transform(crs_use) %>%
  st_intersection(study_boundary) %>%
  st_collection_extract("POLYGON") %>%   # <<--- drop lines/points
  dplyr::filter(PROVINCE_S %in% c("ONTARIO", "ON", "Ontario") | is.na(PROVINCE_S)) %>%
  dplyr::select(BCR, BCRNAME, PROVINCE_S) %>%
  group_by(BCR, BCRNAME) %>%
  summarise(geometry = st_union(geometry), .groups = "drop") %>%
  st_make_valid()

# Water
water_sf <- NULL
if (file.exists(in_water)) {
  water_sf <- st_read(in_water, quiet = TRUE) %>%
    st_make_valid() %>%
    st_transform(crs_use)
}

# ------------------------------------------------------------
# Custom polygon (given as lat, lon)
# ------------------------------------------------------------

# vertices as given (lat, lon)
verts_latlon <- tribble(
  ~lat,        ~lon,
  57, -88.144358,
  55.6, -81,
  50, -87.5,
  50, -81
)

# points in EPSG:4326 (x=lon, y=lat)
verts_pts_ll <- st_as_sf(verts_latlon, coords = c("lon", "lat"), crs = 4326)

# transform to your analysis CRS and compute convex hull
custom_poly <- verts_pts_ll |>
  st_transform(crs_use) |>
  st_union() |>
  st_convex_hull() |>
  st_make_valid()

# ensure it's an sfc (often fine as-is, but this is explicit)
custom_poly <- st_sfc(custom_poly, crs = crs_use) %>%
  st_intersection(study_boundary)

# Plot focal area
ggplot()+
  geom_sf(data = study_boundary, fill = "transparent")+
  geom_sf(data = bcr_sf, fill = "transparent", col = "gray70")+
  geom_sf(data = custom_poly, col = "red", fill = "transparent", linewidth = 2)+
  theme_bw()

# ---- water layer
water_sf_sub <- st_transform(water_sf, crs_use) |> 
  dplyr::filter(lengths(st_within(geometry, custom_poly)) > 0)

# ------------------------------------------------------------
# Loop through species and identify which ones have the largest differences
# in mean counts between atlases in that custom polygon
# ------------------------------------------------------------

all_results <- vector("list", length(dat$species_to_model$species_id))

for (i in seq_along(dat$species_to_model$species_id)){

  sp_english <- dat$species_to_model$english_name[i]
  sp_code <- dat$species_to_model$species_id[i]
  if (nrow(subset(NC_names, species_id == sp_code))==0){
    message("Missing species ",sp_code)
    next
  }

  sp_dat <- all_surveys %>%
    mutate(
      count = counts[[sp_code]],
      days_since_june15 = DayOfYear - 166,
      BCR_factor = as.numeric(factor(BCR)),
      Atlas3 = ifelse(Atlas == "OBBA2", 0, 1)
    ) %>%
    st_transform(crs_use) %>%
    dplyr::filter(lengths(st_within(geometry, custom_poly)) > 0)

  rel_counts <- sp_dat %>%
    st_drop_geometry() %>%
    group_by(Atlas) %>%
    summarize(
      n_survey = n(),
      n_det    = sum(count > 0, na.rm = TRUE),
      mean_count = mean(count, na.rm = TRUE),
      prop_det   = mean(count > 0, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(sp_english = sp_english,
           sp_code    = sp_code)

  all_results[[i]] <- rel_counts

}

rel_counts_long <- bind_rows(all_results)

rel_counts_wide <- rel_counts_long %>%
  pivot_wider(
    names_from = Atlas,
    values_from = c(n_survey, n_det, mean_count, prop_det),
    names_sep = "__"
  ) %>%
  mutate(
    d_mean_count = mean_count__OBBA3 - mean_count__OBBA2,
    d_prop_det   = prop_det__OBBA3   - prop_det__OBBA2,
    abs_d_mean_count = abs(d_mean_count),
    abs_d_prop_det   = abs(d_prop_det),
    det20_any = pmax(n_det__OBBA2, n_det__OBBA3, na.rm = TRUE) >= 20
  )

# Focus on species with at least 20 detections in each atlas, where there is a greater than 50% difference
# in mean counts between atlases
species_to_examine <- rel_counts_wide %>%
  subset(!is.na(sp_english)) %>%
  mutate(
    # >=20 detections in either atlas
    det20_any = pmax(n_det__OBBA2, n_det__OBBA3, na.rm = TRUE) >= 20,

    # relative difference in mean counts (symmetric % change)
    # = |x3 - x2| / mean(x2, x3)
    rel_diff_mean = log(mean_count__OBBA3/mean_count__OBBA2)
  ) %>%
  filter(det20_any, rel_diff_mean > 0.50) %>%
  arrange(desc(rel_diff_mean))

species_to_examine %>% as.data.frame()

# ------------------------------------------------------------
# ------------------------------------------------------------
# Examine available data within focal area for study species
# ------------------------------------------------------------
# ------------------------------------------------------------

sp_english <- "Connecticut Warbler"
sp_file <- sp_english %>%
  str_to_lower() %>%
  str_replace_all("[^a-z0-9]+", "_") %>%
  str_replace_all("^_|_$", "")

sp_code <- subset(dat$species_to_model, english_name == sp_english)$species_id

sp_dat <- all_surveys %>%
  mutate(
    count = counts[[sp_code]],
    days_since_june15 = DayOfYear - 166,
    BCR_factor = as.numeric(factor(BCR)),
    Atlas3 = ifelse(Atlas == "OBBA2", 0, 1)
  )

# ---- surveys
# If sp_dat is already sf, just transform it.
if (inherits(sp_dat, "sf")) {
  
  sp_dat_sub <- st_transform(sp_dat, crs_use) |>
    dplyr::filter(lengths(st_within(geometry, custom_poly)) > 0)
  
} else {
  
  sp_dat_sf <- st_as_sf(
    sp_dat,
    coords = c("Longitude", "Latitude"),  # change if needed
    crs = 4326,
    remove = FALSE
  ) |> st_transform(crs_use)
  
  sp_dat_sub <- sp_dat_sf |>
    dplyr::filter(lengths(st_within(geometry, custom_poly)) > 0)
}

# ------------------------------------------------------------
# Summarize differences between Atlas 2 and 3 in terms of survey counts / timing
# ------------------------------------------------------------

sampling_summary <- sp_dat_sub %>%
  st_drop_geometry() %>%
  mutate(
    year = year(Date_Time),
    hour = hour(Date_Time) + minute(Date_Time)/60
  ) %>%
  group_by(Atlas) %>%
  summarise(
    n_surveys = n(),
    n_squares = n_distinct(square_id),
    frac_ARU  = mean(ARU == 1, na.rm = TRUE),
    n_svy_detected = sum(count>0,na.rm = TRUE),
    mean_svy_count = mean(count, na.rm = TRUE),
    prop_svy_detected = mean(count > 0, na.rm = TRUE),
    mean_duration_min = mean(Survey_Duration_Minutes, na.rm = TRUE),
    median_DOY = median(DayOfYear, na.rm = TRUE),
    q10_DOY = quantile(DayOfYear, 0.10, na.rm = TRUE),
    q90_DOY = quantile(DayOfYear, 0.90, na.rm = TRUE),
    median_HSS = median(Hours_Since_Sunrise, na.rm = TRUE),
    q10_HSS = quantile(Hours_Since_Sunrise, 0.10, na.rm = TRUE),
    q90_HSS = quantile(Hours_Since_Sunrise, 0.90, na.rm = TRUE),
    .groups = "drop"
  )

sampling_summary %>% as.data.frame()

# DOY
ggplot(sp_dat_sub %>% st_drop_geometry(), aes(x = DayOfYear)) +
  geom_density() +
  facet_grid(Atlas ~ .) +
  theme_bw()

# Hours since sunrise
ggplot(sp_dat_sub %>% st_drop_geometry(), aes(x = Hours_Since_Sunrise)) +
  geom_density() +
  facet_grid(Atlas ~ .) +
  theme_bw()

# ------------------------------------------------------------
# Plot *all* data for that species
# ------------------------------------------------------------

ggplot() +
  geom_sf(data = custom_poly, fill = NA, linewidth = 0.6, col = "red") +
  
  geom_sf(data = water_sf_sub %>% st_buffer(100),
          fill = "#d6eaf7",
          col = NA) +
  
  # Plot all survey locations
  geom_sf(data = sp_dat_sub,
          col = "gray50",
          size = 2,
          shape = 18,
          stroke = 1) +
  
  # Overlay detection locations
  geom_sf(
    data = sp_dat_sub %>% filter(count > 0),
    shape = 4,
    stroke = 1,
    size = 2,
    aes(color = factor(ARU))
  ) +
  scale_color_manual(values = c("black","blue"), name = "ARU")+
  coord_sf(expand = FALSE) +
  facet_grid(. ~ Atlas,
             labeller = labeller(
               Atlas = c(
                 OBBA2 = "Atlas 2",
                 OBBA3 = "Atlas 3"
               )
             )) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(size = 12, face = "bold")
  )+
  ggtitle(sp_english)


# ------------------------------------------------------------
# Summarize differences between Atlas 2 and 3 in terms of covariate sampling
# ------------------------------------------------------------

covars_to_check <- c(
  "on_river",
  "on_road",
  "urban_3",
  "lc_1S","lc_1N",
  "lc_4S","lc_4N",
  "lc_5S","lc_5N",
  "lc_8S","lc_8N",
  "lc_9S","lc_9N",
  "lc_10S","lc_10N",
  "lc_11","lc_12","lc_14","lc_17",
  "insect_broadleaf","insect_needleleaf"
)
covars_to_check <- intersect(covars_to_check, names(sp_dat_sub))

# Covariate differences (standardized mean difference)
smd_one <- function(x2, x3) {
  x2 <- x2[is.finite(x2)]
  x3 <- x3[is.finite(x3)]
  if (length(x2) < 5 || length(x3) < 5) return(NA_real_)
  m2 <- mean(x2); m3 <- mean(x3)
  s2 <- sd(x2);   s3 <- sd(x3)
  sp <- sqrt((s2^2 + s3^2) / 2) # mean of variances
  if (!is.finite(sp) || sp == 0) return(NA_real_)
  (m3 - m2) / sp   # OBBA3 minus OBBA2 divided by mean variance
}

covar_compare <- covars_to_check %>%
  purrr::map_dfr(function(cv) {
    x2 <- sp_dat_sub[[cv]][sp_dat_sub$Atlas == "OBBA2"]
    x3 <- sp_dat_sub[[cv]][sp_dat_sub$Atlas == "OBBA3"]
    tibble(
      covariate  = cv,
      OBBA2_mean = mean(x2, na.rm = TRUE),
      OBBA3_mean = mean(x3, na.rm = TRUE),
      SMD        = smd_one(x2, x3)
    )
  }) %>%
  arrange(desc(abs(SMD)))

covar_compare %>% as.data.frame() %>% na.omit()

# ------------------------------------------------------------
# RESTRICT COMPARISON TO SHARED FOOTPRINT
# ------------------------------------------------------------

buffer_m <- 500  # Distance

# Create buffer around atlas 2 points
a2_pts <- sp_dat_sub %>% filter(Atlas == "OBBA2")
a3_pts <- sp_dat_sub %>% filter(Atlas == "OBBA3")

# For each A2 point: which A3 points are within buffer?
a2_to_a3 <- st_is_within_distance(a2_pts, a3_pts, dist = buffer_m)

# Keep only A2 points that have at least one nearby A3
a2_keep <- lengths(a2_to_a3) > 0
a2_matched <- a2_pts[a2_keep, ]

# Now keep A3 points that are within buffer of at least one kept A2 point
a3_to_a2keep <- st_is_within_distance(a3_pts, a2_matched, dist = buffer_m)
a3_keep <- lengths(a3_to_a2keep) > 0

a3_matched <- a3_pts[a3_keep, ]

# Combined "shared footprint" dataset (matched A2 + nearby A3)
shared_sf <- bind_rows(a2_matched, a3_matched) %>%
  mutate(shared_pairwise = TRUE)

# 4) Evidence split: "shared footprint" (near_A2) vs not
evidence_by_shared <- shared_sf %>%
  st_drop_geometry() %>%
  group_by(Atlas, ARU) %>%
  summarise(
    n_surveys = n(),
    n_squares = n_distinct(square_id),
    frac_ARU  = mean(ARU == 1, na.rm = TRUE),
    n_svy_detected = sum(count>0,na.rm = TRUE),
    mean_svy_count = mean(count, na.rm = TRUE),
    prop_svy_detected = mean(count > 0, na.rm = TRUE),
    mean_duration_min = mean(Survey_Duration_Minutes, na.rm = TRUE),
    median_DOY = median(DayOfYear, na.rm = TRUE),
    q10_DOY = quantile(DayOfYear, 0.10, na.rm = TRUE),
    q90_DOY = quantile(DayOfYear, 0.90, na.rm = TRUE),
    median_HSS = median(Hours_Since_Sunrise, na.rm = TRUE),
    q10_HSS = quantile(Hours_Since_Sunrise, 0.10, na.rm = TRUE),
    q90_HSS = quantile(Hours_Since_Sunrise, 0.90, na.rm = TRUE),
    .groups = "drop"
  )

evidence_by_shared %>% as.data.frame() 

ggplot() +
  geom_sf(data = custom_poly, fill = NA, linewidth = 0.6, col = "red") +
  
  geom_sf(data = water_sf_sub %>% st_buffer(100),
          fill = "#d6eaf7",
          col = NA) +
  
  # Plot all survey locations
  geom_sf(data = shared_sf,
          col = "gray50",
          size = 2,
          shape = 18,
          stroke = 2.2) +
  
  # Overlay detection locations
  geom_sf(
    data = shared_sf %>% filter(count > 0),
    shape = 4,
    stroke = 2.2,
    size = 2,
    aes(color = factor(ARU))
  ) +
  scale_color_manual(values = c("black","blue"), name = "ARU")+
  coord_sf(expand = FALSE) +
  facet_grid(. ~ Atlas,
             labeller = labeller(
               Atlas = c(
                 OBBA2 = "Atlas 2",
                 OBBA3 = "Atlas 3"
               )
             )) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(size = 12, face = "bold")
  )+
  ggtitle(sp_english)

# --- Analysis of how strongly "Atlas effect" changes when incorporating survey covariates
gam0 <- gam(
  count ~ Atlas,
  data = shared_sf %>% subset(ARU == 0),
  family = "poisson"
)

gam1 <- gam(
  count ~ Atlas + s(DayOfYear) + s(Hours_Since_Sunrise),
  data = shared_sf %>% subset(ARU == 0),
  family = "poisson"
)

summary(gam0)
summary(gam1)
plot(gam1)

# Atlas 3 effect without/with accounting for timing
summary(gam0)$p.coeff[2] %>% exp()
summary(gam1)$p.coeff[2] %>% exp()
