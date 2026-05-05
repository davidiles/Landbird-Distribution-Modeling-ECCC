rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(sf)
  library(terra)
  library(ggplot2)
  library(INLA)
  library(inlabru)
  library(fmesher)
  library(here)
})

# ------------------------------------------------------------
# Centralized paths
# ------------------------------------------------------------

source(here::here("R", "00_config_paths.R"))
source(file.path(paths$functions, "inla_model_utils2.R"))

# ------------------------------------------------------------
# Load data
# ------------------------------------------------------------

in_file <- file.path(paths$data_clean, "birds", "data_ready_for_analysis.rds")
dat <- readRDS(in_file)

all_surveys      <- dat$all_surveys
counts           <- dat$counts
grid_OBBA2       <- dat$grid_OBBA2
grid_OBBA3       <- dat$grid_OBBA3
study_boundary   <- dat$study_boundary %>% st_as_sf()
species_to_model <- dat$species_to_model
hex_grid         <- dat$hex_grid
safe_dates_bcr   <- dat$safe_dates_bcr

water_path <- file.path(
  paths$data,
  "Spatial",
  "Ontario_Hydro_Network_(OHN)_-_Waterbody",
  "Ontario_Hydro_Network_(OHN)_-_Waterbody.shp"
)

# ------------------------------------------------------------
# Helper: make AOI buffer
# ------------------------------------------------------------

make_latlon_buffer <- function(lon,
                               lat,
                               buffer_m = 5000,
                               crs_lonlat = 4326,
                               crs_projected = 3978) {
  
  pt <- st_as_sf(
    data.frame(lon = lon, lat = lat),
    coords = c("lon", "lat"),
    crs = crs_lonlat
  )
  
  pt %>%
    st_transform(crs_projected) %>%
    st_buffer(dist = buffer_m) %>%
    st_make_valid()
}

# ------------------------------------------------------------
# Choose species / AOI
# ------------------------------------------------------------

coords <- c(51.553226, -80.986479)
lat <- coords[1]
lon <- coords[2]

buffer_m <- 50000
covariate_to_plot <- "on_river"

i <- which(species_to_model$english_name == "Black-and-white Warbler")

sp_name <- species_to_model$english_name[i]
sp_code <- as.character(species_to_model$species_id[i])
sp_file <- sp_filename(sp_name)

sp_dat <- all_surveys %>%
  mutate(count = counts[[sp_code]])

# ------------------------------------------------------------
# Filter data to species-safe dates
# ------------------------------------------------------------

sp_safe_dates <- safe_dates_bcr %>%
  filter(sp_english == sp_name)

if (nrow(sp_safe_dates) == 0) {
  
  pred_doy <- NA_real_
  warning(
    paste0("Species '", sp_name, "' has no BCR safe dates at all."),
    call. = FALSE
  )
  
} else {
  
  fallback_start <- max(sp_safe_dates$start_doy, na.rm = TRUE)
  fallback_end   <- min(sp_safe_dates$end_doy,   na.rm = TRUE)
  
  if (fallback_start <= fallback_end) {
    
    all_bcrs <- sp_dat %>% distinct(BCR)
    safe_bcrs <- sp_safe_dates %>% distinct(BCR)
    missing_bcrs <- all_bcrs %>% anti_join(safe_bcrs, by = "BCR")
    
    if (nrow(missing_bcrs) > 0) {
      sp_safe_dates <- bind_rows(
        sp_safe_dates,
        missing_bcrs %>%
          mutate(
            sp_english = sp_name,
            start_doy = fallback_start,
            end_doy   = fallback_end,
            midpoint  = floor((fallback_start + fallback_end) / 2)
          ) %>%
          select(sp_english, BCR, start_doy, end_doy, midpoint)
      )
    }
    
    pred_doy <- floor((fallback_start + fallback_end) / 2)
    
  } else {
    pred_doy <- NA_real_
  }
}

sp_dat <- sp_dat %>%
  left_join(sp_safe_dates, by = "BCR") %>%
  filter(
    !is.na(start_doy),
    !is.na(end_doy),
    DayOfYear >= start_doy,
    DayOfYear <= end_doy
  ) %>%
  mutate(days_midpoint = DayOfYear - pred_doy)

# ------------------------------------------------------------
# Prepare AOI and spatial objects
# ------------------------------------------------------------

target_crs <- st_crs(grid_OBBA3)

aoi <- make_latlon_buffer(
  lon = lon,
  lat = lat,
  buffer_m = buffer_m,
  crs_projected = target_crs
)

grid_OBBA3_tr <- st_transform(grid_OBBA3, target_crs)
sp_dat_tr     <- st_transform(sp_dat, target_crs)

grid_sub <- st_filter(grid_OBBA3_tr, aoi, .predicate = st_intersects)
sp_sub   <- st_filter(sp_dat_tr,     aoi, .predicate = st_intersects)

if (nrow(grid_sub) == 0) {
  stop("No grid cells intersect the AOI. Check lon/lat, buffer_m, or grid CRS.")
}

if (!covariate_to_plot %in% names(grid_sub)) {
  stop("covariate_to_plot is not a column in grid_OBBA3.")
}

# ------------------------------------------------------------
# Read only water polygons near AOI
# ------------------------------------------------------------

water_prj <- sub("\\.shp$", ".prj", water_path)

water_crs <- st_crs(
  paste(readLines(water_prj), collapse = " ")
)

aoi_water_crs <- aoi %>%
  st_transform(water_crs) %>%
  st_make_valid()

water_sub <- st_read(
  water_path,
  quiet = TRUE,
  wkt_filter = st_as_text(st_geometry(aoi_water_crs))
) %>%
  st_make_valid() %>%
  st_transform(target_crs) %>%
  st_filter(aoi, .predicate = st_intersects)

water_river <- water_sub %>%
  filter(WATERBODY_ == "River")

# ------------------------------------------------------------
# Convert filtered grid centroids + covariate to raster
# ------------------------------------------------------------

grid_xy <- st_coordinates(grid_sub)

grid_df <- grid_sub %>%
  st_drop_geometry() %>%
  mutate(
    x = grid_xy[, 1],
    y = grid_xy[, 2]
  )

x_unique <- sort(unique(grid_df$x))
y_unique <- sort(unique(grid_df$y))

dx <- min(diff(x_unique))
dy <- min(diff(y_unique))

r_cov <- terra::rast(
  xmin = min(grid_df$x) - dx / 2,
  xmax = max(grid_df$x) + dx / 2,
  ymin = min(grid_df$y) - dy / 2,
  ymax = max(grid_df$y) + dy / 2,
  resolution = c(dx, dy),
  crs = target_crs$wkt
)

v_cov <- terra::vect(
  grid_df[, c("x", "y", covariate_to_plot)],
  geom = c("x", "y"),
  crs = target_crs$wkt
)

r_cov <- terra::rasterize(
  v_cov,
  r_cov,
  field = covariate_to_plot
)

cov_df <- as.data.frame(r_cov, xy = TRUE, na.rm = FALSE)
names(cov_df)[3] <- "value"

# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------

ggplot() +
  geom_raster(
    data = cov_df,
    aes(x = x, y = y, fill = value)
  ) +
  geom_sf(
    data = water_river,
    fill = "black",
    color = "black",
    linewidth = 2,
    alpha = 0.5
  ) +
  geom_sf(
    data = aoi,
    fill = NA,
    color = "black",
    linewidth = 0.8
  ) +
  geom_sf(
    data = sp_sub,
    color = "black",
    size = 3,
    alpha = 0.5
  ) +
  geom_sf(
    data = sp_sub %>% filter(count > 0),
    color = "yellow",
    size = 2
  ) +
  coord_sf(
    xlim = st_bbox(aoi)[c("xmin", "xmax")],
    ylim = st_bbox(aoi)[c("ymin", "ymax")],
    expand = FALSE
  ) +
  scale_fill_gradientn(
    na.value = "gray90",
    colours = c("white", "dodgerblue")
  ) +
  theme_bw() +
  labs(fill = covariate_to_plot)

# ------------------------------------------------------------
# Optional summary
# ------------------------------------------------------------

sp_dat %>%
  as.data.frame() %>%
  mutate(presence = as.numeric(count > 0)) %>%
  group_by(presence, BCR) %>%
  summarize(
    water_open = mean(water_open, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  arrange(BCR, presence)