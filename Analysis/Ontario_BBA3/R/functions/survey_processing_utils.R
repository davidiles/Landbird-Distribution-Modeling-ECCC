# ============================================================
# survey_processing_utils.R
#
# Utility functions for standardizing survey metadata, building count tables,
# adding time covariates, thinning duplicate surveys, and creating spatial
# aggregation lookups.
# ============================================================

# Infer standardized OBBA3 survey type labels from metadata and free-text fields.
infer_survey_type_OBBA3 <- function(Remarks, Remarks2, EffortMeasurement1, SurveyAreaIdentifier) {
  
  to_chr <- function(x) {
    x <- as.character(x)
    x[is.na(x)] <- ""
    x
  }
  
  fix_enc <- function(x) {
    iconv(x, from = "", to = "UTF-8", sub = "")
  }
  
  r1  <- fix_enc(to_chr(Remarks))
  r2  <- fix_enc(to_chr(Remarks2))
  em1 <- fix_enc(to_chr(EffortMeasurement1))
  sai <- fix_enc(to_chr(SurveyAreaIdentifier))
  
  txt <- tolower(paste(r1, r2, em1, sai, sep = " | "))
  
  aru_pat <- paste0(
    "(",
    paste(c(
      "\\baru\\b",
      "aru_",
      "\\bsm4\\b", "\\bsm3\\b", "\\bsm2\\b",
      "sm\\s*-?\\s*4\\b",
      "song\\s*meter", "songmeter",
      "wildlife\\s*acoustics",
      "audiomoth",
      "autonomous\\s*(record(ing)?|recorder)",
      "\\bacoustic\\b",
      "\\baudio\\b",
      "\\brec(order|ording)?\\b",
      "\\bwav\\b", "\\bmp3\\b"
    ), collapse = "|"),
    ")"
  )
  
  is_aru <- grepl(aru_pat, txt, perl = TRUE)
  
  em1_trim <- trimws(em1)
  out <- rep("Point_Count", length(em1_trim))
  
  out[em1_trim %in% c("Special")]                <- "Special"
  out[em1_trim %in% c("SWIFT")]                  <- "SWIFT"
  out[em1_trim %in% c("IN_PERSON")]              <- "Point_Count"
  out[em1_trim %in% c("ARU_SM4", "ARU_UNKNOWN")] <- "ARU"
  
  out[is_aru] <- "ARU"
  
  out
}

# Convert a Date plus fractional hour value into a POSIXct datetime.
make_datetime_from_frac_hours <- function(date_ymd, frac_hours, tz = "UTC") {
  stopifnot(inherits(date_ymd, "Date"))
  
  secs <- round(frac_hours * 3600)
  as.POSIXct(date_ymd, tz = tz) + secs
}

# Create a deterministic survey ID from project, location, datetime, and survey type.
make_survey_id <- function(project_name, lat, lon, date_time, survey_type,
                           digits = 6) {
  
  lat_r <- round(as.numeric(lat), digits)
  lon_r <- round(as.numeric(lon), digits)
  
  paste(
    project_name,
    lat_r,
    lon_r,
    date_time,
    survey_type,
    sep = "_"
  )
}

# Convert long-format detections into a survey-by-species count matrix.
build_count_matrix <- function(
    long_data,
    survey_id_col,
    species_col,
    count_col
) {
  sp_ids <- sort(unique(long_data[[species_col]]))
  survey_ids <- sort(unique(long_data[[survey_id_col]]))
  
  mat <- matrix(
    0,
    nrow = length(survey_ids),
    ncol = length(sp_ids),
    dimnames = list(survey_ids, sp_ids)
  )
  
  summary <- long_data %>%
    group_by(
      .data[[survey_id_col]],
      .data[[species_col]]
    ) %>%
    summarise(
      total_count = sum(.data[[count_col]]),
      .groups = "drop"
    )
  
  row_idx <- match(summary[[survey_id_col]], rownames(mat))
  col_idx <- match(summary[[species_col]], colnames(mat))
  valid <- !is.na(row_idx) & !is.na(col_idx)
  
  mat[cbind(row_idx[valid], col_idx[valid])] <- summary$total_count[valid]
  
  mat[, colSums(mat) > 0, drop = FALSE]
}

# Convert a count matrix to a tibble with a stable observation index.
as_counts_tbl <- function(count_matrix, obs_idx_name = "obs_idx") {
  if (inherits(count_matrix, "matrix")) count_matrix <- as.data.frame(count_matrix)
  counts <- as_tibble(count_matrix)
  counts <- counts %>% mutate(!!obs_idx_name := row_number(), .before = 1)
  counts
}

# Add timezone and hours-since-sunrise fields to survey records.
add_hours_since_sunrise <- function(surveys_sf, coord_digits = 6) {
  stopifnot(inherits(surveys_sf, "sf"))
  stopifnot("Date_Time" %in% names(surveys_sf))
  stopifnot(all(c("Latitude", "Longitude") %in% names(surveys_sf)))
  
  wgs84 <- sf::st_transform(surveys_sf, 4326)
  tz <- lutz::tz_lookup(wgs84, method = "fast")
  
  df <- surveys_sf %>%
    sf::st_drop_geometry() %>%
    mutate(
      timezone = tz,
      date = as.Date(Date_Time),
      lat_key = round(as.numeric(Latitude), coord_digits),
      lon_key = round(as.numeric(Longitude), coord_digits)
    )
  
  unique_q <- df %>%
    distinct(lat_key, lon_key, date, timezone)
  
  sunrise <- purrr::pmap_dfr(unique_q, function(lat_key, lon_key, date, timezone) {
    suncalc::getSunlightTimes(
      date = date,
      lat  = lat_key,
      lon  = lon_key,
      keep = "sunrise",
      tz   = timezone
    ) %>%
      transmute(
        lat_key = lat_key,
        lon_key = lon_key,
        date = date,
        timezone = timezone,
        sunrise = sunrise
      )
  })
  
  sunrise <- sunrise %>%
    distinct(lat_key, lon_key, date, timezone, .keep_all = TRUE)
  
  df2 <- df %>%
    left_join(sunrise, by = c("lat_key", "lon_key", "date", "timezone")) %>%
    mutate(
      Hours_Since_Sunrise = as.numeric(difftime(
        lubridate::force_tz(Date_Time, timezone),
        lubridate::force_tz(sunrise, timezone),
        units = "hours"
      ))
    )
  
  if (nrow(df2) != nrow(df)) {
    stop("Row count changed during sunrise join; sunrise keys are not unique.")
  }
  
  surveys_sf$Hours_Since_Sunrise <- df2$Hours_Since_Sunrise
  surveys_sf$timezone <- df2$timezone
  surveys_sf
}

# Select one survey per rounded location-year, optionally closest to a target day-of-year.
dedupe_by_location <- function(surveys_sf, digits = 5, target_doy = NULL, date_col = "date") {
  
  xy <- st_coordinates(surveys_sf)
  
  surveys_sf <- surveys_sf %>%
    mutate(
      .x     = round(xy[, 1], digits = digits),
      .y     = round(xy[, 2], digits = digits),
      .year  = lubridate::year(.data[[date_col]]),
      .doy   = lubridate::yday(.data[[date_col]])
    )
  
  if (!is.null(target_doy)) {
    surveys_sf <- surveys_sf %>%
      mutate(.doy_dist = abs(.doy - target_doy))
    
    surveys_sf <- surveys_sf %>%
      group_by(.x, .y, .year) %>%
      slice_min(order_by = .doy_dist, n = 1, with_ties = FALSE) %>%
      ungroup()
  } else {
    surveys_sf <- surveys_sf %>%
      group_by(.x, .y, .year) %>%
      slice_sample(n = 1) %>%
      ungroup()
  }
  
  surveys_sf %>%
    select(-starts_with("."))
}

# Create a hexagon grid over the study boundary and keep overlapping cells.
make_hex_grid <- function(study_boundary, width_km = 25, seed = 123) {
  
  set.seed(seed)
  
  hex <- sf::st_make_grid(
    study_boundary,
    cellsize = width_km ,
    square = FALSE
  )
  
  hex_sf <- sf::st_sf(
    hex_id = seq_along(hex),
    geometry = hex
  )
  
  hex_sf <- hex_sf[sf::st_intersects(hex_sf, study_boundary, sparse = FALSE), , drop = FALSE]
  
  hex_sf
}

# Build a reusable pixel-to-polygon lookup for posterior aggregation.
build_pixel_polygon_index <- function(grid_sf,
                                      polygons_sf,
                                      poly_id_col,
                                      join = c("within", "intersects")) {
  join <- match.arg(join)
  
  stopifnot(inherits(grid_sf, "sf"))
  stopifnot(inherits(polygons_sf, "sf"))
  stopifnot(poly_id_col %in% names(polygons_sf))
  
  if (sf::st_crs(grid_sf) != sf::st_crs(polygons_sf)) {
    polygons_sf <- sf::st_transform(polygons_sf, sf::st_crs(grid_sf))
  }
  
  polygons_sf <- polygons_sf[sf::st_intersects(polygons_sf, sf::st_union(sf::st_geometry(grid_sf)),
                                               sparse = FALSE), , drop = FALSE]
  
  pix_sf <- sf::st_sf(pixel_id = seq_len(nrow(grid_sf)),
                      geometry = sf::st_geometry(grid_sf))
  
  join_fun <- switch(
    join,
    within     = sf::st_within,
    intersects = sf::st_intersects
  )
  
  pix_joined <- sf::st_join(
    pix_sf,
    polygons_sf[, poly_id_col, drop = FALSE],
    join = join_fun,
    left = TRUE
  )
  
  pix_poly_id <- pix_joined[[poly_id_col]]
  poly_ids <- sort(unique(pix_poly_id[!is.na(pix_poly_id)]))
  
  list(
    pix_poly_id = pix_poly_id,
    poly_ids = poly_ids,
    polygons_sf = polygons_sf,
    crs = sf::st_crs(grid_sf)
  )
}
