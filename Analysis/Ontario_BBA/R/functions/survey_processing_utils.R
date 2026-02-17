# Create numeric survey_id for each survey
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

# Build a matrix (n_survey x n_species) of counts of individuals from each species for each survey
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

# Calculate hours since sunrise for each survey, based on date_time and geographic location
add_hours_since_sunrise <- function(surveys_sf, coord_digits = 6) {
  stopifnot(inherits(surveys_sf, "sf"))
  stopifnot("Date_Time" %in% names(surveys_sf))
  stopifnot(all(c("Latitude", "Longitude") %in% names(surveys_sf)))
  
  # Work in WGS84 for tz lookup / suncalc
  wgs84 <- sf::st_transform(surveys_sf, 4326)
  tz <- lutz::tz_lookup(wgs84, method = "fast")
  
  df <- surveys_sf %>%
    sf::st_drop_geometry() %>%
    mutate(
      timezone = tz,
      date = as.Date(Date_Time),
      # Round coords to stabilize joins (optional but recommended)
      lat_key = round(as.numeric(Latitude), coord_digits),
      lon_key = round(as.numeric(Longitude), coord_digits)
    )
  
  # Build unique query table (stable keys)
  unique_q <- df %>%
    distinct(lat_key, lon_key, date, timezone)
  
  # Compute sunrise for each unique combination
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
  
  # Enforce uniqueness to prevent join fan-out
  sunrise <- sunrise %>%
    distinct(lat_key, lon_key, date, timezone, .keep_all = TRUE)
  
  # Join back and compute hours since sunrise
  df2 <- df %>%
    left_join(sunrise, by = c("lat_key", "lon_key", "date", "timezone")) %>%
    mutate(
      Hours_Since_Sunrise = as.numeric(difftime(
        lubridate::force_tz(Date_Time, timezone),
        lubridate::force_tz(sunrise, timezone),
        units = "hours"
      ))
    )
  
  # Assert row count unchanged
  if (nrow(df2) != nrow(df)) {
    stop("Row count changed during sunrise join; sunrise keys are not unique.")
  }
  
  surveys_sf$Hours_Since_Sunrise <- df2$Hours_Since_Sunrise
  surveys_sf$timezone <- df2$timezone
  surveys_sf
}

