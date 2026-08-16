# ============================================================
# covariate_processing_utils.R
#
# Utility functions for rasterizing vector covariates, matching sf/raster
# coordinate reference systems, and extracting raster summaries to grid cells.
# CRS transformations are intentionally handled inside helper functions; do not
# simplify or replace CRS definitions without checking downstream assumptions.
# ============================================================

# Simplify, buffer, rasterize, and optionally save an sf layer as a binary raster.
buf_rasterize_save <- function(sf_obj, buffer_dist_km, out_path, label,
                               template     = template_30m,
                               simplify_tol = SIMPLIFY_KM) {
  
  message(sprintf("[%s] %d features | simplify %.0f m | buffer %.0f m | rasterizing...",
                  label, nrow(sf_obj), simplify_tol * 1000, buffer_dist_km * 1000))
  
  v <- sf_obj %>%
    sf::st_simplify(preserveTopology = TRUE, dTolerance = simplify_tol) %>%
    sf::st_make_valid() %>%
    sf::st_buffer(dist = buffer_dist_km) %>%
    sf::st_make_valid() %>%
    terra::vect()
  
  r <- terra::rasterize(v, template, field = 1L,
                        background = 0L, touches = TRUE)
  
  if (!is.null(out_path)) {
    terra::writeRaster(r, out_path, datatype = "INT1U",
                       overwrite = TRUE, gdal = GTIFF_OPTS)
    message(sprintf("[%s] Saved -> %s", label, basename(out_path)))
  }
  
  invisible(r)
}

# Load a raster when the file exists; return NULL otherwise.
try_rast <- function(path) {
  if (!file.exists(path)) return(NULL)
  terra::rast(path)
}

# Create square polygons centred on point geometries.
make_square_buffer <- function(points, half_width = 500) {
  coords <- st_coordinates(points)
  
  squares <- lapply(seq_len(nrow(coords)), function(i) {
    x <- coords[i, "X"]
    y <- coords[i, "Y"]
    st_polygon(list(matrix(c(
      x - half_width, y - half_width,
      x + half_width, y - half_width,
      x + half_width, y + half_width,
      x - half_width, y + half_width,
      x - half_width, y - half_width   # close the polygon
    ), ncol = 2, byrow = TRUE)))
  })
  
  st_sfc(squares, crs = st_crs(points))
}

# Assign each point the ID of its containing polygon, optionally using nearest fallback.
assign_poly_id <- function(pts, polys, id_col, nearest_fallback = TRUE) {
  stopifnot(inherits(pts, "sf"), inherits(polys, "sf"))
  stopifnot(id_col %in% names(polys))
  
  polys <- st_transform(polys, st_crs(pts))
  
  hit <- st_within(pts, polys)               # list column of indices
  idx <- rep(NA_integer_, nrow(pts))
  
  has_hit <- lengths(hit) > 0L
  if (any(has_hit)) {
    idx[has_hit] <- vapply(hit[has_hit], function(x) x[1], integer(1))
  }
  
  if (nearest_fallback) {
    missing <- is.na(idx)
    if (any(missing)) {
      idx[missing] <- st_nearest_feature(pts[missing, , drop = FALSE], polys)
    }
  }
  
  polys[[id_col]][idx]
}

# Extract mean raster values within grid cells, optionally treating raster NA as zero.
extract_mean <- function(r, grid_cells, boundary_sf, na_as_zero = FALSE) {
  if (is.null(r)) return(rep(NA_real_, nrow(grid_cells)))
  
  boundary_sf <- match_sf_to_raster_crs(boundary_sf, r)
  grid_cells  <- match_sf_to_raster_crs(grid_cells, r)
  
  r <- crop_if(r, boundary_sf)
  
  if (na_as_zero) {
    r[is.na(r)] <- 0
  }
  
  exactextractr::exact_extract(r, grid_cells, "mean")
}

# Extract fractional coverage of categorical raster classes within grid cells.
extract_frac <- function(r, grid_cells, boundary_sf,
                         prefix = NULL,
                         drop0 = TRUE,
                         clean_names = FALSE) {
  if (is.null(r)) return(NULL)
  
  boundary_sf <- match_sf_to_raster_crs(boundary_sf, r)
  grid_cells  <- match_sf_to_raster_crs(grid_cells, r)
  
  r <- crop_if(r, boundary_sf)
  
  df <- suppressWarnings(exactextractr::exact_extract(r, grid_cells, "frac"))
  if (is.null(df) || nrow(df) == 0) return(NULL)
  
  if (drop0) df <- drop_frac0_cols(df)
  
  if (!is.null(prefix)) {
    names(df) <- gsub("^frac_?", paste0(prefix, "_"), names(df))
  }
  
  if (clean_names) {
    df <- clean_frac_names(df)
  }
  
  df
}

# Transform sf data to match a raster CRS after checking that both CRS values exist.
match_sf_to_raster_crs <- function(x, r) {
  r_crs <- sf::st_crs(terra::crs(r, proj = TRUE))
  
  if (is.na(r_crs)) {
    stop("Raster has missing CRS.")
  }
  if (is.na(sf::st_crs(x))) {
    stop("sf object has missing CRS.")
  }
  
  if (sf::st_crs(x) != r_crs) {
    x <- sf::st_transform(x, r_crs)
  }
  
  x
}

# Clean class-fraction column names returned by exactextractr.
clean_frac_names <- function(df) {
  names(df) <- gsub("^frac_1\\.", "", names(df))
  names(df) <- gsub("^frac1\\.", "", names(df))
  df
}

# Drop fraction columns for the zero/background class.
drop_frac0_cols <- function(df) {
  df %>% dplyr::select(-dplyr::matches("^frac_0\\.?"))
}

crop_if <- function(r, boundary_sf) {
  if (is.null(r)) return(NULL)
  v <- terra::vect(boundary_sf)
  r <- terra::crop(r, v)
  terra::mask(r, v)
}

add_solar_time_covariates_reference <- function(
    surveys_sf,
    reference_hours = -3,
    coord_digits = 6,
    adjacent_days = 1
) {
  
  # ==========================================================================
  # 1. Validate inputs
  # ==========================================================================
  
  if (!inherits(surveys_sf, "sf")) {
    stop("surveys_sf must be an sf object.")
  }
  
  required_cols <- c(
    "Date_Time",
    "Latitude",
    "Longitude"
  )
  
  missing_cols <- setdiff(
    required_cols,
    names(surveys_sf)
  )
  
  if (length(missing_cols) > 0) {
    stop(
      "surveys_sf is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  if (!inherits(surveys_sf$Date_Time, "POSIXct")) {
    stop("Date_Time must be a POSIXct vector.")
  }
  
  if (
    length(reference_hours) != 1 ||
    !is.numeric(reference_hours) ||
    !is.finite(reference_hours)
  ) {
    stop("reference_hours must be one finite numeric value.")
  }
  
  if (
    length(coord_digits) != 1 ||
    !is.numeric(coord_digits) ||
    !is.finite(coord_digits) ||
    coord_digits < 0 ||
    coord_digits != as.integer(coord_digits)
  ) {
    stop("coord_digits must be a non-negative integer.")
  }
  
  if (
    length(adjacent_days) != 1 ||
    !is.numeric(adjacent_days) ||
    !is.finite(adjacent_days) ||
    adjacent_days < 1 ||
    adjacent_days != as.integer(adjacent_days)
  ) {
    stop("adjacent_days must be a positive integer.")
  }
  
  if (
    any(!is.finite(surveys_sf$Latitude)) ||
    any(!is.finite(surveys_sf$Longitude))
  ) {
    stop(
      "Latitude and Longitude must contain finite numeric values."
    )
  }
  
  original_n <- nrow(surveys_sf)
  
  if (original_n == 0) {
    stop("surveys_sf contains no rows.")
  }
  
  # Preserve the original values, whose clock components represent local
  # survey time but are incorrectly labelled as UTC.
  original_date_time <- surveys_sf$Date_Time
  
  # ==========================================================================
  # 2. Determine the timezone of each survey
  # ==========================================================================
  
  surveys_wgs84 <- sf::st_transform(
    surveys_sf,
    crs = 4326
  )
  
  coords_wgs84 <- sf::st_coordinates(
    surveys_wgs84
  )
  
  coord_key <- paste(
    coords_wgs84[, 1],
    coords_wgs84[, 2],
    sep = "_"
  )
  
  unique_coord_idx <- !duplicated(coord_key)
  
  unique_timezone <- lutz::tz_lookup(
    surveys_wgs84[unique_coord_idx, ],
    method = "fast"
  )
  
  names(unique_timezone) <- coord_key[unique_coord_idx]
  
  timezone <- unname(
    unique_timezone[coord_key]
  )
  
  if (length(timezone) != original_n) {
    stop(
      "Timezone lookup returned an unexpected number of values."
    )
  }
  
  if (anyNA(timezone)) {
    stop(
      "Timezone lookup returned NA for ",
      sum(is.na(timezone)),
      " survey(s)."
    )
  }
  
  survey_df <- surveys_sf |>
    sf::st_drop_geometry() |>
    dplyr::mutate(
      .survey_row = dplyr::row_number(),
      
      timezone = timezone,
      
      lat_key = round(
        as.numeric(Latitude),
        digits = coord_digits
      ),
      
      lon_key = round(
        as.numeric(Longitude),
        digits = coord_digits
      )
    )
  
  # ==========================================================================
  # 3. Reinterpret the falsely labelled UTC times as local clock times
  #
  # Example:
  #
  # Stored value:
  #   2005-06-30 07:00:00 UTC
  #
  # Intended interpretation:
  #   2005-06-30 07:00:00 America/Toronto
  #
  # Correct UTC instant:
  #   2005-06-30 11:00:00 UTC
  # ==========================================================================
  
  recorded_clock_text <- format(
    survey_df$Date_Time,
    format = "%Y-%m-%d %H:%M:%S",
    tz = "UTC",
    usetz = FALSE
  )
  
  corrected_time_numeric <- rep(
    NA_real_,
    original_n
  )
  
  date_time_local_chr <- rep(
    NA_character_,
    original_n
  )
  
  local_date_chr <- rep(
    NA_character_,
    original_n
  )
  
  timezone_groups <- split(
    seq_len(original_n),
    survey_df$timezone
  )
  
  for (timezone_i in names(timezone_groups)) {
    
    rows_i <- timezone_groups[[timezone_i]]
    
    corrected_i <- as.POSIXct(
      recorded_clock_text[rows_i],
      format = "%Y-%m-%d %H:%M:%S",
      tz = timezone_i
    )
    
    corrected_time_numeric[rows_i] <- as.numeric(
      corrected_i
    )
    
    date_time_local_chr[rows_i] <- format(
      corrected_i,
      format = "%Y-%m-%d %H:%M:%S %Z",
      tz = timezone_i
    )
    
    local_date_chr[rows_i] <- format(
      corrected_i,
      format = "%Y-%m-%d",
      tz = timezone_i
    )
  }
  
  invalid_time <- is.na(
    corrected_time_numeric
  )
  
  if (any(invalid_time)) {
    
    invalid_rows <- which(
      invalid_time
    )
    
    invalid_details <- paste0(
      "row ",
      survey_df$.survey_row[invalid_rows],
      ": Date_Time=",
      ifelse(
        is.na(survey_df$Date_Time[invalid_rows]),
        "NA",
        recorded_clock_text[invalid_rows]
      ),
      ", timezone=",
      survey_df$timezone[invalid_rows],
      ", lat=",
      survey_df$Latitude[invalid_rows],
      ", lon=",
      survey_df$Longitude[invalid_rows]
    )
    
    stop(
      "Could not reinterpret ",
      length(invalid_rows),
      " Date_Time value(s) as local clock times.\n",
      paste(
        invalid_details,
        collapse = "\n"
      ),
      "\nThis usually indicates a missing Date_Time or a nonexistent ",
      "local clock time during the spring daylight-saving transition."
    )
  }
  
  survey_df$Date_Time <- as.POSIXct(
    corrected_time_numeric,
    origin = "1970-01-01",
    tz = "UTC"
  )
  
  survey_df$Date_Time_Local <- date_time_local_chr
  
  survey_df$local_date <- as.Date(
    local_date_chr
  )
  
  # ==========================================================================
  # 4. Build unique solar-event queries
  # ==========================================================================
  
  unique_queries <- survey_df |>
    dplyr::distinct(
      lat_key,
      lon_key,
      local_date,
      timezone
    ) |>
    dplyr::mutate(
      .solar_query_id = dplyr::row_number()
    )
  
  survey_df <- survey_df |>
    dplyr::left_join(
      unique_queries,
      by = c(
        "lat_key",
        "lon_key",
        "local_date",
        "timezone"
      )
    )
  
  if (anyNA(survey_df$.solar_query_id)) {
    stop(
      "Some surveys could not be matched to a solar-event query."
    )
  }
  
  # Include surrounding dates so that we can identify:
  #
  # 1. The nearest sunrise;
  # 2. The nearest sunset; and
  # 3. The most recent occurrence of the sunrise-based reference point.
  
  day_offsets <- seq.int(
    from = -adjacent_days,
    to = adjacent_days
  )
  
  solar_queries <- tidyr::expand_grid(
    unique_queries,
    day_offset = day_offsets
  ) |>
    dplyr::mutate(
      event_date = local_date + day_offset,
      
      offset_label = dplyr::case_when(
        day_offset < 0 ~ paste0(
          "minus_",
          abs(day_offset)
        ),
        
        day_offset > 0 ~ paste0(
          "plus_",
          day_offset
        ),
        
        TRUE ~ "current"
      )
    )
  
  # ==========================================================================
  # 5. Calculate sunrise and sunset
  # ==========================================================================
  
  solar_events <- solar_queries |>
    dplyr::group_split(
      timezone
    ) |>
    purrr::map_dfr(
      function(query_group) {
        
        timezone_i <- query_group$timezone[[1]]
        
        sunlight_input <- query_group |>
          dplyr::transmute(
            date = event_date,
            lat = lat_key,
            lon = lon_key
          )
        
        sunlight <- suncalc::getSunlightTimes(
          data = sunlight_input,
          keep = c(
            "sunrise",
            "sunset"
          ),
          tz = timezone_i
        )
        
        if (nrow(sunlight) != nrow(query_group)) {
          stop(
            "getSunlightTimes() returned ",
            nrow(sunlight),
            " rows for timezone ",
            timezone_i,
            ", but ",
            nrow(query_group),
            " rows were expected."
          )
        }
        
        query_group |>
          dplyr::transmute(
            .solar_query_id,
            offset_label,
            sunrise_time = sunlight$sunrise,
            sunset_time = sunlight$sunset
          )
      }
    )
  
  solar_events_wide <- solar_events |>
    tidyr::pivot_wider(
      id_cols = .solar_query_id,
      names_from = offset_label,
      values_from = c(
        sunrise_time,
        sunset_time
      ),
      names_glue = "{.value}_{offset_label}"
    )
  
  survey_solar <- survey_df |>
    dplyr::left_join(
      solar_events_wide,
      by = ".solar_query_id"
    ) |>
    dplyr::arrange(
      .survey_row
    )
  
  if (nrow(survey_solar) != original_n) {
    stop(
      "Row count changed during the solar-event join."
    )
  }
  
  if (
    !identical(
      survey_solar$.survey_row,
      seq_len(original_n)
    )
  ) {
    stop(
      "Survey row order changed during solar-event processing."
    )
  }
  
  # ==========================================================================
  # 6. Helper functions
  # ==========================================================================
  
  event_columns_to_matrix <- function(
    data,
    columns
  ) {
    
    if (length(columns) == 0) {
      stop("No solar-event columns were found.")
    }
    
    output <- do.call(
      cbind,
      lapply(
        data[columns],
        as.numeric
      )
    )
    
    if (is.null(dim(output))) {
      output <- matrix(
        output,
        ncol = 1
      )
    }
    
    output
  }
  
  select_nearest_event <- function(
    event_matrix,
    survey_time
  ) {
    
    survey_numeric <- as.numeric(
      survey_time
    )
    
    difference_hours <- sweep(
      event_matrix,
      MARGIN = 1,
      STATS = survey_numeric,
      FUN = function(event_time, survey_time) {
        (survey_time - event_time) / 3600
      }
    )
    
    absolute_difference <- abs(
      difference_hours
    )
    
    absolute_difference[
      !is.finite(absolute_difference)
    ] <- Inf
    
    has_event <- rowSums(
      is.finite(event_matrix)
    ) > 0
    
    nearest_index <- max.col(
      -absolute_difference,
      ties.method = "first"
    )
    
    row_index <- seq_len(
      nrow(event_matrix)
    )
    
    relative_hours <- difference_hours[
      cbind(
        row_index,
        nearest_index
      )
    ]
    
    relative_hours[!has_event] <- NA_real_
    
    relative_hours
  }
  
  select_most_recent_event <- function(
    event_matrix,
    survey_time
  ) {
    
    survey_numeric <- as.numeric(
      survey_time
    )
    
    is_past <- sweep(
      event_matrix,
      MARGIN = 1,
      STATS = survey_numeric,
      FUN = "<="
    )
    
    is_past[is.na(is_past)] <- FALSE
    
    past_events <- event_matrix
    past_events[!is_past] <- -Inf
    
    output <- apply(
      past_events,
      MARGIN = 1,
      FUN = max
    )
    
    output[!is.finite(output)] <- NA_real_
    
    output
  }
  
  # ==========================================================================
  # 7. Calculate relative time to the nearest sunrise and sunset
  # ==========================================================================
  
  sunrise_columns <- grep(
    pattern = "^sunrise_time_",
    x = names(survey_solar),
    value = TRUE
  )
  
  sunset_columns <- grep(
    pattern = "^sunset_time_",
    x = names(survey_solar),
    value = TRUE
  )
  
  sunrise_matrix <- event_columns_to_matrix(
    survey_solar,
    sunrise_columns
  )
  
  sunset_matrix <- event_columns_to_matrix(
    survey_solar,
    sunset_columns
  )
  
  sunrise_relative_hours <- select_nearest_event(
    event_matrix = sunrise_matrix,
    survey_time = survey_solar$Date_Time
  )
  
  sunset_relative_hours <- select_nearest_event(
    event_matrix = sunset_matrix,
    survey_time = survey_solar$Date_Time
  )
  
  # ==========================================================================
  # 8. Calculate hours after the most recent reference point
  #
  # The reference point is defined relative to each sunrise:
  #
  #   reference_time = sunrise_time + reference_hours
  #
  # For example, reference_hours = -3 means that the reference point occurs
  # three hours before sunrise.
  #
  # We select the most recent reference point rather than the nearest one.
  # This creates a continuous daily time axis:
  #
  #   reference point             =  0 hours
  #   sunrise                     =  3 hours, when reference_hours = -3
  #   morning survey              =  approximately 6 hours
  #   evening survey              =  approximately 14 hours
  #   midnight survey             =  approximately 20 hours
  #
  # Values can differ from exactly 0-24 around daylight-saving transitions
  # because the elapsed duration between successive local reference points
  # can be 23 or 25 hours.
  # ==========================================================================
  
  reference_matrix <- sunrise_matrix +
    reference_hours * 3600
  
  most_recent_reference_numeric <- select_most_recent_event(
    event_matrix = reference_matrix,
    survey_time = survey_solar$Date_Time
  )
  
  hours_after_reference <- (
    as.numeric(survey_solar$Date_Time) -
      most_recent_reference_numeric
  ) / 3600
  
  hours_after_reference[
    !is.finite(hours_after_reference)
  ] <- NA_real_
  
  # ==========================================================================
  # 9. Add covariates back to the original sf object
  # ==========================================================================
  
  surveys_sf$Date_Time_Original <-
    original_date_time
  
  # Corrected universal instant, represented in UTC.
  surveys_sf$Date_Time <-
    survey_solar$Date_Time
  
  # Readable local clock time with local timezone abbreviation.
  surveys_sf$Date_Time_Local <-
    survey_solar$Date_Time_Local
  
  surveys_sf$timezone <-
    survey_solar$timezone
  
  surveys_sf$local_date <-
    survey_solar$local_date
  
  surveys_sf$Sunrise_Relative_Hours <-
    sunrise_relative_hours
  
  surveys_sf$Sunset_Relative_Hours <-
    sunset_relative_hours
  
  surveys_sf$Hours_After_Reference <-
    hours_after_reference
  
  # Store the selected reference value so that the resulting object remains
  # self-documenting.
  surveys_sf$Reference_Hours_Since_Sunrise <-
    reference_hours
  
  # ==========================================================================
  # 10. Final consistency checks
  # ==========================================================================
  
  if (nrow(surveys_sf) != original_n) {
    stop(
      "Final row count differs from the input row count."
    )
  }
  
  if (
    any(
      surveys_sf$Hours_After_Reference < -1e-8,
      na.rm = TRUE
    )
  ) {
    stop(
      "At least one survey occurs before its selected most recent ",
      "reference point. This indicates an error in reference selection."
    )
  }
  
  surveys_sf
}
