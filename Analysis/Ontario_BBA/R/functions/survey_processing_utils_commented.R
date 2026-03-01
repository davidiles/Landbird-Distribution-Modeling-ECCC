# Create numeric survey_id for each survey

#' Construct a stable survey identifier
#'
#' Builds a deterministic ID string from key survey attributes (atlas,
#' location, date/time, survey type, and optional identifiers). This is used
#' to:
#' - link records across tables,
#' - deduplicate repeated entries, and
#' - keep joins stable across reruns.
#'
#' @param project_name Atlas/project label (e.g., "OBBA2", "OBBA3").
#' @param lat Latitude.
#' @param lon Longitude.
#' @param date_time Survey datetime (character/POSIXct).
#' @param survey_type Survey type label.
#' @param square_id Optional atlas square identifier.
#' @param point_id Optional point identifier within square.
#' @return Character vector of survey IDs.
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

#' Build an observation-level count matrix from long-format detections
#'
#' Converts detection records (one row per species observation) into a
#' matrix-like wide representation with one row per survey and one column
#' per species, plus survey metadata. This is a common intermediate format
#' for modelling and for computing per-survey summaries.
#'
#' @param det_long Long-format detection table.
#' @param species_col Column name with species identifier.
#' @param count_col Column name with counts.
#' @param survey_id_col Column name identifying surveys.
#' @param keep_survey_cols Survey-level metadata columns to keep.
#' @return A list with count matrix and metadata (see function body).
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

#' Add hours-since-sunrise (HSS) covariate to surveys
#'
#' Computes the time difference between the survey start and local sunrise
#' on the same date, producing a key detectability covariate used in the
#' modelling scripts. The function also typically standardizes datetime
#' parsing and may add auxiliary time fields.
#'
#' @param surveys_sf Survey `sf` object with date/time and coordinates.
#' @param coord_digits Digits used to round coordinates when needed (for joins/caching).
#' @return `surveys_sf` with an added `Hours_Since_Sunrise` (or similar) column.
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

#' Infer OBBA3 survey type from metadata fields
#'
#' OBBA3 records can come from ARUs, point counts, canoe routes, etc.
#' This helper parses relevant text/fields and returns a standardized
#' survey type label used downstream (e.g., "ARU" vs "Point_Count").
#'
#' @param Remarks Free-text remarks field.
#' @param Remarks2 Secondary remarks field.
#' @param EffortMeasurement1 Effort measurement field (often encodes protocol).
#' @param SurveyAreaIdentifier Field that can encode ARU deployment context.
#' @return Character vector of inferred survey types.
infer_survey_type_OBBA3 <- function(Remarks, Remarks2, EffortMeasurement1, SurveyAreaIdentifier) {
  # Returns one of: "ARU", "Point_Count", "Special", "SWIFT"
  # ARU keywords anywhere override everything else.
  # Encoding-safe: repairs invalid multibyte strings before tolower().
  
  #' Coerce to character safely
  #'
  #' Small helper to ensure fields used in string logic are character vectors.
  #'
  #' @param x Input vector.
  #' @return Character vector.
  to_chr <- function(x) {
    x <- as.character(x)
    x[is.na(x)] <- ""
    x
  }
  
  # Fix/normalize encoding (replace invalid bytes rather than erroring)

  #' Fix common text encoding issues in imported data
  #'
  #' Handles invalid UTF-8 sequences and common Windows-1252 artifacts so
  #' string matching and CSV writing behave predictably.
  #'
  #' @param x Character vector.
  #' @return Cleaned character vector.
  fix_enc <- function(x) {
    # Try to convert to UTF-8; any invalid sequences become ""
    # If you prefer to keep something, you can set sub = "?" instead of ""
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
  
  # Override by keyword scan
  out[is_aru] <- "ARU"
  
  out
}

#' Convert a count matrix into a tidy (long) counts table
#'
#' Takes the wide count matrix and returns a long table with one row per
#' (survey × species) non-zero entry (or all entries, depending on the
#' implementation). Useful for joining predictions, summarizing totals,
#' and plotting.
#'
#' @param count_matrix Wide matrix/data.frame of counts.
#' @param obs_idx_name Name of the observation index column to create/keep.
#' @return A tibble in long format.
as_counts_tbl <- function(count_matrix, obs_idx_name = "obs_idx") {
  # count_matrix is expected to have rows aligned to all_surveys_with_covs
  # Convert to tibble and add obs_idx so we can subset by obs_idx safely.
  if (inherits(count_matrix, "matrix")) count_matrix <- as.data.frame(count_matrix)
  counts <- as_tibble(count_matrix)
  counts <- counts %>% mutate(!!obs_idx_name := row_number(), .before = 1)
  counts
}

# For selecting a single survey per site

#' Deduplicate surveys by rounded coordinates
#'
#' Sometimes the same survey appears multiple times or has tiny coordinate
#' jitter. This helper rounds coordinates to `digits` and keeps one record
#' per unique location (optionally with additional rules in the body).
#'
#' @param surveys_sf Survey `sf` points.
#' @param digits Number of decimal digits for rounding.
#' @return Deduplicated `sf` object.
dedupe_by_location <- function(surveys_sf, digits = 5) {
  xy <- st_coordinates(surveys_sf)
  surveys_sf %>%
    mutate(
      x = round(xy[, 1], digits = digits),
      y = round(xy[, 2], digits = digits)
    ) %>%
    group_by(x, y) %>%
    slice_sample(n = 1) %>%
    ungroup() %>%
    select(-x, -y)
}

# Standardizing covariates

#' Standardize (z-score) covariates using grid statistics
#'
#' For comparability between Atlas 2 and 3, this helper standardizes a set
#' of covariates using reference means/SDs computed from the prediction
#' grids (grid2/grid3). This avoids differences caused purely by different
#' covariate scaling choices.
#'
#' @param surveys_sf Survey data with extracted covariates.
#' @param grid2 Prediction grid for Atlas 2.
#' @param grid3 Prediction grid for Atlas 3.
#' @param covars_for_stats Covariates used to compute mean/SD.
#' @param covars_to_standardize Covariates to z-score.
#' @return Updated `surveys_sf` with standardized covariate columns.
standardize_covars <- function(surveys_sf, grid2, grid3, covars_for_stats, covars_to_standardize) {
  present_stats <- intersect(covars_for_stats, names(surveys_sf))
  present_std   <- intersect(covars_to_standardize, names(surveys_sf))
  
  if (length(present_std) == 0) {
    return(list(surveys = surveys_sf, grid2 = grid2, grid3 = grid3, stats = NULL))
  }
  
  stats <- surveys_sf %>%
    st_drop_geometry() %>%
    select(all_of(present_stats)) %>%
    summarise(across(
      everything(),
      list(mean = \(x) mean(x, na.rm = TRUE),
           sd   = \(x) sd(x, na.rm = TRUE))
    ))
  
  for (col in present_std) {
    mu <- stats[[paste0(col, "_mean")]]
    sd <- stats[[paste0(col, "_sd")]]
    
    if (is.na(sd) || sd == 0) next
    
    surveys_sf[[col]] <- (surveys_sf[[col]] - mu) / sd
    if (col %in% names(grid2)) grid2[[col]] <- (grid2[[col]] - mu) / sd
    if (col %in% names(grid3)) grid3[[col]] <- (grid3[[col]] - mu) / sd
  }
  
  list(surveys = surveys_sf, grid2 = grid2, grid3 = grid3, stats = stats)
}
