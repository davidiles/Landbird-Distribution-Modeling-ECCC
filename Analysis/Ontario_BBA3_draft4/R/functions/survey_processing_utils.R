# ============================================================
# survey_processing_utils.R
#
# Utility functions for standardizing survey metadata, building count tables,
# thinning duplicate surveys, adding site identifiers, and creating spatial
# aggregation lookups.
#
# Contents, in the order the pipeline uses them:
#
#   1. Survey metadata        infer_survey_type_OBBA3_PC, infer_survey_type_OBBA3_CL,
#                             make_datetime_from_frac_hours, make_survey_id
#   2. Count tables           build_count_matrix, as_counts_tbl
#   3. Site identifiers       add_site_ids, summarise_site_structure
#   4. Spatial aggregation    make_hex_grid, build_pixel_polygon_index
#
# Sourced by 02-06.
# ============================================================

# ============================================================
# 1. Survey metadata
# ============================================================

# Infer standardized OBBA3 survey type labels from metadata and free-text fields.
infer_survey_type_OBBA3_PC <- function(Remarks, Remarks2, EffortMeasurement1, SurveyAreaIdentifier) {
  
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

# Modified version of infer_survey_type function for checklists
infer_survey_type_OBBA3_CL <- function(
    Survey_Type,
    Remarks,
    Remarks2,
    EffortMeasurement1,
    SurveyAreaIdentifier
) {
  
  to_chr <- function(x) {
    x <- as.character(x)
    x[is.na(x)] <- ""
    x
  }
  
  fix_enc <- function(x) {
    iconv(x, from = "", to = "UTF-8", sub = "")
  }
  
  # Existing survey type is the default
  survey_type <- fix_enc(to_chr(Survey_Type))
  
  # Fields searched for evidence of an ARU survey
  r1  <- fix_enc(to_chr(Remarks))
  r2  <- fix_enc(to_chr(Remarks2))
  em1 <- fix_enc(to_chr(EffortMeasurement1))
  sai <- fix_enc(to_chr(SurveyAreaIdentifier))
  
  txt <- tolower(
    paste(r1, r2, em1, sai, sep = " | ")
  )
  
  aru_pat <- paste0(
    "(",
    paste(
      c(
        "\\baru\\b",
        "aru_",
        "\\bsm4\\b",
        "\\bsm3\\b",
        "\\bsm2\\b",
        "sm\\s*-?\\s*4\\b",
        "song\\s*meter",
        "songmeter",
        "wildlife\\s*acoustics",
        "audiomoth",
        "autonomous\\s*(record(ing)?|recorder)",
        "\\bacoustic\\b",
        "\\baudio\\b",
        "\\brec(order|ording)?\\b",
        "\\bwav\\b",
        "\\bmp3\\b"
      ),
      collapse = "|"
    ),
    ")"
  )
  
  is_aru <- grepl(
    aru_pat,
    txt,
    perl = TRUE
  )
  
  # Retain the original Survey_Type unless ARU evidence is found
  out <- survey_type
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


# ============================================================
# 2. Count tables
# ============================================================

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

# ============================================================
# 3. Site identifiers
# ============================================================

# ==============================================================================
# add_site_ids.R
#
# Assign stable, atlas-specific site identifiers to survey locations, plus an
# optional site-by-year identifier.
#
# Rationale
#   Repeated surveys at one physical location (especially ARU stations with many
#   short-duration recordings) are pseudoreplicates: they share individuals and
#   habitat, and inflate the likelihood weight of that location. A site-level
#   iid random effect quarantines that within-site correlation, but requires a
#   site key, which the survey table does not carry.
#
#   Sites are indexed SEPARATELY PER ATLAS. A location visited in OBBA2 and
#   again in OBBA3 must NOT share a site level, or genuine between-atlas change
#   leaks into the nuisance term.
#
#   `site_year` goes one step further and separates visits to the same station
#   in different years WITHIN an atlas. Clustering is still done once per atlas
#   on the pooled coordinates, and the year is crossed in afterwards, so the
#   spatial footprint of a site is identical in every year it was visited. This
#   is deliberate: clustering within atlas-year would let the same station
#   acquire a different spatial extent (and a different label) depending on
#   which of its neighbours happened to be surveyed that year.
#
# UNITS
#   tolerance_m and snap_m are always in METRES. The CRS unit is detected from
#   sf::st_crs() and arguments are converted internally, so this works whether
#   your projected CRS is in m or km. Override with `crs_units_per_m` if
#   detection fails.
#
#   Snapping too coarsely is actively harmful: if the snap distance exceeds the
#   spacing between genuinely distinct stations, unrelated surveys merge into
#   one site level and the random effect absorbs real habitat variation, which
#   is then discarded at prediction. Keep snap_m well below station spacing.
#
# Method
#   1. Coordinates are snapped to `snap_m` and combined with the atlas label to
#      form an exact-match key, collapsing surveys to distinct locations.
#   2. If `tolerance_m > 0`, distinct locations are single-linkage clustered
#      within each atlas (dbscan::dbscan, minPts = 1; kd-tree, near-instant).
#   3. Cluster membership is mapped back to surveys by index.
#   4. If `year_col` is supplied, the survey year is crossed with the site label
#      to produce site_year_id / site_year.
#
#   Only distinct locations are ever clustered, so cost scales with the number
#   of stations, not the number of surveys.
#
#   Single-linkage is transitive: A-B at 90 m and B-C at 90 m merge A and C even
#   though they are 180 m apart. Inspect `site_radius_m` after any run with
#   tolerance_m > 0.
#
# Returns `surveys_sf` with these new columns:
#   site_id        character; carries the atlas and the cluster centroid in
#                  METRES, so it is stable across reruns, row orderings, and
#                  CRS units
#   site_atlas     integer 1..n_sites, ready for f(site_atlas, model = "iid")
#   site_radius_m  numeric; distance in metres from this survey's snapped
#                  location to its site centroid. max() within a site is that
#                  site's radius. Zero when tolerance_m = 0.
#   site_year_id   character; site_id with the survey year appended
#                  (only when year_col is not NULL)
#   site_year      integer 1..n_site_years, ready for f(site_year, model = "iid")
#                  (only when year_col is not NULL)
#
# Usage
#   surveys_f <- add_site_ids(surveys_f, tolerance_m = 100)
#   print(summarise_site_structure(surveys_f, by_type = FALSE))
#
#   # inspect chaining
#   surveys_f |>
#     sf::st_drop_geometry() |>
#     dplyr::distinct(site_id, site_radius_m) |>
#     dplyr::arrange(dplyr::desc(site_radius_m)) |>
#     head(10)
# ==============================================================================


# --- Internal: coerce a column to an integer year -----------------------------
# Accepts POSIXct/POSIXlt/Date (uses the stored tzone), a numeric year, or a
# character/factor that is either a bare year or a leading-ISO datetime string.

.resolve_survey_year <- function(x, col) {
  
  if (inherits(x, "POSIXt") || inherits(x, "Date")) {
    return(as.integer(format(x, "%Y")))
  }
  
  if (is.numeric(x)) {
    yr <- as.integer(x)
    ok <- is.na(yr) | (yr >= 1800 & yr <= 2200)
    if (!all(ok)) {
      stop("Column '", col, "' is numeric but contains values that are not ",
           "plausible years (e.g. ", x[which(!ok)[1]], ").")
    }
    return(yr)
  }
  
  if (is.character(x) || is.factor(x)) {
    xc <- as.character(x)
    yr <- suppressWarnings(as.integer(xc))
    if (anyNA(yr) && !all(is.na(xc))) {
      # not bare years; try a leading ISO date
      yr <- suppressWarnings(as.integer(substr(xc, 1, 4)))
    }
    ok <- is.na(yr) | (yr >= 1800 & yr <= 2200)
    if (!all(ok)) {
      stop("Could not parse a year from character column '", col,
           "'. Convert it to Date/POSIXct or to an integer year first.")
    }
    return(yr)
  }
  
  stop("Column '", col, "' is of class '", paste(class(x), collapse = "/"),
       "'; expected Date, POSIXct, numeric year, or character.")
}


add_site_ids <- function(surveys_sf,
                         atlas_col       = "Atlas",
                         year_col        = "Date_Time_Local",
                         tolerance_m     = 0,
                         snap_m          = 1,
                         crs_units_per_m = NULL,
                         verbose         = TRUE) {
  
  stopifnot(inherits(surveys_sf, "sf"))
  
  if (!atlas_col %in% names(surveys_sf)) {
    stop("Column '", atlas_col, "' not found in surveys_sf.")
  }
  if (!is.null(year_col) && !year_col %in% names(surveys_sf)) {
    stop("Column '", year_col, "' not found in surveys_sf. ",
         "Pass year_col = NULL to skip the site-by-year index.")
  }
  if (!all(sf::st_geometry_type(surveys_sf) == "POINT")) {
    stop("add_site_ids() expects POINT geometries.")
  }
  if (isTRUE(sf::st_is_longlat(surveys_sf))) {
    stop("surveys_sf is in geographic coordinates; project to a metric CRS first.")
  }
  if (tolerance_m > 0 && !requireNamespace("dbscan", quietly = TRUE)) {
    stop("tolerance_m > 0 requires the 'dbscan' package. install.packages('dbscan')")
  }
  if (tolerance_m > 0 && tolerance_m < snap_m) {
    stop("tolerance_m (", tolerance_m, ") is smaller than snap_m (", snap_m, ").")
  }
  
  # --- 0. Resolve CRS units ---------------------------------------------------
  # m_to_crs converts a distance in metres into CRS units:
  #   CRS in metres     -> 1
  #   CRS in kilometres -> 0.001
  
  if (is.null(crs_units_per_m)) {
    crs_unit <- sf::st_crs(surveys_sf)$ud_unit
    if (is.null(crs_unit)) {
      stop(
        "Could not determine CRS units. Pass crs_units_per_m explicitly ",
        "(1 if the CRS is in metres, 0.001 if in kilometres)."
      )
    }
    m_to_crs <- as.numeric(
      units::set_units(units::set_units(1, "m"), crs_unit, mode = "standard")
    )
  } else {
    m_to_crs <- as.numeric(crs_units_per_m)
  }
  
  if (!is.finite(m_to_crs) || m_to_crs <= 0) {
    stop("Resolved an invalid CRS unit conversion factor: ", m_to_crs)
  }
  
  snap_crs <- snap_m      * m_to_crs
  tol_crs  <- tolerance_m * m_to_crs
  
  if (verbose) {
    message(
      "add_site_ids(): CRS units = ",
      format(1 / m_to_crs), " m per unit; snapping at ", snap_m, " m (= ",
      format(snap_crs), " CRS units)"
    )
  }
  
  n_surveys <- nrow(surveys_sf)
  xy        <- sf::st_coordinates(surveys_sf)[, 1:2, drop = FALSE]
  
  if (anyNA(xy)) {
    stop("Some survey coordinates are NA; resolve these before assigning sites.")
  }
  
  dat   <- sf::st_drop_geometry(surveys_sf)
  atlas <- as.character(dat[[atlas_col]])
  if (anyNA(atlas)) stop("Some values of '", atlas_col, "' are NA.")
  
  # --- 0b. Resolve survey year ------------------------------------------------
  
  if (!is.null(year_col)) {
    survey_year <- .resolve_survey_year(dat[[year_col]], year_col)
    if (anyNA(survey_year)) {
      stop(sum(is.na(survey_year)), " surveys have a missing year in '",
           year_col, "'; resolve these before assigning site-year ids.")
    }
  } else {
    survey_year <- NULL
  }
  
  # --- 1. Collapse surveys to distinct snapped locations ----------------------
  
  xs <- round(xy[, 1] / snap_crs) * snap_crs
  ys <- round(xy[, 2] / snap_crs) * snap_crs
  
  key      <- paste(atlas, xs, ys, sep = "\r")
  key_uniq <- unique(key)
  loc_idx  <- match(key, key_uniq)
  
  first_row <- match(key_uniq, key)
  loc_x     <- xs[first_row]
  loc_y     <- ys[first_row]
  loc_atlas <- atlas[first_row]
  n_loc     <- length(key_uniq)
  
  if (verbose) {
    message(
      "add_site_ids(): ", format(n_surveys, big.mark = ","), " surveys -> ",
      format(n_loc, big.mark = ","), " distinct locations"
    )
  }
  
  # --- 2. Cluster distinct locations within each atlas ------------------------
  
  if (tolerance_m > 0) {
    
    clust  <- integer(n_loc)
    offset <- 0L
    
    for (a in unique(loc_atlas)) {
      in_a <- which(loc_atlas == a)
      cl <- dbscan::dbscan(
        cbind(loc_x[in_a], loc_y[in_a]),
        eps    = tol_crs,
        minPts = 1
      )$cluster
      clust[in_a] <- cl + offset
      offset      <- offset + max(cl)
    }
    
  } else {
    clust <- seq_len(n_loc)
  }
  
  # --- 3. Stable site label, centroid expressed in METRES ---------------------
  
  cx <- as.numeric(tapply(loc_x, clust, mean))
  cy <- as.numeric(tapply(loc_y, clust, mean))
  ca <- as.character(tapply(loc_atlas, clust, function(z) z[1]))
  
  clust_lab <- paste0(ca, "_", round(cx / m_to_crs), "_", round(cy / m_to_crs))
  
  if (anyDuplicated(clust_lab)) {
    warning(
      "Distinct site clusters share a rounded centroid label; ",
      "disambiguating with a numeric suffix."
    )
    clust_lab <- make.unique(clust_lab, sep = "_dup")
  }
  
  clust_key  <- as.integer(names(tapply(loc_x, clust, mean)))
  clust_pos  <- match(clust, clust_key)     # location -> row in cx/cy/clust_lab
  lab_by_loc <- clust_lab[clust_pos]
  site_id    <- lab_by_loc[loc_idx]
  
  # Distance from each distinct location to its site centroid, in metres
  loc_r       <- sqrt((loc_x - cx[clust_pos])^2 + (loc_y - cy[clust_pos])^2) / m_to_crs
  site_radius <- as.numeric(tapply(loc_r, clust, max))   # per cluster
  
  surveys_sf$site_id       <- site_id
  surveys_sf$site_atlas    <- as.integer(factor(site_id))
  surveys_sf$site_radius_m <- loc_r[loc_idx]
  
  # --- 3b. Cross the site with the survey year --------------------------------
  # site_year nests inside site_atlas: every site_year level belongs to exactly
  # one site_atlas level, so the two can be fitted together as a two-level
  # nuisance structure, or site_year can be used alone.
  
  if (!is.null(year_col)) {
    site_year_id <- paste0(site_id, "_", survey_year)
    surveys_sf$site_year_id <- site_year_id
    surveys_sf$site_year    <- as.integer(factor(site_year_id))
  }
  
  # --- 4. Diagnostics ---------------------------------------------------------
  
  if (verbose) {
    n_sites <- length(unique(site_id))
    per     <- tabulate(surveys_sf$site_atlas)
    
    message(
      "add_site_ids(): ", format(n_sites, big.mark = ","), " sites",
      if (tolerance_m > 0) paste0(" (tolerance_m = ", tolerance_m, ")") else "",
      "; surveys per site: median ", stats::median(per),
      ", mean ", round(mean(per), 1), ", max ", max(per)
    )
    message(
      "add_site_ids(): ", sum(per == 1), " single-visit sites, ",
      sum(per >= 10), " sites with >= 10 surveys"
    )
    
    if (!is.null(year_col)) {
      n_sy   <- length(unique(site_year_id))
      per_sy <- tabulate(surveys_sf$site_year)
      message(
        "add_site_ids(): ", format(n_sy, big.mark = ","), " site-years",
        " (", round(n_sy / n_sites, 2), " years per site on average)",
        "; surveys per site-year: median ", stats::median(per_sy),
        ", mean ", round(mean(per_sy), 1), ", max ", max(per_sy)
      )
    }
    
    if (tolerance_m > 0) {
      qs <- stats::quantile(site_radius, c(0.5, 0.95, 1), names = FALSE)
      message(
        "add_site_ids(): site radius (m): median ", round(qs[1]),
        ", 95th ", round(qs[2]), ", max ", round(qs[3]),
        "; ", sum(site_radius > 2 * tolerance_m),
        " sites exceed 2x tolerance_m"
      )
      
      if (qs[3] > 3 * tolerance_m) {
        worst <- order(site_radius, decreasing = TRUE)[seq_len(min(5, length(site_radius)))]
        warning(
          "Largest site extends ", round(qs[3]), " m from its centroid, well ",
          "beyond tolerance_m = ", tolerance_m, " m. Single-linkage chaining ",
          "has likely merged a transect of stations. Widest sites: ",
          paste0(clust_lab[worst], " (", round(site_radius[worst]), " m)",
                 collapse = "; "),
          ". Inspect the site_radius_m column and consider a smaller tolerance_m."
        )
      }
    }
  }
  
  surveys_sf
}
# ==============================================================================
# summarise_site_structure()
#
# by_type = FALSE gives the distribution the site random effect actually sees
# (site_id pools survey types). by_type = TRUE attributes replication to types.
# ==============================================================================

summarise_site_structure <- function(surveys_sf,
                                     atlas_col = "Atlas",
                                     type_col  = "Survey_Type",
                                     by_type   = TRUE) {
  
  stopifnot("site_id" %in% names(surveys_sf))
  
  d   <- sf::st_drop_geometry(surveys_sf)
  grp <- if (by_type) c(atlas_col, type_col) else atlas_col
  
  # Summary column is n_svy, NOT n_surveys: reusing the grouped column's name
  # inside summarise() overwrites the per-site vector with the group total
  # before the median/quantile/max are evaluated.
  per_site <- d %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c(grp, "site_id")))) %>%
    dplyr::summarise(n_svy = dplyr::n(), .groups = "drop")
  
  per_site %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(grp))) %>%
    dplyr::summarise(
      n_sites       = dplyr::n(),
      total_surveys = sum(n_svy),
      per_site_med  = stats::median(n_svy),
      per_site_q90  = as.numeric(stats::quantile(n_svy, 0.9)),
      per_site_max  = max(n_svy),
      pct_single    = round(100 * mean(n_svy == 1), 1),
      pct_svy_in_sites_ge10 = round(100 * sum(n_svy[n_svy >= 10]) / sum(n_svy), 1),
      .groups = "drop"
    )
}


# ============================================================
# 4. Spatial aggregation
# ============================================================

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

