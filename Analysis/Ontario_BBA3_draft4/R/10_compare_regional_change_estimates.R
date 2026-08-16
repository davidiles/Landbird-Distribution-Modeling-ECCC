# ============================================================
# 10_compare_regional_change_estimates.R
#
# Purpose
#   Load the regional (province + BCR) population-change estimates for TWO
#   sources and compare them, region by region and species by species.
#
#   A "source" is one of three kinds, chosen by the model_A / model_B strings:
#     * an atlas model  -> the full-model tidy table written by script 08,
#                          change_estimates_<model>/regional_change_estimates_<model>.rds
#                          Specify with the bare model name, e.g. "PC_ARU_CL".
#     * "paired"        -> the repeated-survey (shared-footprint) table written by
#                          script 08, change_estimates/paired_change_estimates.rds.
#                          The paired analysis (07b) is run SEPARATELY from any
#                          PC_ARU / PC_ARU_CL model and is not a variant of one, so
#                          it is model-independent and referenced simply as
#                          "paired" -- never as "paired:<model>".
#     * "BBS"           -> an external Breeding Bird Survey trends CSV.
#
#   Use cases:
#     model vs model   (PC_ARU_CL_nosite vs PC_ARU_CL)
#     model vs paired  (PC_ARU vs paired)
#     paired vs BBS    (paired vs BBS)
#     model vs BBS     (PC_ARU vs BBS)
#
# What the paired source is / how it differs
#   The paired analysis estimates change only from stations surveyed in BOTH
#   atlases (07b's shared_change_summary), so it is:
#     * per-BCR only (no province-wide row) -> when compared against a full model
#       the join keeps the shared BCRs; against BBS it keeps BCR 12 and 13;
#     * on a fixed 90% interval (q05/q95), so paired comparisons run at ci 0.90;
#     * total percent-change only (no annual column) unless you annualise it via
#       paired_interval_years below.
#
# The region-matching problem
#   BBS strata are not the atlas BCRs. Only three are directly comparable:
#     BBS "ON"        (region_type prov_state) <-> atlas province-wide
#     BBS "CA-ON-12"  (region_type stratum)    <-> atlas BCR 12
#     BBS "CA-ON-13"  (region_type stratum)    <-> atlas BCR 13
#   Each loader maps its regions onto a canonical region_key
#   (province / BCR12 / BCR13 / BCR74 / BCR76 / BCR77). The comparison is an
#   INNER join on region_key, so:
#     model vs model -> all six shared regions are compared;
#     model vs BBS   -> only the three comparable regions survive the join.
#
# Metrics
#   Two commensurate quantities are available from both sides:
#     "percent_change" : total % change over the interval
#                        (atlas: 100*prop_change ; BBS: percent_change)
#     "annual_trend"   : geometric %/year
#                        (atlas: annual_trend_pct_* ; BBS: trend)
#   Pick with compare_metric below. annual_trend is only available if script 08
#   was run with atlas reference years set.
#
# Credible intervals
#   The atlas file's qlow/qhigh sit at its stored ci_level (0.90 by default:
#   5th/95th percentiles). BBS quantile columns are chosen to match. Supported
#   levels: 0.90 and 0.95. Direction (increase / decrease / uncertain) is derived
#   uniformly from whether that interval excludes zero.
#
# Caveats worth remembering
#   * The BBS "ON" estimate pools its own strata (here CA-ON-8/12/13), which are
#     not identical to the atlas study-area BCRs; treat province<->ON as the
#     intended-but-imperfect whole-Ontario comparison.
#   * The BBS interval (e.g. 2002-2022) and the atlas interval (~2003-2023) are
#     close but not identical; total percent_change is more sensitive to that
#     offset than the annualised trend.
#   * Species are matched on common name (atlas sp_english <-> BBS english),
#     normalised for case/whitespace. Unmatched species are reported and saved;
#     supply species_crosswalk_csv to reconcile naming differences.
#
# Outputs (in change_comparisons/)
#   comparison_<A>_vs_<B>.csv/.rds          per species x region, side by side
#   comparison_summary_<A>_vs_<B>.csv       per-region + overall agreement stats
#   comparison_unmatched_<A>_vs_<B>.csv     species present in one source only
#   comparison_<A>_vs_<B>.pdf               scatter with 1:1 line (if make_plot)
# ============================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(ggplot2)
})

source(here::here("R", "00_config_paths.R"))

# ============================================================
# CONFIG
# ============================================================

# Each of model_A / model_B is one of:
#   "<model>"   full atlas model      (e.g. "PC_ARU_CL")
#   "paired"    the (model-independent) repeated-survey analysis from 07b/08
#   "BBS"       external BBS trends CSV
model_A <- "paired"
model_B <- "BBS"             # try "BBS" or another atlas model

# Quantity to compare.
#   "percent_change" -> total % change over the interval (always available)
#   "annual_trend"   -> geometric %/year (needs annual columns on both sides;
#                       paired needs paired_interval_years set)
compare_metric <- "percent_change"

# Credible-interval level. NULL infers it: from the atlas source(s), or forced
# to 0.90 whenever a paired source is involved (paired carries only q05/q95).
# Falls back to 0.90 otherwise. Only 0.90 and 0.95 are supported for BBS.
ci_level <- NULL

# Only used to annualise the paired source (which is stored as total % change).
# Must equal the atlas reference interval used in script 08 so annual trends are
# comparable. NULL leaves paired as percent-change only.
paired_interval_years <- NULL

# BBS trends CSV (only read when a source is "BBS"). Auto-detects comma or tab.
bbs_csv_path <- file.path(paths$data, "bbs_trend_estimates",
                          "Ontario_Atlas_trends_means_by_atlas_period_2002_2022.csv")

# Optional: restrict BBS to a single interval when the file holds several.
# c(start = 2002, end = 2022) or NULL to use every row (duplicates are de-duped
# with a warning).
bbs_years <- NULL

# Optional species-name crosswalk applied to source B before joining, for
# reconciling atlas vs BBS common-name differences. CSV with columns `from`,`to`
# (from = name as it appears in source B; to = the matching source-A name).
# NULL to skip.
species_crosswalk_csv <- NULL

# ---- Detection-based data-sufficiency filter -------------------------------
# Drop species x region combinations that were poorly sampled in the atlas, using
# the per-region squares-with-detection counts (n_sq_det_OBBA2 / n_sq_det_OBBA3)
# that script 08 writes to the FULL-model export. These counts come from the
# survey data, so any model's export gives authoritative counts for the province
# and every BCR; a paired or BBS comparison borrows them from a chosen model.
#
#   min_n_sq_det : squares-with-detection threshold. NULL disables the filter
#                  (counts are still attached to the output when available).
#   detection_rule : "either" keeps a row if the threshold is met in EITHER
#                    atlas; "both" requires it in both.
#   detection_counts_model : atlas model whose full export supplies the counts.
#                    NULL auto-resolves to an atlas source. When NEITHER source is
#                    an atlas model (e.g. paired vs BBS), set this explicitly.
min_n_sq_det           <- 20
detection_rule         <- "either"
detection_counts_model <- "PC_ARU_nosite"

# Draw (and save) the comparison scatter (needs ggplot2).
make_plot <- TRUE

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------

est_root <- paths$model_output
out_dir  <- file.path(paths$model_output, "change_comparisons")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(compare_metric %in% c("percent_change", "annual_trend"))
stopifnot(detection_rule %in% c("either", "both"))
stopifnot(is.null(min_n_sq_det) ||
            (is.numeric(min_n_sq_det) && length(min_n_sq_det) == 1 &&
               min_n_sq_det >= 0))
if (identical(model_A, model_B)) {
  stop("model_A and model_B are the same source; nothing to compare.")
}

safe_name <- function(x) {
  stringr::str_replace_all(x, "[^A-Za-z0-9]+", "_") |>
    stringr::str_replace_all("^_|_$", "")
}

# Canonical species key for joining (case / whitespace insensitive).
norm_sp <- function(x) stringr::str_squish(tolower(x))

# ============================================================
# Loaders -> common schema
#
# Every loader returns one row per (sp_english, region_key) with:
#   source, sp_english, sp_key, region_key, region_label, ci_level,
#   pc_median/pc_lower/pc_upper   (total % change),
#   at_median/at_lower/at_upper   (annual %/yr; NA if unavailable),
#   has_annual
# ============================================================

load_atlas <- function(name) {
  rds <- file.path(est_root, paste0("change_estimates_", name),
                   paste0("regional_change_estimates_", name, ".rds"))
  csv <- sub("\\.rds$", ".csv", rds)

  if (file.exists(rds)) {
    df <- readRDS(rds)
  } else if (file.exists(csv)) {
    df <- readr::read_csv(csv, show_col_types = FALSE)
  } else {
    stop("No script-08 estimates for model '", name, "'. Looked for:\n  ",
         rds, "\n  ", csv, call. = FALSE)
  }

  req <- c("sp_english", "Region_Number", "Region_Name",
           "pct_change_median", "pct_change_qlow", "pct_change_qhigh", "ci_level")
  miss <- setdiff(req, names(df))
  if (length(miss)) {
    stop("Model '", name, "' estimates are missing columns: ",
         paste(miss, collapse = ", "), call. = FALSE)
  }

  has_annual <- all(c("annual_trend_pct_median", "annual_trend_pct_qlow",
                      "annual_trend_pct_qhigh") %in% names(df))

  df %>%
    dplyr::mutate(Region_Number = as.character(Region_Number)) %>%
    dplyr::transmute(
      source       = name,
      sp_english   = sp_english,
      sp_key       = norm_sp(sp_english),
      region_key   = dplyr::if_else(Region_Number == "province",
                                    "province", paste0("BCR", Region_Number)),
      region_label = Region_Name,
      ci_level     = ci_level,
      pc_median    = pct_change_median,
      pc_lower     = pct_change_qlow,
      pc_upper     = pct_change_qhigh,
      at_median    = if (has_annual) annual_trend_pct_median else NA_real_,
      at_lower     = if (has_annual) annual_trend_pct_qlow   else NA_real_,
      at_upper     = if (has_annual) annual_trend_pct_qhigh  else NA_real_,
      has_annual   = has_annual
    )
}

# Per-region squares-with-detection counts, read from a model's FULL script-08
# export (the paired export does not carry them). Returns one row per
# (sp_key, region_key). Used only for the data-sufficiency filter.
load_detection_counts <- function(name) {
  rds <- file.path(est_root, paste0("change_estimates_", name),
                   paste0("regional_change_estimates_", name, ".rds"))
  csv <- sub("\\.rds$", ".csv", rds)

  if (file.exists(rds)) {
    df <- readRDS(rds)
  } else if (file.exists(csv)) {
    df <- readr::read_csv(csv, show_col_types = FALSE)
  } else {
    stop("Detection filter needs a full-model export for '", name,
         "', but none was found. Looked for:\n  ", rds, "\n  ", csv,
         call. = FALSE)
  }

  req <- c("sp_english", "Region_Number", "n_sq_det_OBBA2", "n_sq_det_OBBA3")
  miss <- setdiff(req, names(df))
  if (length(miss)) {
    stop("Full-model export for '", name, "' lacks detection columns: ",
         paste(miss, collapse = ", "),
         ". Re-run script 08, which writes n_sq_det_OBBA2 / n_sq_det_OBBA3.",
         call. = FALSE)
  }

  df %>%
    dplyr::mutate(Region_Number = as.character(Region_Number)) %>%
    dplyr::transmute(
      sp_key         = norm_sp(sp_english),
      region_key     = dplyr::if_else(Region_Number == "province",
                                      "province", paste0("BCR", Region_Number)),
      n_sq_det_OBBA2 = as.integer(n_sq_det_OBBA2),
      n_sq_det_OBBA3 = as.integer(n_sq_det_OBBA3)
    ) %>%
    dplyr::distinct(sp_key, region_key, .keep_all = TRUE)
}

# BBS quantile-column suffixes for a supported ci_level.
bbs_quantile_suffix <- function(ci) {
  if (isTRUE(all.equal(ci, 0.90))) return(list(lo = "0.05",  hi = "0.95"))
  if (isTRUE(all.equal(ci, 0.95))) return(list(lo = "0.025", hi = "0.975"))
  stop("BBS comparison supports ci_level 0.90 or 0.95 only (got ", ci, "). ",
       "Re-run script 08 at one of those levels, or set ci_level accordingly.",
       call. = FALSE)
}

load_bbs <- function(path, ci) {
  if (!file.exists(path)) {
    stop("BBS trends file not found at: ", path, call. = FALSE)
  }

  # Sniff comma vs tab from the header line.
  header <- readLines(path, n = 1L, warn = FALSE)
  delim  <- if (grepl("\t", header) && !grepl(",", header)) "\t" else ","
  bbs <- readr::read_delim(path, delim = delim, show_col_types = FALSE,
                           na = c("", "NA"), guess_max = 100000)

  q <- bbs_quantile_suffix(ci)
  pc_lo <- paste0("percent_change_q_", q$lo)
  pc_hi <- paste0("percent_change_q_", q$hi)
  tr_lo <- paste0("trend_q_", q$lo)
  tr_hi <- paste0("trend_q_", q$hi)

  req <- c("region", "region_type", "english_name",
           "percent_change", pc_lo, pc_hi, "trend", tr_lo, tr_hi)
  miss <- setdiff(req, names(bbs))
  if (length(miss)) {
    stop("BBS file is missing columns: ", paste(miss, collapse = ", "),
         call. = FALSE)
  }

  if (!is.null(bbs_years)) {
    stopifnot(all(c("start_year", "end_year") %in% names(bbs)))
    bbs <- bbs %>%
      dplyr::filter(start_year == bbs_years[["start"]],
                    end_year   == bbs_years[["end"]])
  }

  # region / region_type -> canonical key. Everything else is dropped.
  region_map <- tibble::tribble(
    ~region,     ~region_type,  ~region_key, ~region_label,
    "ON",        "prov_state",  "province",  "Province-wide (BBS ON)",
    "CA-ON-12",  "stratum",     "BCR12",     "BCR 12",
    "CA-ON-13",  "stratum",     "BCR13",     "BCR 13"
  )

  out <- bbs %>%
    dplyr::inner_join(region_map, by = c("region", "region_type")) %>%
    dplyr::transmute(
      source       = "BBS",
      sp_english   = english_name,
      sp_key       = norm_sp(english_name),
      region_key   = region_key,
      region_label = region_label,
      ci_level     = ci,
      pc_median    = percent_change,
      pc_lower     = .data[[pc_lo]],
      pc_upper     = .data[[pc_hi]],
      at_median    = trend,
      at_lower     = .data[[tr_lo]],
      at_upper     = .data[[tr_hi]],
      has_annual   = TRUE
    )

  # Guard against multiple rows per species x region (e.g. several intervals).
  dup <- out %>% dplyr::count(sp_key, region_key) %>% dplyr::filter(n > 1)
  if (nrow(dup) > 0) {
    warning(nrow(dup), " BBS species x region combinations had multiple rows; ",
            "keeping the first of each. Set bbs_years to disambiguate.",
            call. = FALSE)
    out <- out %>% dplyr::distinct(sp_key, region_key, .keep_all = TRUE)
  }

  out
}

# Paired (shared-footprint) estimates from script 08. These are per-BCR only
# (no province row) and carry a fixed 90% interval (pct_change_q05/q50/q95), so
# `ci` is not chosen here -- it is always 0.90. Optionally annualised when
# paired_interval_years is set (monotone transform of the % quantiles).
#
# The paired analysis is model-independent, so this reads the single
# change_estimates/paired_change_estimates.rds -- there is no per-model variant.
load_paired <- function() {
  rds <- file.path(est_root, "change_estimates", "paired_change_estimates.rds")
  csv <- sub("\\.rds$", ".csv", rds)

  if (file.exists(rds)) {
    df <- readRDS(rds)
  } else if (file.exists(csv)) {
    df <- readr::read_csv(csv, show_col_types = FALSE)
  } else {
    stop("No script-08 paired estimates found. Looked for:\n  ",
         rds, "\n  ", csv, "\nRun script 08 with INCLUDE_PAIRED = TRUE ",
         "(after running 07b).", call. = FALSE)
  }

  req <- c("sp_english", "Region_Number", "Region_Name",
           "pct_change_q50", "pct_change_q05", "pct_change_q95")
  miss <- setdiff(req, names(df))
  if (length(miss)) {
    stop("Paired estimates are missing columns: ",
         paste(miss, collapse = ", "),
         ". (Expected the shared_change_summary q05/q50/q95 columns.)",
         call. = FALSE)
  }

  has_annual <- !is.null(paired_interval_years)
  if (has_annual) {
    stopifnot(is.numeric(paired_interval_years), paired_interval_years > 0)
    annualise_pct <- function(pct) {
      100 * ((1 + pct / 100)^(1 / paired_interval_years) - 1)
    }
  }

  df %>%
    dplyr::mutate(Region_Number = as.character(Region_Number)) %>%
    dplyr::transmute(
      source       = "paired",
      sp_english   = sp_english,
      sp_key       = norm_sp(sp_english),
      region_key   = dplyr::if_else(Region_Number == "province",
                                    "province", paste0("BCR", Region_Number)),
      region_label = Region_Name,
      ci_level     = 0.90,
      pc_median    = pct_change_q50,
      pc_lower     = pct_change_q05,
      pc_upper     = pct_change_q95,
      at_median    = if (has_annual) annualise_pct(pct_change_q50) else NA_real_,
      at_lower     = if (has_annual) annualise_pct(pct_change_q05) else NA_real_,
      at_upper     = if (has_annual) annualise_pct(pct_change_q95) else NA_real_,
      has_annual   = has_annual
    )
}

# ============================================================
# Resolve sources and load
# ============================================================

# Parse a model_A / model_B string into a source spec. "paired" is a single
# model-independent source; there is no "paired:<model>" form any more.
parse_spec <- function(x) {
  if (toupper(x) == "BBS")    return(list(kind = "bbs",    model = NA_character_))
  if (tolower(x) == "paired") return(list(kind = "paired", model = NA_character_))
  list(kind = "atlas", model = x)
}

spec_A <- parse_spec(model_A)
spec_B <- parse_spec(model_B)

# Human-readable source labels (used for columns, plot axes, and file names).
label_of <- function(spec) {
  switch(spec$kind,
         bbs    = "BBS",
         atlas  = spec$model,
         paired = "paired")
}
label_A <- label_of(spec_A)
label_B <- label_of(spec_B)
if (identical(label_A, label_B)) {
  stop("Both sources resolve to '", label_A, "'; nothing to compare.",
       call. = FALSE)
}

# Cache atlas loads so ci inference doesn't reread.
atlas_cache <- new.env(parent = emptyenv())
get_atlas <- function(m) {
  if (is.null(atlas_cache[[m]])) atlas_cache[[m]] <- load_atlas(m)
  atlas_cache[[m]]
}

# ---- Target ci_level -------------------------------------------------------
# Priority: explicit ci_level > forced 0.90 if a paired source is present >
# the atlas source ci_level > 0.90.
atlas_models <- unique(c(
  if (spec_A$kind == "atlas") spec_A$model,
  if (spec_B$kind == "atlas") spec_B$model
))
paired_present <- spec_A$kind == "paired" || spec_B$kind == "paired"

if (!is.null(ci_level)) {
  target_ci <- ci_level
  if (paired_present && !isTRUE(all.equal(target_ci, 0.90))) {
    warning("A paired source is present but ci_level = ", target_ci,
            "; paired carries only a 90% interval, so its bounds are used as-is ",
            "and are NOT at the ", target_ci, " level.", call. = FALSE)
  }
} else if (paired_present) {
  target_ci <- 0.90            # paired only has q05/q95
  if (length(atlas_models) > 0) {
    a_ci <- get_atlas(atlas_models[[1]])$ci_level[1]
    if (!isTRUE(all.equal(a_ci, 0.90))) {
      warning("Atlas model was summarised at ci_level ", a_ci, ", but a paired ",
              "source forces the comparison to 0.90. Re-run script 08 at 0.90 ",
              "for a like-for-like interval.", call. = FALSE)
    }
  }
} else if (length(atlas_models) > 0) {
  cis <- vapply(atlas_models, function(m) get_atlas(m)$ci_level[1], numeric(1))
  if (length(unique(round(cis, 6))) > 1) {
    warning("The two atlas models were summarised at different ci_levels (",
            paste(unique(cis), collapse = ", "), "); their intervals are not on ",
            "the same nominal level. Medians remain comparable.", call. = FALSE)
  }
  target_ci <- cis[[1]]
} else {
  target_ci <- 0.90
}

# ---- Load each side --------------------------------------------------------
load_spec <- function(spec) {
  switch(spec$kind,
         atlas  = get_atlas(spec$model),
         paired = load_paired(),
         bbs    = load_bbs(bbs_csv_path, target_ci))
}

A <- load_spec(spec_A)
B <- load_spec(spec_B)

# Optional crosswalk on source-B names.
if (!is.null(species_crosswalk_csv)) {
  xwalk <- readr::read_csv(species_crosswalk_csv, show_col_types = FALSE)
  stopifnot(all(c("from", "to") %in% names(xwalk)))
  xwalk <- xwalk %>% dplyr::mutate(from_key = norm_sp(from), to_key = norm_sp(to))
  B <- B %>%
    dplyr::left_join(xwalk[, c("from_key", "to_key")],
                     by = c("sp_key" = "from_key")) %>%
    dplyr::mutate(sp_key = dplyr::coalesce(to_key, sp_key)) %>%
    dplyr::select(-to_key)
}

# Pick the metric columns.
metric_cols <- if (compare_metric == "percent_change") {
  c(med = "pc_median", lo = "pc_lower", hi = "pc_upper")
} else {
  if (!all(A$has_annual) || !all(B$has_annual)) {
    stop("compare_metric = 'annual_trend' but annual columns are absent on at ",
         "least one source. For atlas models, re-run script 08 with atlas ",
         "reference years set; for a paired source, set paired_interval_years.",
         call. = FALSE)
  }
  c(med = "at_median", lo = "at_lower", hi = "at_upper")
}

metric_label <- if (compare_metric == "percent_change") {
  "Total change (%)"
} else {
  "Annual trend (%/yr)"
}

prep_side <- function(d, tag) {
  d %>%
    dplyr::transmute(
      sp_key, region_key,
      sp_english_side = sp_english,
      region_label_side = region_label,
      !!paste0("est_", tag)  := .data[[metric_cols[["med"]]]],
      !!paste0("lo_", tag)   := .data[[metric_cols[["lo"]]]],
      !!paste0("hi_", tag)   := .data[[metric_cols[["hi"]]]]
    )
}

A2 <- prep_side(A, "A")
B2 <- prep_side(B, "B")

# ============================================================
# Report unmatched species (within comparable regions)
# ============================================================

comparable_regions <- intersect(unique(A2$region_key), unique(B2$region_key))
message("Comparable regions: ", paste(comparable_regions, collapse = ", "))

A_sp <- A2 %>% dplyr::filter(region_key %in% comparable_regions) %>%
  dplyr::distinct(sp_key, sp_english_side)
B_sp <- B2 %>% dplyr::filter(region_key %in% comparable_regions) %>%
  dplyr::distinct(sp_key, sp_english_side)

only_A <- dplyr::anti_join(A_sp, B_sp, by = "sp_key") %>%
  dplyr::transmute(present_in = label_A, sp_english = sp_english_side)
only_B <- dplyr::anti_join(B_sp, A_sp, by = "sp_key") %>%
  dplyr::transmute(present_in = label_B, sp_english = sp_english_side)
unmatched <- dplyr::bind_rows(only_A, only_B)

message("Species only in ", label_A, ": ", nrow(only_A))
message("Species only in ", label_B, ": ", nrow(only_B))

# ============================================================
# Join and derive comparison quantities
# ============================================================

classify_dir <- function(lo, hi) {
  dplyr::case_when(
    is.na(lo) | is.na(hi) ~ NA_character_,
    lo > 0 ~ "increase",
    hi < 0 ~ "decrease",
    TRUE   ~ "uncertain"
  )
}

comparison <- A2 %>%
  dplyr::inner_join(
    B2 %>% dplyr::select(sp_key, region_key, est_B, lo_B, hi_B,
                         region_label_B = region_label_side),
    by = c("sp_key", "region_key")
  ) %>%
  dplyr::transmute(
    sp_english   = sp_english_side,
    region_key,
    region_label = region_label_side,
    source_A = label_A, source_B = label_B,
    metric   = compare_metric,
    ci_level = target_ci,
    est_A, lo_A, hi_A,
    est_B, lo_B, hi_B,
    dir_A = classify_dir(lo_A, hi_A),
    dir_B = classify_dir(lo_B, hi_B),
    diff_median      = est_A - est_B,
    abs_diff_median  = abs(est_A - est_B),
    ci_overlap       = (lo_A <= hi_B) & (lo_B <= hi_A),
    direction_agree  = dir_A == dir_B,
    agreement = dplyr::case_when(
      dir_A == "uncertain" | dir_B == "uncertain" ~ "One/both uncertain",
      dir_A == dir_B ~ "Agree",
      TRUE ~ "Disagree"
    )
  ) %>%
  dplyr::mutate(
    region_key = factor(region_key,
                        levels = c("province", "BCR12", "BCR13",
                                   "BCR74", "BCR76", "BCR77"))
  ) %>%
  dplyr::arrange(region_key, sp_english) %>%
  dplyr::mutate(region_key = as.character(region_key))

# ============================================================
# Detection-based data-sufficiency filter
#
# Attach the per-region squares-with-detection counts and (optionally) drop
# species x region combinations that fall below min_n_sq_det. Counts come from a
# model's full script-08 export, keyed to comparison rows by normalised species
# name + region_key. They are added to the output for transparency even when
# min_n_sq_det is NULL.
# ============================================================

# Which model supplies the counts: explicit config, else an atlas source. A
# paired source is model-independent and carries no counts, so it cannot supply
# them; BBS carries none either.
counts_model <- if (!is.null(detection_counts_model)) {
  detection_counts_model
} else {
  cand <- c(
    if (spec_A$kind == "atlas") spec_A$model,
    if (spec_B$kind == "atlas") spec_B$model
  )
  if (length(cand)) cand[[1]] else NA_character_
}

if (is.na(counts_model)) {
  if (!is.null(min_n_sq_det)) {
    stop("min_n_sq_det is set, but neither source is an atlas model to supply ",
         "detection counts (e.g. paired vs BBS). Set detection_counts_model.",
         call. = FALSE)
  }
} else {
  det_counts <- load_detection_counts(counts_model)

  comparison <- comparison %>%
    dplyr::mutate(sp_key = norm_sp(sp_english)) %>%
    dplyr::left_join(det_counts, by = c("sp_key", "region_key")) %>%
    dplyr::mutate(
      counts_matched = !is.na(n_sq_det_OBBA2) | !is.na(n_sq_det_OBBA3),
      n_sq_det_OBBA2 = dplyr::coalesce(n_sq_det_OBBA2, 0L),
      n_sq_det_OBBA3 = dplyr::coalesce(n_sq_det_OBBA3, 0L)
    ) %>%
    dplyr::select(-sp_key) %>%
    dplyr::relocate(n_sq_det_OBBA2, n_sq_det_OBBA3, counts_matched,
                    .after = region_label)

  if (!is.null(min_n_sq_det)) {
    keep <- if (identical(detection_rule, "either")) {
      comparison$n_sq_det_OBBA2 >= min_n_sq_det |
        comparison$n_sq_det_OBBA3 >= min_n_sq_det
    } else {
      comparison$n_sq_det_OBBA2 >= min_n_sq_det &
        comparison$n_sq_det_OBBA3 >= min_n_sq_det
    }

    n_drop <- sum(!keep)
    n_miss <- sum(!keep & !comparison$counts_matched)
    rule_phrase <- if (identical(detection_rule, "either")) "either atlas" else "both atlases"
    message(sprintf(
      "Detection filter (>= %d squares detected in %s; counts from '%s'): dropped %d of %d rows%s.",
      as.integer(min_n_sq_det), rule_phrase, counts_model,
      n_drop, nrow(comparison),
      if (n_miss > 0)
        sprintf(" (%d had no detection-count match and were treated as 0)", n_miss)
      else ""
    ))

    comparison <- comparison[keep, , drop = FALSE]
  }
}

if (nrow(comparison) == 0) {
  stop("No species x region rows remain to compare between ", label_A, " and ",
       label_B, ". Check species naming (see the unmatched report), region ",
       "overlap, and the min_n_sq_det threshold.", call. = FALSE)
}

message("Matched ", dplyr::n_distinct(comparison$sp_english), " species across ",
        nrow(comparison), " species x region rows.")

# ============================================================
# Per-region and overall agreement summary
# ============================================================

summarise_block <- function(d, label) {
  safe_cor <- function(x, y, method) {
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 3) return(NA_real_)
    tryCatch(stats::cor(x[ok], y[ok], method = method), error = function(e) NA_real_)
  }
  tibble::tibble(
    region_key           = label,
    n                    = nrow(d),
    pearson_r            = safe_cor(d$est_A, d$est_B, "pearson"),
    spearman_rho         = safe_cor(d$est_A, d$est_B, "spearman"),
    bias_A_minus_B       = mean(d$diff_median, na.rm = TRUE),
    median_abs_diff      = stats::median(d$abs_diff_median, na.rm = TRUE),
    rmsd                 = sqrt(mean(d$diff_median^2, na.rm = TRUE)),
    prop_ci_overlap      = mean(d$ci_overlap, na.rm = TRUE),
    prop_direction_agree = mean(d$direction_agree, na.rm = TRUE)
  )
}

summary_tbl <- dplyr::bind_rows(
  comparison %>%
    dplyr::group_split(region_key) %>%
    lapply(function(d) summarise_block(d, unique(d$region_key))) %>%
    dplyr::bind_rows(),
  summarise_block(comparison, "ALL")
)

print(summary_tbl, n = nrow(summary_tbl))

# ============================================================
# Save
# ============================================================

tag       <- paste0(safe_name(label_A), "_vs_", safe_name(label_B))
csv_path  <- file.path(out_dir, paste0("comparison_", tag, ".csv"))
rds_path  <- file.path(out_dir, paste0("comparison_", tag, ".rds"))
sum_path  <- file.path(out_dir, paste0("comparison_summary_", tag, ".csv"))
un_path   <- file.path(out_dir, paste0("comparison_unmatched_", tag, ".csv"))

readr::write_csv(comparison, csv_path)
saveRDS(comparison, rds_path)
readr::write_csv(summary_tbl, sum_path)
readr::write_csv(unmatched, un_path)
message("Wrote: ", csv_path)
message("Wrote: ", rds_path)
message("Wrote: ", sum_path)
message("Wrote: ", un_path)

# ============================================================
# Comparison scatter (optional)
# ============================================================

if (make_plot) {

  # Okabe-Ito, CVD-safe: Agree = blue, Disagree = golden, Uncertain = grey.
  agree_cols <- c(
    "Agree"              = "#0072B2",
    "Disagree"           = "#E69F00",
    "One/both uncertain" = "grey65"
  )

  # Per-region facet labels carrying rho + n.
  facet_lab <- summary_tbl %>%
    dplyr::filter(region_key != "ALL") %>%
    dplyr::mutate(
      facet = sprintf("%s  (rho = %.2f, n = %d)", region_key, spearman_rho, n)
    ) %>%
    dplyr::select(region_key, facet)

  plot_df <- comparison %>%
    dplyr::left_join(facet_lab, by = "region_key") %>%
    dplyr::mutate(
      agreement = factor(agreement,
                         levels = c("Agree", "Disagree", "One/both uncertain")),
      region_key = factor(region_key,
                          levels = c("province", "BCR12", "BCR13",
                                     "BCR74", "BCR76", "BCR77"))
    ) %>%
    dplyr::arrange(region_key)
  plot_df$facet <- factor(plot_df$facet, levels = unique(plot_df$facet))

  # Symmetric axis: warp both axes by log2 fold-change, log2(1 + x/100), so a
  # doubling (+100%) and a halving (-50%) are equidistant from 0. The transform
  # is identical on both axes, so the 1:1 line and the zero lines stay correct.
  # Tick labels remain in percent. Values <= -100% have no finite fold change and
  # are dropped (with a warning) rather than plotted at -Inf.
  sym_on <- TRUE
  to_sym <- function(v) log2(1 + v / 100)

  plot_df <- plot_df %>%
    dplyr::mutate(
      x    = if (sym_on) to_sym(est_B) else est_B,
      y    = if (sym_on) to_sym(est_A) else est_A,
      xmin = if (sym_on) to_sym(lo_B)  else lo_B,
      xmax = if (sym_on) to_sym(hi_B)  else hi_B,
      ymin = if (sym_on) to_sym(lo_A)  else lo_A,
      ymax = if (sym_on) to_sym(hi_A)  else hi_A
    )

  # Percent-labelled breaks placed at their log2 fold-change positions.
  fold_pct <- 100 * (2^(-6:6) - 1)          # ..., -75, -50, 0, +100, +300, ...
  fine_pct <- c(-25, -10, -5, 5, 10, 25)    # extra detail near 0
  nice_pct <- sort(unique(round(c(0, fold_pct, fine_pct))))
  brk_pos  <- to_sym(nice_pct)
  brk_lab  <- ifelse(nice_pct == 0, "0%", sprintf("%+.0f%%", nice_pct))
  sym_scales <- list(
    ggplot2::scale_x_continuous(breaks = brk_pos, labels = brk_lab),
    ggplot2::scale_y_continuous(breaks = brk_pos, labels = brk_lab)
  )

  x_lab <- sprintf("%s  -  %s", label_B, metric_label)
  y_lab <- sprintf("%s  -  %s", label_A, metric_label)

  scale_note <- if (sym_on) "; symmetric log2 fold-change axis" else ""

  # Per-facet SQUARE axis range: the same [min, max] on x and y, spanning every
  # point and error-bar end (finite values only; a symmetric-axis -Inf from a
  # -100% change is ignored here). Invisible corner points at (lo, lo) and
  # (hi, hi) force each free-scaled panel to share one range across both axes;
  # combined with square panels (aspect.ratio = 1) this makes the 1:1 line a
  # true 45-degree diagonal.
  finite_rng <- function(...) {
    v <- c(...)
    v <- v[is.finite(v)]
    if (length(v) == 0) c(NA_real_, NA_real_) else range(v)
  }
  axis_rng <- plot_df %>%
    dplyr::group_by(facet) %>%
    dplyr::summarise(
      lo = finite_rng(x, y, xmin, xmax, ymin, ymax)[1],
      hi = finite_rng(x, y, xmin, xmax, ymin, ymax)[2],
      .groups = "drop"
    )
  square_df <- dplyr::bind_rows(
    dplyr::transmute(axis_rng, facet, x = lo, y = lo),
    dplyr::transmute(axis_rng, facet, x = hi, y = hi)
  )

  # One label per species. ggrepel keeps labels off each other and off the
  # points; if it is not installed, fall back to overlap-thinned geom_text.
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    label_layer <- ggrepel::geom_text_repel(
      ggplot2::aes(label = sp_english),
      size = 2, colour = "grey20", show.legend = FALSE,
      max.overlaps = Inf, min.segment.length = 0,
      segment.size = 0.2, segment.colour = "grey75"
    )
  } else {
    message("Package 'ggrepel' not installed; using geom_text with overlap ",
            "thinning. Install ggrepel for non-overlapping species labels.")
    label_layer <- ggplot2::geom_text(
      ggplot2::aes(label = sp_english),
      size = 2, colour = "grey20", show.legend = FALSE,
      vjust = -0.6, check_overlap = TRUE
    )
  }

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_blank(data = square_df, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.3) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey40") +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = ymin, ymax = ymax, colour = agreement),
      width = 0, alpha = 0.35, linewidth = 0.3
    ) +
    ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = xmin, xmax = xmax, colour = agreement),
      height = 0, alpha = 0.35, linewidth = 0.3
    ) +
    ggplot2::geom_point(ggplot2::aes(colour = agreement), size = 1.4, alpha = 0.9) +
    label_layer +
    sym_scales +
    ggplot2::scale_colour_manual(values = agree_cols, drop = FALSE, name = "Direction") +
    ggplot2::facet_wrap(~ facet, scales = "free") +
    ggplot2::labs(
      title    = sprintf("Regional population change: %s vs %s", label_A, label_B),
      subtitle = sprintf("%s; %.0f%% credible intervals; dashed line is 1:1%s",
                         metric_label, 100 * target_ci, scale_note),
      x = x_lab, y = y_lab
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid      = ggplot2::element_blank(),  # omit grid lines
      aspect.ratio    = 1,                          # square panels -> 1:1 at 45 deg
      plot.title      = ggplot2::element_text(face = "bold"),
      strip.text      = ggplot2::element_text(size = 9),
      legend.position = "bottom"
    )

  n_facet <- dplyr::n_distinct(plot_df$facet)
  n_col   <- min(3, n_facet)
  n_row   <- ceiling(n_facet / n_col)

  pdf_path <- file.path(out_dir, paste0("comparison_", tag, ".pdf"))
  ggplot2::ggsave(pdf_path, p, width = 4.5 * n_col, height = 4.2 * n_row,
                  limitsize = FALSE)
  message("Wrote: ", pdf_path)

  print(p)
}

message("10_compare_regional_change_estimates.R complete")
