# ============================================================
# extract_spde_hyperpars.R   (rename to fit your numbering)
#
# Purpose
#   Read the per-species model summaries, pull the SPDE range + marginal SD out
#   of each species' hyperparameter table, and assemble one tidy table for
#   exploring the distribution of ranges across species -- so you can spot the
#   species whose spatial field is far shorter than the rest (and far from your
#   range prior).
#
#   Why range diagnoses "patchy"/"blobby" maps: for a Matern field the practical
#   RANGE is ~the distance at which spatial correlation decays to ~0.13. A SHORT
#   posterior range means neighbouring cells decorrelate almost immediately, so
#   the surface is isolated hotspots rather than a smooth cline -- the blobby
#   look. Most species pull toward a long range (smooth); the short-range
#   minority are the ones you're seeing. Marginal SD sets feature amplitude, so
#   short range + high SD is the most extreme patchy corner.
#
#   Two spatial fields appear per species: "spde_mean" (baseline ABUNDANCE
#   surface -- the one that drives abundance-map smoothness) and "spde_diff"
#   (the between-atlas change field). Both land in the long table; the flag
#   table below targets SPATIAL_COMPONENT (default spde_mean).
# ============================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(here)
})

source(here::here("R", "00_config_paths.R"))

# ============================================================
# CONFIG
# ============================================================
model_name <- "PC_ARU"

# Which spatial field to diagnose. "spde_mean" = baseline abundance surface;
# "spde_diff" = change field. Everything with a "Range for X" row is a candidate.
SPATIAL_COMPONENT <- "spde_mean"

# pc-matern prior anchor for THAT component (native mesh units; your ranges look
# like km). prior.range = c(PRIOR_RANGE, p) => P(range < PRIOR_RANGE) = p, so it
# is a lower-tail yardstick for "did the data pull the range far below where I
# set it". NA => skip the prior-based flag (mirrors 08's atlas-year NA-skip).
PRIOR_RANGE <- 300
PRIOR_SIGMA <- 1.5

# "Short/unexpected" thresholds.
LOWER_Q            <- 0.10   # flag species in the bottom 10% of ranges (self-referential)
PRIOR_RATIO_THRESH <- 0.25   # flag median range < 25% of PRIOR_RANGE (if set)
RANGE_FLOOR        <- NA_real_  # optional absolute floor in native units (NA = skip)
TOP_N_REPORT       <- 25

MAKE_PLOTS <- TRUE

# ------------------------------------------------------------
out_dir  <- file.path(paths$model_output, paste0("spde_hyperpars_", model_name))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

model_summaries_path <- file.path(
  paths$model_output, paste0("summaries_", model_name), "model_summaries.rds"
)
stopifnot(file.exists(model_summaries_path))

okabe_ito <- c(orange = "#E69F00", skyblue = "#56B4E9", green = "#009E73",
               yellow = "#F0E442", blue = "#0072B2", vermillion = "#D55E00",
               purple = "#CC79A7", black = "#000000")

# ============================================================
# Extract Range + Stdev rows from one hyperparameter table
# ============================================================
# INLA quantile columns -> tidy names; other columns untouched.
std_names <- function(d) {
  dplyr::rename_with(d, .fn = function(nm) dplyr::case_when(
    nm == "0.025quant" ~ "q025",
    nm == "0.5quant"   ~ "median",
    nm == "0.975quant" ~ "q975",
    TRUE               ~ nm
  ))
}

extract_range_sd <- function(sh) {
  if (is.null(sh) || nrow(sh) == 0) return(tibble::tibble())
  tibble::as_tibble(as.data.frame(sh), rownames = "param_row") |>
    std_names() |>
    dplyr::filter(stringr::str_detect(param_row, "^(Range|Stdev) for ")) |>
    dplyr::mutate(
      param     = dplyr::if_else(stringr::str_starts(param_row, "Range"), "range", "sd"),
      component = stringr::str_remove(param_row, "^(Range|Stdev) for ")
    )
}

# ============================================================
# Iterate species (the list is keyed by English name)
# ============================================================
model_summaries <- readRDS(model_summaries_path)
message("Species in model_summaries: ", length(model_summaries))

spde_long <- purrr::imap_dfr(model_summaries, function(ms, sp_english) {
  extract_range_sd(ms$summary_hyperpar) |>
    dplyr::mutate(sp_english = sp_english, model_name = model_name, .before = 1)
})
if (nrow(spde_long) == 0) stop("No Range/Stdev rows found in any summary_hyperpar.")

message("Components seen: ", paste(sort(unique(spde_long$component)), collapse = ", "))

long_csv <- file.path(out_dir, paste0("spde_hyperpars_long_", model_name, ".csv"))
utils::write.csv(spde_long, long_csv, row.names = FALSE)
# saveRDS(spde_long, file.path(out_dir, paste0("spde_hyperpars_long_", model_name, ".rds")))
message("Wrote: ", long_csv)

# ============================================================
# Spatial-component wide table + short-range flags
# ============================================================
range_components <- unique(spde_long$component[spde_long$param == "range"])
if (!SPATIAL_COMPONENT %in% range_components) {
  stop("SPATIAL_COMPONENT '", SPATIAL_COMPONENT, "' not found. Available: ",
       paste(range_components, collapse = ", "))
}

spatial_wide <- spde_long |>
  dplyr::filter(component == SPATIAL_COMPONENT, param %in% c("range", "sd")) |>
  dplyr::select(sp_english, param, mean, median, q025, q975) |>
  tidyr::pivot_wider(names_from = param,
                     values_from = c(mean, median, q025, q975),
                     names_glue = "{param}_{.value}")

range_cut <- 100 # OR use.... stats::quantile(spatial_wide$range_median, LOWER_Q, na.rm = TRUE)

spatial_wide <- spatial_wide |>
  dplyr::mutate(
    range_vs_prior   = if (is.finite(PRIOR_RANGE)) range_median / PRIOR_RANGE else NA_real_,
    sd_vs_prior      = if (is.finite(PRIOR_SIGMA)) sd_median    / PRIOR_SIGMA else NA_real_,
    flag_short_pack  = range_median < range_cut,
    flag_short_prior = if (is.finite(PRIOR_RANGE)) range_vs_prior < PRIOR_RATIO_THRESH else NA,
    flag_short_abs   = if (is.finite(RANGE_FLOOR)) range_median < RANGE_FLOOR else NA,
    flag_patchy = dplyr::coalesce(flag_short_pack,  FALSE) |
      dplyr::coalesce(flag_short_prior, FALSE) |
      dplyr::coalesce(flag_short_abs,   FALSE)
  ) |>
  dplyr::arrange(range_median)

wide_csv <- file.path(out_dir, paste0("spde_spatial_range_sd_", model_name, ".csv"))
utils::write.csv(spatial_wide, wide_csv, row.names = FALSE)
# saveRDS(spatial_wide, file.path(out_dir, paste0("spde_spatial_range_sd_", model_name, ".rds")))
message("Wrote: ", wide_csv)

# ============================================================
# Console exploration
# ============================================================
message("\n--- posterior median RANGE across species (", SPATIAL_COMPONENT, ") ---")
print(round(stats::quantile(spatial_wide$range_median,
                            c(.05,.25,.5,.75,.95), na.rm = TRUE)))
if (is.finite(PRIOR_RANGE)) message("prior range anchor: ", PRIOR_RANGE)
message("bottom-", round(100*LOWER_Q), "% cutoff: ", round(range_cut))
message("flagged short/patchy: ", sum(spatial_wide$flag_patchy, na.rm = TRUE),
        " of ", nrow(spatial_wide))

message("\n--- ", min(TOP_N_REPORT, nrow(spatial_wide)), " shortest-range species ---")
spatial_wide |>
  dplyr::transmute(sp_english,
                   range_median = round(range_median),
                   sd_median    = round(sd_median, 2),
                   vs_prior     = round(range_vs_prior, 2),
                   patchy       = flag_patchy) |>
  utils::head(TOP_N_REPORT) |>
  print(n = TOP_N_REPORT)

# ============================================================
# Optional Okabe-Ito figures
# ============================================================
if (isTRUE(MAKE_PLOTS) && requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  fig_dir <- if (!is.null(paths$figures)) paths$figures else out_dir
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  
  p_dist <- ggplot(spatial_wide, aes(x = range_median)) +
    geom_histogram(bins = 40, fill = okabe_ito[["skyblue"]], colour = "white") +
    geom_rug(aes(colour = flag_patchy), alpha = 0.7, show.legend = FALSE) +
    scale_x_log10() +
    scale_colour_manual(values = c(`FALSE` = okabe_ito[["black"]],
                                   `TRUE`  = okabe_ito[["vermillion"]])) +
    labs(title = paste0("Spatial SPDE range across species (", SPATIAL_COMPONENT, ")"),
         subtitle = "Left tail = patchy/blobby maps; red rug = flagged",
         x = "Posterior median range (native units, log10)", y = "Species") +
    theme_minimal(base_size = 12)
  if (is.finite(PRIOR_RANGE)) p_dist <- p_dist +
    geom_vline(xintercept = PRIOR_RANGE, linetype = "dashed", colour = okabe_ito[["orange"]])
  
  p_scatter <- ggplot(spatial_wide,
                      aes(x = range_median, y = sd_median, colour = flag_patchy)) +
    geom_point(alpha = 0.8, size = 2) +
    scale_x_log10() + scale_y_log10() +
    scale_colour_manual(values = c(`FALSE` = okabe_ito[["skyblue"]],
                                   `TRUE`  = okabe_ito[["vermillion"]]), name = "patchy") +
    labs(title = "Spatial field: range vs marginal SD",
         subtitle = "Left (short range) = blobby; high SD adds contrast",
         x = "Median range (log10)", y = "Median SD (log10)") +
    theme_minimal(base_size = 12)
  
  # ggsave(file.path(fig_dir, paste0("spde_range_distribution_", model_name, ".png")),
  #        p_dist, width = 8, height = 5, dpi = 150)
  # ggsave(file.path(fig_dir, paste0("spde_range_vs_sd_", model_name, ".png")),
  #        p_scatter, width = 7, height = 6, dpi = 150)
  # message("Wrote figures to: ", fig_dir)
} else if (isTRUE(MAKE_PLOTS)) {
  message("ggplot2 not available; skipped plots.")
}


print(p_scatter)
subset(spatial_wide, flag_patchy == TRUE) %>% as.data.frame()

message("\nextract_spde_hyperpars.R complete")