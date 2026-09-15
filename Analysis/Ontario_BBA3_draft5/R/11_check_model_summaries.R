library(dplyr)
library(ggplot2)
library(ggrepel)

rm(list = ls())
source(here::here("R", "00_config_paths.R"))
out_dir  <- paths$model_output


model_name = "PC_ARU"
model_summaries_path <- file.path(
  out_dir, paste0("summaries_", model_name), "model_summaries.rds"
)

# INLA's awkward quantile column names -> syntactic names
quant_rename <- c("0.025quant" = "q0.025",
                  "0.5quant"   = "q0.5",
                  "0.975quant" = "q0.975")

# ---- load (skip if model_summaries is already in scope) --------------------
model_summaries <- readRDS(model_summaries_path)

# ---- flatten every species' summary_hyperpar into one long table ----------
hyper_df <- purrr::imap(
  .x = model_summaries,
  .f = function(sp, species_name) {
    if (!is.list(sp)) return(NULL)
    hp <- sp$summary_hyperpar
    if (is.null(hp) || NROW(hp) == 0) return(NULL)
    
    hp |>
      as.data.frame(check.names = FALSE) |>              # matrix/df -> df, keep names
      tibble::rownames_to_column(var = "parameter") |>
      dplyr::mutate(species = species_name, .before = 1)
  }
) |>
  purrr::compact() |>          # drop species with no hyperpar (NULL / failed fits)
  dplyr::bind_rows()           # aligns by column name; missing cols -> NA

if (nrow(hyper_df) == 0) {
  stop("No summary_hyperpar tables found in model_summaries.")
}

# ---- tidy column names -----------------------------------------------------
hit <- names(hyper_df) %in% names(quant_rename)
names(hyper_df)[hit] <- quant_rename[names(hyper_df)[hit]]

# ---- parse the parameter name + derive RE variances ------------------------
# stat_type is the leading word ("Precision" / "Stdev" / "Range" / "size").
# For iid random effects INLA reports PRECISION, so add SD/variance. These use
# the MEDIAN (q0.5) because the median is invariant under the monotone
# precision -> sd/variance map, so 1/sqrt(median_prec) is exactly the median SD.
# (Inverting the mean would be wrong: 1/mean(prec) != mean(1/prec).)
# The 95% SD interval flips the precision tails: high precision = low SD.
hyper_df <- hyper_df |>
  dplyr::mutate(
    stat_type = sub(pattern = " .*$", replacement = "", x = parameter),
    component = dplyr::case_when(
      grepl("^Precision for ", parameter) ~ sub("^Precision for ", "", parameter),
      grepl("^Stdev for ",     parameter) ~ sub("^Stdev for ",     "", parameter),
      grepl("^Range for ",     parameter) ~ sub("^Range for ",     "", parameter),
      TRUE                                 ~ parameter
    ),
    is_precision = stat_type == "Precision",
    sd_median  = ifelse(is_precision, 1 / sqrt(q0.5),    NA_real_),
    var_median = ifelse(is_precision, 1 / q0.5,          NA_real_),
    sd_q0.025  = ifelse(is_precision, 1 / sqrt(q0.975),  NA_real_),  # tails flip
    sd_q0.975  = ifelse(is_precision, 1 / sqrt(q0.025),  NA_real_)
  )

# ---- roster of species that fit -------------------------------------------
species_fit <- sort(unique(hyper_df$species))
n_fit       <- length(species_fit)
message(n_fit, " species with a fitted PC_ARU model:\n  ",
        paste(species_fit, collapse = ", "))

# ---- ensure native-scale median + scale label exist -----------------------
# precision rows -> SD (1/sqrt(median prec)); size/range/stdev -> the value.
if (!"estimate_median" %in% names(hyper_df)) {
  hyper_df <- hyper_df |>
    dplyr::mutate(
      estimate_median = dplyr::if_else(is_precision, 1 / sqrt(q0.5), q0.5),
      estimate_scale  = dplyr::case_when(
        is_precision         ~ "sd",
        stat_type == "size"  ~ "size",
        stat_type == "Range" ~ "range",
        TRUE                 ~ "stdev"
      )
    )
}

# ---- per-parameter cross-species summary (native scale) -------------------
param_summary <- hyper_df |>
  dplyr::group_by(parameter) |>
  dplyr::summarise(
    scale       = dplyr::first(estimate_scale),
    n_species   = dplyr::n(),                        # fitted species carrying this row
    min         = min(estimate_median,    na.rm = TRUE),
    p10         = quantile(estimate_median, 0.10, na.rm = TRUE),
    median      = median(estimate_median, na.rm = TRUE),
    p90         = quantile(estimate_median, 0.90, na.rm = TRUE),
    max         = max(estimate_median,    na.rm = TRUE),
    post_cv_med = median(sd / abs(mean),  na.rm = TRUE),  # relative posterior width
    .groups = "drop"
  ) |>
  dplyr::mutate(coverage = n_species / n_fit) |>
  dplyr::arrange(dplyr::desc(coverage), parameter)

# ---- readable copy to paste back ------------------------------------------
param_summary_disp <- param_summary |>
  dplyr::mutate(dplyr::across(c(min, p10, median, p90, max, post_cv_med),
                              ~ signif(.x, 3)))
print(param_summary_disp, n = Inf, width = Inf)


per_species_wide <- hyper_df |>
  dplyr::select(species, parameter, estimate_median) |>
  tidyr::pivot_wider(names_from = parameter, values_from = estimate_median)
print(per_species_wide, n = Inf, width = Inf)

hist(per_species_wide$`Range for spde_mean`)

library(ggplot2)
ggplot(per_species_wide, 
       aes(x = `Range for spde_mean`, y = `Stdev for spde_mean`,
           label = species))+
  geom_vline(xintercept = 200, linetype = 2)+
  geom_text_repel(col = "gray50", size = 2) +
  geom_point()+
  theme_bw()+
  scale_x_continuous(trans = "log10")

ggplot(per_species_wide, aes(x = `Range for spde_diff`, y = `Stdev for spde_diff`,
                             label = species))+
  geom_vline(xintercept = 200, linetype = 2)+
  geom_text_repel(col = "gray50", size = 2) +
  geom_point()+
  theme_bw()+
  scale_x_continuous(trans = "log10")


subset(per_species_wide, `Range for spde_mean` <= 150 | `Range for spde_diff` <= 150)
subset(per_species_wide, `Range for spde_diff` <= 150)