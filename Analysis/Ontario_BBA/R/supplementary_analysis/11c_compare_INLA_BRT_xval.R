suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(pROC)
  library(ggplot2)
})

cv_dir_INLA <- "data_clean/model_output/xval_INLA"
cv_dir_BRT  <- "data_clean/model_output/xval_BRT"

# Choose one species or many:
fitted_species <- c("Bobolink",
                    "Blue Jay",
                    "Canada Jay",
                    "Olive-sided Flycatcher",
                    "Winter Wren",
                    "Lesser Yellowlegs",
                    "Blackpoll Warbler",
                    "Connecticut Warbler",
                    "Palm Warbler",
                    "Lincoln's Sparrow",
                    "Fox Sparrow",
                    "Common Nighthawk",
                    "Long-eared Owl",
                    "American Tree Sparrow",
                    "LeConte's Sparrow",
                    "Nelson's Sparrow",
                    "Boreal Chickadee",
                    "Rusty Blackbird",
                    "Yellow-bellied Flycatcher",
                    "Greater Yellowlegs",
                    "Hudsonian Godwit",
                    "Canada Warbler",
                    "Eastern Wood-Peewee",
                    "Grasshopper Sparrow",
                    "Solitary Sandpiper",
                    "White-throated Sparrow",
                    "Bay-breasted Warbler")

# ---- storage objects
cv_comp_summaries <- vector("list", length(fitted_species))  # each element = 1-row tibble
names(cv_comp_summaries) <- fitted_species

# Optional: store the joined block-level comparison table per species for debugging
cv_comp_details <- vector("list", length(fitted_species))
names(cv_comp_details) <- fitted_species

auc_safe <- function(y, p) {
  y <- as.integer(y)
  ok <- is.finite(p) & !is.na(y)
  y <- y[ok]; p <- p[ok]
  if (length(unique(y)) < 2) return(NA_real_)
  as.numeric(pROC::auc(pROC::roc(y, p, quiet = TRUE)))
}

for (i in seq_along(fitted_species)) {
  
  sp_english <- fitted_species[i]
  message("\n====================\n", i, "/", length(fitted_species), ": ", sp_english, "\n====================")
  
  # ------------------------------------------------------------
  # Build filenames
  # ------------------------------------------------------------
  sp_file <- sp_english %>%
    str_to_lower() %>%
    str_replace_all("[^a-z0-9]+", "_") %>%
    str_replace_all("^_|_$", "")
  
  cv_path_block_INLA <- file.path(cv_dir_INLA, "block_summaries", paste0(sp_file, ".rds"))
  cv_path_block_BRT  <- file.path(cv_dir_BRT,  "block_summaries", paste0(sp_file, ".rds"))
  
  if (!file.exists(cv_path_block_INLA)) {
    message("Missing INLA block summary: ", cv_path_block_INLA)
    next
  }
  if (!file.exists(cv_path_block_BRT)) {
    message("Missing BRT block summary: ", cv_path_block_BRT)
    next
  }
  
  # ------------------------------------------------------------
  # Load block summaries
  # ------------------------------------------------------------
  cv_block_summ_INLA <- readRDS(cv_path_block_INLA)
  cv_block_summ_BRT  <- readRDS(cv_path_block_BRT)
  
  # ------------------------------------------------------------
  # Join (fold/block level paired comparisons)
  # ------------------------------------------------------------
  comp <- cv_block_summ_INLA %>%
    st_drop_geometry() %>%
    select(sp_english, sp_code, block_size_km, rep, fold, Atlas, block_id,
           n = n,
           mean_obs,
           rmse_inla = rmse,
           mae_inla  = mae,
           bias_inla = bias,
           pred_mu_inla = mean_pred_mean,
           p_any_detect_inla = p_any_detect,
           obs_any_detect_inla = obs_any_detect) %>%
    inner_join(
      cv_block_summ_BRT %>%
        st_drop_geometry() %>%
        select(sp_english, sp_code, block_size_km, rep, fold, Atlas, block_id,
               rmse_brt = rmse,
               mae_brt  = mae,
               bias_brt = bias,
               pred_mu_brt = mean_pred_mean,
               p_any_detect_brt = p_any_detect),
      by = c("sp_english","sp_code","block_size_km","rep","fold","Atlas","block_id")
    ) %>%
    mutate(
      d_rmse = rmse_brt - rmse_inla,  # >0 => INLA better (lower rmse)
      d_mae  = mae_brt  - mae_inla,   # >0 => INLA better (lower mae)
      d_bias = bias_brt - bias_inla
    )
  
  if (nrow(comp) == 0) {
    message("No overlapping rows between INLA and BRT for this species (check keys).")
    next
  }
  
  # ------------------------------------------------------------
  # Summarize metrics (overall)
  # ------------------------------------------------------------
  summ_overall <- comp %>%
    summarize(
      sp_english = first(sp_english),
      sp_code = first(sp_code),
      block_size_km = first(block_size_km),
      n_blocks_eval = n(),
      
      # Observed data
      mean_count_per_block = mean(mean_obs,na.rm = TRUE),
      prop_blocks_detetected = mean(mean_obs>0,na.rm = TRUE),
      
      # Crossvalidation metrics
      mean_rmse_inla = mean(rmse_inla, na.rm = TRUE),
      mean_rmse_brt  = mean(rmse_brt,  na.rm = TRUE),
      mean_d_rmse    = mean(d_rmse,    na.rm = TRUE),
      brt_win_rate_rmse = mean(d_rmse < 0, na.rm = TRUE),
      
      mean_mae_inla  = mean(mae_inla, na.rm = TRUE),
      mean_mae_brt   = mean(mae_brt,  na.rm = TRUE),
      mean_d_mae     = mean(d_mae,    na.rm = TRUE),
      brt_win_rate_mae = mean(d_mae < 0, na.rm = TRUE),
      
      cor_inla = suppressWarnings(cor(mean_obs, pred_mu_inla, use = "complete.obs")),
      cor_brt  = suppressWarnings(cor(mean_obs, pred_mu_brt,  use = "complete.obs")),
      
      slope_inla = coef(lm(mean_obs ~ pred_mu_inla, data = comp))[2],
      slope_brt  = coef(lm(mean_obs ~ pred_mu_brt,  data = comp))[2],
      
      auc_inla = auc_safe(obs_any_detect_inla, p_any_detect_inla),
      auc_brt  = auc_safe(obs_any_detect_inla, p_any_detect_brt)
    )
  
  # Optional: also summarize separately by Atlas
  summ_by_atlas <- comp %>%
    group_by(Atlas) %>%
    summarize(
      sp_english = first(sp_english),
      sp_code = first(sp_code),
      block_size_km = first(block_size_km),
      n_blocks_eval = n(),
      
      # Observed data
      mean_count_per_block = mean(mean_obs,na.rm = TRUE),
      prop_blocks_detetected = mean(mean_obs>0,na.rm = TRUE),
      
      # Crossvalidation metrics
      mean_rmse_inla = mean(rmse_inla, na.rm = TRUE),
      mean_rmse_brt  = mean(rmse_brt,  na.rm = TRUE),
      mean_d_rmse    = mean(d_rmse,    na.rm = TRUE),
      brt_win_rate_rmse = mean(d_rmse < 0, na.rm = TRUE),
      
      mean_mae_inla  = mean(mae_inla, na.rm = TRUE),
      mean_mae_brt   = mean(mae_brt,  na.rm = TRUE),
      mean_d_mae     = mean(d_mae,    na.rm = TRUE),
      brt_win_rate_mae = mean(d_mae < 0, na.rm = TRUE),
      
      cor_inla = suppressWarnings(cor(mean_obs, pred_mu_inla, use = "complete.obs")),
      cor_brt  = suppressWarnings(cor(mean_obs, pred_mu_brt,  use = "complete.obs")),
      
      slope_inla = coef(lm(mean_obs ~ pred_mu_inla, data = cur_data()))[2],
      slope_brt  = coef(lm(mean_obs ~ pred_mu_brt,  data = cur_data()))[2],
      
      auc_inla = auc_safe(obs_any_detect_inla, p_any_detect_inla),
      auc_brt  = auc_safe(obs_any_detect_inla, p_any_detect_brt),
      .groups = "drop"
    )
  
  # ------------------------------------------------------------
  # Store results
  # ------------------------------------------------------------
  cv_comp_summaries[[sp_english]] <- list(
    overall = summ_overall,
    by_atlas = summ_by_atlas
  )
  cv_comp_details[[sp_english]] <- comp
  
  # For single species “walk-through”, print key summaries:
  print(summ_overall)
  print(summ_by_atlas)
}

# Combine into data frames after loop:
summary_overall_tbl <- bind_rows(lapply(cv_comp_summaries, `[[`, "overall"))
summary_by_atlas_tbl <- bind_rows(lapply(cv_comp_summaries, `[[`, "by_atlas"))

summary_overall_tbl
summary_by_atlas_tbl %>% as.data.frame()


# View xval results for a single species
dat <- cv_comp_details$`Canada Warbler`

ggplot(dat, aes(x = n, y = d_rmse)) +
  geom_point(alpha = 0.4) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_smooth()+
  scale_x_log10() +
  labs(x = "Surveys per block (log scale)", y = "d_rmse = RMSE_BRT - RMSE_INLA")

lim <- range(dat[,c("mean_obs","pred_mu_brt","pred_mu_inla")])
ggplot(data =  dat)+
  geom_point(aes(x = mean_obs, y = pred_mu_inla, size = n))+
  geom_abline(intercept = 0, slope = 1)+
  coord_cartesian(xlim = lim, ylim = lim)

ggplot(data = dat)+
  geom_point(aes(x = mean_obs, y = pred_mu_brt, size = n))+
  geom_abline(intercept = 0, slope = 1)+
  coord_cartesian(xlim = lim, ylim = lim)


