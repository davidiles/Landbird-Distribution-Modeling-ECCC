rm(list=ls())

library(sf)
library(dplyr)
library(ggplot2)
library(terra)
library(ggrepel)
library(naturecounts)

# ------------------------------------------------------------------------------
# Define study region
# ------------------------------------------------------------------------------

# Load study region
bcr <- st_read("../../Data/Spatial/BCR/BCR_Terrestrial_master.shp", quiet = TRUE)

study_region <- bcr %>%
  filter(PROVINCE_S == "ONTARIO") %>%
  st_make_valid() %>%
  mutate(region = paste0("CA-ON-", BCR)) %>%
  group_by(region) %>%
  summarise(do_union = TRUE, .groups = "drop") %>%
  st_make_valid()

# ------------------------------------------------------------------------------
# Load Atlas 20-year change estimates
# NOTE: these are currently produced by Ontario_BBA/Analysis/R/11_generate_BCR_change_estimates.R
# ------------------------------------------------------------------------------

atlas_model <- "joint2"
atlas_change_summary <- readRDS(paste0("../Ontario_BBA/data_clean/model_output/summaries_",atlas_model,"/atlas_BCR_change_estimates.rds"))

# ------------------------------------------------------------------------------
# Load BBS trend estimates (2002 - 2022)
# ------------------------------------------------------------------------------

BBS_trend_summary <- read.csv(
  "../../Data/bbs_trend_estimates/Ontario_Atlas_trends_means_by_atlas_period_2002_2022.csv"
) %>%
  rename(
    BBS_q025 = percent_change_q_0.025,
    BBS_mean = percent_change,
    BBS_q975 = percent_change_q_0.975,
    n_BBS_routes_mean = mean_n_routes
  ) %>%
  select(
    english_name, region, n_routes, n_BBS_routes_mean,
    BBS_mean, BBS_q025, BBS_q975
  ) %>%
  filter(
    n_BBS_routes_mean >= 10,
    region %in% unique(study_region$region)
  )


# ------------------------------------------------------------------------------
# Load BAM-based change estimates (processed by 1_process_BAM_rasters.R)
# ------------------------------------------------------------------------------

bam_summary <- readRDS("output/bam_change_summary_2000_2020.rds")

# ------------------------------------------------------------------------------
# Join all estimates together into a single wide dataframe
# ------------------------------------------------------------------------------

# Convert percent change to log (and vice versa)
pct_to_log <- function(pct) log1p(pct / 100)
log_to_pct <- function(x) 100 * (exp(x) - 1)

# Join dataframes
comparison <- bam_summary %>%
  dplyr::rename(english_name = "bam_commonName") %>%
  full_join(BBS_trend_summary) %>%
  
  full_join(
    atlas_change_summary %>% 
      dplyr::select(english_name,region,n_sq_OBBA2,n_sq_OBBA3,atlas_median,atlas_q025,atlas_q975)
  ) %>%
  
  mutate(region = factor(region, levels = c("CA-ON-13","CA-ON-12","CA-ON-8","CA-ON-7")),
         
         atlas_median_log = pct_to_log(atlas_median),
         atlas_q025_log = pct_to_log(atlas_q025),
         atlas_q975_log = pct_to_log(atlas_q975),
         
         bam_median_log = pct_to_log(bam_percent_change),
         
         BBS_mean_log = pct_to_log(BBS_mean),
         BBS_q025_log = pct_to_log(BBS_q025),
         BBS_q975_log = pct_to_log(BBS_q975),
         BBS_CIW      = BBS_q975 - BBS_q025) %>%
  arrange(BBS_mean)


# Optionally, limit to only species that are present in all 3 datasets
shared_spp <- bam_summary$bam_commonName %>% 
  intersect(BBS_trend_summary$english_name) %>%
  intersect(atlas_change_summary$english_name)

comparison <- comparison %>% subset(english_name %in% shared_spp)

# ------------------------------------------------------------
# Compare BBS and BAM
# ------------------------------------------------------------

n_bbs_required <- 10

cmpr1 <- comparison %>%
  filter(n_BBS_routes_mean >= n_bbs_required)

lim <- range(cmpr1[, c("BBS_mean_log","bam_median_log","atlas_median_log")],na.rm = TRUE)

# correlation label per region
cor_labels <- cmpr1 %>%
  group_by(region) %>%
  summarise(
    cor_val = cor(BBS_mean_log, bam_median_log, use = "complete.obs", method = "spearman"),
    label = paste0("Cor = ", round(cor_val, 2)),
    .groups = "drop"
  )

fig1 <- ggplot(cmpr1) +
  geom_abline(slope = 1, col = "black") +
  
  geom_smooth(aes(x = BBS_mean_log,
                  y = bam_median_log),
              method = "lm",
              col = "dodgerblue",
              fill = "dodgerblue",
              alpha = 0.2)+
  
  geom_text_repel(
    aes(
      x = BBS_mean_log,
      y = bam_median_log,
      label = english_name
    ),
    size = 2,
    max.overlaps = 20,
    col = "dodgerblue",
    alpha = 0.5
  ) +
  
  geom_text(
    data = cor_labels,
    aes(x = -Inf, y = Inf, label = label),
    hjust = -0.1,
    vjust = 1.2,
    inherit.aes = FALSE
  ) +
  
  geom_point(
    aes(
      x = BBS_mean_log,
      y = bam_median_log
    ),
    col = "dodgerblue"
  ) +
  scale_x_continuous(
    name = "BBS overall percent change from 2002 to 2022",
    labels = function(x) sprintf("%+d%%", round(log_to_pct(x)))
  ) +
  scale_y_continuous(
    name = "BAM overall percent change from 2000 to 2020",
    labels = function(x) sprintf("%+d%%", round(log_to_pct(x)))
  ) +
  coord_cartesian() +
  theme_bw() +
  coord_cartesian(xlim = lim, ylim = lim) +
  facet_grid(.~region, scales = "free") +
  ggtitle("BAM vs BBS")

fig1


cmpr1_bcr13 <- cmpr1 %>% subset(region == "CA-ON-13")
lm1 <- lm(bam_median_log ~ BBS_mean_log, data = cmpr1_bcr13)
sd(cmpr1_bcr13$bam_median_log)/sd(cmpr1_bcr13$BBS_mean_log)
mean(abs(cmpr1_bcr13$bam_median_log)) / mean(abs(cmpr1_bcr13$BBS_mean_log))

# ------------------------------------------------------------
# Compare Atlas and BBS
# ------------------------------------------------------------

n_atlas_required <- 20
cmpr2 <- comparison %>%
  filter(n_BBS_routes_mean >= n_bbs_required &
        (n_sq_OBBA2 >= n_atlas_required | n_sq_OBBA3 >= n_atlas_required))
  
# correlation label per region
cor_labels <- cmpr2 %>%
  group_by(region) %>%
  summarise(
    cor_val = cor(BBS_mean_log, atlas_median_log, use = "complete.obs", method = "spearman"),
    label = paste0("Cor = ", round(cor_val, 2)),
    .groups = "drop"
  )

fig2 <- ggplot(cmpr2) +
  geom_abline(slope = 1, col = "black") +
  
  geom_smooth(aes(x = BBS_mean_log,
                  y = atlas_median_log),
              method = "lm",
              col = "dodgerblue",
              fill = "dodgerblue",
              alpha = 0.2)+
  
  geom_text_repel(
    aes(
      x = BBS_mean_log,
      y = atlas_median_log,
      label = english_name
    ),
    size = 2,
    max.overlaps = 20,
    col = "dodgerblue",
    alpha = 0.5
  ) +
  
  geom_text(
    data = cor_labels,
    aes(x = -Inf, y = Inf, label = label),
    hjust = -0.1,
    vjust = 1.2,
    inherit.aes = FALSE
  ) +
  
  geom_point(
    aes(
      x = BBS_mean_log,
      y = atlas_median_log
    ),
    col = "dodgerblue"
  ) +
  scale_x_continuous(
    name = "BBS overall percent change from 2002 to 2022",
    labels = function(x) sprintf("%+d%%", round(log_to_pct(x)))
  ) +
  scale_y_continuous(
    name = "Atlas overall percent change from Atlas 2 to 3",
    labels = function(x) sprintf("%+d%%", round(log_to_pct(x)))
  ) +
  coord_cartesian() +
  theme_bw() +
  coord_cartesian(xlim = lim, ylim = lim) +
  facet_grid(.~region, scales = "free") +
  ggtitle("Atlas vs BBS")

fig2

cmpr2_bcr13 <- cmpr2 %>% subset(region == "CA-ON-13")
lm2 <- lm(atlas_median_log ~ BBS_mean_log, data = cmpr2_bcr13)
summary(lm2)
sd(cmpr2_bcr13$atlas_median_log)/sd(cmpr2_bcr13$BBS_mean_log)
mean(abs(cmpr2_bcr13$atlas_median_log)) / mean(abs(cmpr2_bcr13$BBS_mean_log))

# ------------------------------------------------------------
# Compare Atlas and BAM
# ------------------------------------------------------------

cmpr3 <- comparison %>%
  filter(n_sq_OBBA2 >= n_atlas_required | n_sq_OBBA3 >= n_atlas_required)

# correlation label per region
cor_labels <- cmpr3 %>%
  group_by(region) %>%
  summarise(
    cor_val = cor(bam_median_log, atlas_median_log, use = "complete.obs", method = "spearman"),
    label = paste0("Cor = ", round(cor_val, 2)),
    .groups = "drop"
  )

fig3 <- ggplot(cmpr3) +
  geom_abline(slope = 1, col = "black") +
  
  geom_smooth(aes(x = bam_median_log,
                  y = atlas_median_log),
              method = "lm",
              col = "dodgerblue",
              fill = "dodgerblue",
              alpha = 0.2)+
  
  geom_text_repel(
    aes(
      x = bam_median_log,
      y = atlas_median_log,
      label = english_name
    ),
    size = 2,
    max.overlaps = 20,
    col = "dodgerblue",
    alpha = 0.5
  ) +
  
  geom_text(
    data = cor_labels,
    aes(x = -Inf, y = Inf, label = label),
    hjust = -0.1,
    vjust = 1.2,
    inherit.aes = FALSE
  ) +
  
  geom_point(
    aes(
      x = bam_median_log,
      y = atlas_median_log
    ),
    col = "dodgerblue"
  ) +
  scale_x_continuous(
    name = "BAM overall percent change from 2000 to 2020",
    labels = function(x) sprintf("%+d%%", round(log_to_pct(x)))
  ) +
  scale_y_continuous(
    name = "Atlas overall percent change from Atlas 2 to 3",
    labels = function(x) sprintf("%+d%%", round(log_to_pct(x)))
  ) +
  coord_cartesian() +
  theme_bw() +
  coord_cartesian(xlim = lim, ylim = lim) +
  facet_grid(.~region, scales = "free") +
  ggtitle("Atlas vs BAM")

fig3

# ------------------------------------------------------------
# Compare all 3 programs in a single figure
# ------------------------------------------------------------

head(comparison)

