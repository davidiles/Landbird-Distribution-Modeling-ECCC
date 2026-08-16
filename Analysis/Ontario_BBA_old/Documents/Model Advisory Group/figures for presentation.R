
setwd("C:/Users/IlesD/OneDrive - EC-EC/Iles/Projects/Landbirds/Landbird-Distribution-Modeling-ECCC/Analysis/Ontario_BBA")

library(ggplot2)

set.seed(123)

n <- 200
draws <- rnorm(n, mean = 0.1, sd = 0.1)
med <- median(draws)
ci <- quantile(draws, c(0.025, 0.975))

dat <- data.frame(
  x = draws,
  y = runif(n, -0.07, 0.07)
)

ggplot(dat, aes(x, y)) +
  annotate(
    "rect",
    xmin = ci[1], xmax = ci[2],
    ymin = -0.11, ymax = 0.11,
    fill = "gray50", alpha = 0.25
  ) +
  geom_point(size = 2.2, col = "gray50") +
  geom_vline(xintercept = med, linewidth = 1.1, col = "gray50") +
  geom_vline(xintercept = 0, linetype = 2) +
  
  labs(
    x = "Posterior draws for one pixel",
    y = NULL,
    title = "Single estimate vs full range of plausible values"
  ) +
  coord_cartesian(ylim = c(-0.16, 0.16)) +
  theme_classic(base_size = 14) +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank()
  )

# -------------------------------------------------------------------------------
# Comparison of different rasters

# Load a species relative abundance raster
library(terra)
library(sf)

r <- rast("data_clean/model_output/rasters_PC_ARU/bay_breasted_warbler_a3.tif")

# Layer to plot
abund <- r$mu_q50


# -----------------------------
# 1. Continuous map
# -----------------------------
# 99th percentile
q99 <- quantile(values(abund), probs = 0.99, na.rm = TRUE)

abund_cap <- clamp(abund, lower = 0, upper = q99)

# Colour ramp
colramp <- colorRampPalette(c(
  "#FEFEFE", "#FBF7E2", "#FCF8D0", "#EEF7C2", "#CEF2B0",
  "#94E5A0", "#51C987", "#18A065", "#008C59", "#007F53", "#006344"
))(100)


plot(
  abund_cap,
  col = colramp,
  type = "continuous",
  range = c(0, q99),   # cap colour scale at 99th percentile
  main = "Relative abundance (continuous)",
  colNA = "grey90"
)

# -----------------------------
# 2. Discrete map: 5 equal-quantile bins
# -----------------------------

# Colour ramp (4 bins)
cols4 <- colorRampPalette(c(
  "#FEFEFE", "#FBF7E2", "#FCF8D0", "#EEF7C2", "#CEF2B0",
  "#94E5A0", "#51C987", "#18A065", "#008C59", "#007F53", "#006344"
))(4)

# 99th percentile
q99 <- quantile(values(abund), probs = 0.99, na.rm = TRUE)

# Lower threshold
low <- 0.01

# One midpoint between low and q99
mid <- (low + q99) / 2

# Classification matrix: 4 bins total
rcl <- rbind(
  c(0,   low, 1),
  c(low, mid, 2),
  c(mid, q99, 3),
  c(q99, Inf, 4)
)

abund_binned <- classify(
  abund,
  rcl,
  include.lowest = TRUE,
  right = FALSE
)

# Labels
fmt <- function(x) format(round(x, 3), nsmall = 3)

labels <- c(
  paste0("0 – ", fmt(low)),
  paste0(fmt(low), " – ", fmt(mid)),
  paste0(fmt(mid), " – ", fmt(q99)),
  paste0("> ", fmt(q99))
)

# Plot
plot(
  abund_binned,
  col = cols4,
  type = "classes",
  levels = labels,
  main = "Relative abundance (4 bins)",
  colNA = "grey90"
)
