# ============================================================
# 00_config_paths.R
#
# Purpose
#   Centralize project and shared-data paths for the Ontario_BBA workflow.
#   Every numbered script sources this file first and then uses the `paths` list.
# ============================================================

suppressPackageStartupMessages({
  library(here)
})

paths <- list()

# Project root (.../Analysis/Ontario_BBA when Ontario_BBA.Rproj is open).
paths$project <- here::here()

# Repo root, two levels up (.../Landbird-Distribution-Modeling-ECCC).
paths$repo <- normalizePath(here::here("..", ".."), winslash = "/", mustWork = FALSE)

# Shared data folder (lives at <repo>/Data).
paths$data <- normalizePath(file.path(paths$repo, "Data"), winslash = "/", mustWork = FALSE)

# Analysis outputs (kept within the Ontario_BBA project).
paths$data_clean   <- here::here("data_clean")
paths$model_output <- here::here("data_clean", "model_output")
paths$figures      <- here::here("figures")
paths$tables       <- here::here("tables")

# Code locations inside the Ontario project.
paths$r_dir     <- here::here("R")
paths$functions <- here::here("R", "functions")

# ------------------------------------------------------------
# Safety check (fail fast with a helpful message)
# ------------------------------------------------------------

if (!dir.exists(paths$project)) {
  stop("Project root not found: ", paths$project)
}
