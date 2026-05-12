# ============================================================
# R/00_config_paths.R
#
# Purpose:
#   Centralize project + shared-data paths for the Ontario_BBA workflow.
#   Scripts should source this file and then use the `paths` list.
# ============================================================

suppressPackageStartupMessages({
  library(here)
})

# Project root (should be .../Analysis/Ontario_BBA when Ontario_BBA.Rproj is open)
paths <- list()

paths$project <- here::here()

# Repo root is two levels up from Ontario_BBA:
#   .../Landbird-Distribution-Modeling-ECCC
paths$repo <- normalizePath(here::here("..", ".."), winslash = "/", mustWork = FALSE)

# Shared data folder (lives at repo root /Data)
paths$data <- normalizePath(file.path(paths$repo, "Data"), winslash = "/", mustWork = FALSE)

# Ontario analysis outputs (kept within the Ontario_BBA project)
paths$data_clean <- here::here("data_clean")
paths$model_output <- here::here("data_clean", "model_output")
paths$figures <- here::here("figures")
paths$tables  <- here::here("tables")

# Common code locations inside the Ontario project
paths$r_dir <- here::here("R")
paths$functions <- here::here("R", "functions")

# ------------------------------------------------------------
# Safety checks (fail fast with helpful messages)
# ------------------------------------------------------------

if (!dir.exists(paths$project)) {
  stop("Project root not found: ", paths$project)
}

if (!dir.exists(paths$data)) {
  stop(
    "Shared Data folder not found at: ", paths$data,
    "\nExpected repo layout: <repo>/Data and <repo>/Analysis/Ontario_BBA",
    "\nIf your layout differs, update paths$repo / paths$data in R/config_paths.R"
  )
}