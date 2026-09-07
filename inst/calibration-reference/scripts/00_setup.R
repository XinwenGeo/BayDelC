# BayDelC v0.3 clean project setup -----------------------------------------
# Run this file from the BayDelC project root.
# It loads the active dependencies and public PSM interface only. Plotting
# files, legacy user interfaces, and model-fitting scripts are never sourced.

if (!isTRUE(getOption("baydelc.setup.loaded", FALSE))) {

  required_project_files <- c(
    "R/utils.R",
    "R/model_bundle.R",
    "R/forward_psm.R",
    "R/inverse_psm.R",
    "config/model_registry.R"
  )

  missing_project_files <- required_project_files[
    !file.exists(required_project_files)
  ]

  if (length(missing_project_files) > 0L) {
    stop(
      paste0(
        "BayDelC setup must be run from the project root. Missing: ",
        paste(missing_project_files, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  required_packages <- c(
    # Core data/model dependencies
    "dplyr", "readr", "readxl", "tibble", "purrr", "tidyr",
    "brms", "rstan", "loo", "posterior",
    "sessioninfo", "jsonlite"
  )

  package_available <- vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )

  missing_packages <- required_packages[!package_available]

  if (length(missing_packages) > 0L) {
    stop(
      paste0(
        "Missing required R package(s): ",
        paste(missing_packages, collapse = ", "),
        ". Install them before running BayDelC."
      ),
      call. = FALSE
    )
  }

  suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(readxl)
    library(tibble)
    library(purrr)
    library(tidyr)
    library(brms)
    library(rstan)
    library(loo)
    library(posterior)
    library(sessioninfo)
    library(jsonlite)
  })

  # Active core helpers. Order matters because model_bundle.R defines %||%.
  source("R/utils.R")
  source("R/model_bundle.R")
  source("config/model_registry.R")

  # Public user-facing PSM interface. Legacy forward_user.R/inverse_user.R are
  # intentionally not loaded in the clean project.
  source("R/forward_psm.R")
  source("R/inverse_psm.R")

  ensure_dirs()

  # Preserve a reproducibility marker without changing its timestamp every
  # time an individual figure script sources this setup file.
  if (!file.exists("VERSION.json")) {
    write_version_file("VERSION.json")
  }

  rstan::rstan_options(auto_write = TRUE)
  options(mc.cores = max(1L, parallel::detectCores(logical = TRUE)))
  options(baydelc.setup.loaded = TRUE)

  message("BayDelC v0.3 clean setup loaded from: ", normalizePath(".", winslash = "/"))
}

invisible(TRUE)
