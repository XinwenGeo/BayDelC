# Bundled prediction-model registry ---------------------------------------

.baydelc_registry_data <- function() {
  data.frame(
    model_id = c("BLR_Sub", "BLR_All", "GPR_Sub", "GPR_All"),
    algorithm = c("BLR", "BLR", "GPR", "GPR"),
    calibration = c("Sub", "All", "Sub", "All"),
    display_name = c(
      "BLR restricted-taxon calibration",
      "BLR mixed-taxon calibration",
      "GPR restricted-taxon calibration",
      "GPR mixed-taxon calibration"
    ),
    file = c(
      "Sub_Linear_EIV_Gaussian_v0.3_portable.rds",
      "All_Linear_EIV_Gaussian_v0.3_portable.rds",
      "Sub_GPR_Baseline_v0.3_portable.rds",
      "All_GPR_Baseline_v0.3_portable.rds"
    ),
    default = c(TRUE, FALSE, FALSE, FALSE),
    stringsAsFactors = FALSE
  )
}
.baydelc_model_dir <- function() {
  installed <- system.file("models", "main", package = "BayDelC")
  if (nzchar(installed)) return(installed)

  # Development fallback for source-tree examples and tests.
  candidate <- file.path("inst", "models", "main")
  if (dir.exists(candidate)) return(candidate)
  ""
}

#' List the calibrated models bundled with BayDelC
#'
#' @return A data frame describing the available models and calibration sets.
#' @export
baydelc_models <- function() {
  registry <- .baydelc_registry_data()
  root <- .baydelc_model_dir()
  registry$available <- nzchar(root) & file.exists(file.path(root, registry$file))
  registry[c("model_id", "algorithm", "calibration", "display_name", "default", "available")]
}

.baydelc_model_path <- function(algorithm, calibration, must_exist = TRUE) {
  algorithm <- match.arg(toupper(algorithm), c("BLR", "GPR"))
  calibration <- match.arg(calibration, c("Sub", "All"))
  registry <- .baydelc_registry_data()
  hit <- registry[
    registry$algorithm == algorithm & registry$calibration == calibration,
    , drop = FALSE
  ]
  if (nrow(hit) != 1L) stop("Internal model registry is inconsistent.", call. = FALSE)
  path <- file.path(.baydelc_model_dir(), hit$file[[1L]])
  if (must_exist && (!nzchar(path) || !file.exists(path))) {
    stop(
      sprintf("Bundled %s/%s model asset was not found.", algorithm, calibration),
      call. = FALSE
    )
  }
  path
}

#' Load a bundled BayDelC prediction model
#'
#' @param algorithm Either `"BLR"` or `"GPR"`.
#' @param calibration Either the restricted-taxon `"Sub"` calibration or the
#'   mixed-taxon `"All"` calibration.
#' @return A lightweight calibrated model bundle.
#' @export
baydelc_load_model <- function(algorithm = c("BLR", "GPR"),
                               calibration = c("Sub", "All")) {
  algorithm <- match.arg(algorithm)
  calibration <- match.arg(calibration)
  readRDS(.baydelc_model_path(algorithm, calibration))
}

#' Inspect metadata for a bundled BayDelC model
#'
#' @inheritParams baydelc_load_model
#' @return A named list containing model identity, calibration domain, version,
#'   draw count, and prediction backend.
#' @export
baydelc_model_info <- function(algorithm = c("BLR", "GPR"),
                               calibration = c("Sub", "All")) {
  model <- baydelc_load_model(match.arg(algorithm), match.arg(calibration))
  list(
    model_name = model$model_name,
    model_type = model$model_type,
    model_version = model$bundle_version %||% NA_character_,
    data_version = model$data_version %||% NA_character_,
    calibration_range = model$calibration_range,
    posterior_draws = nrow(as.data.frame(model$posterior_draws)),
    prediction_backend = model$prediction_backend %||% "posterior parameters"
  )
}
