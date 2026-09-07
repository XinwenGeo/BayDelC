# Bundled prediction-model registry ---------------------------------------

.baydelc_registry_data <- function() {
  data.frame(
    model_id = c(
      "BLR_Sub_main", "BLR_All_main", "GPR_Sub_main", "GPR_All_main",
      "BLR_Sub_StudentT", "BLR_All_StudentT",
      "GPR_Sub_shrinkage", "GPR_All_shrinkage",
      "GPR_Sub_xjitter", "GPR_All_xjitter", "GPR_Sub_latent"
    ),
    model_set = c(rep("main", 4), rep("sensitivity", 7)),
    algorithm = c("BLR", "BLR", "GPR", "GPR", "BLR", "BLR",
                  "GPR", "GPR", "GPR", "GPR", "GPR"),
    calibration = c("Sub", "All", "Sub", "All", "Sub", "All",
                    "Sub", "All", "Sub", "All", "Sub"),
    variant = c("baseline", "baseline", "baseline", "baseline",
                "student_t", "student_t", "shrinkage", "shrinkage",
                "x_jitter", "x_jitter", "latent_input"),
    display_name = c(
      "BLR Gaussian main calibration (restricted taxa)",
      "BLR Gaussian main calibration (mixed taxa)",
      "GPR baseline main calibration (restricted taxa)",
      "GPR baseline main calibration (mixed taxa)",
      "BLR Student-t sensitivity calibration (restricted taxa)",
      "BLR Student-t sensitivity calibration (mixed taxa)",
      "GPR shrinkage sensitivity calibration (restricted taxa)",
      "GPR shrinkage sensitivity calibration (mixed taxa)",
      "GPR oxygen-jitter ensemble diagnostic (restricted taxa)",
      "GPR oxygen-jitter ensemble diagnostic (mixed taxa)",
      "GPR latent-input sensitivity calibration (restricted taxa)"
    ),
    file = c(
      "Sub_Linear_EIV_Gaussian_v0.3_portable.rds",
      "All_Linear_EIV_Gaussian_v0.3_portable.rds",
      "Sub_GPR_Baseline_v0.3_portable.rds",
      "All_GPR_Baseline_v0.3_portable.rds",
      "Sub_Linear_EIV_StudentT_v0.3_portable.rds",
      "All_Linear_EIV_StudentT_v0.3_portable.rds",
      "Sub_GPR_Shrinkage_v0.3_portable.rds",
      "All_GPR_Shrinkage_v0.3_portable.rds",
      "Sub_GPR_XJitter_Baseline_v0.3_diagnostic.rds",
      "All_GPR_XJitter_Baseline_v0.3_diagnostic.rds",
      "Sub_GPR_LatentInput_v0.3_portable.rds"
    ),
    default = c(TRUE, rep(FALSE, 10)),
    prediction_ready = c(rep(TRUE, 8), FALSE, FALSE, TRUE),
    stringsAsFactors = FALSE
  )
}

.baydelc_model_dir <- function(model_set) {
  installed <- system.file("models", model_set, package = "BayDelC")
  if (nzchar(installed)) return(installed)
  candidate <- file.path("inst", "models", model_set)
  if (dir.exists(candidate)) return(candidate)
  ""
}

#' List the calibrated models bundled with BayDelC
#'
#' @param model_set Which model collection to list: `"all"`, `"main"`, or
#'   `"sensitivity"`.
#' @param prediction_only If `TRUE`, omit diagnostic-only assets that cannot be
#'   passed to `forward_psm()` or `inverse_psm()`.
#' @return A data frame describing the available model and calibration choices.
#' @export
baydelc_models <- function(model_set = c("all", "main", "sensitivity"),
                           prediction_only = FALSE) {
  model_set <- match.arg(model_set)
  registry <- .baydelc_registry_data()
  if (model_set != "all") registry <- registry[registry$model_set == model_set, , drop = FALSE]
  if (isTRUE(prediction_only)) registry <- registry[registry$prediction_ready, , drop = FALSE]
  registry$available <- vapply(seq_len(nrow(registry)), function(i) {
    root <- .baydelc_model_dir(registry$model_set[[i]])
    nzchar(root) && file.exists(file.path(root, registry$file[[i]]))
  }, logical(1))
  registry[c("model_id", "model_set", "algorithm", "calibration", "variant",
             "display_name", "default", "prediction_ready", "available")]
}

.baydelc_resolve_registry_row <- function(algorithm, calibration,
                                           model_set = "main", variant = NULL) {
  algorithm <- match.arg(toupper(algorithm), c("BLR", "GPR"))
  calibration <- match.arg(calibration, c("Sub", "All"))
  model_set <- match.arg(tolower(model_set), c("main", "sensitivity"))
  if (is.null(variant)) {
    variant <- if (model_set == "main") "baseline" else if (algorithm == "BLR") "student_t" else "shrinkage"
  }
  variant <- tolower(variant)
  registry <- .baydelc_registry_data()
  hit <- registry[
    registry$algorithm == algorithm & registry$calibration == calibration &
      registry$model_set == model_set & registry$variant == variant,
    , drop = FALSE
  ]
  if (nrow(hit) != 1L) {
    available <- registry[
      registry$algorithm == algorithm & registry$calibration == calibration &
        registry$model_set == model_set,
      "variant", drop = TRUE
    ]
    stop(sprintf(
      "No %s/%s/%s model has variant `%s`. Available variant(s): %s.",
      model_set, algorithm, calibration, variant,
      if (length(available)) paste(available, collapse = ", ") else "none"
    ), call. = FALSE)
  }
  hit
}

.baydelc_model_path <- function(algorithm, calibration, model_set = "main",
                                 variant = NULL, must_exist = TRUE) {
  hit <- .baydelc_resolve_registry_row(algorithm, calibration, model_set, variant)
  path <- file.path(.baydelc_model_dir(hit$model_set[[1L]]), hit$file[[1L]])
  if (must_exist && (!nzchar(path) || !file.exists(path))) {
    stop(sprintf("Bundled model asset `%s` was not found.", hit$model_id[[1L]]), call. = FALSE)
  }
  path
}

#' Load a bundled BayDelC model or diagnostic asset
#'
#' @param algorithm Either `"BLR"` or `"GPR"`.
#' @param calibration Either restricted-taxon `"Sub"` or mixed-taxon `"All"`.
#' @param model_set Either the primary `"main"` collection or
#'   `"sensitivity"` collection.
#' @param variant Model variant. Defaults to `"baseline"` for main models,
#'   `"student_t"` for sensitivity BLR, and `"shrinkage"` for sensitivity GPR.
#' @return A lightweight calibrated model or diagnostic bundle.
#' @export
baydelc_load_model <- function(algorithm = c("BLR", "GPR"),
                               calibration = c("Sub", "All"),
                               model_set = c("main", "sensitivity"),
                               variant = NULL) {
  algorithm <- match.arg(algorithm)
  calibration <- match.arg(calibration)
  model_set <- match.arg(model_set)
  readRDS(.baydelc_model_path(algorithm, calibration, model_set, variant))
}

#' Inspect metadata for a bundled BayDelC model
#'
#' @inheritParams baydelc_load_model
#' @return A named list containing model identity, collection, variant,
#'   calibration domain, version, draw count, and prediction backend.
#' @export
baydelc_model_info <- function(algorithm = c("BLR", "GPR"),
                               calibration = c("Sub", "All"),
                               model_set = c("main", "sensitivity"),
                               variant = NULL) {
  algorithm <- match.arg(algorithm)
  calibration <- match.arg(calibration)
  model_set <- match.arg(model_set)
  row <- .baydelc_resolve_registry_row(algorithm, calibration, model_set, variant)
  model <- readRDS(.baydelc_model_path(algorithm, calibration, model_set, variant))
  list(
    model_id = row$model_id[[1L]],
    model_name = model$model_name,
    model_type = model$model_type,
    model_set = row$model_set[[1L]],
    variant = row$variant[[1L]],
    prediction_ready = row$prediction_ready[[1L]],
    model_version = model$bundle_version %||% NA_character_,
    data_version = model$data_version %||% NA_character_,
    calibration_range = model$calibration_range %||% range(model$summary_curve$O2),
    posterior_draws = if (is.null(model$posterior_draws)) NA_integer_ else nrow(as.data.frame(model$posterior_draws)),
    prediction_backend = model$prediction_backend %||% "diagnostic curve ensemble"
  )
}
