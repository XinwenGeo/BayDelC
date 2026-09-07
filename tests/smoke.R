library(BayDelC)

stopifnot(
  identical(baydelc_version(), as.character(packageVersion("BayDelC"))),
  nrow(baydelc_models()) == 11L,
  all(baydelc_models()$available)
)

registry <- baydelc_models(prediction_only = TRUE)
stopifnot(nrow(registry) == 9L)
for (i in seq_len(nrow(registry))) {
    algorithm <- registry$algorithm[[i]]
    calibration <- registry$calibration[[i]]
    model_set <- registry$model_set[[i]]
    variant <- registry$variant[[i]]
    info <- baydelc_model_info(algorithm, calibration, model_set, variant)
    stopifnot(isTRUE(info$prediction_ready))

    f <- forward_psm(
      bwo = c(80, 150, 220), algorithm = algorithm,
      calibration = calibration, model_set = model_set, variant = variant,
      n_draw = 100, seed = 11,
      return_draws = TRUE
    )
    stopifnot(inherits(f, "forward_psm_result"), nrow(f$summary) == 3L,
              all(is.finite(f$summary$median)), all(dim(f$draws) == c(3L, 100L)))

    inv <- inverse_psm(
      dd13c = 1.5, dd13c_sd = 0.08, algorithm = algorithm,
      calibration = calibration, model_set = model_set, variant = variant,
      inversion = if (algorithm == "BLR") "strict" else "grid",
      prior = if (algorithm == "BLR") "none" else "uniform_physical",
      n_draw = 100, integration_draws = 50, seed = 12,
      return_draws = TRUE
    )
    stopifnot(inherits(inv, "inverse_psm_result"),
              is.finite(inv$summary$median), all(dim(inv$draws) == c(1L, 100L)))
}

xj <- baydelc_load_model("GPR", "Sub", "sensitivity", "x_jitter")
stopifnot(identical(xj$prediction_backend, "diagnostic_curve_ensemble"))
xj_error <- try(forward_psm(bwo = 100, algorithm = "GPR", calibration = "Sub",
                            model_set = "sensitivity", variant = "x_jitter"), silent = TRUE)
stopifnot(inherits(xj_error, "try-error"))

component <- psm_error_components(
  epi_measurement_sd = 0.04, infa_measurement_sd = 0.04
)
z <- inverse_psm(dd13c = 1.5, error_components = component,
                 error_method = "components", n_draw = 100)
stopifnot(is.finite(z$summary$median))

stopifnot(
  nzchar(system.file("calibration-reference", "README.md", package = "BayDelC")),
  nzchar(system.file("examples", "paper", "Figure02_calibration_curves.R", package = "BayDelC"))
)
