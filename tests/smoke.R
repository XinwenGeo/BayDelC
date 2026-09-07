library(BayDelC)

stopifnot(
  identical(baydelc_version(), as.character(packageVersion("BayDelC"))),
  nrow(baydelc_models()) == 4L,
  all(baydelc_models()$available)
)

for (algorithm in c("BLR", "GPR")) {
  for (calibration in c("Sub", "All")) {
    info <- baydelc_model_info(algorithm, calibration)
    stopifnot(identical(info$prediction_backend,
                        if (algorithm == "BLR") "posterior_parameter_draws"
                        else "interpolated_posterior_mean_curves"))

    f <- forward_psm(
      bwo = c(80, 150, 220), algorithm = algorithm,
      calibration = calibration, n_draw = 100, seed = 11,
      return_draws = TRUE
    )
    stopifnot(inherits(f, "forward_psm_result"), nrow(f$summary) == 3L,
              all(is.finite(f$summary$median)), all(dim(f$draws) == c(3L, 100L)))

    inv <- inverse_psm(
      dd13c = 1.5, dd13c_sd = 0.08, algorithm = algorithm,
      calibration = calibration,
      inversion = if (algorithm == "BLR") "strict" else "grid",
      prior = if (algorithm == "BLR") "none" else "uniform_physical",
      n_draw = 100, integration_draws = 50, seed = 12,
      return_draws = TRUE
    )
    stopifnot(inherits(inv, "inverse_psm_result"),
              is.finite(inv$summary$median), all(dim(inv$draws) == c(1L, 100L)))
  }
}

component <- psm_error_components(
  epi_measurement_sd = 0.04, infa_measurement_sd = 0.04
)
z <- inverse_psm(dd13c = 1.5, error_components = component,
                 error_method = "components", n_draw = 100)
stopifnot(is.finite(z$summary$median))
