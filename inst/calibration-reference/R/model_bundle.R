# Model bundle utilities ----------------------------------------------------

make_model_bundle <- function(model_name, model_type, data, fit, posterior_draws,
                              log_lik = NULL, y_rep = NULL, mu = NULL, loo = NULL,
                              metrics = list(), calibration_range = NULL,
                              settings = list(), notes = NULL) {
  list(
    model_name = model_name,
    model_type = model_type,
    bundle_version = baydelc_version(),
    data_version = "O2DD13C_NEW_v0.2",
    training_data = data,
    data_filter_rule = "foram-only; outlier removed; O2 >= 50 µmol/kg; complete O2/O2_sd/dd13C/dd13C_sd",
    fit_object = fit,
    posterior_draws = posterior_draws,
    log_lik_matrix = log_lik,
    posterior_predictive_draws = y_rep,
    posterior_epred_draws = mu,
    loo_object = loo,
    pareto_k = if (!is.null(loo)) tryCatch(loo$diagnostics$pareto_k, error = function(e) NA_real_) else NA_real_,
    metrics = metrics,
    calibration_range = calibration_range %||% range(data$O2, na.rm = TRUE),
    settings = settings,
    notes = notes,
    created_time = Sys.time(),
    session_info = sessioninfo::session_info()
  )
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

save_model_bundle <- function(bundle, path, legacy_rdata = TRUE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(bundle, path)
  if (legacy_rdata) {
    legacy_path <- sub("\\.rds$", ".RData", path)
    model_bundle <- bundle
    save(model_bundle, file = legacy_path)
  }
  invisible(path)
}

load_model_bundles <- function(paths) {
  lapply(paths, readRDS)
}

posterior_summary_table <- function(bundles, pars = NULL) {
  purrr::map_dfr(bundles, function(b) {
    draws <- as.data.frame(b$posterior_draws)
    keep <- pars %||% intersect(c("alpha", "beta", "sigma_model", "nu", "eta", "rho"), names(draws))
    purrr::map_dfr(keep, function(p) {
      x <- draws[[p]]
      tibble::tibble(
        Model = b$model_name,
        Type = b$model_type,
        Parameter = p,
        Mean = mean(x, na.rm = TRUE),
        Median = stats::median(x, na.rm = TRUE),
        SD = stats::sd(x, na.rm = TRUE),
        Q2.5 = stats::quantile(x, 0.025, na.rm = TRUE),
        Q97.5 = stats::quantile(x, 0.975, na.rm = TRUE)
      )
    })
  })
}
