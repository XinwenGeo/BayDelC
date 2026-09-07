# BayDelC utility helpers ---------------------------------------------------

baydelc_version <- function() "0.3.0-calibration"

ensure_dirs <- function(paths = c(
  "data/processed",
  "inst/models/main",
  "inst/models/sensitivity",
  "outputs/tables",
  "outputs/figures/main",
  "outputs/figures/supplement",
  "outputs/logs"
)) {
  invisible(lapply(paths, dir.create, recursive = TRUE, showWarnings = FALSE))
}

safe_loo_from_loglik <- function(log_lik) {
  tryCatch(
    loo::loo(log_lik),
    error = function(e) {
      warning("LOO failed: ", conditionMessage(e))
      NULL
    }
  )
}

safe_loo_brms <- function(fit, moment_match = TRUE) {
  out <- tryCatch(
    brms::loo(fit, moment_match = moment_match),
    error = function(e) {
      warning("brms::loo(moment_match = ", moment_match, ") failed: ", conditionMessage(e),
              "; retrying with moment_match = FALSE.")
      tryCatch(brms::loo(fit, moment_match = FALSE), error = function(e2) {
        warning("brms::loo failed: ", conditionMessage(e2)); NULL
      })
    }
  )
  out
}

extract_pareto_k <- function(loo_obj) {
  if (is.null(loo_obj)) return(NA_real_)
  pk <- tryCatch(loo_obj$diagnostics$pareto_k, error = function(e) NA_real_)
  suppressWarnings(max(pk, na.rm = TRUE))
}

extract_elpd <- function(loo_obj, field = c("Estimate", "SE")) {
  field <- match.arg(field)
  if (is.null(loo_obj)) return(NA_real_)
  est <- loo_obj$estimates
  if (!"elpd_loo" %in% rownames(est)) return(NA_real_)
  est["elpd_loo", field]
}

posterior_coverage <- function(y, y_rep, prob = 0.95) {
  alpha <- (1 - prob) / 2
  lo <- apply(y_rep, 2, stats::quantile, alpha, na.rm = TRUE)
  hi <- apply(y_rep, 2, stats::quantile, 1 - alpha, na.rm = TRUE)
  mean(y >= lo & y <= hi, na.rm = TRUE)
}

rmse <- function(obs, pred) sqrt(mean((obs - pred)^2, na.rm = TRUE))

bayes_r2_simple <- function(yhat, obs) {
  var_yhat <- stats::var(yhat, na.rm = TRUE)
  var_res <- stats::var(obs - yhat, na.rm = TRUE)
  var_yhat / (var_yhat + var_res)
}

sample_draws_df <- function(draws, n_draw = 2000, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(draws)
  draws[sample(seq_len(n), min(n_draw, n), replace = FALSE), , drop = FALSE]
}

write_version_file <- function(path = "VERSION.json") {
  ensure_dirs()
  v <- list(
    version = baydelc_version(),
    data_version = "O2DD13C_NEW_v0.2",
    created_time = as.character(Sys.time()),
    r_version = as.character(getRversion())
  )
  jsonlite::write_json(v, path, auto_unbox = TRUE, pretty = TRUE)
}
