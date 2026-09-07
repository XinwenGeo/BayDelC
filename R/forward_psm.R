# BayDelC v0.3 unified forward interface -----------------------------------
# The user-facing package loads lightweight calibrated assets; no model is
# retrained by this file.

.psm_issue <- function(message, validate = c("strict", "warn", "none"), fatal = TRUE) {
  validate <- match.arg(validate)
  if (validate == "strict" && fatal) stop(message, call. = FALSE)
  if (validate == "warn") warning(message, call. = FALSE)
  invisible(FALSE)
}

.psm_recycle_numeric <- function(x, n, name, allow_null = FALSE) {
  if (is.null(x)) {
    if (allow_null) return(NULL)
    stop(sprintf("`%s` is required.", name), call. = FALSE)
  }
  x <- suppressWarnings(as.numeric(x))
  if (!length(x)) stop(sprintf("`%s` is empty.", name), call. = FALSE)
  if (!(length(x) %in% c(1L, n))) {
    stop(sprintf("`%s` must have length 1 or %d.", name, n), call. = FALSE)
  }
  rep(x, length.out = n)
}

.psm_validate_bounds <- function(bounds, name = "bounds") {
  if (is.null(bounds) || length(bounds) != 2L || any(is.na(bounds)) || bounds[1] >= bounds[2]) {
    stop(sprintf("`%s` must be an increasing numeric vector of length 2.", name), call. = FALSE)
  }
  as.numeric(bounds)
}

.psm_detect_column <- function(data, requested, candidates, required = FALSE, label = "column") {
  if (!is.null(requested)) {
    if (!requested %in% names(data)) stop(sprintf("%s `%s` was not found.", label, requested), call. = FALSE)
    return(requested)
  }
  hit <- candidates[candidates %in% names(data)]
  if (length(hit)) return(hit[[1L]])
  if (required) stop(sprintf("Could not detect the %s. Supply its column name explicitly.", label), call. = FALSE)
  NULL
}

.psm_site_labels <- function(data, site_col, n) {
  if (!is.null(data)) {
    sc <- .psm_detect_column(data, site_col, c("Site", "site", "Core", "core", "Core name"), FALSE, "site column")
    if (!is.null(sc)) return(as.character(data[[sc]]))
  }
  as.character(seq_len(n))
}

.psm_as_draw_matrix <- function(x, n_obs, name, validate = "strict") {
  if (is.null(x)) return(NULL)
  if (is.vector(x) && !is.list(x)) {
    if (n_obs != 1L) stop(sprintf("Vector `%s` is only valid for one observation.", name), call. = FALSE)
    mat <- matrix(as.numeric(x), nrow = 1L)
  } else {
    mat <- as.matrix(x)
    storage.mode(mat) <- "double"
  }
  if (nrow(mat) != n_obs) {
    if (ncol(mat) == n_obs && validate != "strict") {
      .psm_issue(sprintf("`%s` was transposed to enforce rows = observations and columns = realizations.", name), validate, fatal = FALSE)
      mat <- t(mat)
    } else {
      stop(sprintf("`%s` must have %d rows (one row per observation).", name, n_obs), call. = FALSE)
    }
  }
  if (!ncol(mat)) stop(sprintf("`%s` contains no realizations.", name), call. = FALSE)
  if (any(!is.finite(mat))) stop(sprintf("`%s` contains non-finite values.", name), call. = FALSE)
  mat
}

.psm_normalize_weights <- function(w, n, name = "weights") {
  if (is.null(w)) return(rep(1 / n, n))
  w <- as.numeric(w)
  if (length(w) != n || any(!is.finite(w)) || any(w < 0) || sum(w) <= 0) {
    stop(sprintf("`%s` must contain %d finite non-negative values with a positive sum.", name, n), call. = FALSE)
  }
  w / sum(w)
}

.psm_truncated_normal <- function(n, mean, sd, lower = 0, upper = Inf) {
  if (!is.finite(mean) || !is.finite(sd) || sd < 0) stop("Invalid normal input parameters.", call. = FALSE)
  if (sd == 0) {
    if (mean < lower || mean > upper) stop("A fixed input lies outside `input_bounds`.", call. = FALSE)
    return(rep(mean, n))
  }
  pl <- if (is.infinite(lower) && lower < 0) 0 else stats::pnorm((lower - mean) / sd)
  pu <- if (is.infinite(upper) && upper > 0) 1 else stats::pnorm((upper - mean) / sd)
  if (!is.finite(pl) || !is.finite(pu) || pu <= pl) stop("The truncated normal has negligible probability inside the requested bounds.", call. = FALSE)
  u <- stats::runif(n, min = max(pl, .Machine$double.eps), max = min(pu, 1 - .Machine$double.eps))
  mean + sd * stats::qnorm(u)
}

.psm_model_registry <- function(algorithm, calibration) {
  if (calibration == "custom") return(NULL)
  .baydelc_model_path(algorithm, calibration, must_exist = FALSE)
}

.psm_resolve_model <- function(model = NULL, algorithm = "BLR", calibration = "Sub") {
  algorithm <- match.arg(toupper(algorithm), c("BLR", "GPR"))
  calibration <- match.arg(calibration, c("Sub", "All", "custom"))
  if (is.null(model)) {
    if (calibration == "custom") stop("`calibration = 'custom'` requires `model`.", call. = FALSE)
    path <- .psm_model_registry(algorithm, calibration)
    if (is.null(path) || !file.exists(path)) stop(sprintf("Model file not found: %s", path %||% "<unknown>"), call. = FALSE)
    model <- readRDS(path)
    attr(model, "psm_model_path") <- path
  } else if (is.character(model) && length(model) == 1L) {
    if (!file.exists(model)) stop(sprintf("Model file not found: %s", model), call. = FALSE)
    path <- model
    model <- readRDS(model)
    attr(model, "psm_model_path") <- path
  }
  if (!is.list(model)) stop("`model` must be a v0.3 bundle, custom model list, or RDS path.", call. = FALSE)
  inferred <- if (grepl("GPR", paste(model$model_name %||% "", model$model_type %||% ""), ignore.case = TRUE)) "GPR" else "BLR"
  if (!is.null(model$algorithm)) inferred <- toupper(model$algorithm)
  if (algorithm != inferred && calibration != "custom") {
    stop(sprintf("Requested algorithm `%s`, but the supplied model appears to be `%s`.", algorithm, inferred), call. = FALSE)
  }
  list(bundle = model, algorithm = inferred, calibration = calibration)
}

.psm_calibration_range <- function(model) {
  r <- model$calibration_range
  if (is.numeric(r) && length(r) >= 2L && all(is.finite(r[1:2]))) return(range(r[1:2]))
  if (!is.null(model$training_data) && "O2" %in% names(model$training_data)) {
    r <- range(as.numeric(model$training_data$O2), na.rm = TRUE)
    if (all(is.finite(r))) return(r)
  }
  c(50, 300)
}

.psm_sample_rows <- function(n_available, n_draw, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  sample.int(n_available, size = n_draw, replace = n_draw > n_available)
}

.psm_mvn2 <- function(n, mean, covariance) {
  covariance <- as.matrix(covariance)
  if (!all(dim(covariance) == c(2L, 2L))) stop("The alpha/beta covariance matrix must be 2 x 2.", call. = FALSE)
  ev <- eigen(covariance, symmetric = TRUE, only.values = TRUE)$values
  if (any(ev < -1e-10)) stop("The alpha/beta covariance matrix is not positive semidefinite.", call. = FALSE)
  L <- chol(covariance + diag(1e-12, 2L))
  z <- matrix(stats::rnorm(n * 2L), nrow = n, ncol = 2L)
  sweep(z %*% L, 2L, mean, "+")
}

.psm_blr_parameter_draws <- function(model, n_draw, seed = 123) {
  dr <- model$posterior_draws
  if (!is.null(dr)) {
    dr <- as.data.frame(dr)
    alpha_name <- intersect(c("alpha", "b_Intercept", "Intercept"), names(dr))
    beta_name <- intersect(c("beta", "b_O2", "b_BWO"), names(dr))
    sigma_name <- intersect(c("sigma_model", "sigma", "residual_sd"), names(dr))
    nu_name <- intersect(c("nu", "student_df"), names(dr))
    if (!length(alpha_name) || !length(beta_name)) stop("The BLR posterior draws lack `alpha` and/or `beta`.", call. = FALSE)
    idx <- .psm_sample_rows(nrow(dr), n_draw, seed)
    out <- data.frame(
      alpha = as.numeric(dr[[alpha_name[[1L]]]][idx]),
      beta = as.numeric(dr[[beta_name[[1L]]]][idx]),
      sigma_model = if (length(sigma_name)) abs(as.numeric(dr[[sigma_name[[1L]]]][idx])) else 0,
      nu = if (length(nu_name)) as.numeric(dr[[nu_name[[1L]]]][idx]) else NA_real_
    )
    return(out)
  }
  if (is.null(model$alpha) || is.null(model$beta)) stop("A custom BLR model requires posterior draws or `alpha` and `beta`.", call. = FALSE)
  set.seed(seed)
  alpha <- as.numeric(model$alpha); beta <- as.numeric(model$beta)
  if (length(alpha) > 1L || length(beta) > 1L) {
    m <- max(length(alpha), length(beta))
    raw <- data.frame(alpha = rep(alpha, length.out = m), beta = rep(beta, length.out = m))
    idx <- .psm_sample_rows(m, n_draw, seed)
    a <- raw$alpha[idx]; b <- raw$beta[idx]
  } else if (!is.null(model$alpha_sd) || !is.null(model$beta_sd) || !is.null(model$alpha_beta_cov)) {
    a_sd <- as.numeric(model$alpha_sd %||% 0); b_sd <- as.numeric(model$beta_sd %||% 0)
    cov_ab <- as.numeric(model$alpha_beta_cov %||% 0)
    ab <- .psm_mvn2(n_draw, c(alpha, beta), matrix(c(a_sd^2, cov_ab, cov_ab, b_sd^2), 2L, 2L))
    a <- ab[, 1L]; b <- ab[, 2L]
  } else {
    a <- rep(alpha, n_draw); b <- rep(beta, n_draw)
  }
  sigma <- rep(abs(as.numeric(model$sigma %||% model$sigma_model %||% 0)), length.out = n_draw)
  nu <- rep(as.numeric(model$nu %||% NA_real_), length.out = n_draw)
  data.frame(alpha = a, beta = b, sigma_model = sigma, nu = nu)
}

.psm_is_student_t <- function(model, pars = NULL) {
  txt <- paste(model$model_name %||% "", model$model_type %||% "", model$residual_distribution %||% "")
  grepl("STUDENT|student|_T", txt) || (!is.null(pars) && any(is.finite(pars$nu)))
}

.psm_gpr_transform <- function(model) {
  tr <- model$settings$x_transform %||% list()
  m <- tr$mean %||% if (!is.null(model$training_data$O2_mean_train)) unique(model$training_data$O2_mean_train)[1L] else mean(model$training_data$O2, na.rm = TRUE)
  s <- tr$sd %||% if (!is.null(model$training_data$O2_sd_train_scale)) unique(model$training_data$O2_sd_train_scale)[1L] else stats::sd(model$training_data$O2, na.rm = TRUE)
  m <- as.numeric(m[1L]); s <- as.numeric(s[1L])
  if (!is.finite(m) || !is.finite(s) || s <= 0) stop("Could not recover the GPR O2 standardization parameters.", call. = FALSE)
  list(mean = m, sd = s)
}

.psm_gpr_available_draws <- function(model) {
  if (!is.null(model$portable_mean_draws)) return(nrow(as.matrix(model$portable_mean_draws)))
  if (!is.null(model$posterior_draws)) return(nrow(as.data.frame(model$posterior_draws)))
  if (!is.null(model$fit_object) && requireNamespace("posterior", quietly = TRUE)) return(posterior::ndraws(model$fit_object))
  stop("Could not determine the number of GPR posterior draws.", call. = FALSE)
}

.psm_is_portable_gpr <- function(model) {
  is.numeric(model$portable_bwo_grid) &&
    length(model$portable_bwo_grid) >= 2L &&
    !is.null(model$portable_mean_draws)
}

.psm_portable_gpr_mean <- function(model, x, draw_ids) {
  grid <- as.numeric(model$portable_bwo_grid)
  curves <- as.matrix(model$portable_mean_draws)
  if (ncol(curves) != length(grid)) {
    stop("Portable GPR grid and posterior curves have incompatible dimensions.", call. = FALSE)
  }
  if (any(!is.finite(x)) || any(x < min(grid) | x > max(grid))) {
    stop(
      sprintf(
        "Portable GPR predictions support BWO from %g to %g umol kg^-1.",
        min(grid), max(grid)
      ),
      call. = FALSE
    )
  }
  if (any(draw_ids < 1L | draw_ids > nrow(curves))) {
    stop("Portable GPR draw index is out of range.", call. = FALSE)
  }
  out <- matrix(NA_real_, nrow = length(draw_ids), ncol = length(x))
  for (k in seq_along(draw_ids)) {
    out[k, ] <- stats::approx(
      grid, curves[draw_ids[k], ], xout = x,
      method = "linear", ties = "ordered", rule = 1
    )$y
  }
  out
}

.psm_gpr_newdata <- function(model, x) {
  tr <- .psm_gpr_transform(model)
  data.frame(
    O2 = as.numeric(x),
    O2_std = (as.numeric(x) - tr$mean) / tr$sd,
    dd13C_sd = 0,
    dd13c_sd = 0
  )
}

.psm_gpr_sigma_draws <- function(model, draw_ids) {
  dr <- model$posterior_draws
  if (!is.null(dr)) {
    dr <- as.data.frame(dr)
    nm <- intersect(c("sigma", "sigma_model", "residual_sd"), names(dr))
    if (length(nm)) return(abs(as.numeric(dr[[nm[[1L]]]][draw_ids])))
  }
  rep(0, length(draw_ids))
}


.psm_is_latent_gpr <- function(model) {
  grepl("LATENT", paste(model$model_name %||% "", model$model_type %||% ""), ignore.case = TRUE) &&
    !is.null(model$fit_object) && inherits(model$fit_object, "stanfit")
}

.psm_latent_extract <- function(model) {
  if (!requireNamespace("rstan", quietly = TRUE)) stop("Package `rstan` is required for the latent-input GPR bundle.", call. = FALSE)
  ext <- rstan::extract(model$fit_object)
  needed <- c("alpha", "eta", "rho", "x_true", "f_latent")
  miss <- setdiff(needed, names(ext))
  if (length(miss)) stop(sprintf("Latent-input GPR is missing posterior arrays: %s", paste(miss, collapse = ", ")), call. = FALSE)
  ext
}

.psm_latent_mean_one <- function(ext, draw_id, x_std) {
  xt <- as.numeric(ext$x_true[draw_id, ])
  f <- as.numeric(ext$f_latent[draw_id, ])
  eta <- as.numeric(ext$eta[draw_id])
  rho <- as.numeric(ext$rho[draw_id])
  alpha <- as.numeric(ext$alpha[draw_id])
  K <- outer(xt, xt, function(a, b) eta^2 * exp(-0.5 * ((a - b) / rho)^2))
  diag(K) <- diag(K) + 1e-6
  Ks <- outer(x_std, xt, function(a, b) eta^2 * exp(-0.5 * ((a - b) / rho)^2))
  sol <- try(solve(K, f), silent = TRUE)
  if (inherits(sol, "try-error")) return(rep(NA_real_, length(x_std)))
  as.numeric(alpha + Ks %*% sol)
}

.psm_latent_mean_grid <- function(model, x, draw_ids) {
  ext <- .psm_latent_extract(model)
  tr <- .psm_gpr_transform(model)
  x_std <- (as.numeric(x) - tr$mean) / tr$sd
  out <- matrix(NA_real_, nrow = length(draw_ids), ncol = length(x))
  for (k in seq_along(draw_ids)) out[k, ] <- .psm_latent_mean_one(ext, draw_ids[k], x_std)
  out
}

.psm_latent_diagonal_predict <- function(model, x_draws, draw_ids, predictive = TRUE) {
  ext <- .psm_latent_extract(model)
  tr <- .psm_gpr_transform(model)
  n_obs <- nrow(x_draws); n_draw <- ncol(x_draws)
  out <- matrix(NA_real_, nrow = n_obs, ncol = n_draw)
  for (k in seq_len(n_draw)) {
    x_std <- (x_draws[, k] - tr$mean) / tr$sd
    out[, k] <- .psm_latent_mean_one(ext, draw_ids[k], x_std)
  }
  if (predictive) {
    sig <- .psm_gpr_sigma_draws(model, draw_ids)
    out <- out + matrix(stats::rnorm(n_obs * n_draw), nrow = n_obs) * rep(sig, each = n_obs)
  }
  out
}

.psm_gpr_mean_grid <- function(model, x, draw_ids) {
  if (.psm_is_portable_gpr(model)) return(.psm_portable_gpr_mean(model, x, draw_ids))
  if (.psm_is_latent_gpr(model)) return(.psm_latent_mean_grid(model, x, draw_ids))
  if (!is.null(model$predict_mean_draws) && is.function(model$predict_mean_draws)) {
    ans <- model$predict_mean_draws(model = model, x = x, draw_ids = draw_ids)
    ans <- as.matrix(ans)
    if (!all(dim(ans) == c(length(draw_ids), length(x)))) stop("Custom GPR `predict_mean_draws` returned incorrect dimensions.", call. = FALSE)
    return(ans)
  }
  if (is.null(model$fit_object) || !inherits(model$fit_object, "brmsfit")) stop("GPR prediction requires a brms `fit_object` or custom `predict_mean_draws` function.", call. = FALSE)
  if (!requireNamespace("brms", quietly = TRUE)) stop("Package `brms` is required for GPR prediction.", call. = FALSE)
  as.matrix(brms::posterior_epred(model$fit_object, newdata = .psm_gpr_newdata(model, x), draw_ids = draw_ids, re_formula = NA))
}

.psm_gpr_diagonal_predict <- function(model, x_draws, draw_ids, predictive = TRUE, chunk_size = 200L) {
  if (.psm_is_portable_gpr(model)) {
    if (length(draw_ids) != ncol(x_draws)) {
      stop("Portable GPR requires one posterior draw index per input draw.", call. = FALSE)
    }
    grid <- as.numeric(model$portable_bwo_grid)
    curves <- as.matrix(model$portable_mean_draws)
    if (any(x_draws < min(grid) | x_draws > max(grid))) {
      stop(
        sprintf(
          "Portable GPR predictions support BWO from %g to %g umol kg^-1.",
          min(grid), max(grid)
        ),
        call. = FALSE
      )
    }
    out <- matrix(NA_real_, nrow = nrow(x_draws), ncol = ncol(x_draws))
    for (k in seq_len(ncol(x_draws))) {
      out[, k] <- stats::approx(
        grid, curves[draw_ids[k], ], xout = x_draws[, k],
        method = "linear", ties = "ordered", rule = 1
      )$y
    }
    if (predictive) {
      sig <- .psm_gpr_sigma_draws(model, draw_ids)
      out <- out + matrix(stats::rnorm(length(out)), nrow = nrow(out)) *
        rep(sig, each = nrow(out))
    }
    return(out)
  }
  if (.psm_is_latent_gpr(model)) return(.psm_latent_diagonal_predict(model, x_draws, draw_ids, predictive))
  if (!is.null(model$predict_draws) && is.function(model$predict_draws)) {
    ans <- model$predict_draws(model = model, bwo_draws = x_draws, draw_ids = draw_ids, predictive = predictive)
    ans <- as.matrix(ans)
    if (!all(dim(ans) == dim(x_draws))) stop("Custom GPR `predict_draws` must return observations x draws.", call. = FALSE)
    return(ans)
  }
  if (is.null(model$fit_object) || !inherits(model$fit_object, "brmsfit")) stop("GPR prediction requires a brms `fit_object` or custom `predict_draws` function.", call. = FALSE)
  if (!requireNamespace("brms", quietly = TRUE)) stop("Package `brms` is required for GPR prediction.", call. = FALSE)
  n_obs <- nrow(x_draws); n_draw <- ncol(x_draws)
  out <- matrix(NA_real_, nrow = n_obs, ncol = n_draw)
  chunk_size <- max(1L, as.integer(chunk_size))
  for (i in seq_len(n_obs)) {
    starts <- seq.int(1L, n_draw, by = chunk_size)
    for (st in starts) {
      en <- min(n_draw, st + chunk_size - 1L)
      k <- st:en
      nd <- .psm_gpr_newdata(model, x_draws[i, k])
      mat <- if (predictive) {
        brms::posterior_predict(model$fit_object, newdata = nd, draw_ids = draw_ids[k], re_formula = NA)
      } else {
        brms::posterior_epred(model$fit_object, newdata = nd, draw_ids = draw_ids[k], re_formula = NA)
      }
      mat <- as.matrix(mat)
      out[i, k] <- diag(mat)
    }
  }
  out
}

.psm_summary_matrix <- function(draws, sites, input_value, input_sd = NULL) {
  data.frame(
    row_id = seq_len(nrow(draws)),
    Site = sites,
    input_value = as.numeric(input_value),
    input_sd = if (is.null(input_sd)) NA_real_ else as.numeric(input_sd),
    mean = rowMeans(draws, na.rm = TRUE),
    median = apply(draws, 1L, stats::median, na.rm = TRUE),
    q025 = apply(draws, 1L, stats::quantile, probs = 0.025, na.rm = TRUE),
    q975 = apply(draws, 1L, stats::quantile, probs = 0.975, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

#' Propagate bottom-water oxygen through a BayDelC calibration
#'
#' `forward_psm()` estimates the distribution of benthic foraminiferal
#' Delta-delta-13-C expected from one or more bottom-water oxygen (BWO) values.
#' It propagates calibration uncertainty, optional BWO input uncertainty, and
#' optionally residual model variability.
#'
#' @param data Optional data frame containing BWO values and, optionally,
#'   uncertainty and site columns.
#' @param bwo Numeric BWO values in micromoles per kilogram.
#' @param bwo_sd Optional one-standard-deviation uncertainty for `bwo`.
#' @param bwo_draws Optional matrix of BWO draws with observations in rows.
#' @param bwo_draw_weights Optional weights for the columns of `bwo_draws`.
#' @param bwo_col,bwo_sd_col,site_col Column names used when `data` is supplied.
#' @param input_error_method Input uncertainty representation. `"auto"`
#'   selects from the supplied arguments.
#' @param input_distribution Optional custom random generator for BWO inputs.
#' @param input_bounds Two-element physical bounds for sampled BWO inputs.
#' @param algorithm Calibration family: `"BLR"` or `"GPR"`.
#' @param calibration Calibration dataset: restricted-taxon `"Sub"`,
#'   mixed-taxon `"All"`, or `"custom"`.
#' @param model Optional custom model list or path to an RDS model bundle.
#' @param prediction `"predictive"` includes residual variability;
#'   `"latent"` returns the latent calibration relationship.
#' @param include_model_residual Logical override controlling residual
#'   variability.
#' @param n_draw Number of Monte Carlo draws.
#' @param seed Random seed.
#' @param return_draws Return the observation-by-draw matrix when `TRUE`.
#' @param validate Validation policy: `"strict"`, `"warn"`, or `"none"`.
#' @param gpr_chunk_size Chunk size used only with a custom full `brms` model.
#'
#' @return An object of class `forward_psm_result` containing `summary`,
#'   optional `draws`, diagnostics, and model metadata.
#' @export
#' @examples
#' forward_psm(bwo = c(80, 150, 220), algorithm = "BLR", n_draw = 500)
#' forward_psm(bwo = 120, bwo_sd = 15, algorithm = "GPR", n_draw = 500)
forward_psm <- function(
    data = NULL,
    bwo = NULL,
    bwo_sd = NULL,
    bwo_draws = NULL,
    bwo_draw_weights = NULL,
    bwo_col = "O2",
    bwo_sd_col = NULL,
    site_col = "Site",
    input_error_method = c("auto", "fixed", "normal", "draws", "custom"),
    input_distribution = NULL,
    input_bounds = c(0, Inf),
    algorithm = c("BLR", "GPR"),
    calibration = c("Sub", "All", "custom"),
    model = NULL,
    prediction = c("predictive", "latent"),
    include_model_residual = TRUE,
    n_draw = 2000,
    seed = 123,
    return_draws = FALSE,
    validate = c("strict", "warn", "none"),
    gpr_chunk_size = 200L) {

  input_error_method <- match.arg(input_error_method)
  algorithm <- match.arg(algorithm)
  calibration <- match.arg(calibration)
  prediction <- match.arg(prediction)
  validate <- match.arg(validate)
  input_bounds <- .psm_validate_bounds(input_bounds, "input_bounds")
  n_draw <- as.integer(n_draw)
  if (!is.finite(n_draw) || n_draw < 1L) stop("`n_draw` must be a positive integer.", call. = FALSE)

  if (!is.null(data)) {
    data <- as.data.frame(data)
    if (!is.null(bwo)) .psm_issue("Supply BWO through either `data` or `bwo`, not both.", validate)
    bc <- .psm_detect_column(data, bwo_col, c("O2", "BWO", "bwo", "Bottom_Water_Oxygen"), TRUE, "BWO column")
    bwo <- as.numeric(data[[bc]])
    if (!is.null(bwo_sd_col)) {
      sc <- .psm_detect_column(data, bwo_sd_col, c("O2_sd", "BWO_sd", "bwo_sd"), TRUE, "BWO SD column")
      if (!is.null(bwo_sd)) .psm_issue("Supply BWO SD through either `data` or `bwo_sd`, not both.", validate)
      bwo_sd <- as.numeric(data[[sc]])
    } else if (input_error_method == "auto" && is.null(bwo_sd)) {
      auto_sc <- .psm_detect_column(data, NULL, c("O2_sd", "BWO_sd", "bwo_sd"), FALSE, "BWO SD column")
      if (!is.null(auto_sc)) bwo_sd <- as.numeric(data[[auto_sc]])
    }
  }
  bwo <- as.numeric(bwo)
  if (!length(bwo) || any(!is.finite(bwo))) stop("`bwo` must contain finite values.", call. = FALSE)
  n_obs <- length(bwo)
  sites <- .psm_site_labels(data, site_col, n_obs)
  if (length(sites) != n_obs) stop("The site column length does not match the BWO input.", call. = FALSE)

  if (any(bwo < input_bounds[1] | bwo > input_bounds[2])) {
    .psm_issue("At least one central BWO input lies outside `input_bounds`.", validate)
  }
  if (!is.null(bwo_sd)) {
    bwo_sd <- .psm_recycle_numeric(bwo_sd, n_obs, "bwo_sd")
    if (any(!is.finite(bwo_sd) | bwo_sd < 0)) stop("`bwo_sd` must be finite and non-negative.", call. = FALSE)
  }
  draw_mat <- .psm_as_draw_matrix(bwo_draws, n_obs, "bwo_draws", validate)
  representations <- c(sd = !is.null(bwo_sd), draws = !is.null(draw_mat), custom = !is.null(input_distribution))
  if (input_error_method == "auto") {
    if (sum(representations) > 1L) .psm_issue("Conflicting BWO uncertainty representations were supplied.", validate)
    input_error_method <- if (!is.null(draw_mat)) "draws" else if (!is.null(bwo_sd)) "normal" else if (!is.null(input_distribution)) "custom" else "fixed"
  }
  if (input_error_method == "draws" && is.null(draw_mat)) stop("`bwo_draws` is required for `input_error_method = 'draws'`.", call. = FALSE)
  if (input_error_method == "normal" && is.null(bwo_sd)) stop("`bwo_sd` is required for `input_error_method = 'normal'`.", call. = FALSE)
  if (input_error_method == "custom" && !is.function(input_distribution)) stop("`input_distribution` must be a function for the custom method.", call. = FALSE)
  if (input_error_method == "fixed" && any(representations)) .psm_issue("Uncertainty inputs are incompatible with `input_error_method = 'fixed'`.", validate)

  set.seed(seed)
  x_draws <- switch(input_error_method,
    fixed = matrix(rep(bwo, n_draw), nrow = n_obs),
    normal = {
      m <- matrix(NA_real_, nrow = n_obs, ncol = n_draw)
      for (i in seq_len(n_obs)) m[i, ] <- .psm_truncated_normal(n_draw, bwo[i], bwo_sd[i], input_bounds[1], input_bounds[2])
      m
    },
    draws = {
      if (any(draw_mat < input_bounds[1] | draw_mat > input_bounds[2])) .psm_issue("`bwo_draws` contains values outside `input_bounds`.", validate)
      w <- .psm_normalize_weights(bwo_draw_weights, ncol(draw_mat), "bwo_draw_weights")
      idx <- sample.int(ncol(draw_mat), n_draw, replace = n_draw > ncol(draw_mat), prob = w)
      draw_mat[, idx, drop = FALSE]
    },
    custom = {
      m <- matrix(NA_real_, nrow = n_obs, ncol = n_draw)
      for (i in seq_len(n_obs)) {
        ans <- try(do.call(input_distribution, list(n = n_draw, mean = bwo[i], sd = (bwo_sd %||% rep(0, n_obs))[i], lower = input_bounds[1], upper = input_bounds[2], row = i, data = data)), silent = TRUE)
        if (inherits(ans, "try-error")) ans <- input_distribution(n_draw, bwo[i], (bwo_sd %||% rep(0, n_obs))[i])
        ans <- as.numeric(ans)
        if (length(ans) != n_draw || any(!is.finite(ans))) stop("`input_distribution` must return `n_draw` finite values per observation.", call. = FALSE)
        if (any(ans < input_bounds[1] | ans > input_bounds[2])) .psm_issue("Custom BWO draws fall outside `input_bounds`.", validate)
        m[i, ] <- ans
      }
      m
    }
  )

  resolved <- .psm_resolve_model(model, algorithm, calibration)
  bundle <- resolved$bundle
  algorithm_used <- resolved$algorithm
  include_residual_effective <- isTRUE(include_model_residual) && prediction == "predictive"

  if (algorithm_used == "BLR") {
    pars <- .psm_blr_parameter_draws(bundle, n_draw, seed + 1L)
    y_draws <- sweep(x_draws, 2L, pars$beta, "*")
    y_draws <- sweep(y_draws, 2L, pars$alpha, "+")
    if (include_residual_effective) {
      if (.psm_is_student_t(bundle, pars)) {
        nu <- ifelse(is.finite(pars$nu) & pars$nu > 2, pars$nu, 5)
        eps <- matrix(stats::rt(n_obs * n_draw, df = rep(nu, each = n_obs)), nrow = n_obs)
        eps <- sweep(eps, 2L, pars$sigma_model, "*")
      } else {
        eps <- matrix(stats::rnorm(n_obs * n_draw), nrow = n_obs)
        eps <- sweep(eps, 2L, pars$sigma_model, "*")
      }
      y_draws <- y_draws + eps
    }
  } else {
    n_avail <- .psm_gpr_available_draws(bundle)
    draw_ids <- .psm_sample_rows(n_avail, n_draw, seed + 1L)
    y_draws <- .psm_gpr_diagonal_predict(bundle, x_draws, draw_ids, include_residual_effective, gpr_chunk_size)
  }

  calib_range <- .psm_calibration_range(bundle)
  summary <- .psm_summary_matrix(y_draws, sites, bwo, bwo_sd)
  summary$BWO <- bwo
  summary$BWO_sd <- if (is.null(bwo_sd)) NA_real_ else bwo_sd
  summary$P_BWO_below_50 <- rowMeans(x_draws < 50)
  summary$P_BWO_outside_calibration <- rowMeans(x_draws < calib_range[1] | x_draws > calib_range[2])
  summary$algorithm <- algorithm_used
  summary$calibration <- calibration
  summary$prediction <- prediction
  summary$input_error_method <- input_error_method

  input_out <- if (is.null(data)) data.frame(row_id = seq_len(n_obs), Site = sites, BWO = bwo, stringsAsFactors = FALSE) else data
  diagnostics <- summary[, c("row_id", "Site", "P_BWO_below_50", "P_BWO_outside_calibration"), drop = FALSE]
  result <- list(
    input = input_out,
    summary = summary,
    draws = if (return_draws) y_draws else NULL,
    density = NULL,
    metadata = list(
      direction = "forward",
      model_name = bundle$model_name %||% "custom",
      model_type = bundle$model_type %||% algorithm_used,
      model_path = attr(bundle, "psm_model_path"),
      algorithm = algorithm_used,
      calibration = calibration,
      calibration_range = calib_range,
      prediction_requested = prediction,
      include_model_residual = include_residual_effective,
      input_error_method = input_error_method,
      input_bounds = input_bounds,
      draw_orientation = "rows = observations; columns = joint realizations",
      n_draw = n_draw,
      seed = seed
    ),
    diagnostics = diagnostics
  )
  class(result) <- c("forward_psm_result", "psm_result", "list")
  result
}
