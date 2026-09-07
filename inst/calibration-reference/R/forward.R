# BayDelC forward model -----------------------------------------------------
# Forward maps BWO/O2 to Δδ13C distributions.

forward_blr <- function(bundle, O2_grid, O2_sd = 0,
                        mode = c("strict", "practical", "deterministic"),
                        include_residual = TRUE, n_draw = 2000, seed = 9001,
                        return_draws = FALSE) {
  mode <- match.arg(mode)
  stopifnot(grepl("BLR_EIV", bundle$model_type))
  draws <- as.data.frame(bundle$posterior_draws)
  draws <- sample_draws_df(draws, n_draw, seed)

  if (mode %in% c("practical", "deterministic")) {
    draws <- data.frame(
      alpha = stats::median(draws$alpha, na.rm = TRUE),
      beta = stats::median(draws$beta, na.rm = TRUE),
      sigma_model = stats::median(draws$sigma_model, na.rm = TRUE)
    )
  }

  O2_grid <- as.numeric(O2_grid)
  O2_sd <- rep(O2_sd, length.out = length(O2_grid))
  nd <- nrow(draws)
  ng <- length(O2_grid)
  pred <- matrix(NA_real_, nrow = nd, ncol = ng)

  for (j in seq_len(nd)) {
    x_eff <- if (mode == "deterministic") O2_grid else stats::rnorm(ng, O2_grid, O2_sd)
    mu <- draws$alpha[j] + draws$beta[j] * x_eff
    if (include_residual && mode != "deterministic") {
      pred[j, ] <- stats::rnorm(ng, mu, draws$sigma_model[j])
    } else {
      pred[j, ] <- mu
    }
  }
  if (mode == "deterministic") pred <- matrix(pred[1, ], nrow = 1)

  out <- tibble::tibble(
    O2 = O2_grid,
    median = apply(pred, 2, stats::median, na.rm = TRUE),
    q025 = apply(pred, 2, stats::quantile, 0.025, na.rm = TRUE),
    q975 = apply(pred, 2, stats::quantile, 0.975, na.rm = TRUE)
  )
  if (return_draws) attr(out, "draws") <- pred
  out
}

forward_gpr <- function(bundle, O2_grid, mode = c("predict", "epred"),
                        ndraws = NULL, return_draws = FALSE) {
  mode <- match.arg(mode)
  if (!inherits(bundle$fit_object, "brmsfit")) stop("forward_gpr requires a brms GPR bundle.")
  tr <- bundle$settings$x_transform
  newdata <- tibble::tibble(
    O2 = as.numeric(O2_grid),
    O2_std = (as.numeric(O2_grid) - tr$mean) / tr$sd,
    dd13C_sd = 0
  )
  if (mode == "predict") {
    pred <- brms::posterior_predict(bundle$fit_object, newdata = newdata, ndraws = ndraws)
  } else {
    pred <- brms::posterior_epred(bundle$fit_object, newdata = newdata, ndraws = ndraws)
  }
  out <- tibble::tibble(
    O2 = as.numeric(O2_grid),
    median = apply(pred, 2, stats::median, na.rm = TRUE),
    q025 = apply(pred, 2, stats::quantile, 0.025, na.rm = TRUE),
    q975 = apply(pred, 2, stats::quantile, 0.975, na.rm = TRUE)
  )
  if (return_draws) attr(out, "draws") <- pred
  out
}

forward_baydelc <- function(bundle, O2_grid, ...) {
  if (grepl("BLR_EIV", bundle$model_type)) {
    forward_blr(bundle, O2_grid, ...)
  } else if (grepl("GPR", bundle$model_type) && inherits(bundle$fit_object, "brmsfit")) {
    forward_gpr(bundle, O2_grid, ...)
  } else {
    stop("Unsupported bundle type for forward prediction: ", bundle$model_type)
  }
}
