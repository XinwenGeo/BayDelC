# Numerical validation of portable assets against the archived full bundles.
# Requires the project renv library because the full GPR fits require brms.
# This script reads models and predicts from existing posterior draws; it does
# not sample, update, or refit any model.

source_root <- Sys.getenv(
  "BAYDELC_SOURCE_ROOT",
  unset = file.path("..", "BayDelC_v0.3_publishable")
)
full_dir <- file.path(normalizePath(source_root), "inst", "models", "main")
portable_dir <- file.path("inst", "models", "main")

checks <- list()

for (calibration in c("Sub", "All")) {
  full_name <- paste0(calibration, "_Linear_EIV_Gaussian_v0.3.rds")
  portable_name <- paste0(calibration, "_Linear_EIV_Gaussian_v0.3_portable.rds")
  full <- readRDS(file.path(full_dir, full_name))
  portable <- readRDS(file.path(portable_dir, portable_name))
  columns <- names(portable$posterior_draws)
  delta <- max(abs(
    as.matrix(as.data.frame(full$posterior_draws)[columns]) -
      as.matrix(portable$posterior_draws)
  ))
  checks[[paste0("BLR_", calibration)]] <- data.frame(
    model = paste0("BLR_", calibration), comparison = "posterior parameters",
    max_abs_difference = delta, full_repeat_difference = 0, tolerance = 0
  )
}

if (!requireNamespace("brms", quietly = TRUE)) {
  stop("Validation of archived full GPR fits requires brms.")
}

validation_grid <- c(60, 80, 120, 160, 200, 240, 280)
interpolation_grid <- c(60.5, 119.5, 199.5, 279.5)

for (calibration in c("Sub", "All")) {
  full_name <- paste0(calibration, "_GPR_Baseline_v0.3.rds")
  portable_name <- paste0(calibration, "_GPR_Baseline_v0.3_portable.rds")
  full <- readRDS(file.path(full_dir, full_name))
  portable <- readRDS(file.path(portable_dir, portable_name))
  tr <- full$settings$x_transform

  predict_full <- function(x, seed) {
    nd <- data.frame(
      O2 = x,
      O2_std = (x - as.numeric(tr$mean)) / as.numeric(tr$sd),
      dd13C_sd = 0
    )
    set.seed(seed)
    as.matrix(brms::posterior_epred(
      full$fit_object, newdata = nd,
      draw_ids = portable$source_draw_ids, re_formula = NA
    ))
  }
  predict_portable <- function(x) {
    t(vapply(seq_len(nrow(portable$portable_mean_draws)), function(i) {
      stats::approx(portable$portable_bwo_grid,
                    portable$portable_mean_draws[i, ], xout = x,
                    ties = "ordered")$y
    }, numeric(length(x))))
  }
  summarize_draws <- function(z) {
    rbind(
      mean = colMeans(z),
      q025 = apply(z, 2L, stats::quantile, 0.025),
      median = apply(z, 2L, stats::median),
      q975 = apply(z, 2L, stats::quantile, 0.975)
    )
  }
  compare_distribution <- function(x) {
    full_a <- summarize_draws(predict_full(x, 901L))
    full_b <- summarize_draws(predict_full(x, 902L))
    portable_summary <- summarize_draws(predict_portable(x))
    repeat_difference <- max(abs(full_a - full_b))
    difference <- max(abs(full_a - portable_summary))
    # New-point GP predictions conditionally draw latent function values, so
    # even two calls using identical posterior IDs are not drawwise identical.
    # Require portable summaries to agree within a small absolute limit or
    # twice the observed repeat-prediction Monte Carlo variation.
    tolerance <- max(0.005, 2 * repeat_difference)
    c(difference = difference, repeat_difference = repeat_difference,
      tolerance = tolerance)
  }
  exact <- compare_distribution(validation_grid)
  interpolated <- compare_distribution(interpolation_grid)
  checks[[paste0("GPR_", calibration, "_grid")]] <- data.frame(
    model = paste0("GPR_", calibration), comparison = "stored grid points",
    max_abs_difference = exact[["difference"]],
    full_repeat_difference = exact[["repeat_difference"]],
    tolerance = exact[["tolerance"]]
  )
  checks[[paste0("GPR_", calibration, "_interpolation")]] <- data.frame(
    model = paste0("GPR_", calibration), comparison = "half-step interpolation",
    max_abs_difference = interpolated[["difference"]],
    full_repeat_difference = interpolated[["repeat_difference"]],
    tolerance = interpolated[["tolerance"]]
  )
}

results <- do.call(rbind, checks)
results$passed <- results$max_abs_difference <= results$tolerance
print(results, row.names = FALSE)
if (!all(results$passed)) stop("At least one portable-model validation failed.")
utils::write.csv(results, "PORTABLE_MODEL_VALIDATION.csv", row.names = FALSE)
