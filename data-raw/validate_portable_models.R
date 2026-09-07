# Validate portable main and sensitivity assets against archived full bundles.
# This script predicts from existing posterior draws and never refits a model.

source_root <- normalizePath(Sys.getenv(
  "BAYDELC_SOURCE_ROOT", unset = file.path("..", "BayDelC_v0.3_publishable")
))
full_root <- file.path(source_root, "inst", "models")
portable_root <- file.path("inst", "models")
checks <- list()

add_check <- function(id, comparison, difference, repeat_difference = 0,
                      tolerance = 0) {
  checks[[length(checks) + 1L]] <<- data.frame(
    model = id, comparison = comparison,
    max_abs_difference = difference,
    full_repeat_difference = repeat_difference,
    tolerance = tolerance,
    stringsAsFactors = FALSE
  )
}

blr_jobs <- data.frame(
  id = c("BLR_Sub_main", "BLR_All_main", "BLR_Sub_StudentT", "BLR_All_StudentT"),
  set = c("main", "main", "sensitivity", "sensitivity"),
  source = c("Sub_Linear_EIV_Gaussian_v0.3.rds", "All_Linear_EIV_Gaussian_v0.3.rds",
             "Sub_Linear_EIV_StudentT_v0.3.rds", "All_Linear_EIV_StudentT_v0.3.rds"),
  portable = c("Sub_Linear_EIV_Gaussian_v0.3_portable.rds", "All_Linear_EIV_Gaussian_v0.3_portable.rds",
               "Sub_Linear_EIV_StudentT_v0.3_portable.rds", "All_Linear_EIV_StudentT_v0.3_portable.rds")
)
for (i in seq_len(nrow(blr_jobs))) {
  job <- blr_jobs[i, ]
  full <- readRDS(file.path(full_root, job$set, job$source))
  portable <- readRDS(file.path(portable_root, job$set, job$portable))
  columns <- names(portable$posterior_draws)
  difference <- max(abs(as.matrix(as.data.frame(full$posterior_draws)[columns]) -
                        as.matrix(portable$posterior_draws)))
  add_check(job$id, "posterior parameters", difference)
}

if (!requireNamespace("brms", quietly = TRUE)) stop("Full GPR validation requires brms.")
validation_grid <- c(60, 80, 120, 160, 200, 240, 280)
interpolation_grid <- c(60.5, 119.5, 199.5, 279.5)
summarize_draws <- function(z) {
  rbind(mean = colMeans(z), q025 = apply(z, 2L, stats::quantile, 0.025),
        median = apply(z, 2L, stats::median),
        q975 = apply(z, 2L, stats::quantile, 0.975))
}
portable_predict <- function(model, x) {
  t(vapply(seq_len(nrow(model$portable_mean_draws)), function(i) {
    stats::approx(model$portable_bwo_grid, model$portable_mean_draws[i, ],
                  xout = x, ties = "ordered")$y
  }, numeric(length(x))))
}

gpr_jobs <- data.frame(
  id = c("GPR_Sub_main", "GPR_All_main", "GPR_Sub_shrinkage", "GPR_All_shrinkage"),
  set = c("main", "main", "sensitivity", "sensitivity"),
  source = c("Sub_GPR_Baseline_v0.3.rds", "All_GPR_Baseline_v0.3.rds",
             "Sub_GPR_Shrinkage_v0.3.rds", "All_GPR_Shrinkage_v0.3.rds"),
  portable = c("Sub_GPR_Baseline_v0.3_portable.rds", "All_GPR_Baseline_v0.3_portable.rds",
               "Sub_GPR_Shrinkage_v0.3_portable.rds", "All_GPR_Shrinkage_v0.3_portable.rds")
)
for (i in seq_len(nrow(gpr_jobs))) {
  job <- gpr_jobs[i, ]
  full <- readRDS(file.path(full_root, job$set, job$source))
  portable <- readRDS(file.path(portable_root, job$set, job$portable))
  tr <- full$settings$x_transform
  predict_full <- function(x, seed) {
    nd <- data.frame(O2 = x, O2_std = (x - tr$mean) / tr$sd, dd13C_sd = 0)
    set.seed(seed)
    as.matrix(brms::posterior_epred(full$fit_object, newdata = nd,
                                    draw_ids = portable$source_draw_ids,
                                    re_formula = NA))
  }
  compare_distribution <- function(x) {
    a <- summarize_draws(predict_full(x, 901L))
    b <- summarize_draws(predict_full(x, 902L))
    p <- summarize_draws(portable_predict(portable, x))
    repeat_difference <- max(abs(a - b))
    difference <- max(abs(a - p))
    c(difference, repeat_difference, max(0.005, 2 * repeat_difference))
  }
  exact <- compare_distribution(validation_grid)
  interpolated <- compare_distribution(interpolation_grid)
  add_check(job$id, "stored-grid posterior summaries", exact[1], exact[2], exact[3])
  add_check(job$id, "half-step posterior summaries", interpolated[1], interpolated[2], interpolated[3])
}

# The latent-input forward function is deterministic for each stored posterior
# draw, so its compact curves can be compared draw by draw.
if (!requireNamespace("rstan", quietly = TRUE)) stop("Latent GPR validation requires rstan.")
latent_full <- readRDS(file.path(full_root, "sensitivity", "Sub_GPR_LatentInput_v0.3.rds"))
latent_portable <- readRDS(file.path(portable_root, "sensitivity", "Sub_GPR_LatentInput_v0.3_portable.rds"))
ext <- rstan::extract(latent_full$fit_object)
tr <- latent_full$settings$x_transform
latent_predict <- function(x) {
  x_std <- (x - tr$mean) / tr$sd
  out <- matrix(NA_real_, nrow = length(latent_portable$source_draw_ids), ncol = length(x))
  for (k in seq_along(latent_portable$source_draw_ids)) {
    i <- latent_portable$source_draw_ids[[k]]
    xt <- as.numeric(ext$x_true[i, ]); f <- as.numeric(ext$f_latent[i, ])
    eta <- ext$eta[i]; rho <- ext$rho[i]
    covariance <- outer(xt, xt, function(a, b) eta^2 * exp(-0.5 * ((a - b) / rho)^2))
    diag(covariance) <- diag(covariance) + 1e-6
    cross <- outer(x_std, xt, function(a, b) eta^2 * exp(-0.5 * ((a - b) / rho)^2))
    out[k, ] <- as.numeric(ext$alpha[i] + cross %*% solve(covariance, f))
  }
  out
}
latent_exact <- max(abs(latent_predict(validation_grid) - portable_predict(latent_portable, validation_grid)))
latent_interp <- max(abs(latent_predict(interpolation_grid) - portable_predict(latent_portable, interpolation_grid)))
add_check("GPR_Sub_latent", "stored grid points", latent_exact, tolerance = 1e-10)
add_check("GPR_Sub_latent", "half-step interpolation", latent_interp, tolerance = 0.005)

for (calibration in c("Sub", "All")) {
  original <- readRDS(file.path(full_root, "sensitivity", paste0(calibration, "_GPR_XJitter_Baseline_v0.3.rds")))
  diagnostic <- readRDS(file.path(portable_root, "sensitivity", paste0(calibration, "_GPR_XJitter_Baseline_v0.3_diagnostic.rds")))
  difference <- max(abs(as.matrix(original$member_curves[c("O2", "median", "q025", "q975")]) -
                        as.matrix(diagnostic$member_curves[c("O2", "median", "q025", "q975")])))
  add_check(paste0("GPR_", calibration, "_xjitter"), "diagnostic curves", difference)
}

results <- do.call(rbind, checks)
results$passed <- results$max_abs_difference <= results$tolerance
print(results, row.names = FALSE)
if (!all(results$passed)) stop("At least one portable-model validation failed.")
utils::write.csv(results, "PORTABLE_MODEL_VALIDATION.csv", row.names = FALSE)
