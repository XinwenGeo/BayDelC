library(BayDelC)

result <- inverse_psm(
  dd13c = 1.50,
  dd13c_sd = 0.08,
  algorithm = "GPR",
  calibration = "Sub",
  inversion = "grid",
  prior = "normal",
  prior_args = list(mean = 140, sd = 60, bounds = c(0, 300)),
  o2_step = 1,
  integration_draws = 500,
  n_draw = 2000,
  seed = 2026,
  return_draws = TRUE,
  return_density = TRUE
)

print(result)
print(result$diagnostics)
