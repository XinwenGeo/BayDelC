library(BayDelC)

print(baydelc_models())

forward_result <- forward_psm(
  bwo = c(75, 125, 200),
  bwo_sd = c(5, 10, 15),
  algorithm = "BLR",
  calibration = "Sub",
  n_draw = 2000,
  seed = 2026,
  return_draws = TRUE
)
print(forward_result)

inverse_result <- inverse_psm(
  dd13c = c(1.10, 1.50, 2.10),
  dd13c_sd = c(0.08, 0.08, 0.10),
  algorithm = "BLR",
  calibration = "Sub",
  n_draw = 2000,
  seed = 2026,
  return_draws = TRUE
)
print(inverse_result)
print(inverse_result$diagnostics)
