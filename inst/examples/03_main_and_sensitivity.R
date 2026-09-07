library(BayDelC)

# The same observations under the primary Gaussian BLR and Student-t
# sensitivity calibration.
main <- inverse_psm(
  dd13c = c(1.1, 1.5, 2.1), dd13c_sd = 0.08,
  algorithm = "BLR", calibration = "Sub",
  model_set = "main", n_draw = 2000, seed = 2026
)
sensitivity <- inverse_psm(
  dd13c = c(1.1, 1.5, 2.1), dd13c_sd = 0.08,
  algorithm = "BLR", calibration = "Sub",
  model_set = "sensitivity", variant = "student_t",
  n_draw = 2000, seed = 2026
)

comparison <- rbind(
  transform(summary(main), selected_model = "main Gaussian BLR"),
  transform(summary(sensitivity), selected_model = "sensitivity Student-t BLR")
)
print(comparison[c("input_value", "median", "q025", "q975", "selected_model")])

# GPR sensitivity variants are selected the same way.
shrinkage <- forward_psm(
  bwo = c(80, 150, 220), algorithm = "GPR", calibration = "Sub",
  model_set = "sensitivity", variant = "shrinkage", n_draw = 1000
)
latent_input <- forward_psm(
  bwo = c(80, 150, 220), algorithm = "GPR", calibration = "Sub",
  model_set = "sensitivity", variant = "latent_input", n_draw = 1000
)
print(summary(shrinkage))
print(summary(latent_input))
