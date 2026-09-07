# Package-native reproduction of the BayDelC article Figure 2.
# Requires the optional packages ggplot2 and patchwork.

library(BayDelC)
if (!requireNamespace("ggplot2", quietly = TRUE) ||
    !requireNamespace("patchwork", quietly = TRUE)) {
  stop("Install ggplot2 and patchwork to run this figure example.")
}

output_dir <- Sys.getenv("BAYDELC_EXAMPLE_OUTPUT", "BayDelC_paper_examples")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
n_draw <- 1000L

h25 <- data.frame(calibration = c("Sub", "All"),
                  intercept = c(0.093, 0.232), slope = c(0.0110, 0.0094))

curve_data <- function(algorithm, calibration, seed) {
  model <- baydelc_load_model(algorithm, calibration, model_set = "main")
  grid <- seq(model$calibration_range[1], model$calibration_range[2], length.out = 151)
  latent <- forward_psm(bwo = grid, algorithm = algorithm,
                        calibration = calibration, model_set = "main",
                        prediction = "latent", n_draw = n_draw, seed = seed)
  predictive <- forward_psm(bwo = grid, algorithm = algorithm,
                            calibration = calibration, model_set = "main",
                            prediction = "predictive", n_draw = n_draw, seed = seed)
  out <- latent$summary[c("BWO", "median", "q025", "q975")]
  names(out)[2:4] <- c("fit", "ci_low", "ci_high")
  out$pi_low <- predictive$summary$q025
  out$pi_high <- predictive$summary$q975
  h <- h25[h25$calibration == calibration, ]
  out$h25 <- h$intercept + h$slope * out$BWO
  list(curve = out, training = model$training_data, model = model)
}

make_panel <- function(algorithm, calibration, seed) {
  z <- curve_data(algorithm, calibration, seed)
  ggplot2::ggplot(z$curve, ggplot2::aes(BWO, fit)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = pi_low, ymax = pi_high),
                         fill = "grey75", alpha = 0.35) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = ci_low, ymax = ci_high),
                         fill = if (algorithm == "BLR") "#F8766D" else "#619CFF",
                         alpha = 0.30) +
    ggplot2::geom_line(color = if (algorithm == "BLR") "red3" else "blue3",
                       linewidth = 0.8) +
    ggplot2::geom_line(ggplot2::aes(y = h25), color = "black",
                       linetype = "dashed", linewidth = 0.6) +
    ggplot2::geom_errorbar(data = z$training,
      ggplot2::aes(x = O2, y = dd13C, ymin = dd13C - dd13C_sd,
                   ymax = dd13C + dd13C_sd), inherit.aes = FALSE,
      width = 3, color = "grey40") +
    ggplot2::geom_point(data = z$training, ggplot2::aes(O2, dd13C),
                        inherit.aes = FALSE, shape = 21, fill = "skyblue",
                        color = "black", size = 2) +
    ggplot2::coord_cartesian(xlim = c(0, 300), ylim = c(-0.5, 4.5)) +
    ggplot2::labs(title = paste(calibration, algorithm),
      x = expression("Bottom-water O"[2] * " (" * mu * "mol kg"^-1 * ")"),
      y = expression(Delta * delta^13 * C ~ ("‰"))) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
}

figure <- patchwork::wrap_plots(
  make_panel("BLR", "Sub", 1201), make_panel("BLR", "All", 1202),
  make_panel("GPR", "Sub", 2201), make_panel("GPR", "All", 2202),
  ncol = 2
) + patchwork::plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")

ggplot2::ggsave(file.path(output_dir, "Figure02_package_reproduction.pdf"),
                figure, width = 180, height = 160, units = "mm")
ggplot2::ggsave(file.path(output_dir, "Figure02_package_reproduction.png"),
                figure, width = 180, height = 160, units = "mm", dpi = 300)
