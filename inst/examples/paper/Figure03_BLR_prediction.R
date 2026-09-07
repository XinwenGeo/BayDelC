# Package-native reproduction of the BayDelC article Figure 3.
# Requires the optional packages ggplot2 and patchwork.

library(BayDelC)
if (!requireNamespace("ggplot2", quietly = TRUE) ||
    !requireNamespace("patchwork", quietly = TRUE)) {
  stop("Install ggplot2 and patchwork to run this figure example.")
}

output_dir <- Sys.getenv("BAYDELC_EXAMPLE_OUTPUT", "BayDelC_paper_examples")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

make_panel <- function(calibration, fill, seed) {
  model <- baydelc_load_model("BLR", calibration, model_set = "main")
  dat <- model$training_data
  prediction <- forward_psm(
    data = dat, bwo_col = "O2", bwo_sd_col = "O2_sd",
    algorithm = "BLR", calibration = calibration, model_set = "main",
    prediction = "latent", n_draw = 2000, seed = seed
  )$summary
  dat$predicted <- prediction$median
  dat$low <- prediction$q025
  dat$high <- prediction$q975
  rmse <- sqrt(mean((dat$dd13C - dat$predicted)^2))
  r2 <- 1 - stats::var(dat$dd13C - dat$predicted) / stats::var(dat$dd13C)
  ggplot2::ggplot(dat, ggplot2::aes(dd13C, predicted)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = low, ymax = high), width = 0.04,
                           color = "grey55") +
    ggplot2::geom_point(shape = 21, fill = fill, size = 2.5) +
    ggplot2::annotate("text", x = 0.1, y = 3.35,
      label = sprintf("RMSE: %.2f ‰\nBayesian R²: %.2f", rmse, r2),
      hjust = 0, vjust = 1, size = 3.2) +
    ggplot2::coord_cartesian(xlim = c(0, 3.5), ylim = c(0, 3.5)) +
    ggplot2::labs(title = calibration,
      x = expression("Observed " * Delta * delta^13 * C ~ ("‰")),
      y = expression("Predicted " * Delta * delta^13 * C ~ ("‰"))) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), aspect.ratio = 1)
}

figure <- make_panel("Sub", "skyblue", 3101) +
  make_panel("All", "#0072B2", 3102) +
  patchwork::plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")

ggplot2::ggsave(file.path(output_dir, "Figure03_package_reproduction.pdf"),
                figure, width = 180, height = 90, units = "mm")
ggplot2::ggsave(file.path(output_dir, "Figure03_package_reproduction.png"),
                figure, width = 180, height = 90, units = "mm", dpi = 300)
