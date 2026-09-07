#' @export
print.forward_psm_result <- function(x, ...) {
  cat("BayDelC forward prediction\n")
  cat(sprintf("Model: %s/%s/%s (%s calibration)\n", x$metadata$model_set,
              x$metadata$algorithm, x$metadata$variant, x$metadata$calibration))
  print(x$summary, row.names = FALSE)
  invisible(x)
}

#' @export
print.inverse_psm_result <- function(x, ...) {
  cat("BayDelC inverse reconstruction\n")
  cat(sprintf("Model: %s/%s/%s (%s calibration); inversion: %s\n",
              x$metadata$model_set, x$metadata$algorithm, x$metadata$variant,
              x$metadata$calibration, x$metadata$inversion))
  print(x$summary, row.names = FALSE)
  invisible(x)
}

#' @export
summary.forward_psm_result <- function(object, ...) object$summary

#' @export
summary.inverse_psm_result <- function(object, ...) object$summary

.baydelc_interval_plot <- function(x, xlab, ylab, ...) {
  tab <- x$summary
  graphics::plot(tab$input_value, tab$median,
                 ylim = range(c(tab$q025, tab$q975), finite = TRUE),
                 xlab = xlab, ylab = ylab, pch = 19, ...)
  graphics::segments(tab$input_value, tab$q025,
                     tab$input_value, tab$q975)
  graphics::points(tab$input_value, tab$median, pch = 19)
  invisible(x)
}

#' @export
plot.forward_psm_result <- function(x, ...) {
  .baydelc_interval_plot(
    x, "Bottom-water oxygen (micromol kg^-1)",
    expression(Delta * delta^{13} * C ~ (per~mille)), ...
  )
}

#' @export
plot.inverse_psm_result <- function(x, ...) {
  .baydelc_interval_plot(
    x, expression(Delta * delta^{13} * C ~ (per~mille)),
    "Bottom-water oxygen (micromol kg^-1)", ...
  )
}
