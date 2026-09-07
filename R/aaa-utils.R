# Internal utilities -------------------------------------------------------

`%||%` <- function(a, b) if (!is.null(a)) a else b

#' Report the installed BayDelC version
#'
#' @return A character scalar.
#' @export
baydelc_version <- function() {
  as.character(utils::packageVersion("BayDelC"))
}
.baydelc_assert_scalar_string <- function(x, name) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop(sprintf("`%s` must be one non-empty character string.", name), call. = FALSE)
  }
  invisible(x)
}
