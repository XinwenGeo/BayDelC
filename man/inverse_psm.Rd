\name{inverse_psm}
\alias{inverse_psm}
\title{Reconstruct Bottom-Water Oxygen from Delta-delta-13-C}
\usage{
inverse_psm(
  data = NULL, dd13c = NULL, dd13c_sd = NULL, glob_site_sd = NULL,
  error_components = NULL, dd13c_draws = NULL,
  dd13c_draw_weights = NULL, dd13c_col = "dd13C",
  dd13c_sd_col = NULL, glob_site_sd_col = NULL, site_col = "Site",
  error_method = c("auto", "provided_total", "h15_proxy_default",
    "lu_site", "components", "draws"),
  algorithm = c("BLR", "GPR"),
  calibration = c("Sub", "All", "custom"),
  model_set = c("main", "sensitivity"), variant = NULL, model = NULL,
  inversion = c("auto", "strict", "grid"),
  prior = c("none", "uniform_physical", "normal", "draws", "density", "custom"),
  prior_args = NULL, o2_bounds = c(0, 300), o2_step = 1, o2_grid = NULL,
  dd13c_convention = "epifaunal_minus_infaunal",
  include_model_residual = TRUE, n_draw = 2000, integration_draws = NULL,
  seed = 123, return_draws = FALSE, return_density = FALSE,
  validate = c("strict", "warn", "none")
)
}
\arguments{
\item{data}{Optional data frame containing observations and uncertainty/site columns.}
\item{dd13c}{Numeric Delta-delta-13-C observations in per mille.}
\item{dd13c_sd}{Total one-standard-deviation observation uncertainty.}
\item{glob_site_sd}{Site-specific infaunal uncertainty.}
\item{error_components}{Output of \code{psm_error_components()} or a compatible list.}
\item{dd13c_draws}{Optional observations-by-realizations matrix.}
\item{dd13c_draw_weights}{Optional weights for draw columns.}
\item{dd13c_col, dd13c_sd_col, glob_site_sd_col, site_col}{Column names.}
\item{error_method}{Observation-error representation.}
\item{algorithm}{Model family: \code{"BLR"} or \code{"GPR"}.}
\item{calibration}{Calibration dataset or \code{"custom"}.}
\item{model_set}{Primary \code{"main"} or alternative \code{"sensitivity"} models.}
\item{variant}{Optional model variant; see \code{baydelc_models()}.}
\item{model}{Optional custom model list or RDS path.}
\item{inversion}{Automatic, analytic BLR, or grid inversion.}
\item{prior}{Prior family for grid inversion.}
\item{prior_args}{Named list configuring the prior.}
\item{o2_bounds, o2_step, o2_grid}{Numerical BWO-grid controls.}
\item{dd13c_convention}{Proxy sign convention.}
\item{include_model_residual}{Whether residual calibration variability is included.}
\item{n_draw}{Number of output Monte Carlo draws.}
\item{integration_draws}{Number of model draws used in likelihood integration.}
\item{seed}{Random seed.}
\item{return_draws, return_density}{Whether to return draws or density values.}
\item{validate}{Validation policy.}
}
\value{An object of class \code{inverse_psm_result} containing a summary,
optional draws/densities, diagnostics, and model metadata.}
\description{Propagates calibration and observation uncertainty to estimate BWO.}
\examples{
inverse_psm(dd13c = 1.5, dd13c_sd = 0.08, n_draw = 100)
}
