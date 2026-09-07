\name{forward_psm}
\alias{forward_psm}
\title{Propagate Bottom-Water Oxygen Through a BayDelC Calibration}
\usage{
forward_psm(
  data = NULL, bwo = NULL, bwo_sd = NULL, bwo_draws = NULL,
  bwo_draw_weights = NULL, bwo_col = "O2", bwo_sd_col = NULL,
  site_col = "Site",
  input_error_method = c("auto", "fixed", "normal", "draws", "custom"),
  input_distribution = NULL, input_bounds = c(0, Inf),
  algorithm = c("BLR", "GPR"),
  calibration = c("Sub", "All", "custom"), model = NULL,
  prediction = c("predictive", "latent"), include_model_residual = TRUE,
  n_draw = 2000, seed = 123, return_draws = FALSE,
  validate = c("strict", "warn", "none"), gpr_chunk_size = 200L
)
}
\arguments{
\item{data}{Optional data frame containing BWO and optional uncertainty/site columns.}
\item{bwo}{Numeric BWO values in micromoles per kilogram.}
\item{bwo_sd}{Optional one-standard-deviation BWO uncertainty.}
\item{bwo_draws}{Optional observations-by-realizations matrix of BWO draws.}
\item{bwo_draw_weights}{Optional weights for columns of \code{bwo_draws}.}
\item{bwo_col, bwo_sd_col, site_col}{Column names used with \code{data}.}
\item{input_error_method}{Input uncertainty representation.}
\item{input_distribution}{Optional custom input random generator.}
\item{input_bounds}{Two-element bounds for sampled BWO.}
\item{algorithm}{Model family: \code{"BLR"} or \code{"GPR"}.}
\item{calibration}{Calibration dataset or \code{"custom"}.}
\item{model}{Optional custom model list or RDS path.}
\item{prediction}{Predictive or latent-mean output.}
\item{include_model_residual}{Whether residual model variability is included.}
\item{n_draw}{Number of Monte Carlo draws.}
\item{seed}{Random seed.}
\item{return_draws}{Whether to return the observations-by-draws matrix.}
\item{validate}{Validation policy.}
\item{gpr_chunk_size}{Chunk size for custom full \pkg{brms} GPR models.}
}
\value{An object of class \code{forward_psm_result} containing a summary,
optional draws, diagnostics, and model metadata.}
\description{Estimates Delta-delta-13-C expected from one or more BWO values
while propagating calibration and input uncertainty.}
\examples{
forward_psm(bwo = c(80, 150, 220), algorithm = "BLR", n_draw = 100)
}
