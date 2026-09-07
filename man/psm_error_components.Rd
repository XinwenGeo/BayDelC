\name{psm_error_components}
\alias{psm_error_components}
\title{Describe Delta-delta-13-C Uncertainty Components}
\usage{
psm_error_components(
  epi_total_sd = NULL, infa_total_sd = NULL,
  epi_measurement_sd = NULL, epi_sample_sd = NULL,
  infa_measurement_sd = NULL, infa_sample_sd = NULL,
  pairing_sd = 0, other_sd = 0, rho_epi_infa = 0,
  scale = c("sd", "se")
)
}
\arguments{
\item{epi_total_sd, infa_total_sd}{Total epifaunal or infaunal standard deviation.}
\item{epi_measurement_sd, epi_sample_sd}{Epifaunal error components.}
\item{infa_measurement_sd, infa_sample_sd}{Infaunal error components.}
\item{pairing_sd, other_sd}{Additional independent error terms.}
\item{rho_epi_infa}{Correlation between epifaunal and infaunal errors.}
\item{scale}{Input scale, \code{"sd"} or \code{"se"}.}
}
\value{An object for the \code{error_components} argument of \code{inverse_psm()}.}
\description{Records the components used to calculate proxy-observation uncertainty.}
