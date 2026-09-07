\name{psm-result-methods}
\alias{print.forward_psm_result}
\alias{print.inverse_psm_result}
\alias{summary.forward_psm_result}
\alias{summary.inverse_psm_result}
\alias{plot.forward_psm_result}
\alias{plot.inverse_psm_result}
\title{Methods for BayDelC Results}
\usage{
\method{print}{forward_psm_result}(x, ...)
\method{print}{inverse_psm_result}(x, ...)
\method{summary}{forward_psm_result}(object, ...)
\method{summary}{inverse_psm_result}(object, ...)
\method{plot}{forward_psm_result}(x, ...)
\method{plot}{inverse_psm_result}(x, ...)
}
\arguments{
\item{x, object}{A BayDelC forward or inverse result.}
\item{...}{Additional arguments; plotting arguments are passed to base graphics.}
}
\value{Print and plot methods return the input invisibly. Summary methods return
the result's summary data frame.}
\description{Print, summarize, or plot the median and 95-percent interval of a BayDelC result.}
