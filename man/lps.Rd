% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lps.R
\name{lps}
\alias{lps}
\title{Local polynomial smoothing}
\usage{
lps(xin, yin, xout, kty = "epan", bw, npoly = 1)
}
\arguments{
\item{xin}{A vector of length \eqn{n} with measurement points.}

\item{yin}{A vector of length \eqn{n} with measurement values.}

\item{xout}{A vector of length \eqn{m} with output measurement points.}

\item{kty}{A character holding the kernel type; default: \code{'epan'}.}

\item{bw}{A scalar holding the bandwidth.}

\item{npoly}{The order of polynomials used; e.g., \code{1} for local linear smoothing.}
}
\value{
An \eqn{m}-by-\eqn{p} matrix of which columns hold the estimates of the function 
and its derivatives up to order \eqn{p} evaluated on \code{xout}, where \eqn{p} equals (\code{npoly-1}).
}
\description{
Local polynomial smoothing
}
