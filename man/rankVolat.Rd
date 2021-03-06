% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rankVolat.R
\name{rankVolat}
\alias{rankVolat}
\title{Compute subject-specific rank volatility}
\usage{
rankVolat(rank, tgrid)
}
\arguments{
\item{rank}{A vector of length \eqn{l} or an \eqn{n}-by-\eqn{l} matrix holding the ranks evaluated on \code{tgrid}.}

\item{tgrid}{A vector of length \eqn{l} holding the time grid where the ranks are evaluated.}
}
\value{
A vector holding the subject-specific rank volatilities for each subject.
}
\description{
Compute subject-specific rank volatility
}
