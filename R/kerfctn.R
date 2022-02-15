#' Kernel function
#' @param kernel_type Character holding the kernel type: \code{"rect"}, \code{"gaus"}, \code{"epan"} (default), \code{"tgaus"}.
#' @param fctn_type Character holding the function type: \code{"pdf"} (default) , \code{"cdf"}, \code{"dpdf"} (derivatives of pdfs).
#' @import stats
#' @return A function computing the kernel function.
#' @export
kerfctn <- function ( kernel_type = "epan", fctn_type = "pdf" ) {
  if ( ! kernel_type %in% c( "epan", "gaus", "rect", "tgaus" ) ) {
    stop('Unavailable kernel type.')
  }
  if ( ! fctn_type %in% c( "pdf", "cdf", "dpdf" ) ) {
    stop('Unavailable function type.')
  }
  if ( fctn_type == "pdf" ) {
    if ( kernel_type == 'gaus' ) {
      K <- function(x) dnorm(x)
    } else if ( kernel_type == 'rect' ) {
      K <- function(x) dunif (x, -1, 1)
    } else if ( kernel_type == 'epan' ) {
      K <- function(x) 0.75 * (1-x^2) * (abs(x)<=1)
    } else if ( kernel_type == 'tgaus' ) {
      K <- function(x) dnorm(x) * (abs(x)<=4) / (pnorm(4) - pnorm(-4))
    }
  } else if ( fctn_type == "cdf" ) {
    if (kernel_type=='gaus') {
      K <- function(x) pnorm(x)
    } else if (kernel_type=='rect') {
      K <- function(x) punif (x,-1,1)
    } else if (kernel_type=='epan') {
      K <- function(x) (0.75 * x - x^3 / 4 + 1/2) * (abs(x)<=1) + (x>1)
    } else if (kernel_type=='tgaus') {
      K <- function(x) (pnorm(x) - pnorm(-4)) * (abs(x)<=4) / (pnorm(4) - pnorm(-4)) + (x>4)
    }
  } else if ( fctn_type == "dpdf" ) {
    if (kernel_type=='gaus') {
      K <- function(x) -x * exp(-x^2 / 2) / sqrt(2*pi)
    } else if (kernel_type=='rect') {
      K <- function(x) numeric(length(x))
    } else if (kernel_type=='epan') {
      K <- function(x) { - 1.5 * x * (abs(x)<=1) }
    } else if (kernel_type=='tgaus') {
      K <- function(x) { - x * exp(-x^2 / 2) / sqrt(2*pi) * (abs(x)<=4) / (pnorm(4) - pnorm(-4)) }
    }
  }
  return (K)
}
