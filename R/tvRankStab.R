#' Compute time-dependent population rank stability
#' @param rankDeriv An \eqn{n}-by-\eqn{l} matrix holding the rank derivatives evaluated on \code{tgrid}.
#' @return A vector of length \eqn{l} holding the time-dependent population rank stability.
#' @export
tvRankStab <- function (rankDeriv) {
  if ( is.vector(rankDeriv) ) rankDeriv <- matrix( rankDeriv, nrow = 1 )
  colMeans( rankDeriv^2 )
}
