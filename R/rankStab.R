#' Compute overall population rank stability
#' @param rankDeriv An \eqn{n}-by-\eqn{l} matrix holding the rank derivatives evaluated on \code{tgrid}.
#' @param tgrid A vector of length \eqn{l} holding the time grid where the rank derivatives are evaluated.
#' @return A scalar holding the overall population rank stability.
#' @export
rankStab <- function (rankDeriv, tgrid) {
  tgrid <- normalizeGrid(tgrid)
  exp( - pracma::trapz( x = tgrid, y = tvRankStab( rankDeriv = rankDeriv ) ) )
}

