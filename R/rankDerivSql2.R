#' Compute subject-specific squared \eqn{L^2} norm of the rank derivatives
#' @param rankDeriv A vector of length \eqn{l} or an \eqn{n}-by-\eqn{l} matrix holding the rank derivatives evaluated on \code{tgrid}.
#' @param tgrid A vector of length \eqn{l} holding the time grid where the rank derivatives are evaluated.
#' @return A vector holding the squared \eqn{L^2} norm of the rank derivatives for each subject.
#' @export
rankDerivSql2 <- function (rankDeriv, tgrid) {
  tgrid <- normalizeGrid(tgrid)
  if ( is.vector(rankDeriv) ) rankDeriv <- matrix( rankDeriv, nrow = 1 )
  apply(rankDeriv^2, 1, pracma::trapz, x = tgrid)
}
