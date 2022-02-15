#' Compute subject-specific rank volatility
#' @param rank A vector of length \eqn{l} or an \eqn{n}-by-\eqn{l} matrix holding the ranks evaluated on \code{tgrid}.
#' @param tgrid A vector of length \eqn{l} holding the time grid where the ranks are evaluated.
#' @return A vector holding the subject-specific rank volatilities for each subject.
#' @export
rankVolat <- function (rank, tgrid) {
  if ( is.vector(rank) ) rank <- matrix( rank, nrow = 1 )
  tgrid <- normalizeGrid(tgrid)
  intrank <- intRank(rank = rank, tgrid = tgrid )
  apply( (rank - intrank)^2, 1, pracma::trapz, x = tgrid )
  # apply( (rank - matrix(intrank, nrow = nrow(rank), ncol = ncol(rank)))^2, 1, pracma::trapz, x = tgrid )
}

