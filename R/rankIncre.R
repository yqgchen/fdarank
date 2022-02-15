#' Compute subject-specific rank increments.
#' @param rank A vector of length \eqn{l} or an \eqn{n}-by-\eqn{l} matrix holding the ranks evaluated on \code{tgrid}.
#' @param tgrid A vector of length \eqn{l} holding the time grid where the ranks are evaluated.
#' @return A vector holding the difference between the ending and starting ranks for each subject.
#' @export
rankIncre <- function (rank, tgrid) {
  if ( is.vector(rank) ) rank <- matrix( rank, nrow = 1 )
  rank[, which.max(tgrid)] - rank[, which.min(tgrid)]
}
