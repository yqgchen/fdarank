#' Compute subject-specific integrated ranks
#' @param rank A vector of length \eqn{l} or an \eqn{n}-by-\eqn{l} matrix holding the ranks evaluated on \code{tgrid}.
#' @param tgrid A vector of length \eqn{l} holding the time grid where the ranks are evaluated.
#' @return A vector holding the subject-specific integrated ranks for each subject.
#' @export
intRank <- function (rank, tgrid) {
  if ( is.vector(rank) ) rank <- matrix( rank, nrow = 1 )
  apply(rank, 1, pracma::trapz, x = normalizeGrid(tgrid) )
}