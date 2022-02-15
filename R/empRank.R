#' Compute empirical ranks
#' @param Ly A list of \eqn{n} vectors holding the observed values for each subject.
#' @param Lt A list of \eqn{n} vectors holding the observation time points for each subject. 
#' Each vector should be in ascending order and of the same length as the corresponding element in \code{Ly}.
#' @param workGrid A vector holding the time grid where the empirical ranks are computed. Default: set up by \code{\link{setWorkGrid}}.
#' @return An \eqn{n}-by-\eqn{l} matrix of which each column holds the empirical ranks for one subject evaluated on \code{workGrid},
#' where \eqn{l} is the length of \code{workGrid}.
#' @importFrom dplyr near
#' @importFrom stats approx
#' @export
empRank <- function(Ly, Lt, workGrid = NULL) {
  if ( is.null(workGrid) ) workGrid <- setWorkGrid(Lt)
  data <- sapply(seq_along(Ly), function(i) {
    noInterpol <- dplyr::near( length(Lt[[i]]), length(workGrid) )
    if ( noInterpol ) {
      noInterpol <- all( near( Lt[[i]], workGrid ) )
    }
    if ( noInterpol ) {
      return ( Ly[[i]] )
    } else {
      approx(x = Lt[[i]], y = Ly[[i]], xout = workGrid, rule = 2)$y
    }
  })
  if ( is.vector(data) ) data <- matrix( data, nrow = 1 )
  n <- ncol(data)
  apply(data, 1, function(x) {
    sapply(1:n,function(i) {
      (sum(x<=x[i]) - 1)
    })/n
  })
}
