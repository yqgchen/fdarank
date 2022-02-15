#' Set working grid for a given list of time grids.
#' @param Lt A list of vectors holding the observed time points.
#' @param nGrid The number of equidistant grid points of the output grid. Default: 101.
#' @return A vector holding an equidistant time grid.
#' @export
setWorkGrid <- function(Lt, nGrid = 101) {
  workGrid <- range(unlist(Lt))
  seq(workGrid[1], workGrid[2], length.out = nGrid)
}