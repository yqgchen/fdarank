#' Normalize time grid
#' @param tgrid A vector holding the time grid to be normalized.
#' @return A grid on \eqn{[0,1]}.
normalizeGrid <- function (tgrid) {
  range <- range(tgrid)
  (tgrid - range[1]) / diff(range)
}
