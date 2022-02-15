#' Compute rank partial derivatives
#' @description Compute rank partial derivatives \eqn{D_1(Y_i(t_j),t_j)} and \eqn{D_2(Y_i(t_j),t_j)}, 
#' where \eqn{D_1(y,t) = \partial F_t(y)/\partial t} and \eqn{D_2(y,t) = \partial F_t(y)/\partial y}.
#' @param Ly A list of \eqn{n} vectors holding the observed values for each subject.
#' @param Lt A list of \eqn{n} vectors holding the observation time points for each subject. 
#' Each vector should be in ascending order and of the same length as the corresponding element in \code{Ly}.
#' @param optns A list of options control parameters specified by \code{list(name=value)}. 
#' Available control options are \code{tWorkGrid}, \code{yWorkMat}, \code{tbw}, \code{ybw}, 
#' \code{Kty}, \code{Hty}, \code{ncores}; see \code{\link{smoRank}} for details.
#' @import stats
#' @return A list of the following elements:
#' \describe{
#' \item{D1, D2}{Matrices of dimension \eqn{n}-by-\eqn{l} of which the \eqn{(i,j)}-th entry holds the corresponding derivatives for the \eqn{i}-th subject evaluated at \code{tWorkGrid[j]} and \code{yWorkMat[i,j]}.}
#' \item{optns}{A list of control options used in computing \code{D1} and \code{D2}.}
#' }
#' @export
rankPtDeriv <- function(Ly, Lt, optns = list()) {
  tWorkGrid <- optns$tWorkGrid
  yWorkMat <- optns$yWorkMat
  tbw <- optns$tbw
  ybw <- optns$ybw
  Kty <- optns$Kty
  Hty <- optns$Hty
  
  yWorkGrid <- setWorkGrid(Ly)
  
  if (is.null(tWorkGrid))
    tWorkGrid <- setWorkGrid(Lt)
  if (is.null(yWorkMat))
    yWorkMat <- t(sapply(seq_along(Ly), function(i) approx(x = Lt[[i]], y = Ly[[i]], xout = tWorkGrid, rule = 2)$y))
  
  if (is.null(tbw) | is.null(ybw)) {
    optns1 <- optns
    optns1$bwchoice <- TRUE
    ybw <- smoRank(Ly = Ly, Lt = Lt, optns = optns1)
    tbw <- ybw$tbw
    ybw <- ybw$ybw
  }
  
  if (is.null(Kty))
    Kty <- 'epan'
  if (is.null(Hty))
    Hty <- 'epan'
  
  yin <- unlist(Ly)
  tin <- unlist(Lt)
  
  K <- kerfctn( kernel_type = Kty, fctn_type = "pdf" )
  H <- kerfctn( kernel_type = Hty, fctn_type = "cdf" )
  dK <- kerfctn( kernel_type = Kty, fctn_type = "dpdf" )
  
  d1 <- sapply(tWorkGrid, function(t) {
    sapply(yWorkGrid, function(y) {
      mean(H((y-yin)/ybw)*dK((t-tin)/tbw))/ tbw^2 / mean(K((t-tin)/tbw)) * tbw - mean(H((y-yin)/ybw)*K((t-tin)/tbw)) * mean(dK((t-tin)/tbw)) /tbw / mean(K((t-tin)/tbw))^2
    })
  })
  d2 <- sapply(tWorkGrid, function(t) {
    sapply(yWorkGrid, function(y) {
      mean(K((y-yin)/ybw)*K((t-tin)/tbw))/ (ybw*tbw) / mean(K((t-tin)/tbw)) * tbw
    })
  })
  
  d1 <- sapply(seq_along(tWorkGrid), function(j) approx(x = yWorkGrid, y = d1[,j], xout = yWorkMat[,j])$y)
  d2 <- sapply(seq_along(tWorkGrid), function(j) approx(x = yWorkGrid, y = d2[,j], xout = yWorkMat[,j])$y)
  optns <- list(tWorkGrid = tWorkGrid, yWorkMat = yWorkMat, tbw = tbw, ybw = ybw, Kty = Kty, Hty = Hty)
  return(list(D1 = d1, D2 = d2, optns = optns))
}
