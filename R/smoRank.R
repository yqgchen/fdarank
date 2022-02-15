#' Compute smooth ranks
#' @param Ly A list of \eqn{n} vectors holding the observed values for each subject.
#' @param Lt A list of \eqn{n} vectors holding the observation time points for each subject. 
#' Each vector should be in ascending order and of the same length as the corresponding element in \code{Ly}.
#' @param optns A list of options control parameters specified by \code{list(name=value)}. See `Details'.
#' @details Available control options are 
#' \describe{
#' \item{tWorkGrid}{A vector of length \eqn{l} holding the time grid where the smooth ranks are computed.}
#' \item{yWorkMat}{An \eqn{n}-by-\eqn{l} matrix; each row holding the values of the process for one subject observed on \code{tWorkGrid} where the smooth ranks are computed.}
#' \item{tbw, ybw}{Positive scalars holding the bandwidths.}
#' \item{Kty, Hty}{A character string holding the type of kernel functions; \code{"rect"}, \code{"gaus"}, \code{"epan"}, \code{"tgaus"}; default: \code{"epan"}.}
#' \item{training}{A scalar between 0 and 1 holding the proportion of the subjects used for choosing the bandwidths; default: 1.}
#' \item{ncores}{Number of cores used in parallel; default: 1.}
#' \item{bwchoice}{logical. If \code{TRUE}, only return the choices for \code{ybw} and \code{tbw}. Default: \code{FALSE}.}
#' }
#' @import stats
#' @import doParallel
#' @import parallel
#' @import foreach
#' @importFrom pracma trapz
#' @return A list of the following elements if \code{bwchoice} is \code{FALSE} and a list including \code{ybw} and \code{tbw} otherwise.
#' \describe{
#' \item{rank}{An \eqn{n}-by-\eqn{l} matrix of which the \eqn{(i,j)}-th entry holds the smooth ranks of the \eqn{i}-th subject evaluated at \code{tWorkGrid[j]} and \code{yWorkMat[i,j]}.}
#' \item{optns}{A list of control options used in computing \code{rank}.}
#' }
#' @export
smoRank <- function(Ly, Lt, optns = list()) {
  tWorkGrid <- optns$tWorkGrid
  yWorkMat <- optns$yWorkMat
  tbw <- optns$tbw
  ybw <- optns$ybw
  Kty <- optns$Kty
  Hty <- optns$Hty
  training <- optns$training
  ncores <- optns$ncores
  bwchoice <- optns$bwchoice
  
  yWorkGrid <- setWorkGrid(Ly)
  
  if (is.null(tWorkGrid))
    tWorkGrid <- setWorkGrid(Lt)
  if (is.null(yWorkMat))
    yWorkMat <- t(sapply(seq_along(Ly), function(i) approx(x = Lt[[i]], y = Ly[[i]], xout = tWorkGrid, rule = 2)$y))
  if (is.null(Kty))
    Kty <- 'epan'
  if (is.null(Hty))
    Hty <- 'epan'
  if (is.null(training))
    training <- 1
  if (is.null(bwchoice))
    bwchoice <- FALSE
  
  K <- kerfctn( kernel_type = Kty, fctn_type = "pdf" )
  H <- kerfctn( kernel_type = Hty, fctn_type = "cdf" )
  
  cdfsm <- function(tin, yin, tout, yout, tbw, ybw, K, H) {
    res <- sapply(tout, function(t) {
      sapply(yout, function(y) {
        mean(H((y-yin)/ybw)*K((t-tin)/tbw))/tbw / mean(K((t-tin)/tbw)) *tbw
      })
    })
    if (ncol(res) == 1) res <- c(res)
    return(res)
  }
  
  if (is.null(tbw) | is.null(ybw)) {
    
    # require(foreach)
    # require(doParallel)
    #if (!require(foreach)) stop("The package foreach was not installed.")
    #if (!require(doParallel)) stop("The package doParallel was not installed.")
    
    setBwCand <- function(xin, xout) {
      x <- unique(xin)
      xbwCand <- max(min(diff(sort(x))), sort(x - min(xout))[2], sort(max(xout) - x)[2]) * 1.2
      xbwCand <- seq(xbwCand^(1/2), (diff(range(x))/4)^(1/2), length.out = 10)^2
      return(xbwCand)
    }
    n <- length(Ly)
    trIdx <- sample(1:n, round(n*training))
    yin <- unlist(Ly[trIdx])
    tin <- unlist(Lt[trIdx])
    if (is.null(ybw)) {
      ybwCand <- setBwCand(yin, yWorkGrid)
    } else ybwCand <- ybw
    if (is.null(tbw)) {
      tbwCand <- setBwCand(tin, tWorkGrid)
    } else tbwCand <- tbw
    
    AISE <- array(dim = c(length(ybwCand), length(tbwCand)))
    
    cores <- detectCores()
    if (is.null(ncores)) {
      # ncores <- cores[1] - 1
      ncores <- 1
    } else {
      if (!is.integer(ncores)) {
        ncores <- round(ncores)
      }
      if (ncores > cores[1] - 1) {
        ncores <- cores[1] - 1
      }
    }
    if ( ncores > 1 ) {
      cl <- makeCluster(ncores)
      registerDoParallel(cl)
      `%doforeach%` <- `%dopar%`
    } else {
      `%doforeach%` <- `%do%`
    }
    idxInnerPts <- which( tin > min(tWorkGrid) + max(tbwCand) & tin < max(tWorkGrid) - max(tbwCand) )
    for (tbwIdx in seq_along(ybwCand)) {
      for (ybwIdx in seq_along(tbwCand)) {
        ise <- foreach(j = idxInnerPts, .combine=c ) %doforeach% {
          est <- cdfsm(tin=tin[-j], yin=yin[-j], tout=tin[j], yout=yWorkGrid,
                       ybw=ybwCand[ybwIdx], tbw=tbwCand[tbwIdx], H=H, K=K)
          err <- ((yWorkGrid >= yin[j]) - est)^2
          pracma::trapz(x = yWorkGrid, y = err)
        }
        AISE[ybwIdx,tbwIdx] <- mean(ise)
      }
    }
    if ( ncores > 1 ) {
      stopCluster(cl)
    }
    optIdx <- which(AISE == min(AISE), arr.ind = TRUE)
    ybw <- ybwCand[optIdx[1]]
    tbw <- tbwCand[optIdx[2]]
  }
  
  if (bwchoice) {
    return(list(ybw = ybw, tbw = tbw))
  } else {
    yin <- unlist(Ly)
    tin <- unlist(Lt)
    res <- cdfsm(tin = tin, yin = yin, tout = tWorkGrid, yout = yWorkGrid, tbw = tbw, ybw = ybw, K = K, H = H)
    res <- sapply(seq_along(tWorkGrid), function(j) approx(x = yWorkGrid, y = res[,j], xout = yWorkMat[,j])$y)
    optns <- list(tWorkGrid = tWorkGrid, yWorkMat = yWorkMat, tbw = tbw, ybw = ybw, Kty = Kty, Hty = Hty,
                  training = training, ncores = ncores)
    return(list(rank = res, optns = optns))
  }
}
