#' Local polynomial smoothing
#' @param xin A vector of length \eqn{n} with measurement points.
#' @param yin A vector of length \eqn{n} with measurement values.
#' @param xout A vector of length \eqn{m} with output measurement points.
#' @param kty A character holding the kernel type; default: \code{'epan'}.
#' @param bw A scalar holding the bandwidth.
#' @param npoly The order of polynomials used; e.g., \code{1} for local linear smoothing.
#' @return An \eqn{m}-by-\eqn{p} matrix of which columns hold the estimates of the function 
#' and its derivatives up to order \eqn{p} evaluated on \code{xout}, where \eqn{p} equals (\code{npoly-1}).
#' @export
lps <- function(xin,yin,xout,kty='epan',bw,npoly=1) {
  ker <- kerfctn( kernel_type = kty, fctn_type = "pdf" )
  
  n <- length(xin)
  if (length(bw) > 1) {
    bw <- bw[1]
    warning("bw holds more than one elements; only the first value is used.")
  }
  yout <- matrix(numeric(length(xout)*(npoly+1)),ncol=npoly+1)
  for(k in 1:length(xout)) {
    xstar <- sapply(0:npoly, function(j) {
      return((xin-xout[k])^j)
    })
    wsr <- sqrt(ker((xout[k]-xin)/bw))
    xstar <- wsr * xstar
    temp <- svd(t(xstar) %*% xstar)
    lambda <- temp$d[abs(temp$d)>1e-5]
    inv <- temp$v %*% diag(c(1/lambda, rep(0, npoly+1-length(lambda)))) %*% t(temp$u)
    yout[k,] <- inv %*% t(xstar) %*% (wsr * yin)
  }
  return(yout[,1:npoly])
}

