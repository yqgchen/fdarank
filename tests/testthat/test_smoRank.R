
require(testthat)

test_that('Check it works if no bandwidths are pre-specified', {
  n <- 10 #sample size
  m <- 51 # number of time points
  K <- 20 # number of eigenfunctions
  
  tgrid <- seq( 0, 1, length.out = m )
  egnfctn <- function ( t, k ) {
    j <- k %/% 2
    if ( k %% 2 > 0 ) {
      sqrt(2) * sin( j * pi * t )
    } else {
      sqrt(2) * cos( j * pi * t )
    }
  }
  egnval <- 1 / ( seq_len(K)^2 )
  
  set.seed(1)
  FPCs <- matrix( rnorm( n * K ), nrow = n ) * matrix( sqrt(egnval), nrow = n, ncol = K, byrow = TRUE )
  egnfctnMat <- t( sapply( seq_len(K), egnfctn, t = tgrid ) )
  data <- FPCs %*% egnfctnMat
  
  sd <- sqrt( colSums( egnfctnMat^2 * egnval ) )
  true <- pnorm( data / matrix( sd, nrow = n, ncol = m, byrow = TRUE ) )
  
  res <- smoRank(
    Ly = lapply( seq_len(n), function (i) data[i,] ),
    Lt = rep( list(tgrid), n ),
    optns = list( tWorkGrid = tgrid )
  )
  est <- res$rank
  
  expect_lt( mean( abs( true - est ) ), 0.08 )
})

test_that('Check it works if bandwidths are given', {
  n <- 10 #sample size
  m <- 51 # number of time points
  K <- 20 # number of eigenfunctions
  
  tgrid <- seq( 0, 1, length.out = m )
  egnfctn <- function ( t, k ) {
    j <- k %/% 2
    if ( k %% 2 > 0 ) {
      sqrt(2) * sin( j * pi * t )
    } else {
      sqrt(2) * cos( j * pi * t )
    }
  }
  egnval <- 1 / ( seq_len(K)^2 )
  
  set.seed(1)
  FPCs <- matrix( rnorm( n * K ), nrow = n ) * matrix( sqrt(egnval), nrow = n, ncol = K, byrow = TRUE )
  egnfctnMat <- t( sapply( seq_len(K), egnfctn, t = tgrid ) )
  data <- FPCs %*% egnfctnMat
  
  sd <- sqrt( colSums( egnfctnMat^2 * egnval ) )
  true <- pnorm( data / matrix( sd, nrow = n, ncol = m, byrow = TRUE ) )
  
  res <- smoRank(
    Ly = lapply( seq_len(n), function (i) data[i,] ),
    Lt = rep( list(tgrid), n ),
    optns = list( tWorkGrid = tgrid, tbw = 2.5/(m-1), ybw = 1 )
  )
  est <- res$rank
  
  expect_lt( mean( abs( true - est ) ), 0.08 )
})
