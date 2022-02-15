
require(testthat)

test_that('Check computed empirical ranks on target', {
  n <- 10 #sample size
  m <- 51 # number of time points
  K <- 10 # number of eigenfunctions
  
  tgrid <- seq( 0, 1, length.out = m )
  egnfctn <- function ( t, k ) {
    j <- k %/% 2
    if ( k %% 2 > 0 ) {
      sqrt(2) * sin( j * pi * t )
    } else {
      sqrt(2) * cos( j * pi * t )
    }
  }
  FPCs <- matrix( rnorm( n * K ), nrow = n )
  data <- FPCs %*% t( sapply( seq_len(K), egnfctn, t = tgrid ) )
  est <- empRank( Ly = lapply( seq_len(n), function (i) data[i,] ), Lt = rep( list(tgrid), n ), workGrid = tgrid )
  expected <- ( apply( data, 2, rank ) - 1 ) / n
  expect_equal( all( abs( expected - est ) == 0 ), TRUE )
})
