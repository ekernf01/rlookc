
#' Create Gaussian knockoffs with dense covariance in a large p, small n setting.
#'
#' @param X input features
#' @param rho tuning parameter with range 0<rho<1. Larger is expected to work better.
#' @param lambda degree of shrinkage: 0 returns the sample correlations and 1
#' returns correlations of all 0.
#' Set this optimally using corpcor::estimate_lambda.
highDimensionalKnockoffs = function(X, rho = 0.9, lambda = 0.1){
  stopifnot("Use this only when p>>n"=nrow(X)>ncol(X))
  # center to zero mean and scale to unit variance (it goes back at the end of the function)
  mu = apply(X, 2, mean)
  X = sweep(X, 2, mu, FUN = "-")
  standard_deviations = apply(X, 2, sd)
  X = sweep(X, 2, standard_deviations, FUN = "/")

  # svd of reweighted data matrix
  n = nrow(X)
  W = X
  svdW = svd(W, nv = nrow(X), nu = 0)
  # basic_transform goes singular values to eigs of Sigma
  basic_transform = function(ti) { (1 - lambda)*ti^2/(n-1)  + lambda }
  # additional_transform goes from eigs of Sigma to user defined
  # 1/x for Sigma to Sigmainv or the ugly thing for 2S-S\Sigma_inv S
  multiplywithCustomEigs = function(Z, additional_transform ){
    Zx    = Z  %*% svdW$v %*% Matrix::Diagonal( x = additional_transform( basic_transform( svdW$d ) ) )
    Zx    = Zx %*% t(svdW$v)
    Zperp = Z - ( Z %*% svdW$v ) %*% t( svdW$v )
    Zperp = Zperp * additional_transform( basic_transform( 0 ) )
    return( Zx + Zperp )
  }
  # mean
  knockoff_mean = X - multiplywithCustomEigs( X, additional_transform = function(x) 1/x ) * (2*rho*lambda)
  knockoff_random = matrix(rnorm(prod(dim(X))), nrow = nrow(X))
  go_from_sigma_to_sqrt_c = function(x) sqrt(4*rho*lambda - 4*rho^2*lambda^2/x)
  knockoff_random = multiplywithCustomEigs(knockoff_random, additional_transform = go_from_sigma_to_sqrt_c)
  knockoffs = knockoff_mean + knockoff_random
  # Restore previous mean and SD
  knockoffs = sweep(knockoffs, 2, standard_deviations, FUN = "*")
  knockoffs = sweep(knockoffs, 2, mu, FUN = "+")
  return(knockoffs)
}
