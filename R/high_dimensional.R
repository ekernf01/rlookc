
#' Create Gaussian knockoffs with dense covariance in a large p, small n setting.
#'
#' @param X input features
#' @param rho tuning parameter with range 0<rho<1. Larger is expected to work better.
#' @param lambda degree of shrinkage: 0 returns the sample correlations and 1
#' @param output_type What to return: either knockoffs themselves or parameters
#' for downstream use in e.g. unit tests or leave-one-outs.
#' @param silent Do not issue warnings.

#' returns correlations of all 0.
#' You can set this for optimal MSE using corpcor::estimate_lambda.
#' @export
#'
createHighDimensionalKnockoffs = function(X, rho = 0.9, lambda = NULL, silent =F,
                                          output_type = c("knockoffs", "parameters"), seed = NULL ){
  if(!is.null(seed)){
    set.seed(seed)
  }
  # Handle too many or too few observations
  stopifnot("Use this only when p>>n"=nrow(X)<ncol(X))
  if(nrow(X)<=3){
    if(!silent){
      warning("Too few observations to estimate mean and covariance. Returning data verbatim.")
    }
    return(X)
  }
  # center to zero mean and scale to unit variance (it goes back at the end of the function)
  mu = apply(X, 2, mean)
  X = sweep(X, 2, mu, FUN = "-")
  standard_deviations = apply(X, 2, sd)
  nzsd = standard_deviations>0
  # Handle zero-variance features
  X[,nzsd] = sweep(X[,nzsd], 2, standard_deviations[nzsd], FUN = "/")
  n = nrow(X)
  svdX = svd(X, nv = nrow(X), nu = 0)
  # Optimal shrinkage
  if(is.null(lambda)){
    lambda = corpcor::estimate.lambda(X, verbose = F)
  }
  # basic_transform goes singular values to eigs of Sigma
  basic_transform = function(ti) { (1 - lambda)*ti^2/(n-1)  + lambda }
  # additional_transform goes from eigs of Sigma to user defined. Use
  # 1/x for Sigmainv or the ugly thing provided below for the conditional
  # covariance C = 2S-S\Sigma_inv S .
  multiplyBySigma = function(Z, additional_transform = function(x) x ){
    Zx    = (Z  %*% svdX$v) %*% Matrix::Diagonal( x = additional_transform( basic_transform( svdX$d ) ) )
    Zx    = Zx %*% t(svdX$v)
    Zperp = Z - ( Z %*% svdX$v ) %*% t( svdX$v )
    Zperp = Zperp * additional_transform( basic_transform( 0 ) )
    return( Zx + Zperp )
  }
  # Conditional mean assuming marginal mean is 0
  knockoff_mean = X - multiplyBySigma( X, additional_transform = function(x) 1/x ) * (2*rho*lambda)
  # Undo centering
  knockoff_mean = sweep(knockoff_mean, 2, mu, FUN = "+")
  # The covariance calculations are finalized inside of the
  # conditional depending on the type of output desired.
  transform_eigenvectors_from_sigma_to_sqrt_c = function(x) sqrt(4*rho*lambda - 4*rho^2*lambda^2/x)
  output_type = match.arg(output_type)
  if(output_type == "knockoffs"){
    knockoff_random = matrix(rnorm(prod(dim(X))), nrow = nrow(X))
    knockoff_random = multiplyBySigma(knockoff_random,
                                      additional_transform = transform_eigenvectors_from_sigma_to_sqrt_c)
    # Restore previous SD
    knockoff_random = sweep(knockoff_random, 2, standard_deviations, FUN = "*")
    return(as.matrix(knockoff_mean + knockoff_random))
  } else if(output_type == "parameters"){
    return( list(
      knockoff_mean=knockoff_mean,
      multiplyBySigma=multiplyBySigma,
      transform_eigenvectors_from_sigma_to_sqrt_c = transform_eigenvectors_from_sigma_to_sqrt_c,
      standard_deviations = standard_deviations
    ) )
  }
}
