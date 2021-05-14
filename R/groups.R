#' Make sure variable groups match data matrix (or any other matrix with one variable per column) 
#'
checkGroups = function(X, groups){
  vars_expected = seq( ncol( X ) )
  vars_in_groups = sort( Reduce( c, groups ) )
  if( !all( vars_in_groups == vars_expected ) ){
    browser()
    stop("Each variable must be in exactly one group. Use numbers, not names.\n")
  }
}
  
  
#' Given one knockoff test stat per variable, return one per group.
#'
aggregateStats = function(stats, groups, FUN = mean){
  if(is.null(dim(stats))){
    stats = matrix(stats, nrow = 1)
  }
  checkGroups(stats, groups)
  stats_agg = matrix(0, ncol = length(groups), nrow = nrow(stats))
  for(gi in seq_along(groups)){
    stats_agg[,gi] = apply(stats[ , groups[[gi]], drop = F ], 1, FUN)
  }
  stats_agg
}

#' Compute the block-diagonal matrix S for grouped model-X knockoffs using "equi" extension in Dai & Barber 2016.
#'
#' Dai, R., & Barber, R. (2016, June). The knockoff filter for FDR control in
#' group-sparse and multitask regression. In International Conference on
#' Machine Learning (pp. 1851-1859). PMLR.
#'
#' See equation 11 and material immediately following it.
#'
solveGroupEqui = function(Sigma, groups, do_fast = T){
  D = S = Matrix::Matrix(0, nrow = nrow(Sigma), ncol = ncol(Sigma))
  sqrt_inv = function(X, power = -0.5){
    eigenstuff = eigen(X, symmetric = T)
    if(any(eigenstuff$values <= 0)){
      stop(paste0(
        "Sorry, this grouped knockoff implementation cannot handle singular covariance matrices within groups.\n",
        "You'll need to prune one or more of these linearly redundant features: \n", paste0(g, collapse = " "), "\n")
      )
    }
    eigenstuff$vectors %*% diag(eigenstuff$values ^ power) %*% t(eigenstuff$vectors)
  }
  # Test of sqrt inv closure
  # M = matrix(rnorm(160), ncol = 4)
  # M = t(M) %*% M
  # round(sqrt_inv(M, power = -1) %*% M)
  # round(sqrt_inv(M, power = -0.5) %*% M %*% sqrt_inv(M, power = -0.5))
  # round(sqrt_inv(M, power = 1) - M)
  for(g in groups){
    S[g,g] = Sigma[g,g]
    D[g,g] = sqrt_inv(Sigma[g,g])
  }
  M = D %*% Sigma %*% D
  trace_M = sum(Matrix::diag(M))
  if(do_fast){
    # Get the smallest eigenvalue of M in a hopefully fast and non-horrible way
    # Via largest eigenvalue of tr(M)*I - M.
    min_eigenvalue = trace_M - irlba::partial_eigen(n = 1, trace_M * diag(nrow(M)) - M, symmetric = T )$values
  } else {
    min_eigenvalue = min(eigen(M)$values)
  }
  S = min(min_eigenvalue, 1)*S
  return(S)
}
