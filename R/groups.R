#' Make sure variable groups match data matrix (or any other matrix with one variable per column)
#'
#' @export
#'
checkGroups = function(X, groups){
  vars_expected = seq( ncol( X ) )
  vars_in_groups = sort( Reduce( c, groups ) )
  if( !all( vars_in_groups == vars_expected ) ){
    stop("Each variable must be in exactly one group. Use numbers, not names.\n")
  }
}

#' Remove one variable from a list of groups, then shift everything down.
#'
removeKFromGroups = function( groups, k ){
  if(length(groups)==0){
    groups = NULL
  } else {
    # Remove index k wherever it occurs
    groups = lapply(groups, setdiff, k)
    # If that leaves an empty group, remove it
    groups = groups[sapply(groups, length)>0]
    # Any variable >k now appears one slot to the left
    groups = lapply(groups, function(x) x - ifelse(x>k, 1, 0))
  }
  groups
}

#' Deprecated. See stats__aggregateStats
#'
#' @export
aggregateStats = function(...){ stats__aggregateStats(...) }

#' Given one knockoff test stat per variable, return one per group.
#'
#' @export
#'
stats__aggregateStats = function(stats, groups, FUN = mean){
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

#' Get a diagonal matrix from input of unknown shape
#'
safe_diag = function(x) {
  stopifnot(is.vector(x))
  if(length(x)>1){
    return(diag(x))
  } else {
    return(x)
  }
}

#' Get a symmetric inverse square root of an spd matrix via eigendecomposition
#'
sqrt_inv = function(X, power = -0.5){
  eigenstuff = eigen(X, symmetric = T)
  if(any(eigenstuff$values <= 0)){
    stop(paste0(
      "Sorry, this grouped knockoff implementation cannot handle singular covariance matrices within groups.\n",
      "You'll need to prune one or more of these linearly redundant features: \n", paste0(g, collapse = " "), "\n")
    )
  }
  eigenstuff$vectors %*% safe_diag(eigenstuff$values ^ power) %*% t(eigenstuff$vectors)
}

# Test of sqrt inv
# M = matrix(rnorm(160), ncol = 4)
# M = t(M) %*% M
# round(sqrt_inv(M, power = -1) %*% M)
# round(sqrt_inv(M, power = -0.5) %*% M %*% sqrt_inv(M, power = -0.5))
# round(sqrt_inv(M, power = 1) - M)

#' Deprecated. See stats__aggregateStats
#'
#' @export
solveGroupEqui = function(...){ create__solveGroupEqui(...) }


#' Compute the block-diagonal matrix S for grouped model-X knockoffs using "equi" extension in Dai & Barber 2016.
#'
#' Dai, R., & Barber, R. (2016, June). The knockoff filter for FDR control in
#' group-sparse and multitask regression. In International Conference on
#' Machine Learning (pp. 1851-1859). PMLR.
#'
#' See equation 11 and material immediately following it.
#'
create__solveGroupEqui = function(Sigma, groups, do_fast = T){
  stopifnot(isSymmetric(Sigma))
  # I'm not sure whether this assumes unit-variance data, so I scale it as such and rescale it at the end.
  original_variances = diag(Sigma)
  Sigma = diag(1/sqrt(original_variances)) %*% Sigma %*% diag(1/sqrt(original_variances))

  # The actual calculation of S
  D = S = Matrix::Matrix(0, nrow = nrow(Sigma), ncol = ncol(Sigma))
  for(g in groups){
    S[g,g] = Sigma[g,g]
    D[g,g] = sqrt_inv(Sigma[g,g])
  }
  M = D %*% Sigma %*% D
  trace_M = sum(Matrix::diag(M))
  min_eigenvalue = RSpectra::eigs_sym(M, 1, which = "SA", opts = list(retvec = FALSE, maxitr = 1e+05, tol = 1e-08))$values
  S = min(2*min_eigenvalue, 1)*S

  # Make sure it's not too close to causing a 0 eigenvalue in G
  S = S*0.99
  S = as.matrix(S)

  # Restore the original scaling per coordinate
  S = diag(sqrt(original_variances)) %*% S %*% diag(sqrt(original_variances))
  return(S)
}
