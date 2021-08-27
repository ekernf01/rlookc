
#' Compute S, mean, covariance, and realizations of Gaussian knockoffs.
#'
#' This function is copied and lightly modified from the function \code{create.gaussian}
#' in Matteo Sesia's R package 'knockoff', licensed
#' under GPLv3 and published at https://github.com/msesia/knockoff-filter/tree/master/R .
#' I needed to modify it to return some internals and handle groups of variables.
#'
#' Public use is not recommended when "output_type" is anything but "knockoffs",
#' as results are poorly labeled.
#' This was done in order to make the interface similar to generateLooks
#' while still providing the correct functionality for internal use.
#' Users should typically prefer generateLooks or Matteo Sesia's \code{knockoff::create.gaussian}.
#'
#' @param X @param mu  @param sigma  @param method @param groups  @param diag_s  @param output_type See \code{?generateLooks}.
#' @param num_realizations Since the knockoff procedure is random, it can be useful to generate many instances,
#' for instance in simulation studies. If >1, this returns a list of multiple realizations.
#' @param seeds Random seed. Set this for repeatable results.
#' @export
#'
computeGaussianKnockoffs = function (
  X,
  mu = colMeans(X),
  Sigma = cov(X),
  groups = NULL,
  method = c("asdp", "sdp", "equi", "group"),
  diag_s = NULL,
  output_type = c("knockoffs", "knockoffs_compact", "statistics", "parameters"),
  num_realizations = 1
){
  # Input checking for S solver (Same as Sesia's code)
  method = match.arg(method)
  if ((nrow(Sigma) <= 500) && method == "asdp") {
    method = "sdp"
  }
  # Input checking for dimensions
  stopifnot("Input dimension mismatch for mu"=length(mu) %in% c(1, ncol(X)))
  stopifnot("Input dimension mismatch for Sigma"=nrow(Sigma)==ncol(X))
  stopifnot("Input dimension mismatch for Sigma"=ncol(Sigma)==ncol(X))
  # Input checking for groups
  if( method == "group" ){
    if( is.null( groups ) ){
      stop("If method is 'group', variable groups must be specified.\n")
    }
    checkGroups(X, groups)
  }
  if( method != "group" & !is.null( groups ) ){
    stop("Input unclear. If variable groups are specified, method must be 'group'.\n")
  }
  # Input checking for output type
  if( match.arg(output_type) == "knockoffs_compact" ){
    stop("This function cannot produce a compact representation.\n")
  }
  if( num_realizations != 1 & match.arg(output_type) != "knockoffs" ){
    stop("For multiple realizations, only `output_type='knockoffs'` is supported.\n")
  }


  # This part is identical to Sesia's original except the added "group" option
  if (is.null(diag_s)) {
    diag_s =
      switch(
        match.arg(method),
        equi = knockoff::create.solve_equi(Sigma),
        sdp  = knockoff::create.solve_sdp(Sigma),
        asdp = knockoff::create.solve_asdp(Sigma),
        group = solveGroupEqui(Sigma, groups)
      )
  }
  if( is.vector( diag_s ) ){
    diag_s = diag( diag_s, length( diag_s ) )
  }
  diag_s = as.matrix(diag_s)
  if( all( diag_s == 0 ) ) {
    warning( "Knockoffs cannot be constructed except as copies of X. Knockoffs will have no power." )
    return( X )
  }
  # Mostly the same as the original code with a few more named intermediates that will be needed later
  SigmaInv = solve(Sigma)
  SigmaInv_s = SigmaInv %*% diag_s
  mu_ko = X - sweep(X, 2, mu, "-") %*% SigmaInv_s
  Sigma_ko = 2 * diag_s - diag_s %*% SigmaInv_s
  R = chol(Sigma_ko)
  # For simulation studies, make multiple realizations without redoing all the computation
  if(num_realizations > 1){
    X_ko = list()
    for(i in seq(num_realizations)){
      X_ko[[i]] = mu_ko + matrix(rnorm(ncol(X) * nrow(X)), nrow(X)) %*% R
    }
  } else {
    X_ko = mu_ko + matrix(rnorm(ncol(X) * nrow(X)), nrow(X)) %*% R
  }
  switch(
    match.arg(output_type),
    knockoffs = X_ko,
    statistics = X_ko, # Not what the user wants (yet) but it's what the wrapper needs to get there
    pearson = X_ko,    # ^ Same
    parameters = list(ko_mean = mu_ko, ko_sqrt_covariance = R, S = diag_s, SigmaInv = SigmaInv, knockoffs = X_ko)
  )
}

#' Generate leave-one-out knockoffs for each variable in X slowly and simply. Suitable for small problems or software tests.
#'
#' @details see ?generateLooks .
generateLooksSlow = function(
  X, mu, Sigma,
  method = c("asdp", "sdp", "equi", "group"),
  diag_s = NULL,
  groups = NULL,
  statistic = knockoff::stat.glmnet_coefdiff,
  vars_to_omit = 1:ncol(X),
  output_type = c("knockoffs", "knockoffs_compact", "statistics", "parameters", "pearson"),
  ...
){
  if(length(mu)==1){
    mu = rep(mu, ncol(X))
  }
  if(is.null(dim(diag_s))){
    diag_s = diag(diag_s)
  }
  do_one = function(k){
    ko = computeGaussianKnockoffs(X[,-k],
                                  mu[-k],
                                  Sigma[-k,-k],
                                  method = method,
                                  groups = groups %>% removeKFromGroups(k),
                                  output_type = output_type,
                                  diag_s = diag_s[-k,-k])
    if(output_type=="pearson"){
      return(marginalScreen(X[,-k], ko, y = X[,k], ...))
    } else if(output_type=="statistics"){
      if(is.null(statistic)){
        stop("'statistic' must be provided if output_type = 'statistics'. Consult options from the 'knockoff' package.")
      }
      return(statistic(X[,-k], ko, y = X[,k], ...))
    } else {
      return( ko )
    }
  }
  lapply(vars_to_omit, do_one)
}

#' Generate and use leave-one-out knockoffs for each variable in X.
#'
#' @param X n-by-p matrix of original variables.
#' @param mu vector of length p, indicating the mean parameter of the Gaussian model for \eqn{X}.
#' @param Sigma p-by-p covariance matrix for the Gaussian model of \eqn{X}.
#' @param method either "equi", "sdp" or "asdp" (default: "asdp").
#' This determines the method that will be used to minimize the correlation between the original variables and the knockoffs.
#' @param diag_s Square PxP matrix or a vector of length P, containing the pre-computed matrix S where
#' \deqn{s_{jk} = cov(x_k, x_j) - cov(x_k, \tilde x_j)}.
#' This will be computed according to \code{method}, if not supplied.
#' @param statistic function used to assess variable importance. In addition to the options from
#' the \code{knockoff} package, you can enter the string \code{"pearson"}. This will take advantage
#' of the low-rank knockoff updates computed by generateLooks. This is only used if \code{output_type}
#' is "statistics".
#' @param output_type
#' If "knockoffs" (default), result is a list of knockoff realizations. List is of the same length as vars_to_omit.
#' This can still be faster for generateLooks than for generateLooksSlow, as low-rank updates are used internally.
#' If "knockoffs_compact", a low-rank, low-memory update for each knockoff is returned, along with the precomputed
#' knockoffs that they modify. This option is not available for generateLooksSlow.
#' If "parameters", mean and covariance for each is returned.
#' This is for internal use.
#' If "statistics", knockoffs themselves are not returned; rather, test statistics are returned.
#' See parameter \code{statistic.}
#' This will err if \code{statistic} is set to NULL.
#' If "pearson", statistics are returned like with \code{output_type="statistic"}.
#' It's fast but rigid: the \code{statistic} arg is ignored and Pearson correlations are used.
#' @param ... Passed to statistic
#' @return See parameter \code{return_type}.
#' @export
generateLooks = function(
  X, mu, Sigma,
  method = c("asdp", "sdp", "equi", "group"),
  groups = NULL,
  diag_s = NULL,
  vars_to_omit = 1:ncol(X),
  statistic = knockoff::stat.glmnet_coefdiff,
  output_type = c("knockoffs", "knockoffs_compact", "statistics", "parameters", "pearson"),
  ...
){
  # Check and clean up input
  if(length(mu)==1){
    mu = rep(mu, ncol(X))
  }
  if(is.vector(diag_s)){
    diag_s = diag(diag_s)
  }
  output_type = match.arg(output_type)
  # This calls a modification of Matteo Sesia's code that reveals some internals for downstream use.
  precomputed_quantities = computeGaussianKnockoffs(
    X = X,
    mu = mu,
    Sigma = Sigma,
    method = method,
    groups = groups,
    diag_s = diag_s,
    output_type = "parameters"
  )

  # Prepare to do the rest of the work via rank-1 updates.
  # I used drop=F to make sure I don't do an inner product when an outer product is needed.
  get_compact_updates = function(k){
    if(length(k)!=1){
      stop("Current implementation can only omit one variable at a time, sorry.\n")
    }
    g = precomputed_quantities$SigmaInv[+k, -k, drop = F]
    h = precomputed_quantities$SigmaInv[+k, +k]
    sqrt_h = sqrt(h)
    list(
      # E(ko[,-k] | X[,-k]) - E(ko[,-k] | X) -- line 14 in derivation
      mean_update_left =  X[,-k] %*% t(g) / h + X[,k, drop = F],
      mean_update_right = g %*% precomputed_quantities$S[-k,-k] + h*precomputed_quantities$S[k, -k, drop = F],
      # Sqrt of cov(ko[,-k] | X[,-k]) - cov(ko[,-k] | X)  -- line 24 in derivation
      sqrt_cov_update =
        (g/sqrt_h) %*% precomputed_quantities$S[-k, -k, drop = F] + # always used
        sqrt_h * precomputed_quantities$S[+k, -k, drop = F], # 0 if S diagonal (no grouping)
      # Derivation relies on
      #
      #    cov(X + RZ) = cov(X) + RR^T
      #
      # where Z is IID standard normal.
      # It helps to have 'Z' stored with the fixed updates, not generated de novo upon each use, so that
      # the exact same knockoffs are always re-assembled.
      random_update = matrix( rnorm( nrow( X ) ), ncol = 1 )
    )
  }
  updates = lapply(vars_to_omit, get_compact_updates)
  # Send'em out on the cheap if possible
  if ( output_type == "pearson"){
    stop("Sorry, this is not implemented yet.")
  } else if ( output_type == "knockoffs_compact" ){
    return(
      list(
        knockoffs = precomputed_quantities$knockoffs,
        updates = updates,
        groups = groups,
        vars_to_omit = vars_to_omit
      )
    )
  } else if ( output_type == "parameters"){
    return( lapply(vars_to_omit, function(k) formLookParameters(precomputed_quantities, vars_to_omit, updates, k) ) )
  } else if(output_type == "knockoffs"){
    return( formAllLooks(precomputed_quantities$knockoffs, vars_to_omit, updates) )
  } else if(output_type == "statistics"){
    if( is.null( statistic ) ){
      stop("'statistic' must be provided if output_type = 'statistics'. Consult options from the 'knockoff' package.")
    }
    return( formAllLooks(precomputed_quantities$knockoffs, vars_to_omit, updates, statistic, X, ...) )
  }
}

#' Extract one update from the low-rank representation, while handling a tricky case properly.
#'
getUpdateK = function(k, vars_to_omit, updates){
  correct_index = which(k==vars_to_omit) # Handle case when vars_to_omit is not 1, 2, ... P
  updates[[correct_index]]
}

#' Given the low-rank representations, update knockoffs to omit each variable.
#'
#' @param k variable to omit.
#' @param statistic Optional but useful for memory efficiency: instead of returning all knockoffs, compute statistics using this function, and return those instead.
#' @param updates @param knockoffs @param vars_to_omit
#' Inputs should be from \code{loadCompactLooks} or from \code{generateLooks(..., output_type = 'knockoffs_compact'}.)
#' Those functions return a list with the same names as the necessary args.
#' @export
#'
formAllLooks = function(knockoffs, vars_to_omit, updates, statistic = NULL, X = NULL, ...){
  if(is.null(statistic)){
    return( lapply( vars_to_omit, function(k) formOneLook(knockoffs, vars_to_omit, updates, k) ) )
  } else {
    stopifnot("Pass original data to formAllLooks if you want test statistics as output.\n"=!is.null(X))
    return( lapply( vars_to_omit, function(k) statistic( X[,-k], formOneLook(knockoffs, vars_to_omit, updates, k), y = X[,k], ... ) ) )
  }
}

#' Given the low-rank representations, update knockoffs to omit one variable.
#'
#' @param k variable to omit.
#' @param updates @param knockoffs @param vars_to_omit
#' Inputs should be from \code{loadCompactLooks} or from \code{generateLooks(..., output_type = 'knockoffs_compact'}.)
#' Those functions return a list with the same names as the necessary args.
#' @export
#'
formOneLook = function(knockoffs, vars_to_omit, updates, k){
  one_update = getUpdateK(k, vars_to_omit, updates)
  with(one_update,
       knockoffs[,-k] +
         getMeanUpdate(one_update) +
         getRandomUpdate(one_update, n_obs = nrow(knockoffs))
  )
}

#' Return mean and covariance used to sample a look. This is only used for checking the math.
#'
formLookParameters = function(precomputed_quantities,  vars_to_omit, updates, k){
  one_update = getUpdateK(k, vars_to_omit, updates)
  list(
    ko_mean                 = precomputed_quantities$ko_mean[           ,-k] + getMeanUpdate( one_update ),
    ko_covariance = crossprod(precomputed_quantities$ko_sqrt_covariance[,-k]) + getCovarianceUpdate( one_update ),
    ko_mean_naive           = precomputed_quantities$ko_mean[           ,-k],
    ko_cov_naive  = crossprod(precomputed_quantities$ko_sqrt_covariance[,-k]),
    S        = precomputed_quantities$S[-k,-k],
    SigmaInv = precomputed_quantities$SigmaInv[-k,-k]
  )
}

#' Output can be added to knockoffs (or their covariance) to correct for removal of a variable.
#'
getMeanUpdate = function(updates){
  with(updates,
       mean_update_left %*% mean_update_right
  )
}
getRandomUpdate = function(updates, n_obs){
  updates$random_update %*% updates$sqrt_cov_update
}

# Used only for checking the math -- not needed for knockoff generation
getCovarianceUpdate = function(updates){
  crossprod(updates$sqrt_cov_update)
}
