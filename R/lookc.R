

#' Compute S, mean, covariance, and realizations of Gaussian knockoffs.
#'
#' This function is copied and lightly modified from the function \code{create.gaussian}
#' in Matteo Sesia's R package 'knockoff', licensed
#' under GPLv3 and published at https://github.com/msesia/knockoff-filter/tree/master/R .
#' I needed to modify it to return some internals and handle groups of variables.
#'
#' Results of "output_type" are downright mislabeled -- that's just a holdover from generateLooks and it's not meant for public use.
#'
computeGaussianKnockoffs = function (
  X,
  mu,
  Sigma,
  groups = NULL,
  method = c("asdp", "sdp", "equi", "group"),
  diag_s = NULL,
  output_type = c("knockoffs", "knockoffs_compact", "statistics", "parameters")
){
  # Input checking for S solver (Same as Sesia's code)
  method = match.arg(method)
  if ((nrow(Sigma) <= 500) && method == "asdp") {
    method = "sdp"
  }
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
  diag_s %<>% as.matrix()
  if( all( diag_s == 0 ) ) {
    warning( "The conditional knockoff covariance matrix is not positive definite. Knockoffs will have no power." )
    return( X )
  }
  # Mostly the same as the original with a few more named intermediates that will be needed later
  SigmaInv = solve(Sigma)
  SigmaInv_s = SigmaInv %*% diag_s
  mu_ko = X - sweep(X, 2, mu, "-") %*% SigmaInv_s
  Sigma_ko = 2 * diag_s - diag_s %*% SigmaInv_s
  R = chol(Sigma_ko)
  X_ko = mu_ko + matrix(rnorm(ncol(X) * nrow(X)), nrow(X)) %*% R
  if( match.arg(output_type) == "knockoffs_compact" ){
    stop("This function cannot produce a compact representation.")
  }
  switch(
    match.arg(output_type),
    knockoffs = X_ko,
    statistics = X_ko, # Not what the user wants (yet) but it's what the wrapper needs to get there
    pearson = X_ko,    # ^ Same
    parameters = list(ko_mean = mu_ko, ko_sqrt_covariance = R, S = diag_s, SigmaInv = SigmaInv, knockoffs = X_ko)
  )
}

#' Helper function for simple marginal screening via the pearson correlation.
#'
do_pearson_screen = function(X, ko, y){
  rbind(as.vector(cor(X, y)), as.vector(cor(ko, y))) %>% set_rownames(c("X", "knockoff"))
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
      return(do_pearson_screen(X[,-k], ko, y = X[,k], ...))
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
  if(length(mu)==1){
    mu = rep(mu, ncol(X))
  }
  if(is.vector(diag_s)){
    diag_s = diag(diag_s)
  }
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

  # Do the rest of the work via rank-1 updates.
  get_compact_updates = function(k){
    g = precomputed_quantities$SigmaInv[+k, -k, drop = F]
    h = precomputed_quantities$SigmaInv[+k, +k]
    sqrt_h = sqrt(h)
    list(
      # Mean of P(ko | X) -- line 3 in derivation
      mean_update1_left = X[,+k, drop = F] + X[,-k, drop = F] %*% t(g/h),
      mean_update1_right = g %*% precomputed_quantities$S[-k,-k, drop = F],
      # These extra terms are used only for grouped variables (when S is not diagonal)
      mean_update2_left = -( X[,-k] %*% t(g) + X[,+k]*h ),
      mean_update2_right = precomputed_quantities$S[+k,-k, drop = F],
      # Sqrt of cov(ko | X) -- lines 8,18 in derivation
      sqrt_cov_update =
        (g/sqrt_h) %*% precomputed_quantities$S[-k, -k, drop = F] + # always used
        sqrt_h * precomputed_quantities$S[+k, -k, drop = F] # 0 if S diagonal (no grouping)
    )
  }
  all_updates = lapply(vars_to_omit, get_compact_updates)

  # Send'em out on the cheap if possible
  if ( match.arg(output_type) == "pearson"){
    stop("Sorry, this is not implemented yet.")
  }
  if ( match.arg(output_type) == "knockoffs_compact"){
    return(
      list(
        knockoffs = precomputed_quantities$knockoffs,
        updates = all_updates
      )
    )
  }

  # Otherwise, return nicer but more expensive stuff as requested
  assemble_output = function(k){
    ko = updateKnockoffs(precomputed_quantities$knockoffs, all_updates[[which(k==vars_to_omit)]], k)
    if(output_type == "knockoffs"){
      return( ko )
    } else if(output_type == "statistics"){
      if(is.null(statistic)){
        stop("'statistic' must be provided if output_type = 'statistics'. Consult options from the 'knockoff' package.")
      }
      return( statistic(X[,-k], ko, y = X[,k], ...) )
    } else if ( output_type == "parameters"){
      updates = all_updates[[which(vars_to_omit==k)]]
      return(
        list(
          ko_mean                 = precomputed_quantities$ko_mean[           ,-k] + getMeanUpdate( updates ),
          ko_covariance = crossprod(precomputed_quantities$ko_sqrt_covariance[,-k]) + getCovarianceUpdate( updates ),
          ko_mean_naive           = precomputed_quantities$ko_mean[           ,-k],
          ko_cov_naive  = crossprod(precomputed_quantities$ko_sqrt_covariance[,-k]),
          S        = precomputed_quantities$S,
          SigmaInv = precomputed_quantities$SigmaInv
        )
      )
    }
  }
  return( lapply( vars_to_omit, assemble_output ) )
}

#' Given the low-rank updates from generateLooks, update knockoffs to omit variable k.
#'
#' @param knockoffs Knockoffs computed with variable k included.
#' @param k variable to omit.
#' @param updates updates from generateLooks with output_type = "knockoffs_compact"
#' @export
updateKnockoffs = function(knockoffs, updates, k){
  with(updates,
       knockoffs[,-k] +
         getMeanUpdate(updates) +
         getRandomUpdate(updates, n_obs = nrow(knockoffs))
  )
}


#' Output can be added to knockoffs (or their covariance) to correct for removal of a variable.
#'
getMeanUpdate = function(updates){
  with(updates,
       mean_update1_left %*% mean_update1_right +
         mean_update2_left %*% mean_update2_right
  )
}
getRandomUpdate = function(updates, n_obs){
  matrix( rnorm( n_obs ), ncol = 1 ) %*% (updates$sqrt_cov_update)
}
getCovarianceUpdate = function(updates){
  crossprod(updates$sqrt_cov_update)
}
