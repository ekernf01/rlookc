#' Deprecated. See calibrate__chooseDiverseY.
#'
#' @export
chooseDiverseY = function(...){ calibrate__chooseDiverseY(...) }

#' Pick P(Y|X) to span a wide range of possibilities.
#'
#' @param X  Covariates from which you want to select an active set, but you are
#' afraid your modeling assumptions might be wrong.
#' @param n_quantiles Each variable is split at this many quantiles and
#' Y is set to a collection of step functions indicating their locations.
#' Their outputs are smashed together into a matrix.
#' @return A list with values of Y, one matrix of step functions per variable.
#'
#' @export
#'
calibrate__chooseDiverseY = function(X, n_quantiles = 10 ){
  y = list()
  for(k in seq(ncol(X))){
    var_k_quantiles = quantile(X[,k], probs = (1:n_quantiles) / ( n_quantiles+1 ) )
    idx = (k-1)*n_quantiles + (1:n_quantiles)
    y[idx] = lapply( var_k_quantiles, function(q) { as.numeric( X[,k] > q ) } )
  }
  return( list( y = y, ground_truth = rep( seq( ncol( X ) ), each = n_quantiles ) ) )
}

#' Deprecated. See calibrate__simulateY
#'
#' @export
simulateY = function(...){ calibrate__simulateY(...) }

#' Given X, simulate Y|X and check calibration.
#'
#' @export
#' @param X @param knockoffs A real dataset and a corresponding model-X knockoff realization.
#' To average out over multiple knockoff realizations, pass in a list of matrices. Any number is fine
#' if your computer can handle it; the code will cycle through if it's less than n_sim.
#' @param X_observed if there is error in the covariates, you can study it in simulations
#' by passing in the true X above and the observed X (with error) here.
#' The true X will be used to simulate Y and the observed X will be used for subset selection.
#' @param statistic Function used to compute variable importance, e.g. knockoff::stat.glmnet_lambdasmax.
#' @param ... Passed to statistic, e.g. n_lambda=100, or alternate_method (these are never both used).
#' @param groups Groups for composite hypothesis testing.
#' @param active_set_size how many groups to include in the active set for each simulation.
#' @param n_sim How many simulations to perform.
#' @param FUN How to simulate Y given the active variables from a single observation, e.g. function(x) 1*all(x>0) + rnorm(1).
#' Should produce numeric output (not logical, even if you use a step function).
#' @param plot_savepath Passed to check.calibration
#' @param cores Number of cores to use.
#' @param alternate_method An alternative method of controlling FDR. Should
#' accept args y (vector length n, covariates) and X (matrix n by p). Should return
#' a vector of q-values of length p.
#'
#' @importFrom magrittr %>%
#'
#' @returns List with elements:
#' - groups: list of groups of variables, one per group unless you are testing grouped subset selection (e.g. LD blocks)
#' - ground_truth: list where each element is a list of groups of variables used in constructing the corresponding y.
#' - stats: symmetric stats from
#' - calibration: dataframe with expected and observed FDR
#' - statistic: function used inside the knockoff filter as variable or group importance
#' - FUN: function used to simulate y
#' - active_set_size: size of each active set (number of groups of variables). Provide a positive integer or a vector of
#'     positive integers of length n_sim. Each integer refers to a group of variables (in the "groups" arg above).
#' - n_sim: number of simulations
#'
calibrate__simulateY = function(X,
                                knockoffs = NULL,
                                X_observed = X,
                                statistic = knockoff::stat.glmnet_lambdasmax,
                                groups = seq(ncol(X)),
                                active_set_size = 1,
                                n_sim = 500,
                                FUN = function(x) as.numeric(all(x>0)) + rnorm(n=1,sd=0.1),
                                plot_savepath = NULL,
                                shuddup = F,
                                alternative_method = NULL,
                                rng_seed = 0,
                                cores = parallel::detectCores(),
                                ... ){
  set.seed(rng_seed)
  # You can pass in just one set of knockoffs, or a bunch of knockoff realizations.
  if(is.null(alternative_method)){
    stopifnot("Knockoffs must be a matrix or a list of matrices equal in size to X" = !is.null(knockoffs))
    if(!is.list(knockoffs)){
      knockoffs = list(knockoffs)
    }
    stopifnot("Knockoffs must be a matrix or a list of matrices equal in size to X" = dim(knockoffs[[1                ]])==dim(X))
    stopifnot("Knockoffs must be a matrix or a list of matrices equal in size to X" = dim(knockoffs[[length(knockoffs)]])==dim(X))
    if(length(knockoffs)<n_sim & !shuddup){
      warning("More simulations than knockoffs. Knockoffs will be recycled.\n")
    }
  } else {
    stopifnot("Knockoffs are ignored if alternative_method is specified." = is.null(knockoffs))
  }

  # Set up the active set
  if(length(active_set_size)==1){
    active_set_size = rep(active_set_size, n_sim)
  }
  active_group_idx = as.list(rep(NA, n_sim))
  for( i in 1:n_sim ){
    active_group_idx[[i]] = sample(
      seq_along(groups),
      replace = F,
      size = active_set_size[[i]]
    )
  }
  active_variable_idx = active_group_idx %>% lapply(function(g) Reduce(union, groups[g]) )
  # Generate the targets
  y = lapply(active_variable_idx, function(s) apply(X[,s, drop = F], 1, FUN))
  # Run the knockoff filter using different knockoffs and y's each time.
  if(is.null(alternative_method)){
    do_one = function( i ) {
      if(i%%100==1){cat("\n",i)}
      try(
        {
          w = statistic( X = X_observed,
                         # Recycle if needed. This recycling is dumb and it works best if length(y) and length(knockoffs) are
                         # relatively prime.
                         X_k = knockoffs[[1 + magrittr::mod(i-1, length(knockoffs))]],
                         y = y[[i]], ... )
          return(w)
        }
      )
      return(rep(0, ncol(X)))
    }
    W = parallel::mclapply( seq(n_sim), do_one, mc.cores = cores) %>%
      simplify2array %>%
      t %>%
      aggregateStats(groups)
    calibration = checkCalibration(ground_truth = active_group_idx,
                                   W = W,
                                   plot_savepath = plot_savepath)
  } else {
    do_one = function( i ) {
      if(i%%100==1){cat("\n",i)}
      q = tryCatch(
        expr = {alternative_method( y = y[[i]], X = X_observed, ... )},
        error = function(e) {
          warning("On iteration ", i, " alternative_method hit an error with text:\n ",
                  e, " \n. We will set q=1 and proceed.")
          rep(1, ncol(X))
        }
      )
      stopifnot("alternative_method must return a vector of qvals of length ncol(X)" = length(q)==ncol(X))
      return(q)
    }
    W = NULL
    qvals = lapply( seq(n_sim), do_one)
    calibration = checkCalibration(ground_truth = active_group_idx,
                                   qvals = qvals,
                                   plot_savepath = plot_savepath)
  }

  return(list(
    groups = groups,
    ground_truth = active_group_idx,
    stats = W,
    calibration = calibration,
    statistic = statistic,
    FUN = FUN,
    active_set_size = active_set_size,
    n_sim = n_sim
  ))
}

#' Deprecated. See calibrate__findWorstY
#'
#' @export
findWorstY = function(...){ calibrate__findWorstY(...) }


#' Given X and various simulated Y|X, find the worst-calibrated Y|X.
#'
#' @param X @param X_k A real dataset and a corresponding model-X knockoff realization.
#' @param y @param ground_truth Two lists of the same length, one containing y's that you simulated
#' and another containing the indices of the variables in X used to simulate them.
#' @param split to get an unbiased estimate of the calibration for the worst-calibrated P(Y|X),
#' half the data are used for choosing the worst and the other half for estimating its calibration.
#' @param statistic Function used to compute variable importance, e.g. knockoff::stat.glmnet_lambdasmax.
#' @param ... Passed to statistic, e.g. n_lambda=100.
#' @param plot_savepath Passed to check.calibration
#'
#' @importFrom magrittr %>%
#'
#' @export
#'
calibrate__findWorstY = function(X,
                                 X_k,
                                 y,
                                 ground_truth,
                                 split = rep(c(F, T), length.out = nrow(X)),
                                 statistic = knockoff::stat.glmnet_lambdasmax,
                                 plot_savepath = NULL,
                                 ... ){
  stopifnot("Knockoffs must be a matrix or a list of matrices equal in size to X" = dim(X_k)==dim(X))
  # Run the knockoff filter across all Y and all knockoffs.
  W = sapply(
    y,
    function(yi) {
      statistic(X[split,], X_k[split,], yi[split])
    }
  ) %>% t
  fdr = checkCalibration(ground_truth = ground_truth,
                         W = W,
                         plot_savepath = plot_savepath)$fdr %>%
    as.data.frame %>%
    dplyr::mutate(y_index = seq_along(y))

  # Identify worst-case P(Y|X)
  worst_y = fdr %>%
    dplyr::group_by(y_index) %>%
    dplyr::summarise_all(.funs = mean) %>%
    tidyr::pivot_longer(!y_index, names_to = "nominal_fdr", values_to = "empirical_fdr") %>%
    dplyr::mutate(nominal_fdr = as.numeric(nominal_fdr)) %>%
    dplyr::group_by(y_index) %>%
    dplyr::summarise(max_miscalibration = max(empirical_fdr - nominal_fdr)) %>%
    dplyr::arrange(desc(max_miscalibration)) %>%
    head(1)

  # Re-estimate miscalibration to avoid selection bias
  W = statistic(X[!split,], X_k[!split,], y[[worst_y$y_index]][!split]) %>% matrix(nrow = 1)
  calibration = checkCalibration( ground_truth = list(worst_y$y_index),
                                  W,
                                  plot_savepath = plot_savepath)
  return(list(calibration = calibration, worst_y = worst_y$y_index))
}



#' Deprecated. See calibrate__checkCalibration
#'
#' @export
checkCalibration = function(...){ calibrate__checkCalibration(...) }

#' Given a set of knockoff stats and ground-truth nonzero entries, compute
#' q-values and compare them against actual false discoveries.
#'
#' @export
#' @param ground_truth List of active sets of length n_sim
#' @param W Knockoff statistics from simulations, matrix of size n_sim by n_vars. Ignored if qvals are given.
#' @param qvals q-values from simulations. List of length n_sim containing vectors of length n_vars.
#' @param targeted_fdrs Abcissae for the resulting plots.
#' @param plot_savepath Where to save plots.
#'
calibrate__checkCalibration = function(ground_truth, W = NULL, targeted_fdrs = (1:10)/10, qvals = list(), plot_savepath = NULL, verbose = F){
  if(length(qvals)>0){
    if(!is.null(W)){
      stop("Please specify either 'W' or 'qvals' but not both.")
    }
    stopifnot(length(qvals)==length(ground_truth))
    n_var = length(qvals[[1]])
  } else {
    if(is.null(W)){
      stop("Please specify 'W' or 'qvals' (but not both).")
    }
    n_var = ncol(W)
    stopifnot(nrow(W)==length(ground_truth))
  }
  p_true_discoveries = fdr = matrix(NA, nrow = length(ground_truth), ncol = length(targeted_fdrs)) %>% magrittr::set_colnames(targeted_fdrs)
  for(i in seq_along(ground_truth)){
    if(i%%25==0 & verbose){
      cat(i, "")
    }
    idx = ground_truth[[i]]
    if(!is.null(W)){
      qvals[[i]] = knockoffQvals(W[i, ])
    }
    na2zero =function(q) { q[is.na(q)]=0; q }
    for(j in seq_along(targeted_fdrs)){
      selected = qvals[[i]] %>% magrittr::is_weakly_less_than(targeted_fdrs[[j]]) %>% which
      fdr[i,j]                = selected %>% is.element(idx) %>% magrittr::not() %>% mean %>% na2zero
      p_true_discoveries[i,j] = selected %>% intersect(idx) %>% length %>% magrittr::divide_by(length(idx))
    }
  }
  if( !is.null(plot_savepath)){
    pdf(plot_savepath, width = 5, height = 5)
  }
  errorplot = function(Y, ...){
    moe = 1.96*sqrt( apply(Y, 2, var) / nrow( Y ) )
    my_estimate = colMeans(Y)
    plot(   y = my_estimate, ...)
    points( y = my_estimate + moe, pch = "-", ...)
    points( y = my_estimate - moe, pch = "-", ...)
  }
  errorplot( x = colMeans(fdr), Y = p_true_discoveries, main = "Power", ylab = "True discoveries / active set size", xlab = "Empirical FDR")
  errorplot( x = targeted_fdrs, Y = fdr, main = "Error control", ylab = "Empirical FDR", xlab = "Targeted FDR"); abline(a = 0, b = 1)
  if( !is.null(plot_savepath)){
    dev.off()
  }
  list(targeted_fdrs = targeted_fdrs, fdr = fdr, p_true_discoveries = p_true_discoveries)
}


#' Deprecated. See calibrate__checkCalibration
#'
#' @export
knockoffQvals = function(...){ calibrate__getQvals(...) }


#' Given knockoff statistics, return q-values.
#'
#' modified from knockoff::knockoff.threshold
#' @param W Knockoff statistics, symmetric under H0
#' @param offset Control modified (offset=0) vs traditional (offset=1) FDR. Same as knockoff::knockoff.threshold.
#' @export
#'
calibrate__getQvals = function (W, offset = 1){
  if (offset != 1 && offset != 0) {
    stop("Input offset must be either 0 or 1")
  }
  q = rep(NA, length(W))
  rankByW = ecdf(W)
  rankByNegW = ecdf(-W)
  n = length(W)
  ts = sort(c(0, abs(W)))
  ratio =
    (offset + n*rankByW(-ts)) /
    pmax(1, n*rankByNegW(-ts))
  rankByThresholds = ecdf(ts)
  cum_min_ratio = cummin(ratio)
  Smax = (n+1)*rankByThresholds(W) - 2
  is_empty = Smax<1
  min_eligible_ratio = ifelse(is_empty, Inf, cum_min_ratio[pmax(Smax, 1)])
  q = min_eligible_ratio %>% pmin(1)
  q
}

#' Compute correlations, but constant input gives 0 instead of causing a panic
#'
safeCorrelation = function(...) {r = cor(...); r[is.na(r)] = 0; return(r)}

#' Deprecated. See stat.marginalScreen.
#'
#' @export
marginalScreen = function(...){ stat.marginalScreen(...) }

#' A simple, fast test statistic, symmetric under the null, useful for testing knockoff software
#'
#' @export
stat.marginalScreen =   function(X, X_k, y) as.vector(safeCorrelation(X, y) - safeCorrelation(X_k, y))

#' A knockoff statistic capable of handling multivariate Y.
#'
#' For a single variable or knockoff x, the scheme is to regress x (as the target)
#' on Y (as the covariates) and compute cor(x, hat x) where hat x are the fitted values.
stat.CCA = function(X, X_k, y){
  hat_matrix = solve(t(y) %*% y) %*% t(y)
  do_one = function(x){
    cor( y %*% ( hat_matrix %*% x ), x )
  }
  apply(X, 2, do_one) - apply(X_k, 2, do_one)
}

#' Deprecated. See calibrate__KNNTest
#'
#' @export
KNNTest = function(...){ calibrate__KNNTest(...) }

#' Run the "KNN diagnostics" from Section 5.1 (page 14) of the Romano et al. Deep Knockoffs paper.
#'
#' @param X Features
#' @param X_k knockoffs
#' @param swap_type Full compares X, X_k to X_k, X. Partial swaps each column with probability 0.5.
#'
#' @export
#'
calibrate__KNNTest = function(X, X_k, n_neighbors = min(20, nrow(X)), swap_type = c("full", "partial")){
  if (!requireNamespace("FNN", quietly=TRUE)) {
    stop("Running the KNN diagnostic requires the FNN package.\n")
  }
  stopifnot(dim(X)==dim(X_k))
  # Do swaps
  p = ncol(X)
  Z = cbind(X, X_k)
  Zswap = Z
  zero_p = c(0, p)
  swap_type = match.arg(swap_type)
  if(swap_type == "full"){
    vars_to_swap = 1:p
  } else if(swap_type == "partial"){
    vars_to_swap = which(1==rbinom(prob = 0.5, size = 1, n = p))
  } else {
    stop("Swap type not recognized\n.")
  }
  for(j in vars_to_swap){
    Zswap[,j + zero_p] = Z[,j + rev(zero_p)]
  }
  # get neighbors
  neighbors = FNN::get.knn( data = rbind(Z, Zswap),
                            k = n_neighbors + 1)$nn.index
  proportion_not_swapped = c()
  for(i in seq(2*nrow(X))){
    current_side = ifelse(i>nrow(X), "Z_swap", "Z")
    is_neighbor_on_current_side = function(nn){
      if(current_side=="Z_swap"){
        return(nn > nrow(X))
      } else {
        return(nn <= nrow(X))
      }
    }
    proportion_not_swapped[[i]] =
      neighbors[i,] %>%
      setdiff(i + c(0, nrow(X), -nrow(X))) %>% # same row of Z or Zswap is not eligible
      magrittr::extract(1:n_neighbors) %>% #If the last line did nothing, we still want just 20
      is_neighbor_on_current_side() %>%
      mean
  }
  proportion_not_swapped = unlist(proportion_not_swapped)
  # Null distribution taken from:
  #
  #  Schilling, M. F. (1986). Multivariate two-sample tests based on nearest neighbors.
  # Journal of the American Statistical Association, 81(395), 799-806.
  #
  # Theorem 3.1 with lambdas=1/2 and using the limiting value from 3.2.
  null_mean = 0.5
  null_variance = 0.5 / (2*nrow(X)*n_neighbors)
  list(
    prop_not_swapped_per_observation = proportion_not_swapped,
    prop_not_swapped = mean(proportion_not_swapped),
    p_value = pnorm(q = mean(proportion_not_swapped), mean = null_mean, sd = sqrt(null_variance), lower.tail = F)
  )
}


