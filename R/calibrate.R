#' Pick P(Y|X) to span a wide range of possibilities.
#'
#' @param X  Covariates from which you want to select an active set, but you are
#' afraid your modeling assumptions might be wrong.
#' @param n_quantiles Each variable is split at this many quantiles and
#' Y is set to a collection of step functions indicating their locations.
#' Their outputs are smashed together into a matrix.
#' @return A list with values of Y, one matrix of step functions per variable.
#'
chooseDiverseY = function(X, n_quantiles = 10 ){
  y = list()
  for(k in seq(ncol(X))){
    var_k_quantiles = quantile(X[,k], probs = (1:n_quantiles) / ( n_quantiles+1 ) )
    y[[k]] = sapply( var_k_quantiles, function(q) { X[,k] > q } )
  }
  return(y)
}

#' Pick P(Y|X) purposefully to reveal bad calibration.
#'
#' @param X Covariates from which you want to select an active set, but you are
#' afraid your modeling assumptions might be wrong.
#' @param kmeans_centers The diagnostic looks for situations where the variance of
#' a given variable differs sharply across clusters. You have to manually
#' specify how many clusters
#'
#' @return A list with metadata and values of Y.
#' Each Y is generated from a univariate step function.
#' metadata includes the true active variable,
#' the cutoff used in the step function,
#' and the variable likely to be found as a false positive
#' if homoskedastic knockoffs are used.
#'
chooseAdversarialY = function(X, kmeans_centers ){
  # Cluster the dataset and look for heteroskedasticity
  clusters = stats::kmeans(X, centers = kmeans_centers)$cluster
  variances = matrix(NA, nrow = kmeans_centers, ncol = ncol(X))
  for(j in unique(clusters)){
    variances[j,] = apply(X[j==clusters,], 2, var)
  }
  # Cause each variable in turn to appear as a false positive
  # Strategy: find an indicator function for a low-variance
  # region, but based on a different variable.
  output = data.frame(active_variable = NA,
                      likely_false_hit = seq(ncol(X)),
                      cutoff = NA)
  y =list()
  for(k in seq(ncol(X))){
    low_variance_indicator = which.min(variances[,k])==clusters
    output[k, "active_variable"] = which.max(abs(cor(X, low_variance_indicator)))
    # Set the cutoff above or below to separate high from low variance
    direction_of_association = sign(cor(low_variance_indicator, X[, output[k, "active_variable"]]))
    output[k, "cutoff"] =
      X[low_variance_indicator, output[k, "active_variable"]] %>%
      quantile(0.5 - direction_of_association*0.4)
    y[[k]] = as.numeric(X[,output[k, "active_variable"]] > output[k, "cutoff"])
  }
  return(list(metadata = output, y = y))
}

#' Given X, simulate Y|X and check calibration.
#'
#' @export
#' @param X @param knockoffs A real dataset and a corresponding model-X knockoff realization.
#' To average out over multiple knockoff realizations, pass in a list of matrices. Any number is fine
#' if your computer can handle it; the code will cycle through if it's less than n_sim.
#' @param statistic Function used to compute variable importance, e.g. knockoff::stat.glmnet_lambdasmax.
#' @param ... Passed to statistic, e.g. n_lambda=100.
#' @param groups Groups for composite hypothesis testing.
#' @param active_set_size how many grups to include in the active set for each simulation.
#' @param n_sim How many simulations to perform.
#' @param FUN How to simulate Y given some x's, e.g. function(x) all(x>0) + rnorm(1).
#' Should produce numeric output (not logical, even if you use a step function).
#' Pass the string "adversarial" or "diverse" for custom options designed to reveal
#' hidden harms of heteroskedasticity and other unexpected pitfalls.
#' @param kmeans_centers passed to chooseAdversarialY
#' @param plot_savepath Passed to check.calibration
#'
#' @importFrom magrittr %>%
#'
simulateY = function(X,
                     knockoffs,
                     statistic = knockoff::stat.glmnet_lambdasmax,
                     groups = seq(ncol(X)),
                     active_set_size = 1,
                     n_sim = 500,
                     FUN = function(x) as.numeric(x>0),
                     plot_savepath = NULL,
                     shuddup = F,
                     kmeans_centers = kmeans_centers,
                     ... ){
  # You can pass in just one set, or a bunch of knockoff realizations.
  if(!is.list(knockoffs)){
    knockoffs = list(knockoffs)
  }
  stopifnot("Knockoffs must be a matrix or a list of matrices equal in size to X" = dim(knockoffs[[1                ]])==dim(X))
  stopifnot("Knockoffs must be a matrix or a list of matrices equal in size to X" = dim(knockoffs[[length(knockoffs)]])==dim(X))
  if(length(knockoffs)<n_sim & !shuddup){
    warning("More simulations than knockoffs. Knockoffs will be recycled.\n")
  }
  # You can specify Y|X and choose X's randomly, or have it done
  # in a way designed to increase sensitivity.
  {
    if( identical( FUN, "adversarial" ) )
    {
      if( !missing( active_set_size ) ){
        warning("Adversarial Y|X is only implemented for active_set_size=1 for now.\n")
        active_set_size = 1
      }
      if( !missing( groups ) ){
        warning("Adversarial Y|X is only implemented for non-grouped tests for now.\n")
        # groups = seq(length.out = n_sim)
      }
      adversarial = chooseAdversarialY(X, kmeans_centers = kmeans_centers)
      y = adversarial$y
      active_group_idx = active_variable_idx = adversarial$metadata[,"active_variable"]
      if(n_sim > length(y)){
        warning("Adversarial y diagnostic code only generates one y per input variable. Some will be recycled.\n")
        active_group_idx %<>% rep(length.out = n_sim) %>% as.list()
      }
    }
    else if( identical( FUN, "diverse" ) )
    {
      if( !missing( active_set_size ) ){
        warning("The 'diverse' scheme for Y|X is only implemented for active_set_size=1 for now.\n")
      }
      if( !missing( groups ) ){
        warning("The 'diverse' scheme for Y|X is only implemented for non-grouped tests for now.\n")
      }
      if( !missing( statistic ) ){
        warning("The 'diverse' scheme for Y|X has multivariate 'y' so 'statistic' will be overridden by a special default.\n")
      }
      active_set_size = 1
      # groups = seq(length.out = n_sim)
      statistic = stat.CCA
      y = chooseDiverseY(X)
      active_group_idx = active_variable_idx = seq_along(y)
      if(n_sim != length(y)){
        warning("The 'diverse' scheme for Y|X generates one y per input variable. If n_sim is different, some Y's will be recycled or ignored.\n")
        active_group_idx %<>% rep(length.out = n_sim) %>% as.list()
      }
    }
    else
    {
      active_group_idx    =
        replicate(n = n_sim,
                  simplify = F,
                  sample(seq_along(groups), replace = F, size = active_set_size)
        )
      active_variable_idx = active_group_idx %>% lapply(function(g) Reduce(union, groups[g]) )
      y = lapply(active_variable_idx, function(s) apply(X[,s, drop = F], 1, FUN))
    }
  }
  # Run the knockoff filter using different knockoffs and y's each time.
  # Recycle if needed. This recycling works best if length(y) and length(knockoffs) are
  # relatively prime... sorryyyyyy...
  do_one = function( i ) {
    statistic( X = X,
               X_k = knockoffs[[1 + magrittr::mod(i-1, length(knockoffs))]],
               y = y[[1 + magrittr::mod(i-1, length(y))]], ... )
  }
  W =
    sapply( seq(n_sim), do_one) %>%
    t %>%
    aggregateStats(groups)
  # How did it do?
  calibration_grouped = checkCalibration(ground_truth = active_group_idx,
                                         W = W,
                                         plot_savepath = plot_savepath)
  return(list(
    groups = groups,
    ground_truth = active_group_idx,
    stats = W,
    calibration = calibration_grouped,
    statistic = statistic,
    FUN = FUN,
    active_set_size = active_set_size,
    n_sim = n_sim
  ))
}

#' Given a set of knockoff stats and ground-truth nonzero entries, compute
#' q-values and compare them against actual false discoveries.
#'
#' @export
#' @param ground_truth List of active sets of length n_sim
#' @param W Knockoff statistics from simulations, matrix of size n_sim by n_vars
#' @param targeted_fdrs Abcissae for the resulting plots.
#' @param plot_savepath Where to save plots.
#'
checkCalibration = function(ground_truth, W, targeted_fdrs = (1:10)/10, plot_savepath = NULL, verbose = F, n_var = ncol(W)){
  qvals = list()
  p_true_discoveries = fdr = matrix(NA, nrow = length(ground_truth), ncol = 10) %>% magrittr::set_colnames(targeted_fdrs)
  stopifnot(nrow(W)==length(ground_truth))
  stopifnot(ncol(W)==n_var)
  for(i in seq_along(ground_truth)){
    if(i%%25==0 & verbose){
      cat(i, "")
    }
    idx = ground_truth[[i]]
    qvals[[i]] = knockoffQvals(W[i, ])
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

#' Given knockoff statistics, return q-values.
#'
#' modified from knockoff::knockoff.threshold
#' @param W Knockoff statistics, symmetric under H0
#' @export
#'
knockoffQvals = function (W, offset = 1) {
  q = rep(NA, length(W))
  if (offset != 1 && offset != 0) {
    stop("Input offset must be either 0 or 1")
  }
  ts = sort(c(0, abs(W)))
  ratio = sapply(ts, function(t) (offset + sum(W <= -t))/max(1, sum(W >= t)))
  # This step is quadratic in length(W) and could probably be done faster if needed
  for(k in seq_along(W)){
    q[[k]] = min(c(ratio[W[[k]]>ts], 1))
  }
  q
}

#' Compute correlations, but constant input gives 0 instead of causing a panic
#'
safeCorrelation = function(...) {r = cor(...); r[is.na(r)] = 0; return(r)}

#' A simple, fast test statistic, symmetric under the null, useful for testing knockoff software
#'
marginalScreen =   function(X, X_k, y) as.vector(safeCorrelation(X, y) - safeCorrelation(X_k, y))

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
