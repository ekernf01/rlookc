#' Given X, simulate Y|X and check calibration.
#'
#' @export
#' @param X @param knockoffs A real dataset and a corresponding model-X knockoff realization
#' @param statistic Function used to compute variable importance, e.g. knockoff::stat.glmnet_lambdasmax.
#' @param ... Passed to statistic, e.g. n_lambda=100.
#' @param groups Groups for composite hypothesis testing.
#' @param active_set_size how many grups to include in the active set for each simulation.
#' @param n_sim How many simulations to perform.
#' @param FUN How to simulate Y given some x's, e.g. function(x) all(x>0) + rnorm(1)
#' @param plot_savepath Passed to check.calibration
#'
simulateY = function(X,
                     knockoffs,
                     statistic = marginalScreen,
                     groups = seq(ncol(X)),
                     active_set_size = 5,
                     n_sim = 500,
                     FUN = function(x) mean(x) + 0.1*rnorm(1),
                     plot_savepath = NULL,
                     ... ){
    active_group_idx    =
      replicate(n = n_sim,
                simplify = F,
                sample(seq_along(groups), replace = F, size = n_groups)
      )
    active_variable_idx = active_group_idx %>% lapply(function(g) Reduce(union, groups[g]) )
    y = lapply(active_variable_idx, function(s) apply(X[,s, drop = F], 2, FUN))
    W =
      sapply( seq(n_sim), function( i ) statistic( X = X, X_k = knockoffs, y = y[[i]], ... ) ) %>%
      t %>%
      aggregateStats(groups)
    calibration_grouped = checkCalibration(ground_truth = active_group_idx,
                                            W = W,
                                            plot_savepath = plot_savepath)
    return(list(groups = groups, ground_truth = active_group_idx, stats = W, calibration = calibration))
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
  p_true_discoveries = fdr = matrix(NA, nrow = length(ground_truth), ncol = 10) %>% set_colnames(targeted_fdrs)
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
      selected = qvals[[i]] %>% is_weakly_less_than(targeted_fdrs[[j]]) %>% which
      fdr[i,j]                = selected %>% is.element(idx) %>% magrittr::not() %>% mean %>% na2zero
      p_true_discoveries[i,j] = selected %>% intersect(idx) %>% length %>% divide_by(length(idx))
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

