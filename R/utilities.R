#' Given a set of knockoff stats and ground-truth nonzero entries, compute
#' q-values and compare them against actual false discoveries.
#'
#' @export
#'
check.calibration = function(ground_truth, W, targeted_fdrs = (1:10)/10, plot_savepath = NULL, verbose = F, n_var = ncol(W)){
  qvals = list()
  p_true_discoveries = fdr = matrix(NA, nrow = length(ground_truth), ncol = 10) %>% set_colnames(targeted_fdrs)
  stopifnot(nrow(W)==length(ground_truth))
  stopifnot(ncol(W)==n_var)
  for(i in seq_along(ground_truth)){
    if(i%%25==0 & verbose){
      cat(i, "")
    }
    idx = ground_truth[[i]]
    qvals[[i]] = knockoff.qvals(W[i, ])
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
knockoff.qvals = function (W, offset = 1) {
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
safecor = function(...) {r = cor(...); r[is.na(r)] = 0; return(r)}

#' A simple, fast test statistic, symmetric under the null, useful for testing knockoff software
#'
marginal_screen =   function(X, X_k, y) as.vector(safecor(X, y) - safecor(X_k, y))

