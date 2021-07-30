#' Gaussian mixture analog to computeGaussianKnockoffs.
#'
#' @param ... passed to computeGaussianKnockoffs. See computeGaussianKnockoffs for parameter documentation.
#' @param seed Meant for internal use only (in a unit test). This seed is set after sampling cluster assignments and before sampling knockoffs per cluster.
#' @param output_type If "knockoffs", a matrix matching X in shape.
#' If "parameters", a list with one element per cluster, containing parameters taken directly from the
#' "parameters" output of computeGaussianKnockoffs.

#' @value See output_type parameter.
#'
computeGaussianMixtureKnockoffs = function(
  X,
  mus ,
  sigmas,
  posterior_probs,
  diag_s = NULL,
  output_type = c("knockoffs", "statistics", "parameters"),
  num_realizations = 1,
  seed = NULL,
  ...
){
  stopifnot("Number of clusters must match for all inputs.\n"=length(mus)==length(sigmas))
  stopifnot("Number of clusters must match for all inputs.\n"=length(mus)==ncol(posterior_probs))
  stopifnot("Posterior probs must be n_obs by n_clusters.\n"=nrow(X)==nrow(posterior_probs))
  if( num_realizations != 1  ){
    stop("Sorry, the ability to get multiple realizations is not yet implemented.\n")
  }

  assignments = apply(posterior_probs, 1, function(p) sample(seq_along(p), prob = p, size = 1, replace = F))
  knockoffs_by_cluster = list()
  for(cluster_idx in unique(assignments)){
    if(!is.null(seed)){
      set.seed(seed)
    }
    knockoffs_by_cluster[[cluster_idx]] =
      computeGaussianKnockoffs(
        X[assignments==cluster_idx,],
        mu = mus[[cluster_idx]],
        Sigma =sigmas[[cluster_idx]],
        diag_s = diag_s[[cluster_idx]],
        ...
      )
  }
  put_back_in_order = function(data_by_cluster){
    data_in_order = X
    for(cluster_idx in unique(assignments)){
      data_in_order[assignments==cluster_idx,] = data_by_cluster[[cluster_idx]]
    }
    data_in_order
  }
  switch(
    match.arg(output_type),
    knockoffs  = knockoffs_by_cluster %>% put_back_in_order,
    parameters = knockoffs_by_cluster
  )
}

#' Gaussian mixture analog to GenerateLooks.
#'
generateGaussianMixtureLooks = function(
  X, mus, Sigmas,
  posterior_probs,
  method = c("asdp", "sdp", "equi", "group"),
  groups = NULL,
  diag_s = NULL,
  vars_to_omit = 1:ncol(X),
  statistic = knockoff::stat.glmnet_coefdiff,
  output_type = c("knockoffs", "knockoffs_compact", "statistics", "parameters", "pearson"),
  ...){
  stop("Not implemented yet.\n")
}
