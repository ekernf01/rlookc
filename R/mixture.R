#' Gaussian mixture analog to computeGaussianKnockoffs.
#'
#' @param ... passed to computeGaussianKnockoffs or highDimensionalKnockoffs.
#' @param seed Meant for internal use only (in a unit test). This seed is set after sampling cluster assignments and before sampling knockoffs per cluster.
#' @param do_high_dimensional If T, use highDimensionalKnockoffs instead of computeGaussianKnockoffs.
#' @param posterior_probs @param hard_assignments Soft or hard cluster assignments. Specify one but not both.
#' @param mus  @param sigmas  @param diag_s Lists of means, covariances, and valid choices of S, with one element per cluster.
#' @param lambdas List of shrinkage parameters with one element per cluster. Set with corpcor::estimate.lambda for optimal MSE.
#' Used only in the high-dimensional case.
#' @param output_type If "knockoffs", a matrix matching X in shape.
#' If "parameters", a list with one element per cluster, containing parameters taken directly from the
#' "parameters" output of computeGaussianKnockoffs or highDimensionalKnockoffs.

#' @return See output_type parameter.
#' @export
#'
computeGaussianMixtureKnockoffs = function(
  X,
  do_high_dimensional = F,
  mus = NULL,
  sigmas = NULL,
  posterior_probs = NULL,
  hard_assignments = NULL,
  diag_s = NULL,
  lambdas = NULL,
  output_type = c("knockoffs", "statistics", "parameters"),
  num_realizations = 1,
  seed = NULL,
  ...
){
  # Check cluster membership input
  stopifnot("Input only posterior probs or hard assignments.\n"=(!is.null(posterior_probs) & !is.null(posterior_probs)))
  if(!is.null(posterior_probs)){
    stopifnot("Posterior probs must be n_obs by n_clusters.\n"=nrow(X)==nrow(posterior_probs))
    n_cluster = ncol(posterior_probs)
  }
  if(!is.null(hard_assignments)){
    stopifnot("hard_assignments must be n_obs in length \n"=nrow(X)==length(hard_assignments))
    stopifnot("hard_assignments must contain only numbers 1:K"=is.integer(hard_assignments))
    stopifnot("hard_assignments must contain only numbers 1:K"=all(hard_assignments>0))
    n_cluster = max(hard_assignments)
  }
  # Check other input
  if(do_high_dimensional){
    stopifnot("Cannot handle user-specified parameters in the high-dimensional case."=is.null(sigmas))
    stopifnot("Cannot handle user-specified parameters in the high-dimensional case."=is.null(mus))
    stopifnot("Cannot handle user-specified parameters in the high-dimensional case."=is.null(diag_s))
  } else {
    stopifnot()
  }
  if(!is.null(sigmas)){
    stopifnot("Number of clusters must match for all inputs.\n"=length(sigmas)==n_cluster)
  }
  if(!is.null(mus)){
    stopifnot("Number of clusters must match for all inputs.\n"=length(mus)==n_cluster)
  }
  if( num_realizations != 1  ){
    stop("Sorry, the ability to get multiple realizations is not yet implemented.\n")
  }
  # Assign obs for knockoffs
  if(!is.null(hard_assignments)){
    assignments = hard_assignments
  } else {
    assignments = apply(posterior_probs, 1, function(p) sample(seq_along(p), prob = p, size = 1, replace = F))
  }
  # Generate knockoffs by cluster
  knockoffs_by_cluster = list()
  for(cluster_idx in unique(assignments)){
    if(!is.null(seed)){
      set.seed(seed)
    }
    if(do_high_dimensional){
      knockoffs_by_cluster[[cluster_idx]] =
        highDimensionalKnockoffs(
          X[assignments==cluster_idx,,drop = F],
          lambda = lambdas[[cluster_idx]]
        )
    } else {
      knockoffs_by_cluster[[cluster_idx]] =
        computeGaussianKnockoffs(
          X[assignments==cluster_idx,,drop = F],
          mu = mus[[cluster_idx]],
          Sigma =sigmas[[cluster_idx]],
          diag_s = diag_s[[cluster_idx]],
          ...
        )
    }
  }
  # Return'em in the original order
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
#' @export
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
