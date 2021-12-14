# This vignette demonstrates the new features of our software that
# allow construction of knockoffs:
#
# - Gaussian knockoffs with p>n
# - Gaussian mixture model knockoffs
# - fast leave-one-out knockoffs (LOOKs)
# - grouped Gaussian model-X knockoffs

library("magrittr")
set.seed(0)

# Generate leave-one-out knockoffs fast
looks = rlookc::create.looks(matrix(rnorm(1e4), nrow = 1e2), mu = 0, Sigma = diag(100), output_type = "knockoffs_compact")
names(looks)
dim(looks$knockoffs)
lapply(looks$updates[[1]], dim)
# Save them to disk in a compact form
rlookc::saveload.saveCompactLooks(looks, "demo_compact_looks")
# Load them from disk
looks = rlookc::saveload.loadCompactLooks("demo_compact_looks")
# Clean up the files so the vignette has no side effects
unlink("demo_compact_looks", recursive = T)
# Reconstitute the LOOKs in a more immediately usable form
# To save memory, you can do this for one variable at a time.
looks_one_at_a_time = rlookc::saveload.formOneLook(looks$knockoffs,
                                                   vars_to_omit = looks$vars_to_omit,
                                                   updates = looks$updates,
                                                   k = 17)
dim(looks_one_at_a_time)
# Or just do them all at once, if you have the RAM for it.
looks_explicit = rlookc::saveload.formAllLooks(looks$knockoffs,
                                               vars_to_omit = looks$vars_to_omit,
                                               updates = looks$updates)
length(looks_explicit)
sapply(looks_explicit, dim)
# We also have code to load LOOKs into Python and Julia and reconstitute them there.
# It is minimal but it is unit-tested.

# You can also output explicit knockoffs or feature-importance statistics immediately upon creating LOOKs.
# For the latter, it assumes you're doing structure learning, so Y is set equal to the omitted feature.
looks_explicit  = rlookc::create.looks(matrix(rnorm(1e4), nrow = 1e2), mu = 0, Sigma = diag(100), output_type = "knockoffs")
symmetric_stats = rlookc::create.looks(matrix(rnorm(1e4), nrow = 1e2), mu = 0, Sigma = diag(100), output_type = "statistics")
sapply(symmetric_stats, length)

# We also have a couple of other features for scalable knockoff construction.
# Generate high-dimensional Gaussian knockoffs
rather_wide_knockoffs = rlookc::createHighDimensionalKnockoffs( matrix(rnorm(1e7), nrow = 1e2) )
dim(rather_wide_knockoffs)
# Generate mixture-model knockoffs
X = matrix(rnorm(1e4), nrow = 1e3)
X = rbind(X, matrix(rnorm(1e4), nrow = 1e3) + 10)
library("mclust")
fitted_mixture_model = mclust::Mclust(X)
mixture_model_knockoffs = rlookc::create.gaussianMixtureKnockoffs(
  X = X,
  do_high_dimensional = F,
  mus = fitted_mixture_model$parameters$mean %>%
    as.data.frame() %>%
    as.list(),
  sigmas = list(
    fitted_mixture_model$parameters$variance$Sigma,
    fitted_mixture_model$parameters$variance$Sigma
  ),
  posterior_probs = fitted_mixture_model$z
)
plot(X[,1:2], col = "black", main = "Quick exchangeability demo")
points(mixture_model_knockoffs[,1:2], col = "red", pch = ".")
points(X[,1], mixture_model_knockoffs[,2], col = "blue", pch = "+")
points(mixture_model_knockoffs[,1], X[,2], col = "green", pch = "*")

# Generate grouped knockoffs (experimental)
# Start with nearly-duplicated features
my_groups = lapply(1:10, function(k) ((k-1)*10) + 1:10)
my_duplicated_features = matrix(0, nrow = 1e3, ncol = 1e2)
most_of_the_info = 0
for(i in 1:10){
  most_of_the_info = rnorm(1000) + 0.5*most_of_the_info
  for(j in 1:10){
    my_duplicated_features[,my_groups[[i]][[j]]] = most_of_the_info + 0.001*rnorm(1000)
  }
}
# Use the same groups for knockoff construction
grouped_knockoffs = rlookc::create.gaussianKnockoffs(
  my_duplicated_features,
  groups = my_groups,
  method = "group"
)
# Make non-grouped knockoffs for comparison
plain_knockoffs =  knockoff::create.second_order(
  my_duplicated_features
)
# Note the diagonal of the off-diagonal block.
# Grouped knockoffs are much less correlated with the original features.
# Power will be correspondingly higher.
image(cor(cbind(my_duplicated_features, grouped_knockoffs)),
      main = "Joint correlation with grouped knockoffs")
image(cor(cbind(my_duplicated_features, plain_knockoffs)),
      main = "Joint correlation with plain knockoffs")


