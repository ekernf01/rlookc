
# ========= ========= ========= ========= =========
#                      Purpose
# ========= ========= ========= ========= =========
#
# This vignette demonstrates the new features of our software that
# allow construction of knockoffs:
#
# - fast and low-memory leave-one-out knockoffs (LOOKs)
# - fast and low-memory Gaussian knockoffs with p>n
# - Gaussian mixture model knockoffs
# - grouped Gaussian model-X knockoffs
#

# ========= ========= ========= ========= =========
#                    Organization
# ========= ========= ========= ========= =========
#
# Most functions fall into one of three categories, each marked by a prefix.
#
# - "create__" functions generate knockoffs. This vignette covers these functions.
# - "saveLoad__" functions save knockoffs to disk, or load them, or convert them from 
#   a compact form into an immediately useful form. This vignette covers some of these functions.
# - "calibrate__" functions help check if the knockoffs are good enough to control false discoveries.
#   More accurately, they will tell you if the knockoffs are obviously/detectably bad. There is a
#   separate vignette about these functions.
#

# ========= ========= ========= ========= =========
#                       Demos
# ========= ========= ========= ========= =========
library("magrittr")
set.seed(0)

# The namesake of this package is the ability to generate leave-one-out knockoffs fast.
# The flagship function create__looks gives a list of matrices where element i contains knockoffs for all but the ith feature, i.e. for X[,-i].
looks_explicit  = rlookc::create__looks(matrix(rnorm(1e4), nrow = 1e2), mu = 0, Sigma = diag(100), output_type = "knockoffs")
dim(looks_explicit[[6]])

# You can also generate leave-one-out knockoffs with a low memory footprint.
# These are stored as a set of knockoffs generated with no variables left out, plus
# a set of low-rank updates to correct them as if they had been done with certain variables left out.
# The eventual result is mathematically equivalent to the above, but the internals are more efficient.
looks = rlookc::create__looks(matrix(rnorm(1e4), nrow = 1e2), mu = 0, Sigma = diag(100), output_type = "knockoffs_compact")
# You can explore the structure of the compact output if you wish.
names(looks)
dim(looks$knockoffs) # This shows the full-data knockoffs, with no variables left out.
lapply(looks$updates[[1]], dim) # This shows the size of the low-rank updates.

# You can save the LOOKs to disk in this compact form and load them later.
# Our Python and Julia code can also read these saved files.
rlookc::saveLoad__saveCompactLooks(looks, "demo_compact_looks")
looks = rlookc::saveLoad__loadCompactLooks("demo_compact_looks")
# Clean up the files so the vignette has no lasting side effects.
unlink("demo_compact_looks", recursive = T)

# For the compact knockoffs, they are not immediately usable. 
# We provide functions that will reconstitute the LOOKs into an immediately usable form, like putting a raisin in water. 
# To save memory, you can do this for one variable at a time. (Not the whole box of raisins.)
# For example, this will give you Gaussian knockoffs for X[,-17].
looks_just_one = rlookc::saveLoad__formOneLook(looks$knockoffs,
                                              vars_to_omit = looks$vars_to_omit,
                                              updates = looks$updates,
                                              k = 17)
dim(looks_just_one)
# Or you can do the whole box of raisins all at once, if you have the RAM for it.
looks_explicit = rlookc::saveLoad__formAllLooks(looks$knockoffs,
                                               vars_to_omit = looks$vars_to_omit,
                                               updates = looks$updates)
length(looks_explicit)
sapply(looks_explicit, dim)

# Another way to save memory is to output feature-importance statistics immediately upon creating LOOKs.
# This will make one set of knockoffs at a time instead of making them all at once, saving memory.
# This setting assumes you're doing structure learning, so no "Y" is needed.
# The argument vars_to_omit indicates which variables to leave out.
symmetric_stats = rlookc::create__looks(matrix(rnorm(1e4), nrow = 1e2), mu = 0,
                                        Sigma = diag(100),
                                        vars_to_omit = 1:10,
                                        output_type = "statistics")
sapply(symmetric_stats, length)

# We also have a couple of other features for scalable knockoff construction.
# You can generate high-dimensional Gaussian knockoffs (p >> n) very efficiently.
rather_wide_knockoffs = rlookc::create__highDimensionalKnockoffs( matrix(rnorm(1e7), nrow = 1e2) )
dim(rather_wide_knockoffs)
# You can generate mixture-model knockoffs from a Gaussian mixture model fitted using mclust.
X = matrix(rnorm(1e4), nrow = 1e3)
X = rbind(X, matrix(rnorm(1e4), nrow = 1e3) + 10)
library("mclust")
fitted_mixture_model = mclust::Mclust(X)
mixture_model_knockoffs = rlookc::create__gaussianMixtureKnockoffs(
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

# You can even generate grouped knockoffs (but this capability is experimental).
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
grouped_knockoffs = rlookc::create__gaussianKnockoffs(
  my_duplicated_features,
  groups = my_groups,
  method = "group"
)
# Make non-grouped knockoffs for comparison
plain_knockoffs =  knockoff::create.second_order(
  my_duplicated_features
)
# Grouped knockoffs are much less correlated with the original features.
# Power will be correspondingly higher.
# You can see it in the joint correlation matrix. Check the values that are off the diagonal but parallel to it.
image(cor(cbind(my_duplicated_features, grouped_knockoffs)),
      main = "Joint correlation with grouped knockoffs")
image(cor(cbind(my_duplicated_features, plain_knockoffs)),
      main = "Joint correlation with plain knockoffs")


