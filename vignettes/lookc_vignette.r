# Generate data from a simple chain-like hierarchical model
set.seed(0)
p = 500
X = matrix(NA, ncol = p, nrow = 1e4)
X[,1] = rnorm(1e4)
for(k in 2:p){
  X[,k] = rnorm(1e4)
  if(k-100 < 4){
    X[,k] = X[,k] + X[,k-1]
  }
}
Sigma = cov(X)
# Leave one out using Dr. Sesia's original code
filter_100 = knockoff::knockoff.filter(X[,-100], y = X[,100], statistic = knockoff::stat.glmnet_coefdiff, fdr = 0.5)
# Leave one out using a wrapper to Dr. Sesia's original code
stats_from_wrapper      = rlookc::generateLooksSlow(X, mu = 0, Sigma = Sigma, vars_to_omit = 100, output_type = "statistics")
# Leave one out using a reimplementation that I hope will eventually be faster for our use case
stats_from_reimplementation = rlookc::generateLooks(X, mu = 0, Sigma = Sigma, vars_to_omit = 100, output_type = "statistics")

# All three show good error control on this toy example
filter_100
plot(stats_from_wrapper[[1]])
plot(stats_from_reimplementation[[1]])

# Regarding scaling, here's 10 looks by brute force and by low-rank updating.
microbenchmark::microbenchmark(
  times = 1,
  looks_from_reimplementation = rlookc::generateLooks(X, mu = 0, Sigma = Sigma, vars_to_omit = 96:100, output_type = "knockoffs")
)
microbenchmark::microbenchmark(
  times = 1,
  looks_brute_force =  lapply(96:100, function(k) {knockoff::create.gaussian(X[,-k], mu = 0, Sigma = Sigma)})
)
