# This vignette shows off the scaling of our leave-one-out construction.

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
mu = rep(0, ncol(X))
# Leave one out using Dr. Sesia's original code
filter_100 = knockoff::knockoff.filter(X[,-100], y = X[,100], statistic = knockoff::stat.glmnet_coefdiff, fdr = 0.5)
# Leave one out using a wrapper to a minor modification of Dr. Sesia's original code
stats_from_wrapper          = rlookc::create__looksSlow(X, mu = mu, Sigma = Sigma, vars_to_omit = 100, output_type = "statistics")
# Leave one out using lookc
stats_from_reimplementation = rlookc::create__looks(    X, mu = mu, Sigma = Sigma, vars_to_omit = 100, output_type = "statistics")

# All three show a very obvious correct signal on this toy example
filter_100
plot(stats_from_wrapper[[1]])
plot(stats_from_reimplementation[[1]])

# Regarding scaling, here are 5 looks done by brute force and then by low-rank updating.
# On my 2016 Dell XPS13 with Dual-Core Intel® Core™ i5-6200U CPU @ 2.30GHz,
# the low-rank updates win by a factor of 4, 46s to 11s -- nearly ideal scaling.
five_looks_sesia = microbenchmark::microbenchmark(
  times = 1,
  looks_brute_force =  lapply(96:100, function(k) {knockoff::create.gaussian(X[,-k], mu = mu[-k], Sigma = Sigma[-k,-k])})
)
five_looks_lookc = microbenchmark::microbenchmark(
  times = 1,
  looks_from_reimplementation = rlookc::create__looks(X, mu = mu, Sigma = Sigma, vars_to_omit = 96:100, output_type = "knockoffs")
)

# In fact, we can we generate all 100 of the LOOKs and compute statistics in a couple of minutes.
# Runs in 105s on my computer.
do_pearson_screen = function(X, ko, y){
  rbind(cor(X, y), cor(ko, y)) %>% set_rownames(c("X", "knockoff"))
}
five_hundred_looks_lookc = microbenchmark::microbenchmark(
  times = 1,
  looks_from_reimplementation = rlookc::create__looks(X,
                                                      mu = mu,
                                                      Sigma = Sigma,
                                                      vars_to_omit = 1:500,
                                                      statistic = do_pearson_screen,
                                                      output_type = "statistic")
)
knitr::kable(t(t(c(five_looks_sesia = five_looks_sesia$time, five_looks_lookc = five_looks_lookc$time, five_hundred_looks_lookc = five_hundred_looks_lookc$time ))) / 1e9)

