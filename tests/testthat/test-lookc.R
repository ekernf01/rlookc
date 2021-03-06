# Generate data from a simple chain-like hierarchical model
p = 300
X = matrix(NA, ncol = p, nrow = 1e4)
X[,1] = rnorm(1e4)
for(k in 2:p){
  X[,k] = rnorm(1e4)
  X[,k] = X[,k] + X[,k-1]
}
X = X %*% diag(1/apply(X, 2, sd))
max(cov(X))

test_that("mean and covariance match", {
  looks_from_wrapper     = rlookc::generateLooksSlow(X, mu = 0, Sigma = cov(X), vars_to_omit = 100, output_type = "parameters")
  looks_from_reimplementation = rlookc::generateLooks(X, mu = 0, Sigma = cov(X), vars_to_omit = 100, output_type = "parameters")
  testthat::expect_equal( looks_from_wrapper[[1]][["ko_mean"]],
                          looks_from_reimplementation[[1]][["ko_mean"]], tolerance = 1E-2 )
  testthat::expect_equal( crossprod(looks_from_wrapper[[1]][["ko_sqrt_covariance"]]),
                          crossprod(looks_from_reimplementation[[1]][["ko_sqrt_covariance"]]), tolerance = 1E-2 )
})

