# Sim data
test_that("high-dimensional knockoffs are mathematically correct", {
  set.seed(0)
  X = matrix(rnorm(300), nrow = 10)
  mu = apply(X, 2, mean)
  X = sweep(X, 2, mu, FUN = "-")
  standard_deviations = apply(X, 2, sd)
  X = sweep(X, 2, standard_deviations, FUN = "/")
  lambda = 0.1
  rho = 0.9
  params = createHighDimensionalKnockoffs(X, output_type = "parameters", lambda = lambda)
  Sigma = cor(X)*(1-lambda) + lambda*diag(ncol(X))
  I = diag(ncol(X))
  set_dimnames = function(X, nm) {
    dimnames(X) = nm
    X
  }
  # Check Sigma and its inverse
  I %>% (params$multiplyBySigma) %>% as.matrix %>% set_dimnames(NULL) %>% expect_equal(Sigma)
  params$multiplyBySigma(I, additional_transform = function(x) 1/x) %>% as.matrix %>% set_dimnames(NULL) %>% expect_equal(solve(Sigma))

  # Check mean
  S = 2*rho*lambda*I
  knockoff_mean_other = X - X %*% solve(Sigma) %*% S
  # plot(knockoff_mean_other, params$knockoff_mean)
  # abline(a = 0, b = 1, col = "red")
  expect_equal(knockoff_mean_other, params$knockoff_mean %>% as.matrix %>% set_dimnames(NULL))

  # Check covariance
  actual_covariance = 2*S - ( S %*% solve(Sigma) %*% S )
  what_i_think_it_is = params$multiplyBySigma(I, additional_transform = params$transform_eigenvectors_from_sigma_to_sqrt_c)
  what_i_think_it_is = what_i_think_it_is %*% what_i_think_it_is
  # plot(actual_covariance,  what_i_think_it_is      )
  # abline(a = 0, b = 1, col = "red")
  what_i_think_it_is %>% as.matrix %>% set_dimnames(NULL) %>% expect_equal(actual_covariance)
})
