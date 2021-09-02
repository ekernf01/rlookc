

test_that("Single-cluster results match the reference", {
  X =
    rbind(
      matrix( rnorm( 1e3 ), ncol = 10 )
    )
  knockoffs = computeGaussianMixtureKnockoffs(
    X,
    mus = list(0),
    sigmas = list(diag(10)),
    posterior_probs = matrix(1, ncol = 1, nrow = 100),
    seed = 0
  )
  set.seed(0)
  knockoffs_reference = computeGaussianKnockoffs(
    X,
    mu = 0,
    Sigma = diag(10)
  )
  testthat::expect_equal(knockoffs, knockoffs_reference, tolerance = 1e-8)
})

set.seed(1)
test_that("Code runs when some clusters are unused", {
  X =
    rbind(
      matrix( rnorm( 1e3 ), ncol = 10 ) + 10,
      matrix( rnorm( 1e3 ), ncol = 10 ) - 10
    )
  knockoffs = computeGaussianMixtureKnockoffs(
    X,
    mus = list(10, -10),
    sigmas = list(diag(10), diag(10)),
    posterior_probs = rep(1, 200) %o% 1:0
  )
})

set.seed(1)
test_that("Error control works", {
  X =
    rbind(
      matrix( rnorm( 1e3 ), ncol = 10 ) + 10,
      matrix( rnorm( 1e3 ), ncol = 10 ) - 10
    )
  knockoffs = computeGaussianMixtureKnockoffs(
    X,
    mus = list(10, -10),
    sigmas = list(diag(10), diag(10)),
    posterior_probs =
      data.frame(cluster = rep(c("a", "b"), each = 100)) %>%
      model.matrix(~cluster + 0, .)
  )
  bad_knockoffs = computeGaussianKnockoffs(X)
  calibration = simulateY(X, knockoffs, shuddup = T)
  calibration_bad = simulateY(X, bad_knockoffs, shuddup = T)
  expect_true(all( calibration$calibration$targeted_fdrs > calibration$calibration$fdr %>% colMeans))
  expect_true(mean( calibration_bad$calibration$targeted_fdrs < calibration_bad$calibration$fdr %>% colMeans) > 0.7)
})

