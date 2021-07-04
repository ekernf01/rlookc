library(ggplot2)
test_that("q-values match reference", {
  set.seed(0)
  x = rnorm(1000) + 1
  for(threshold in (1:9)/10){
    expect_equal(
      knockoffQvals(x) < threshold,
      x > knockoff::knockoff.threshold(x, fdr = threshold)
    )
  }
})


test_that("Simulations run", {
  expect_silent({
    calibration = simulateY(X = matrix(rnorm(10000), ncol = 50),
                            knockoffs = replicate(50, matrix(rnorm(10000), ncol = 50), simplify = F),
                            n_sim = 50,
                            shuddup = T)
  })
  X = matrix(rnorm(10000), ncol = 50)
  expect_silent({
    calibration = findWorstY(
      X,
      knockoffs = replicate(10, matrix(rnorm(10000), ncol = 50), simplify = F),
      y = lapply( seq(ncol(X)), function(k) X[,k] ),
      ground_truth = seq(ncol(X))
    )
  })
})

test_that("Diagnostic catches heteroskedasticity", {
  set.seed(0)
  # Generate heteroskedastic data
  X = rbind(
    matrix(rnorm(1e4), ncol = 50),
    matrix(rnorm(1e4), ncol = 50)*0.1 + 5
  )
  for(k in seq(ncol(X))){
    X[,k] = X[,k] - mean(X[,k])
    X[,k] = X[,k] / sd(X[,k])
  }
  knockoffs = rlookc::computeGaussianKnockoffs(X, mu = 0, Sigma = cor(X), num_realizations = 24)
  reserved_knockoffs = knockoffs[1:5]
  knockoffs = knockoffs[-(1:5)]
  calibration_regular = rlookc::simulateY(X, knockoffs, reserved_knockoffs, FUN = function(x) x)
  calibration_nayusty = rlookc::simulateY(X, knockoffs, reserved_knockoffs, FUN = "adversarial", kmeans_centers = 2)
  calibration_diverse = rlookc::simulateY(X, knockoffs, reserved_knockoffs, FUN = "diverse")
  expect_lt(
    calibration_regular$calibration$fdr %>% colMeans %>% sum,
    calibration_nayusty$calibration$fdr %>% colMeans %>% sum
  )
  expect_lt(
    calibration_regular$calibration$fdr %>% colMeans %>% sum,
    calibration_diverse$calibration$fdr %>% colMeans %>% sum
  )
  data.frame(
    targeted_fdr = calibration_regular$calibration$targeted_fdrs,
    linear = calibration_regular$calibration$fdr %>% colMeans,
    adversarial = calibration_nayusty$calibration$fdr %>% colMeans,
    diverse_step = calibration_diverse$calibration$fdr %>% colMeans
  ) %>%
    tidyr::pivot_longer(cols = !targeted_fdr,
                        names_to = "Scheme",
                        values_to = "empirical_fdr") %>%
    ggplot() +
    geom_point(aes(x = targeted_fdr, y = empirical_fdr, colour = Scheme, shape = Scheme)) +
    geom_abline(aes(intercept=0, slope=1)) +
    ggtitle("Towards a more sensitive diagnostic")
})
