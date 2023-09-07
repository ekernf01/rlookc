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


test_that("Simulations handle failures of user-provided alternative_method", {
  expect_warning({
    calibration = calibrate__simulateY(X = matrix(rnorm(1000), ncol = 20),
                            alternative_method = function(y, X) stop("nyah nyah"),
                            n_sim = 5,
                            shuddup = T)
  })
})

test_that("Simulations run", {
  expect_invisible({
    calibration = calibrate__simulateY(X = matrix(rnorm(1000), ncol = 20),
                            knockoffs = replicate(5, matrix(rnorm(1000), ncol = 20), simplify = F),
                            n_sim = 50,
                            shuddup = T)
  })
  expect_invisible({
    calibration = calibrate__simulateY(X = matrix(rnorm(1000), ncol = 20),
                            knockoffs = NULL,
                            alternative_method = function(y, X) runif(ncol(X)),
                            n_sim = 50,
                            shuddup = T)
  })
  X = matrix(rnorm(1000), ncol = 20)
  expect_invisible({
    calibration = findWorstY(
      X,
      X_k = matrix(rnorm(1000), ncol = 20),
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
  X_k = rlookc::computeGaussianKnockoffs(X, mu = 0, Sigma = cor(X))
  diverse_y = chooseDiverseY(X)
  worst_y_calibration = findWorstY(
    X,
    X_k,
    y = diverse_y$y,
    ground_truth = diverse_y$ground_truth
  )
  plot(worst_y_calibration$calibration$targeted_fdr, colMeans(worst_y_calibration$calibration$fdr))
})

test_that("KNN test has correct null distribution", {
  set.seed(0)
  outs = list()
  for(i in 1:1000){
    X   = rnorm(1000) %>% matrix(ncol = 10)
    X_k = rnorm(1000) %>% matrix(ncol = 10)
    outs[[i]] = KNNTest(X, X_k)
  }
  outs %>% sapply(magrittr::extract2, "prop_not_swapped") %>% summary
  outs %>% sapply(magrittr::extract2, "prop_not_swapped") %>% hist(40)
  outs %>%
    sapply(magrittr::extract2, "p_value") %>%
    sort %>%
    plot(main ="KNN test calibration",
         xlab = "rank",
         ylab = "null p")
  abline(a = 0, b = 1e-3, col = "red")
})

