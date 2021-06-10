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
    calibration = simulateY(X = matrix(rnorm(1000), ncol = 10),
                            knockoffs = matrix(rnorm(1000), ncol = 10))
  })

})
