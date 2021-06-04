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

