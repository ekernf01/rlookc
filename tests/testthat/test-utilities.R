test_that("q-values match reference works", {
  set.seed(0)
  x = rnorm(1000) + 1
  for(threshold in (1:9)/10){
    expect_equal(
      knockoff.qvals(x) < 0.5,
      x > knockoff::knockoff.threshold(x, fdr = 0.5)
    )
  }
})

