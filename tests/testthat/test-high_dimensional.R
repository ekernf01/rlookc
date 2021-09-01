

test_that("jacobi works", {
  # Set up a linear system with an SPD matrix
  set.seed(0)
  U = rnorm(5e4) %>% matrix(ncol = 100)
  A = cov(U)
  x = rnorm(100)
  b = A %*% x
  # Solve it, hopefully
  plot(x, jacobi(A, b, 20), main = "Testing my hack-job weighted Jacobi iteration")
  abline(a = 0, b = 1, col = "red")
  expect_equal(x, as.numeric(jacobi(A, b, 20)), tolerance = 1e-3)
})

test_that("multiplySqrtA works", {
  # Set up an SPD matrix
  set.seed(0)
  U = rnorm(5e2) %>% matrix(ncol = 10)
  A = cov(U)
  eigenA = eigen(A)
  sqrtA = eigenA$vectors %*% diag(sqrt(eigenA$values)) %*% t(eigenA$vectors)
  max(abs(sqrtA %*% sqrtA - A))
  I = diag(nrow(A))
  supposedlySqrtA = multiplySqrtA(A, I, 1000)
  plot(sqrtA, supposedlySqrtA)
  plot(A, supposedlySqrtA %*% supposedlySqrtA )
  expect_equal(sqrtA, as.matrix(supposedlySqrtA), tolerance = 1e-3)
})
