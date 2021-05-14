library("magrittr")

# Tests related to leave-one-out knockoff construction.


# Example where knockoffs look really different without v1
Z = matrix(NA, ncol = 10, nrow = 100)
Z[,1] = rnorm(1e2)
for(k in 2:10){
  Z[,k] = rnorm(1e2)
}
Z[,1] = 3*rowMeans(Z[,-1]) + rnorm(1e2)
Z = Z %*% diag(1/apply(Z, 2, sd))
Sigma = cov(Z)
Sigma %>% image
Sigma %>% solve %>% heatmap(symm = T, col = colorRampPalette(c("red", "white", "blue"))(25))
Sigma_inv = Sigma %>% solve

# First test theoretic quantities
# This is not the best software test because these are not used directly.
# But it's a good check on the derivation.
test_that("mean and covariance match in theory", {
  params_naive = rlookc::computeGaussianKnockoffs(Z, mu = 0, Sigma = Sigma,                   output_type = "parameters")
  same_s =  diag(params_naive$S)
  params_lookc = rlookc::generateLooks(           Z, mu = 0, Sigma = Sigma, vars_to_omit = 1, output_type = "parameters", diag_s = same_s)
  params_sesia = rlookc::generateLooksSlow(       Z, mu = 0, Sigma = Sigma, vars_to_omit = 1, output_type = "parameters", diag_s = same_s)

  i = cor(params_sesia[[1]][["ko_mean"]],
          params_lookc[[1]][["ko_mean"]]) %>% diag %>% which.min
  plot(params_sesia[[1]][["ko_mean"]][,i],
      params_lookc[[1]][["ko_mean"]][,i], main = "Worst match between means of \n reference and looks implementations")
  abline(a=0, b=1)
  testthat::expect_equal( params_sesia[[1]][["ko_mean"]],
                          params_lookc[[1]][["ko_mean"]], tolerance = 1E-2 )
  testthat::expect_equal( params_sesia[[1]][["ko_sqrt_covariance"]] %>% crossprod(),
                          params_lookc[[1]][["ko_covariance"]], tolerance = 1E-2 )
})
# Monte carlo estimate of LOOK covariance and mean
# It's random, so gotta run lots of trials
look_by_sesia =list()
look_by_lookc =list()
look_by_naive =list()
for(i in 1:1000){
  if(i%%10==0){cat(".")}
  if(i%%100==0){cat("meow\n")}
  look_by_naive[[i]] = knockoff::create.gaussian(Z,      mu = 0, Sigma = Sigma)[, -1] #Mathematically incorrect
  look_by_sesia[[i]] = knockoff::create.gaussian(Z[,-1], mu = 0, Sigma = Sigma[-1, -1]) #Mathematically correct
  look_by_lookc[[i]] = rlookc::generateLooks(Z, mu = 0, Sigma = Sigma, vars_to_omit = 1, output_type = "knockoffs")[[1]] #Mathematically correct unless I fucked up
}

test_that("lookc covariance matches reference better than naive covariance", {
  cov_naive = look_by_naive %>% lapply(cov) %>% Reduce("+", .) %>% divide_by(1000)
  cov_sesia = look_by_sesia %>% lapply(cov) %>% Reduce("+", .) %>% divide_by(1000)
  cov_lookc = look_by_lookc %>% lapply(cov) %>% Reduce("+", .) %>% divide_by(1000)
  plot(cov_sesia, cov_naive, main = "Naive vs actual")
  plot(cov_naive, cov_lookc, main = "Naive vs lookc")
  plot(cov_sesia, cov_lookc, main = "Actual vs lookc")
  cor(cov_sesia, cov_naive) %>% diag %>% sum
  cor(cov_naive, cov_lookc) %>% diag %>% sum
  cor(cov_sesia, cov_lookc) %>% diag %>% sum
  testthat::expect_lt(
    cor(cov_naive, cov_lookc) %>% diag %>% sum,
    cor(cov_sesia, cov_lookc) %>% diag %>% sum
    )
})

test_that("lookc mean matches reference better than naive mean", {
  m_naive = look_by_naive %>% Reduce("+", .) %>% divide_by(1000)
  m_sesia = look_by_sesia %>% Reduce("+", .) %>% divide_by(1000)
  m_lookc = look_by_lookc %>% Reduce("+", .) %>% divide_by(1000)
  plot(m_sesia, m_naive, main = "Naive vs actual")
  plot(m_naive, m_lookc, main = "Naive vs lookc")
  plot(m_sesia, m_lookc, main = "Actual vs lookc")
  cor(m_sesia, m_naive) %>% diag %>% sum
  cor(m_naive, m_lookc) %>% diag %>% sum
  cor(m_sesia, m_lookc) %>% diag %>% sum
  testthat::expect_lt(
    cor(m_naive, m_lookc) %>% diag %>% sum,
    cor(m_sesia, m_lookc) %>% diag %>% sum
  )
})

test_that("lookc predicts left-out variable way worse than naive does", {
  # given the simulation, the best predictor should be the sum of the knockoffs.
  corr_naive = look_by_naive %>% lapply(rowSums) %>% sapply(cor, Z[, 1])
  corr_sesia = look_by_sesia %>% lapply(rowSums) %>% sapply(cor, Z[, 1])
  corr_lookc = look_by_lookc %>% lapply(rowSums) %>% sapply(cor, Z[, 1])
  boxplot(list("lookc" = corr_lookc, "naive" = corr_naive, "sesia" = corr_sesia),
          main = "Correlation with the held-out variable \n across 1000 simulations")
  testthat::expect_lt(median(corr_sesia), median(corr_naive))
})


