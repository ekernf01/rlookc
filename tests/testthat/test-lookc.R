library("magrittr")
library("Matrix")

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
example_groups = c(list(1:2), as.list(3:10))

for(do_grouping in c(T, F)){
  if(do_grouping){
    example_groups = c(list(1:2), as.list(3:10))
    s_method = "group"
    reference_knockoff_maker = rlookc::computeGaussianKnockoffs
    plot(1,1, main = "Next plots: leave-one-out honesty \nchecks with grouped variables")
  } else {
    example_groups = NULL
    s_method = "sdp"
    reference_knockoff_maker = function(groups, ...) knockoff::create.gaussian(...)
    plot(1,1, main = "Next plots: leave-one-out honesty \nchecks without grouped variables")
  }
  # First test theoretical quantities
  # This is not the best software test because these are not used directly.
  # But it's a good check on the derivation.
  test_that("mean and covariance match in theory", {
    params_naive = rlookc::computeGaussianKnockoffs(Z, mu = 0, Sigma = Sigma,
                                                    output_type = "parameters",
                                                    groups = example_groups,
                                                    method = s_method)
    # For matching results, need to use the same S matrix for lookc and reference.
    # Don't let the reference code recompute it every time.
    same_s =  params_naive$S
    params_lookc = rlookc::generateLooks(           Z,
                                                    mu = 0,
                                                    Sigma = Sigma,
                                                    vars_to_omit = 1,
                                                    output_type = "parameters",
                                                    diag_s = same_s,
                                                    method = s_method,
                                                    groups = example_groups)
    params_sesia = rlookc::generateLooksSlow(       Z,
                                                    mu = 0,
                                                    Sigma = Sigma,
                                                    vars_to_omit = 1,
                                                    output_type = "parameters",
                                                    diag_s = same_s,
                                                    method = s_method,
                                                    groups = example_groups)
    i = cor(params_sesia[[1]][["ko_mean"]],
            params_lookc[[1]][["ko_mean"]]) %>% diag %>% which.min
    plot(params_sesia[[1]][["ko_mean"]][,i],
         params_lookc[[1]][["ko_mean"]][,i], main = "Worst match between means from \n reference and rlookc implementations")
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
  for(i in 1:500){
    if(i%%10==0){cat(".")}
    if(i%%100==0){cat("meow\n")}
    # Mathematically incorrect
    look_by_naive[[i]] = reference_knockoff_maker(Z,      mu = 0, Sigma = Sigma,
                                                   groups = example_groups, method = s_method)[, -1]
    # Mathematically correct, but slow
    look_by_sesia[[i]] = reference_knockoff_maker(Z[,-1], mu = 0, Sigma = Sigma[-1, -1],
                                                   groups = example_groups %>% removeKFromGroups(1), method = s_method)
    # Mathematically correct and fast, unless I fucked up
    look_by_lookc[[i]] = rlookc::generateLooks(Z, mu = 0, Sigma = Sigma, vars_to_omit = 1,
                                               groups = example_groups, method = s_method, output_type = "knockoffs")[[1]]
  }
  label_whether_grouped = function( plot_title ) paste0( plot_title, ifelse( do_grouping, " (grouped)", "" ) )
  test_that("lookc covariance matches reference better than naive covariance", {
    cov_naive = look_by_naive %>% lapply(cov) %>% Reduce("+", .) %>% divide_by(1000)
    cov_sesia = look_by_sesia %>% lapply(cov) %>% Reduce("+", .) %>% divide_by(1000)
    cov_lookc = look_by_lookc %>% lapply(cov) %>% Reduce("+", .) %>% divide_by(1000)
    plot(cov_sesia, cov_naive, main = "Cov: Naive vs actual" %>% label_whether_grouped)
    plot(cov_naive, cov_lookc, main = "Cov: Naive vs lookc"  %>% label_whether_grouped)
    plot(cov_sesia, cov_lookc, main = "Cov: Actual vs lookc" %>% label_whether_grouped)
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
    plot(m_sesia, m_naive, main = "Mean: Naive vs actual" %>% label_whether_grouped)
    plot(m_naive, m_lookc, main = "Mean: Naive vs lookc"  %>% label_whether_grouped)
    plot(m_sesia, m_lookc, main = "Mean: Actual vs lookc" %>% label_whether_grouped)
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
    boxplot(list("lookc" = corr_lookc, "naive" = corr_naive, "reference" = corr_sesia),
            main = "Correlation with the held-out variable \n across 1000 simulations" %>% label_whether_grouped)
    testthat::expect_lt(median(corr_sesia), median(corr_naive))
  })

}
