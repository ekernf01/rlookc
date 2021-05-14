library("magrittr")
library("Matrix")

# Tests related to grouped knockoff selection

# Normal small example
easy_example = list(
  groups = list(1:2, 3:4, 5:6),
  X = matrix(rnorm(40*6), ncol = 6)
)
easy_example$Sigma = cor(easy_example$X)

# Hard example: nearly exact dupes of features
hard_example = list(
  groups =  10 %>% multiply_by(0:49) %>% lapply(add, 1:10),
  X = matrix(0, nrow = 2e3, ncol = 500)
)
for(g in hard_example$groups){
  x = rnorm(2e3)
  for(j in g){
    hard_example$X[,j] = x + rnorm(2e3)*0.1
  }
}
hard_example$Sigma = cor(hard_example$X)

# Make sure choice of S is valid
testthat::test_that("Joint cov of X and knockoffs is SPD", {
  for( e in list(easy_example, hard_example)){
    testthat::expect_silent({
      S = solveGroupEqui(e$Sigma, e$groups)
      G = rbind(
        cbind( e$Sigma , e$Sigma - S ),
        cbind( e$Sigma - S , e$Sigma )
      )
      R = chol(G)
    })
  }
})

# Test if group-level FDR is honest 
# Also, does it beat the singly-generated knockoffs in power?
easy_example$active_group_idx = sample(seq_along(easy_example$groups), replace = T, size = 2500)
hard_example$active_group_idx = sample(seq_along(hard_example$groups), replace = T, size = 2500)
easy_example$y = lapply(easy_example$active_group_idx, function(g) easy_example$X[,easy_example$groups[[g]][[1]]] + 0.1*rnorm(40))
hard_example$y = lapply(hard_example$active_group_idx, function(g) hard_example$X[,hard_example$groups[[g]][[1]]] + 0.1*rnorm(2e3))
testthat::test_that("Group selection controls the FDR", {
  for( e in list(easy_example, hard_example)){
    X_ko_eachone  = computeGaussianKnockoffs(X = e$X, mu = 0, Sigma = e$Sigma, method = "equi")
    X_ko_grouped  = computeGaussianKnockoffs(X = e$X, mu = 0, Sigma = e$Sigma, groups = e$groups, method = "group") %>% as.matrix
    stats_eachone = sapply( e$y, function( y ) marginal_screen( X = e$X, X_k = X_ko_eachone, y = y ) ) %>% t
    stats_grouped = sapply( e$y, function( y ) marginal_screen( X = e$X, X_k = X_ko_grouped, y = y ) ) %>% t %>%
      aggregateStats(e$groups)
    calibration_eachone = check.calibration(ground_truth = e$groups[e$active_group_idx] %>% lapply(extract2, 1),  
                                            W = stats_eachone, 
                                            plot_savepath = "tests/individual_calibration.pdf")
    calibration_grouped = check.calibration(ground_truth = e$active_group_idx,           
                                            W = stats_grouped,
                                            plot_savepath = "tests/group_calibration.pdf")
  }
})



# TODO: finish the math for grouped LOOKs
# TODO: Test if LOOKs still control the error when grouping is present, if the derivation works out
