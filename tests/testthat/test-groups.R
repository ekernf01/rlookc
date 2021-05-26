library("magrittr")
library("Matrix")
set.seed(0)

# Tests related to grouped knockoff selection

# Normal small example
easy_example = list(
  groups = list(1:2, 3:4, 5:6),
  X = matrix(rnorm(40*6), ncol = 6)
)
easy_example$Sigma = cor(easy_example$X)

# Hard example: nearly exact dupes of features
group_size = 5
num_groups = 50
num_obs = 1e3
hard_example = list(
  groups =  group_size %>% multiply_by(seq(num_groups)-1) %>% lapply(add, 1:group_size),
  X = matrix(0, nrow = num_obs, ncol = num_groups*group_size)
)
for(g in hard_example$groups){
  x = rnorm(num_obs)
  for(j in g){
    hard_example$X[,j] = x + rnorm(num_obs)*0.1
  }
}
hard_example$Sigma = matrix(0, nrow = num_groups*group_size, ncol = num_groups*group_size)
for(g in hard_example$groups){
  hard_example$Sigma[g,g] = 1/1.01
}
diag(hard_example$Sigma) = 1

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
testthat::test_that("Group selection controls the FDR", {
  num_knockoffs = 500
  X_ko_eachone  = lapply(1:num_knockoffs, function(i) knockoff::create.gaussian(X = hard_example$X, mu = 0, Sigma = hard_example$Sigma, method = "equi"))
  X_ko_grouped  = computeGaussianKnockoffs(X = hard_example$X, mu = 0, Sigma = hard_example$Sigma, num_realizations = num_knockoffs,
                                           groups = hard_example$groups, method = "group") %>% as.matrix
  for(signal in c("concentrated", "spread_out")){
    num_y = 500
    hard_example$active_group_idx    = sample(seq_along(hard_example$groups), replace = T, size = num_y)
    hard_example$active_variable_idx = hard_example$active_group_idx %>% lapply(function(g) hard_example$groups[[g]] )
    if(signal=="concentrated"){
      hard_example$active_variable_idx %<>% lapply(extract2, 1)
    }
    hard_example$y = lapply(hard_example$active_variable_idx, function(s) rowMeans(hard_example$X[,s, drop = F]) + 0.1*rnorm(num_obs))
    stats_eachone = sapply( seq(num_y), function( i ) marginal_screen( X = hard_example$X, X_k = X_ko_eachone[[i%%num_knockoffs + 1]], y = hard_example$y[[i]] ) ) %>% t
    stats_grouped = sapply( seq(num_y), function( i ) marginal_screen( X = hard_example$X, X_k = X_ko_grouped[[i%%num_knockoffs + 1]], y = hard_example$y[[i]] ) ) %>% t %>%
      aggregateStats(hard_example$groups)
    dir.create("~/Desktop/jhu/research/projects/rlookc/tests/grouping", recursive = T, showWarnings = F)
    calibration_eachone = check.calibration(ground_truth = hard_example$active_variable_idx,
                                            W = stats_eachone,
                                            plot_savepath = "~/Desktop/jhu/research/projects/rlookc/tests/grouping/individual_calibration.pdf")
    calibration_grouped = check.calibration(ground_truth = hard_example$active_group_idx,
                                            W = stats_grouped,
                                            plot_savepath = "~/Desktop/jhu/research/projects/rlookc/tests/grouping/grouped_calibration.pdf")
    # power comparison
    {
      pdf(paste0("~/Desktop/jhu/research/projects/rlookc/tests/grouping/power_grouped_vs_not__signal_type_", signal, ".pdf"))
      with(calibration_eachone,
           plot( x = colMeans(fdr),
                 y = colMeans(p_true_discoveries),
                 main = paste0("Power comparison\nsignal type: ", signal),
                 xlab = "empirical fdr", col ="red",
                 ylab = "True discoveries / active set size")
      )
      with(calibration_grouped,
           points( x = colMeans(fdr),
                   y = colMeans(p_true_discoveries), col = "blue")
      )
      legend(x = 0, y = 1, legend = c("individual", "grouped"), fill = c("red", "blue"))
      dev.off()
    }
  }
})

