library("rlookc")
library("magrittr")
library("knockoff")

# This vignette tests knockoff error control in subset selection for a vector autoregressive model.
# Since there is weird overlap between "Y" and "X" in an autoregressive model, it is unclear (to me) if knockoffs
# retain their ability to control FDR in subset selection. I chose to investigate this via simulation.

# VAR model parameters
n_time = 1000
n_var = 100
nnz = 10
n_rep = 100
A = matrix(0, n_var, n_var)
active_sets = list()
for(i in seq(n_var)){
  active_sets[[i]] = sample(seq(n_var), size = nnz, replace = F) %>% union(i)
  A[i, active_sets[[i]]] = 1
}
# For each experiment, make data, make knockoffs, use knockoffs.
X = knockoffs = stats = list()
for(j in seq(n_rep)){
  X[[j]] = matrix(0, nrow = n_time, ncol = n_var)
  X[[j]][1, ] = rnorm(n_var)
  for(i in seq(2,n_time)){
    X[[j]][i, ] = (0.5/n_var)*A %*% X[[j]][i-1, ] + 0.5*X[[j]][i-1, ] + rnorm( n_var )
  }
  knockoffs[[j]] = knockoff::create.second_order(X[[j]])
  # I test only 99 variables. This is solely to check the dimensions in my code.
  # Otherwise it's 100x100 and I risk using the transpose by accident.
  stats[[j]] =
    cor( X[[j]][2:n_time,2:n_var],         X[[j]][1:(n_time-1),] ) -
    cor( X[[j]][2:n_time,2:n_var], knockoffs[[j]][1:(n_time-1),] )
}

dir.create("tests/autoregressive", recursive = T, showWarnings = F)
# This plot has formally correct CI's.
# Each datum is draw from a separate trial.
calibration_results = checkCalibration(
  ground_truth = active_sets[2] %>% rep(n_rep),
  W = stats %>% lapply(t) %>% sapply(extract, 2, ) %>% t,
  n_var = n_var - 1,
  plot_savepath = "tests/autoregressive/FDR control one variable.pdf"
)
# This plot uses all the data but CI's might be off because the knockoffs are re-used.
calibration_results2 = checkCalibration(
  ground_truth = active_sets %>% rep(n_rep),
  W = stats %>% lapply(t) %>% lapply(as.data.frame) %>% data.table::rbindlist() %>% as.matrix(),
  plot_savepath = "tests/autoregressive/FDR control all variables.pdf"
)
