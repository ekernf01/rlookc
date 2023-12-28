## ----setup, include=FALSE----------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ----------------------------------------------------------------------------------------------------------------------------------
library("magrittr")
library("knockoff")
library("rlookc")
library("ggplot2")


## ----------------------------------------------------------------------------------------------------------------------------------
set.seed(0)
p = 25
X = matrix(NA, ncol = p, nrow = 1e4)
X[,1] = rnorm(1e4)
for(k in 2:p){
  X[,k] = rnorm(1e4)
  if(k-100 < 4){
    X[,k] = (X[,k] + X[,k-1]) / sqrt(2)
  }
}


## ----------------------------------------------------------------------------------------------------------------------------------
Sigma = cov(X)
mu = rep(0, ncol(X))
stats = rlookc::create__looks( X, mu = mu, Sigma = Sigma, output_type = "statistics")


## ----------------------------------------------------------------------------------------------------------------------------------
stats_tidy = matrix(0, p,p)
for(k in 1:p){
  stats_tidy[-k, k] = stats[[k]]
}
stats_tidy = stats_tidy %>% 
  as.data.frame %>% 
  set_colnames(1:p) %>%
  set_rownames(1:p) %>%
  tibble::rownames_to_column("rowname") %>%
  tidyr::pivot_longer(cols = !matches("rowname")) %>%
  dplyr::mutate(rowname = as.integer(rowname), name = as.integer(name)) %>%
  dplyr::mutate(true_signal = abs(rowname - name) == 1)
ggplot(stats_tidy) + geom_histogram(aes(x = value, fill = true_signal), bins = 100)


## ----------------------------------------------------------------------------------------------------------------------------------
stats_tidy$q = rlookc::calibrate__getQvals(stats_tidy$value)
stats_tidy %<>% dplyr::arrange(q)
stats_tidy$observed_fdr = cumsum(!stats_tidy$true_signal) / (1:(p^2))
ggplot(stats_tidy) + geom_tile(aes(x = name, y=  rowname, fill = q))
ggplot(stats_tidy) + 
  geom_point(aes(x = q, y=  observed_fdr), color = "blue") + 
  geom_abline(intercept=0, slope=1)

