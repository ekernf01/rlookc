% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calibrate.R
\name{calibrate__findWorstY}
\alias{calibrate__findWorstY}
\title{Given X and various simulated Y|X, find the worst-calibrated Y|X.}
\usage{
calibrate__findWorstY(
  X,
  X_k,
  y,
  ground_truth,
  split = rep(c(F, T), length.out = nrow(X)),
  statistic = knockoff::stat.glmnet_lambdasmax,
  plot_savepath = NULL,
  ...
)
}
\arguments{
\item{X}{@param X_k A real dataset and a corresponding model-X knockoff realization.}

\item{y}{@param ground_truth Two lists of the same length, one containing y's that you simulated
and another containing the indices of the variables in X used to simulate them.}

\item{split}{to get an unbiased estimate of the calibration for the worst-calibrated P(Y|X),
half the data are used for choosing the worst and the other half for estimating its calibration.}

\item{statistic}{Function used to compute variable importance, e.g. knockoff::stat.glmnet_lambdasmax.}

\item{plot_savepath}{Passed to check.calibration}

\item{...}{Passed to statistic, e.g. n_lambda=100.}
}
\description{
Given X and various simulated Y|X, find the worst-calibrated Y|X.
}
