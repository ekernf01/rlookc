# This vignette demonstrates the new features of our package that
# are related to checking FDR calibration and model assumptions.
library("magrittr")

# Using calibrate__simulateY, you can simulate from P(Y|X) and check any set of knockoffs you generate.
# The choice of P(Y|X) has a big effect on the results.
# When E[Y|X] is linear, this diagnostic has almost no power, as you might expect
# if you know a lot about fixed-X knockoffs. This is shown in "calibration_check3" below.
# We try to choose a sensible family.
# The default is step functions of individual features, indicating when they cross 0.
# But you can choose anything and provide it as the FUN argument.
set.seed(0)
make_bimodal = function(X) rbind(X, X + 10)
X = rnorm(1e4) %>% matrix(ncol = 10)
good_knockoffs = knockoff::create.second_order(X)
X %<>% make_bimodal
good_knockoffs %<>% make_bimodal
bad_knockoffs = knockoff::create.second_order(X)
calibration_check = rlookc::calibrate__simulateY(X, good_knockoffs)
calibration_check2 = rlookc::calibrate__simulateY(X, bad_knockoffs)
calibration_check3 = rlookc::calibrate__simulateY(X, bad_knockoffs, FUN = function(x) 2*x)
plot(   calibration_check$calibration$targeted_fdrs, calibration_check$calibration$fdr %>% colMeans, main = "Calibration")
points(calibration_check2$calibration$targeted_fdrs, calibration_check2$calibration$fdr %>% colMeans, col = "blue", pch = "+")
points(calibration_check3$calibration$targeted_fdrs, calibration_check3$calibration$fdr %>% colMeans, col = "red", pch = "*")
abline(a = 0, b = 1)

# You can use calibrate__chooseDiverseY and rlookc::calibrate__findWorstY to try
# to search more thoroughly for a worst-case P(Y|X).
# rlookc::calibrate__findWorstY will split the data, searching for the worst P(Y|X)
# in one half and checking calibration in the other half.
collection_of_y = rlookc::calibrate__chooseDiverseY(X, n_quantiles = 5)
worst_y = rlookc::calibrate__findWorstY(X, good_knockoffs,
                                       y = collection_of_y$y,
                                       ground_truth = collection_of_y$ground_truth)
worst_y = rlookc::calibrate__findWorstY(X, bad_knockoffs,
                                       y = collection_of_y$y,
                                       ground_truth = collection_of_y$ground_truth)
# Once these have run, you can compare the worst-case Y to another Y or otherwise
# try to determine informative properties of the worst-case Y.
collection_of_y$ground_truth[[worst_y$worst_y]]
cor( X, collection_of_y$y[[worst_y$worst_y]])
cor( bad_knockoffs, collection_of_y$y[[worst_y$worst_y]])
plot(X[,2], collection_of_y$y[[worst_y$worst_y]], main = "Worst Y")
plot(X[,2], collection_of_y$y[[6]], main = "Aparently not as bad")

# Using calibrate__KNNTest, you can check exchangeability for any set of knockoffs you generate.
# This test is less customizable than the above but we find it's much more sensitive.
# It is sensitive to any deviation from exchangeability without the need
# for a careful choice of P(Y|X).
set.seed(0)
X = rnorm(1e4) %>% matrix(ncol = 10)
A = rnorm(1e2) %>% matrix(ncol = 10) %>% add(1) # to induce known correlation across variables
X = X %*% A
good_knockoffs = knockoff::create.gaussian(X, mu = 0, Sigma = t(A) %*% A)
bad_knockoffs = apply(X, 2, sample)
rlookc::calibrate__KNNTest(X, good_knockoffs)[2:3]
rlookc::calibrate__KNNTest(X, bad_knockoffs)[2:3]
# You can even pinpoint outliers. (This trick can have weak signal, though.)
X[1:50,1:5] = -100
rlookc::calibrate__KNNTest(X, good_knockoffs, n_neighbors = 100)$prop_not_swapped_per_observation %>%
  head(100) %>%
  plot
# The conservative action to take is to use exact dupes of the data for
# any observations you think don't fit the model for P(X).
# This will restore exchangeability, possibly with a loss of power.
good_knockoffs_outlier_dupe = good_knockoffs
good_knockoffs_outlier_dupe[1:50,] = X[1:50,]
rlookc::calibrate__KNNTest(X, good_knockoffs_outlier_dupe, n_neighbors = 100)$prop_not_swapped_per_observation %>%
  head(100) %>%
  plot


