# Sim data
X = matrix(rnorm(300), nrow = 10)
Sigma = cor(X)*(1-lambda) + lambda*diag(ncol(X))

# Check mean
S = 2*rho*lambda
knockoff_mean_other = X - X %*% solve(Sigma)*S
plot(knockoff_mean_other, knockoff_mean); abline(a = 0, b = 1)

# Check covariance
S = 2*rho*lambda
knockoff_mean_other = X - X %*% solve(Sigma)*S
plot(knockoff_mean_other, knockoff_mean); abline(a = 0, b = 1)
