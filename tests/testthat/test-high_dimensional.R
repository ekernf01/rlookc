# Check mean
S = 2*rho*lambda
knockoff_mean_other = X - X %*% solve(Sigma)*S
plot(knockoff_mean_other, knockoff_mean); abline(a = 0, b = 1)
