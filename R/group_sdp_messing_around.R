# Sigma = cor(matrix(rnorm(3e3), ncol = 3))
# stopifnot(isSymmetric(Sigma))
# G = cov2cor(Sigma)
# p = dim(G)[1]
# if (!is_posdef(G)) {
#   warning("The covariance matrix is not positive-definite: knockoffs may not have power.",
#           immediate. = T)
# }
# Cl1 = rep(0, p)
# Al1 = -Matrix::Diagonal(p)
# Cl2 = rep(1, p)
# Al2 = Matrix::Diagonal(p)
# d_As = c(diag(p))
# As = Matrix::Diagonal(length(d_As), x = d_As)
# As = As[which(Matrix::rowSums(As) > 0), ]
# Cs = c(2 * G)
# A = cbind(Al1, Al2, As)
# C = matrix(c(Cl1, Cl2, Cs), 1)
# K = NULL
# K$s = p
# K$l = 2 * p
# b = rep(1, p)
# OPTIONS = NULL
# OPTIONS$gaptol = gaptol
# OPTIONS$maxit = maxit
# OPTIONS$logsummary = 0
# OPTIONS$outputstats = 0
# OPTIONS$print = 0
# if (verbose)
#   cat("Solving SDP ... ")
# sol = Rdsdp::dsdp(A, b, C, K, OPTIONS)
