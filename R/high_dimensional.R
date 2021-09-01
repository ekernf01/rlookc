#' Solve Ax=b fast.
#'
jacobi = function( A, b, n_step, w = 2/3 ){
  dA = diag(A)
  D = Matrix::Diagonal(x = dA)
  Dinv = Matrix::Diagonal(x = 1/dA)
  x = b
  for(k in seq(n_step)){
    current_b = A %*% x
    x = ( 1 - w ) * x + w * Dinv %*% (b - ( current_b - D %*% x ) )
  }
  x
}


#' Multiply by sqrt A fast.
#'
#' This implements the Euler method in eq. 3.6 of:
#'
#' Allen, E. J., Baglama, J., & Boyd, S. K. (2000). Numerical
#' approximation of the product of the square root of a matrix
#' with a vector. Linear Algebra and its Applications, 310(1-3),
#' 167-181.
#'
multiplySqrtA = function(A, x, n_step){
  dt = 1/n_step
  I = diag(nrow(A))
  for(k in seq(0, n_step-1)){
    kt = k*dt
    r = jacobi(
      I * ( 1 - kt ) + A * kt,
      ( I - A ) %*% x * 0.5,
      20)
    x = x + dt*r
  }
  x
}

