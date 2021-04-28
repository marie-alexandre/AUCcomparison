#' @title Spline Interpolation Method - Matrix of the zero order derivative coefficients
#' @description \loadmathjax In the area under the curve calculation using the spline interpolation method, the vector of the second derivative of the outcome of interest \mjseqn{Y} is expressed as \mjseqn{A Y^{''} = B Y  + F}. This function calculate calculate the matrix B.
#' 
#' @param time a numerical vector of time points of length m (x-axis coordinates).
#' 
#' @details The tridiagonal matrix \mjteqn{B}{B}{B} is defined as (for the "not-a-knot boundary conditions):
#' The \mjteqn{j}{j}{j}th line of the matrix, \mjteqn{B_{\[j,\ :\]}}{B_{\[j,\ :\]}}{B_{\[j,\ :\]}} is given by
#' \mjtdeqn{B_{\[j,\ :\]} = \left(0, \cdots, 0\right) \ if \ j=1}{B_{\[j,\ :\]} = \left(0, \cdots, 0\right) \ if \ j=1}{B_{\[j,\ :\]} = \left(0, \cdots, 0\right) \ if \ j=1}  
#' \mjtdeqn{B_{\[j,\ :\]} = \left(0, \cdots, 0\right) \ if \ j=m}{B_{\[j,\ :\]} = \left(0, \cdots, 0\right) \ if \ j=m}{B_{\[j,\ :\]} = \left(0, \cdots, 0\right) \ if \ j=m}
#' \mjtdeqn{B_{\[j,\ :\]} = \left(0_1, \cdots, 0_{j-2}, \frac{1}{h_j},-\left\[\frac{1}{h_j} + \frac{1}{h_{j+1}}\right\], \frac{1}{h_{j+1}}, 0_{j+2}, \cdots, 0_{m} \right) \ otherwise }{B_{\[j,\ :\]} = \left(0_1, \cdots, 0_{j-2}, \frac{1}{h_j},-\left\[\frac{1}{h_j} + \frac{1}{h_{j+1}}\right\], \frac{1}{h_{j+1}}, 0_{j+2}, \cdots, 0_{m} \right) \ otherwise }{B_{\[j,\ :\]} = \left(0_1, \cdots, 0_{j-2}, \frac{1}{h_j},-\left\[\frac{1}{h_j} + \frac{1}{h_{j+1}}\right\], \frac{1}{h_{j+1}}, 0_{j+2}, \cdots, 0_{m} \right) \ otherwise }

#' @return  a tridiagonal matrix corresponding to the weights of the variable of interest in the spline interpolation method. In this version, the matrix is build considering the "not-a-knot" spline boundary conditions.
#' @rdname AUC_Spline_matrix_B
#' @export 


AUC_Spline_matrix_B <- function(time){
  m <- length(time)
  hj <- NULL  # Warning: length(hj)=m-1
  for(j in 2:m){
    hj <- c(hj,(time[j]-time[j-1]))
  }
  Matrix_B <- matrix(0,ncol=m,nrow=m)
  
  for(j in 2:(m-1)){
    Matrix_B[j,j-1] <- 1/hj[j-1] 
    Matrix_B[j,j] <- -(1/hj[j-1]+1/hj[j])
    Matrix_B[j,j+1] <- 1/hj[j]
  }
  return(Matrix_B)
}
