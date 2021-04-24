#' @title Weights for AUC Matrix Formulation
#' @description \loadmathjax In matrix formulation, the area under a curve of interest, named \emph{Y}, can be expressed as matrix product of a vector of weights \emph{W} and the vector of the values of \emph{Y}. This function calculates the weights \emph{W} when AUC is calculated either by the trapezoid, the Lagrange or the Spline interpolation methods.
#' 
#' @param time a numerical vector of time points of length m (x-axis coordinates for AUC calculation).
#' @param method a character scalar indicating the interpolation method of interest. Options are 'trapezoid', 'lagrange' and 'spline'. In this version the 'spline' interpolation method is implemented with the "not-a-knot" spline boundary conditions. 
#' 
#' @details In matrix formulation, the AUC of the outcome \mjseqn{Y} can be expressed as \mjseqn{\text{AUC} = W \cdot Y}, with \mjseqn{W} defined by the following expressions for the trapezoid, the Lagrange and the spline interpolation methods.
#' 
#' **Trapezoid method:** 
#' \mjsdeqn{W_j = \frac{t_{j+1} - t_j}{2} \text{ if } j=1} 
#' \mjsdeqn{W_j = \frac{t_{j} - t_{j-1}}{2} \text{ if } j=m} 
#' \mjsdeqn{W_j = \frac{t_{j+1} - t_{j-1}}{2} \text{ otherwise}} 
#'
#' **Lagrange method:** (see \code{\link[AUCcomparison]{AUC_Lagrange_Cjp_coefficients}} for the definition of the Cjp coefficients)
#' \mjsdeqn{W_j = \frac{C_{\[2\]\[j-1\]}}{\prod_{l=0 ;\ l\neq (j-1)}^{P=2} (t_j-t_{j+1})} + \sum_{p=0}^{P=3} \frac{C_{\[j-1+p\]\[3-p\]}}{\prod_{l=0 ;\ l\neq (3-p)}^{P=3} (t_j-t_{j-3+p+l})} \text{ if } j=1,2,3} 
#' \mjsdeqn{W_j = \frac{C_{\[m\]\[j-(m-2)\]}}{\prod_{l=0 ;\ l\neq (j-(m-2))}^{P=2} (t_j-t_{j-2+l})} + \sum_{p=0}^{m-j} \frac{C_{\[j-1+p\]\[3-p\]}}{\prod_{l=0 ;\ l\neq (3-p)}^{P=3} (t_j-t_{j-3+p+l})} \text{ if } j=m-2,m-1,m} 
#' \mjsdeqn{W_j = \sum_{p=0}^{m-j} \frac{C_{\[j-1+p\]\[3-p\]}}{\prod_{l=0 ;\ l\neq (3-p)}^{P=3} (t_j-t_{j-3+p+l})} \text{ othermise}} 
#' 
#' **Spline method:**  (see \code{\link[AUCcomparison]{AUC_Spline_matrix_A}} and \code{\link[AUCcomparison]{AUC_Spline_matrix_B}} for the definition of Matrices A and B)
#'  \mjsdeqn{W_j = \sum_{p=2}^m -\frac{(t_p-t_{p-1})^3}{24}(u_{pj}+u_{p-1j}) + W_j^{trap.}}
#'  where \mjseqn{(u_{pj})} is the element \mjseqn{U(p,j)} with \mjseqn{U} a matrix defined as \mjseqn{U = A^{-1}B}.
#'  
#' 
#' 
#' @return A numerical scalar with same length than the vector \code{time} corresponding to the weights \emph{W}.
#' 
#' @rdname AUC_time_weights_estimation
#' @export 
#' @importFrom matlib


AUC_time_weights_estimation <- function(time,method){
  Weights <- NULL
  mg <- length(time)
  if(method == "trapezoid"){
    Weights <- sapply(seq(1,mg),function(j){
      if(j==1){
        res <- (time[j+1]-time[j])/2
      }else if(j == mg){
        res <- (time[j]-time[j-1])/2
      }else{
        res <- (time[j+1]-time[j-1])/2
      }
      return(res)
    },simplify=TRUE)
  }else if(method == "lagrange"){
    Weights <- sapply(seq(1,mg),function(j){
      if(j %in% c(1,2,3)){
        term1 <- AUC_Lagrange_Cjp_coefficients(ind_j=2,ind_p=j-1,t=time)/prod(sapply(seq(0,2),function(l){
          if(l!=(j-1)){
            return(time[j]-time[1+l])
          }else{
            return(1)
          }
        },simplify=TRUE))
        term2 <- sum(sapply(seq(3-(j-1),3),function(p){
          coeff <- AUC_Lagrange_Cjp_coefficients(ind_j=j-1+p,ind_p=3-p,t=time)
          den_prod <- prod(sapply(seq(0,3),function(l){
            if(l!=(3-p)){
              return(time[j]-time[j-3+p+l])
            }else{
              return(1)
            }
          },simplify=TRUE))
          return(coeff/den_prod)
        },simplify=TRUE))
        res <- term1 + term2
      }else if(j %in% c(mg-2,mg-1,mg)){
        term1 <- AUC_Lagrange_Cjp_coefficients(ind_j=mg,ind_p=j-(mg-2),t=time)/prod(sapply(seq(0,2),function(l){
          if(l!=(j-(mg-2))){
            return(time[j]-time[mg-2+l])
          }else{
            return(1)
          }
        },simplify=TRUE))
        term2 <- sum(sapply(seq(0,mg-j),function(p){
          coeff <- AUC_Lagrange_Cjp_coefficients(ind_j=j-1+p,ind_p=3-p,t=time)
          den_prod <- prod(sapply(seq(0,3),function(l){
            if(l!=(3-p)){
              return(time[j]-time[j-3+p+l])
            }else{
              return(1)
            }
          },simplify=TRUE))
          return(coeff/den_prod)
        },simplify=TRUE))
        res <- term1 + term2
      }else{
        res <- sum(sapply(seq(0,3),function(p){
          coeff <- AUC_Lagrange_Cjp_coefficients(ind_j=j-1+p,ind_p=3-p,t=time)
          den_prod <- prod(sapply(seq(0,3),function(l){
            if(l!=(3-p)){
              return(time[j]-time[j-3+p+l])
            }else{
              return(1)
            }
          },simplify=TRUE))
          return(coeff/den_prod)
        },simplify=TRUE))
      }
      return(res)
    },simplify=TRUE)
  }else if(method == "spline"){
    w_trap <- AUC_time_weights_estimation(time=time,method="trapezoid")
    mat_A <- AUC_Spline_matrix_A(time=time)
    mat_B <- AUC_Spline_matrix_B(time=time)
    mat_U <- inv(mat_A)%*%mat_B
    Weights <- sapply(seq(1,mg),function(j){
      res <- sum(sapply(seq(2,mg),function(p){
        tmp <- -(time[p]-time[p-1])^3/24*(mat_U[p,j]+mat_U[p-1,j])
        return(tmp)
      },simplify=TRUE))
    },simplify=TRUE) + w_trap 
  }
  return(Weights)
}
