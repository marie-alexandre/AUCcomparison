#' @title Time Dependent Coefficients Cjp for AUC Lagrange Interpolation Method
#' @description \loadmathjax This function calculate the time-dependent coefficients Cjp involved in the calculation of the are under the curve when the Lagrange interpolation method is used.
#' 
#' @param ind_j a numerical scalar indicating the value of the index j.
#' @param ind_p a numerical scalar indicating the value of the index p
#' @param t a numerical vector of time points (x-axis coordinates) to consider for the AUC calculation.
#' 
#' @details The coefficients \mjteqn{C_{jp}}{C_{jp}}{C_{jp}} involved in the calculation of the AUC are defined as
#' \mjtdeqn{C_{2p} = (t_2-t_1)\prod_{l=0 ;\ l\neq p}^{P=2} t_{1+l} - \frac{(t_2^2-t_1^2)}{2}\sum_{l=0 ;\ l\neq p}^{P=2} t_{1+l} + \frac{(t_2^3-t_1^3)}{3} }{C_{2p} = (t_2-t_1)\prod_{l=0 ;\ l\neq p}^{P=2} t_{1+l} - \frac{(t_2^2-t_1^2)}{2}\sum_{l=0 ;\ l\neq p}^{P=2} t_{1+l} + \frac{(t_2^3-t_1^3)}{3} }{C_{2p} = (t_2-t_1)\prod_{l=0 ;\ l\neq p}^{P=2} t_{1+l} - \frac{(t_2^2-t_1^2)}{2}\sum_{l=0 ;\ l\neq p}^{P=2} t_{1+l} + \frac{(t_2^3-t_1^3)}{3} }
#' \mjtdeqn{C_{mp} = (t_m-t_{m-1})\prod_{l=0 ;\ l\neq p}^{P=2} t_{m-2+l} - \frac{(t_m^2-t_{m-1}^2)}{2}\sum_{l=0 ;\ l\neq p}^{P=2} t_{m-2+l} + \frac{(t_m^3-t_{m-1}^3)}{3} }{C_{mp} = (t_m-t_{m-1})\prod_{l=0 ;\ l\neq p}^{P=2} t_{m-2+l} - \frac{(t_m^2-t_{m-1}^2)}{2}\sum_{l=0 ;\ l\neq p}^{P=2} t_{m-2+l} + \frac{(t_m^3-t_{m-1}^3)}{3} }{C_{mp} = (t_m-t_{m-1})\prod_{l=0 ;\ l\neq p}^{P=2} t_{m-2+l} - \frac{(t_m^2-t_{m-1}^2)}{2}\sum_{l=0 ;\ l\neq p}^{P=2} t_{m-2+l} + \frac{(t_m^3-t_{m-1}^3)}{3} }
#' \mjtdeqn{C_{jp} = -(t_j-t_{j-1})\prod_{l=0 ;\ l\neq p}^{P=3} t_{j-2+l} + \frac{(t_j^2-t_{j-1}^2)}{2}\sum_{l_1=0 ;\ l_1\neq p}^{P-1=2}\sum_{l_2=l_1+1 ;\ l_2\neq p}^{P=3}t_{j-2+l_1}\cdot t_{j-2+l_2} - }{C_{jp} = -(t_j-t_{j-1})\prod_{l=0 ;\ l\neq p}^{P=3} t_{j-2+l} + \frac{(t_j^2-t_{j-1}^2)}{2}\sum_{l_1=0 ;\ l_1\neq p}^{P-1=2}\sum_{l_2=l_1+1 ;\ l_2\neq p}^{P=3}t_{j-2+l_1}\cdot t_{j-2+l_2} - }{C_{jp} = -(t_j-t_{j-1})\prod_{l=0 ;\ l\neq p}^{P=3} t_{j-2+l} + \frac{(t_j^2-t_{j-1}^2)}{2}\sum_{l_1=0 ;\ l_1\neq p}^{P-1=2}\sum_{l_2=l_1+1 ;\ l_2\neq p}^{P=3}t_{j-2+l_1}\cdot t_{j-2+l_2} - }
#' \mjtdeqn{\frac{(t_j^3-t_{j-1}^3)}{3} \sum_{l=0 ;\ l\neq p}^{P=3} t_{j-2+l} + \frac{(t_j^4-t_{j-1}^4)}{4}}{\frac{(t_j^3-t_{j-1}^3)}{3} \sum_{l=0 ;\ l\neq p}^{P=3} t_{j-2+l} + \frac{(t_j^4-t_{j-1}^4)}{4}}{\frac{(t_j^3-t_{j-1}^3)}{3} \sum_{l=0 ;\ l\neq p}^{P=3} t_{j-2+l} + \frac{(t_j^4-t_{j-1}^4)}{4}}
#' where \mjteqn{m}{m}{m} is the number of time points in the vector \code{t}.
#' 
#' @return a numerical scalar corresponding to the coefficient \mjteqn{C_{jp}}{C_{jp}}{C_{jp}} evaluated for j = \code{ind_j} and p = \code{ind_p}.
#' 
#' @rdname AUC_Lagrange_Cjp_coefficients
#' @export 

AUC_Lagrange_Cjp_coefficients <- function(ind_j,ind_p,t){
  if(ind_j == 2){
    term1 <- (t[ind_j] - t[ind_j-1])*prod(sapply(seq(0,2),function(l){
      if(l!=ind_p){
        return(t[1+l])
      }else{
        return(1)
      }
    },simplify=TRUE))
    term2 <- (t[ind_j]^2 - t[ind_j-1]^2)/2*sum(sapply(seq(0,2),function(l){
      if(l!=ind_p){
        return(t[1+l])
      }else{
        return(0)
      }
    },simplify=TRUE))
    term3 <- (t[ind_j]^3 - t[ind_j-1]^3)/3
    Cjp <- term1 - term2 + term3
  }else if(ind_j == length(t)){
    term1 <- (t[ind_j] - t[ind_j-1])*prod(sapply(seq(0,2),function(l){
      if(l!=ind_p){
        return(t[ind_j-2+l])
      }else{
        return(1)
      }
    },simplify=TRUE))
    term2 <- (t[ind_j]^2 - t[ind_j-1]^2)/2*sum(sapply(seq(0,2),function(l){
      if(l!=ind_p){
        return(t[ind_j-2+l])
      }else{
        return(0)
      }
    },simplify=TRUE))
    term3 <- (t[ind_j]^3 - t[ind_j-1]^3)/3
    Cjp <- term1 - term2 + term3
  }else{
    term1 <- (t[ind_j] - t[ind_j-1])*prod(sapply(seq(0,3),function(l){
      if(l!=ind_p){
        return(t[ind_j-2+l])
      }else{
        return(1)
      }
    },simplify=TRUE))
    term2 <- (t[ind_j]^2 - t[ind_j-1]^2)/2*sum(sapply(seq(0,2),function(l1){
      if(l1!=ind_p){
        res <- 0
        for(l2 in seq(l1+1,3)){
          if(l2!=ind_p){
            res <- res + t[ind_j-2+l1]*t[ind_j-2+l2]
          }
        }
        return(res)
      }else{
        return(0)
      }
    },simplify=TRUE))
    term3 <- (t[ind_j]^3 - t[ind_j-1]^3)/3*sum(sapply(seq(0,3),function(l){
      if(l!=ind_p){
        return(t[ind_j-2+l])
      }else{
        return(0)
      }
    },simplify=TRUE))
    term4 <- (t[ind_j]^4 - t[ind_j-1]^4)/4
    Cjp <- -term1 + term2 - term3 + term4
  }
  return(Cjp)
}
