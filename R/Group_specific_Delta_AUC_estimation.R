#' @title Difference of AUC of Two Group-Specific Polynomial Marginal Dynamics
#' @description \loadmathjax This function estimates the difference of area under the curve of marginal dynamics from two groups when marginal dynamics are modeled by group-structured polynomials or B-spline curves.

#' @param MEM_Pol_group A list with similar structure than the output provided by the function \link[AUCcomparison]{MEM_Polynomial_Group_structure}. 
#' 
#' A list containing: \tabular{ll}{
#' \code{Model_estimation} \tab either the vector of the marginal parameter estimates (at least for the groups \code{Group1} and \code{Group2}), or a list containing this vector labeled _'beta'_ among other arguments (see \link[AUCcomparison]{MEM_Polynomial_Group_structure} for details about the parameter order). \cr
#' \code{Model_features} \tab a list of at least 2 elements: \cr
#' \tab 1. \code{Groups}  -  a vector indicating the names of the groups whose fixed parameters are given.  \cr
#' \tab 2. \code{Marginal.dyn.feature}  -  a list summarizing the features of the marginal dynamics defined in the model:  \cr
#' \tab \itemize{
#' \item \code{dynamic.type} - a character scalar indicating the chosen type of marginal dynamics. Options are 'polynomial' or 'spline'
#' \item \code{intercept} -  a logical vector summarizing choices about global and group-specific intercepts (Number of groups + 1) elements whose elements are named as ('global.intercept','group.intercept1', ..., 'group.interceptG') if G Groups are defined in \code{MEM_Pol_group}. For each element of the vector, if TRUE, the considered intercept is considered as included in the model. 
#'  
#' If \code{dynamic.type} is defined as 'polynomial': 
#' \item \code{polynomial.degree} - an integer vector indicating the degree of polynomial functions, one value for each group. 
#' 
#' If \code{dynamic.type} is defined as 'spline':
#' \item \code{spline.degree} - an integer vector indicating the degree of B-spline curves, one for each group. 
#' \item \code{knots} - a list of group-specific internal knots used to build B-spline basis (one numerical vector for each group) (see \link[splines]{bs} for more details).
#' \item \code{df} - a numerical vector of group-specific degrees of freedom used to build B-spline basis, (one for each group). 
#' \item \code{boundary.knots} - a list of group-specific boundary knots used to build B-spline  basis (one vector for each group) (see \link[splines]{bs} for more details).
#' } \cr
#' }
#' 
#' @param Group1 a character scalar indicating the name of the first group whose marginal dynamics must be considered. This group name must belong to the set of groups involved in the MEM (see \code{Groups} vector in \code{MEM_Pol_group}). 
#' @param Group2 a character scalar indicating the name of the second group whose marginal dynamics must be considered. This group name must belong to the set of groups involved in the MEM (see \code{Groups} vector in \code{MEM_Pol_group}).
#' @param time.G1 a numerical vector of time points (x-axis coordinates) to use for the Group1 AUC calculation.
#' @param time.G2 a numerical vector of time points (x-axis coordinates) to use for the Group2 AUC calculation.
#' @param method a character scalar indicating the interpolation method to use to estimate the AUC. Options are 'trapezoid' (default), 'lagrange' and 'spline'. In this version, the 'spline' interpolation is implemented with "not-a-knot" spline boundary conditions.
#' @param Averaged a logical scalar. If TRUE, the function return the difference of normalized AUC (nAUC) where nAUC is computated as the AUC divided by the range of time of calculation. If FALSE (default), the classic AUC is calculated.
#' @return A numerical scalar defined as \mjseqn{\Delta \text{AUC} = \text{AUC}_2 - \text{AUC}_1} (or \mjseqn{\Delta\text{nAUC} = \text{nAUC}_2 - \text{nAUC}_1})  with \mjseqn{\text{AUC}_1} (or \mjseqn{\text{nAUC}_1}) and  \mjseqn{\text{AUC}_2} (or \mjseqn{\text{nAUC}_2}) being respectively estimated as the AUC (or nAUC) for the Group1 and for the Group2.
#' 
#' @examples 
#' # Download of data
#' data("HIV_Simu_Dataset_Delta01_cens")
#' data <- HIV_Simu_Dataset_Delta01_cens
#' 
#' # Change factors in character vectors
#' data$id <- as.character(data$id) ; data$Group <- as.character(data$Group)
#' 
#' # Example 1: We consider the variable \code{MEM_Pol_Group} as the output of our function \link[AUCcomparison]{MEM_Polynomial_Group_structure}
#' MEM_estimation_1 <- MEM_Polynomial_Group_structure(y=data$VL,x=data$time,Group=data$Group,Id=data$id,Cens=data$cens)
#' Delta_AUC_1 <- Group_specific_Delta_AUC_estimation(MEM_Pol_group=MEM_estimation_1,Group1="Group1",Group2="Group2",
#'                                                    time.G1=unique(data$time[which(data$Group=="Group1")]),
#'                                                    time.G2=unique(data$time[which(data$Group=="Group2")]))
#'
#' Example 2: We consider results of MEM estimation from another source. We have to give build the variable 'MEM_Pol_group' with the good structure
#'  # We build the variable 'MEM_estimation_2' with the results of MEM estimation obtained for two groups 
#' MEM_estimation_2 <- list(Model_estimation=c(1.077,0.858,-0.061,0.0013,0.887,-0.066,0.0014), # c(global.intercept,beta1.G1,beta2.G1,beta2.G1,beta1.G2,beta2.G2,beta3.G2)
#'                         Model_features=list(Groups=c("Group1","Group2"),
#'                                             Marginal.dyn.feature=list(dynamic.type="polynomial",intercept=c(global.intercept=TRUE,group.intercept1=FALSE,group.intercept2=FALSE),polynomial.degree=c(3,3))))
#' Delta_AUC_2 <- Group_specific_Delta_AUC_estimation(MEM_Pol_group=MEM_estimation_2,Group1="Group1",Group2="Group2",
#'                                                    time.G1=unique(data$time[which(data$Group=="Group1")]),
#'                                                    time.G2=unique(data$time[which(data$Group=="Group2")]))

#' @seealso 
#'  \code{\link[AUCcomparison]{MEM_Polynomial_Group_structure}},
#'  \code{\link[AUCcomparison]{Group_specific_AUC_estimation}}
#' @rdname Group_specific_Delta_AUC_estimation
#' @export 

Group_specific_Delta_AUC_estimation <- function(MEM_Pol_group,Group1,Group2,time.G1,time.G2,method="trapezoid",Averaged=FALSE){
  AUC_Group1 <- Group_specific_AUC_estimation(MEM_Pol_group=MEM_Pol_group,time=time.G1,Groups=Group1,method=method,Averaged=Averaged) 
  AUC_Group2 <- Group_specific_AUC_estimation(MEM_Pol_group=MEM_Pol_group,time=time.G2,Groups=Group2,method=method,Averaged=Averaged) 
  Delta_AUC <- as.numeric(AUC_Group2 - AUC_Group1)
  return(Delta_AUC)
}
