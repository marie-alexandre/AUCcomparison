#' @title T-statistic of the Difference of AUC of Two Group-Specific Polynomial Marginal Dynamics
#' @description This function performs the t-test evaluating whether the difference of area under the curve of two marginal dynamics, modeled by group-structured polynomials or B-spline curve in Mixed-Effects model, is null. 
#' 
#' @param MEM_Pol_group A list with similar structure than the output provided by the function \link[AUCcomparison]{MEM_Polynomial_Group_structure}. 
#' 
#' A list containing:
#' \itemize{
#' \item \code{Model_estimation}: a list containing at least 2 elements: \enumerate{
#'     \item the vector of the marginal (fixed) parameters estimates (at least for the groups whose AUC is to estimate), labeled _'beta'_.
#'     \item the variance-covariance matrix of these parameters, labeled _'varFix'_ (see \link[AUCcomparison]{MEM_Polynomial_Group_structure} for details about the parameter order).
#'     }
#' \item \code{Model_features}: a list of at least 2 elements: \enumerate{
#'     \item \code{Groups}: a vector indicating the names of the groups whose fixed parameters are given.
#'     \item \code{Marginal.dyn.feature}: a list summarizing the features of the marginal dynamics defined in the model: 
#'     \itemize{
#'         \item \code{dynamic.type}: a character scalar indicating the chosen type of marginal dynamics. Options are 'polynomial' or 'spline'.
#'         \item \code{intercept}: a logical vector summarizing choices about global and group-specific intercepts (Number of groups + 1) elements whose elements are named as ('global.intercept','group.intercept1', ..., 'group.interceptG') if G Groups are defined in \code{MEM_Pol_group}. For each element of the vector, if TRUE, the considered intercept is considered as included in the model (see \emph{Examples}).
#'         }
#'      If \code{dynamic.type} is defined as 'polynomial':\itemize{
#'          \item \code{polynomial.degree}: an integer vector indicating the degree of polynomial functions, one value for each group.
#'          }
#'      If \code{dynamic.type} is defined as 'spline':\itemize{
#'          \item \code{spline.degree}: an integer vector indicating the degree of B-spline curves, one for each group. 
#'          \item \code{knots}: a list of group-specific internal knots used to build B-spline basis (one numerical vector for each group) (see \link[splines]{bs} for more details).
#'          \item \code{df}: a numerical vector of group-specific degrees of freedom used to build B-spline basis, (one for each group).
#'          \item \code{boundary.knots}: a list of group-specific boundary knots used to build B-spline  basis (one vector for each group) (see \link[splines]{bs} for more details).
#'          }
#'     }
#' }
#' @param Group1 a character scalar indicating the name of the first group whose marginal dynamics must be considered. This group name must belong to the set of groups involved in the MEM (see \code{Groups} vector in \code{MEM_Pol_group}). 
#' @param Group2 a character scalar indicating the name of the second group whose marginal dynamics must be considered. This group name must belong to the set of groups involved in the MEM (see \code{Groups} vector in \code{MEM_Pol_group}).
#' @param time.G1 a numerical vector of time points (x-axis coordinates) to use for the Group1 AUC calculation.
#' @param time.G2 a numerical vector of time points (x-axis coordinates) to use for the Group2 AUC calculation.
#' @param common.interval a logical scalar. If FALSE, the difference of AUC is calculated as the difference of AUCs where the AUC of each group is calculated on its specific interval of time. If TRUE (default), the difference of AUC is estimated on a common interval of time defined as the intersect of the two group-specific interval (see \link[AUCcomparison]{Group_specific_Delta_AUC_estimation} for more details about calculation).
#' @param method a character scalar indicating the interpolation method to use to estimate the AUC. Options are 'trapezoid' (default), 'lagrange' and 'spline'. In this version, the 'spline' interpolation is implemented with "not-a-knot" spline boundary conditions.
#' @param Group.dependence a logical scalar indicating whether the two groups, whose the difference of AUC (\eqn{\Delta}AUC) is studied, are considered as dependent. By default, this variable is defined as TRUE.
#' @param Averaged a logical scalar. If TRUE, the function return the difference of normalized AUC (nAUC) where nAUC is computated as the AUC divided by the range of time of calculation. If FALSE (default), the classic AUC is calculated.
#' @param conf_level a numerical value (between 0 and 1) indicating the confidence level of the interval. By default, this variable is fixed at 0.95
#' @param alternative a character scalar specifying the alternative hypothesis. Options are 'two.sided' (default), 'greater' or 'less'. 
#' 
#' 
#' @return A list containing:\itemize{
#' \item \code{Tstat}: the value of the t-statistic.
#' \item \code{Pvalue}:  the P-value.
#' \item \code{Conf.int}: the confidence interval.
#' \item \code{Delta_AUC}: the estimated value of the difference of AUC between the two groups (nAUC2 - nAUC1) (see \code{\link[AUCcomparison]{Group_specific_Delta_AUC_estimation}} for more details).
#' \item \code{AUCs}: the estimated values of the Group-specific AUC (AUC1 and AUC 2)  
#' (see \code{\link[AUCcomparison]{Group_specific_AUC_estimation}} for more details).
#' }

#' @examples 
#' \donttest{# Download of data
#' data("HIV_Simu_Dataset_Delta01_cens")
#' data <- HIV_Simu_Dataset_Delta01_cens
#' 
#' # Change factors in character vectors
#' data$id <- as.character(data$id) ; data$Group <- as.character(data$Group)
#' 
#' MEM_estimation <- MEM_Polynomial_Group_structure(y=data$VL,x=data$time,Group=data$Group,
#'                                                  Id=data$id,Cens=data$cens)
#'                                                  
#' time_group1 <- unique(data$time[which(data$Group=="Group1")])
#' time_group2 <- unique(data$time[which(data$Group=="Group2")])
#' Test <- Stat_test_Delta_AUC_Group_Specific(MEM_Pol_group=MEM_estimation,
#'                                            Group1="Group1",Group2="Group2",
#'                                            time.G1=time_group1,time.G2=time_group2)
#'}
#' 
#' @seealso 
#'  \code{\link[AUCcomparison]{MEM_Polynomial_Group_structure}},
#'  \code{\link[AUCcomparison]{Group_specific_Delta_AUC_estimation}},
#'  
#'  \code{\link[AUCcomparison]{Group_specific_AUC_estimation}}
#' @rdname Stat_test_Delta_AUC_Group_Specific
#' @export 
#' @importFrom stats pnorm qnorm

Stat_test_Delta_AUC_Group_Specific <- function(MEM_Pol_group,Group1,Group2,time.G1,time.G2,common.interval=TRUE,method="trapezoid",Group.dependence=TRUE,Averaged=FALSE,conf_level=0.95,alternative="two.sided"){
  # Estimation of the Difference of AUC
  Delta_AUC <- Group_specific_Delta_AUC_estimation(MEM_Pol_group=MEM_Pol_group,Group1=Group1,Group2=Group2,time.G1=time.G1,time.G2=time.G2,common.interval=common.interval,method=method,Averaged=Averaged)
  # Estimation of its variance
  Var_Delta_AUC <- Group_specific_Var_Delta_AUC_estimation(MEM_Pol_group=MEM_Pol_group,Group1=Group1,Group2=Group2,time.G1=time.G1,time.G2=time.G2,method=method,Group.dependence=Group.dependence,Averaged=Averaged)
  # Estimation of the AUC
  if(common.interval == TRUE){
    min.time.interval <- max(min(time.G1,na.rm=TRUE),min(time.G2,na.rm=TRUE),na.rm=TRUE)
    max.time.interval <- min(max(time.G1,na.rm=TRUE),max(time.G2,na.rm=TRUE),na.rm=TRUE)
    time.G1 <- time.G1[which(time.G1 >= min.time.interval & time.G1 <= max.time.interval)]
    time.G2 <- time.G2[which(time.G2 >= min.time.interval & time.G2 <= max.time.interval)]
  }
  AUC_G1 <- Group_specific_AUC_estimation(MEM_Pol_group=MEM_Pol_group,Groups=Group1,time=time.G1,method=method,Averaged=Averaged)
  AUC_G2 <- Group_specific_AUC_estimation(MEM_Pol_group=MEM_Pol_group,Groups=Group2,time=time.G2,method=method,Averaged=Averaged)
  AUCs <- c(AUC_G1,AUC_G2)
  names(AUCs) <- c(Group1,Group2)
  
  # Estimation of the statistic
  Tstat <- (Delta_AUC)/sqrt(Var_Delta_AUC)
  # Estimation of the Pvalue
  if(alternative == "two.sided"){
    Pvalue <- 2*stats::pnorm(abs(Tstat),mean=0,sd=1,lower.tail = FALSE)
    alpha <- (1-conf_level)/2
  }else if(alternative == "greater"){
    Pvalue <- 1-stats::pnorm(abs(Tstat),mean=0,sd=1,lower.tail = FALSE)
    alpha <- 1-conf_level
  }else{
    Pvalue <- stats::pnorm(abs(Tstat),mean=0,sd=1,lower.tail = FALSE)
    alpha <- 1-conf_level
  }
  
  # Estimation of the confidence interval
  t.theoric <- stats::qnorm(alpha,mean=0,sd=1,lower.tail = FALSE)*(alternative%in%c("two.sided","greater")) - stats::qnorm(alpha,mean=0,sd=1,lower.tail = FALSE)*(alternative%in%c("less"))
  
  Confidence.Interval <- sort(c(Delta_AUC-t.theoric*sqrt(Var_Delta_AUC), Delta_AUC+t.theoric*sqrt(Var_Delta_AUC)))
  attr(Confidence.Interval,'conf.level') <- conf_level
  
  Results <- list(Tstat=Tstat,Pvalue=Pvalue,Conf.int=Confidence.Interval,Delta_AUC=Delta_AUC,AUCs=AUCs) 
  return(Results)
}
