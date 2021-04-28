#' @title Variance of the Difference of AUC of Two Group-Specific Polynomial Marginal Dynamics
#' @description \loadmathjax This function calculates the variance of the difference of area under the curve of two marginal dynamics modeled by group-structured polynomials or B-spline curve in Mixed-Effects models.
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
#' @param time.G1 a numerical vector of time points (x-axis coordinates) to use for the variance of the Group1 AUC calculation.
#' @param time.G2 a numerical vector of time points (x-axis coordinates) to use for the variance of the Group2 AUC calculation.
#' @param common.interval a logical scalar. If FALSE, the variance of difference of AUC is calculated as the variance of the difference of AUCs where the AUC of each group is calculated on its specific interval of time. If TRUE (default), the variance is estimated on a common interval of time defined as the intersect of the two group-specific interval (see \link[AUCcomparison]{Group_specific_Delta_AUC_estimation} for more details about calculation.).
#' @param method  a character scalar indicating the interpolation method to use to estimate the AUC. Options are 'trapezoid' (default), 'lagrange' and 'spline'. In this version, the 'spline' interpolation is implemented with "not-a-knot" spline boundary conditions.
#' @param Group.dependence a logical scalar indicating whether the two groups, whose the difference of AUC (\mjteqn{\Delta AUC}{\Delta AUC}{\Delta AUC}) is studied, are considered as dependent. By default, this variable is defined as TRUE.
#' @param Averaged a logical scalar. If TRUE, the function return the difference of normalized AUC (nAUC) where nAUC is computed as the AUC divided by the range of time of calculation. If FALSE (default), the classic AUC is calculated.
#' 
#' @return A numerical scalar corresponding to the variance of the difference of AUC (\mjteqn{\Delta AUC}{\Delta AUC}{\Delta AUC}) between the Group1 and the Group2. If the two groups are considered as dependent (\code{Group.dependence}=TRUE), the variance of \mjteqn{\Delta AUC}{\Delta AUC}{\Delta AUC} is calculated as \mjteqn{Var(AUC_1) + Var(AUC_2) - 2Cov(AUC_1,AUC_2)}{Var(AUC_1) + Var(AUC_2) - 2Cov(AUC_1,AUC_2)}{Var(AUC_1) + Var(AUC_2) - 2Cov(AUC_1,AUC_2)}. Otherwise, only the sum of the two variance is used.
#'
#' @seealso 
#'  \code{\link[splines]{bs}}, 
#'  \code{\link[AUCcomparison]{Group_specific_Var_AUC_estimation}}, 
#'  \code{\link[AUCcomparison]{MEM_Polynomial_Group_structure}}
#'  
#' @examples 
#' \donttest{# Download of data
#' data("HIV_Simu_Dataset_Delta01_cens")
#' data <- HIV_Simu_Dataset_Delta01_cens
#' 
#' # Change factors in character vectors
#' data$id <- as.character(data$id) ; data$Group <- as.character(data$Group)
#' 
#' # Example 1: We consider the variable \code{MEM_Pol_Group} as the output of our function 
#' # \link[AUCcomparison]{MEM_Polynomial_Group_structure}
#' MEM_estimation <- MEM_Polynomial_Group_structure(y=data$VL,x=data$time,Group=data$Group,
#'                                                  Id=data$id,Cens=data$cens)
#'
#' time_group1 <- unique(data$time[which(data$Group=="Group1")])
#' time_group2 <- unique(data$time[which(data$Group=="Group2")])
#' 
#' Var_Delta_AUC_estimation <- Group_specific_Var_Delta_AUC_estimation(
#'                                               MEM_Pol_group=MEM_estimation,
#'                                               Group1="Group1",Group2="Group2",
#'                                               time.G1=time_group1,time.G2=time_group2)
#'                                                                     
#' # Example 2: We consider results of MEM estimation from another source. 
#' # We have to give build the variable 'MEM_Pol_group' with the good structure
#' # We build the variable 'MEM_Pol_group.1' with the results of MEM estimation obtained for 2 groups 
#' # Generation of random matrix
#' Covariance_Matrix_1 <- matrix(rnorm(7*7,mean=0,sd=0.01),ncol=7,nrow=7) 
#' # Transform the matrix into symmetric one
#' Covariance_Matrix_1 <- Covariance_Matrix_1 %*% t(Covariance_Matrix_1) 
#' MEM_Pol_group.1 <- list(Model_estimation=Covariance_Matrix_1,
#'                        Model_features=list(Groups=c("Group1","Group2"),
#'                               Marginal.dyn.feature=list(dynamic.type="polynomial",
#'                                               intercept=c(global.intercept=TRUE,
#'                                               group.intercept1=FALSE,group.intercept2=FALSE),
#'                                               polynomial.degree=c(3,3))))
#'                                            
#' Var_Delta_AUC_estimation_2 <- Group_specific_Var_Delta_AUC_estimation(
#'                                               MEM_Pol_group=MEM_Pol_group.1,
#'                                               Group1="Group1",Group2="Group2",
#'                                               time.G1=time_group1,time.G2=time_group2)
#'                                                                     
#'}
#' 
#' @rdname Group_specific_Var_Delta_AUC_estimation
#' @export 
#' @importFrom ArgumentCheck newArgCheck addError finishArgCheck addWarning addMessage
#' @importFrom splines bs

Group_specific_Var_Delta_AUC_estimation <- function(MEM_Pol_group,Group1,Group2,time.G1,time.G2,common.interval=TRUE,method="trapezoid",Group.dependence=TRUE,Averaged=FALSE){
  
  # Step 1: Verification of the type of arguments
  # ----- #
  # Verification of 'common.interval'
  Check_commonTime <- ArgumentCheck::newArgCheck()
  if(isFALSE(is.logical(common.interval))){
    ArgumentCheck::addError(
      msg = "Error - The variable 'common.interval' must be a boolean",
      argcheck = Check_commonTime
    )
  }
  ArgumentCheck::finishArgCheck(Check_commonTime)
  # Verification of 'Group.dependence'
  Check_group.dep <- ArgumentCheck::newArgCheck()
  if(isFALSE(is.logical(Group.dependence))){
    ArgumentCheck::addError(
      msg = "The argument 'Group.dependence' must be defined as boolean variable",
      argcheck = Check_group.dep
    )
  }
  ArgumentCheck::finishArgCheck(Check_group.dep)
  # Verification of Group1 and Group2 not null
  Check_groups_null <-  ArgumentCheck::newArgCheck()
  if(is.null(Group1) || is.null(Group2)){
    ArgumentCheck::addError(
      msg = "One of the two Groups to compared has been assigned to 'NULL' value",
      argcheck = Check_groups_null
    )
  }
  ArgumentCheck::finishArgCheck(Check_groups_null)
  
  Groups <- c(Group1,Group2)
  
  if(common.interval == TRUE){
    min.time.interval <- max(min(time.G1,na.rm=TRUE),min(time.G2,na.rm=TRUE),na.rm=TRUE)
    max.time.interval <- min(max(time.G1,na.rm=TRUE),max(time.G2,na.rm=TRUE),na.rm=TRUE)
    if(min.time.interval > max.time.interval){
      stop("Impossible to estimate the difference of AUC on the common interval of time for the two groups: the inrestection of the two intervals is equal to an empty interval.")
    }else if(min.time.interval == max.time.interval){
      stop(paste("Impossible to estimate the difference of AUC on the common interval of time for the two groups: the inrestection of the two intervals is equal to a single time point:",min.time.interval))
    }else{
      time.G1 <- time.G1[which(time.G1 >= min.time.interval & time.G1 <= max.time.interval)]
      time.G2 <- time.G2[which(time.G2 >= min.time.interval & time.G2 <= max.time.interval)]
    }
  }
  
  time <- list(time.G1,time.G2)
  Check_argument_Group_specific_Var_AUC(MEM_Pol_group,time,Groups,method,Averaged)
  
  Model_features <- MEM_Pol_group$Model_features
  Marginal_dynamics <- Model_features$Marginal.dyn.feature
  
  # Extraction of population parameters according to their groups
  if(is.list(MEM_Pol_group$Model_estimation)){
    Population_variance <- MEM_Pol_group$Model_estimation$varFix
  }else{
    Population_variance <- MEM_Pol_group$Model_estimation
  }
  MEM_groups <- as.vector(Model_features$Groups)
  
  global_intercept <- Marginal_dynamics$intercept["global.intercept"]
  ind_params <- 0
  Group_index_params <- list()
  for(g in 1:length(MEM_groups)){
    params <- NULL
    index_params <- NULL
    if(global_intercept){
      index_params <- c(index_params,1)  
      if(g == 1){
        ind_params <- ind_params + 1
      }
    }
    if(Marginal_dynamics$dynamic.type == "spline"){
      Nb_group_params <- as.numeric(1*Marginal_dynamics$intercept[paste("group.intercept",g,sep="")] + 
                                      length(Marginal_dynamics$knots[[MEM_groups[g]]]) + Marginal_dynamics$spline.degree[g])
    }else if(Marginal_dynamics$dynamic.type == "polynomial"){
      Nb_group_params <- as.numeric(1*Marginal_dynamics$intercept[paste("group.intercept",g,sep="")] + 
                                      Marginal_dynamics$polynomial.degree[g])
    }
    index_params <- c(index_params,seq(ind_params+1,ind_params+Nb_group_params))
    ind_params <- ind_params + Nb_group_params
    Group_index_params[[MEM_groups[g]]] <- index_params
  }
  Groups_index_params <- unique(sort(as.numeric(unlist(Group_index_params))))
  Groups_Variance <- Population_variance[Groups_index_params,Groups_index_params]
  
  # Step 2: Calculation of the variance of Delta AUC
  # ----- #
  if(Group.dependence){
    Combined_Pop_Covariate <- NULL
    Pop_Covariate <- list()
    for(g in 1:length(Groups)){
      time_group <- time[[g]]
      
      if(Marginal_dynamics$dynamic.type == "polynomial"){
        # Creation of covariate matrix
        Covariate_group <- do.call(cbind,lapply(1*isFALSE(Marginal_dynamics$intercept[paste("group.intercept",g,sep="")]):Marginal_dynamics$polynomial.degree[g],function(d) time_group^d))
      }else if(Marginal_dynamics$dynamic.type == "spline"){
        # Creation of covariate matrix
        if(is.null(Marginal_dynamics$boundary.knots[[Groups[g]]])){
          Covariate_group <- splines::bs(x=time_group,knots=Marginal_dynamics$knots[[Groups[g]]],df=Marginal_dynamics$df[g],degree=Marginal_dynamics$spline.degree[g])
        }else{
          Covariate_group <- splines::bs(x=time_group,knots=Marginal_dynamics$knots[[Groups[g]]],df=Marginal_dynamics$df[g],degree=Marginal_dynamics$spline.degree[g],Boundary.knots=Marginal_dynamics$boundary.knots[[Groups[g]]])
        }
        if(Marginal_dynamics$intercept[paste("group.intercept",g,sep="")]){
          Covariate_group <- cbind(rep(1,length(time_group)),Covariate_group)
        }
      } # End spline covariate
      Pop_Covariate[[Groups[g]]] <- Covariate_group
    }
    Combined_Pop_Covariate <- matrix(0,ncol=ncol(Pop_Covariate[[1]]) + ncol(Pop_Covariate[[2]]),
                                     nrow=nrow(Pop_Covariate[[1]]) + nrow(Pop_Covariate[[2]]))
    Combined_Pop_Covariate[1:nrow(Pop_Covariate[[1]]),1:ncol(Pop_Covariate[[1]])] <- Pop_Covariate[[1]]
    Combined_Pop_Covariate[(nrow(Pop_Covariate[[1]])+1):nrow(Combined_Pop_Covariate),(ncol(Pop_Covariate[[1]])+1):ncol(Combined_Pop_Covariate)] <- Pop_Covariate[[2]]
    if(global_intercept){
      Combined_Pop_Covariate <- cbind(rep(1,nrow(Combined_Pop_Covariate)),Combined_Pop_Covariate)
    }
    
    time_weights.G1 <- AUC_time_weights_estimation(time=time[[1]],method)
    time_weights.G2 <- AUC_time_weights_estimation(time=time[[2]],method)
    if(Averaged){
      Global_time_weights <- c(rep(0,length(time_weights.G1)),time_weights.G2)/(diff(range(time.G2))) -c(time_weights.G1,rep(0,length(time_weights.G2)))/(diff(range(time.G1)))
    }else{
      Global_time_weights <- c(rep(0,length(time_weights.G1)),time_weights.G2) -c(time_weights.G1,rep(0,length(time_weights.G2)))
    }
    Var_Delta_AUC <- as.numeric(Global_time_weights %*% Combined_Pop_Covariate %*% Groups_Variance %*% t(Combined_Pop_Covariate) %*% Global_time_weights)
  }else{
    Var_AUC_Group1 <- Group_specific_Var_AUC_estimation(MEM_Pol_group=MEM_Pol_group,time=time.G1,Groups=Group1,method=method,Averaged=Averaged)
    Var_AUC_Group2 <- Group_specific_Var_AUC_estimation(MEM_Pol_group=MEM_Pol_group,time=time.G2,Groups=Group2,method=method,Averaged=Averaged)
    Var_Delta_AUC <- Var_AUC_Group1 + Var_AUC_Group2
  }
  return(Var_Delta_AUC)
}
