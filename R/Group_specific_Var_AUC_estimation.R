#' @title Variance of the Area Under The Curve of Group-Specific Polynomial Marginal Dynamics
#' @description This function calculates the variance of the area under the curve of marginal dynamics modeled by group-structured polynomials or B-spline curves in Mixed-Effects models
#' 
#' @param MEM_Pol_group A list with similar structure than the output provided by the function \link[AUCcomparison]{MEM_Polynomial_Group_structure}. 
#' 
#' A list containing: \tabular{ll}{
#' \code{Model_estimation} \tab either the variance-covariance matrix of the marginal (fixed) parameters (at least for the groups whose Variance of AUC is to estimate), or a list containing at least this matrix labeled _'varFix'_ (see \link[AUCcomparison]{MEM_Polynomial_Group_structure} for details about the parameter order). \cr
#' \code{Model_features} \tab a list of at least 2 elements: \cr
#' \tab 1. \code{Groups}  -  a vector indicating the names of the groups whose variance of fixed parameters are given.\cr
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
#' @param time a numerical vector of time points (x-axis coordinates) or a list of numerical vectors (with as much elements than the number of groups in \code{Groups}).
#' @param Groups  a vector indicating the names of the groups belonging to the set of groups involved in \code{MEM_Pol_group} for which we want to estimate the AUC  (a subset or the entire set of groups involved in the model can be considered). If NULL (default), the AUC for all the groups involved the MEM is calculated.
#' @param method a character scalar indicating the interpolation method to use to estimate the AUC. Options are 'trapezoid' (default), 'lagrange' and 'spline'. In this version, the 'spline' interpolation is implemented with the "not-a-knot" spline boundary conditions. 
#' @param Averaged a logical scalar. If TRUE, the function return the normalized AUC (nAUC) computed as the AUC divided by the range of the time calculation. If FALSE (default), the classic AUC is calculated.
#' 
#' @return A numerical vector containing the estimation of the variance of the AUC (or nAUC) for each group defined in the \code{Groups} vector.
#' 
#' @seealso 
#'  \code{\link[splines]{bs}}, 
#'  \code{\link[AUCcomparison]{MEM_Polynomial_Group_structure}}
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
#' MEM_estimation <- MEM_Polynomial_Group_structure(y=data$VL,x=data$time,Group=data$Group,Id=data$id,Cens=data$cens)
#' Var_AUC_estimation <- Group_specific_Var_AUC_estimation(MEM_Pol_group=MEM_estimation,time=list(unique(data$time[which(data$Group == "Group1")]),
#'                                                          unique(data$time[which(data$Group == "Group2")])))
#'                                                          
#' Example 2: We consider results of MEM estimation from another source. We have to give build the variable 'MEM_Pol_group' with the good structure
#'  # We build the variable 'MEM_Pol_group.1' with the results of MEM estimation obtained for two groups 
#' Covariance_Matrix_1 <- matrix(rnorm(7*7,mean=0,sd=0.01),ncol=7,nrow=7) # Generation of random matrix
#' Covariance_Matrix_1 <- Covariance_Matrix_1 %*% t(Covariance_Matrix_1) # Transform the matrix into symmetric one
#' MEM_Pol_group.1 <- list(Model_estimation=Covariance_Matrix_1, # Covariance matrix of fixed effects for all parameters
#'                        Model_features=list(Groups=c("Group1","Group2"),
#'                                            Marginal.dyn.feature=list(dynamic.type="polynomial",intercept=c(global.intercept=TRUE,group.intercept1=FALSE,group.intercept2=FALSE),polynomial.degree=c(3,3))))
#'                                            
#' Var_AUC_estimation_G1.1 <- Group_specific_Var_AUC_estimation(MEM_Pol_group.1,time=unique(data$time[which(data$Group == "Group1")]),Groups=c("Group1"))
#' 
#' # We build the variable 'MEM_Pol_group.2' with the results of MEM estimation obtained only for the group of interest (extraction)
#' Covariance_Matrix_2 <-  matrix(rnorm(4*4,mean=0,sd=0.01),ncol=4,nrow=4) # Generation of random matrix
#' Covariance_Matrix_2 <- Covariance_Matrix_2 %*% t(Covariance_Matrix_2) # Transform the matrix into a symetric one
#' MEM_Pol_group.2 <- list(Model_estimation=Covariance_Matrix_2, # Covariance matrix of fixed effects, only for the parameters from Group1
#'                        Model_features=list(Groups=c("Group1"),
#'                                            Marginal.dyn.feature=list(dynamic.type="polynomial",intercept=c(global.intercept=TRUE,group.intercept1=FALSE),polynomial.degree=c(3))))
#'                                            
#' Var_AUC_estimation_G1.2 <- Group_specific_Var_AUC_estimation(MEM_Pol_group=MEM_Pol_group.2,time=unique(data$time[which(data$Group == "Group1")]))  

#'  
#' @rdname Group_specific_Var_AUC_estimation
#' @export 
#' @importFrom splines bs


Group_specific_Var_AUC_estimation <- function(MEM_Pol_group,time,Groups=NULL,method="trapezoid",Averaged=FALSE){
  '%notin%' <- Negate('%in%') 
  
  # Step 1: Verification of the type of arguments
  # ----- #
  Check_argument_Group_specific_Var_AUC(MEM_Pol_group,time,Groups,method,Averaged)
  
  Model_features <- MEM_Pol_group$Model_features
  Marginal_dynamics <- Model_features$Marginal.dyn.feature
  if(is.null(Groups)){
    Groups <- Model_features$Groups
  }
  if(is.numeric(time)){
    time <- lapply(seq(1,length(Groups)),function(g) return(time))
  }
  
  
  # Extraction of population parameters according to their groups
  if(is.list(MEM_Pol_group$Model_estimation)){
    Population_variance <- MEM_Pol_group$Model_estimation$varFix
  }else{
    Population_variance <- MEM_Pol_group$Model_estimation
  }
  MEM_groups <- as.vector(Model_features$Groups)
  
  global_intercept <- Marginal_dynamics$intercept["global.intercept"]
  ind_params <- 0
  Group_variance <- list()
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
    variance_matrix <- Population_variance[index_params,index_params]
    Group_variance[[MEM_groups[g]]] <- variance_matrix
  }
  
  
  # Step 2: Calculation of the variance of AUC
  # ----- #
  Estimated_var_AUC <- NULL
  for(g in 1:length(Groups)){
    time_group <- time[[g]]
    variance_mat_group <- Group_variance[[Groups[g]]]
    Pop_Covariate <- NULL
    
    if(global_intercept){
      Pop_Covariate <- cbind(Pop_Covariate,rep(1,length(time_group)))
    }
    
    # Extraction of information about model
    if(Marginal_dynamics$dynamic.type == "polynomial"){
      Covariate_poly_group <- do.call(cbind,lapply(1*isFALSE(Marginal_dynamics$intercept[paste("group.intercept",g,sep="")]):Marginal_dynamics$polynomial.degree[g],function(d) time_group^d))
      Pop_Covariate <- cbind(Pop_Covariate,Covariate_poly_group)
    }else if(Marginal_dynamics$dynamic.type == "spline"){
      # Creation of covariate matrix
      if(is.null(Marginal_dynamics$boundary.knots[[Groups[g]]])){
        Covariate_spline_group <- splines::bs(x=time_group,knots=Marginal_dynamics$knots[[Groups[g]]],df=Marginal_dynamics$df[g],degree=Marginal_dynamics$spline.degree[g])
      }else{
        Covariate_spline_group <- splines::bs(x=time_group,knots=Marginal_dynamics$knots[[Groups[g]]],df=Marginal_dynamics$df[g],degree=Marginal_dynamics$spline.degree[g],Boundary.knots=Marginal_dynamics$boundary.knots[[Groups[g]]])
      }
      if(Marginal_dynamics$intercept[paste("group.intercept",g,sep="")]){
        Covariate_spline_group <- cbind(rep(1,length(time_group)),Covariate_spline_group)
      }
      Pop_Covariate <- cbind(Pop_Covariate,Covariate_spline_group)
    }# End spline covariate
    
    # Estimation of the dynamics variance
    group_dyn_variance <- Pop_Covariate%*%variance_mat_group%*%t(Pop_Covariate)
    # Creation of method time weights (W) vector
    time_weights <- AUC_time_weights_estimation(time=time_group,method)
    Var_AUC_group <- time_weights %*% group_dyn_variance %*% time_weights
    Estimated_var_AUC <- c(Estimated_var_AUC,Var_AUC_group)
  }
  
  names(Estimated_var_AUC) <- Groups
  if(Averaged){
    Estimated_var_nAUC <- sapply(seq(1,length(Groups)),function(g) Estimated_var_AUC[g]/(diff(range(time[[g]]))^2),simplify=TRUE)
    Results <- Estimated_var_nAUC
  }else{
    Results <- Estimated_var_AUC
  }
  return(Results)
}
