#' @title Area Under The Curve of Group-Specific Polynomial Marginal Dynamics
#' @description This function estimates the area under the curve of marginal dynamics modeled by group-structured polynomials or B-spline curves.
#' 
#' @param MEM_Pol_group A list with similar structure than the output provided by the function \link[AUCcomparison]{MEM_Polynomial_Group_structure}. 
#' 
#' A list containing: \tabular{ll}{
#' \code{Model_estimation} \tab either the vector of the marginal parameter estimates (at least for the groups whose AUC is to estimate), or a list containing this vector labeled _'beta'_ among other arguments (see \link[AUCcomparison]{MEM_Polynomial_Group_structure} for details about the parameter order). \cr
#' \code{Model_features} \tab a list of at least 2 elements: \cr
#' \tab 1. \code{Groups}  -  a vector indicating the names of the groups whose fixed parameters are given. \cr
#' \tab 2. \code{Marginal.dyn.feature}  -  a list summarizing the features of the marginal dynamics defined in the model:  \cr
#' \tab \itemize{
#' \item \code{dynamic.type} - a character scalar indicating the chosen type of marginal dynamics. Options are 'polynomial' or 'spline'
#' \item \code{intercept} -  a logical vector summarizing choices about global and group-specific intercepts (Number of groups + 1) elements whose elements are named as ('global.intercept','group.intercept1', ..., 'group.interceptG') if G Groups are defined in \code{MEM_Pol_group}. For each element of the vector, if TRUE, the considered intercept is considered as included in the model (see \emph{Examples}). 
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
#' @param Groups a vector indicating the names of the groups belonging to the set of groups involved in \code{MEM_Pol_group} for which we want to estimate the AUC  (a subset or the entire set of groups involved in the model can be considered). If NULL (default), the AUC for all the groups involved the MEM is calculated.
#' @param method a character scalar indicating the interpolation method to use to estimate the AUC. Options are 'trapezoid' (default), 'lagrange' and 'spline'. In this version, the 'spline' interpolation is implemented with the "not-a-knot" spline boundary conditions. 
#' @param Averaged a logical scalar. If TRUE, the function return the normalized AUC (nAUC) computed as the AUC divided by the range of the time calculation. If FALSE (default), the classic AUC is calculated.
#' @return A numerical vector containing the estimation of the AUC (or nAUC) for each group defined in the \code{Groups} vector.

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
#' 
#' # Estimation of the AUC for the two groups defined in the dataset
#' AUC_estimation <- AUC_estimation <- Group_specific_AUC_estimation(MEM_Pol_group=MEM_estimation,time=list(unique(data$time[which(data$Group == "Group1")]),unique(data$time[which(data$Group == "Group2")])),Groups=unique(data$Group))
#' 
#' Estimation of the AUC only for the group "Group1"
#' AUC_estimation_G1 <- Group_specific_AUC_estimation(MEM_Pol_group=MEM_estimation,time=unique(data$time[which(data$Group == "Group1")]),Groups=c("Group1"))
#'
#' Example 2: We consider results of MEM estimation from another source. We have to give build the variable 'MEM_Pol_group' with the good structure
#'  # We build the variable 'MEM_Pol_group.1' with the results of MEM estimation obtained for two groups (even if only "Group1" is called in AUC estimation function)
#' MEM_Pol_group.1 <- list(Model_estimation=c(1.077,0.858,-0.061,0.0013,0.887,-0.066,0.0014), # c(global.intercept,beta1.G1,beta2.G1,beta2.G1,beta1.G2,beta2.G2,beta3.G2)
#'                         Model_features=list(Groups=c("Group1","Group2"),
#'                                             Marginal.dyn.feature=list(dynamic.type="polynomial",intercept=c(global.intercept=TRUE,group.intercept1=FALSE,group.intercept2=FALSE),polynomial.degree=c(3,3))))
#' 
#'  # We build the variable 'MEM_Pol_group.2' with the results of MEM estimation obtained only for the group of interest (extraction)
#' MEM_Pol_group.2 <- list(Model_estimation=c(1.077,0.858,-0.061,0.0013), # c(global.intercept,beta1.G1,beta2.G1,beta3.G1)
#'                         Model_features=list(Groups=c("Group1"),
#'                                             Marginal.dyn.feature=list(dynamic.type="polynomial",intercept=c(global.intercept=TRUE,group.intercept1=FALSE),polynomial.degree=c(3))))
#'                       
#'# Estimation of the AUC for the group "Group1"
#' AUC_estimation_G1.1 <- Group_specific_AUC_estimation(MEM_Pol_group=MEM_Pol_group.1,time=unique(data$time[which(data$Group == "Group1")]),Groups=c("Group1"))
#' AUC_estimation_G1.2 <- Group_specific_AUC_estimation(MEM_Pol_group=MEM_Pol_group.2,time=unique(data$time[which(data$Group == "Group1")]))



#' @seealso 
#'  \code{\link[splines]{bs}}, 
#'  \code{\link[AUCcomparison]{MEM_Polynomial_Group_structure}}
#' @rdname Group_specific_AUC_estimation
#' @export 
#' @importFrom splines bs
#' 

Group_specific_AUC_estimation <- function(MEM_Pol_group,time,Groups=NULL,method="trapezoid",Averaged=FALSE){
  '%notin%' <- Negate('%in%') 
  
  # Step 1: Verification of the type of arguments
  # ----- #
  Check_argument_Group_specific_AUC(MEM_Pol_group,time,Groups,method,Averaged)
  
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
    Population_params <- MEM_Pol_group$Model_estimation$beta
  }else{
    Population_params <- MEM_Pol_group$Model_estimation
  }
  MEM_groups <- as.vector(Model_features$Groups)
  
  global_intercept <- Marginal_dynamics$intercept["global.intercept"]
  ind_params <- 0
  Group_parameters <- list()
  for(g in 1:length(MEM_groups)){
    params <- NULL
    if(global_intercept){
      params <- c(params,Population_params[1])
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
    params <- c(params,Population_params[(ind_params+1):(ind_params+Nb_group_params)])
    ind_params <- ind_params + Nb_group_params
    Group_parameters[[MEM_groups[g]]] <- params
  }
  
  # Step 2: Calculation of AUC
  # ----- #
  Estimated_AUC <- NULL
  for(g in 1:length(Groups)){
    time_group <- time[[g]]
    beta_group <- Group_parameters[[Groups[g]]]
    Pop_Covariate <- NULL
    
    if(global_intercept){
      Pop_Covariate <- cbind(Pop_Covariate,rep(1,length(time_group)))
    }
    
    # Extraction of information about model
    if(Marginal_dynamics$dynamic.type == "polynomial"){
      # Creation of covariate matrix
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
    
    # Estimation of the marginal dynamics
    Group_dynamics <- as.numeric(Pop_Covariate %*% beta_group)
    # Creation of method time weights (W) vector
    time_weights <- AUC_time_weights_estimation(time=time_group,method)
    AUC_group <- as.numeric(Group_dynamics %*% time_weights)
    Estimated_AUC <- c(Estimated_AUC,AUC_group)
  }
  
  names(Estimated_AUC) <- Groups
  if(Averaged){
    Estimated_nAUC <- sapply(seq(1,length(Groups)),function(g) Estimated_AUC[g]/diff(range(time[[g]])),simplify=TRUE)
    Results <- Estimated_nAUC
  }else{
    Results <- Estimated_AUC
  }
  return(Results)
}
