#' @title One-way analysis of variance of group specific AUC
#' @description This function performs a one-way ANOVA to compare the the area under the curves of multiple groups marginal dynamics, modeled by group-structured polynomials or B-spline curve in Mixed-Effects model.
#' Before performing the  ANOVA, this function can parform a Bartlett's test to evaluate homoscedasticity. In addition to ANOVA, users can decide to evaluate all the 2 by 2 comparisons.
#' @param MEM_Pol_group A list with similar structure than the output provided by the function \link[AUCcomparison]{MEM_Polynomial_Group_structure}. 
#' 
#' A list containing: \tabular{ll}{
#' \code{Model_estimation} \tab a list containing at least two elements: (1) the vector of the marginal (fixed) parameters estimates (at least for the groups whose AUC is to estimate), labeled _'beta'_ ;  (2) the variance-covariance matrix of these parameters, labeled _'varFix'_ (see \link[AUCcomparison]{MEM_Polynomial_Group_structure} for details about the parameter order). \cr
#' \code{Model_features} \tab a list of at least 2 elements: \cr
#' \tab 1. \code{Groups}  -  a vector indicating the names of the groups whose fixed parameters are given. \cr
#' \tab 2. \code{Marginal.dyn.feature}  -  a list summarizing the features of the marginal dynamics defined in the model:  \cr
#' \tab \itemize{
#' \item \code{dynamic.type} - a character scalar indicating the chosen type of marginal dynamics. Options are 'polynomial' or 'spline'
#' \item \code{intercept} - a logical vector summarizing choices about global and group-specific intercepts (Number of groups + 1) elements whose elements are named as ('global.intercept','group.intercept1', ..., 'group.interceptG') if G Groups are defined in \code{MEM_Pol_group}. For each element of the vector, if TRUE, the considered intercept is considered as included in the model (see \emph{Examples}). 
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

#' @param Groups a vector indicating the names of the groups belonging to the set of groups involved in \code{MEM_Pol_group} we want to include in the ANOVA (a subset or the entire set of groups involved in the model can be considered).
#' @param Time_groups a list of numerical vectors indicating the time points to consider in AUC calculation for each group (as much elements than the number of groups in \code{Groups}).
#' @param Nb_id_group a numerical vector indicating the number of individuals belonging to each group (as much elements than the number of groups in \code{Groups}).
#' @param common.interval a logical scalar. If TRUE (default) AUCs of all the compared groups are calculated on the same time interval defined as the intersect of all the time interval defined in \code{Time_groups}. If FALSE, AUC specific to each group is evaluated on its own interval.
#' @param method a character scalar indicating the interpolation method to use to estimate the AUC. Options are 'trapezoid' (default), 'lagrange' and 'spline'. In this version, the 'spline' interpolation is implemented with "not-a-knot" spline boundary conditions.
#' @param Averaged a logical scalar. If TRUE, AUC are evaluated as normalized AUC (nAUC) where nAUC is computated as the AUC divided by the range of time of calculation. If FALSE (default), the classic AUC is calculated (see \link[AUCcomparison]{Group_specific_AUC_estimation} for more details).
#' @param conf_level a numerical value (between 0 and 1) indicating the confidence level of the interval. By default, this variable is fixed at 0.95
#' @param bartlettTest PARAM_DESCRIPTION, Default: FALSE
#' @param data a dataframe gathering data for the groups involved in the ANOVA that have been fitted by the MEM model summarized in \code{MEM_Pol_group}. This dataframe has to contain at least 4 columns:
#' \itemize{
#'\item A column labelled 'Group' containing the information of the group for each observation (the same Groups than those defined in \code{Groups}).
#'\item A column labelled 'id' containing the information of individual identifier.
#'\item A column labelled 'time' containing information about the time of observations.
#'\item A column labelled 'value' containing the longitudinal observations.
#' }
#' @param twobytwo.comp a logical scalar indicating whether all the 2 by 2 comparisons must be evaluated after the ANOVA. Default: TRUE.
#' @param alternative a character scalar specifying the alternative hypothesis for the 2 by 2 comparisons. Options are 'two.sided' (default), 'greater' or 'less'. 
#'
#' @return A list of three elements elements:\tabular{ll}{
#' \code{bartlettTest} \tab a list of class "htest" corresponding to the Bartlett's test results (see \link[stats]{bartlett.test} for more details). If the test is not performed, a character 'Not performed' is returned. \cr
#' \code{ANOVA_F} \tab a list containg: \cr
#' \tab 1. \code{Followup} - a vector of two numerical values indicating the time interval when \code{Common.interval} = TRUE. If \code{common.interval}=FALSE, if variable is not included in the list.\cr
#' \tab 2. \code{Estimated_AUCs} - the vector of AUC estimated for each group involved in the ANOVA. \cr
#' \tab 3. \code{Estimated_VarAUCs} - the vector of the intra group variance estimated by the MEM (equiv. squared Standard error). \cr
#' \tab 4. \code{Between} - a vector gathering the between groups sum of squares, degree of freedom and variance with var=SS/df. \cr
#' \tab 5. \code{Within} - a vector gathering the within groups sum of squares, degree of dreedom and variance with var=SS/df. \cr
#' \tab 6. \code{ANOVA_F} - the value of the ANOVA statistic F. \cr
#' \tab 7. \code{Pvalue} - the value of the Pvalue. \cr
#' \tab 8. \code{Reject_H0} - a boolean indicating whether the null hypothesis (all AUC are equals) is rejected. \cr
#' \code{TwobyTwo_Comparison} \tab a list of K sublists where each sublist k gathers results of the kth comparison. If the 2 by 2 comparison is nor performed, a character 'Not performed' is returned. \cr
#' \tab 1. \code{Groups} - a vector indicating the 2 compared groups. \cr
#' \tab 2. \code{Estimated.AUC} - a vector of the 2 estimated values of AUC. \cr
#' \tab 3. \code{Delta_AUC} - the value of the difference of AUC between the 2 compared groups. \cr
#' \tab 4. \code{Tstat} - the value of the t-statistic \cr
#' \tab 5. \code{Pvalue} - the P-value (without any adjustment on multiple tests) \cr
#' \tab 6. \code{Conf.int} - the confidence interval \cr
#' \tab 7. \code{Adjusted.Pvalue} - the P-value adjusted on multiple testing by "bonferroni" method (see \link[stats]{p.adjust} for more details). \cr
#' }
#'  
#' @examples 
#' # Download of data
#' data("HIV_Simu_Dataset_Delta01_cens")
#' data <- HIV_Simu_Dataset_Delta01_cens
#' 
#' colnames(data) <- c("id","time","Group","value","cens") 
#' # Change factors in character vectors
#' data$id <- as.character(data$id) ; data$Group <- as.character(data$Group)
#' 
#' MEM_estimation <- MEM_Polynomial_Group_structure(y=data$value,x=data$time,Group=data$Group,Id=data$id,Cens=data$cens)
#' 
#' Groups <- unique(data$Group)
#' Time_groups <- lapply(seq(1,length(Groups)),function(g) return(unique(data$time[which(data$Group == Groups[g])])))
#' Nb_id_group <- sapply(seq(1,length(Groups)),function(g) return(length(unique(data$id[which(data$Group == Groups[g])]))))
#' 
#' ANOVA_test <- ANOVA_AUC_Group_Specific(MEM_Pol_group=MEM_estimation,Groups=Groups,Time_groups=Time_groups,
#'                                        Nb_id_group=Nb_id_group,common.interval=TRUE,method="trapezoid",Averaged=FALSE,conf_level=0.95,
#'                                        bartlettTest=TRUE,data=data,twobytwo.comp=TRUE,alternative="two.sided")
#'                                                    
#'                                                    
#' @seealso 
#'  \code{\link[ArgumentCheck]{addError}}
#'  \code{\link[dplyr]{reexports}}
#'  \code{\link[plyr]{ddply}}
#'  \code{\link[AUCcomparison]{AUC_time_weights_estimation}}
#'  \code{\link[stats]{bartlett.test}},\code{\link[stats]{FDist}}
#'  \code{\link[gtools]{combinations}}
#' @rdname ANOVA_AUC_Group_Specific
#' @export 
#' @importFrom ArgumentCheck newArgCheck addError finishArgCheck
#' @importFrom dplyr setdiff
#' @importFrom stats bartlett.test qf pf p.adjust
#' @importFrom gtools combinations
ANOVA_AUC_Group_Specific <- function(MEM_Pol_group,Groups,Time_groups,Nb_id_group,common.interval=TRUE,method="trapezoid",Averaged=FALSE,conf_level=0.95,
                                     bartlettTest=FALSE,data=NULL,twobytwo.comp=TRUE,alternative="two.sided"){
  
  '%notin%' <- Negate('%in%') 
  
  # Verification of bartlett 
  Check_bartlett <- ArgumentCheck::newArgCheck()
  if(isFALSE(is.logical(bartlettTest))){
    ArgumentCheck::addError(
      msg = "'bartlettTest' must be a boolean",
      argcheck = Check_bartlett
    )
  }else{
    # Verification of 'data' according to the value of bartlettTest
    if(isTRUE(bartlettTest)){
      Check_data <- ArgumentCheck::newArgCheck()
      if(is.null(data)){
        ArgumentCheck::addError(
          msg = "In order to realize the bartlett test, the data must be given in order to evaluate AUC at individual level",
          argcheck = Check_data
        )
      }else if(isFALSE(is.data.frame(data)) | length(dplyr::setdiff(c("Group","time","id","value"),names(data))) != 0){
        ArgumentCheck::addError(
          msg = "The variable 'data' must be a dataframe with at least 4 columns: 'Group','id','time' and 'value'",
          argcheck = Check_data
        )
      }#else{} Add the verification for the column of data
    }
    ArgumentCheck::finishArgCheck(Check_data)
  }
  ArgumentCheck::finishArgCheck(Check_bartlett)
  
  # Verification of Nb_id_group
  
  # Verification Groups
  
  # Verification Time_groups: a list of times (one element for each Group)
  
  # Verification twobytwo.comp
  Check_2by2comp <- ArgumentCheck::newArgCheck()
  if(isFALSE(is.logical(twobytwo.comp))){
    ArgumentCheck::addError(
      msg = "'twobytwo.comp' must be a boolean",
      argcheck = Check_2by2comp
    )
  }
  
  
  
  Multiple_comparison <- list()
  
  data <- data[which(!is.na(data$value)),]
  data <- data[order(data$time),]
  
  data.with.method <- data
  data.with.method$method <- method
  
  # Bartlett test 
  if(bartlettTest == TRUE){
    # Estimation of the individual AUC
    individuals <- unique(data$id)
    individual_AUC <- NULL
    for(id in 1:length(individuals)){
      data_id <- data[which(data$id == individuals[id]),]
      tmp_AUC <- as.numeric(data_id$value%*%AUC_time_weights_estimation(data_id$time,method=method))/diff(range(data_id$time))
      tmp_nAUC <- as.numeric(data_id$value%*%AUC_time_weights_estimation(data_id$time,method=method))
      individual_AUC <- rbind(individual_AUC,data.frame(id=individuals[id],Group=unique(data_id$Group),AUC=tmp_AUC,nAUC=tmp_nAUC,stringsAsFactors = FALSE))
    }
    # individual_AUC <- plyr::ddply(.data=data.with.method,.variables = .(id,Group),here(summarise),
    #                               AUC=as.numeric(value%*%AUC_time_weights_estimation(time,method=unique(method)))/diff(range(time)),
    #                               nAUC=as.numeric(value%*%AUC_time_weights_estimation(time,method=unique(method))))
    if(Averaged){
      bartlett_test <- stats::bartlett.test(nAUC~Group,data=individual_AUC)
    }else{
      bartlett_test <- stats::bartlett.test(AUC~Group,data=individual_AUC)
    }
  }else{
    bartlett_test <- "Not performed"
  }
  Multiple_comparison$bartlettTest <- bartlett_test
  
  ## Analyse of variance
  if(common.interval){
    Min_Timefollowp <- sapply(seq(1,length(Time_groups)),function(g) return(min(Time_groups[[g]],na.rm=TRUE)))
    Max_Timefollowp <- sapply(seq(1,length(Time_groups)),function(g) return(max(Time_groups[[g]],na.rm=TRUE)))
    min.time.interval <- max(Min_Timefollowp,na.rm=TRUE)
    max.time.interval <- min(Max_Timefollowp,na.rm=TRUE)
    
    if(min.time.interval > max.time.interval){
      stop("Impossible to estimate the difference of AUC on the common interval of time for all the groups: the inrestection of the time intervals is equal to an empty interval.")
    }else if(min.time.interval == max.time.interval){
      stop(paste("Impossible to estimate the difference of AUC on the common interval of time for the two groups: the inrestection of the time intervals is equal to a single time point:",min.time.interval))
    }else{
      Time_groups <- lapply(seq(1,length(Time_groups)),function(g,minT,maxT){
        time_group <- Time_groups[[g]]
        time_group <- time_group[which(time_group >= minT & time_group <= maxT)]
        return(time_group)
      },minT=min.time.interval,maxT=max.time.interval)
    }
  }
  

  # Estimation of AUC and Var AUC for all the groups 
  Estimated_AUCs <- Group_specific_AUC_estimation(MEM_Pol_group=MEM_Pol_group,time=Time_groups,Groups=Groups,method=method,Averaged=Averaged)
  Estimated_VarAUC <- Group_specific_Var_AUC_estimation(MEM_Pol_group=MEM_Pol_group,time=Time_groups,Groups=Groups,method=method,Averaged=Averaged)
  
  # Estimation of the inter-group variability (SS group) ####
  Mean_AUCs <- mean(Estimated_AUCs,na.rm=TRUE)  # Mean of the estimated restricted nAUCs 
  SSgroup <- sum((Estimated_AUCs-Mean_AUCs)^2*Nb_id_group,na.rm=TRUE)
  df_group <- length(Groups) - 1
  S2group <- SSgroup/df_group
  
  # Step 4. Estimation of the intra-group variability (SS individual) ####
  SSind <- sum(Estimated_VarAUC*(Nb_id_group)^2,na.rm = TRUE)
  df_ind <- sum(Nb_id_group) - length(Groups)
  S2ind <- SSind/df_ind
  
  # Step 5. Estimation of the ANOVA statistic
  Anova_stat_F <- S2group/S2ind
  F_critic <- stats::qf(conf_level,df1=df_group,df2=df_ind)
  P_value <- round(stats::pf(abs(Anova_stat_F),df1=df_group,df2=df_ind,lower.tail=FALSE),digits = 5)
  Reject_H0 <- Anova_stat_F > F_critic 
  
  if(common.interval){
    Results_ANOVA <- list(Followup=c(min.time.interval,max.time.interval),
                          Estimated_AUCs=Estimated_AUCs,Estimated_VarAUCs=Estimated_VarAUC,
                          Between=c(SSbetween=SSgroup,dfbetween=df_group,Varbetween=S2group),
                          Within=c(SSwithin=SSind,dfwithin=df_ind,Varwithin=S2ind),
                          ANOVA_F=Anova_stat_F,Pvalue=P_value,Reject_H0=Reject_H0)
  }else{
    Results_ANOVA <- list(Estimated_AUCs=Estimated_AUCs,Estimated_VarAUCs=Estimated_VarAUC,
                          Between=c(SSbetween=SSgroup,dfbetween=df_group,Varbetween=S2group),
                          Within=c(SSwithin=SSind,dfwithin=df_ind,Varwithin=S2ind),
                          ANOVA_F=Anova_stat_F,Pvalue=P_value,Reject_H0=Reject_H0)
  }
  Multiple_comparison$ANOVA <- Results_ANOVA
  
  
  ## 2 by 2 comparisons
  if(twobytwo.comp){
    TwobyTwo_Comparison <- list()
    Pvalues_2by2 <- NULL
    if(!Reject_H0){
      warning(paste("As requested, the 2 by 2 comparisons are performed but the ANOVA concluded to the null hypothesis with pvalue of",P_value))
    }
    Comparisons <- data.frame(gtools::combinations(n=length(Groups),r=2,repeats.allowed = FALSE,v=Groups),stringsAsFactors = FALSE)
    names(Comparisons) <- c("GroupA","GroupB")
    
    for(comp in 1:nrow(Comparisons)){
      GroupA <- Comparisons$GroupA[comp]
      time_point_groupA <- Time_groups[[which(Groups == GroupA)]]
      
      GroupB <- Comparisons$GroupB[comp]
      time_point_groupB <- Time_groups[[which(Groups == GroupB)]]
      
      tstat <- Stat_test_Delta_AUC_Group_Specific(MEM_Pol_group=MEM_Pol_group,Group1=GroupA,Group2=GroupB,time.G1=time_point_groupA,time.G2=time_point_groupB,
                                                  common.interval=common.interval,method=method,Group.dependence=TRUE,Averaged=Averaged,conf_level=conf_level,alternative=alternative)
      
      Pvalues_2by2 <- c(Pvalues_2by2,tstat$Pvalue)
      TwobyTwo_Comparison[[comp]] <- list(Groups=c(GroupA,GroupB),Estimated.AUC=tstat$AUCs,Delta_AUC=tstat$Delta_AUC,Tstat=tstat$Tstat,Pvalue=tstat$Pvalue,
                                          Conf.int=tstat$Conf.int) 
    }
    Adjusted.Pvalues <- round(stats::p.adjust(p=Pvalues_2by2,method = "bonferroni"),digits=5)
    for(comp in 1:nrow(Comparisons)){
      TwobyTwo_Comparison[[comp]]$Adjusted.Pvalue <- Adjusted.Pvalues[comp]
    }
    Multiple_comparison$TwobyTwo_Comparison <- TwobyTwo_Comparison
  }else{
    Multiple_comparison$TwobyTwo_Comparison <- "Not performed"
  }
  return(Multiple_comparison)
}


