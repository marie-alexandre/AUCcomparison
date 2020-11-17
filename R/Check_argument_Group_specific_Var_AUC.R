#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param MEM_Pol_group PARAM_DESCRIPTION
#' @param time PARAM_DESCRIPTION
#' @param Groups PARAM_DESCRIPTION
#' @param method PARAM_DESCRIPTION
#' @param Averaged PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#'
#' @seealso 
#'  \code{\link[ArgumentCheck]{addError}}
#' @rdname Check_argument_Group_specific_Var_AUC
#' @importFrom ArgumentCheck newArgCheck addError finishArgCheck addMessage
#' @export
#' @keywords internal


Check_argument_Group_specific_Var_AUC <- function(MEM_Pol_group,time,Groups,method,Averaged){
  '%notin%' <- Negate('%in%') 
  # Verification of 'MEM_Pol_group'
  Check_MEM <- ArgumentCheck::newArgCheck()
  if(isFALSE(is.list(MEM_Pol_group))){
    ArgumentCheck::addError(
      msg = "'MEM_Pol_group' must be a list with similar structure than results provided by the function 'MEM_Polynomial_Group_structure' (only information about individual individual dynamics (random effects) can be missing)",
      argcheck = Check_MEM
    )
  }else{
    if(length(MEM_Pol_group)<2 || names(MEM_Pol_group) %notin% c("Model_estimation","Model_features")){
      ArgumentCheck::addError(
        msg = "The list 'MEM_Pol_group' must contain at least two lists named 'Model_estimation' and 'Model_features' as defined by the function 'MEM_Polynomial_Group_structure'",
        argcheck = Check_MEM
      )
    }else{
      # Verification of elements in 'Model_estimation'
      Check_Model_estimation <- ArgumentCheck::newArgCheck()
      if(isFALSE(is.list(MEM_Pol_group$Model_estimation) || is.matrix(MEM_Pol_group$Model_estimation))){
        ArgumentCheck::addError(
          msg = "'Model_estimation' in 'MEM_Pol_group' must be a list or a matrix of numerical values",
          argcheck = Check_Model_estimation
        )
      }else if(is.list(MEM_Pol_group$Model_estimation)){
        if("varFix" %notin% names(MEM_Pol_group$Model_estimation)){
          ArgumentCheck::addError(
            msg = "The list 'Model_estimation' in 'MEM_Pol_group' does not contain the covariance matrix of Fixed effects named 'varFix' ",
            argcheck = Check_Model_estimation
          )
        }else{
          if(isFALSE(is.matrix(MEM_Pol_group$Model_estimation$varFix) & is.numeric(MEM_Pol_group$Model_estimation$varFix))){
            ArgumentCheck::addError(
              msg = "The variable 'varFix' in 'Model_estimation' must be a matrix of numerical values",
              argcheck = Check_Model_estimation
            )
          }
        }
      }else{
        if(isFALSE(is.matrix(MEM_Pol_group$Model_estimation) & is.numeric(MEM_Pol_group$Model_estimation))){
          ArgumentCheck::addError(
            msg = "The variable 'Model_estimation' must contain matrix of numerical values corresponding to the covariance matrix of fixed effects",
            argcheck = Check_Model_estimation
          )
        }
      }
      ArgumentCheck::finishArgCheck(Check_Model_estimation)
      
      # Verification of elements in 'Model_features'
      Check_Model_features <- ArgumentCheck::newArgCheck()
      Model_features <- MEM_Pol_group$Model_features
      if(isFALSE(is.list(MEM_Pol_group$Model_features))){
        ArgumentCheck::addError(
          msg = "'Model_features' in 'MEM_Pol_group' must be a list as defined in by the function 'MEM_Polynomial_Group_structure'",
          argcheck = Check_Model_features
        )
      }else{
        # The list 'Model_features' must contain at least a vector 'Groups' and a list 'Marginal.dyn.feature'
        # Groups
        Check_Model_features_Groups <- ArgumentCheck::newArgCheck()
        if("Groups" %notin% names(Model_features)){
          ArgumentCheck::addError(
            msg = "'Model_features' does not contain the vector 'Groups'",
            argcheck = Check_Model_features_Groups
          )
        }
        ArgumentCheck::finishArgCheck(Check_Model_features_Groups)
        
        # # Marginal dynamic features
        Check_Model_features_dyn <- ArgumentCheck::newArgCheck()
        if("Marginal.dyn.feature" %notin% names(Model_features)){
          ArgumentCheck::addError(
            msg = "'Model_features' does not contain the list 'Marginal.dyn.feature' with information of marginal dynamics",
            argcheck = Check_Model_features_dyn
          )
        }else{
          Marginal_dynamics <- Model_features$Marginal.dyn.feature
          if(isFALSE(is.list(Marginal_dynamics))){
            ArgumentCheck::addError(
              msg = "'Marginal.dyn.feature' in 'Model_features' must be a list",
              argcheck = Check_Model_features_dyn
            )
          }else{
            # Dynamic type
            Check_Model_features_dyn_type <- ArgumentCheck::newArgCheck()
            if("dynamic.type" %notin% names(Marginal_dynamics)){
              ArgumentCheck::addError(
                msg = "'Marginal.dyn.feature' does not contain the argument 'dynamic.type'",
                argcheck = Check_Model_features_dyn_type
              )
            }else{
              if(isFALSE(is.character(Marginal_dynamics$dynamic.type)) || Marginal_dynamics$dynamic.type %notin% c("polynomial","spline")){
                ArgumentCheck::addError(
                  msg = "The variable 'dynamic.type' in 'Marginal.dyn.feature' must be defined as character type ('polynomial' or 'spline)",
                  argcheck = Check_Model_features_dyn_type
                )
              }
            }
            ArgumentCheck::finishArgCheck(Check_Model_features_dyn_type)
            
            # intercept
            Check_Model_features_dyn_intercept <- ArgumentCheck::newArgCheck()
            if("intercept" %notin% names(Marginal_dynamics)){
              ArgumentCheck::addError(
                msg = "'Marginal.dyn.feature' does not contain the argument 'intercept'",
                argcheck = Check_Model_features_dyn
              )
            }else{
              if(isFALSE(is.logical(Marginal_dynamics$intercept))){
                ArgumentCheck::addError(
                  msg = "The variable 'intercept' in 'Marginal.dyn.feature' must be a vector of booleans",
                  argcheck = Check_Model_features_dyn_intercept
                )
              }else if(length(Marginal_dynamics$intercept) != c(length(Model_features$Groups)+1)){
                ArgumentCheck::addError(
                  msg = paste("The variable 'intercept' in 'Marginal.dyn.feature' must contain ",length(Model_features$Groups)+1," booleans (global.intercept and group intercept for each group)",sep=""),
                  argcheck = Check_Model_features_dyn_intercept
                )
              }
            }
            ArgumentCheck::finishArgCheck(Check_Model_features_dyn_intercept)
            
            # Other arguments according to the type of dynamics
            if(Marginal_dynamics$dynamic.type == "polynomial"){
              # polynomial degree
              Check_Model_features_dyn_degree <- ArgumentCheck::newArgCheck()
              if("polynomial.degree" %notin% names(Marginal_dynamics)){
                ArgumentCheck::addError(
                  msg = "'Marginal.dyn.feature' does not contain the argument 'polynomial.degree'",
                  argcheck = Check_Model_features_dyn_degree
                )
              }else{
                if(isFALSE(is.numeric(Marginal_dynamics$polynomial.degree))){
                  ArgumentCheck::addError(
                    msg = "The variable 'polynomial.degree' in 'Marginal_dynamics' must be a vector of numeric values",
                    argcheck = Check_Model_features_dyn_degree
                  )
                }else if(length(Marginal_dynamics$polynomial.degree) != length(Model_features$Groups)){
                  ArgumentCheck::addError(
                    msg = paste("The variable 'polynomial.degree' in 'Marginal.dyn.feature' must contain ",length(Model_features$Groups)," values (one for each group involved in the model)",sep=""),
                    argcheck = Check_Model_features_dyn_degree
                  )
                }
              }
              ArgumentCheck::finishArgCheck(Check_Model_features_dyn_degree)
            }else{
              # spline degree
              Check_Model_features_dyn_degree <- ArgumentCheck::newArgCheck()
              if("spline.degree" %notin% names(Marginal_dynamics)){
                ArgumentCheck::addError(
                  msg = "'Marginal.dyn.feature' does not contain the argument 'spline.degree'",
                  argcheck = Check_Model_features_dyn_degree
                )
              }else{
                if(isFALSE(is.numeric(Marginal_dynamics$spline.degree))){
                  ArgumentCheck::addError(
                    msg = "The variable 'spline.degree' in 'Marginal_dynamics' must be a vector of numeric values",
                    argcheck = Check_Model_features_dyn_degree
                  )
                }else if(length(Marginal_dynamics$spline.degree) != length(Model_features$Groups)){
                  ArgumentCheck::addError(
                    msg = paste("The variable 'spline.degree' in 'Marginal.dyn.feature' must contain ",length(Model_features$Groups)," values (one for each group involved in the model)",sep=""),
                    argcheck = Check_Model_features_dyn_degree
                  )
                }
              }
              ArgumentCheck::finishArgCheck(Check_Model_features_dyn_degree)
              
              # knots
              Check_Model_features_dyn_knots <- ArgumentCheck::newArgCheck()
              if("knots" %notin% names(Marginal_dynamics)){
                ArgumentCheck::addError(
                  msg = "'Marginal.dyn.feature' does not contain the argument 'knots'",
                  argcheck = Check_Model_features_dyn_knots
                )
              }else{
                if(isFALSE(is.list(Marginal_dynamics$knots))){
                  ArgumentCheck::addError(
                    msg = paste("The variable 'knots' in 'Marginal.dyn.feature' must be a list of ",length(Model_features$Groups)," elements",sep=""),
                    argcheck = Check_Model_features_dyn_knots
                  )
                }else if(length(Marginal_dynamics$knots) != length(Model_features$Groups)){
                  ArgumentCheck::addError(
                    msg = paste("The variable 'knots' in 'Marginal.dyn.feature' must contain ",length(Model_features$Groups)," elements (one for each group involved in the model)",sep=""),
                    argcheck = Check_Model_features_dyn_degree
                  )
                }else{
                  for(g in 1:length(Marginal_dynamics$knots)){
                    if(isFALSE(is.null(Marginal_dynamics$knots[[g]])|| is.numeric(Marginal_dynamics$knots[[g]]))){
                      ArgumentCheck::addError(
                        msg = paste("The variable 'knots' in 'Marginal.dyn.feature' for Group ",Model_features$Groups[g]," must be NULL or a vector of numerical values",sep=""),
                        argcheck = Check_Model_features_dyn_degree
                      )
                    }
                  }
                }
              }
              ArgumentCheck::finishArgCheck(Check_Model_features_dyn_knots)
              
              # df
              Check_Model_features_dyn_df <- ArgumentCheck::newArgCheck()
              if("df" %notin% names(Marginal_dynamics)){
                ArgumentCheck::addError(
                  msg = "'Marginal.dyn.feature' does not contain the argument 'df'",
                  argcheck = Check_Model_features_dyn_df
                )
              }else{
                if(isFALSE(is.vector(Marginal_dynamics$df) || length(Marginal_dynamics$df) == length(Model_features$Groups))){
                  ArgumentCheck::addError(
                    msg = paste("The variable 'df' in 'Marginal.dyn.feature' must be a vector of ",length(Model_features$Groups),sep=""),
                    argcheck = Check_Model_features_dyn_df
                  )
                }else if(isFALSE(is.numeric(Marginal_dynamics$df))){
                  ArgumentCheck::addError(
                    msg = "The variable 'df' in 'Marginal.dyn.feature' must be a vector of numerical values",
                    argcheck = Check_Model_features_dyn_df
                  )
                }
              }
              ArgumentCheck::finishArgCheck(Check_Model_features_dyn_df)
              
              # Boundary knots
              Check_Model_features_dyn_boundary <- ArgumentCheck::newArgCheck()
              if("boundary.knots" %notin% names(Marginal_dynamics)){
                ArgumentCheck::addError(
                  msg = "'Marginal.dyn.feature' does not contain the argument 'boundary.knots'",
                  argcheck = Check_Model_features_dyn_boundary
                )
              }else{
                if(isFALSE(is.list(Marginal_dynamics$boundary.knots))){
                  ArgumentCheck::addError(
                    msg = "The variable 'boundary.knots' in 'Marginal.dyn.feature' must be a list",
                    argcheck = Check_Model_features_dyn_boundary
                  )
                }else if(length(Marginal_dynamics$boundary.knots) != length(Model_features$Groups)){
                  ArgumentCheck::addError(
                    msg = paste("The variable 'boundary.knots' in 'Marginal.dyn.feature' must contain ",length(Model_features$Groups)," elements",sep=""),
                    argcheck = Check_Model_features_dyn_boundary
                  )
                }else{
                  for(g in 1:length(Model_features$Groups)){
                    if(isFALSE(is.numeric(Marginal_dynamics$boundary.knots[[g]]) || is.null(Marginal_dynamics$boundary.knots[[g]]) || length(Marginal_dynamics$boundary.knots[[g]] != 2))){
                      ArgumentCheck::addError(
                        msg = paste("The variable 'boundary.knots' in 'Marginal.dyn.feature' in Group",Model_features$Groups[[g]]," must be either NULL or a vector of 2 numerical values",sep=""),
                        argcheck = Check_Model_features_dyn_boundary
                      )
                    }
                  }
                }
              }
              ArgumentCheck::finishArgCheck(Check_Model_features_dyn_boundary)
            }
          }
        }
        ArgumentCheck::finishArgCheck(Check_Model_features_dyn)
      }
      ArgumentCheck::finishArgCheck(Check_Model_features)
    }
  }
  ArgumentCheck::finishArgCheck(Check_MEM)
  
  # Verification of 'Group'
  Check_Group <- ArgumentCheck::newArgCheck()
  if(isFALSE(is.null(Groups) || is.vector(Groups))){
    ArgumentCheck::addError(
      msg = "'Groups' must be defined as NULL or a vector",
      argcheck = Check_Group
    )
  }else{
    if(is.null(Groups)){
      ArgumentCheck::addMessage(
        msg = "'Groups' being not defined, all Groups defined in 'MEM_Pol_group' are considered",
        argcheck = Check_Group
      )
    }else{
      # we verify that group names defined in Groups belongs are contained in MEM_Pol_group
      for(g in 1:length(Groups)){
        if(Groups[g] %notin% Model_features$Groups){
          ArgumentCheck::addError(
            msg = paste("The group name ",Groups[g]," defined in 'Groups' does not belong to Groups defined in 'MEM_Pol_group'",sep=""),
            argcheck = Check_Group
          )
        }
      }
    }
  }
  ArgumentCheck::finishArgCheck(Check_Group)
  if(is.null(Groups)){
    Groups <- Model_features$Groups
  }
  
  # Verification of 'time'
  Check_time <- ArgumentCheck::newArgCheck()
  if(isFALSE(is.numeric(time) || is.list(time))){
    ArgumentCheck::addError(
      msg = "'time' must be defined as a vector or a list of vectors (one for each group to estimate)",
      argcheck = Check_time
    )
  }else{
    if(is.numeric(time)){
      if(length(Groups) >1){
        ArgumentCheck::addMessage(
          msg = "The same 'time' vector has been used for all Groups",
          argcheck = Check_time
        )
      }
      if(length(time)<2){
        ArgumentCheck::addMessage(
          msg = "The vector 'time' must contain at least two values",
          argcheck = Check_time
        )
      }
    }else if(is.list(time)){
      if(length(time) != length(Groups)){
        ArgumentCheck::addError(
          msg = paste("The list 'time' must contain as much vector of time as the number of Group",sep=""),
          argcheck = Check_time
        )
      }else{
        for(g in 1:length(time)){
          if(isFALSE(is.numeric(time[[g]]))){
            ArgumentCheck::addError(
              msg = paste("The vector 'time' for Group", Groups[g],"must be a vector of numerical values",sep=""),
              argcheck = Check_time
            )
          }else if(length(time[[g]])<2){
            ArgumentCheck::addError(
              msg = paste("The vector 'time' for Group", Groups[g],"must contain at least two values",sep=""),
              argcheck = Check_time
            )
          }
        }
      }
    }
  }
  ArgumentCheck::finishArgCheck(Check_time)
  if(is.numeric(time)){
    time <- lapply(seq(1,length(Groups)),function(g) return(time))
  }
  
  # Verification of 'method'
  Check_method <- ArgumentCheck::newArgCheck()
  if(isFALSE(is.character(method))){
    ArgumentCheck::addError(
      msg = "The variable 'method' must be a string type",
      argcheck = Check_method
    )
  }else{
    if(method %notin% c("trapezoid","lagrange","spline")){
      ArgumentCheck::addError(
        msg = "The variable 'method' must take a value between 'trapezoid','lagrange' or 'spline'",
        argcheck = Check_method
      )
    }
  }  
  ArgumentCheck::finishArgCheck(Check_method)
  
  # Verification of 'Averaged'
  Check_Averaged <- ArgumentCheck::newArgCheck()
  if(isFALSE(is.logical(Averaged))){
    ArgumentCheck::addError(
      msg = "The variable 'Averaged' must be a boolean",
      argcheck = Check_Averaged
    )
  }
  ArgumentCheck::finishArgCheck(Check_Averaged)
  
}
