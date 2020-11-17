MEM_Polynomial_Group_structure <- function(y,x=NULL,Group=NULL,Id=NULL,Cens=NULL,
                                           marginal_dyn_type="polynomial",ind_dyn_type="polynomial",
                                           global_intercept=TRUE,group_intercept=FALSE,degree_group=3,
                                           Adaptive="none",min_knots_group=2,max_knots_group=2,
                                           knots_group = NULL,df_group = NULL,Boundary.knots_group=NULL,
                                           ind_intercept=FALSE,degree_ind=3,min_knots_ind=2,max_knots_ind=2,
                                           same_base_group_ind=FALSE,knots_ind=NULL,df_ind=NULL,Boundary.knots_ind=NULL,...){
  
  '%notin%' <- Negate('%in%') 
  
  # Step 1: Creation of the dataframe gathering y, x, cens and Group: #####
  # ------ #
  Check_y <- ArgumentCheck::newArgCheck()
  if(isFALSE(is.data.frame(y) || is.numeric(y))){
    ArgumentCheck::addError(
      msg = "'y' must be a dataframe or a vector of numerical values",
      argcheck = Check_y
    )
  }else if(is.data.frame(y)){
    
    # Verification of 'y'
    Check_dataframe_y <- ArgumentCheck::newArgCheck()
    if("y" %notin% colnames(y)){
      ArgumentCheck::addError(
        msg = "The variable 'y' must be specified in the dataframe 'y'.",
        argcheck = Check_dataframe_y
      )
    }
    ArgumentCheck::finishArgCheck(Check_dataframe_y)
    
    # Verification about 'x'
    Check_dataframe_x <- ArgumentCheck::newArgCheck()
    if("x" %notin% colnames(y)){
      if(is.null(x)){
        ArgumentCheck::addError(
          msg = "The variable 'x' must be provided either in y dataframe or x vector. Be sure to use the right variable name in 'y'.",
          argcheck = Check_y
        )
      }else{
        if(isFALSE(is.numeric(x) || is.integer(x))){
          ArgumentCheck::addError(
            msg = "The variable 'x' must be a vector of numerical values",
            argcheck = Check_x
          )
        }else{
          if(nrow(y) != length(x)){
            ArgumentCheck::addError(
              msg = paste("The variables 'y' and 'x' must have the same size:",
                          nrow(y),"rows in y and",length(x),"elements in x",sep=" "),
              argcheck = Check_x
            )
          }else{
            y <- cbind(y,x=x)
          }
        }
      }
    }else{
      if(isFALSE(is.numeric(y$x) || is.integer(y$x))){
        ArgumentCheck::addError(
          msg = "The variable 'x' provided in the dataframe 'y' must be a vector of numerics or integers.",
          argcheck = Check_x
        )
      }
      if(!is.null(x)){
        ArgumentCheck::addWarning(
          msg = "The variable 'x' has been provided twice (in 'y' and as argument 'x'). Values given in 'y' are used. Be sure to use the right values.",
          argcheck = Check_x
        )
      }
    }
    ArgumentCheck::finishArgCheck(Check_dataframe_x)
    
    # Verification about 'Group'
    Check_dataframe_group <- ArgumentCheck::newArgCheck()
    if("Group" %notin% colnames(y)){
      if(is.null(Group)){
        ArgumentCheck::addMessage(
          msg = "The variable 'Group' has not been specified in the dataframe y or as independent argument 'Group'. Consequently a single group is considered.",
          argcheck = Check_y
        )
        Group <- rep("Group1",nrow(y))
        y <- cbind(y,Group=Group)
      }else{
        if(isFALSE(is.numeric(Group) || is.character(Group))){
          ArgumentCheck::addError(
            msg = "The variable 'Group' must be a vector of numerical values or a vector of characters.",
            argcheck = Check_group
          )
        }else{
          if(nrow(y) != length(Group)){
            ArgumentCheck::addError(
              msg = paste("The variables 'y' and 'Group' must have the same size:",
                          nrow(y),"rows in y and",length(Group),"elements in Group",sep=" "),
              argcheck = Check_group
            )
          }else{
            y <- cbind(y,Group=Group)
          }
        }
      }
    }else{
      if(isFALSE(is.numeric(y$Group) || is.integer(y$Group) || is.character(y$Group))){
        ArgumentCheck::addError(
          msg = "The variable 'Group' provided in the dataframe 'y' must be a vector of numerics or integers or characters.",
          argcheck = Check_Group
        )
      }
      if(isFALSE(is.null(Group))){
        ArgumentCheck::addWarning(
          msg = "The variable 'Group' has been provided twice (in 'y' and as argument 'Group'). Values given in 'y' are used. Be sure to use the right values.",
          argcheck = Check_group
        )
      }
    }
    ArgumentCheck::finishArgCheck(Check_dataframe_group)
    
    # Verification about 'Id'
    Check_dataframe_Id <- ArgumentCheck::newArgCheck()
    if("Id" %notin% colnames(y)){
      if(is.null(Id)){
        ArgumentCheck::addError(
          msg = "The variable 'Id' must be provided either in y dataframe or Id vector. Be sure to use the right variable name in 'y'.",
          argcheck = Check_y
        )
      }else{
        if(isFALSE(is.numeric(Id) || is.character(Id))){
          ArgumentCheck::addError(
            msg = "The variable 'Id' must be a vector of numerical values or a vector of characters.",
            argcheck = Check_Id
          )
        }else{
          if(nrow(y) != length(Id)){
            ArgumentCheck::addError(
              msg = paste("The variable 'y' and 'Id' must have the same size:",
                          nrow(y),"rows in y and",length(Id),"elements in Id",sep=" "),
              argcheck = Check_Id
            )
          }else{
            y <- cbind(y,Id=Id)
          }
        }
      }
    }else{
      if(isFALSE(is.numeric(y$Id) || is.integer(y$Id) || is.character(y$Id))){
        ArgumentCheck::addError(
          msg = "The variable 'Id' provided in the dataframe 'y' must be a vector of numerics, integers or characters.",
          argcheck = Check_Id
        )
      }
      if(!is.null(Id)){
        ArgumentCheck::addWarning(
          msg = "The variable 'Id' has been provided twice (in 'y' and as argument 'Id'). Values given in 'y' are used. Be sure to use the right values.",
          argcheck = Check_Id
        )
      }
    }
    ArgumentCheck::finishArgCheck(Check_dataframe_Id)
    
    # Verification about 'Cens'
    Check_dataframe_Cens <- ArgumentCheck::newArgCheck()
    if("Cens" %notin% colnames(y)){
      if(is.null(Cens)){
        ArgumentCheck::addMessage(
          msg = "The variable 'Cens' has not been specified. Data will be treated as uncensored",
          argcheck = Check_y
        )
        y <- cbind(y,Cens=rep(0,nrow(y)))
      }else{
        if(isFALSE(is.numeric(Cens))){
          ArgumentCheck::addError(
            msg = "The variable 'Cens' must be a vector of numerical values",
            argcheck = Check_cens
          )
        }else{
          if(nrow(y) != length(Cens)){
            ArgumentCheck::addError(
              msg = paste("The variable 'y' and 'Cens' must have the same size:",
                          nrow(y),"rows in y and",length(Cens),"elements in Cens",sep=" "),
              argcheck = Check_cens
            )
          }else{
            y <- cbind(y,Cens=Cens)
          }
        }
      }
    }else{
      if(isFALSE(is.numeric(y$Cens) || is.integer(y$Cens))){
        ArgumentCheck::addError(
          msg = "The variable 'Cens' provided in the dataframe 'y' must be a vector of numerics or integers.",
          argcheck = Check_cens
        )
      }
      if(!is.null(Cens)){
        ArgumentCheck::addWarning(
          msg = "The variable 'Cens' has been provided twice (in 'y' and as arguement 'Cens'). Values given in 'y' are used. Be sure to use the right values.",
          argcheck = Check_cens
        )
      }
    }
    ArgumentCheck::finishArgCheck(Check_dataframe_Cens)
    
    data <- y[,c("Id","x","Group","y","Cens")]
  }else{
    # Verification about 'x'
    Check_x <- ArgumentCheck::newArgCheck()
    if(is.null(x)){
      ArgumentCheck::addError(
        msg = "The variable 'x' must be specified",
        argcheck = Check_x
      )
    }else{
      if(isFALSE(is.numeric(x) || is.integer(x))){
        ArgumentCheck::addError(
          msg = "The variable 'x' must be a vector of numerical values",
          argcheck = Check_x
        )
      }else{
        if(length(y) != length(x)){
          ArgumentCheck::addError(
            msg = paste("The variables 'y' and 'x' must have the same size:",
                        length(y),"elements in 'y' and",length(x),"elements in 'x'.",sep=" "),
            argcheck = Check_x
          )
        }
      }
    }
    ArgumentCheck::finishArgCheck(Check_x)
    
    # Verification about 'Group'
    Check_group <- ArgumentCheck::newArgCheck()
    if(is.null(Group)){
      ArgumentCheck::addMessage(
        msg = "The variable 'Group' has not been specified. Consequently, a single group is considered",
        argcheck = Check_group
      )
      Group <- rep("Group1",length(y))
    }else{
      if(isFALSE(is.numeric(Group) || is.character(Group))){
        ArgumentCheck::addError(
          msg = "The variable 'Group' must be a vector of numerical values or a vector of characters.",
          argcheck = Check_group
        )
      }else{
        if(length(y) != length(Group)){
          ArgumentCheck::addError(
            msg = paste("The variables 'y' and 'Group' must have the same size:",
                        lenght(y),"rows in y and",length(Group),"elements in Group",sep=" "),
            argcheck = Check_group
          )
        }
      }
    }
    ArgumentCheck::finishArgCheck(Check_group)
    
    # Verification of 'Id'
    Check_Id <- ArgumentCheck::newArgCheck()
    if(is.null(Id)){
      ArgumentCheck::addError(
        msg = "The variable 'Id' must be specified.",
        argcheck = Check_Id
      )
    }else{
      if(isFALSE(is.numeric(Id) || is.character(Id) || is.integer(Id))){
        ArgumentCheck::addError(
          msg = "The variable 'Id' must be a vector of numerical values or a vector of characters.",
          argcheck = Check_Id
        )
      }else{
        if(length(y) != length(Id)){
          ArgumentCheck::addError(
            msg = paste("The variable 'y' and 'Id' must have the same size:",
                        length(y),"rows in y and",length(Id),"elements in Id",sep=" "),
            argcheck = Check_Id
          )
        }
      }
    }
    ArgumentCheck::finishArgCheck(Check_Id)
    
    # Verification of 'Cens'
    Check_cens <- ArgumentCheck::newArgCheck()
    if(is.null(Cens)){
      ArgumentCheck::addMessage(
        msg = "The variable 'Cens'has not been specified. Data will be treated as uncensored.",
        argcheck = Check_cens
      )
      Cens <- rep(0,length(y))
    }else{
      if(isFALSE(is.numeric(Cens))){
        ArgumentCheck::addError(
          msg = "The variable 'Cens' must be a vector of numerical values",
          argcheck = Check_cens
        )
      }else{
        if(length(y) != length(Cens)){
          ArgumentCheck::addError(
            msg = paste("The variable 'y' and 'Cens' must have the same size:",
                        length(y),"rows in y and",length(Cens),"elements in Cens",sep=" "),
            argcheck = Check_cens
          )
        }
      }
    }
    ArgumentCheck::finishArgCheck(Check_cens)
    
    data <- data.frame(Id=Id,x=x,Group=Group,y=y,Cens=Cens)
  } 
  ArgumentCheck::finishArgCheck(Check_y)
  
  data <- data[with(data,order(Group,Id,x)),]
  rownames(data) <- seq(1,nrow(data))
  
  
  # Estimation of the number of Group based on data
  Groups <- as.vector(unique(data$Group))
  Nb_groups <- length(Groups)
  
  
  # Step 2: Choice of the type of model according to given arguments ####
  # ----- #
  
  # 1. Marginal dynamics 
  Check_marg_typ <- ArgumentCheck::newArgCheck()
  check_degree_group <- ArgumentCheck::newArgCheck()
  check_global_intercept <- ArgumentCheck::newArgCheck()
  check_group_intercept <- ArgumentCheck::newArgCheck()
  
  # Verification of 'marginal_dyn_type'
  if(isFALSE(is.character(marginal_dyn_type))){
    ArgumentCheck::addError(
      msg = "The variable 'marginal_dyn_type' must be a string type",
      argcheck = Check_marg_typ
    )
  }else{
    if(marginal_dyn_type %notin% c("polynomial","spline")){
      ArgumentCheck::addError(
        msg = "The variable 'marginal_dyn_type' has been wrongly chosen. Please choose between 'polynomial' or 'spline'.",
        argcheck = Check_marg_typ
      )
    }
  }
  ArgumentCheck::finishArgCheck(Check_marg_typ)
  
  # Verification of 'degree_group'
  if(isFALSE(is.numeric(degree_group))){
    ArgumentCheck::addError(
      msg = "The variable 'degree_group' must be a numerical value",
      argcheck = check_degree_group
    )
  }else{
    if(degree_group<=0){
      ArgumentCheck::addError(
        msg = "The variable 'degree_group' must be >0",
        argcheck = check_degree_group
      )
    }
  }
  ArgumentCheck::finishArgCheck(check_degree_group)
  
  # Verification of 'global_intercept'
  if(isFALSE(is.logical(global_intercept))){
    ArgumentCheck::addError(
      msg = "The variable 'global_intecept' must be a boolean",
      argcheck = check_global_intercept
    )
  }
  ArgumentCheck::finishArgCheck(check_degree_group)
  
  # Verification of 'group_intercept'
  if(isFALSE(is.logical(group_intercept))){
    ArgumentCheck::addError(
      msg = "The variable 'group_intercept' must be a boolean or a vector of boolean",
      argcheck = check_group_intercept
    )
  }else{
    if(length(group_intercept) %notin% c(1,Nb_groups)){
      ArgumentCheck::addError(
        msg = paste("The variable 'group_intercept' must have a length of 1 or the number of group defined in the data (",Nb_groups,")",sep=""),
        argcheck = check_group_intercept
      )
    }
  }
  ArgumentCheck::finishArgCheck(check_group_intercept)
  
  if(marginal_dyn_type == "spline"){
    # Verification of 'Adaptive"
    Check_Adaptive <- ArgumentCheck::newArgCheck()
    if(isFALSE(is.character(Adaptive))){
      ArgumentCheck::addError(
        msg = "The variable 'Adaptive' must be a character.",
        argcheck = Check_Adaptive
      )
    }else{
      if(length(Adaptive) != 1){
        ArgumentCheck::addError(
          msg ="The variable 'Adaptive' must have a length of 1",
          argcheck = Check_Adaptive
        )
      }else{
        if(Adaptive %notin% c("none","group","individual","both")){
          ArgumentCheck::addError(
            msg = "The variable 'Adaptive' must take a value between 'none', 'group', 'individual' or 'both'.",
            argcheck = Check_Adaptive
          )
        }
      }
    }
    ArgumentCheck::finishArgCheck(Check_Adaptive)
    
    if(Adaptive %in% c("group","both")){
      # Verification of 'min_knots_group' and 'max_knots_group'
      Check_min_knots_group <- ArgumentCheck::newArgCheck()
      Check_max_knots_group <- ArgumentCheck::newArgCheck()
      
      if(isFALSE(is.numeric(min_knots_group))){
        ArgumentCheck::addError(
          msg = "The variable 'min_knots_group' must be an integer.",
          argcheck = Check_min_knots_group
        )
      }else{
        min_knots_group <- as.integer(min_knots_group)
        if(min_knots_group<=0){
          ArgumentCheck::addError(
            msg = "The variable 'min_knots_group' must be positive.",
            argcheck = Check_min_knots_group
          )
        }
      }
      if(isFALSE(is.numeric(max_knots_group))){
        ArgumentCheck::addError(
          msg = "The variable 'max_knots_group' must be an integer.",
          argcheck = Check_max_knots_group
        )
      }else{
        max_knots_group <- as.integer(max_knots_group)
        if(max_knots_group < min_knots_group){
          ArgumentCheck::addError(
            msg = "The variable 'max_knots_group' must be bigger than'min_knots_group'.",
            argcheck = Check_min_knots_group
          )
        }
      }
      ArgumentCheck::finishArgCheck(Check_min_knots_group)
      ArgumentCheck::finishArgCheck(Check_max_knots_group)
    }else{
      # Verification of 'knots_group', 'df_group' and 'Boundary.knots_group'
      Check_knot_group <- ArgumentCheck::newArgCheck()
      Check_df_group <- ArgumentCheck::newArgCheck()
      Check_Bound.knots_group <- ArgumentCheck::newArgCheck()
      
      if(isFALSE(is.list(knots_group) || is.numeric(knots_group) || is.null(knots_group))){
        ArgumentCheck::addError(
          msg = "The variable 'knots_group' must either be NULL, a numeric vector or a list.",
          argcheck = Check_knot_group
        )
      }else if(is.list(knots_group)){
        if(length(knots_group) %notin% c(1,Nb_groups)){
          ArgumentCheck::addError(
            msg = paste("The variable 'knots_group' must be a list of length 1 or equal to the number of groups (",Nb_groups,")",sep=""),
            argcheck = Check_knot_group
          )
        }else{
          for(g in 1:Nb_groups){
            if(isFALSE(is.numeric(knots_group[[g]]))){
              ArgumentCheck::addError(
                msg = "The variable 'knots_group' must be a list of vectors of numerical values, for all groups",
                argcheck = Check_knot_group
              )
            }
          }
        }
      }
      
      if(isFALSE(is.null(df_group) || is.numeric(df_group))){
        ArgumentCheck::addError(
          msg = "The variable 'df_group' must either be NULL or a vector of numerical values.",
          argcheck = Check_df_group
        )
      }else if(isFALSE(is.null(df_group))){
        if(length(df_group) %notin% c(1,Nb_groups)){
          ArgumentCheck::addError(
            msg = paste("The variable 'df_group' must be a vector of length 1 or equal to the number of groups (",Nb_groups,")",sep=""),
            argcheck = Check_df_group
          )
        }
      }
      
      if(isFALSE(is.null(Boundary.knots_group) || is.numeric(Boundary.knots_group) || is.list(Boundary.knots_group))){
        ArgumentCheck::addError(
          msg = "The variable 'Boundary.knots_group' must either be NULL or a list.",
          argcheck = Check_Bound.knots_group
        )
      }else if(is.list(Boundary.knots_group)){
        if(length(Boundary.knots_group) %notin% c(1,Nb_groups)){
          ArgumentCheck::addError(
            msg = paste("The variable 'Boundary.knots_group' must be a list of length 1 or equal to the number of groups (",Nb_groups,")",sep=""),
            argcheck = Check_Bound.knots_group
          )
        }else{
          for(g in 1:Nb_groups){
            if(isFALSE(is.numeric(Boundary.knots_group[[g]]))){
              ArgumentCheck::addError(
                msg = "The variable 'Boundary.knots_group' must be a list of vectors of numerical values, for all groups",
                argcheck = Check_Bound.knots_group
              )
            }
          }
        }
      }
      ArgumentCheck::finishArgCheck(Check_knot_group)
      ArgumentCheck::finishArgCheck(Check_df_group)
      ArgumentCheck::finishArgCheck(Check_Bound.knots_group)
    }
  }
  
  
  # 2. Individual dynamics
  Check_ind_typ <- ArgumentCheck::newArgCheck()
  check_degree_ind <- ArgumentCheck::newArgCheck()
  check_ind_intercept <- ArgumentCheck::newArgCheck()
  
  # Verification of 'ind_dyn_type'
  if(isFALSE(is.character(ind_dyn_type))){
    ArgumentCheck::addError(
      msg = "The variable 'ind_dyn_type' must be a string type",
      argcheck = Check_ind_typ
    )
  }else{
    if(ind_dyn_type %notin% c("polynomial","spline")){
      ArgumentCheck::addError(
        msg = "The variable 'ind_dyn_type' has been wrongly chosen. Please choose between 'polynomial' or 'spline'.",
        argcheck = Check_ind_typ
      )
    }
  }
  ArgumentCheck::finishArgCheck(Check_ind_typ)
  
  # Verification of 'degree_ind'
  if(isFALSE(is.numeric(degree_ind))){
    ArgumentCheck::addError(
      msg = "The variable 'degree_ind' must be a numerical value",
      argcheck = check_degree_ind
    )
  }else{
    if(degree_ind<=0){
      ArgumentCheck::addError(
        msg = "The variable 'degree_ind' must be >0",
        argcheck = check_degree_ind
      )
    }
  }
  ArgumentCheck::finishArgCheck(check_degree_ind)
  
  # Verification of 'ind_intercept'
  if(isFALSE(is.logical(ind_intercept))){
    ArgumentCheck::addError(
      msg = "The variable 'ind_intercept' must be a boolean or a vector of boolean",
      argcheck = check_ind_intercept
    )
  }else{
    if(length(ind_intercept) != 1){
      ArgumentCheck::addError(
        msg = "The variable 'ind_intercept' must have a length of 1",
        argcheck = check_ind_intercept
      )
    }
  }
  ArgumentCheck::finishArgCheck(check_ind_intercept)
  
  if(ind_dyn_type == "spline"){
    # Verification of 'Adaptive"
    Check_Adaptive <- ArgumentCheck::newArgCheck()
    if(isFALSE(is.character(Adaptive))){
      ArgumentCheck::addError(
        msg = "The variable 'Adaptive' must be a character.",
        argcheck = Check_Adaptive
      )
    }else{
      if(length(Adaptive) != 1){
        ArgumentCheck::addError(
          msg ="The variable 'Adaptive' must have a length of 1",
          argcheck = Check_Adaptive
        )
      }else{
        if(Adaptive %notin% c("none","group","individual","both")){
          ArgumentCheck::addError(
            msg = "The variable 'Adaptive' must take a value between 'none', 'group', 'individual' or 'both'.",
            argcheck = Check_Adaptive
          )
        }
      }
    }
    ArgumentCheck::finishArgCheck(Check_Adaptive)
    
    if(Adaptive %in% c("individual","both")){
      # Verification of 'min_knots_ind' and 'max_knots_ind'
      Check_min_knots_ind <- ArgumentCheck::newArgCheck()
      Check_max_knots_ind <- ArgumentCheck::newArgCheck()
      
      if(isFALSE(is.numeric(min_knots_ind))){
        ArgumentCheck::addError(
          msg = "The variable 'min_knots_ind' must be an integer.",
          argcheck = Check_min_knots_ind
        )
      }else{
        min_knots_group <- as.integer(min_knots_ind)
        if(min_knots_ind<=0){
          ArgumentCheck::addError(
            msg = "The variable 'min_knots_ind' must be positive.",
            argcheck = Check_min_knots_ind
          )
        }
      }
      if(isFALSE(is.numeric(max_knots_ind))){
        ArgumentCheck::addError(
          msg = "The variable 'max_knots_ind' must be an integer.",
          argcheck = Check_max_knots_ind
        )
      }else{
        max_knots_group <- as.integer(max_knots_ind)
        if(max_knots_ind < min_knots_ind){
          ArgumentCheck::addError(
            msg = "The variable 'max_knots_ind' must be bigger than 'min_knots_ind'.",
            argcheck = Check_max_knots_ind
          )
        }
      }
      ArgumentCheck::finishArgCheck(Check_min_knots_ind)
      ArgumentCheck::finishArgCheck(Check_max_knots_ind)
      
      if(Adaptive == "both"){
        # Verification of 'same_base_group_ind'
        Check_same_base_group_ind <- ArgumentCheck::newArgCheck()
        if(isFALSE(is.logical(same_base_group_ind))){
          ArgumentCheck::addError(
            msg = "The variable 'same_base_group_ind' must be bigger a boolean",
            argcheck = Check_same_base_group_ind
          )
        }
      }
    }else{
      # Verification of 'knots_ind', 'df_ind' and 'Boundary.knots_ind'
      Check_knot_ind <- ArgumentCheck::newArgCheck()
      Check_df_ind <- ArgumentCheck::newArgCheck()
      Check_Bound.knots_ind <- ArgumentCheck::newArgCheck()
      
      if(isFALSE(is.vector(knots_ind) || is.null(knots_ind) || is.list(knots_ind))){
        ArgumentCheck::addError(
          msg = "The variable 'knots_ind' must either be NULL, a vector or a list.",
          argcheck = Check_knot_ind
        )
      }else if(is.list(knots_ind)){
        if(length(knots_ind) %notin% c(Nb_groups,length(unique(data$Id)))){
          ArgumentCheck::addError(
            msg = "The variable 'knots_ind' as list must contain as much elements as the number of groups or individuals.",
            argcheck = Check_knot_ind
          )
        }
      }
      
      if(isFALSE(is.null(df_ind) || is.numeric(df_ind))){
        ArgumentCheck::addError(
          msg = "The variable 'df_ind' must either be NULL or numerical values.",
          argcheck = Check_df_ind
        )
      }else if(isFALSE(is.null(df_ind))){
        if(length(df_ind) %notin% c(1,Nb_groups,length(unique(data$Id)))){
          ArgumentCheck::addError(
            msg = paste("The variable 'df_ind' must be a single value (same df for all individuals)",
                        "\n","or a vector containing as much elements as the number of groups or individuals",sep=""),
            argcheck = Check_df_ind
          )
        }
      }
      
      if(isFALSE(is.null(Boundary.knots_ind) || is.vector(Boundary.knots_ind))){
        ArgumentCheck::addError(
          msg = "The variable 'Boundary.knots_ind' must either be NULL or a vector of 2 elements",
          argcheck = Check_Bound.knots_ind
        )
      }else if(is.vector(Boundary.knots_ind)){
        if(length(Boundary.knots_ind) != 2){
          ArgumentCheck::addError(
            msg = "The variable 'Boundary.knots_ind' must be a vector of 2 elements",
            argcheck = Check_Bound.knots_ind
          )
        }
      }
      
      ArgumentCheck::finishArgCheck(Check_knot_ind)
      ArgumentCheck::finishArgCheck(Check_df_ind)
      ArgumentCheck::finishArgCheck(Check_Bound.knots_ind)
    }
  }
  
  
  # Step 3: Build of the design matrix of the marginal dynamics ####
  # ----- #
  Pop_Covariate <- NULL
  Bsplines_groups <- list()
  
  if(marginal_dyn_type == "polynomial"){
    if(global_intercept){
      Pop_Covariate <- cbind(Pop_Covariate,intercept=rep(1,nrow(data)))
    }
    if(length(degree_group) == 1){
      degree_group <- rep(degree_group,Nb_groups)
    }
    if(length(group_intercept) == 1){
      group_intercept <- rep(group_intercept,Nb_groups)
    }
    for(g in 1:Nb_groups){
      data_group <- data[which(data$Group == Groups[g]),]
      tmp_pop_covariate <- do.call(cbind,lapply(1*isFALSE(group_intercept[g]):degree_group[g],function(d){
        res <- rbind(matrix(0,nrow=nrow(data[which(data$Group != Groups[g] & as.numeric(rownames(data)) < as.numeric(rownames(data_group))[1]),]),ncol=1),
                     matrix(data_group$x^d,ncol=1),
                     matrix(0,nrow=nrow(data[which(data$Group != Groups[g] & as.numeric(rownames(data)) > as.numeric(rownames(data_group))[nrow(data_group)]),]),ncol=1))
        return(res)
      }))
      colnames(tmp_pop_covariate) <- paste("X",1*isFALSE(group_intercept[g]):degree_group[g],".G",g,sep="")
      Pop_Covariate <- cbind(Pop_Covariate,tmp_pop_covariate)
    }
  }else if(marginal_dyn_type == "spline"){
    if(global_intercept){
      Pop_Covariate <- cbind(Pop_Covariate,intercept=rep(1,nrow(data)))
    }
    
    if(length(degree_group) == 1){
      degree_group <- rep(degree_group,Nb_groups)
    }
    if(length(group_intercept) == 1){
      group_intercept <- rep(group_intercept,Nb_groups)
    }
    if(length(Boundary.knots_group) == 1 || is.null(Boundary.knots_group)){
      Boundary.knots_group <- setNames(lapply(seq(1,Nb_groups),function(g) return(Boundary.knots_group)),Groups)
    }
    if(length(df_group) == 1 || is.null(df_group)){
      df_group <- setNames(lapply(seq(1,Nb_groups),function(g) return(df_group)),Groups)
    }
    
    if(Adaptive %in% c("group","both")){
      # Research of optimal knots for each group
      knots_group <- setNames(lapply(seq(1,Nb_groups), function(g) Optimal_knot_research(data=data[which(data$Group == Groups[g]),],degree=degree_group[g],minknot=min_knots_group,maxknot=max_knots_group)),Groups)
    }else if(is.null(knots_group) || isTRUE(is.list(knots_group) & length(knots_group) == 1 )){
      knots_group <- setNames(lapply(seq(1,Nb_groups), function(g) knots_group),Groups)
    }
    
    for(g in 1:Nb_groups){
      data_group <- data[which(data$Group == Groups[g]),]
      if(is.null(Boundary.knots_group[[g]])){
        Covariate_spline_group <- splines::bs(x=data_group$x,knots=knots_group[[g]],df=df_group[[g]],degree=degree_group[g])
      }else{
        Covariate_spline_group <- splines::bs(x=data_group$x,knots=Adaptive_Knot_groups[[g]],df=df_group[[g]],degree=degree_group[g],
                                              Boundary.knots=Boundary.knots_group[[g]])
      }
      tmp_pop_covariate <- rbind(matrix(0,nrow=nrow(data[which(data$Group != Groups[g] & as.numeric(rownames(data)) < as.numeric(rownames(data_group))[1]),]),ncol=ncol(Covariate_spline_group)),
                                 Covariate_spline_group,
                                 matrix(0,nrow=nrow(data[which(data$Group != Groups[g] & as.numeric(rownames(data)) > as.numeric(rownames(data_group))[nrow(data_group)]),]),ncol=ncol(Covariate_spline_group)))
      colnames(tmp_pop_covariate) <- paste("Spl",seq(1,ncol(Covariate_spline_group)),".G",g,sep="")
      if(group_intercept[g]){
        add_group_intercept <- c(rep(0,nrow(data[which(data$Group != Groups[g] & as.numeric(rownames(data)) < as.numeric(rownames(data_group))[1]),])),
                                 rep(1,nrow(Covariate_spline_group)),
                                 rep(0,nrow(data[which(data$Group != Groups[g] & as.numeric(rownames(data)) > as.numeric(rownames(data_group))[nrow(data_group)]),])))
        tmp_pop_covariate <- cbind(add_group_intercept,tmp_pop_covariate)
        colnames(tmp_pop_covariate)[1] <- paste("intercept.G",g,sep="")
      }
      Pop_Covariate <- cbind(Pop_Covariate,tmp_pop_covariate)
      Bsplines_groups[[Groups[g]]] <- Covariate_spline_group
    }
  }
  
  
  # Step 4: Build of the design matrix of the random effects ####
  # ----- #
  Rand_covariate <-  NULL 
  IDs <- unique(data$Id)
  
  if(ind_intercept){
    Rand_covariate <- cbind(Rand_covariate,intercept=rep(1,nrow(data)))
  }
  
  if(ind_dyn_type == "polynomial"){
    tmp_rnd_covariate <- NULL
    for(i in 1:length(IDs)){
      data_pat <- data[which(data$Id == IDs[i]),]
      tmp_rnd_covariate <- rbind(tmp_rnd_covariate,do.call(cbind,lapply(1:degree_ind,function(d) data_pat$x^d)))
      colnames(tmp_rnd_covariate) <- paste("Z",1:degree_ind,sep="")
    }
    Rand_covariate <- cbind(Rand_covariate,tmp_rnd_covariate)
  }else if(ind_dyn_type == "spline"){
    if(Adaptive %in% c("individual","both")){
      # Research of optimal knots for each individual
      knots_ind <- setNames(lapply(seq(1,length(IDs)),function(i) Optimal_knot_research(data=data[which(data$Id == IDs[i]),],degree=degree_ind,minknot=min_knots_ind,maxknot=max_knots_ind)),IDs)
    }else if(length(knots_ind)==1 || is.null(knots_ind)){
      knots_ind <- setNames(lapply(seq(1,length(IDs)),function(i) return(knots_ind)),IDs)
    }else if(length(knots_ind) == Nb_groups){
      knots_ind <- setNames(lapply(seq(1,length(IDs)),function(i) return(knots_ind[[unique(data$Group[which(data$Id == IDs[i])])]])),IDs) 
    }
    tmp_rnd_covariate <- NULL
    for(i in 1:length(IDs)){
      data_pat <- data[which(data$Id == IDs[i]),]
      if(same_base_group_ind & marginal_dyn_type == "spline"){
        covariate_spline_pat <- predict(Bsplines_groups[[unique(data_pat$Group)]],newx=data_pat$x)
      }else{
        if(is.null(Boundary.knots_ind)){
          covariate_spline_pat <- splines::bs(data_pat$x,df=df_ind,knots=knots_ind[[i]],degree=degree_ind)
        }else{
          covariate_spline_pat <- splines::bs(data_pat$x,df=df_ind,knots=knots_ind[[i]],degree=degree_ind,Boundary.knots=Boundary.knots_ind)
        }
      }
      tmp_rnd_covariate <- rbind(tmp_rnd_covariate,covariate_spline_pat) 
      colnames(tmp_rnd_covariate) <- paste("Ind.Spl",seq(1,ncol(covariate_spline_pat)),sep="")
    }
    Rand_covariate <- cbind(Rand_covariate,tmp_rnd_covariate)
  }
  
  
  # Step 5:  Model Estimation ####
  # ----- #
  cluster <- unlist(lapply(seq(1,length(IDs)),function(x) rep(x,nrow(data[which(data$Id == IDs[x]),]))))
  model <- lmec::lmec(yL=data$y,cens=data$Cens,X=Pop_Covariate,Z=Rand_covariate,cluster=cluster,...)
  
  # OUTPUT: ####
  # ----- #
  marginal.dyn.feature <- list(dynamic.type=marginal_dyn_type,
                               intercept=c(global.intercept=global_intercept,group.intercept=group_intercept))
  names(degree_group) <- Groups
  if(marginal_dyn_type == "polynomial"){
    marginal.dyn.feature[["polynomial.degree"]] <- degree_group
  }else{
    marginal.dyn.feature[["spline.degree"]] <- degree_group
    marginal.dyn.feature[["adaptive.splines"]] <- Adaptive %in% c("group","both")
    marginal.dyn.feature[["knots"]] <- setNames(lapply(seq(1,Nb_groups),function(g){
      if(is.null(knots_group[[g]])){
        if(is.null(df_group[[g]]) || df_group[[g]] <= degree_group[g]){
          return(knots_group[[g]])
        }else{
          return(attr(Bsplines_groups[[Groups[g]]],"knots"))
        }
      }else{
        return(knots_group[[g]])
      } 
    }),Groups)
    marginal.dyn.feature[["df"]] <- sapply(seq(1,Nb_groups),function(g){
      res <- NULL
      if(is.null(df_group[[g]])){
        if(is.null(knots_group[[g]])){
          res <- degree_group[g]
        }else{
          res <- length(knots_group[[g]])+degree_group[g]
        }
      }else if(df_group[[g]] > degree_group[g]){
        res <- df_group[[g]]
      }else{
        res <- degree_group[g]
      }
      names(res) <- Groups[g]
      return(res)
    })
    marginal.dyn.feature[["boundary.knots"]] <- setNames(lapply(seq(1,Nb_groups),function(g){
      res <- Boundary.knots_group[[g]]
      if(is.null(res)){
        res <- range(data$x[which(data$Group == Groups[g])])
      }
      return(res)
    }),Groups)
  }
  
  individual.dyn.feature <- list(dynamic.type=ind_dyn_type,
                                 intercept=ind_intercept)
  if(ind_dyn_type == "polynomial"){
    individual.dyn.feature[["polynomial.degree"]] <- degree_ind
  }else{
    if(same_base_group_ind & marginal_dyn_type == "spline"){
      individual.dyn.feature[["splines"]] <- "Same spline basis than defined in group for each individual"
    }else{
      individual.dyn.feature[["spline.degree"]] <- degree_ind
      individual.dyn.feature[["adaptive.splines"]] <-  Adaptive %in% c("individual","both")
      if(Adaptive %in% c("individual","both")){
        individual.dyn.feature[["knots"]] <- do.call("rbind",knots_ind)
      }else{
        if(is.null(unique(knots_ind)[[1]])){
          individual.dyn.feature[["knots"]] <- "defined by bs function"
        }else if(length(unique(knots_ind)) == Nb_groups){
          individual.dyn.feature[["knots"]] <- setNames(unique(knots_ind),Groups)
        }else{
          individual.dyn.feature[["knots"]] <- unique(knots_ind)[[1]]
        }
      }
      if(is.null(df_ind)){
        individual.dyn.feature[["df"]] <- "defined by bs function"
      }else{
        individual.dyn.feature[["df"]] <- df_ind
      }
      if(is.null(Boundary.knots_ind)){
        individual.dyn.feature[["boundary.knots"]] <- "defined by bs function"
      }else{
        individual.dyn.feature[["boundary.knots"]] <- Boundary.knots_ind
      }
    }
  }
  
  Results <- list(Model_estimation=model,
                  Model_features=list(Groups=Groups,
                                      Marginal.dyn.feature = marginal.dyn.feature,
                                      Individual.dyn.feature = individual.dyn.feature))
  return(Results)
}
