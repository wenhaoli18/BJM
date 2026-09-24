#' conditional distribution of Y|T, if no competing risk; 
#' 
#' @description This function computes the conditional probability density function of 
#' longitudinal variable Y, given the survival time T without competing risk D.
#' 
#' @param data.predict.all This involves a collection of \code{data.frame} objects for 
#' dynamic prediction, each corresponding to a distinct longitudinal outcome. 
#' These data frames should contain the variables specified in \code{LongSubFixed} 
#' and \code{LongSubRandom}. Utilizing a list structure 
#' allows for the incorporation of multiple longitudinal outcomes, 
#' each potentially following different measurement protocols. 
#' In instances where all longitudinal outcomes are recorded at identical 
#' time points across patients, a singular \code{data.frame} object may 
#' be used in a \code{list}. It is presumed that each data frame is 
#' structured in a long format.
#' 
#' @param long_fit_all Outputs from the model fitting process using the \code{nlme} package, 
#' encompassing the results and parameters obtained from the analysis.
#' @param l_i A vector of time points to calculate the conditional probability.
#' @param survival_variable Time-to-event outcomes variable name.
#' @param time_variable The name of time variable in linear mixed model.
#' @param survivalVariableAll The name of the transformed time-to-event outcomes variable.
#' @param survivalTransFunction The transformation function used for time-to-event outcomes, 
#' in the order of \code{survivalVariableAll}.
#' 
#' @return The output is a list containing probability matrices. In the presence of 
#' competing risks, this list includes two elements; otherwise, 
#' it contains only one element. Each element within the list is a probability matrix, 
#' with the number of rows (l_i) corresponding to specific time points and 
#' columns representing different patients. Every matrix element represents 
#' the conditional probability derived from the conditional distribution 
#' of longitudinal variable Y given the survival time T without competing risk D
#' for a particular patient at a specific time point.
#' @keywords internal
conditionalYT = function(data.predict.all, long_fit_all, l_i, survival_variable, 
                         time_variable, survivalVariableAll, survivalTransFunction){
  #LME model fitting
  lfit = long_fit_all[[1]]
  #variance-covariance matrix
  Sigma = long_fit_all[[2]]
  #patient ID
  num <- as.character(nlme::splitFormula(long_fit_all[[4]][[1]], "|")[[2]])[2]
  
  #number of longitudinal biomarkers
  n_longitudinal <- length(lfit)  #length(data_num_i_list)
  #residual variance
  sigma.longitudinal = c()
  for(i in 1:n_longitudinal){
    sigma.longitudinal[i] = lfit[[i]]$sigma
  }

  # data.long is a list containing all your input data frames
  # data.long is an input list from user
  # data.long must be a list, it can contain n data frames and each element contains one biomarker
  # or data.long can be a list and only contain one data matrix, all biomarkers are contained
  data.long <- data.predict.all
  
  # Convert 'data.long' to a list if it is not a list
  if (!is.list(data.long)) {
    data.long <- list(data.long)
    data.long <- rep(data.long, each = n_longitudinal)
  }
  
  # A probability matrix, 
  # with the number of rows (l_i) corresponding to specific time points and 
  # columns representing different patients. 
  f_Y_T_D_w1 = matrix(NA, length(l_i), length(unlist(data.long[[1]][!duplicated(data.long[[1]][num]), ][num])))
  
  iii = 0
  for(num_i in unlist(data.long[[1]][!duplicated(data.long[[1]][num]), ][num])){
    iii = iii + 1
    #if(data.raw.sim.1[data.raw.sim.1[num] == num_i,]$status[1] == 1) next
    ### each rows represent intercept, slope, covariates numbers, time to event/l_i
    ### each columns represent different repeated measurements with different times.

    # Initialize lists to store the results for each subject 'num_i'
    # rep_num_i_list will store the repeated ones for each data frame.
    # data_num_i_list will store the filtered data for each patient num_i.
    rep_num_i_list <- list()
    data_num_i_list <- list()
    
    # Iterate over each data frame
    for (i in 1:n_longitudinal) {#
      df <- data.long[[i]]
      # Extract data for patient ID of 'num_i',  where 'num' equals 'num_i'
      selected_data <- df[df[num] == num_i, ]
      # Store the row length of patient ID of 'num_i', in a vector of repeated 1, 
      # Store in the list for different biomarkers
      rep_num_i_list[[i]] <- rep(1, length(unlist(selected_data[time_variable])))
      # Store the filtered patient ID of 'num_i' data with all variables in the list
      data_num_i_list[[i]] <- selected_data
    }
    #if all biomarkers contained in one data frame
    if(length(data.long) == 1){
      for (i in 1:n_longitudinal) {#
        rep_num_i_list[[i]] <- rep_num_i_list[[1]]
        data_num_i_list[[i]] <- data_num_i_list[[1]]
      }
    }
    
    # Use lapply to check the length of each element, 
    # and then use any to determine whether there is an element with a length of 0
    if(any(sapply(data_num_i_list, nrow) == 0)) next
    
    ####Initialize longitudinal matrix for all biomarkers
    longitudinal_all_matrix <- matrix(0, nrow = sum(sapply(data_num_i_list, nrow)), ncol = n_longitudinal) #n_data_num_i_list
    ####Constructing the longitudinal matrix for all biomarkers
    length_y = rep(0, n_longitudinal + 1)
    for (i in 1:n_longitudinal) {
      longname = as.character(formula(lfit[[i]]))[2]
      
      #length_y[i + 1] = length_y[i] + length(c(data_num_i_list[[i]][,longname])) #paste('longitudinal',i, sep = "")
      # change Mar 24 add unlist
      length_y[i + 1] = length_y[i] + length(c(unlist(data_num_i_list[[i]][,longname]))) #paste('longitudinal',i, sep = "")
      # change Mar 24 add unlist
      longitudinal_all_matrix[c((length_y[i] + 1) : length_y[i + 1]), i]  <- 
        unlist(data_num_i_list[[i]][,longname]) #paste('longitudinal',i, sep = "")
    }

    ####constructing the regression parameters' matrix
    n_lfit_total = 0 #total number of rows
    for(i in seq_len(length(lfit))){
      n_lfit_total = n_lfit_total + length(lfit[[i]]$coefficients$fixed)
    }
    parameter_matrix <- matrix(0, nrow = n_lfit_total, ncol = length(lfit))
    length_p = rep(0, n_longitudinal + 1)
    for (i in 1:n_longitudinal) {
      length_p[i + 1] = length_p[i] + length(lfit[[i]]$coefficients$fixed)
      parameter_matrix[c((length_p[i] + 1) : length_p[i + 1]), i]  <- lfit[[i]]$coefficients$fixed
    }

    A_i_ = list()
    ###random intercept or slope, depend on variance-covariance matrix
    if(dim(Sigma)[1] == n_longitudinal){
      ###random intercept
      for(i in 1:n_longitudinal){
        A_i_[[i]] = rbind(rep_num_i_list[[i]])
      }
    }else{
      ###random slope
      for(i in 1:n_longitudinal){
        A_i_[[i]] = rbind(rep_num_i_list[[i]], unlist(data_num_i_list[[i]][time_variable]))
      }
    }

    A_i <- matrix(0, nrow = sum(sapply(A_i_, ncol)), ncol = sum(sapply(A_i_, nrow)))
    length_A = rep(0, n_longitudinal + 1)
    for (i in 1:n_longitudinal) {
      length_A[i + 1] = length_A[i] + dim(A_i_[[i]])[2]
      if(dim(Sigma)[1] == n_longitudinal){
        ###random intercept
        A_i[c((length_A[i] + 1) : length_A[i + 1]), i]  <- t(A_i_[[i]])
      }else{
        ###random slope
        A_i[c((length_A[i] + 1) : length_A[i + 1]), (2*i-1):(2*i)]  <- t(A_i_[[i]])
      }
    }

    Sigma_vector = c()
    for(i in 1:n_longitudinal){
      Sigma_vector = c(Sigma_vector, rep(sigma.longitudinal[i]^2, dim(data_num_i_list[[i]])[1]))
    }
    Sigma_all =  A_i %*% Sigma %*% t(A_i) + diag(Sigma_vector)

    det_Var_cov_estep = det(2 * pi* Sigma_all)
    Sigma_all_solve = solve(Sigma_all)
    
    #A_matrix_1_1_loop = LME_all_matrix %*% Sigma_all_solve %*% t(LME_all_matrix)
    #A_matrix_2_1_loop = LME_all_matrix %*% Sigma_all_solve %*% longitudinal_all_matrix
    
    long_sigma_long = t(longitudinal_all_matrix) %*% Sigma_all_solve %*%  longitudinal_all_matrix
    
    ### for loop and make prediction probability for all time points in l_i
    for(it in 1: length(l_i)){
      
      ### get prediction data matrix for each biomarker for patient 'num_i'
      LME_indi_matrix = list()
      #interaction_orders = list()
      for(i in 1:n_longitudinal){
        model_formula = formula(lfit[[i]]) #lfit[[1]]
        terms_model <- terms(model_formula)
        variable_names <- attr(terms_model, "term.labels")
        all_variables <- all.vars(model_formula)
        
        #survival variable replaced by l_i[it]
        data_num_i_list[[i]][survival_variable] = l_i[it]
        
        #transformed survival variable/basis function of survival variable
        #replaced by trans_function(l_i[it])
        if(length(survivalVariableAll) != 0){
          for(surv_i in 1 : length(survivalVariableAll)){
            data_num_i_list[[i]][survivalVariableAll[[surv_i]]] = survivalTransFunction[[surv_i]](l_i[it])
          }
        }

        target_covariate = survival_variable 
        ### if fuyrs exists or not
        if(target_covariate %in% all_variables != TRUE)
          stop("Error: Condition is false. Please add survival variable to linear mixed model.")
        else
          
          ### NA in nlme outcome (longitudinal biomarkers), replace with 999
          data_num_i_list[[i]][as.character(formula(long_fit_all[[3]][[i]])[[2]])][is.na(data_num_i_list[[i]][as.character(formula(long_fit_all[[3]][[i]])[[2]])])] <- 999
          LME_indi_matrix[[i]] = t(model.matrix(long_fit_all[[3]][[i]], data_num_i_list[[i]]))
         if(dim( LME_indi_matrix[[i]] )[2] != dim(data_num_i_list[[i]])[1]){
           LME_indi_matrix[[i]] = cbind(LME_indi_matrix[[i]], matrix(NA, 
                dim(LME_indi_matrix[[i]] )[1], dim(data_num_i_list[[i]])[1] - dim( LME_indi_matrix[[i]] )[2]))
         }
           
      }
      
      # Initialize empty lists for rows and columns
      rows = list()
      columns = list()
      # Loop over the LME_indi_matrix
      # get prediction data matrix for all biomarkers for patient 'num_i'
      for (i in seq_len(length(LME_indi_matrix))) {
        for (j in seq_len(length(LME_indi_matrix))) {
          if (i == j) {
            # Add the matrix itself when row and column index are the same
            columns[[j]] = LME_indi_matrix[[i]]
          } else {
            # Add a zero matrix otherwise
            columns[[j]] = matrix(0, nrow=nrow(LME_indi_matrix[[i]]), ncol=ncol(LME_indi_matrix[[j]]))
          }
        }
        # Combine the columns for this row
        rows[[i]] = do.call(cbind, columns)
      }
      # Combine all the rows, combine all individual longitudinal matrix
      LME_all_matrix = do.call(rbind, rows)
      
      A_matrix_1_1_loop = LME_all_matrix %*% Sigma_all_solve %*% t(LME_all_matrix)
      A_matrix_2_1_loop = LME_all_matrix %*% Sigma_all_solve %*% longitudinal_all_matrix
      
      para_matrix_A_21 = t(parameter_matrix) %*% A_matrix_2_1_loop
      f_Y_T_D_w1[it, iii] = det_Var_cov_estep^{-0.5} * 
        exp(sum(diag(-0.5*( long_sigma_long + 
                              t(parameter_matrix) %*% A_matrix_1_1_loop %*% parameter_matrix - 
                              para_matrix_A_21 - 
                              t(para_matrix_A_21) ) ))) 
    }
    
  }
  return(f_Y_T_D = list(f_Y_T_D_w1))
}
