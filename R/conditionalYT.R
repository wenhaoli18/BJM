#' conditional distribution of Y|T, if no competing risk; 
#' 
#' @description This function computes the conditional probability density function of 
#' longitudinal variable Y, given the survival time T without competing risk D.
#' 
#' @param data_predict_all This involves a collection of \code{data.frame} objects for 
#' dynamic prediction, each corresponding to a distinct longitudinal outcome. 
#' These data frames should contain the variables specified in \code{long_sub_fixed} 
#' and \code{long_sub_random}. Utilizing a list structure 
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
#' @param survival_variable_all The name of the transformed time-to-event outcomes variable.
#' @param survival_trans_function The transformation function used for time-to-event outcomes, 
#' in the order of \code{survival_variable_all}.
#' 
#' @return The output is a list containing probability matrices of \strong{log} densities
#' (so that they neither overflow nor underflow; see \code{predictRisk()} for how they
#' are exponentiated). In the presence of 
#' competing risks, this list includes two elements; otherwise, 
#' it contains only one element. Each element within the list is a probability matrix, 
#' with the number of rows (l_i) corresponding to specific time points and 
#' columns representing different patients. Every matrix element represents 
#' the conditional probability derived from the conditional distribution 
#' of longitudinal variable Y given the survival time T without competing risk D
#' for a particular patient at a specific time point.
#' @keywords internal
conditionalYT = function(data_predict_all, long_fit_all, l_i, survival_variable, 
                         time_variable, survival_variable_all, survival_trans_function){
  #LME model fitting
  lfit = long_fit_all$lfit
  #variance-covariance matrix
  Sigma = long_fit_all$Sigma_fit
  #patient ID
  num <- as.character(nlme::splitFormula(long_fit_all$long_sub_random[[1]], "|")[[2]])[2]
  
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
  data.long <- data_predict_all
  
  # Convert 'data.long' to a list if it is not a list
  if (!is.list(data.long) || is.data.frame(data.long)) {
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

    patient_data <- select_patient_longitudinal_data(data.long, num, num_i, n_longitudinal, time_variable)
    if(is.null(patient_data)) next
    rep_num_i_list <- patient_data$rep_num_i_list
    data_num_i_list <- patient_data$data_num_i_list

    design <- build_conditional_design(rep_num_i_list, data_num_i_list, lfit, Sigma,
                                        sigma.longitudinal, time_variable, n_longitudinal,
                                        long_fit_all$long_sub_random)
    longitudinal_all_matrix <- design$longitudinal_all_matrix
    parameter_matrix <- design$parameter_matrix
    log_det_Var_cov_estep <- design$log_det_Var_cov_estep
    Sigma_all_solve <- design$Sigma_all_solve

    ### for loop and make prediction probability for all time points in l_i
    for(it in 1: length(l_i)){
      
      ### get prediction data matrix for each biomarker for patient 'num_i'
      LME_indi_matrix = list()
      #interaction_orders = list()
      for(i in 1:n_longitudinal){
        model_formula = formula(lfit[[i]]) #lfit[[1]]
        ### terms cached from the fit on the *full* training data, not
        ### re-derived from this patient's own small slice -- see
        ### conditionalYDT.R for why this matters for poly()/
        ### splines::ns()/splines::bs()/factor() terms.
        terms_model <- lfit[[i]]$terms
        variable_names <- attr(terms_model, "term.labels")
        all_variables <- all.vars(model_formula)
        
        #survival variable replaced by l_i[it]
        data_num_i_list[[i]][survival_variable] = l_i[it]
        
        #transformed survival variable/basis function of survival variable
        #replaced by trans_function(l_i[it])
        if(length(survival_variable_all) != 0){
          for(surv_i in 1 : length(survival_variable_all)){
            data_num_i_list[[i]][survival_variable_all[[surv_i]]] = apply_survival_trans(survival_trans_function[[surv_i]], l_i[it], surv_i)
          }
        }

        target_covariate = survival_variable 
        ### if fuyrs exists or not
        if(target_covariate %in% all_variables != TRUE)
          stop("Error: Condition is false. Please add survival variable to linear mixed model.")
        else
          
          ### NA in nlme outcome (longitudinal biomarkers), replace with 999
          data_num_i_list[[i]][as.character(formula(long_fit_all$long_sub_fixed[[i]])[[2]])][is.na(data_num_i_list[[i]][as.character(formula(long_fit_all$long_sub_fixed[[i]])[[2]])])] <- 999
          xlev_i = if (!is.null(long_fit_all$xlevels)) long_fit_all$xlevels[[i]] else NULL
          mf_i = model.frame(terms_model, data_num_i_list[[i]], xlev = xlev_i)
          LME_indi_matrix[[i]] = t(model.matrix(terms_model, mf_i, contrasts.arg = lfit[[i]]$contrasts))
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
      
      ### Quadratic form (Y - mu)' Sigma_all_solve (Y - mu) for the FULL
      ### stacked (across all markers) observation vector.
      ###
      ### mean_all_matrix mirrors longitudinal_all_matrix's layout exactly:
      ### t(LME_all_matrix) %*% parameter_matrix is block-diagonal in the
      ### same way, so its column i holds mu_i = X_i %*% beta_i at marker
      ### i's own block of rows and zero elsewhere. rowSums() of each then
      ### collapses the M columns down to the single stacked Y/mu vector
      ### (safe because each row has a nonzero entry in exactly one column).
      ###
      ### NOTE: an earlier version of this function computed this quantity
      ### via a "trace trick" (sum(diag(...)) of an M x M matrix built from
      ### parameter_matrix/longitudinal_all_matrix sandwiched through
      ### Sigma_all_solve). That trace only sums the *diagonal* (i == j)
      ### blocks of the reduced M x M matrix, which silently discards every
      ### cross-marker (i != j) contribution of Sigma_all_solve -- i.e. it
      ### implicitly assumed the biomarkers were conditionally independent
      ### given the random effects. That contradicts the whole point of
      ### fitting a joint (possibly correlated) Sigma_fit across markers,
      ### and was numerically verified to diverge from the true joint
      ### quadratic form whenever markers are correlated (see NEWS.md).
      mean_all_matrix = t(LME_all_matrix) %*% parameter_matrix
      resid_full = rowSums(longitudinal_all_matrix) - rowSums(mean_all_matrix)
      quad_form = as.numeric(t(resid_full) %*% Sigma_all_solve %*% resid_full)

      ### log density (see predictRisk() for how it is exponentiated safely)
      f_Y_T_D_w1[it, iii] = -0.5 * log_det_Var_cov_estep - 0.5 * quad_form
    }
    
  }
  return(f_Y_T_D = list(f_Y_T_D_w1))
}
