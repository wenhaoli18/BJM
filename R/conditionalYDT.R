#' conditional distribution of Y|D, T, if with competing risk
#' 
#' @description This function computes the conditional probability density function of 
#' longitudinal variable Y, given the survival time T with competing risk D.
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
#' @return The output is a list containing probability matrices. In the presence of 
#' competing risks, this list includes two elements; otherwise, 
#' it contains only one element. Each element within the list is a probability matrix, 
#' with the number of rows (l_i) corresponding to specific time points and 
#' columns representing different patients. Every matrix element represents 
#' the conditional probability derived from the conditional distribution 
#' of longitudinal variable Y given the survival time T with competing risk D
#' for a particular patient at a specific time point.
#' @keywords internal
conditionalYDT = function(data_predict_all, long_fit_all, survival_fit_all, 
                          l_i, survival_variable, 
                          time_variable, survival_variable_all, survival_trans_function){
  #LME model fitting
  lfit = long_fit_all$lfit
  #variance-covariance matrix
  Sigma = long_fit_all$Sigma_fit
  #patient ID
  num <- as.character(nlme::splitFormula(long_fit_all$long_sub_random[[1]], "|")[[2]])[2]
  ### event type variable name
  event_type_variable = as.character(formula(survival_fit_all$form_conditional_cr)[[2]])
  
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
  f_Y_T_D_w0 = matrix(NA, length(l_i), length(unlist(data.long[[1]][!duplicated(data.long[[1]][num]), ][num])))

  ### Each biomarker's longitudinal model formula (and the variables/outcome
  ### name derived from it) is fixed for the whole function call -- it does
  ### not depend on num_i or l_i[it] at all -- so it is derived once here,
  ### rather than being recomputed inside the patient/l_i loops below (where
  ### it previously ran length(l_i) times per patient per biomarker).
  model_formula_all = list()
  all_variables_all = list()
  outcome_var_all = list()
  for(i in 1:n_longitudinal){
    model_formula_all[[i]] = formula(lfit[[i]])
    all_variables_all[[i]] = all.vars(model_formula_all[[i]])
    outcome_var_all[[i]] = as.character(model_formula_all[[i]][[2]])
  }

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
                                        sigma.longitudinal, time_variable, n_longitudinal)
    longitudinal_all_matrix <- design$longitudinal_all_matrix
    parameter_matrix <- design$parameter_matrix
    det_Var_cov_estep <- design$det_Var_cov_estep
    Sigma_all_solve <- design$Sigma_all_solve
    long_sigma_long <- design$long_sigma_long

    ### Build a reusable model.frame "template" per biomarker/event-type,
    ### once per patient, instead of inside the l_i loop below. Only the
    ### survival_variable/survival_variable_all cells change across l_i --
    ### everything else about this patient's data (rows, other covariates,
    ### NA pattern) is fixed, so it is derived once here and then just has
    ### those cells overwritten in place below. This is also what makes
    ### model.matrix(terms_i, mf) below cheap: since mf already carries a
    ### "terms" attribute, model.matrix() skips re-running model.frame()/
    ### terms() from scratch -- which is what dominated profiling time here,
    ### especially with large l_i grids (i.e. large bandcount1/bandcount2,
    ### see checkBandcountConvergence()).
    ###
    ### terms_i is `lfit[[i]]$terms`, the terms object saved by
    ### longitudinalSub() when the model was fit on the *full* training
    ### data -- not re-derived from this patient's own (small) data. For
    ### ordinary terms this makes no difference, but poly()/splines::ns()/
    ### splines::bs() embed their basis parameters (knots, centering) in the
    ### terms object's `predvars` attribute, and factor() needs the full set
    ### of observed levels; re-deriving terms() from a small per-patient
    ### slice would silently recompute a different (wrong) basis, or error
    ### outright with too few distinct values/levels. Reusing the cached
    ### terms/xlevels/contrasts from the fit, via model.frame(..., xlev=)
    ### and model.matrix(..., contrasts.arg=), makes model.matrix() build
    ### exactly the design matrix the model was fit with, matching
    ### predict.lme()'s own approach.
    terms_i_list = list()
    mf_1_list = list()
    mf_0_list = list()
    n_expected_1 = integer(n_longitudinal)
    n_expected_0 = integer(n_longitudinal)
    for(i in 1:n_longitudinal){
      data_i_1 = data_num_i_list[[i]]
      data_i_0 = data_num_i_list[[i]]

      ### placeholder value for the l_i-dependent columns -- overwritten on
      ### every l_i grid point in the loop below, so its value here is
      ### irrelevant to the result; it only fixes the template's column type.
      data_i_1[survival_variable] = l_i[1]
      data_i_0[survival_variable] = l_i[1]
      if(length(survival_variable_all) != 0){
        for(surv_i in 1 : length(survival_variable_all)){
          data_i_1[survival_variable_all[[surv_i]]] = apply_survival_trans(survival_trans_function[[surv_i]], l_i[1], surv_i)
          data_i_0[survival_variable_all[[surv_i]]] = apply_survival_trans(survival_trans_function[[surv_i]], l_i[1], surv_i)
        }
      }
      data_i_1[event_type_variable] = 1
      data_i_0[event_type_variable] = 0

      if(!(survival_variable %in% all_variables_all[[i]]))
        stop("Error: Condition is false. Please add survival variable to linear mixed model.")

      ### NA in nlme outcome (longitudinal biomarkers), replace with 999
      data_i_1[[outcome_var_all[[i]]]][is.na(data_i_1[[outcome_var_all[[i]]]])] <- 999
      data_i_0[[outcome_var_all[[i]]]][is.na(data_i_0[[outcome_var_all[[i]]]])] <- 999

      ### terms/xlevels/contrasts cached from the fit (see comment above),
      ### not re-derived from this patient's own small slice of data.
      terms_i_list[[i]] = lfit[[i]]$terms
      xlev_i = if (!is.null(long_fit_all$xlevels)) long_fit_all$xlevels[[i]] else NULL
      mf_1_list[[i]] = model.frame(terms_i_list[[i]], data_i_1, xlev = xlev_i)
      mf_0_list[[i]] = model.frame(terms_i_list[[i]], data_i_0, xlev = xlev_i)
      n_expected_1[i] = nrow(data_i_1)
      n_expected_0[i] = nrow(data_i_0)
    }

    ### for loop and calculate the prediction probability for all time points in l_i
    for(it in 1: length(l_i)){

      ### get prediction data matrix for each biomarker for patient 'num_i'
      LME_indi_matrix_1 = list()
      LME_indi_matrix_0 = list()
      for(i in 1:n_longitudinal){
        terms_i = terms_i_list[[i]]

        #survival variable replaced by l_i[it]
        if(survival_variable %in% names(mf_1_list[[i]])){
          mf_1_list[[i]][[survival_variable]] = l_i[it]
          mf_0_list[[i]][[survival_variable]] = l_i[it]
        }

        #transformed survival variable/basis function of survival variable
        #replaced by trans_function(l_i[it])
        if(length(survival_variable_all) != 0){
          for(surv_i in 1 : length(survival_variable_all)){
            svar = survival_variable_all[[surv_i]]
            if(svar %in% names(mf_1_list[[i]])){
              trans_val = apply_survival_trans(survival_trans_function[[surv_i]], l_i[it], surv_i)
              mf_1_list[[i]][[svar]] = trans_val
              mf_0_list[[i]][[svar]] = trans_val
            }
          }
        }

        ## extract data matrix to calcuate the probability
        LME_indi_matrix_1[[i]] = t(model.matrix(terms_i, mf_1_list[[i]], contrasts.arg = lfit[[i]]$contrasts))
        LME_indi_matrix_0[[i]] = t(model.matrix(terms_i, mf_0_list[[i]], contrasts.arg = lfit[[i]]$contrasts))

          ### data missing when extrat the data using model.matrix,
          ### model.matrix will automatic delete the missing data
         if(dim( LME_indi_matrix_1[[i]] )[2] != n_expected_1[i]){
           LME_indi_matrix_1[[i]] = cbind(LME_indi_matrix_1[[i]], matrix(NA,
                dim(LME_indi_matrix_1[[i]] )[1], n_expected_1[i] -
                  dim( LME_indi_matrix_1[[i]] )[2]))
         }
         if(dim( LME_indi_matrix_0[[i]] )[2] != n_expected_0[i]){
           LME_indi_matrix_0[[i]] = cbind(LME_indi_matrix_0[[i]], matrix(NA,
                dim(LME_indi_matrix_0[[i]] )[1], n_expected_0[i] -
                  dim( LME_indi_matrix_0[[i]] )[2]))
         }

      }

      # Initialize empty lists for rows and columns
      rows1 = list()
      columns1 = list()
      rows0 = list()
      columns0 = list()
      # Loop over the LME_indi_matrix_1
      # get prediction data matrix for all biomarkers for patient 'num_i'
      for (i in seq_len(length(LME_indi_matrix_1))) {
        for (j in seq_len(length(LME_indi_matrix_1))) {
          if (i == j) {
            # Add the matrix itself when row and column index are the same
            columns1[[j]] = LME_indi_matrix_1[[i]]
            columns0[[j]] = LME_indi_matrix_0[[i]]
          } else {
            # Add a zero matrix otherwise
            columns1[[j]] = matrix(0, nrow=nrow(LME_indi_matrix_1[[i]]), ncol=ncol(LME_indi_matrix_1[[j]]))
            columns0[[j]] = matrix(0, nrow=nrow(LME_indi_matrix_0[[i]]), ncol=ncol(LME_indi_matrix_0[[j]]))
          }
        }
        # Combine the columns for this row
        rows1[[i]] = do.call(cbind, columns1)
        rows0[[i]] = do.call(cbind, columns0)
        
      }
      # Combine all the rows, combine all individual longitudinal matrix 
      # for all patients using to do prediction 
      LME_all_matrix_1 = do.call(rbind, rows1)
      LME_all_matrix_0 = do.call(rbind, rows0)
      
      A_matrix_1_1_loop = LME_all_matrix_1 %*% Sigma_all_solve %*% t(LME_all_matrix_1)
      A_matrix_2_1_loop = LME_all_matrix_1 %*% Sigma_all_solve %*% longitudinal_all_matrix
      
      A_matrix_1_0_loop = LME_all_matrix_0 %*% Sigma_all_solve %*% t(LME_all_matrix_0)
      A_matrix_2_0_loop = LME_all_matrix_0 %*% Sigma_all_solve %*% longitudinal_all_matrix
      
      para_matrix_A_21 = t(parameter_matrix) %*% A_matrix_2_1_loop
      para_matrix_A_20 = t(parameter_matrix) %*% A_matrix_2_0_loop
      f_Y_T_D_w1[it, iii] = det_Var_cov_estep^{-0.5} * 
        exp(sum(diag(-0.5*( long_sigma_long + 
                              t(parameter_matrix) %*% A_matrix_1_1_loop %*% parameter_matrix - 
                              para_matrix_A_21 - 
                              t(para_matrix_A_21) ) ))) 
      f_Y_T_D_w0[it, iii] = det_Var_cov_estep^{-0.5} * 
        exp(sum(diag(-0.5*( long_sigma_long + 
                              t(parameter_matrix) %*% A_matrix_1_0_loop %*% parameter_matrix - 
                              para_matrix_A_20 - 
                              t(para_matrix_A_20) ) ))) 
    }
    
  }
  return(f_Y_T_D = list(f_Y_T_D_w0, f_Y_T_D_w1))
}
