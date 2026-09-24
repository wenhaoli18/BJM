#' conditional distribution of Y|T, if no competing risk; 
#' 
#' @description This function computes the conditional probability density function of 
#' longitudinal variable Y, given the survival time T without competing risk D.
#' 
#' @param time_new Prediction time add horizon
#' @param bio_i Biomarker used to do prediction
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
#' of longitudinal variable Y given the survival time T without competing risk D
#' for a particular patient at a specific time point.
#' @keywords internal
conditionalYTBio = function(Y_all, time_new, bio_i, data_predict_all, 
                            long_fit_all, l_i, survival_variable, 
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
  
  # MVN variance 
  # Apply the function over each unique num using lapply for variance list
  Sigma_all <- lapply(as.numeric(unlist(unique(data.long[[1]][num]))), process_variance, 
                      time_new, bio_i, data_predict_all, long_fit_all, time_variable)
  
  # A probability matrix, 
  # with the number of rows (l_i) corresponding to specific time points and 
  # columns representing different patients. 
  f_Y_T_D_w1 = list(); f_Y_T_D_w0 = list()
  for(Y_i in 1 : length(Y_all)){
    f_Y_T_D_w1[[Y_i]] = matrix(NA, length(l_i), length(unlist(data.long[[1]][!duplicated(data.long[[1]][num]), ][num])))
  }
  
  iii = 0
  for(num_i in unlist(data.long[[1]][!duplicated(data.long[[1]][num]), ][num])){
    iii = iii + 1
    #if(data.raw.sim.1[data.raw.sim.1[num] == num_i,]$status[1] == 1) next
    ### each rows represent intercept, slope, covariates numbers, time to event/l_i
    ### each columns represent different repeated measurements with different times.

    patient_data <- select_patient_longitudinal_data_bio(data.long, num, num_i, n_longitudinal, time_variable,
                                                          bio_i, time_new, Y_all, long_fit_all)
    if(is.null(patient_data)) next
    rep_num_i_list <- patient_data$rep_num_i_list
    data_num_i_list <- patient_data$data_num_i_list
    Y_select_all <- patient_data$Y_select_all

    ####Constructing the longitudinal matrix for all biomarkers, grid search for all Y_all
    longitudinal_all_matrix <- build_longitudinal_matrix_bio(data_num_i_list, lfit, bio_i, Y_select_all,
                                                              n_longitudinal, Y_all)

    #### MVN mean function
    Amean_list1 = list()
    #Amean_list0 = list()
    ### for loop and calculate the prediction probability for all time points in l_i
    for(it in 1: length(l_i)){
      
      ### get prediction data matrix for each biomarker for patient 'num_i'
      LME_indi_matrix_1 = list()
      #LME_indi_matrix_0 = list()
      for(i in 1:n_longitudinal){
        model_formula = formula(lfit[[i]]) #lfit[[1]]
        ### terms cached from the fit on the *full* training data, not
        ### re-derived from this patient's own small slice -- see
        ### conditionalYDT.R for why this matters for poly()/
        ### splines::ns()/splines::bs()/factor() terms.
        terms_model <- lfit[[i]]$terms
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
        
        #event type indicator replaced by 1/0
        data_num_i_list_1 = data_num_i_list_0 = data_num_i_list
        #data_num_i_list_1[[i]][event_type_variable] = 1
        #data_num_i_list_0[[i]][event_type_variable] = 0
        
        target_covariate = survival_variable 
        ### if fuyrs exists or not
        if(target_covariate %in% all_variables != TRUE)
          stop("Error: Condition is false. Please add survival variable to linear mixed model.")
        else
          ### NA in nlme outcome (longitudinal biomarkers), replace with 999
          data_num_i_list_1[[i]][as.character(formula(long_fit_all$long_sub_fixed[[i]])[[2]])][is.na(data_num_i_list_1[[i]][as.character(formula(long_fit_all$long_sub_fixed[[i]])[[2]])])] <- 999
        data_num_i_list_0[[i]][as.character(formula(long_fit_all$long_sub_fixed[[i]])[[2]])][is.na(data_num_i_list_0[[i]][as.character(formula(long_fit_all$long_sub_fixed[[i]])[[2]])])] <- 999
        
        ## extract data matrix to calcuate the probability
        xlev_i = if (!is.null(long_fit_all$xlevels)) long_fit_all$xlevels[[i]] else NULL
        mf_i_1 = model.frame(terms_model, data_num_i_list_1[[i]], xlev = xlev_i)
        LME_indi_matrix_1[[i]] = t(model.matrix(terms_model, mf_i_1, contrasts.arg = lfit[[i]]$contrasts))

        ### data missing when extract the data using model.matrix, 
        ### model.matrix will automatic delete the missing data
        if(dim( LME_indi_matrix_1[[i]] )[2] != dim(data_num_i_list_1[[i]])[1]){
          LME_indi_matrix_1[[i]] = cbind(LME_indi_matrix_1[[i]], matrix(NA, 
                                                                        dim(LME_indi_matrix_1[[i]] )[1], dim(data_num_i_list_1[[i]])[1] - 
                                                                          dim( LME_indi_matrix_1[[i]] )[2]))
        }
        
      }
      
      mean_list1 = c()
      for (i in 1:n_longitudinal) {
        mean_list1 = c(mean_list1, t(LME_indi_matrix_1[[i]]) %*% lfit[[i]]$coefficients$fixed)
      }
      Amean_list1[[it]] = mean_list1
    
    }
    
    results_lapply1 <- lapply(seq_along(Amean_list1), function(it) {
      mvtnorm::dmvnorm(x = longitudinal_all_matrix, mean = c(Amean_list1[[it]]), sigma = Sigma_all[[iii]])
    })
    
    for(Y_i in 1 : length(Y_all)){
      f_Y_T_D_w1[[Y_i]][,iii] = unlist(lapply(results_lapply1, function(x) x[Y_i]))
    }
    
  }
  return(f_Y_T_D = list(f_Y_T_D_w1))
}
