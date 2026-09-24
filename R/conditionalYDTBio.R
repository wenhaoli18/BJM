#' conditional distribution of Y|D, T, if with competing risk
#' 
#' @description This function computes the conditional probability density function of 
#' longitudinal variable Y, given the survival time T with competing risk D.
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
#' of longitudinal variable Y given the survival time T with competing risk D
#' for a particular patient at a specific time point.
#' @keywords internal
conditionalYDTBio = function(Y_all, time_new, bio_i, data.predict.all, 
                             long_fit_all, 
                              survival_fit_all, l_i, survival_variable, 
                              time_variable, survivalVariableAll, survivalTransFunction){
  
  #LME model fitting
  lfit = long_fit_all[[1]]
  #variance-covariance matrix
  Sigma = long_fit_all[[2]]
  #patient ID
  num <- as.character(nlme::splitFormula(long_fit_all[[4]][[1]], "|")[[2]])[2]
  ### event type variable name
  event_type_variable = as.character(formula(survival_fit_all[[4]])[[2]])
  
  #number of longitudinal biomarkers
  n_longitudinal <- length(lfit)  #length(data_num_i_list)
  
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
 
  # MVN variance 
  # Apply the function over each unique num using lapply for variance list
  Sigma_all <- lapply(as.numeric(unlist(unique(data.long[[1]][num]))), process_variance, 
                      time_new, bio_i, data.predict.all, long_fit_all, time_variable)
  
  # A probability matrix, 
  # with the number of rows (l_i) corresponding to specific time points and 
  # columns representing different patients. 
  f_Y_T_D_w1 = list(); f_Y_T_D_w0 = list()
  for(Y_i in 1 : length(Y_all)){
    f_Y_T_D_w1[[Y_i]] = matrix(NA, length(l_i), length(unlist(data.long[[1]][!duplicated(data.long[[1]][num]), ][num])))
    f_Y_T_D_w0[[Y_i]] = matrix(NA, length(l_i), length(unlist(data.long[[1]][!duplicated(data.long[[1]][num]), ][num])))
  }
  
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
    
    # Iterate over each data frame for different biomarkers
    for (i in 1:n_longitudinal) {
      df <- data.long[[i]]
      # Extract data for patient ID of 'num_i',  where 'num' equals 'num_i'
      selected_data <- df[df[num] == num_i, ]
      
      # the biomarker used to predict,
      if(i == bio_i){
        selected_data <- rbind(selected_data, selected_data[nrow(selected_data), ])
        #replace time variable with predict time
        selected_data[time_variable][nrow(selected_data),] = time_new
        #name of biomarker
        bio_i_name = as.character(formula(long_fit_all[[3]][[i]])[[2]]) 
        Y_select_all = c() ### matrix for all Y_all for predicted biomarker
        for(Y_new in Y_all){
          selected_data[bio_i_name][nrow(selected_data),] = Y_new
          Y_select_all = cbind(Y_select_all, unlist(selected_data[bio_i_name]))
        }
      }
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
    longitudinal_all_matrix <- c() 
    ####Constructing the longitudinal matrix for all biomarkers
    ####Grid search for all Y_all
    for(Y_i in 1 : length(Y_all)){
      longitudinal_all_matrix_tran <- c() 
      for (i in 1:n_longitudinal) {
        longname = as.character(formula(lfit[[i]]))[2]
        if(i == bio_i){
          longitudinal_all_matrix_tran = c(longitudinal_all_matrix_tran, Y_select_all[,Y_i])
        }else{
          longitudinal_all_matrix_tran = c(longitudinal_all_matrix_tran, unlist(data_num_i_list[[i]][,longname]) )
        }
      }
      longitudinal_all_matrix = rbind(longitudinal_all_matrix, longitudinal_all_matrix_tran)
    }

    
    #### MVN mean function
    Amean_list1 = list()
    Amean_list0 = list()
    ### for loop and calculate the prediction probability for all time points in l_i
    for(it in 1: length(l_i)){
      
      ### get prediction data matrix for each biomarker for patient 'num_i'
      LME_indi_matrix_1 = list()
      LME_indi_matrix_0 = list()
      for(i in 1:n_longitudinal){
        model_formula = formula(lfit[[i]]) #lfit[[1]]
        #terms_model <- terms(model_formula)
        #variable_names <- attr(terms_model, "term.labels")
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

        #event type indicator replaced by 1/0
        data_num_i_list_1 = data_num_i_list_0 = data_num_i_list
        data_num_i_list_1[[i]][event_type_variable] = 1
        data_num_i_list_0[[i]][event_type_variable] = 0
        
        target_covariate = survival_variable 
        ### if fuyrs exists or not
        if(target_covariate %in% all_variables != TRUE)
          stop("Error: Condition is false. Please add survival variable to linear mixed model.")
        else
          ### NA in nlme outcome (longitudinal biomarkers), replace with 999
          data_num_i_list_1[[i]][as.character(formula(long_fit_all[[3]][[i]])[[2]])][is.na(data_num_i_list_1[[i]][as.character(formula(long_fit_all[[3]][[i]])[[2]])])] <- 999
          data_num_i_list_0[[i]][as.character(formula(long_fit_all[[3]][[i]])[[2]])][is.na(data_num_i_list_0[[i]][as.character(formula(long_fit_all[[3]][[i]])[[2]])])] <- 999
          
          ## extract data matrix to calcuate the probability
          LME_indi_matrix_1[[i]] = t(model.matrix(long_fit_all[[3]][[i]], data_num_i_list_1[[i]]))
          LME_indi_matrix_0[[i]] = t(model.matrix(long_fit_all[[3]][[i]], data_num_i_list_0[[i]]))
        
          ### data missing when extract the data using model.matrix, 
          ### model.matrix will automatic delete the missing data
         if(dim( LME_indi_matrix_1[[i]] )[2] != dim(data_num_i_list_1[[i]])[1]){
           LME_indi_matrix_1[[i]] = cbind(LME_indi_matrix_1[[i]], matrix(NA, 
                dim(LME_indi_matrix_1[[i]] )[1], dim(data_num_i_list_1[[i]])[1] - 
                  dim( LME_indi_matrix_1[[i]] )[2]))
           LME_indi_matrix_0[[i]] = cbind(LME_indi_matrix_0[[i]], matrix(NA, 
                dim(LME_indi_matrix_0[[i]] )[1], dim(data_num_i_list_1[[i]])[1] - 
                  dim( LME_indi_matrix_0[[i]] )[2]))
         }
           
      }
      
      mean_list1 = c()
      for (i in 1:n_longitudinal) {
        mean_list1 = c(mean_list1, t(LME_indi_matrix_1[[i]]) %*% lfit[[i]]$coefficients$fixed)
      }
      Amean_list1[[it]] = mean_list1
      
      mean_list0 = c()
      for (i in 1:n_longitudinal) {
        mean_list0 = c(mean_list0, t(LME_indi_matrix_0[[i]]) %*% lfit[[i]]$coefficients$fixed)
      }
      Amean_list0[[it]] = mean_list0
    }
    
    results_lapply1 <- lapply(seq_along(Amean_list1), function(it) {
      mvtnorm::dmvnorm(x = longitudinal_all_matrix, mean = c(Amean_list1[[it]]), sigma = Sigma_all[[iii]])
    })
    results_lapply0 <- lapply(seq_along(Amean_list0), function(it) {
      mvtnorm::dmvnorm(x = longitudinal_all_matrix, mean = c(Amean_list0[[it]]), sigma = Sigma_all[[iii]])
    })
    
    for(Y_i in 1 : length(Y_all)){
      ### get #Y_i from all elements of a list
      f_Y_T_D_w1[[Y_i]][,iii] = unlist(lapply(results_lapply1, function(x) x[Y_i]))
      f_Y_T_D_w0[[Y_i]][,iii] = unlist(lapply(results_lapply0, function(x) x[Y_i]))
    }

  }
  
  return(f_Y_T_D = list(f_Y_T_D_w0, f_Y_T_D_w1))
}
