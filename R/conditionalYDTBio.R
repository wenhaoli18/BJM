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
#' @return The output is a list containing probability matrices of \strong{log} densities
#' (so that they neither overflow nor underflow; see \code{predictRisk()} for how they
#' are exponentiated). In the presence of 
#' competing risks, this list includes two elements; otherwise, 
#' it contains only one element. Each element within the list is a probability matrix, 
#' with the number of rows (l_i) corresponding to specific time points and 
#' columns representing different patients. Every matrix element represents 
#' the conditional probability derived from the conditional distribution 
#' of longitudinal variable Y given the survival time T with competing risk D
#' for a particular patient at a specific time point.
#' @keywords internal
conditionalYDTBio = function(Y_all, time_new, bio_i, data_predict_all, 
                             long_fit_all, 
                              survival_fit_all, l_i, survival_variable, 
                              time_variable, survival_variable_all, survival_trans_function){
  
  #LME model fitting
  lfit = long_fit_all$lfit
  #patient ID
  num <- as.character(nlme::splitFormula(long_fit_all$long_sub_random[[1]], "|")[[2]])[2]
  ### event type variable name
  event_type_variable = as.character(formula(survival_fit_all$form_conditional_cr)[[2]])
  
  #number of longitudinal biomarkers
  n_longitudinal <- length(lfit)  #length(data_num_i_list)
  
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
  cov_fac_all <- lapply(unique(data.long[[1]][[num]]), process_variance, 
                      time_new, bio_i, data_predict_all, long_fit_all, time_variable)
  
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

    patient_data <- select_patient_longitudinal_data_bio(data.long, num, num_i, n_longitudinal, time_variable,
                                                          bio_i, time_new, Y_all, long_fit_all)
    if(is.null(patient_data)) next
    data_num_i_list <- patient_data$data_num_i_list
    Y_select_all <- patient_data$Y_select_all

    ####Constructing the longitudinal matrix for all biomarkers, grid search for all Y_all
    longitudinal_all_matrix <- build_longitudinal_matrix_bio(data_num_i_list, lfit, bio_i, Y_select_all,
                                                              n_longitudinal, Y_all)

    ### every biomarker's mean at every grid point, per event type (see
    ### fixed_mean_grid()), then the log density of every candidate value at
    ### every grid point; the covariance is the same at every grid point and
    ### for both event types, so it is factored once per patient (by
    ### process_variance())
    for (i in seq_len(n_longitudinal)) {
      assert_uses_event_time(formula(lfit[[i]]), i, survival_variable, survival_variable_all)
    }
    mean_at <- function(event_type) {
      patient_mean_grid(data_num_i_list, long_fit_all, l_i, survival_variable,
                        survival_variable_all, survival_trans_function,
                        fixed_values = stats::setNames(list(event_type), event_type_variable))
    }
    log_dens_1 <- cov_logdens_grid(cov_fac_all[[iii]], longitudinal_all_matrix, mean_at(1))
    log_dens_0 <- cov_logdens_grid(cov_fac_all[[iii]], longitudinal_all_matrix, mean_at(0))

    for(Y_i in 1 : length(Y_all)){
      f_Y_T_D_w1[[Y_i]][,iii] = log_dens_1[Y_i, ]
      f_Y_T_D_w0[[Y_i]][,iii] = log_dens_0[Y_i, ]
    }

  }
  
  return(f_Y_T_D = list(f_Y_T_D_w0, f_Y_T_D_w1))
}
