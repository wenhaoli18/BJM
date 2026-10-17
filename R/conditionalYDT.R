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

  patient_ids <- unlist(data.long[[1]][!duplicated(data.long[[1]][num]), ][num])

  ### each biomarker's mean X(t) beta for every row and grid point at once,
  ### per event type (see fixed_mean_grid()); per patient, only the
  ### quadratic form is left
  if (length(patient_ids) > 0) {
    for (i in seq_len(n_longitudinal)) {
      assert_uses_event_time(formula(lfit[[i]]), i, survival_variable, survival_variable_all)
    }
  }
  data_index <- if (length(data.long) == 1) rep(1, n_longitudinal) else seq_len(n_longitudinal)
  mean_grid <- function(event_type) {
    lapply(seq_len(n_longitudinal), function(i) {
      xlev_i <- if (!is.null(long_fit_all$xlevels)) long_fit_all$xlevels[[i]] else NULL
      fixed_mean_grid(data.long[[data_index[i]]], lfit[[i]], xlev_i, l_i, survival_variable,
                      survival_variable_all, survival_trans_function,
                      fixed_values = stats::setNames(list(event_type), event_type_variable))
    })
  }
  mean_grid_1 <- mean_grid(1)
  mean_grid_0 <- mean_grid(0)

  iii = 0
  for(num_i in patient_ids){
    iii = iii + 1
    patient_data <- select_patient_longitudinal_data(data.long, num, num_i, n_longitudinal, time_variable)
    if(is.null(patient_data)) next

    design <- build_conditional_design(patient_data$rep_num_i_list, patient_data$data_num_i_list,
                                        lfit, Sigma, sigma.longitudinal, time_variable, n_longitudinal,
                                        long_fit_all$long_sub_random)

    ### (Y - mu(t)) for the full stacked observation vector (all biomarkers),
    ### one column per grid point and event type; the log density uses the
    ### quadratic form with the patient's whole covariance, cross-biomarker
    ### terms included (see NEWS.md)
    rows <- lapply(seq_len(n_longitudinal), function(i)
      which(as.character(data.long[[data_index[i]]][[num]]) == as.character(num_i)))
    stack_mean <- function(grid) {
      do.call(rbind, lapply(seq_len(n_longitudinal), function(i) grid[[i]][rows[[i]], , drop = FALSE]))
    }
    y <- rowSums(design$longitudinal_all_matrix)
    ### log densities (see predictRisk() for how they are exponentiated safely)
    f_Y_T_D_w1[, iii] = cov_logdens(design$cov_fac, y - stack_mean(mean_grid_1))
    f_Y_T_D_w0[, iii] = cov_logdens(design$cov_fac, y - stack_mean(mean_grid_0))
  }
  return(f_Y_T_D = list(f_Y_T_D_w0, f_Y_T_D_w1))
}
