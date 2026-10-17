#' conditional distribution of D|T
#' 
#' @description This function computes the conditional probability density function of 
#' competing risk event type D, given the survival time T.
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
#' @param survival_fit_all Results and parameters generated from the model fitting 
#' procedure, utilizing the \code{coxph} function. These outputs include the comprehensive 
#' findings and variables derived from the analysis.
#' @param l_i A vector of time points to calculate the conditional probability.
#' 
#' @return Probability matrices of competing risk event type D conditional on 
#' survival outcome T.
#' @keywords internal
conditionalDT = function(data_predict_all, long_fit_all, survival_fit_all, l_i){
  
  ### extract data to calculate the conditional probability
  num <- as.character(nlme::splitFormula(long_fit_all$long_sub_random[[1]], "|")[[2]])[2]
  data.surv =  data_predict_all[[1]][!duplicated(data_predict_all[[1]][num]), ]
  ### censor variable name
  ### time-to-event variable name
  survival_variable = survival_time_variable(survival_fit_all)
  ### event type variable name
  event_type_variable = as.character(formula(survival_fit_all$form_conditional_cr)[[2]])
  
  ### NA in glm outcome (event type), replace with 999
  data.surv[event_type_variable][is.na( data.surv[event_type_variable])] <- 999
  ## data matrix used to calculate the probability
  ### the survival-time column only contributes through l_i below, so its
  ### value here is irrelevant -- but it may be NA for a new patient, and an
  ### all-NA column is logical, which model.matrix() would treat as a factor
  ### (so the survival-time coefficient could no longer be matched by name).
  data.surv[[survival_variable]] <- 0
  ### na.pass keeps one row per patient: a missing covariate should give
  ### that patient NA rather than drop the row and misalign every later
  ### patient.
  mf_cr = model.frame(survival_fit_all$form_conditional_cr, data.surv,
                      na.action = stats::na.pass, xlev = survival_fit_all$glm_fit$xlevels)
  data_matrix_probability = model.matrix(survival_fit_all$form_conditional_cr, mf_cr)
  
  ## fit glm vertical model
  glm_fit = survival_fit_all$glm_fit

  ### the survival time enters the event-type model linearly (checked by
  ### survivalSub()); its effect is evaluated at each l_i below. If it is
  ### not in the model at all, event type does not depend on the event
  ### time: its coefficient is taken as 0 (previously this indexed with an
  ### empty vector and returned NA for every patient).
  is_time_column <- colnames(data_matrix_probability) == survival_variable
  beta_time <- if (any(is_time_column)) glm_fit$coefficients[is_time_column] else 0

  ## covariates * parameter matrix
  covariate_para_matrix = c(data_matrix_probability[, !is_time_column, drop = FALSE] %*%
                              glm_fit$coefficients[!is_time_column])
  n_patients <- dim(data_matrix_probability)[1]

  ### conditional probability
  surv_med_w1 = 1/ (1 +  exp(matrix(beta_time * l_i, length(l_i), n_patients)  +
                                 t(matrix(covariate_para_matrix, n_patients, length(l_i))) ))

  surv_med_w2 = exp(matrix(beta_time * l_i, length(l_i), n_patients)  +
                      t(matrix(covariate_para_matrix, n_patients, length(l_i))) ) * surv_med_w1

  return(list(surv_med_w1, surv_med_w2))
}
