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
  ### the event type (the glm response) is not needed: only the right-hand
  ### side is used below. It used to be read from the data and its missing
  ### values set to 999, which failed when the column was absent (as it is
  ### for a new patient) and warned when it was a factor.
  ### the terms of the glm fit keep the data-dependent bases (e.g. ns()) it
  ### was fit with
  glm_fit = survival_fit_all$glm_fit
  rhs_cr = stats::delete.response(stats::terms(glm_fit))
  n_patients <- nrow(data.surv)
  design_at <- function(d) {
    ### na.pass keeps one row per patient: a missing covariate should give
    ### that patient NA rather than drop the row and misalign every later
    ### patient.
    mf <- model.frame(rhs_cr, d, na.action = stats::na.pass, xlev = glm_fit$xlevels)
    model.matrix(rhs_cr, mf)
  }

  if (time_enters_linearly(rhs_cr, survival_variable)) {
    ### the survival time enters the event-type model linearly or not at
    ### all: its effect is its coefficient times each l_i. The column's
    ### value here is irrelevant -- but it may be NA for a new patient, and
    ### an all-NA column is logical, which model.matrix() would treat as a
    ### factor (so the coefficient could no longer be matched by name).
    data.surv[[survival_variable]] <- 0
    data_matrix_probability = design_at(data.surv)
    ### if the survival time is not in the model, event type does not
    ### depend on event time: its coefficient is taken as 0
    is_time_column <- colnames(data_matrix_probability) == survival_variable
    beta_time <- if (any(is_time_column)) glm_fit$coefficients[is_time_column] else 0
    covariate_para_matrix = c(data_matrix_probability[, !is_time_column, drop = FALSE] %*%
                                glm_fit$coefficients[!is_time_column])
    eta = matrix(beta_time * l_i, length(l_i), n_patients) +
      t(matrix(covariate_para_matrix, n_patients, length(l_i)))
  } else {
    ### a transformation or interaction of the survival time: rebuild the
    ### model matrix with every patient at every l_i (rows ordered by l_i
    ### within patient)
    d <- data.surv[rep(seq_len(n_patients), each = length(l_i)), , drop = FALSE]
    d[[survival_variable]] <- rep(l_i, n_patients)
    eta = matrix(c(design_at(d) %*% glm_fit$coefficients), length(l_i), n_patients)
  }

  ### conditional probability
  ### plogis() rather than 1 / (1 + exp(eta)) and exp(eta) * that, which
  ### gave Inf * 0 = NaN once exp(eta) overflowed
  surv_med_w1 = stats::plogis(-eta)
  surv_med_w2 = stats::plogis(eta)

  return(list(surv_med_w1, surv_med_w2))
}
