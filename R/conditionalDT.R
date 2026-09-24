#' conditional distribution of D|T
#' 
#' @description This function computes the conditional probability density function of 
#' competing risk event type D, given the survival time T.
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
#' @param survival_fit_all Results and parameters generated from the model fitting 
#' procedure, utilizing the \code{coxph} function. These outputs include the comprehensive 
#' findings and variables derived from the analysis.
#' @param l_i A vector of time points to calculate the conditional probability.
#' 
#' @return Probability matrices of competing risk event type D conditional on 
#' survival outcome T.
#' @keywords internal
conditionalDT = function(data.predict.all, long_fit_all, survival_fit_all, l_i){
  
  coxph_fit = survival_fit_all[[1]]
  ### extract data to calculate the conditional probability
  num <- as.character(nlme::splitFormula(long_fit_all[[4]][[1]], "|")[[2]])[2]
  data.surv =  data.predict.all[[1]][!duplicated(data.predict.all[[1]][num]), ]
  ### censor variable name
  #censor_variable = as.character(formula(coxph_fit)[[2]])[3]
  ### time-to-event variable name
  survival_variable = as.character(formula(coxph_fit)[[2]])[2]
  ### event type variable name
  event_type_variable = as.character(formula(survival_fit_all[[4]])[[2]])
  
  ### NA in glm outcome (event type), replace with 999
  data.surv[event_type_variable][is.na( data.surv[event_type_variable])] <- 999
  ## data matrix used to calculate the probability
  data_matrix_probability = model.matrix(survival_fit_all[[4]], data.surv)
  
  ### is missing in covariates, model.matrix will delete automatically, then we need add NA to data_matrix_probability
  if(dim(data_matrix_probability)[1] != dim(data.surv)[1]){
    data_matrix_probability = rbind(data_matrix_probability, 
                                    matrix(NA, dim(data.surv)[1] - dim(data_matrix_probability)[1],
                                   dim(data_matrix_probability)[2] ))
  }
  
  survival_variable_index <- which(colnames(data_matrix_probability) == survival_variable)
  
  ## fit glm vertical model
  glm_fit = survival_fit_all[[3]]
  
  ## covariates * parameter matrix
  if(is.null(dim(data_matrix_probability[,-survival_variable_index]))){
    ## one sample
    covariate_para_matrix = c(glm_fit$coefficients[-survival_variable_index] %*% data_matrix_probability[,-survival_variable_index])
  }else{
    covariate_para_matrix = c(glm_fit$coefficients[-survival_variable_index] %*% t(data_matrix_probability[,-survival_variable_index]))
  }
  
  ### conditional probability
  surv_med_w1 = 1/ (1 +  exp(matrix(glm_fit$coefficients[survival_variable_index] * l_i, length(l_i), 
                                      dim(data_matrix_probability)[1])  + 
                                 t(matrix(covariate_para_matrix, 
                                          dim(data_matrix_probability)[1], length(l_i))) )) 
  
  surv_med_w2 = exp(matrix(glm_fit$coefficients[survival_variable_index] * l_i, length(l_i), 
                           dim(data_matrix_probability)[1])  + 
                      t(matrix(covariate_para_matrix, 
                               dim(data_matrix_probability)[1], length(l_i))) ) * surv_med_w1
  
  return(list(surv_med_w1, surv_med_w2))
}
