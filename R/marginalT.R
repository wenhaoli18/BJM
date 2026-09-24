#' Marginal distribution of T
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
#' @param upper_bound Upper limit for integration. To manage the hazard function, 
#' extrapolation is required. The \code{upper_bound} parameter specifies the upper time 
#' points at which extrapolation is performed.
#' 
#' @return Marginal density probability of the survival variable T is represented as a 
#' probability matrix. In this matrix, the rows (l_i) are aligned with specific 
#' time points, while the columns correspond to individual patients. 
#' Each entry in the matrix denotes the marginal density probability of 
#' survival for a given patient at a particular time point.
#' @keywords internal
marginalT = function(data.predict.all, long_fit_all, survival_fit_all, l_i, upper_bound){
  
  coxph_fit = survival_fit_all[[1]]
  # survival data frame
  num <- as.character(nlme::splitFormula(long_fit_all[[4]][[1]], "|")[[2]])[2]
  data.surv =  data.predict.all[[1]][!duplicated(data.predict.all[[1]][num]), ]
  
  ## baseline hazard
  cum_basehaz_weibull = basehaz(coxph_fit, centered = FALSE)
  ## extrapolation
  cum_base_lm = lm(hazard ~ time, cum_basehaz_weibull)
  new_data <-  data.frame(time = seq(max(cum_basehaz_weibull), 2 * upper_bound, 0.005))
  new_data_all = data.frame(hazard = c(predict(cum_base_lm, new_data)), time =  new_data)
  cum_basehaz_weibull = rbind(cum_basehaz_weibull, new_data_all)
  #haz_weibull = cbind(diff(cum_basehaz_weibull[,1])/diff(cum_basehaz_weibull[,2]), cum_basehaz_weibull[,2][-1])
  
  ## covariates * parameter matrix
  if(dim(data.surv)[1] == 1){
    ## one sample
    covariate_para_matrix = c(coxph_fit$coefficients  %*%  model.matrix(survival_fit_all[[2]], data.surv)[,-1])
  }else{
    covariate_para_matrix = c(coxph_fit$coefficients  %*%  t(model.matrix(survival_fit_all[[2]], data.surv)[,-1]))
  }
  
  ### cumulative survival
  surv_med = exp(-cum_basehaz_weibull[,1][apply(abs(t(matrix(t(l_i), length(l_i), length(cum_basehaz_weibull[,2]))) -
                                                      matrix(cum_basehaz_weibull[,2], length(cum_basehaz_weibull[,2]), length(l_i))), 2, which.min)] %*% 
                   t(exp(covariate_para_matrix) ))
  ### density
  S_T_all = - diff(surv_med)
  return(S_T_all)
}
