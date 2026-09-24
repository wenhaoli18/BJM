#' Dynamic prediction function
#' 
#' The time values in the prediction data subset must be less than the 
#' specified \code{prediction.time} which is the prediction time. The time points for 
#' longitudinal repeated measurements must not surpass the prediction time.
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
#' @param prediction.time Time used to make the prediction.
#' @param horizon Prediction horizon.
#' @param time_variable The name of time variable in linear mixed model.
#' @param survivalVariableAll The name of the transformed time-to-event outcomes variable.
#' @param survivalTransFunction The transformation function used for time-to-event outcomes, 
#' in the order of \code{survivalVariableAll}.
#' @param bandcount1 The number of points used to perform the numerical integral, 
#' from the prediction time to the prediction time plus the horizon.
#' @param bandcount2 The number of points used to perform the numerical integral,
#'  from the prediction time to infinity.
#' 
#' @return An object of class \code{"dynamicPrediction.BJM"}, a named list with elements:
#' \describe{
#'   \item{risk_prob_1}{A vector of dynamically predicted probabilities, one per patient,
#'   of experiencing the (first) event within the prediction horizon. \code{0} when
#'   \code{horizon <= 0}.}
#'   \item{risk_prob_2}{When \code{survival_fit_all} was fit with competing risks, a vector
#'   of dynamically predicted probabilities, one per patient, of experiencing the competing
#'   event within the prediction horizon. \code{NULL} when there is no competing risk, or
#'   when \code{horizon <= 0}.}
#' }
#'
#' @examples 
#' 
#' \donttest{
#' data(pbc3)
#' 
#' data.survival.fitting =  pbc3[!duplicated(pbc3$id), ]
#' 
#' formMarginalSurv = Surv(years, status3) ~ age + sex
#' formConditionalCR = NULL
#' 
#' survival_fit_all = survivalSub(data.survival.fitting, formMarginalSurv, 
#'                                formConditionalCR)
#' 
#' LongSubFixed = list(
#'   "long1" = serBilir ~ year + age + sex +  (years) + (years) * year,  
#'   "long2" = prothrombin ~ year + age + sex + (years) + (years) * year,  
#'   "long3" = albumin ~ year + age + age * year + sex + (years) + (years) * year,  
#'   "long4" = alkaline ~ year + age + sex + (years) + (years) * year, 
#'   "long5" = SGOT ~ year + age + sex + (years) + (years) * year, 
#'   "long6" = platelets ~ year + age + sex + (years)  + (years) * year)
#' 
#' LongSubRandom =list(
#'   "long1" =  ~ year| id,   
#'   "long2" =  ~ year| id,    
#'   "long3" =  ~ year| id,    
#'   "long4" =  ~ year| id,    
#'   "long5" =  ~ year| id,    
#'   "long6" =  ~ year| id)
#' 
#' survivalVariableAll = list(
#'   "Tyears1",  "Tyears2", "Tyears3", "Tyears4"
#' )
#' 
#' survivalTransFunction = list(
#'   fun1 = function(x){abs(x - 1)}, 
#'   fun2 = function(x){abs(x - 3)}, 
#'   fun3 = function(x){abs(x - 5)}, 
#'   fun4 = function(x){abs(x - 7)}
#' )
#' 
#' # Complete case analysis
#' data.fit.all = list()
#' for(i in seq_len(length(LongSubFixed))){
#'   data.fit.all[[i]] = pbc3[pbc3$status3 == 1, ]
#' }
#' 
#' # fitting longitudinal submodel
#' long_fit_all = longitudinalSub(data.fit.all, LongSubFixed, LongSubRandom)
#' 
#' i_PID = 2
#' data.raw.predict.1 = pbc3[pbc3$id == i_PID, ]
#' 
#' data.predict.all = list()
#' for(i in seq_len(length(LongSubFixed))){
#'   data.predict.all[[i]] = data.raw.predict.1[data.raw.predict.1$year <= 3,]
#' }
#' 
#' # predict risk probability
#' risk.prob = dynamicPrediction(data.predict.all, long_fit_all, survival_fit_all, 
#'                               prediction.time = 3, 
#'                               horizon = 3, time_variable = "year",
#'                               survivalVariableAll, survivalTransFunction,
#'                               bandcount1 = 10, bandcount2 = 10)
#' 
#' }
#' 
#' @export
dynamicPrediction = function(data.predict.all, long_fit_all, survival_fit_all, 
                             prediction.time, horizon, time_variable, 
                             survivalVariableAll, survivalTransFunction, 
                             bandcount1 = 10, bandcount2 = 40){
  
  coxph_fit = survival_fit_all$coxph_fit
  survival_variable = as.character(formula(coxph_fit)[[2]])[2] #survival_variable = "fuyrs"
  ## at risk sample
  ## for loop number of biomarkers
  for(i in seq_len(length(data.predict.all))){
    #data.predict.all[[i]] = data.predict.all[[i]][data.predict.all[[i]]$time <= prediction.time, ]
    data.predict.all[[i]] = data.predict.all[[i]][data.predict.all[[i]][survival_variable] >= prediction.time, ]
  }
  
  upper_bound = 2 * max(data.predict.all[[1]][survival_variable])
  
  #### handle horizon = 0 edge case: probability of event in zero-length window is 0
  if(horizon <= 0){
    out <- list(risk_prob_1 = 0, risk_prob_2 = NULL)
    class(out) <- "dynamicPrediction.BJM"
    return(out)
  }
  
  #### time frame used to do the integral
  bandwidth1 = horizon/bandcount1
  predict.time.horizon = seq(prediction.time, prediction.time + horizon, bandwidth1)
  predict.time.horizon.1 = seq(prediction.time - bandwidth1/2, prediction.time + horizon + bandwidth1/2, bandwidth1)
  
  bandwidth2 = (upper_bound - prediction.time)/bandcount2  
  predict.time.infinity = seq(prediction.time, upper_bound, bandwidth2)
  predict.time.infinity.1 = seq(prediction.time - bandwidth2/2, upper_bound + bandwidth2/2, bandwidth2)
  
  risk.prob.0 = risk.prob.1 = NULL
  ### marginal probability T
  S_T_all_predict = marginalT(data.predict.all, long_fit_all, survival_fit_all, l_i = predict.time.horizon.1, upper_bound)
  S_T_all_infinity = marginalT(data.predict.all, long_fit_all, survival_fit_all, l_i = predict.time.infinity.1, upper_bound)
  
  #conditional probability D|T, survival_fit_all$formConditionalCR == formConditionalCR
  #with competing risk
  if(length(survival_fit_all$formConditionalCR) != 0){
    #conditional probability D|T
    D_T_all_predict = conditionalDT(data.predict.all, long_fit_all, survival_fit_all, 
                                    l_i = predict.time.horizon)
    D_T_all_infinity = conditionalDT(data.predict.all, long_fit_all, survival_fit_all, 
                                     l_i = predict.time.infinity)
    #conditional probability Y|D,T
    f_y_D_all_predict = conditionalYDT(data.predict.all, long_fit_all, survival_fit_all, 
                                       l_i = predict.time.horizon, survival_variable, 
                                       time_variable, survivalVariableAll, 
                                       survivalTransFunction)
    f_y_D_all_infinity = conditionalYDT(data.predict.all, long_fit_all, survival_fit_all, 
                                        l_i = predict.time.infinity, survival_variable, 
                                        time_variable, survivalVariableAll, 
                                        survivalTransFunction)
    T.surv.predict.0 = t(f_y_D_all_predict[[1]] * D_T_all_predict[[1]] * S_T_all_predict)
    T.surv.infinity.0 = t(f_y_D_all_infinity[[1]] * D_T_all_infinity[[1]] * S_T_all_infinity)
    T.surv.predict.1 = t(f_y_D_all_predict[[2]] * D_T_all_predict[[2]] * S_T_all_predict)
    T.surv.infinity.1 = t(f_y_D_all_infinity[[2]] * D_T_all_infinity[[2]] * S_T_all_infinity)
    
    risk.prob.0 =  rowSums(T.surv.predict.0) / rowSums(T.surv.infinity.0 + T.surv.infinity.1)
    risk.prob.0[risk.prob.0 > 1] = 1
    risk.prob.0[risk.prob.0 < 0] = 0
    
    risk.prob.1 =  rowSums(T.surv.predict.1) / rowSums(T.surv.infinity.0 + T.surv.infinity.1)
    risk.prob.1[risk.prob.1 > 1] = 1
    risk.prob.1[risk.prob.1 < 0] = 0
    
  }else{
    #without competing risk
    #conditional probability Y|T
    f_y_D_all_predict = conditionalYT(data.predict.all, long_fit_all, l_i = predict.time.horizon, 
                                      survival_variable, time_variable, survivalVariableAll, survivalTransFunction)
    f_y_D_all_infinity = conditionalYT(data.predict.all, long_fit_all, l_i = predict.time.infinity, 
                                       survival_variable, time_variable, survivalVariableAll, survivalTransFunction)
    
    T.surv.predict.0 = t(f_y_D_all_predict[[1]]  * S_T_all_predict)
    T.surv.infinity.0 = t(f_y_D_all_infinity[[1]]  * S_T_all_infinity)
    
    risk.prob.0 =  rowSums(T.surv.predict.0) / rowSums(T.surv.infinity.0 + 1e-20)
    risk.prob.0[risk.prob.0 > 1] = 1
    risk.prob.0[risk.prob.0 < 0] = 0
  }
  
  out <- list(risk_prob_1 = risk.prob.0, risk_prob_2 = risk.prob.1)
  class(out) <- "dynamicPrediction.BJM"
  return(out)
}


