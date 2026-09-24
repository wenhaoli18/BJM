#' Dynamic prediction function for future biomarker
#' 
#' The time values in the prediction data subset must be less than the 
#' specified \code{prediction_time} which is the prediction time. The time points for 
#' longitudinal repeated measurements must not surpass the prediction time.
#' 
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
#' @param survival_fit_all Results and parameters generated from the model fitting 
#' procedure, utilizing the \code{coxph} function. These outputs include the comprehensive 
#' findings and variables derived from the analysis.
#' @param prediction_time Time used to make the prediction
#' @param horizon Prediction horizon
#' @param time_variable The name of time variable in linear mixed model.
#' @param survival_variable_all The name of the transformed time-to-event outcomes variable.
#' @param survival_trans_function The transformation function used for time-to-event outcomes, 
#' in the order of \code{survival_variable_all}.
#' @param bandcount2 The number of points used to perform the numerical integral,
#'  from the prediction time to infinity.
#' @param bandcount3 The number of points used to calculate the probability density function.
#' 
#' @return An object of class \code{"dynamicPredictionBio.BJM"}, a named list with elements:
#' \describe{
#'   \item{Y_predict}{A vector, one entry per patient, giving the MAP (most likely) predicted
#'   value of biomarker \code{bio_i} at \code{prediction_time + horizon}.}
#'   \item{Y_density}{A probability matrix whose rows correspond to the candidate biomarker
#'   values in \code{Y_all} and whose columns correspond to individual patients; each entry
#'   is the dynamically predicted density of the biomarker taking that value.}
#'   \item{Y_all}{The grid of candidate biomarker values used to build \code{Y_density}.}
#' }
#'
#' @examples 
#' 
#' \donttest{
#' data(pbc3)
#' 
#' data_survival_fitting =  pbc3[!duplicated(pbc3$id), ]
#' 
#' form_marginal_surv = Surv(years, status3) ~ age + sex
#' form_conditional_cr = NULL
#' 
#' survival_fit_all = survivalSub(data_survival_fitting, form_marginal_surv, 
#'                                form_conditional_cr)
#' 
#' long_sub_fixed = list(
#'   "long1" = serBilir ~ year + age + sex +  (years) + (years) * year,  
#'   "long2" = prothrombin ~ year + age + sex + (years) + (years) * year,  
#'   "long3" = albumin ~ year + age + age * year + sex + (years) + (years) * year,  
#'   "long4" = alkaline ~ year + age + sex + (years) + (years) * year, 
#'   "long5" = SGOT ~ year + age + sex + (years) + (years) * year, 
#'   "long6" = platelets ~ year + age + sex + (years)  + (years) * year)
#' 
#' long_sub_random =list(
#'   "long1" =  ~ year| id,   
#'   "long2" =  ~ year| id,    
#'   "long3" =  ~ year| id,    
#'   "long4" =  ~ year| id,    
#'   "long5" =  ~ year| id,    
#'   "long6" =  ~ year| id)
#' 
#' survival_variable_all = list(
#'   "Tyears1",  "Tyears2", "Tyears3", "Tyears4"
#' )
#' 
#' survival_trans_function = list(
#'   fun1 = function(x){abs(x - 1)}, 
#'   fun2 = function(x){abs(x - 3)}, 
#'   fun3 = function(x){abs(x - 5)}, 
#'   fun4 = function(x){abs(x - 7)}
#' )
#' 
#' # Complete case analysis
#' data_fit_all = list()
#' for(i in seq_len(length(long_sub_fixed))){
#'   data_fit_all[[i]] = pbc3[pbc3$status3 == 1, ]
#' }
#' 
#' # fitting longitudinal submodel
#' long_fit_all = longitudinalSub(data_fit_all, long_sub_fixed, long_sub_random)
#' 
#' i_PID = 2
#' data.raw.predict.1 = pbc3[pbc3$id == i_PID, ]
#' 
#' data_predict_all = list()
#' for(i in seq_len(length(long_sub_fixed))){
#'   data_predict_all[[i]] = data.raw.predict.1[data.raw.predict.1$year <= 3,]
#' }
#' 
#' Y_predict = dynamicPredictionBio(bio_i = 1, data_predict_all, long_fit_all, 
#'                                  survival_fit_all, prediction_time = 3, 
#'                                  horizon = 3, time_variable = "year",
#'                                  survival_variable_all, survival_trans_function,
#'                                  bandcount2 = 40, bandcount3 = 400)
#' 
#' }
#' 
#' @export
dynamicPredictionBio = function(bio_i, data_predict_all, long_fit_all, survival_fit_all, 
                                 prediction_time, horizon, time_variable, 
                                 survival_variable_all, survival_trans_function, 
                                 bandcount2 = 40, bandcount3 = 300){
  
  coxph_fit = survival_fit_all$coxph_fit
  survival_variable = as.character(formula(coxph_fit)[[2]])[2] #survival_variable = "fuyrs"
  ## at risk sample
  ## for loop number of biomarkers
  for(i in seq_len(length(data_predict_all))){
    data_predict_all[[i]] = data_predict_all[[i]][data_predict_all[[i]][survival_variable] >= prediction_time, ]
  }
  
  upper_bound = 2 * max(data_predict_all[[1]][survival_variable])
  #### time frame used to do the integral
  bandwidth2 = (upper_bound - prediction_time)/bandcount2   
  predict.time.infinity = seq(prediction_time, upper_bound, bandwidth2)
  predict.time.infinity.1 = seq(prediction_time - bandwidth2/2, upper_bound + bandwidth2/2, bandwidth2)
  
  risk.prob.0 = risk.prob.1 = NULL
  ### marginal probability T
  S_T_all_infinity = marginalT(data_predict_all, long_fit_all, survival_fit_all, l_i = predict.time.infinity.1, upper_bound)
  

  Y_upper = max(data_predict_all[[bio_i]][as.character(formula(long_fit_all$long_sub_fixed[[bio_i]])[[2]])], na.rm = TRUE)
  Y_lower = min(data_predict_all[[bio_i]][as.character(formula(long_fit_all$long_sub_fixed[[bio_i]])[[2]])], na.rm = TRUE)
  Y_all = seq(Y_lower - 5 * (Y_upper - Y_lower), Y_upper + 5 * (Y_upper - Y_lower), 
              11 * (Y_upper - Y_lower)/bandcount3) #seq(0, 100, 5)
  
  Y_density = c()
  #conditional probability D|T, survival_fit_all$form_conditional_cr == form_conditional_cr
  #with competing risk
  if(length(survival_fit_all$form_conditional_cr) != 0){
    #conditional probability D|T
    D_T_all_infinity = conditionalDT(data_predict_all, long_fit_all, survival_fit_all, 
                                     l_i = predict.time.infinity)
    #conditional probability Y|D,T
    f_y_D_all_infinity = conditionalYDT(data_predict_all, long_fit_all, survival_fit_all, 
                                        l_i = predict.time.infinity, survival_variable, 
                                        time_variable, survival_variable_all, 
                                        survival_trans_function)
    
    f_y_D_all_predict = conditionalYDTBio(Y_all, time_new = prediction_time + horizon, 
                                          bio_i, data_predict_all, long_fit_all, 
                                          survival_fit_all, 
                                          l_i = predict.time.infinity, survival_variable, 
                                          time_variable, survival_variable_all, 
                                          survival_trans_function)

    for(Y_i in seq_len(length(Y_all))){
      T.surv.predict.0 = t(f_y_D_all_predict[[1]][[Y_i]] * D_T_all_infinity[[1]] * S_T_all_infinity)
      T.surv.predict.1 = t(f_y_D_all_predict[[2]][[Y_i]] * D_T_all_infinity[[2]] * S_T_all_infinity)
      T.surv.infinity.0 = t(f_y_D_all_infinity[[1]] * D_T_all_infinity[[1]] * S_T_all_infinity)
      T.surv.infinity.1 = t(f_y_D_all_infinity[[2]] * D_T_all_infinity[[2]] * S_T_all_infinity)
      
      risk.prob.1 = rowSums(T.surv.predict.1) / rowSums(T.surv.infinity.1 + T.surv.infinity.0)
      risk.prob.1[risk.prob.1 > 1] = 1
      risk.prob.1[risk.prob.1 < 0] = 0
      risk.prob.0 =  rowSums(T.surv.predict.0) / rowSums(T.surv.infinity.1 + T.surv.infinity.0)
      risk.prob.0[risk.prob.0 > 1] = 1
      risk.prob.0[risk.prob.0 < 0] = 0
      Y_density = rbind(Y_density, risk.prob.1 + risk.prob.0)
    }
    
    Y_predict = c()
    for(i in 1:dim(Y_density)[2] ){
      if(length(Y_all[which.max(Y_density[,i])]) == 0){
        Y_predict = c(Y_predict, NA)
      }else{
        Y_predict = c(Y_predict, Y_all[which.max(Y_density[,i])])
      }
    }
    
  }else{ #without competing risk
    
    #conditional probability Y|T
    f_y_D_all_infinity = conditionalYT(data_predict_all, long_fit_all, l_i = predict.time.infinity, 
                                       survival_variable, time_variable, survival_variable_all, survival_trans_function)
    f_y_D_all_predict = conditionalYTBio(Y_all, time_new = prediction_time + horizon, 
                                          bio_i, data_predict_all, long_fit_all, 
                                          l_i = predict.time.infinity, survival_variable, 
                                          time_variable, survival_variable_all, 
                                          survival_trans_function)
    
    for(Y_i in seq_len(length(Y_all))){
      
      T.surv.predict.0 = t(f_y_D_all_predict[[1]][[Y_i]] * S_T_all_infinity)
      T.surv.infinity.0 = t(f_y_D_all_infinity[[1]] * S_T_all_infinity)

      risk.prob.0 =  rowSums(T.surv.predict.0) / rowSums(T.surv.infinity.0 + 1e-20)
      risk.prob.0[risk.prob.0 > 1] = 1
      risk.prob.0[risk.prob.0 < 0] = 0
      Y_density = rbind(Y_density, risk.prob.0)

    }
    
    Y_predict = c()
    for(i in 1:dim(Y_density)[2] ){
      if(length(Y_all[which.max(Y_density[,i])]) == 0){
        Y_predict = c(Y_predict, NA)
      }else{
        Y_predict = c(Y_predict, Y_all[which.max(Y_density[,i])])
      }
    }
    
  }
  
  out <- list(Y_predict = Y_predict, Y_density = Y_density, Y_all = Y_all)
  class(out) <- "dynamicPredictionBio.BJM"
  return(out)
}

