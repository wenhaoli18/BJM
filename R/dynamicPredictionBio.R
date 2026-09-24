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
#' be used in a \code{list}. Alternatively, a single bare \code{data.frame}
#' (not wrapped in a list) may be supplied directly; it is then reused for
#' every longitudinal outcome. It is presumed that each data frame is
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
#' @param bandcount2 The number of grid points spanning
#' \code{[prediction_time, upper_bound]}, where \code{upper_bound} is set
#' internally to twice the longest observed survival/censoring time among
#' at-risk patients; this approximates integrating out to infinity for the
#' denominator that normalizes the predicted density. A wider follow-up range
#' needs a larger \code{bandcount2} to keep the grid spacing comparable.
#' Defaults to \code{"auto"} (see Details).
#' @param bandcount3 The number of points in the candidate-biomarker-value
#' grid (\code{Y_all}) used to build the predicted density
#' (\code{Y_density}) and locate its mode (\code{Y_predict}). This controls
#' the resolution of the density curve, not a time integral; increase it if
#' the density looks jagged or \code{Y_predict} jumps erratically between
#' nearby grid points. Defaults to \code{"auto"} (see Details).
#'
#' @details
#' There is no universal correct value for \code{bandcount2}/\code{bandcount3}:
#' as a practical check, double both and confirm the results barely change;
#' if they do, keep doubling. By default (\code{bandcount2 = "auto"},
#' \code{bandcount3 = "auto"}), this doubling check is done for you:
#' starting from small built-in values, both are doubled together, and the
#' result is compared to the previous round, until the largest relative
#' change in \code{Y_predict} drops below 1%, or 2 doublings have been
#' tried (so at most 3 calls' worth of work). If it still has not converged
#' by then, a warning reports this and the result at the largest value
#' tried is returned anyway (not an error), so this never silently loops
#' for an unbounded amount of time. Pass an explicit number for either
#' argument to skip auto-tuning it and use a fixed value instead (as in
#' previous package versions), or call \code{checkBandcountConvergence()}
#' directly for more control over the tolerance and doubling count. See
#' also \code{vignette("BJM-intro", package = "BJM")} for a worked example.
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
                                 bandcount2 = "auto", bandcount3 = "auto"){

  assert_class(long_fit_all, "longitudinalSub.BJM", "long_fit_all", "longitudinalSub")
  assert_class(survival_fit_all, "survivalSub.BJM", "survival_fit_all", "survivalSub")
  assert_index(bio_i, length(long_fit_all$lfit), "bio_i", "longitudinal outcomes in long_fit_all")
  assert_data_list(data_predict_all, "data_predict_all", length(long_fit_all$lfit), allow_bare_df = TRUE)
  if (!is.list(data_predict_all) || is.data.frame(data_predict_all)) {
    data_predict_all <- rep(list(data_predict_all), each = length(long_fit_all$lfit))
  }
  assert_scalar_numeric(prediction_time, "prediction_time")
  assert_scalar_numeric(horizon, "horizon")
  assert_string(time_variable, "time_variable")
  assert_bandcount(bandcount2, "bandcount2")
  assert_bandcount(bandcount3, "bandcount3")
  assert_survival_trans(survival_variable_all, survival_trans_function, probe_value = prediction_time)

  # bandcount2/bandcount3 = "auto" (the default): see dynamicPrediction()
  # for the rationale; same doubling-until-stable check, applied here to
  # Y_predict instead of the risk probabilities. `call_args` is captured
  # before `auto_names` is computed, so it is not swept up into `call_args`
  # too (see dynamicPrediction() for the same pattern).
  call_args <- as.list(environment())
  auto_names <- c("bandcount2", "bandcount3")[c(identical(bandcount2, "auto"), identical(bandcount3, "auto"))]
  if (length(auto_names) > 0) {
    return(auto_tune_bandcount(dynamicPredictionBio, call_args, auto_names)$result)
  }

  coxph_fit = survival_fit_all$coxph_fit
  survival_variable = as.character(formula(coxph_fit)[[2]])[2] #survival_variable = "fuyrs"
  for (i in seq_along(data_predict_all)) {
    assert_vars_in_data(time_variable, data_predict_all[[i]],
                         "time_variable", sprintf("data_predict_all[[%d]]", i))
    assert_vars_in_data(survival_variable, data_predict_all[[i]],
                         "the survival-time variable used to fit survival_fit_all",
                         sprintf("data_predict_all[[%d]]", i))
  }

  ## at risk sample
  data_predict_all = subset_at_risk(data_predict_all, survival_variable, prediction_time)

  upper_bound = 2 * max(data_predict_all[[1]][survival_variable])

  infinity_grid <- prepare_infinity_grid(data_predict_all, long_fit_all, survival_fit_all,
                                          prediction_time, upper_bound, bandcount2)
  predict.time.infinity = infinity_grid$predict.time.infinity
  predict.time.infinity.1 = infinity_grid$predict.time.infinity.1
  S_T_all_infinity = infinity_grid$S_T_all_infinity

  risk.prob.0 = risk.prob.1 = NULL


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

    # Y_density is filled row-by-row into a pre-allocated matrix rather than
    # grown with rbind() inside the loop: rbind()-in-a-loop reallocates and
    # copies the whole growing matrix on every iteration (O(bandcount3^2)
    # copies total), which becomes non-negligible at the large end of
    # bandcount3's auto-tuned range. The matrix is allocated on the first
    # iteration once the per-patient row length is known, so this makes no
    # assumption about that length elsewhere.
    Y_density = NULL
    for(Y_i in seq_len(length(Y_all))){
      T.surv.predict.0 = t(f_y_D_all_predict[[1]][[Y_i]] * D_T_all_infinity[[1]] * S_T_all_infinity)
      T.surv.predict.1 = t(f_y_D_all_predict[[2]][[Y_i]] * D_T_all_infinity[[2]] * S_T_all_infinity)
      T.surv.infinity.0 = t(f_y_D_all_infinity[[1]] * D_T_all_infinity[[1]] * S_T_all_infinity)
      T.surv.infinity.1 = t(f_y_D_all_infinity[[2]] * D_T_all_infinity[[2]] * S_T_all_infinity)

      risk.prob.1 = clamp_risk_prob(rowSums(T.surv.predict.1), rowSums(T.surv.infinity.1 + T.surv.infinity.0))
      risk.prob.0 = clamp_risk_prob(rowSums(T.surv.predict.0), rowSums(T.surv.infinity.1 + T.surv.infinity.0))
      Y_density_row = risk.prob.1 + risk.prob.0
      if(is.null(Y_density)) Y_density = matrix(NA_real_, length(Y_all), length(Y_density_row))
      Y_density[Y_i, ] = Y_density_row
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
    
    # See the competing-risk branch above for why Y_density is filled into a
    # pre-allocated matrix instead of grown with rbind() in the loop.
    Y_density = NULL
    for(Y_i in seq_len(length(Y_all))){

      T.surv.predict.0 = t(f_y_D_all_predict[[1]][[Y_i]] * S_T_all_infinity)
      T.surv.infinity.0 = t(f_y_D_all_infinity[[1]] * S_T_all_infinity)

      risk.prob.0 = clamp_risk_prob(rowSums(T.surv.predict.0), rowSums(T.surv.infinity.0 + 1e-20))
      if(is.null(Y_density)) Y_density = matrix(NA_real_, length(Y_all), length(risk.prob.0))
      Y_density[Y_i, ] = risk.prob.0

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

