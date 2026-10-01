#' Dynamic prediction function for future event risk
#'
#' @description
#' Combines a fitted longitudinal sub-model (\code{\link{longitudinalSub}})
#' and survival sub-model (\code{\link{survivalSub}}) into a backward joint
#' model prediction: for each at-risk subject at \code{prediction_time}, the
#' subject's observed longitudinal history up to that time is used to update
#' their individual random effects (empirical Bayes), which are then
#' integrated against the survival sub-model's hazard to give the predicted
#' probability of experiencing the event within
#' \code{(prediction_time, prediction_time + horizon]}, conditional on being
#' event-free at \code{prediction_time}. The integrals in both the numerator
#' (event within the window) and denominator (survival to
#' \code{prediction_time}, used for normalization) are evaluated on
#' numerical grids whose resolution is controlled by \code{bandcount1}/
#' \code{bandcount2}; see Details.
#'
#' The prediction is conditional on the longitudinal history observed up to
#' \code{prediction_time}: rows of \code{data_predict_all} whose
#' \code{time_variable} is later than \code{prediction_time} are dropped,
#' with a warning, before predicting.
#'
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
#' @param prediction_time Time used to make the prediction.
#' @param horizon Prediction horizon: the length of the window after
#' \code{prediction_time}, \code{>= 0}.
#' @param time_variable The name of time variable in linear mixed model.
#' @param survival_variable_all The name of the transformed time-to-event outcomes variable.
#' @param survival_trans_function The transformation function used for time-to-event outcomes, 
#' in the order of \code{survival_variable_all}.
#' @param bandcount1 The number of grid points spanning the prediction window,
#' from \code{prediction_time} to \code{prediction_time + horizon} (the
#' numerator of the risk probability). Larger values give a more accurate but
#' slower estimate. Defaults to \code{"auto"} (see Details).
#' @param bandcount2 The number of grid points spanning
#' \code{[prediction_time, upper_bound]}, where \code{upper_bound} is set
#' internally as the earliest time by which every at-risk patient's
#' model-based probability of still being event-free (given event-free at
#' \code{prediction_time}) has dropped below \code{1e-4}; this approximates
#' integrating out to infinity for the denominator that normalizes the risk probability. A wider follow-up range
#' needs a larger \code{bandcount2} to keep the grid spacing comparable.
#' Defaults to \code{"auto"} (see Details).
#'
#' @details
#' There is no universal correct value for \code{bandcount1}/\code{bandcount2}:
#' as a practical check, double both and confirm the resulting risk
#' probabilities barely change; if they do, keep doubling. By default
#' (\code{bandcount1 = "auto"}, \code{bandcount2 = "auto"}), this doubling
#' check is done for you: starting from small built-in values, both are
#' doubled together, and the result is compared to the previous round,
#' until the largest relative change in the risk probabilities drops below
#' 1%, or 2 doublings have been tried (so at most 3 calls' worth of work).
#' If it still has not converged by then, a warning reports this and the
#' result at the largest value tried is returned anyway (not an error), so
#' this never silently loops for an unbounded amount of time. Pass an
#' explicit number for either argument to skip auto-tuning it and use a
#' fixed value instead (as in previous package versions), or call
#' \code{checkBandcountConvergence()} directly for more control over the
#' tolerance and doubling count. See also \code{vignette("BJM-intro",
#' package = "BJM")} for a worked example.
#'
#' The denominator integrates over every event time after
#' \code{prediction_time}, including times beyond the last follow-up time
#' in the data \code{survivalSub()} was fit on. There, the baseline
#' cumulative hazard is extrapolated linearly (a constant hazard), and the
#' longitudinal sub-model's mean is evaluated at event times it was never
#' fit on. When much of an at-risk patient's survival probability lies
#' beyond the last follow-up time (e.g. a prediction late in follow-up),
#' the prediction depends on this extrapolation; in \code{pbc3} it moved
#' the risks checked by less than one percentage point.
#'
#' For a fit with an ordinal biomarker, each ordinal measurement contributes
#' a multivariate normal probability computed by Monte Carlo
#' (\code{mvtnorm::pmvnorm()}), so repeated calls differ slightly (around
#' the fifth significant digit); call \code{set.seed()} first for exactly
#' reproducible results.
#' 
#' @return An object of class \code{"predictRisk.BJM"}, a named list with elements:
#' \describe{
#'   \item{risk_prob_1}{A vector of dynamically predicted probabilities, one per at-risk
#'   patient, of experiencing the (first) event within the prediction horizon, named by
#'   patient id. All \code{0} when \code{horizon = 0}.}
#'   \item{risk_prob_2}{When \code{survival_fit_all} was fit with competing risks, a vector
#'   of dynamically predicted probabilities, one per at-risk patient, of experiencing the
#'   competing event within the prediction horizon, named by patient id (all \code{0} when
#'   \code{horizon = 0}). \code{NULL} when there is no competing risk.}
#' }
#' Only patients still at risk at \code{prediction_time} are predicted: a patient whose
#' recorded survival time is before \code{prediction_time} is left out (one with a missing
#' survival time is kept). The names show which patients each value belongs to.
#' The numerical integration grids are shared by all patients predicted in one
#' call, so a patient's value can differ slightly (within the grid's
#' discretization error, which shrinks as \code{bandcount1}/\code{bandcount2}
#' grow) depending on which other patients are predicted alongside it.
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
#' # predict risk probability
#' risk.prob = predictRisk(data_predict_all, long_fit_all, survival_fit_all,
#'                               prediction_time = 3,
#'                               horizon = 3, time_variable = "year",
#'                               survival_variable_all, survival_trans_function,
#'                               bandcount1 = 10, bandcount2 = 10)
#'
#' # poly() in its default orthogonal mode, splines::ns()/bs(), and factor()
#' # in a continuous biomarker's long_sub_fixed produce correct dynamic
#' # predictions (see ?longitudinalSub) -- including
#' # when a patient has only a single longitudinal observation to condition
#' # on -- because the basis/contrasts fit on the full training data are
#' # cached (via each biomarker's terms object and long_fit_all$xlevels) and
#' # reused here at prediction time, instead of being recomputed from that
#' # patient's small per-prediction-time slice:
#' long_sub_fixed_nonlinear = list(
#'   "long1" = serBilir ~ poly(year, 2) + age + sex + years,
#'   "long2" = albumin ~ year + age + sex + years)
#' long_sub_random_nonlinear = list("long1" = ~ year | id, "long2" = ~ year | id)
#' long_fit_nonlinear = longitudinalSub(pbc3[pbc3$status3 == 1, ],
#'                                      long_sub_fixed_nonlinear,
#'                                      long_sub_random_nonlinear)
#'
#' data_predict_normal = data.raw.predict.1[data.raw.predict.1$year <= 3, ]
#' data_predict_sparse = data.raw.predict.1[1, ]
#'
#' risk.prob.normal = predictRisk(data_predict_normal, long_fit_nonlinear,
#'                                      survival_fit_all, prediction_time = 3,
#'                                      horizon = 3, time_variable = "year",
#'                                      survival_variable_all, survival_trans_function,
#'                                      bandcount1 = 10, bandcount2 = 10)
#' risk.prob.sparse = predictRisk(data_predict_sparse, long_fit_nonlinear,
#'                                      survival_fit_all, prediction_time = 3,
#'                                      horizon = 3, time_variable = "year",
#'                                      survival_variable_all, survival_trans_function,
#'                                      bandcount1 = 10, bandcount2 = 10)
#' # both give a sane, non-degenerate risk_prob_1 (not 0, no error)
#' risk.prob.normal$risk_prob_1
#' risk.prob.sparse$risk_prob_1
#'
#' }
#'
#' @export
predictRisk = function(data_predict_all, long_fit_all, survival_fit_all,
                             prediction_time, horizon, time_variable,
                             survival_variable_all, survival_trans_function,
                             bandcount1 = "auto", bandcount2 = "auto"){

  assert_class(long_fit_all, "longitudinalSub.BJM", "long_fit_all", "longitudinalSub")
  assert_class(survival_fit_all, "survivalSub.BJM", "survival_fit_all", "survivalSub")
  assert_data_list(data_predict_all, "data_predict_all", length(long_fit_all$lfit), allow_bare_df = TRUE)
  if (!is.list(data_predict_all) || is.data.frame(data_predict_all)) {
    data_predict_all <- rep(list(data_predict_all), each = length(long_fit_all$lfit))
  }
  assert_scalar_numeric(prediction_time, "prediction_time")
  assert_scalar_numeric(horizon, "horizon")
  if (horizon < 0) {
    stop(sprintf("`horizon` must be >= 0 (the length of the window after prediction_time), not %g.", horizon),
         call. = FALSE)
  }
  assert_string(time_variable, "time_variable")
  assert_bandcount(bandcount1, "bandcount1")
  assert_bandcount(bandcount2, "bandcount2")
  assert_survival_trans(survival_variable_all, survival_trans_function, probe_value = prediction_time)
  data_predict_all <- drop_after_prediction_time(data_predict_all, time_variable, prediction_time)
  data_predict_all <- align_ordinal_levels(data_predict_all, long_fit_all)
  data_predict_all <- drop_missing_longitudinal(data_predict_all, long_fit_all,
                                               as.character(formula(survival_fit_all$coxph_fit)[[2]])[2],
                                               survival_variable_all)

  # bandcount1/bandcount2 = "auto" (the default): resolve them by doubling
  # from their built-in starting values until the returned risk
  # probabilities stabilize (see auto_tune_bandcount()), instead of
  # requiring the caller to pick a value. `call_args` is captured here,
  # right after validation and before any other local variables exist
  # (other than `call_args` itself), so it is exactly this call's
  # (possibly-defaulted) arguments -- `auto_names` is deliberately computed
  # afterwards so it is not swept up into `call_args` too.
  call_args <- as.list(environment())
  auto_names <- c("bandcount1", "bandcount2")[c(identical(bandcount1, "auto"), identical(bandcount2, "auto"))]
  if (length(auto_names) > 0) {
    return(auto_tune_bandcount(predictRisk, call_args, auto_names)$result)
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
  patient_ids <- prediction_patient_ids(data_predict_all, long_fit_all)
  has_cr <- length(survival_fit_all$form_conditional_cr) != 0

  #### handle horizon = 0 edge case: probability of event in zero-length window is 0
  if(horizon == 0){
    zeros <- stats::setNames(rep(0, length(patient_ids)), patient_ids)
    out <- list(risk_prob_1 = zeros, risk_prob_2 = if (has_cr) zeros else NULL)
    class(out) <- "predictRisk.BJM"
    return(out)
  }

  upper_bound = integration_upper_bound(data_predict_all, long_fit_all, survival_fit_all,
                                        prediction_time, min_upper = prediction_time + horizon)

  #### time frame used to do the integral
  ### bandcount1 intervals tiling exactly (prediction_time, prediction_time +
  ### horizon], each evaluated at its midpoint. The previous grid's edges ran
  ### from half an interval before prediction_time to half an interval after
  ### the horizon, so the window it integrated over was horizon + 1 interval
  ### long, inflating the risk by roughly 1 / bandcount1 (about 10% at
  ### bandcount1 = 10) -- the main reason "auto" bandcount tuning kept
  ### reporting non-convergence.
  predict.time.horizon.1 = seq(prediction_time, prediction_time + horizon, length.out = bandcount1 + 1)
  predict.time.horizon = (predict.time.horizon.1[-1] + predict.time.horizon.1[-(bandcount1 + 1)]) / 2

  infinity_grid <- prepare_infinity_grid(data_predict_all, long_fit_all, survival_fit_all,
                                          prediction_time, upper_bound, bandcount2)
  predict.time.infinity = infinity_grid$predict.time.infinity
  predict.time.infinity.1 = infinity_grid$predict.time.infinity.1
  S_T_all_infinity = infinity_grid$S_T_all_infinity

  risk.prob.0 = risk.prob.1 = NULL
  ### marginal probability T
  S_T_all_predict = marginalT(data_predict_all, long_fit_all, survival_fit_all, l_i = predict.time.horizon.1, upper_bound)

  ### Mixed continuous/ordinal (Gaussian-copula) fits -- see
  ### longitudinalSubCopula() -- need conditionalYTCopula()/
  ### conditionalYDTCopula() instead of conditionalYT()/conditionalYDT(),
  ### since the latter assume every biomarker is continuous (an exact
  ### Gaussian density evaluated at the observed value). All-continuous
  ### fits never set long_fit_all$biomarker_type at all (see
  ### longitudinalSubGaussian()), so this dispatch is a strict no-op for
  ### every fit that predates the copula extension.
  use_copula <- !is.null(long_fit_all$biomarker_type) && any(long_fit_all$biomarker_type == "ordinal")
  conditionalYT_fun <- if (use_copula) conditionalYTCopula else conditionalYT
  conditionalYDT_fun <- if (use_copula) conditionalYDTCopula else conditionalYDT

  #conditional probability D|T, survival_fit_all$form_conditional_cr == form_conditional_cr
  #with competing risk
  if(length(survival_fit_all$form_conditional_cr) != 0){
    #conditional probability D|T
    D_T_all_predict = conditionalDT(data_predict_all, long_fit_all, survival_fit_all,
                                    l_i = predict.time.horizon)
    D_T_all_infinity = conditionalDT(data_predict_all, long_fit_all, survival_fit_all,
                                     l_i = predict.time.infinity)
    #conditional probability Y|D,T
    f_y_D_all_predict = conditionalYDT_fun(data_predict_all, long_fit_all, survival_fit_all,
                                       l_i = predict.time.horizon, survival_variable,
                                       time_variable, survival_variable_all,
                                       survival_trans_function)
    f_y_D_all_infinity = conditionalYDT_fun(data_predict_all, long_fit_all, survival_fit_all,
                                        l_i = predict.time.infinity, survival_variable,
                                        time_variable, survival_variable_all,
                                        survival_trans_function)
    ### the conditional densities are log densities: shift each patient's
    ### by that patient's largest value before exponentiating, so they
    ### neither underflow to 0 nor overflow (the shift cancels in the ratio)
    shift = patient_log_shift(f_y_D_all_predict[[1]], f_y_D_all_predict[[2]],
                              f_y_D_all_infinity[[1]], f_y_D_all_infinity[[2]])
    T.surv.predict.0 = t(exp_shifted(f_y_D_all_predict[[1]], shift) * D_T_all_predict[[1]] * S_T_all_predict)
    T.surv.infinity.0 = t(exp_shifted(f_y_D_all_infinity[[1]], shift) * D_T_all_infinity[[1]] * S_T_all_infinity)
    T.surv.predict.1 = t(exp_shifted(f_y_D_all_predict[[2]], shift) * D_T_all_predict[[2]] * S_T_all_predict)
    T.surv.infinity.1 = t(exp_shifted(f_y_D_all_infinity[[2]], shift) * D_T_all_infinity[[2]] * S_T_all_infinity)

    risk.prob.0 = clamp_risk_prob(rowSums(T.surv.predict.0), rowSums(T.surv.infinity.0 + T.surv.infinity.1))
    risk.prob.1 = clamp_risk_prob(rowSums(T.surv.predict.1), rowSums(T.surv.infinity.0 + T.surv.infinity.1))

  }else{
    #without competing risk
    #conditional probability Y|T
    f_y_D_all_predict = conditionalYT_fun(data_predict_all, long_fit_all, l_i = predict.time.horizon,
                                      survival_variable, time_variable, survival_variable_all, survival_trans_function)
    f_y_D_all_infinity = conditionalYT_fun(data_predict_all, long_fit_all, l_i = predict.time.infinity,
                                       survival_variable, time_variable, survival_variable_all, survival_trans_function)

    ### see the competing-risk branch: shifted log densities. This also
    ### replaces a "+ 1e-20" added to the denominator to avoid dividing by
    ### 0, which, being absolute, swamped the denominator whenever the
    ### densities were small (biomarkers on a large scale or many
    ### observations) and drove the risk towards 0.
    shift = patient_log_shift(f_y_D_all_predict[[1]], f_y_D_all_infinity[[1]])
    T.surv.predict.0 = t(exp_shifted(f_y_D_all_predict[[1]], shift) * S_T_all_predict)
    T.surv.infinity.0 = t(exp_shifted(f_y_D_all_infinity[[1]], shift) * S_T_all_infinity)

    risk.prob.0 = clamp_risk_prob(rowSums(T.surv.predict.0), rowSums(T.surv.infinity.0))
  }
  
  names(risk.prob.0) <- patient_ids
  if (!is.null(risk.prob.1)) names(risk.prob.1) <- patient_ids
  out <- list(risk_prob_1 = risk.prob.0, risk_prob_2 = risk.prob.1)
  class(out) <- "predictRisk.BJM"
  return(out)
}


