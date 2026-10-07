#' Accessors for the marginal survival sub-model
#'
#' @description Everything the prediction code needs from the marginal
#' survival model of \code{survivalSub()} -- the name of the event-time
#' variable, each patient's linear predictor, the baseline cumulative hazard
#' (per stratum), and the last time it is tabulated at -- goes through these
#' helpers, so the same prediction code serves both a right-censored fit (a
#' \code{\link[survival]{coxph}} model, \code{survival_fit_all$coxph_fit})
#' and an interval-censored fit (a piecewise-constant-hazard proportional
#' hazards model, \code{survival_fit_all$ic_fit}; see \code{icphFit()}).
#' For a right-censored fit each accessor returns exactly what the code used
#' to compute from \code{coxph_fit} inline.
#'
#' @param survival_fit_all Output of \code{survivalSub()}.
#' @param newdata One row per patient.
#' @name survival_accessors
#' @keywords internal
NULL

#' @describeIn survival_accessors \code{TRUE} for an interval-censored fit.
is_interval_censored <- function(survival_fit_all) {
  !is.null(survival_fit_all$ic_fit)
}

#' @describeIn survival_accessors Name of the event-time variable that the
#'   longitudinal sub-model conditions on (for an interval-censored fit, the
#'   \code{event_time} column given to \code{survivalSub()}).
survival_time_variable <- function(survival_fit_all) {
  if (is_interval_censored(survival_fit_all)) return(survival_fit_all$ic_fit$event_time)
  as.character(formula(survival_fit_all$coxph_fit)[[2]])[2]
}

#' @describeIn survival_accessors Uncentered linear predictor, \code{NA} for
#'   a patient with a missing covariate.
survival_lp <- function(survival_fit_all, newdata) {
  if (is_interval_censored(survival_fit_all)) return(icph_lp(survival_fit_all$ic_fit, newdata))
  c(stats::predict(survival_fit_all$coxph_fit, newdata = newdata, type = "lp",
                   reference = "zero", na.action = stats::na.pass))
}

#' @describeIn survival_accessors Baseline cumulative hazard: a data frame
#'   with columns \code{hazard} and \code{time} (and \code{strata} for a
#'   stratified Cox model).
survival_cum_basehaz <- function(survival_fit_all) {
  if (is_interval_censored(survival_fit_all)) return(survival_fit_all$ic_fit$cum_basehaz)
  basehaz(survival_fit_all$coxph_fit, centered = FALSE)
}

#' @describeIn survival_accessors Each patient's stratum (as character), or
#'   \code{NULL} when the model is not stratified.
survival_patient_strata <- function(survival_fit_all, newdata) {
  if (is_interval_censored(survival_fit_all)) return(NULL)
  coxph_fit <- survival_fit_all$coxph_fit
  strata_vars <- survival::untangle.specials(stats::terms(coxph_fit), "strata")$vars
  if (length(strata_vars) == 0) return(NULL)
  rhs_terms <- stats::delete.response(stats::terms(coxph_fit))
  mf_surv <- model.frame(rhs_terms, newdata, na.action = stats::na.pass, xlev = coxph_fit$xlevels)
  as.character(mf_surv[[strata_vars]])
}

#' @describeIn survival_accessors Largest follow-up time in the data the
#'   survival model was fit on (for interval-censored data, the largest
#'   finite interval endpoint).
survival_max_followup <- function(survival_fit_all) {
  if (is_interval_censored(survival_fit_all)) return(survival_fit_all$ic_fit$max_followup)
  surv_y <- survival_fit_all$coxph_fit$y
  max(surv_y[, ncol(surv_y) - 1])
}
