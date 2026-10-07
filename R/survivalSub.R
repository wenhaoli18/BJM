#' Fit the survival sub-model
#'
#' @description
#' Fits the survival building block used by \code{\link{predictRisk}}
#' and \code{\link{dynamicPredictionBio}}: a marginal Cox proportional-hazards
#' model (via \code{\link[survival]{coxph}}, with \code{x = TRUE, y = TRUE} so
#' the fit is self-contained for later prediction) for the time-to-event
#' outcome given \code{form_marginal_surv}, plus, optionally, a logistic
#' event-type/competing-risks model (via \code{\link[stats]{glm}} with
#' \code{family = binomial}) given \code{form_conditional_cr}. The
#' competing-risks model is fit only among subjects who experienced an event
#' (i.e. whose censoring indicator is non-zero), predicting which type of
#' event occurred conditional on an event having occurred. Together with
#' \code{\link{longitudinalSub}}, the object returned here forms the pair of
#' sub-models that dynamic prediction is built on.
#'
#' @param data_survival_fitting One row per subject, containing the
#'   time-to-event outcome, censoring/event-type indicator, and any baseline
#'   covariates referenced in \code{form_marginal_surv} or
#'   \code{form_conditional_cr}.
#' @param form_marginal_surv A survival formula, e.g.
#'   \code{Surv(time, status) ~ covariates}, passed to
#'   \code{\link[survival]{coxph}} to fit the marginal event-time model.
#' @param form_conditional_cr An optional formula for the competing-risks
#'   (event-type) model, e.g. \code{event_type ~ covariates}, passed to
#'   \code{\link[stats]{glm}} with \code{family = binomial}. Set to
#'   \code{NULL} when there is only a single event type (no competing
#'   risks).
#' @param event_time Only for an interval-censored outcome: the name of the
#'   event-time column that the longitudinal sub-model's formulas use (the
#'   exact time is unknown, so it is not part of \code{form_marginal_surv}).
#' @param baseline Only for an interval-censored outcome: the baseline
#'   hazard, \code{"spline"} (default; Royston--Parmar spline for the log
#'   cumulative hazard in log time) or \code{"piecewise"}
#'   (piecewise-constant hazard).
#' @param df Only for an interval-censored outcome: spline degrees of
#'   freedom (default 3) or number of pieces (default 6).
#'
#' @details \strong{Interval censoring (experimental).} When
#'   \code{form_marginal_surv} is \code{Surv(L, R, type = "interval2") ~
#'   covariates} -- the event is only known to lie in \code{(L, R]}, with
#'   \code{R = NA} or \code{Inf} for right-censored subjects -- the marginal
#'   model is a proportional hazards model with a smooth spline (or
#'   piecewise-constant) baseline hazard fit by maximum likelihood (see
#'   \code{icphFit()}), stored as \code{ic_fit} in place of
#'   \code{coxph_fit}. Competing risks and \code{strata()} are not yet
#'   supported in this case. Ordinary right-censored \code{Surv(time,
#'   status)} outcomes are fit with \code{coxph} exactly as before.
#'
#'   Events that happened before a subject's first visit (\code{L = 0}) only
#'   say that the event fell between time 0 and that visit. Before the
#'   earliest visit in the data there are no interval endpoints, so the
#'   baseline hazard there is extrapolated, not estimated, and
#'   \code{\link{fitIntervalBJM}} places these subjects' event times
#'   according to that extrapolation. \code{print()} reports how many events
#'   are of this kind: when they are a large share of all events (in
#'   simulations, around a third or more), results can depend noticeably on
#'   the assumed shape of the early hazard.
#' @return An object of class \code{"survivalSub.BJM"}, a named list with elements:
#' \describe{
#'   \item{coxph_fit}{The fitted \code{\link[survival]{coxph}} marginal survival model
#'   (\code{ic_fit}, an interval-censored PH fit, in its place for an
#'   interval-censored outcome).}
#'   \item{form_marginal_surv}{The \code{form_marginal_surv} formula, as supplied.}
#'   \item{glm_fit}{The fitted \code{\link[stats]{glm}} competing-risks (event type) model,
#'   or \code{NULL} if \code{form_conditional_cr} was not supplied.}
#'   \item{form_conditional_cr}{The \code{form_conditional_cr} formula, as supplied (or \code{NULL}).}
#' }
#'
#' @examples 
#' 
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
#' @export
survivalSub = function(data_survival_fitting, form_marginal_surv, form_conditional_cr,
                       event_time = NULL, baseline = c("spline", "piecewise"), df = NULL){

  assert_data_frame(data_survival_fitting, "data_survival_fitting")
  if (!inherits(form_marginal_surv, "formula")) {
    stop("`form_marginal_surv` must be a formula, e.g. Surv(time, status) ~ covariates.", call. = FALSE)
  }
  assert_vars_in_data(all.vars(form_marginal_surv), data_survival_fitting,
                       "form_marginal_surv", "data_survival_fitting")
  if (is_interval_surv(form_marginal_surv, data_survival_fitting)) {
    return(survivalSubInterval(data_survival_fitting, form_marginal_surv, form_conditional_cr,
                               event_time, baseline, df))
  }
  if (length(form_conditional_cr) != 0) {
    if (!inherits(form_conditional_cr, "formula")) {
      stop("`form_conditional_cr` must be a formula or NULL.", call. = FALSE)
    }
    assert_vars_in_data(all.vars(form_conditional_cr), data_survival_fitting,
                         "form_conditional_cr", "data_survival_fitting")
    assert_linear_time_term(form_conditional_cr, all.vars(form_marginal_surv[[2]])[1])
  }

  ### fit cox weibull model
  ### x = TRUE, y = TRUE make the fit self-contained: survfit.coxph()/basehaz()
  ### (called later, at prediction time, possibly in a different environment)
  ### skip re-deriving the model frame from the captured call + its formula's
  ### environment, which would otherwise fail whenever the caller's data
  ### object is not literally named `data_survival_fitting`.
  coxph_fit = coxph(form_marginal_surv, data = data_survival_fitting, ties = "breslow",
                     x = TRUE, y = TRUE)
  ### censoring indicator name
  censor_variable = as.character(formula(coxph_fit)[[2]])[3]
  
  data.glm = data_survival_fitting[data_survival_fitting[censor_variable] != 0, ]
  ##with or without competing risk
  if(length(form_conditional_cr) != 0){
    ## competing risks outcomes
    ## Vertical model
    glm_fit = glm(form_conditional_cr, data.glm, family = binomial)
  }else{
    glm_fit = NULL
  }
  
  out <- list(coxph_fit = coxph_fit, form_marginal_surv = form_marginal_surv,
              glm_fit = glm_fit, form_conditional_cr = form_conditional_cr)
  class(out) <- "survivalSub.BJM"
  return(out)
}

#' Is the outcome of a survival formula interval-censored?
#' @keywords internal
is_interval_surv <- function(form_marginal_surv, data) {
  y <- eval(form_marginal_surv[[2]], data, environment(form_marginal_surv))
  inherits(y, "Surv") && attr(y, "type") %in% c("interval", "interval2")
}

#' Interval-censored branch of \code{survivalSub()}
#' @keywords internal
survivalSubInterval <- function(data_survival_fitting, form_marginal_surv, form_conditional_cr,
                                event_time, baseline, df) {
  baseline <- match.arg(baseline, c("spline", "piecewise"))
  if (is.null(df)) df <- if (baseline == "spline") 3 else 6
  if (length(form_conditional_cr) != 0) {
    stop("Competing risks are not yet supported with an interval-censored outcome; set form_conditional_cr = NULL.",
         call. = FALSE)
  }
  assert_string(event_time, "event_time")
  if (!is.numeric(df) || length(df) != 1 || df < 1 || df != round(df)) {
    stop("`df` must be a single positive whole number.", call. = FALSE)
  }
  ic_fit <- icphFit(form_marginal_surv, data_survival_fitting, event_time, baseline, df)
  out <- list(ic_fit = ic_fit, form_marginal_surv = form_marginal_surv,
              glm_fit = NULL, form_conditional_cr = NULL)
  class(out) <- "survivalSub.BJM"
  out
}
