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
#' @return An object of class \code{"survivalSub.BJM"}, a named list with elements:
#' \describe{
#'   \item{coxph_fit}{The fitted \code{\link[survival]{coxph}} marginal survival model.}
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
survivalSub = function(data_survival_fitting, form_marginal_surv, form_conditional_cr){

  assert_data_frame(data_survival_fitting, "data_survival_fitting")
  if (!inherits(form_marginal_surv, "formula")) {
    stop("`form_marginal_surv` must be a formula, e.g. Surv(time, status) ~ covariates.", call. = FALSE)
  }
  assert_vars_in_data(all.vars(form_marginal_surv), data_survival_fitting,
                       "form_marginal_surv", "data_survival_fitting")
  if (length(form_conditional_cr) != 0) {
    if (!inherits(form_conditional_cr, "formula")) {
      stop("`form_conditional_cr` must be a formula or NULL.", call. = FALSE)
    }
    assert_vars_in_data(all.vars(form_conditional_cr), data_survival_fitting,
                         "form_conditional_cr", "data_survival_fitting")
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
