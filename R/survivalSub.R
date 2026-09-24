#' Fitting survival sub-model
#' 
#' @param data_survival_fitting Input data containing survival outcomes and baseline covariates.
#' @param form_marginal_surv Survival input formats.
#' @param form_conditional_cr Competing risks input formats.
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
  coxph_fit = coxph(form_marginal_surv, data = data_survival_fitting, ties = "breslow")
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
