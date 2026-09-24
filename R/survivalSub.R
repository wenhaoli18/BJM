#' Fitting survival sub-model
#' 
#' @param data.survival.fitting Input data containing survival outcomes and baseline covariates.
#' @param formMarginalSurv Survival input formats.
#' @param formConditionalCR Competing risks input formats.
#' @return An object of class \code{"survivalSub.BJM"}, a named list with elements:
#' \describe{
#'   \item{coxph_fit}{The fitted \code{\link[survival]{coxph}} marginal survival model.}
#'   \item{formMarginalSurv}{The \code{formMarginalSurv} formula, as supplied.}
#'   \item{glm_fit}{The fitted \code{\link[stats]{glm}} competing-risks (event type) model,
#'   or \code{NULL} if \code{formConditionalCR} was not supplied.}
#'   \item{formConditionalCR}{The \code{formConditionalCR} formula, as supplied (or \code{NULL}).}
#' }
#'
#' @examples 
#' 
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
#' @export
survivalSub = function(data.survival.fitting, formMarginalSurv, formConditionalCR){
  
  ### fit cox weibull model
  coxph_fit = coxph(formMarginalSurv, data = data.survival.fitting, ties = "breslow")
  ### censoring indicator name
  censor_variable = as.character(formula(coxph_fit)[[2]])[3]
  
  data.glm = data.survival.fitting[data.survival.fitting[censor_variable] != 0, ]
  ##with or without competing risk
  if(length(formConditionalCR) != 0){
    ## competing risks outcomes
    ## Vertical model
    glm_fit = glm(formConditionalCR, data.glm, family = binomial)
  }else{
    glm_fit = NULL
  }
  
  out <- list(coxph_fit = coxph_fit, formMarginalSurv = formMarginalSurv,
              glm_fit = glm_fit, formConditionalCR = formConditionalCR)
  class(out) <- "survivalSub.BJM"
  return(out)
}
