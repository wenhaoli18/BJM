#' @importFrom ggplot2 ggplot aes geom_point geom_line geom_smooth
#'   geom_text geom_vline geom_hline geom_ribbon scale_color_manual
#'   scale_y_continuous scale_x_continuous sec_axis guide_legend
#'   ylab xlab theme_bw theme element_blank
#' @importFrom stats binomial formula glm lm model.frame model.matrix
#'   model.response na.omit predict terms
#' @importFrom survival coxph basehaz Surv strata
#' @importFrom nlme lme lmeControl splitFormula getVarCov fixef
#' @importFrom Matrix bdiag
#' @importFrom mvtnorm dmvnorm
NULL

utils::globalVariables(c(
  "time", "longitudinal", "probEvent", "probType1", "probType2",
  "predMode", "predQuan1", "predQuan2", "predQuan3", "predQuan4",
  "predQuan5", "predQuan6", "predQuan7", "predQuan8", "predQuan9",
  "Plot_p1"
))