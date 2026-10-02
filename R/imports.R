#' @importFrom ggplot2 ggplot aes geom_point geom_line geom_smooth
#'   geom_text geom_vline geom_hline geom_ribbon scale_color_manual
#'   scale_y_continuous scale_x_continuous sec_axis guide_legend
#'   ylab xlab theme_bw theme element_blank geom_step labs guides
#'   facet_wrap geom_qq geom_qq_line geom_tile geom_col geom_errorbar
#'   scale_fill_gradient2 scale_x_log10 element_text geom_abline facet_grid
#' @importFrom graphics plot
#' @importFrom stats binomial formula glm lm model.frame model.matrix
#'   model.response na.omit predict terms .getXlevels as.formula sd rnorm
#' @importFrom survival coxph basehaz Surv strata survfit
#' @importFrom nlme lme lmeControl splitFormula getVarCov fixef
#' @importFrom Matrix bdiag
#' @importFrom mvtnorm dmvnorm
NULL

utils::globalVariables(c(
  "time", "longitudinal", "probEvent", "probType1", "probType2",
  "predMode", "predQuan1", "predQuan2", "predQuan3", "predQuan4",
  "predQuan5", "predQuan6", "predQuan7", "predQuan8", "predQuan9",
  "Plot_p1",
  ### spaghettiPlot()/cifPlot() aes() columns
  "x", "y", "id", "group", "cif", "lower", "upper", "event", "curve", "hue",
  ### plot methods for longitudinalSub/survivalSub/predictLongitudinal results
  "fitted", "resid", "biomarker", "panel", "value", "row", "col", "hazard", "strata",
  "estimate", "term", "model", "category", "density", "subject",
  ### performancePlot() aes() columns
  "landmark", "measure", "cause", "predicted", "observed",
  ### `self` is torch::nn_module()'s implicit binding to the module
  ### instance inside `initialize`/`encode`/`decode`, supplied by torch at
  ### call time -- not a real undefined global, just invisible to R CMD
  ### check's static analysis.
  "self"
))