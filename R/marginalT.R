#' Marginal distribution of T
#' 
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
#' @param l_i A vector of time points to calculate the conditional probability.
#' @param upper_bound Upper limit for integration. To manage the hazard function, 
#' extrapolation is required. The \code{upper_bound} parameter specifies the upper time 
#' points at which extrapolation is performed.
#' 
#' @return Marginal density probability of the survival variable T is represented as a 
#' probability matrix. In this matrix, the rows (l_i) are aligned with specific 
#' time points, while the columns correspond to individual patients. 
#' Each entry in the matrix denotes the marginal density probability of 
#' survival for a given patient at a particular time point.
#' @keywords internal
marginalT = function(data_predict_all, long_fit_all, survival_fit_all, l_i, upper_bound){
  
  # survival data frame
  num <- as.character(nlme::splitFormula(long_fit_all$long_sub_random[[1]], "|")[[2]])[2]
  data.surv =  data_predict_all[[1]][!duplicated(data_predict_all[[1]][num]), ]
  
  ## baseline hazard, extrapolated (per stratum, for a stratified Cox model)
  cum_basehaz_all = survival_cum_basehaz(survival_fit_all)

  ## covariates * parameter matrix
  ### uncentered linear predictor (see survival_lp()): needs only the
  ### right-hand side (the outcome columns are unknown for a new patient),
  ### handles factors and strata() terms, and a missing covariate gives that
  ### patient NA instead of dropping the row and misaligning later patients.
  covariate_para_matrix = survival_lp(survival_fit_all, data.surv)
  patient_strata <- survival_patient_strata(survival_fit_all, data.surv)

  if (is.null(patient_strata)) {
    cum_hazard = cumulative_baseline_at(cum_basehaz_all, l_i)
    surv_med = exp(-cum_hazard %*% t(exp(covariate_para_matrix)))
  } else {
    ### each patient uses the baseline hazard of their own stratum
    unknown <- setdiff(stats::na.omit(patient_strata), as.character(cum_basehaz_all$strata))
    if (length(unknown) > 0) {
      stop(sprintf("Stratum %s of the survival sub-model has no baseline hazard (no such subjects in the data survivalSub() was fit on).",
                   paste(unknown, collapse = ", ")), call. = FALSE)
    }
    surv_med = matrix(NA_real_, length(l_i), length(patient_strata))
    for (s in unique(stats::na.omit(patient_strata))) {
      cum_hazard = cumulative_baseline_at(cum_basehaz_all[cum_basehaz_all$strata == s, c("hazard", "time")],
                                          l_i)
      cols = which(patient_strata == s)
      surv_med[, cols] = exp(-cum_hazard %*% t(exp(covariate_para_matrix[cols])))
    }
  }
  ### density
  S_T_all = - diff(surv_med)
  return(S_T_all)
}

#' Baseline cumulative hazard at a set of time points
#'
#' @description Helper for \code{marginalT} and the integration-grid
#' helpers: for each element of \code{l_i}, the Breslow cumulative hazard
#' as the right-continuous step function it is -- the value at the last
#' tabulated time \code{<= l_i}, or 0 before the first. Past the table's
#' last time it continues from the last tabulated value with the slope of
#' the least-squares line through the whole table (a linear extrapolation
#' of the cumulative hazard, i.e. a constant hazard). It used to switch to
#' that line itself, which does not pass through the last tabulated value:
#' the cumulative hazard jumped there, and when it jumped down (a convex
#' cumulative hazard, i.e. an increasing hazard) the interval spanning the
#' last time got a negative probability mass. The line
#' used to be tabulated
#' on a fixed 0.005 grid out to twice the integration upper limit and then
#' looked up, which assumed time was measured in years (with time in days
#' the table ran to millions of rows) and failed with "wrong sign in 'by'"
#' when the upper limit was below the last training time; it is now
#' evaluated directly.
#' The table used to be read at the \emph{nearest} tabulated time, which
#' for a point just before an event time took that event's jump early.
#'
#' @param cum_basehaz A data frame with columns \code{hazard} and \code{time}.
#' @param l_i Time points.
#' @return A numeric vector, one cumulative hazard per element of \code{l_i}.
#' @keywords internal
cumulative_baseline_at <- function(cum_basehaz, l_i) {
  cum_basehaz = cum_basehaz[, c("hazard", "time")]
  cum_basehaz = cum_basehaz[order(cum_basehaz$time), ]
  at_or_before <- findInterval(l_i, cum_basehaz$time)
  out <- c(0, cum_basehaz$hazard)[at_or_before + 1]
  beyond <- l_i > max(cum_basehaz$time)
  if (any(beyond)) {
    last <- nrow(cum_basehaz)
    out[beyond] <- cum_basehaz$hazard[last] +
      extrapolation_slope(cum_basehaz) * (l_i[beyond] - cum_basehaz$time[last])
  }
  out
}

#' Slope of the extrapolated baseline cumulative hazard
#'
#' @description The constant hazard \code{cumulative_baseline_at()} uses
#' past the last tabulated time: the slope of the least-squares line through
#' the whole table, or 0 if that slope is negative or undefined (a table with
#' a single time), so the cumulative hazard never decreases.
#'
#' @param cum_basehaz A data frame with columns \code{hazard} and \code{time},
#'   ordered by time.
#' @return A single non-negative number.
#' @keywords internal
extrapolation_slope <- function(cum_basehaz) {
  slope <- if (nrow(cum_basehaz) < 2) NA_real_ else stats::coef(lm(hazard ~ time, cum_basehaz))[[2]]
  if (is.na(slope) || slope < 0) 0 else slope
}
