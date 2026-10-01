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
  
  coxph_fit = survival_fit_all$coxph_fit
  # survival data frame
  num <- as.character(nlme::splitFormula(long_fit_all$long_sub_random[[1]], "|")[[2]])[2]
  data.surv =  data_predict_all[[1]][!duplicated(data_predict_all[[1]][num]), ]
  
  ## baseline hazard, extrapolated (per stratum, for a stratified Cox model)
  cum_basehaz_all = basehaz(coxph_fit, centered = FALSE)
  strata_vars <- survival::untangle.specials(stats::terms(coxph_fit), "strata")$vars

  ## covariates * parameter matrix
  ### predict(type = "lp", reference = "zero") is exactly the uncentered
  ### sum(coefficients * covariates) used before, but built by survival
  ### itself: it needs only the right-hand side (the time/status columns of
  ### Surv(...) are unknown for a new patient), handles factors and
  ### strata() terms, and with na.pass a missing covariate gives that
  ### patient NA instead of dropping the row and misaligning later patients.
  covariate_para_matrix = c(stats::predict(coxph_fit, newdata = data.surv, type = "lp",
                                           reference = "zero", na.action = stats::na.pass))

  if (length(strata_vars) == 0) {
    cum_hazard = cumulative_baseline_at(cum_basehaz_all, l_i)
    surv_med = exp(-cum_hazard %*% t(exp(covariate_para_matrix)))
  } else {
    ### each patient uses the baseline hazard of their own stratum
    rhs_terms <- stats::delete.response(stats::terms(coxph_fit))
    mf_surv <- model.frame(rhs_terms, data.surv, na.action = stats::na.pass, xlev = coxph_fit$xlevels)
    patient_strata <- as.character(mf_surv[[strata_vars]])
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
#' helpers: for each element of \code{l_i}, the nearest tabulated value of a
#' \code{basehaz()} table, or -- past the table's last time -- the
#' least-squares line through the whole table (a linear extrapolation of the
#' cumulative hazard, i.e. a constant hazard). The line used to be tabulated
#' on a fixed 0.005 grid out to twice the integration upper limit and then
#' looked up, which assumed time was measured in years (with time in days
#' the table ran to millions of rows) and failed with "wrong sign in 'by'"
#' when the upper limit was below the last training time; it is now
#' evaluated directly.
#'
#' @param cum_basehaz A data frame with columns \code{hazard} and \code{time},
#'   sorted by \code{time}.
#' @param l_i Time points.
#' @return A numeric vector, one cumulative hazard per element of \code{l_i}.
#' @keywords internal
cumulative_baseline_at <- function(cum_basehaz, l_i) {
  cum_basehaz = cum_basehaz[, c("hazard", "time")]
  out <- cum_basehaz$hazard[nearest_index(cum_basehaz$time, l_i)]
  beyond <- l_i > max(cum_basehaz$time)
  if (any(beyond)) {
    line <- stats::coef(lm(hazard ~ time, cum_basehaz))
    out[beyond] <- line[1] + line[2] * l_i[beyond]
  }
  out
}

#' Index of the nearest tabulated time
#'
#' @description For each element of \code{query}, the index of the closest
#' element of the non-decreasing vector \code{x}, choosing the earliest index
#' on ties -- the same answer as \code{which.min(abs(x - q))}, but by binary
#' search (\code{findInterval()}) instead of building a
#' \code{length(query) x length(x)} distance matrix, whose cost grew with
#' the integration upper bound.
#'
#' @param x A non-decreasing numeric vector.
#' @param query Numeric values to look up.
#' @return An integer vector of indices into \code{x}.
#' @keywords internal
nearest_index <- function(x, query) {
  n <- length(x)
  lower <- findInterval(query, x)  # x[lower] <= q < x[lower + 1]
  upper <- pmin(lower + 1L, n)
  lower_ok <- lower >= 1L
  d_lower <- ifelse(lower_ok, query - x[pmax(lower, 1L)], Inf)
  d_upper <- ifelse(lower < n, x[upper] - query, Inf)
  chosen <- ifelse(d_lower <= d_upper, pmax(lower, 1L), upper)
  ### map to the first of any run of equal times (which.min's tie rule)
  match(x[chosen], x)
}
