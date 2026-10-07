#' Drop longitudinal measurements taken after the prediction time
#'
#' @description Shared helper for \code{predictRisk},
#' \code{dynamicPredictionBio}, and \code{dynamicPredictionBioAll}: a dynamic
#' prediction at \code{prediction_time} may only condition on the history
#' observed up to \code{prediction_time}, so rows of \code{data_predict_all}
#' whose \code{time_variable} is later than that are removed, with a warning
#' saying how many. (\code{predictPlot}/\code{riskPlot} already truncate per
#' landmark time before calling these, so they never trigger the warning.)
#' Rows with a missing \code{time_variable} are kept, as before.
#'
#' @return \code{data_predict_all}, filtered per element.
#' @keywords internal
drop_after_prediction_time <- function(data_predict_all, time_variable, prediction_time) {
  n_dropped <- integer(length(data_predict_all))
  for (i in seq_along(data_predict_all)) {
    assert_vars_in_data(time_variable, data_predict_all[[i]],
                         "time_variable", sprintf("data_predict_all[[%d]]", i))
    t <- data_predict_all[[i]][[time_variable]]
    after <- !is.na(t) & t > prediction_time + 1e-8
    n_dropped[i] <- sum(after)
    data_predict_all[[i]] <- data_predict_all[[i]][!after, , drop = FALSE]
  }
  if (any(n_dropped > 0)) {
    warning(sprintf(paste0(
      "Dropped measurements taken after prediction_time = %g (%s > %g) from data_predict_all: %s. ",
      "A dynamic prediction conditions only on the history up to prediction_time; ",
      "subset data_predict_all to %s <= prediction_time to silence this warning."),
      prediction_time, time_variable, prediction_time,
      paste(sprintf("%d row(s) in [[%d]]", n_dropped, seq_along(n_dropped))[n_dropped > 0], collapse = ", "),
      time_variable), call. = FALSE)
  }
  data_predict_all
}

#' Outcome variables of the survival sub-model
#'
#' @description The survival-time variable of \code{survival_fit_all} and,
#' with competing risks, its event-type variable (the response of
#' \code{form_conditional_cr}). Both are unknown for a patient still at risk
#' at \code{prediction_time}; the prediction functions integrate over them.
#'
#' @param survival_fit_all Output of \code{survivalSub()}.
#' @return A character vector of variable names.
#' @keywords internal
outcome_variables <- function(survival_fit_all) {
  out <- survival_time_variable(survival_fit_all)
  if (length(survival_fit_all$form_conditional_cr) != 0) {
    out <- c(out, all.vars(survival_fit_all$form_conditional_cr[[2]]))
  }
  out
}

#' Drop prediction rows with a missing biomarker value or covariate
#'
#' @description Shared helper for \code{predictRisk},
#' \code{dynamicPredictionBio}, and \code{dynamicPredictionBioAll}: a
#' patient's conditional density only involves the measurements actually
#' observed, so for each biomarker, rows of its \code{data_predict_all}
#' element with a missing response, or a missing covariate/time/ID used by
#' that biomarker's \code{long_sub_fixed}/\code{long_sub_random} formulas,
#' are removed. Previously such a row put an \code{NA} into the stacked
#' outcome vector (or misaligned the design matrix) and made that patient's
#' whole prediction \code{NA}. The outcome variables
#' (\code{survival_variable}: the survival time and, with competing risks,
#' the event type -- see \code{outcome_variables()}) and
#' \code{survival_variable_all} are not checked: they are overwritten with
#' each integration grid point / event type before use, and are unknown
#' (typically \code{NA}) for a patient still at risk. The event type used
#' to be checked, so with the event type in \code{long_sub_fixed} every
#' at-risk patient was dropped. A patient left with no rows for some
#' biomarker cannot be predicted; a warning names them, and it is an error
#' if that leaves no patient at all.
#'
#' @return \code{data_predict_all}, filtered per element.
#' @keywords internal
drop_missing_longitudinal <- function(data_predict_all, long_fit_all, survival_variable,
                                       survival_variable_all) {
  id <- as.character(nlme::splitFormula(long_fit_all$long_sub_random[[1]], "|")[[2]])[2]
  all_ids <- unique(unlist(lapply(data_predict_all, function(d) as.character(d[[id]]))))
  unpredictable <- character(0)
  for (i in seq_along(data_predict_all)) {
    vars <- unique(c(all.vars(long_fit_all$long_sub_fixed[[i]]),
                     all.vars(long_fit_all$long_sub_random[[i]]), id))
    vars <- setdiff(vars, c(survival_variable, unlist(survival_variable_all)))
    vars <- intersect(vars, names(data_predict_all[[i]]))
    complete <- rowSums(is.na(data_predict_all[[i]][vars])) == 0
    data_predict_all[[i]] <- data_predict_all[[i]][complete, , drop = FALSE]
    unpredictable <- c(unpredictable,
                       setdiff(all_ids, as.character(data_predict_all[[i]][[id]])))
  }
  unpredictable <- unique(unpredictable)
  if (length(all_ids) > 0 && all(all_ids %in% unpredictable)) {
    stop(paste0(
      "No patient in data_predict_all has a non-missing measurement of every biomarker (and of ",
      "the covariates its long_sub_fixed/long_sub_random formulas use) up to prediction_time, ",
      "so no prediction can be made."), call. = FALSE)
  }
  if (length(unpredictable) > 0) {
    warning(sprintf(paste0(
      "Patient(s) %s have no non-missing measurement of at least one biomarker up to ",
      "prediction_time, so no prediction can be made for them."),
      paste(unpredictable, collapse = ", ")), call. = FALSE)
  }
  data_predict_all
}

#' Recode ordinal biomarker responses onto the training categories
#'
#' @description Shared helper for \code{predictRisk},
#' \code{dynamicPredictionBio}, and \code{dynamicPredictionBioAll}. The
#' copula densities turn an ordinal response into its category code with
#' \code{as.numeric()} and index the fitted thresholds with it, and the
#' candidate categories of an ordinal \code{bio_i} are read off
#' \code{levels()}. Both used to come from the prediction data itself, so
#' any difference from the training categories -- unused levels dropped
#' (e.g. by \code{droplevels()} or subsetting and re-creating the factor),
#' a different level order, or a character column -- silently bracketed the
#' latent score between the wrong thresholds (in one check, a risk of 0.53
#' instead of 0.10). Each ordinal response is now rebuilt as an ordered
#' factor on the categories the \code{clmm()} fit was estimated with,
#' matching by label; a value that is not one of those categories is an
#' error. Missing values stay missing (\code{drop_missing_longitudinal()}
#' removes them afterwards).
#'
#' @return \code{data_predict_all}, with each ordinal response recoded.
#' @keywords internal
align_ordinal_levels <- function(data_predict_all, long_fit_all) {
  biomarker_type <- long_fit_all$biomarker_type
  if (is.null(biomarker_type)) return(data_predict_all)
  for (i in which(biomarker_type == "ordinal")) {
    resp_name <- all.vars(long_fit_all$long_sub_fixed[[i]])[1]
    training_levels <- long_fit_all$lfit[[i]]$y.levels
    if (is.null(training_levels) || !resp_name %in% names(data_predict_all[[i]])) next
    values <- as.character(data_predict_all[[i]][[resp_name]])
    unknown <- setdiff(unique(values[!is.na(values)]), training_levels)
    if (length(unknown) > 0) {
      stop(sprintf(paste0(
        "Ordinal biomarker '%s' in data_predict_all[[%d]] has value(s) %s that are not ",
        "categories of the fitted model (%s)."),
        resp_name, i, paste(sQuote(unknown, FALSE), collapse = ", "),
        paste(training_levels, collapse = " < ")), call. = FALSE)
    }
    data_predict_all[[i]][[resp_name]] <- factor(values, levels = training_levels, ordered = TRUE)
  }
  data_predict_all
}

#' Patient ids, in the order predictions are returned
#'
#' @description The density helpers (\code{marginalT()},
#' \code{conditionalYT()}, ...) all produce one column per patient, in order
#' of first appearance in \code{data_predict_all[[1]]}; this returns those
#' ids, as character, to name the returned predictions with.
#'
#' @return A character vector.
#' @keywords internal
prediction_patient_ids <- function(data_predict_all, long_fit_all) {
  num <- as.character(nlme::splitFormula(long_fit_all$long_sub_random[[1]], "|")[[2]])[2]
  unique(as.character(data_predict_all[[1]][[num]]))
}

#' Restrict prediction data to patients still at risk
#'
#' @description Shared helper for \code{predictRisk} and
#' \code{dynamicPredictionBio}: drops rows whose survival-time variable is
#' below \code{prediction_time} from every biomarker's data frame. Rows
#' whose survival time is missing (e.g. a new patient whose event time is
#' not yet known) are kept: they are treated as at risk. It is an error if
#' that leaves no patient (it used to return an empty result from
#' \code{predictRisk()} without comment, and fail with an unrelated error in
#' \code{predictLongitudinal()}).
#'
#' @return \code{data_predict_all}, filtered in place per element.
#' @keywords internal
subset_at_risk <- function(data_predict_all, survival_variable, prediction_time) {
  had_rows <- nrow(data_predict_all[[1]]) > 0
  for (i in seq_len(length(data_predict_all))) {
    s <- data_predict_all[[i]][[survival_variable]]
    data_predict_all[[i]] = data_predict_all[[i]][is.na(s) | s >= prediction_time, , drop = FALSE]
  }
  if (had_rows && nrow(data_predict_all[[1]]) == 0) {
    stop(sprintf(paste0(
      "No patient in data_predict_all is at risk at prediction_time = %g (every %s is below it), ",
      "so no prediction can be made. Set %s to NA for a patient whose event time is not yet known."),
      prediction_time, survival_variable, survival_variable), call. = FALSE)
  }
  data_predict_all
}

#' Upper limit of the prediction-to-infinity integration grid
#'
#' @description Shared helper for \code{predictRisk} and
#' \code{dynamicPredictionBio}. The denominator of a dynamic prediction
#' integrates over every event time after \code{prediction_time}, out to
#' infinity; the grid has to stop somewhere, and it should stop where the
#' probability left beyond it is negligible. So the upper limit is the
#' earliest time by which every at-risk patient's model-based conditional
#' survival \eqn{S(t \mid x) / S(s \mid x)}, \eqn{s} = \code{prediction_time},
#' has dropped below \code{tail_prob}, using the same (linearly
#' extrapolated, per-stratum) baseline cumulative hazard as
#' \code{marginalT()}. It never falls below \code{min_upper}.
#'
#' This replaces twice the largest \emph{observed survival time of the
#' patients being predicted}, which used each patient's own future outcome
#' (not available at \code{prediction_time}), stopped the integral early --
#' and inflated the risk -- for a patient whose event came soon after
#' \code{prediction_time}, and failed when that time was missing.
#'
#' Beyond the last time in the data \code{survivalSub()} was fit on, the
#' baseline hazard is an extrapolation; if it is so flat that the tail
#' criterion is not met by \code{max_multiple} times that last time, the
#' bound is capped there, with a warning if some patient's conditional
#' survival at the cap is still above \code{warn_prob}.
#'
#' @param data_predict_all At-risk prediction data (list of data frames).
#' @param long_fit_all Output of \code{longitudinalSub()} (for the id variable).
#' @param survival_fit_all Output of \code{survivalSub()}.
#' @param prediction_time The prediction (landmark) time.
#' @param min_upper The bound is at least this (e.g. \code{prediction_time +
#'   horizon}, so the denominator grid covers the prediction window).
#' @param tail_prob Remaining conditional survival probability treated as
#'   negligible.
#' @param max_multiple Cap, as a multiple of the last training time.
#' @param warn_prob Warn when the cap leaves more than this conditional
#'   survival probability unintegrated for some patient.
#' @return A single number.
#' @keywords internal
integration_upper_bound <- function(data_predict_all, long_fit_all, survival_fit_all,
                                    prediction_time, min_upper = prediction_time,
                                    tail_prob = 1e-4, max_multiple = 20, warn_prob = 0.01) {
  setup <- conditional_survival_setup(data_predict_all, long_fit_all, survival_fit_all)
  needed <- unlist(lapply(setup$groups, function(g)
    tail_time(g$cum_basehaz, setup$lp[g$patients], prediction_time, tail_prob)))
  needed <- needed[!is.na(needed)]

  cap <- max(max_multiple * setup$last_time, min_upper)
  upper <- max(c(needed, min_upper))
  if (upper > cap) {
    ### only worth a warning when the ignored conditional survival is
    ### material; it used to fire whenever it exceeded tail_prob, e.g. for
    ### an ordinary Weibull-like fit with 0.1% of the mass left at the cap
    left <- unlist(lapply(setup$groups, function(g) {
      dH <- diff(cumulative_baseline_at(g$cum_basehaz, c(prediction_time, cap)))
      exp(-dH * exp(setup$lp[g$patients]))
    }))
    left <- max(c(left[is.finite(left)], 0))
    if (left > warn_prob) {
      warning(sprintf(paste0(
        "The survival sub-model's extrapolated baseline hazard is so flat that up to %.1f%% of ",
        "a patient's conditional survival probability lies beyond %g (%g x the last follow-up ",
        "time), where the integration is capped; that probability mass is ignored."),
        100 * left, cap, max_multiple), call. = FALSE)
    }
    upper <- cap
  }
  upper
}

#' Per-patient pieces of the marginal survival model
#'
#' @description Shared by \code{integration_upper_bound()} and
#' \code{prepare_infinity_grid()}: each at-risk patient's Cox linear
#' predictor (\code{reference = "zero"}, \code{NA} if a covariate is
#' missing), and the patients grouped by the baseline cumulative hazard
#' they use (one group, or one per stratum of a stratified model).
#'
#' @return A list with \code{lp}, \code{groups} (each a list with
#' \code{cum_basehaz} and \code{patients}, indices into \code{lp}), and
#' \code{last_time} (the last time in the data \code{survivalSub()} was fit on).
#' @keywords internal
conditional_survival_setup <- function(data_predict_all, long_fit_all, survival_fit_all) {
  num <- as.character(nlme::splitFormula(long_fit_all$long_sub_random[[1]], "|")[[2]])[2]
  data.surv <- data_predict_all[[1]][!duplicated(data_predict_all[[1]][[num]]), , drop = FALSE]
  cum_basehaz_all <- survival_cum_basehaz(survival_fit_all)

  lp <- survival_lp(survival_fit_all, data.surv)
  patient_strata <- survival_patient_strata(survival_fit_all, data.surv)
  if (is.null(patient_strata)) {
    groups <- list(list(cum_basehaz = cum_basehaz_all[, c("hazard", "time")],
                        patients = seq_along(lp)))
  } else {
    groups <- lapply(intersect(unique(patient_strata), as.character(cum_basehaz_all$strata)), function(s)
      list(cum_basehaz = cum_basehaz_all[cum_basehaz_all$strata == s, c("hazard", "time")],
           patients = which(patient_strata == s)))
  }
  list(lp = lp, groups = groups, last_time = max(cum_basehaz_all$time))
}

#' Time by which conditional survival falls below a threshold
#'
#' @description Helper for \code{integration_upper_bound()}: for each linear
#' predictor in \code{lp}, the earliest time \eqn{t > s} with
#' \eqn{\exp(-(H_0(t) - H_0(s)) e^{lp}) <} \code{tail_prob}, where
#' \eqn{H_0} is the tabulated baseline cumulative hazard and, past its last
#' time, the same linear extrapolation \code{cumulative_baseline_at()} uses.
#' \code{Inf} if the extrapolated hazard is 0.
#'
#' @param cum_basehaz A data frame with columns \code{hazard} and \code{time}.
#' @param lp Linear predictors (\code{reference = "zero"}).
#' @param s The prediction time.
#' @param tail_prob The survival threshold.
#' @return A numeric vector, one time per element of \code{lp}.
#' @keywords internal
tail_time <- function(cum_basehaz, lp, s, tail_prob) {
  cum_basehaz <- cum_basehaz[order(cum_basehaz$time), c("hazard", "time")]
  slope <- extrapolation_slope(cum_basehaz)
  last <- nrow(cum_basehaz)
  H_s <- cumulative_baseline_at(cum_basehaz, s)
  target <- H_s - log(tail_prob) / exp(lp)
  after <- cum_basehaz$time > s
  vapply(target, function(h) {
    if (is.na(h)) return(NA_real_)
    hit <- which(after & cum_basehaz$hazard >= h)
    if (length(hit) > 0) return(cum_basehaz$time[hit[1]])
    if (slope <= 0) return(Inf)
    max(s, cum_basehaz$time[last]) +
      (h - max(H_s[1], cum_basehaz$hazard[last])) / slope
  }, numeric(1))
}

#' Build the prediction-to-infinity integration grid and marginal survival
#'
#' @description Shared helper for \code{predictRisk} and
#' \code{dynamicPredictionBio}: builds the numerical-integration grid from
#' \code{prediction_time} out to \code{upper_bound}, and evaluates the
#' marginal survival function \code{S(T)} over it.
#'
#' The \code{bandcount2 + 1} intervals are spaced by probability, not by
#' time: each holds an equal share of the at-risk patients' (averaged)
#' model-based conditional survival mass beyond \code{prediction_time}. The
#' upper bound sits far out in the tail (see
#' \code{integration_upper_bound()}), so equally spaced intervals would spend
#' most of the grid where there is almost no mass and converge more slowly
#' as \code{bandcount2} grows. The first interval starts exactly at
#' \code{prediction_time}; the previous equally spaced grid started half an
#' interval earlier, counting mass from before the prediction time. Each
#' grid point is its interval's midpoint.
#'
#' @return A list with \code{predict.time.infinity} (the grid points),
#' \code{predict.time.infinity.1} (the interval edges), and
#' \code{S_T_all_infinity}.
#' @keywords internal
prepare_infinity_grid <- function(data_predict_all, long_fit_all, survival_fit_all,
                                   prediction_time, upper_bound, bandcount2) {
  n_intervals <- bandcount2 + 1
  edges <- equal_mass_edges(data_predict_all, long_fit_all, survival_fit_all,
                            prediction_time, upper_bound, n_intervals)
  predict.time.infinity.1 = edges
  predict.time.infinity = (edges[-1] + edges[-length(edges)]) / 2

  S_T_all_infinity = marginalT(data_predict_all, long_fit_all, survival_fit_all,
                                l_i = predict.time.infinity.1, upper_bound)

  list(predict.time.infinity = predict.time.infinity,
       predict.time.infinity.1 = predict.time.infinity.1,
       S_T_all_infinity = S_T_all_infinity)
}

#' Interval edges holding equal conditional survival mass
#'
#' @description Helper for \code{prepare_infinity_grid()}: \code{n_intervals
#' + 1} edges from \code{prediction_time} to \code{upper_bound} such that the
#' average, over at-risk patients, of \eqn{S(t \mid x) / S(s \mid x)} drops
#' by the same amount across every interval. Falls back to equally spaced
#' edges if no patient has a usable linear predictor.
#'
#' @return A strictly increasing numeric vector of length \code{n_intervals + 1}.
#' @keywords internal
equal_mass_edges <- function(data_predict_all, long_fit_all, survival_fit_all,
                             prediction_time, upper_bound, n_intervals) {
  uniform <- seq(prediction_time, upper_bound, length.out = n_intervals + 1)
  setup <- conditional_survival_setup(data_predict_all, long_fit_all, survival_fit_all)
  fine <- seq(prediction_time, upper_bound, length.out = max(2000, 20 * n_intervals))
  surv_sum <- 0
  n_used <- 0
  for (g in setup$groups) {
    lp <- setup$lp[g$patients]
    lp <- lp[is.finite(lp)]
    if (length(lp) == 0) next
    H <- cumulative_baseline_at(g$cum_basehaz, c(prediction_time, fine))
    dH <- pmax(H[-1] - H[1], 0)
    surv_sum <- surv_sum + rowSums(exp(-outer(dH, exp(lp))))
    n_used <- n_used + length(lp)
  }
  if (n_used == 0) return(uniform)
  G <- cummin(surv_sum / n_used)  # conditional survival, forced non-increasing
  if (G[1] - G[length(G)] <= 0) return(uniform)
  levels <- seq(G[1], G[length(G)], length.out = n_intervals + 1)
  keep <- !duplicated(G)  # approx() needs distinct x
  edges <- stats::approx(rev(G[keep]), rev(fine[keep]), xout = levels, ties = "ordered")$y
  edges[1] <- prediction_time
  edges[length(edges)] <- upper_bound
  ### guard against flat stretches of G producing repeated edges
  if (any(diff(edges) <= 0)) {
    edges <- sort(unique(c(edges, uniform)))
    edges <- edges[round(seq(1, length(edges), length.out = n_intervals + 1))]
  }
  edges
}

#' Per-patient shift for exponentiating log densities
#'
#' @description The conditional-density helpers (\code{conditionalYT()} and
#' relatives) return log densities, one column per patient. Every quantity
#' a prediction needs is a ratio of sums of those densities within one
#' patient, so any per-patient constant cancels: this returns, for each
#' patient, the largest finite log density across all the matrices given,
#' to subtract before exponentiating (see \code{exp_shifted()}). Without it
#' the densities -- whose scale is the determinant of a covariance matrix
#' that grows with the number of observations and with the biomarkers'
#' units -- underflowed to 0 or overflowed to \code{Inf}.
#'
#' @param ... Matrices of log densities, rows = grid points, columns =
#'   patients (all with the same columns).
#' @return A numeric vector, one shift per patient (\code{0} for a patient
#'   with no finite value).
#' @keywords internal
patient_log_shift <- function(...) {
  stacked <- do.call(rbind, list(...))
  apply(stacked, 2, function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0) 0 else max(x)
  })
}

#' Exponentiate log densities after a per-patient shift
#'
#' @param log_density A matrix of log densities, columns = patients.
#' @param shift Per-patient shifts, from \code{patient_log_shift()}.
#' @return \code{exp(log_density - shift)}, column-wise.
#' @keywords internal
exp_shifted <- function(log_density, shift) {
  exp(log_density - rep(shift, each = nrow(log_density)))
}

#' Normalize a risk-probability ratio into a valid probability
#'
#' @description Shared helper for \code{predictRisk} and
#' \code{dynamicPredictionBio}: divides summed predicted-event mass by
#' summed total mass and clamps the result to \code{[0, 1]}.
#'
#' @return A numeric vector of risk probabilities in \code{[0, 1]}.
#' @keywords internal
clamp_risk_prob <- function(numerator_sum, denominator_sum) {
  risk.prob <- numerator_sum / denominator_sum
  risk.prob[risk.prob > 1] = 1
  risk.prob[risk.prob < 0] = 0
  risk.prob
}

#' Compare two prediction results' plain numeric-vector fields
#'
#' @description Shared comparison logic for \code{checkBandcountConvergence()}
#' and the \code{"auto"} bandcount support in \code{predictRisk()}/
#' \code{dynamicPredictionBio()}: only plain numeric vectors (no \code{dim})
#' that have the same length in both results are compared. This naturally
#' skips fields whose *size* is itself controlled by the bandcount being
#' varied (e.g. \code{dynamicPredictionBio()}'s \code{Y_density} matrix and
#' \code{Y_all} grid, whose resolution is exactly what \code{bandcount3}
#' sets), while still comparing the actual per-patient estimates derived
#' from them (\code{risk_prob_1}/\code{risk_prob_2}, \code{Y_predict}).
#' @return A list with \code{max} (the largest relative change across all
#' comparable fields, or \code{NA} if none were comparable) and
#' \code{by_field} (a named numeric vector, one entry per comparable
#' field).
#' @keywords internal
max_relative_diff <- function(result_a, result_b) {
  is_comparable <- function(x) is.numeric(x) && is.null(dim(x))
  common_fields <- intersect(
    names(result_a)[vapply(result_a, is_comparable, logical(1))],
    names(result_b)[vapply(result_b, is_comparable, logical(1))]
  )
  by_field <- c()
  for (field in common_fields) {
    a <- result_a[[field]]
    b <- result_b[[field]]
    if (length(a) == 0 || length(a) != length(b)) next
    ### a predicted biomarker value is compared on the scale of its own
    ### predictive distribution: relative to |value| alone, a value near 0
    ### (e.g. a centred or log-scale biomarker) turned a negligible absolute
    ### change into a large "relative" one and auto-tuning never converged.
    scale <- if (field == "Y_predict") predictive_sd(result_b) else 0
    by_field[field] <- max(abs(a - b) / pmax(abs(a), scale, 1e-8), na.rm = TRUE)
  }
  list(max = if (length(by_field) == 0) NA_real_ else max(by_field), by_field = by_field)
}

#' Standard deviation of each patient's predictive distribution
#'
#' @description Helper for \code{max_relative_diff()}: from a
#' \code{dynamicPredictionBio()}-style result (\code{Y_density}, one column
#' per patient, tabulated on \code{Y_all}), the standard deviation of each
#' patient's predicted biomarker distribution. \code{0} if the result has no
#' density.
#'
#' @param result A prediction result.
#' @return A numeric vector, one per patient, or \code{0}.
#' @keywords internal
predictive_sd <- function(result) {
  if (is.null(result$Y_density) || is.null(result$Y_all)) return(0)
  y <- as.numeric(result$Y_all)
  apply(as.matrix(result$Y_density), 2, function(d) {
    w <- d / sum(d)
    if (!all(is.finite(w))) return(0)
    m <- sum(w * y)
    sqrt(sum(w * (y - m)^2))
  })
}

#' Starting values for "auto" bandcount doubling
#'
#' @description The built-in starting point \code{auto_tune_bandcount()}
#' doubles from for each bandcount argument. \code{bandcount1}/
#' \code{bandcount2}/\code{bandcount3} used to default to fixed numbers
#' (\code{10}, \code{40}, and \code{300} respectively, across
#' \code{predictRisk()}/\code{dynamicPredictionBio()}); those same
#' numbers are reused here as the starting point for auto-tuning, so that
#' the first call \code{auto_tune_bandcount()} makes matches what a caller
#' relying on the old fixed defaults would have gotten. This cannot instead
#' be read off \code{formals(predict_fun)}, because that default is now the
#' literal string \code{"auto"} itself. The \code{bandcount3} start was
#' lowered from 300 to 100 when its grid started covering only where the
#' predictive densities are not negligible (see
#' \code{continuous_value_grid()}): 100 points there are finer than 300
#' were over the old fixed, wide range.
#' @keywords internal
bandcount_auto_start <- c(bandcount1 = 10, bandcount2 = 40, bandcount3 = 100)

#' Auto-select "auto" bandcount arguments by doubling until convergence
#'
#' @description Shared implementation backing \code{bandcount1}/
#' \code{bandcount2}/\code{bandcount3 = "auto"} support in
#' \code{predictRisk()}/\code{dynamicPredictionBio()}, and the
#' bandcount pre-resolution done once, up front, by \code{predictPlot()}/
#' \code{riskPlot()} (so their internal horizon/landmark loops do not repeat
#' the auto-tuning search on every iteration).
#'
#' Starts every argument named in \code{auto_names} at its entry in
#' \code{bandcount_auto_start()}, doubles all of them together each round,
#' and compares consecutive results with \code{max_relative_diff()} until
#' the largest relative change drops below \code{tol}, or \code{max_rounds}
#' extra doublings have been tried -- a hard cap, so this never loops
#' indefinitely: at most \code{max_rounds + 1} calls to \code{predict_fun}
#' (the default \code{max_rounds = 2} means at most 3 calls). If the cap is
#' hit without converging, a warning is issued and the result/bandcount at
#' the largest value tried is returned anyway, rather than erroring, so
#' automated pipelines are not interrupted.
#'
#' @param predict_fun \code{predictRisk} or \code{dynamicPredictionBio}.
#' @param args A named list of all of \code{predict_fun}'s arguments
#' (typically \code{as.list(environment())} captured right after argument
#' validation, before any other local variables are created).
#' @param auto_names Character vector naming which element(s) of \code{args}
#' to auto-tune (e.g. \code{"bandcount1"}, or \code{c("bandcount1", "bandcount2")}).
#' @return A list with \code{result} (\code{predict_fun}'s return value at
#' the resolved bandcount) and \code{bandcount} (a named list of the
#' resolved numeric bandcount value(s), one per element of \code{auto_names}).
#' @keywords internal
auto_tune_bandcount <- function(predict_fun, args, auto_names, tol = 0.01, max_rounds = 2, multiplier = 2) {
  for (n in auto_names) {
    args[[n]] <- unname(bandcount_auto_start[n])
  }

  prev_result <- do.call(predict_fun, args)
  round_i <- 0
  repeat {
    scaled_args <- args
    for (n in auto_names) scaled_args[[n]] <- scaled_args[[n]] * multiplier
    new_result <- do.call(predict_fun, scaled_args)
    round_i <- round_i + 1

    comparison <- max_relative_diff(prev_result, new_result)
    converged <- !is.na(comparison$max) && comparison$max < tol

    args <- scaled_args
    prev_result <- new_result

    if (converged || round_i >= max_rounds) {
      if (!converged) {
        warning(sprintf(paste0(
          "Auto-selected bandcount (%s) had not converged (max relative change %s) after ",
          "%d doubling(s) from the default; returning the result at the largest value tried ",
          "(%s). Pass an explicit, larger bandcount if you need tighter convergence, or use ",
          "checkBandcountConvergence() to investigate further."
        ), paste(auto_names, collapse = "/"),
           if (is.na(comparison$max)) "NA" else format(comparison$max, digits = 3),
           round_i,
           paste(sprintf("%s = %s", auto_names, unlist(args[auto_names])), collapse = ", ")
        ), call. = FALSE)
      }
      break
    }
  }
  list(result = prev_result, bandcount = args[auto_names])
}

#' Auto-select "auto" bandcount3 for a single biomarker's per-marker step
#'
#' @description \code{dynamicPredictionBioAll()}'s counterpart to
#' \code{auto_tune_bandcount()}: \code{bandcount3} only controls
#' \code{compute_bio_marker_step()}'s own candidate-value grid (\code{Y_all}),
#' not the shared step, so it is tuned per biomarker by re-calling
#' \code{compute_bio_marker_step()} directly against an already-computed
#' \code{shared} object (from \code{compute_bio_shared_step()}), rather than
#' re-running the whole \code{dynamicPredictionBio()} pipeline -- which would
#' recompute the shared denominator on every doubling round, for every
#' biomarker, exactly the redundant work \code{dynamicPredictionBioAll()} is
#' meant to avoid. Same doubling-until-stable check as
#' \code{auto_tune_bandcount()} (compares \code{Y_predict} via
#' \code{max_relative_diff()}, capped at \code{max_rounds} doublings, warns
#' instead of erroring if not converged by then).
#'
#' @param shared Output of \code{compute_bio_shared_step()}.
#' @inheritParams compute_bio_marker_step
#' @return A list with \code{result} (\code{compute_bio_marker_step()}'s
#' return value at the resolved \code{bandcount3}) and \code{bandcount3}
#' (the resolved numeric value).
#' @keywords internal
auto_tune_marker_bandcount3 <- function(shared, bio_i, long_fit_all, survival_fit_all,
                                         prediction_time, horizon, time_variable,
                                         survival_variable_all, survival_trans_function,
                                         tol = 0.01, max_rounds = 2, multiplier = 2) {
  bandcount3 <- unname(bandcount_auto_start["bandcount3"])
  prev_result <- compute_bio_marker_step(shared, bio_i, long_fit_all, survival_fit_all,
                                          prediction_time, horizon, time_variable,
                                          survival_variable_all, survival_trans_function, bandcount3)
  round_i <- 0
  repeat {
    bandcount3 <- bandcount3 * multiplier
    new_result <- compute_bio_marker_step(shared, bio_i, long_fit_all, survival_fit_all,
                                           prediction_time, horizon, time_variable,
                                           survival_variable_all, survival_trans_function, bandcount3)
    round_i <- round_i + 1

    comparison <- max_relative_diff(prev_result, new_result)
    converged <- !is.na(comparison$max) && comparison$max < tol

    prev_result <- new_result

    if (converged || round_i >= max_rounds) {
      if (!converged) {
        warning(sprintf(paste0(
          "Auto-selected bandcount3 for biomarker %d had not converged (max relative change %s) ",
          "after %d doubling(s) from the default; returning the result at the largest value tried ",
          "(bandcount3 = %s). Pass an explicit, larger bandcount3 if you need tighter convergence."
        ), bio_i,
           if (is.na(comparison$max)) "NA" else format(comparison$max, digits = 3),
           round_i, bandcount3
        ), call. = FALSE)
      }
      break
    }
  }
  list(result = prev_result, bandcount3 = bandcount3)
}
