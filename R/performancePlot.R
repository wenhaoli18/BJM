#' Plot predictive performance across landmark times
#'
#' @description
#' Evaluates how well the backward joint model's dynamic risk predictions
#' discriminate and are calibrated, at each of several landmark times. For
#' each landmark \code{s} in \code{prediction_time}, every subject in
#' \code{data_predict_all} still event-free at \code{s} gets a predicted risk
#' from \code{\link{predictRisk}}, using only their longitudinal history up
#' to \code{s}, of an event in \code{(s, s + horizon]}; this is compared with
#' what actually happened, by
#' \describe{
#'   \item{AUC}{the time-dependent (cumulative/dynamic) area under the ROC
#'   curve: the probability that a subject with an event in the window has a
#'   higher predicted risk than one still event-free at \code{s + horizon}.
#'   0.5 is no better than chance, 1 is perfect.}
#'   \item{Brier score}{the mean squared difference between the predicted
#'   risk and the observed event indicator in the window. Lower is better.}
#' }
#' Subjects censored within the window are handled by inverse probability of
#' censoring weighting (IPCW), with the censoring distribution among the
#' subjects at risk at \code{s} estimated by Kaplan--Meier. With competing
#' risks, both measures are computed separately for each event type: the
#' cases are the subjects with that event type in the window, and for the
#' Brier score a subject with the other event type in the window counts as
#' not having the event.
#'
#' Evaluating on the data the sub-models were fit on gives an optimistic
#' (apparent) performance; pass held-out validation data in
#' \code{data_predict_all} when available.
#'
#' \strong{Interval-censored fits} (experimental). The event is then only
#' known to lie between the last negative visit \eqn{L} and the visit
#' \eqn{R} that detected it. A subject is at risk at \code{s} when followed
#' up past \code{s} without a detected event, and the predicted risk is
#' that of \code{predictRisk()} (conditional on being event-free at the last
#' visit). Following Yang, Rizopoulos, Newcomb and Erler (2026, Biometrical
#' Journal 68:e70108), \code{interval_method} chooses how subjects whose
#' interval straddles \code{s} or \code{s + horizon} are handled:
#' \describe{
#'   \item{\code{"ipcw"} (default)}{only subjects whose status is certain
#'   are used -- cases detected inside the window after a negative visit at
#'   or after \code{s}, controls with a negative visit at or after
#'   \code{s + horizon} -- weighted by the inverse Kaplan--Meier probability
#'   of still being followed up. It does not depend on the model being
#'   evaluated, so it is the fair choice for comparing models, but it is
#'   more variable and loses subjects when visits are sparse; in Yang et
#'   al.'s simulations it underestimated the Brier score.}
#'   \item{\code{"model"}}{every subject at risk counts, as a case or a
#'   control with the probability the fitted model gives their event time
#'   given their interval. Less variable, but optimistic when the model is
#'   misspecified, since the model grades itself.}
#' }
#'
#' @param data_predict_all The evaluation data, in the same format as for
#' \code{\link{predictRisk}}: a list of long-format \code{data.frame}s, one
#' per longitudinal outcome (or a single \code{data.frame} used for all).
#' Unlike for \code{predictRisk()}, it must also contain each subject's
#' observed outcome: the survival time and status variables of
#' \code{survival_fit_all}'s Cox formula and, with competing risks, the event
#' type variable of its \code{form_conditional_cr}. The outcome is taken
#' from the first data frame and is hidden from \code{predictRisk()}.
#' @param long_fit_all A \code{longitudinalSub.BJM} object from
#' \code{\link{longitudinalSub}}.
#' @param survival_fit_all A \code{survivalSub.BJM} object from
#' \code{\link{survivalSub}}.
#' @param prediction_time A vector of landmark times.
#' @param horizon The prediction horizon (a single positive number): the
#' window after each landmark in which events are counted.
#' @param time_variable The name of the time variable in the linear mixed models.
#' @param survival_variable_all,survival_trans_function As for
#' \code{\link{predictRisk}}.
#' @param bandcount1,bandcount2 As for \code{\link{predictRisk}}; passed to
#' it at each landmark.
#' @param interval_method Only for an interval-censored
#' \code{survival_fit_all}: \code{"ipcw"} (default) or \code{"model"}; see
#' Description. For an interval-censored fit, \code{data_predict_all} must
#' contain the interval columns of \code{form_marginal_surv} (and the event
#' type) instead of a survival time and status.
#'
#' @return A \code{ggplot} object with one panel per measure, the landmark
#' time on the horizontal axis. Its \code{data} element is a
#' \code{data.frame} with one row per landmark, measure (and event type), and
#' columns \code{landmark}, \code{measure}, \code{cause}, \code{value},
#' \code{n_at_risk} (subjects predicted) and \code{n_cases} (events of that
#' type observed in the window).
#'
#' @examples
#' \donttest{
#' data(pbc3)
#' survival_fit_all <- survivalSub(pbc3[!duplicated(pbc3$id), ],
#'                                 Surv(years, status3) ~ age + sex, NULL)
#' long_sub_fixed <- list("long1" = serBilir ~ year + age + sex + years)
#' long_sub_random <- list("long1" = ~ year | id)
#' long_fit_all <- longitudinalSub(list(pbc3[pbc3$status3 == 1, ]),
#'                                 long_sub_fixed, long_sub_random)
#' trans <- survivalTrans(c(1, 3, 5, 7))
#'
#' # apparent performance on the fitting data, 2-year window
#' performancePlot(pbc3, long_fit_all, survival_fit_all,
#'   prediction_time = c(1, 3, 5), horizon = 2, time_variable = "year",
#'   trans$survival_variable_all, trans$survival_trans_function)
#' }
#'
#' @export
performancePlot <- function(data_predict_all, long_fit_all, survival_fit_all, prediction_time,
                            horizon, time_variable, survival_variable_all, survival_trans_function,
                            bandcount1 = "auto", bandcount2 = "auto",
                            interval_method = c("ipcw", "model")) {
  interval_method <- match.arg(interval_method)
  interval <- is_interval_censored(survival_fit_all)
  preds <- landmark_predictions(data_predict_all, long_fit_all, survival_fit_all, prediction_time,
                                horizon, time_variable, survival_variable_all, survival_trans_function,
                                bandcount1, bandcount2, posterior = interval && interval_method == "model")
  rows <- list()
  for (lp in preds) {
    for (k in seq_along(lp$risks)) {
      o <- lp$outcome
      m <- if (!interval) {
        performance_metrics(lp$risks[[k]], o$time, o$status, o$cause, k, lp$landmark, horizon)
      } else if (interval_method == "ipcw") {
        performance_metrics_interval_ipcw(lp$risks[[k]], o$L, o$R, o$followup, o$cause, k,
                                          lp$landmark, horizon)
      } else {
        probs <- t(vapply(seq_len(nrow(o)), function(i) {
          window_status_probs(lp$posterior[[o$id[i]]], o$L[i], o$R[i], o$cause[i],
                              lp$landmark, lp$landmark + horizon, k)
        }, numeric(2)))
        performance_metrics_interval_model(lp$risks[[k]], probs[, "case"], probs[, "control"])
      }
      rows[[length(rows) + 1]] <- data.frame(
        landmark = lp$landmark, measure = c("AUC", "Brier score"),
        cause = outcome_cause_label(survival_fit_all, k),
        value = c(m$auc, m$brier), n_at_risk = m$n_at_risk, n_cases = m$n_cases)
    }
  }
  perf <- do.call(rbind, rows)
  perf$measure <- factor(perf$measure, levels = c("AUC", "Brier score"))
  perf$cause <- factor(perf$cause, levels = unique(perf$cause))
  if (anyNA(perf$value[perf$measure == "AUC"])) {
    message("AUC is undefined (no cases or no event-free subjects in the window) at some landmarks; ",
            "those points are left out.")
  }

  reference <- data.frame(measure = factor("AUC", levels = levels(perf$measure)), y = 0.5)
  p <- ggplot(perf, aes(x = landmark, y = value)) +
    geom_hline(data = reference, aes(yintercept = y), linetype = "dashed", color = "grey50")
  if (nlevels(perf$cause) > 1) {
    p <- p + geom_line(aes(color = cause), na.rm = TRUE) + geom_point(aes(color = cause), size = 2.5, na.rm = TRUE) +
      labs(color = "Event type")
  } else {
    p <- p + geom_line(na.rm = TRUE) + geom_point(size = 2.5, na.rm = TRUE)
  }
  p + facet_wrap(~ measure, scales = "free_y") +
    xlab(sprintf("Landmark time (prediction window = %g)", horizon)) + ylab(NULL) +
    theme_bw()
}

#' Predicted risks and observed outcomes at each landmark time
#'
#' @description Shared by \code{performancePlot()} and
#' \code{calibrationPlot()}: validates their common arguments and, for each
#' landmark \code{s}, calls \code{predictRisk()} on the subjects
#' event-free at \code{s}, with their history up to \code{s} and their
#' outcome hidden.
#' @inheritParams performancePlot
#' @return A list with one element per landmark that had subjects to
#' predict, each a list with \code{landmark}, \code{outcome} (see
#' \code{performance_outcome()}, rows matching the predictions) and
#' \code{risks} (a list with one vector of predicted risks per event type).
#' @keywords internal
landmark_predictions <- function(data_predict_all, long_fit_all, survival_fit_all, prediction_time,
                                 horizon, time_variable, survival_variable_all, survival_trans_function,
                                 bandcount1, bandcount2, posterior = FALSE) {
  assert_class(long_fit_all, "longitudinalSub.BJM", "long_fit_all", "longitudinalSub")
  assert_class(survival_fit_all, "survivalSub.BJM", "survival_fit_all", "survivalSub")
  interval <- is_interval_censored(survival_fit_all)
  assert_data_list(data_predict_all, "data_predict_all", length(long_fit_all$lfit), allow_bare_df = TRUE)
  if (!is.list(data_predict_all) || is.data.frame(data_predict_all)) {
    data_predict_all <- rep(list(data_predict_all), each = length(long_fit_all$lfit))
  }
  if (!is.numeric(prediction_time) || length(prediction_time) == 0 || anyNA(prediction_time)) {
    stop("`prediction_time` must be a non-empty numeric vector of landmark times.", call. = FALSE)
  }
  assert_scalar_numeric(horizon, "horizon", positive = TRUE)
  assert_string(time_variable, "time_variable")

  id_variable <- as.character(nlme::splitFormula(long_fit_all$long_sub_random[[1]], "|")[[2]])[2]
  if (interval) {
    ### event known only to lie in (L, R]; at risk at s while followed up
    ### past s without a detected event (follow-up ends at R or L)
    outcome <- performance_outcome_interval(data_predict_all[[1]], survival_fit_all, id_variable)
    outcome$time <- outcome$followup
    hidden <- unique(c(all.vars(survival_fit_all$form_marginal_surv[[2]]),
                       if (!is.null(survival_fit_all$glm_fit)) all.vars(survival_fit_all$form_conditional_cr[[2]])))
  } else {
    outcome <- performance_outcome(data_predict_all[[1]], survival_fit_all, id_variable)
    hidden <- unique(c(all.vars(formula(survival_fit_all$coxph_fit)[[2]]),
                       if (!is.null(survival_fit_all$glm_fit)) all.vars(survival_fit_all$form_conditional_cr[[2]])))
  }
  for (i in seq_along(data_predict_all)) {
    assert_vars_in_data(c(time_variable, id_variable), data_predict_all[[i]],
                        "time_variable/the subject ID", sprintf("data_predict_all[[%d]]", i))
  }

  out <- list()
  for (s in sort(unique(prediction_time))) {
    at_risk <- outcome$id[outcome$time > s]
    if (length(at_risk) == 0) {
      message(sprintf("No subject is event-free at landmark %g; skipping it.", s))
      next
    }
    ### each subject's history up to the landmark, with the outcome hidden
    landmark_data <- lapply(data_predict_all, function(d) {
      d <- d[as.character(d[[id_variable]]) %in% at_risk & d[[time_variable]] <= s, , drop = FALSE]
      d[intersect(hidden, names(d))] <- NA
      ### the event time an interval-censored fit conditions on: unknown
      if (interval) d[[survival_time_variable(survival_fit_all)]] <- rep(NA_real_, nrow(d))
      d
    })
    if (nrow(landmark_data[[1]]) == 0) {
      message(sprintf("No at-risk subject has a measurement by landmark %g; skipping it.", s))
      next
    }
    risk <- predictRisk(landmark_data, long_fit_all, survival_fit_all, prediction_time = s,
                        horizon = horizon, time_variable = time_variable,
                        survival_variable_all = survival_variable_all,
                        survival_trans_function = survival_trans_function,
                        bandcount1 = bandcount1, bandcount2 = bandcount2)
    o <- outcome[match(names(risk$risk_prob_1), outcome$id), , drop = FALSE]
    risks <- if (is.null(risk$risk_prob_2)) list(risk$risk_prob_1) else
      list(risk$risk_prob_1, risk$risk_prob_2)
    lp <- list(landmark = s, outcome = o, risks = risks)
    if (interval) {
      ### each subject's last visit up to s: predictRisk() conditions on T > V
      times <- unlist(lapply(landmark_data, function(d) d[[time_variable]]))
      ids <- unlist(lapply(landmark_data, function(d) as.character(d[[id_variable]])))
      lp$last_visit <- unname(tapply(times, ids, max)[o$id])
      if (posterior) {
        lp$posterior <- event_time_posterior(
          landmark_data, long_fit_all, survival_fit_all, s, time_variable, survival_variable_all,
          survival_trans_function,
          n_inf = if (is.numeric(bandcount2)) bandcount2 else 100,
          n_gap = if (is.numeric(bandcount1)) bandcount1 else 20)
      }
    }
    out[[length(out) + 1]] <- lp
  }
  if (length(out) == 0) {
    stop("No landmark in `prediction_time` had subjects to evaluate.", call. = FALSE)
  }
  out
}

#' One row per subject with the observed outcome, for performancePlot()
#'
#' @param data The first \code{data_predict_all} data frame.
#' @param survival_fit_all A \code{survivalSub.BJM} object.
#' @param id_variable The subject ID column.
#' @return A \code{data.frame} with columns \code{id}, \code{time},
#' \code{status} (1 = event, 0 = censored) and \code{cause} (1 or 2 for the
#' event type matching \code{risk_prob_1}/\code{risk_prob_2} of
#' \code{predictRisk()}; \code{NA} without competing risks or when censored).
#' @keywords internal
performance_outcome <- function(data, survival_fit_all, id_variable) {
  d <- data[!duplicated(data[[id_variable]]), , drop = FALSE]
  surv_lhs <- formula(survival_fit_all$coxph_fit)[[2]]
  assert_vars_in_data(all.vars(surv_lhs), d, "the outcome of survival_fit_all's Cox formula",
                      "data_predict_all")
  y <- eval(surv_lhs, d, environment(formula(survival_fit_all$coxph_fit)))
  out <- data.frame(id = as.character(d[[id_variable]]), time = y[, "time"], status = y[, "status"])

  out$cause <- NA_integer_
  if (!is.null(survival_fit_all$glm_fit)) {
    cr_lhs <- survival_fit_all$form_conditional_cr[[2]]
    assert_vars_in_data(all.vars(cr_lhs), d, "the event type of survival_fit_all's form_conditional_cr",
                        "data_predict_all")
    type <- eval(cr_lhs, d, environment(survival_fit_all$form_conditional_cr))
    ### predictRisk()'s risk_prob_1 is the glm response's "failure" outcome
    ### (0, FALSE, or a factor's first level) and risk_prob_2 its "success"
    first <- if (is.factor(type)) type == levels(type)[1] else as.numeric(type) == 0
    out$cause <- ifelse(out$status == 1, ifelse(first, 1L, 2L), NA_integer_)
    if (any(out$status == 1 & is.na(out$cause))) {
      stop("Some subjects with an event have a missing event type in data_predict_all.", call. = FALSE)
    }
  }
  out[!is.na(out$time) & !is.na(out$status), , drop = FALSE]
}

#' Legend label for event type k in performancePlot()
#'
#' @keywords internal
outcome_cause_label <- function(survival_fit_all, k) {
  if (is.null(survival_fit_all$glm_fit)) {
    return("Event")
  }
  lev <- if (is.factor(survival_fit_all$glm_fit$model[[1]])) levels(survival_fit_all$glm_fit$model[[1]]) else c("0", "1")
  sprintf("Cause %d (%s = %s)", k, deparse(survival_fit_all$form_conditional_cr[[2]]), lev[k])
}

#' IPCW time-dependent AUC and Brier score at one landmark
#'
#' @param risk Predicted risks of event type \code{k} in the window.
#' @param time,status,cause Observed outcome of the same subjects (see
#' \code{performance_outcome()}); all have \code{time > s}.
#' @param k Event type evaluated (ignored when \code{cause} is all
#' \code{NA}, i.e. without competing risks).
#' @param s,horizon Landmark and window length.
#' @return A list with \code{auc}, \code{brier}, \code{n_at_risk},
#' \code{n_cases}.
#' @keywords internal
performance_metrics <- function(risk, time, status, cause, k, s, horizon) {
  end <- s + horizon
  in_window <- time <= end & status == 1
  is_case <- in_window & (is.na(cause) | cause == k)
  is_control <- time > end

  ### Kaplan-Meier of the censoring time among subjects at risk at s;
  ### G(t-) for an event at t, and G(end) for subjects event-free at end
  cens <- survfit(Surv(time, 1 - status) ~ 1)
  G <- stats::stepfun(cens$time, c(1, cens$surv))
  G_minus <- function(t) G(t - sqrt(.Machine$double.eps))
  w <- numeric(length(risk))
  w[in_window] <- 1 / G_minus(time[in_window])
  w[is_control] <- 1 / G(end)
  w[!is.finite(w)] <- 0

  brier <- sum(w * (as.numeric(is_case) - risk)^2) / length(risk)

  auc <- NA_real_
  if (any(is_case) && any(is_control)) {
    rc <- risk[is_case]
    wc <- w[is_case]
    rn <- risk[is_control]
    wn <- w[is_control][1]
    ### for each case, the (equally weighted) controls it outranks, ties half
    concordant <- vapply(rc, function(r) sum(rn < r) + 0.5 * sum(rn == r), numeric(1))
    auc <- sum(wc * wn * concordant) / (sum(wc) * wn * length(rn))
  }
  list(auc = auc, brier = brier, n_at_risk = length(risk), n_cases = sum(is_case))
}
