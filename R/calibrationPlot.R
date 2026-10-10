#' Plot calibration of dynamic risk predictions
#'
#' @description
#' Checks whether the predicted risks are numerically right, not just well
#' ranked. At each landmark \code{s} in \code{prediction_time}, every subject
#' in \code{data_predict_all} still event-free at \code{s} gets a predicted
#' risk from \code{\link{predictRisk}} of an event in
#' \code{(s, s + horizon]}, using only their longitudinal history up to
#' \code{s}. The subjects are split into \code{n_groups} groups of about
#' equal size by predicted risk, and for each group the mean predicted risk
#' (horizontal axis) is plotted against the observed risk (vertical axis):
#' one minus the Kaplan--Meier estimate at \code{s + horizon} or, with
#' competing risks, the Aalen--Johansen cumulative incidence of that event
#' type, both of which account for censoring, with 95\% confidence
#' intervals. Points on the diagonal mean the predictions are well
#' calibrated; points below it mean the model overestimates the risk, above
#' it that it underestimates it.
#'
#' Evaluating on the data the sub-models were fit on gives an optimistic
#' (apparent) calibration; pass held-out validation data in
#' \code{data_predict_all} when available.
#'
#' For an \strong{interval-censored} fit (experimental), the observed risk
#' of a group is a nonparametric (Turnbull-type) estimate from the
#' subjects' intervals, left-truncated at each subject's last visit up to
#' the landmark -- with competing risks, the cumulative incidence estimate
#' of Hudgens, Satten and Longini (2001) -- so it estimates the same
#' quantity \code{predictRisk()} predicts, without using the model. Within
#' an interval the timing of the event is not identified, so the estimate
#' is less precise than with exact times, and no confidence interval is
#' drawn.
#'
#' @inheritParams performancePlot
#' @param n_groups Number of risk groups per landmark. Default \code{10}
#' (deciles); fewer groups give more stable observed risks in small data.
#'
#' @return A \code{ggplot} object with one panel per landmark (and, with
#' competing risks, per event type). Its \code{data} element is a
#' \code{data.frame} with one row per risk group, and columns
#' \code{landmark}, \code{cause}, \code{group}, \code{predicted},
#' \code{observed}, \code{lower}, \code{upper}, \code{n} (subjects in the
#' group) and \code{n_cases} (events of that type observed in the window).
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
#' # apparent calibration on the fitting data: 3-year landmark, 2-year window
#' calibrationPlot(pbc3, long_fit_all, survival_fit_all,
#'   prediction_time = 3, horizon = 2, time_variable = "year",
#'   trans$survival_variable_all, trans$survival_trans_function, n_groups = 5)
#' }
#'
#' @export
calibrationPlot <- function(data_predict_all, long_fit_all, survival_fit_all, prediction_time,
                            horizon, time_variable, survival_variable_all, survival_trans_function,
                            n_groups = 10, bandcount1 = "auto", bandcount2 = "auto") {

  assert_positive_integer(n_groups, "n_groups")
  preds <- landmark_predictions(data_predict_all, long_fit_all, survival_fit_all, prediction_time,
                                horizon, time_variable, survival_variable_all, survival_trans_function,
                                bandcount1, bandcount2)
  rows <- list()
  for (lp in preds) {
    for (k in seq_along(lp$risks)) {
      cal <- if (is_interval_censored(survival_fit_all)) {
        calibration_groups_interval(lp$risks[[k]], lp$outcome$L, lp$outcome$R, lp$outcome$cause,
                                    lp$last_visit, k, lp$landmark, lp$landmark + horizon, n_groups)
      } else {
        calibration_groups(lp$risks[[k]], lp$outcome$time, lp$outcome$status, lp$outcome$cause,
                           k, lp$landmark + horizon, n_groups)
      }
      cal$landmark <- lp$landmark
      cal$cause <- outcome_cause_label(survival_fit_all, k)
      rows[[length(rows) + 1]] <- cal
    }
  }
  cal <- do.call(rbind, rows)
  cal$cause <- factor(cal$cause, levels = unique(cal$cause))
  cal$panel <- factor(sprintf("Landmark %g", cal$landmark),
                      levels = sprintf("Landmark %g", sort(unique(cal$landmark))))
  cal <- cal[, c("landmark", "cause", "group", "predicted", "observed", "lower", "upper", "n",
                 "n_cases", "panel")]

  top <- max(c(cal$predicted, cal$upper), na.rm = TRUE)
  top <- min(1, top * 1.05)
  p <- ggplot(cal, aes(x = predicted, y = observed)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = top / 60, na.rm = TRUE) +
    geom_point(size = 2.5) +
    scale_x_continuous(limits = c(0, top)) +
    scale_y_continuous(limits = c(0, top), oob = function(x, ...) pmin(pmax(x, 0), top))
  p <- if (nlevels(cal$cause) > 1) {
    p + facet_grid(cause ~ panel)
  } else {
    p + facet_wrap(~ panel)
  }
  p + xlab(sprintf("Predicted risk within %g", horizon)) +
    ylab(sprintf("Observed risk within %g", horizon)) +
    theme_bw()
}

#' Predicted and observed risk per risk group at one landmark
#'
#' @param risk Predicted risks of event type \code{k} in the window.
#' @param time,status,cause Observed outcome of the same subjects (see
#' \code{performance_outcome()}); all have \code{time} after the landmark.
#' @param k Event type evaluated (ignored without competing risks, i.e.
#' when \code{cause} is all \code{NA}).
#' @param end End of the prediction window, \code{landmark + horizon}.
#' @param n_groups Number of groups (fewer if there are fewer subjects).
#' @return A \code{data.frame} with one row per group: \code{group},
#' \code{predicted} (mean predicted risk), \code{observed},
#' \code{lower}, \code{upper} (Kaplan--Meier/Aalen--Johansen estimate at
#' \code{end} with 95\% confidence interval), \code{n}, \code{n_cases}.
#' @keywords internal
calibration_groups <- function(risk, time, status, cause, k, end, n_groups) {
  n <- length(risk)
  n_groups <- min(n_groups, n)
  group <- ceiling(rank(risk, ties.method = "first") * n_groups / n)
  competing <- any(!is.na(cause))
  ### 0 = censored, 1 = event type k, 2 = any other event type
  state <- ifelse(status == 0, 0L, ifelse(!competing | cause == k, 1L, 2L))

  do.call(rbind, lapply(seq_len(n_groups), function(g) {
    in_g <- group == g
    obs <- c(est = NA_real_, lower = NA_real_, upper = NA_real_)
    if (any(state[in_g] != 0)) {
      st <- factor(state[in_g], levels = 0:2)
      fit <- survfit(Surv(time[in_g], st) ~ 1)
      at <- summary(fit, times = end, extend = TRUE)
      col <- match("1", fit$states)
      obs <- c(est = unname(at$pstate[1, col]),
               lower = if (is.null(at$lower)) NA else unname(at$lower[1, col]),
               upper = if (is.null(at$upper)) NA else unname(at$upper[1, col]))
    } else {
      ### nobody in the group had any event: observed risk 0
      obs[] <- c(0, NA, NA)
    }
    data.frame(group = g, predicted = mean(risk[in_g]), observed = unname(obs["est"]),
               lower = unname(obs["lower"]), upper = unname(obs["upper"]), n = sum(in_g),
               n_cases = sum(state[in_g] == 1 & time[in_g] <= end))
  }))
}
