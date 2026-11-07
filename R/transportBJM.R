#' Transport a fitted backward joint model to a new population
#'
#' @description
#' \strong{Experimental.} Adapts a survival sub-model fitted in one population
#' (the \emph{source}, e.g. the hospital the model was developed in) to a
#' \emph{target} population (e.g. a new hospital) whose event rates differ,
#' while keeping the longitudinal sub-model unchanged.
#'
#' The backward joint model factors the joint distribution as
#' \eqn{f(T, D \mid X)\, f(Y \mid T, D, X)}: the event-time model and the
#' biomarker trajectory \emph{given} the event time. Between populations it is
#' often the first factor that changes most -- a different case mix, standard
#' of care, or referral pattern moves the event rates -- while how a biomarker
#' evolves in the years before an event is closer to biology and more stable.
#' This is the \emph{label shift} assumption of machine learning, and the
#' backward factorization is exactly the one under which it is easy to act
#' on: only \code{survival_fit_all} needs to change, and the result plugs into
#' \code{\link{predictRisk}}, \code{\link{predictLongitudinal}},
#' \code{\link{simulateTrajectory}}, \code{\link{performancePlot}} and
#' \code{\link{calibrationPlot}} with the source \code{long_fit_all}.
#'
#' Three methods are available, by how much is known about the target:
#' \describe{
#'   \item{\code{"baseline"}}{Target outcomes are available (\code{data_target},
#'   one row per subject). The hazard ratios are kept from the source and only
#'   the baseline cumulative hazard is re-estimated on the target (Breslow
#'   estimator with the source linear predictor as an offset). Needs few
#'   target events.}
#'   \item{\code{"recalibrate"}}{As \code{"baseline"}, but the source linear
#'   predictor is also rescaled by a calibration slope \eqn{\gamma}, estimated
#'   by a Cox model of the target outcome on the source linear predictor
#'   (\eqn{\gamma < 1}: the source hazard ratios are too extreme for the
#'   target).}
#'   \item{\code{"unlabeled"}}{No target outcomes at all, only biomarker
#'   histories (\code{data_predict_all}). The target baseline cumulative hazard
#'   is taken to be the source one times \eqn{e^{\alpha}}, and \eqn{\alpha} is
#'   estimated by maximizing the marginal likelihood of the target histories,
#'   \deqn{\prod_i \int_{s_i}^{\infty} f(Y_i \mid T)\, f_\alpha(T \mid X_i)\, dT,}
#'   where \eqn{s_i} is the time up to which patient \eqn{i} is known to be
#'   event-free (by default their last measurement). Because
#'   \eqn{f(Y \mid T)} is held fixed, the shape of the target trajectories --
#'   how close they look to the trajectories of patients about to have an
#'   event -- tells how far away the events are. This is the same integral
#'   \code{\link{predictRisk}} uses as its denominator, and it is only
#'   possible because the biomarker model conditions on the event time.}
#' }
#'
#' @param survival_fit_all Output of \code{\link{survivalSub}} fit on the
#'   source data (right-censored, not stratified).
#' @param method \code{"baseline"}, \code{"recalibrate"} or \code{"unlabeled"};
#'   see Description.
#' @param data_target For \code{"baseline"} and \code{"recalibrate"}: one row
#'   per target subject, with the outcome and covariates used in
#'   \code{survival_fit_all$form_marginal_surv} (and, with
#'   \code{update_event_type = TRUE}, \code{form_conditional_cr}).
#' @param update_event_type For \code{"baseline"}/\code{"recalibrate"} with
#'   competing risks: also refit the event-type model on the target
#'   (default \code{FALSE}, keep the source one).
#' @param data_predict_all,long_fit_all,time_variable,survival_variable_all,survival_trans_function
#'   For \code{"unlabeled"}: the target biomarker histories and the source
#'   longitudinal sub-model, as for \code{\link{predictRisk}}. The event-time
#'   column may be absent or \code{NA}.
#' @param landmark For \code{"unlabeled"}: the time up to which each target
#'   patient is known to be event-free, a single number or a vector named by
#'   patient id. Defaults to each patient's last measurement time.
#' @param data_source For \code{"unlabeled"}: optionally, biomarker histories
#'   from the source population (in the format of \code{data_predict_all}).
#'   The shift is then estimated on them too, where it should be 0, and
#'   subtracted: whatever misfit of the model biases the estimate in the same
#'   way in both populations cancels. Recommended.
#' @param cap_event_time For \code{"unlabeled"}: if \code{TRUE} (default),
#'   an event time after the source follow-up is given the biomarker
#'   trajectory of an event at the end of follow-up, as in
#'   \code{\link{simulateTrajectory}}, instead of extrapolating the
#'   longitudinal sub-model.
#' @param n_grid,upper For \code{"unlabeled"}: the number of grid intervals
#'   on \code{(0, upper]} used for the integral, and the upper end of the grid
#'   (default three times the source follow-up); the probability of an event
#'   beyond \code{upper} is placed at \code{upper}.
#' @param alpha_range For \code{"unlabeled"}: the search interval for
#'   \eqn{\alpha}, the log ratio of target to source baseline hazard.
#'
#' @return An object of class \code{c("transportBJM", "survivalSub.BJM")}:
#'   \code{survival_fit_all} with an extra element \code{transport} (the
#'   method, the target baseline cumulative hazard, the calibration slope,
#'   and \code{alpha} with its standard error for \code{"unlabeled"}). Use it
#'   in place of \code{survival_fit_all} anywhere in the package.
#'
#' @details Only the marginal event-time model is transported; the
#'   longitudinal sub-model, and by default the event-type model, are the
#'   source ones. If the biomarker trajectories before an event themselves
#'   differ between the populations (different assays, measurement
#'   schedules that depend on health), label shift does not hold and
#'   \code{"unlabeled"} is biased. \code{"unlabeled"} identifies a single
#'   proportional shift; it carries no information from patients without any
#'   measurement, and little from patients with only one.
#'
#' @examples
#' data(pbc3)
#' first <- pbc3[!duplicated(pbc3$id), ]
#' odd <- seq_len(nrow(first)) %% 2 == 1
#' source <- first[!odd, ]
#' target <- first[odd, ]
#' fit_src <- survivalSub(source, Surv(years, status3) ~ age + sex, NULL)
#' fit_tgt <- transportBJM(fit_src, "baseline", data_target = target)
#' fit_tgt
#' @export
transportBJM <- function(survival_fit_all, method = c("baseline", "recalibrate", "unlabeled"),
                         data_target = NULL, update_event_type = FALSE,
                         data_predict_all = NULL, long_fit_all = NULL, time_variable = NULL,
                         survival_variable_all = NULL, survival_trans_function = NULL,
                         landmark = NULL, data_source = NULL, cap_event_time = TRUE,
                         n_grid = 200, upper = NULL, alpha_range = c(-4, 4)) {
  assert_class(survival_fit_all, "survivalSub.BJM", "survival_fit_all", "survivalSub")
  method <- match.arg(method)
  if (is_interval_censored(survival_fit_all)) {
    stop("transportBJM() does not support interval-censored fits yet.", call. = FALSE)
  }
  ### start from the source fit itself, even if it was transported before
  source_fit <- survival_fit_all
  source_fit$transport <- NULL
  class(source_fit) <- "survivalSub.BJM"
  if (length(survival::untangle.specials(stats::terms(source_fit$coxph_fit), "strata")$vars) > 0) {
    stop("transportBJM() does not support stratified survival sub-models yet.", call. = FALSE)
  }

  transport <- if (method == "unlabeled") {
    transport_unlabeled(source_fit, data_predict_all, long_fit_all, time_variable,
                        survival_variable_all, survival_trans_function, landmark,
                        n_grid, upper, alpha_range, data_source, cap_event_time)
  } else {
    transport_labeled(source_fit, data_target, method)
  }

  out <- source_fit
  if (method != "unlabeled" && isTRUE(update_event_type) && length(source_fit$form_conditional_cr) != 0) {
    out$glm_fit <- refit_event_type(source_fit, data_target)
    transport$event_type_updated <- TRUE
  }
  out$transport <- transport
  class(out) <- c("transportBJM", "survivalSub.BJM")
  out
}

#' Breslow baseline cumulative hazard for a fixed linear predictor
#'
#' @param time,status Follow-up time and 0/1 event indicator.
#' @param lp Linear predictor (uncentered).
#' @return A data frame with columns \code{hazard} and \code{time}, at the
#'   distinct follow-up times, as \code{basehaz(centered = FALSE)} would give.
#' @keywords internal
breslow_basehaz <- function(time, status, lp) {
  ok <- !is.na(time) & !is.na(status) & !is.na(lp)
  time <- time[ok]; status <- status[ok]; risk <- exp(lp[ok])
  if (!any(status == 1)) {
    stop("The target data have no events: the baseline hazard cannot be estimated. ",
         "Use method = \"unlabeled\" instead.", call. = FALSE)
  }
  ### tabulated at every distinct follow-up time, as basehaz() does, so the
  ### extrapolation past the last time (see cumulative_baseline_at()) is the
  ### same as for a coxph() fit
  times <- sort(unique(time))
  d <- vapply(times, function(t) sum(time == t & status == 1), numeric(1))
  at_risk <- vapply(times, function(t) sum(risk[time >= t]), numeric(1))
  data.frame(hazard = cumsum(d / at_risk), time = times)
}

#' Target outcome of a survival sub-model formula
#' @keywords internal
target_surv <- function(source_fit, data_target) {
  assert_data_frame(data_target, "data_target")
  form <- source_fit$form_marginal_surv
  assert_vars_in_data(all.vars(form), data_target, "form_marginal_surv", "data_target")
  y <- eval(form[[2]], data_target, environment(form))
  if (ncol(y) != 2) stop("transportBJM() needs a right-censored Surv(time, status) outcome.", call. = FALSE)
  list(time = y[, "time"], status = as.numeric(y[, "status"] != 0))
}

#' The "baseline" and "recalibrate" methods of transportBJM()
#' @keywords internal
transport_labeled <- function(source_fit, data_target, method) {
  if (is.null(data_target)) {
    stop(sprintf("method = \"%s\" needs `data_target` (one row per target subject with the outcome).",
                 method), call. = FALSE)
  }
  y <- target_surv(source_fit, data_target)
  lp <- survival_lp(source_fit, data_target)
  slope <- 1
  slope_se <- NA_real_
  if (method == "recalibrate") {
    recal <- survival::coxph(survival::Surv(y$time, y$status) ~ lp, ties = "breslow")
    slope <- unname(stats::coef(recal))
    slope_se <- sqrt(recal$var[1, 1])
  }
  list(method = method,
       cum_basehaz = breslow_basehaz(y$time, y$status, slope * lp),
       lp_slope = slope, lp_slope_se = slope_se,
       max_followup = max(y$time, na.rm = TRUE),
       n_target = sum(!is.na(y$time)), events_target = sum(y$status == 1, na.rm = TRUE))
}

#' Refit the event-type model on the target
#' @keywords internal
refit_event_type <- function(source_fit, data_target) {
  y <- target_surv(source_fit, data_target)
  assert_vars_in_data(all.vars(source_fit$form_conditional_cr), data_target,
                      "form_conditional_cr", "data_target")
  stats::glm(source_fit$form_conditional_cr, data_target[y$status == 1 & !is.na(y$status), , drop = FALSE],
             family = stats::binomial)
}

#' The "unlabeled" method of transportBJM()
#'
#' @description Estimates the log baseline-hazard shift on the target
#' histories and, when \code{data_source} is given, on the source histories
#' too, and returns the difference (see \code{transportBJM()}).
#' @keywords internal
transport_unlabeled <- function(source_fit, data_predict_all, long_fit_all, time_variable,
                                survival_variable_all, survival_trans_function, landmark,
                                n_grid, upper, alpha_range, data_source, cap_event_time) {
  if (is.null(data_predict_all) || is.null(long_fit_all) || is.null(time_variable) ||
      is.null(survival_variable_all) || is.null(survival_trans_function)) {
    stop("method = \"unlabeled\" needs `data_predict_all`, `long_fit_all`, `time_variable`, ",
         "`survival_variable_all` and `survival_trans_function`.", call. = FALSE)
  }
  assert_class(long_fit_all, "longitudinalSub.BJM", "long_fit_all", "longitudinalSub")
  est <- function(d, lm) unlabeled_alpha(source_fit, d, long_fit_all, time_variable,
                                         survival_variable_all, survival_trans_function, lm,
                                         n_grid, upper, alpha_range, cap_event_time)
  tgt <- est(data_predict_all, landmark)
  ref <- if (!is.null(data_source)) est(data_source, NULL)
  alpha <- tgt$alpha - (if (is.null(ref)) 0 else ref$alpha)
  alpha_se <- sqrt(tgt$alpha_se^2 + (if (is.null(ref)) 0 else ref$alpha_se^2))
  if (min(abs(alpha - alpha_range)) < 1e-3) {
    warning("The estimate of alpha is at the edge of `alpha_range`; widen it.", call. = FALSE)
  }
  cbh_target <- survival_cum_basehaz(source_fit)[, c("hazard", "time")]
  cbh_target$hazard <- cbh_target$hazard * exp(alpha)
  list(method = "unlabeled", cum_basehaz = cbh_target, lp_slope = 1, lp_slope_se = NA_real_,
       alpha = alpha, alpha_se = alpha_se, alpha_target_raw = tgt$alpha,
       alpha_source_raw = if (!is.null(ref)) ref$alpha else NA_real_,
       max_followup = survival_max_followup(source_fit), n_target = tgt$n,
       events_target = NA_integer_)
}

#' Marginal-likelihood estimate of the log baseline-hazard shift
#'
#' @description Computes each patient's log density
#' \eqn{\log \sum_D f(Y \mid T, D) P(D \mid T)} on a grid of event times once
#' -- it does not depend on \eqn{\alpha} -- and then maximizes the marginal
#' log likelihood over \eqn{\alpha}, which only moves the event-time
#' probabilities on the grid.
#' @keywords internal
unlabeled_alpha <- function(source_fit, data_predict_all, long_fit_all, time_variable,
                            survival_variable_all, survival_trans_function, landmark,
                            n_grid, upper, alpha_range, cap_event_time) {
  n_long <- length(long_fit_all$lfit)
  if (is.data.frame(data_predict_all)) data_predict_all <- rep(list(data_predict_all), n_long)
  survival_variable <- survival_time_variable(source_fit)
  ### the event time is unknown: it is integrated over
  data_predict_all <- lapply(data_predict_all, function(d) {
    d[[survival_variable]] <- NA_real_
    d
  })
  id <- as.character(nlme::splitFormula(long_fit_all$long_sub_random[[1]], "|")[[2]])[2]
  ids <- unique(as.character(data_predict_all[[1]][[id]]))

  ### time each patient is known to be event-free
  last_time <- vapply(ids, function(i) max(unlist(lapply(data_predict_all, function(d)
    d[[time_variable]][as.character(d[[id]]) == i])), -Inf, na.rm = TRUE), numeric(1))
  s <- if (is.null(landmark)) last_time else if (length(landmark) == 1) rep(landmark, length(ids)) else {
    if (is.null(names(landmark)) || !all(ids %in% names(landmark))) {
      stop("`landmark` must be a single number or a vector named by every patient id.", call. = FALSE)
    }
    landmark[ids]
  }
  s <- pmax(s, 0)

  max_fu <- survival_max_followup(source_fit)
  if (is.null(upper)) upper <- 3 * max_fu
  edges <- seq(0, upper, length.out = n_grid + 1)
  mids <- (edges[-1] + edges[-(n_grid + 1)]) / 2
  grid_T <- c(mids, upper)  # the last point carries the mass beyond `upper`
  ### as in simulateTrajectory(): past the source follow-up the biomarker
  ### model is extrapolated, so an event there is given the trajectory of an
  ### event at the end of follow-up
  eval_T <- if (cap_event_time) pmin(grid_T, max_fu) else grid_T

  use_copula <- !is.null(long_fit_all$biomarker_type) && any(long_fit_all$biomarker_type == "ordinal")
  has_cr <- length(source_fit$form_conditional_cr) != 0
  if (has_cr) {
    ydt <- if (use_copula) conditionalYDTCopula else conditionalYDT
    logf_d <- ydt(data_predict_all, long_fit_all, source_fit, l_i = eval_T, survival_variable,
                  time_variable, survival_variable_all, survival_trans_function)
    p_d <- conditionalDT(data_predict_all, long_fit_all, source_fit, l_i = eval_T)
    a <- logf_d[[1]] + log(p_d[[1]])
    b <- logf_d[[2]] + log(p_d[[2]])
    m <- pmax(a, b)
    logf <- m + log(exp(a - m) + exp(b - m))
  } else {
    yt <- if (use_copula) conditionalYTCopula else conditionalYT
    logf <- yt(data_predict_all, long_fit_all, l_i = eval_T, survival_variable, time_variable,
               survival_variable_all, survival_trans_function)[[1]]
  }
  ### patients with no usable measurement carry no information on alpha
  usable <- colSums(is.finite(logf)) == nrow(logf)
  if (!any(usable)) stop("No patient has a usable biomarker measurement.", call. = FALSE)
  if (!all(usable)) {
    warning(sprintf("%d patient(s) without a usable measurement were left out.", sum(!usable)),
            call. = FALSE)
  }
  logf <- logf[, usable, drop = FALSE]
  s <- s[usable]
  shift <- apply(logf, 2, max)
  w <- exp(sweep(logf, 2, shift))  # grid x patient, in (0, 1]

  data_surv <- data_predict_all[[1]][!duplicated(data_predict_all[[1]][[id]]), , drop = FALSE]
  lp <- survival_lp(source_fit, data_surv)[usable]
  cbh <- survival_cum_basehaz(source_fit)[, c("hazard", "time")]
  ### per patient, interval (max(lo, s), max(hi, s)]: zero before s
  H_lo <- outer(edges[-(n_grid + 1)], s, pmax)
  H_hi <- outer(edges[-1], s, pmax)
  H_lo[] <- cumulative_baseline_at(cbh, H_lo)
  H_hi[] <- cumulative_baseline_at(cbh, H_hi)
  H_up <- cumulative_baseline_at(cbh, pmax(upper, s))
  H_s <- cumulative_baseline_at(cbh, s)

  loglik <- function(alpha) {
    r <- exp(alpha + lp)
    ### probabilities relative to S(s), times S(s) on the log scale
    rel_lo <- exp(-sweep(sweep(H_lo, 2, H_s), 2, r, `*`))
    rel_hi <- exp(-sweep(sweep(H_hi, 2, H_s), 2, r, `*`))
    mass <- rbind(rel_lo - rel_hi, exp(-(H_up - H_s) * r))
    sum(shift + log(colSums(w * mass)) - H_s * r)
  }
  opt <- stats::optimize(loglik, alpha_range, maximum = TRUE)
  alpha <- opt$maximum
  h <- 1e-3
  d2 <- (loglik(alpha + h) - 2 * opt$objective + loglik(alpha - h)) / h^2
  list(alpha = alpha, alpha_se = if (d2 < 0) sqrt(-1 / d2) else NA_real_, n = sum(usable))
}

#' @export
print.transportBJM <- function(x, digits = 4, ...) {
  tr <- x$transport
  cat("Backward joint model survival sub-model, transported to a target population\n")
  cat(sprintf("  Method                    : %s\n", tr$method))
  cat(sprintf("  Target subjects used      : %d\n", tr$n_target))
  if (!is.na(tr$events_target)) cat(sprintf("  Target events             : %d\n", tr$events_target))
  if (tr$method == "recalibrate") {
    cat(sprintf("  Calibration slope         : %s (SE %s)\n",
                format(tr$lp_slope, digits = digits), format(tr$lp_slope_se, digits = digits)))
  }
  if (tr$method == "unlabeled") {
    cat(sprintf("  log baseline-hazard shift : %s (SE %s), hazard ratio target/source %s\n",
                format(tr$alpha, digits = digits), format(tr$alpha_se, digits = digits),
                format(exp(tr$alpha), digits = digits)))
    if (!is.na(tr$alpha_source_raw)) {
      cat(sprintf("  (target %s minus source reference %s)\n",
                  format(tr$alpha_target_raw, digits = digits), format(tr$alpha_source_raw, digits = digits)))
    }
  }
  if (isTRUE(tr$event_type_updated)) cat("  Event-type model refit on the target\n")
  cat("\nSource fit:\n")
  NextMethod()
  invisible(x)
}
