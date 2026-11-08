#' Predictive distribution of the event time
#'
#' @description \strong{Experimental.} Each patient's whole predictive
#' distribution of the event time \eqn{T} given the biomarker history and
#' \eqn{T > s}, \eqn{p(t \mid Y, T > s) \propto f(Y \mid t) f(t)} for
#' \eqn{t > s}, in one pass: the cumulative distribution function at
#' \code{prediction_time + horizons} and its quantiles. The CDF at a horizon
#' is \code{predictRisk()}'s risk (summed over the event types under
#' competing risks), but computing the whole curve at once is much faster
#' than calling \code{predictRisk()} once per horizon. The quantiles give
#' prediction intervals for \eqn{T} itself, e.g. the 10\% quantile is a
#' time the patient outlives with probability 0.9.
#'
#' Beyond the last follow-up of the survival sub-model the distribution
#' rests on its extrapolated baseline hazard, so upper quantiles that fall
#' there are only as good as that tail (see
#' \code{\link{survivalSubExternal}}). Within follow-up the Breslow baseline
#' puts the probability at the observed event times, so the CDF rises in
#' steps there; quantiles are interpolated within the grid cell where the CDF
#' crosses the probability.
#'
#' @param data_predict_all,long_fit_all,survival_fit_all,time_variable,survival_variable_all,survival_trans_function
#'   As in \code{\link{predictRisk}} (right-censored survival fits only).
#' @param prediction_time The landmark \eqn{s}.
#' @param horizons Times after \code{prediction_time} at which to evaluate
#'   the CDF.
#' @param probs Probabilities of the quantiles.
#' @param n_grid Grid cells between \code{prediction_time} and the last
#'   horizon or the survival sub-model's last follow-up, whichever is later;
#'   \code{n_tail} grid cells beyond it.
#' @param n_tail See \code{n_grid}.
#' @return An object of class \code{"predictEventTime.BJM"}, a list with
#'   \code{cdf} (patients x horizons, \eqn{P(T \le s + h \mid Y, T > s)});
#'   with competing risks, \code{cif}, a list of two such matrices, one per
#'   event type (0, 1; they sum to \code{cdf}); \code{quantile} (patients x
#'   \code{probs}, absolute times); \code{prediction_time}; \code{horizons}.
#' @examples
#' \donttest{
#' data(pbc3)
#' surv <- survivalSub(pbc3[!duplicated(pbc3$id), ], Surv(years, status3) ~ age + sex, NULL)
#' long <- longitudinalSub(pbc3[pbc3$status3 == 1, ], list(serBilir ~ year + age + years),
#'                         list(~ year | id))
#' history <- pbc3[pbc3$id %in% c(2, 7) & pbc3$year <= 3, ]
#' history$years <- NA
#' out <- predictEventTime(history, long, surv, prediction_time = 3, horizons = 1:5,
#'                         time_variable = "year")
#' out$cdf
#' out$quantile
#' }
#' @export
predictEventTime <- function(data_predict_all, long_fit_all, survival_fit_all, prediction_time,
                             horizons, time_variable, survival_variable_all = NULL,
                             survival_trans_function = NULL, probs = c(0.1, 0.5, 0.9),
                             n_grid = 200, n_tail = 100) {
  assert_class(long_fit_all, "longitudinalSub.BJM", "long_fit_all", "longitudinalSub")
  assert_class(survival_fit_all, "survivalSub.BJM", "survival_fit_all", "survivalSub")
  if (is_interval_censored(survival_fit_all)) {
    stop("predictEventTime() supports right-censored survival fits only.", call. = FALSE)
  }
  assert_scalar_numeric(prediction_time, "prediction_time")
  if (!is.numeric(horizons) || length(horizons) == 0 || anyNA(horizons) || any(horizons <= 0)) {
    stop("`horizons` must be positive numbers.", call. = FALSE)
  }
  if (!is.numeric(probs) || anyNA(probs) || any(probs <= 0 | probs >= 1)) {
    stop("`probs` must be in (0, 1).", call. = FALSE)
  }
  assert_string(time_variable, "time_variable")
  horizons <- sort(unique(horizons))
  check_beyond_followup(survival_fit_all, prediction_time + max(horizons))
  if (!is.list(data_predict_all) || is.data.frame(data_predict_all)) {
    data_predict_all <- rep(list(data_predict_all), length(long_fit_all$lfit))
  }
  data_predict_all <- prepare_risk_data(data_predict_all, long_fit_all, survival_fit_all,
                                        prediction_time, time_variable, survival_variable_all)
  survival_variable <- survival_time_variable(survival_fit_all)
  id <- as.character(nlme::splitFormula(long_fit_all$long_sub_random[[1]], "|")[[2]])[2]
  patient_order <- unique(unlist(lapply(data_predict_all, function(d) as.character(d[[id]]))))
  parts <- split_no_measurement(data_predict_all, long_fit_all, outcome_variables(survival_fit_all),
                                survival_variable_all)
  measured <- subset_at_risk(parts$measured, survival_variable, prediction_time, allow_empty = TRUE)
  none <- if (!is.null(parts$no_measurement)) {
    subset_at_risk(parts$no_measurement, survival_variable, prediction_time, allow_empty = TRUE)
  }
  use_copula <- !is.null(long_fit_all$biomarker_type) && any(long_fit_all$biomarker_type == "ordinal")
  has_cr <- !is.null(survival_fit_all$glm_fit)

  s <- prediction_time
  # a fine grid out to the last horizon or the last follow-up, whichever is
  # later, then the survival model's equal-probability tail grid
  t_max <- max(s + max(horizons), survival_max_followup(survival_fit_all))
  edges <- sort(unique(c(seq(s, t_max, length.out = n_grid + 1), s + horizons)))
  one_patient <- function(data_j, use_y) {
    args <- list(data_j, use_y, long_fit_all, survival_fit_all, NA, NA, NA, time_variable,
                 survival_variable_all, survival_trans_function, n_tail,
                 if (use_copula) conditionalYTCopula else conditionalYT,
                 if (use_copula) conditionalYDTCopula else conditionalYDT)
    near <- do.call(subject_time_weights, c(replace(args, 5:6, list(s, t_max)), list(edges = edges)))
    far <- do.call(subject_time_weights, replace(args, 5:6, list(t_max, Inf)))
    log_w <- rbind(near$log_w, far$log_w)
    w <- exp(log_w - max(log_w[is.finite(log_w)]))
    w[!is.finite(w)] <- 0
    w <- w / sum(w)
    all_edges <- c(near$edges, far$edges[-1])
    cum <- rbind(0, apply(w, 2, cumsum))
    total <- rowSums(cum)
    at <- match(s + horizons, all_edges)
    list(cdf = total[at], cif = if (has_cr) lapply(seq_len(ncol(cum)), function(d) cum[at, d]),
         quantile = vapply(probs, function(p) {
           j <- which(total >= p)[1]
           if (is.na(j)) Inf else if (j == 1) all_edges[1] else
             stats::approx(total[(j - 1):j], all_edges[(j - 1):j], xout = p, ties = "ordered")$y
         }, numeric(1)))
  }
  rows <- list()
  for (part in list(list(data = measured, use_y = TRUE), list(data = none, use_y = FALSE))) {
    if (is.null(part$data) || nrow(part$data[[1]]) == 0) next
    for (p in unique(unlist(lapply(part$data, function(d) as.character(d[[id]]))))) {
      data_j <- lapply(part$data, function(d) d[as.character(d[[id]]) == p, , drop = FALSE])
      has_rows <- vapply(data_j, nrow, integer(1)) > 0
      if (!all(has_rows)) data_j <- rep(data_j[which(has_rows)[1]], length(data_j))
      rows[[p]] <- one_patient(data_j, part$use_y && all(has_rows))
    }
  }
  ids <- intersect(patient_order, names(rows))
  hn <- paste0("h", format(horizons, trim = TRUE))
  stack <- function(f) matrix(unlist(lapply(rows[ids], f)), length(ids), byrow = TRUE)
  cdf <- stack(function(r) r$cdf)
  dimnames(cdf) <- list(ids, hn)
  q <- stack(function(r) r$quantile)
  dimnames(q) <- list(ids, paste0(format(100 * probs, trim = TRUE), "%"))
  cif <- if (has_cr) lapply(1:2, function(d) {
    m <- stack(function(r) r$cif[[d]])
    dimnames(m) <- list(ids, hn)
    m
  })
  out <- list(cdf = cdf, cif = cif, quantile = q, prediction_time = s, horizons = horizons)
  class(out) <- "predictEventTime.BJM"
  out
}

#' Print a predictEventTime result
#' @param x A \code{"predictEventTime.BJM"} object.
#' @param digits Number of significant digits.
#' @param ... Not used.
#' @return \code{x}, invisibly.
#' @keywords internal
#' @export
print.predictEventTime.BJM <- function(x, digits = 3, ...) {
  cat(sprintf("Predictive distribution of the event time at prediction_time = %g\n", x$prediction_time))
  cat("\nP(T <= prediction_time + h | history, T > prediction_time):\n")
  print(signif(x$cdf, digits))
  cat("\nQuantiles of T (absolute times):\n")
  print(signif(x$quantile, digits))
  invisible(x)
}
