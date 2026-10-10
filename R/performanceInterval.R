#' Observed outcome per subject for an interval-censored fit
#'
#' @description The interval-censored counterpart of
#' \code{performance_outcome()}: each subject's event is known only to lie in
#' \code{(L, R]} (\code{R = Inf} if no event was detected), with its type
#' when detected.
#' @param data The first \code{data_predict_all} data frame.
#' @param survival_fit_all An interval-censored \code{survivalSub.BJM} object.
#' @param id_variable The subject ID column.
#' @return A \code{data.frame} with columns \code{id}, \code{L}, \code{R},
#' \code{followup} (end of follow-up: \code{R} after an event, \code{L}
#' otherwise) and \code{cause} (1 or 2, matching
#' \code{risk_prob_1}/\code{risk_prob_2}; \code{NA} without competing risks
#' or without an event).
#' @keywords internal
performance_outcome_interval <- function(data, survival_fit_all, id_variable) {
  d <- data[!duplicated(data[[id_variable]]), , drop = FALSE]
  surv_lhs <- survival_fit_all$form_marginal_surv[[2]]
  assert_vars_in_data(all.vars(surv_lhs), d, "the outcome of survival_fit_all's survival formula",
                      "data_predict_all")
  b <- interval_bounds(eval(surv_lhs, d, environment(survival_fit_all$form_marginal_surv)))
  out <- data.frame(id = as.character(d[[id_variable]]), L = b$L, R = b$R, stringsAsFactors = FALSE)
  out$followup <- ifelse(is.finite(out$R), out$R, out$L)
  out$cause <- NA_integer_
  if (!is.null(survival_fit_all$glm_fit)) {
    type_var <- all.vars(survival_fit_all$form_conditional_cr[[2]])[1]
    assert_vars_in_data(type_var, d, "the event type of survival_fit_all's form_conditional_cr",
                        "data_predict_all")
    type <- d[[type_var]]
    out$cause <- ifelse(is.finite(out$R), ifelse(type == 0, 1L, 2L), NA_integer_)
    if (any(is.finite(out$R) & is.na(out$cause))) {
      stop("Some subjects with an event have a missing event type in data_predict_all.", call. = FALSE)
    }
  }
  out[!is.na(out$L), , drop = FALSE]
}

#' Posterior distribution of the event time at a landmark
#'
#' @description For the model-based performance measures: each patient's
#' distribution of \eqn{T} (and event type) given their history up to the
#' landmark \eqn{s} and \eqn{T > V}, their last visit -- the same quantity
#' \code{predictRisk()} integrates -- tabulated on the cells of
#' \code{predictRisk()}'s grids: \eqn{(V, s]} (\code{interval_gap_grid()})
#' and \eqn{(s, upper]} (\code{prepare_infinity_grid()}).
#' @param data_predict_all Landmark data (history up to \code{s}).
#' @param long_fit_all,survival_fit_all The fitted sub-models.
#' @param s The landmark.
#' @param time_variable,survival_variable_all,survival_trans_function As for
#'   \code{predictRisk()}.
#' @param n_inf,n_gap Grid cells after and before \code{s}.
#' @return A named list (patient ids) of lists with \code{edges} (cell
#'   boundaries, from \eqn{V}) and \code{mass} (cells x event types,
#'   summing to 1).
#' @keywords internal
event_time_posterior <- function(data_predict_all, long_fit_all, survival_fit_all, s,
                                 time_variable, survival_variable_all, survival_trans_function,
                                 n_inf = 100, n_gap = 20) {
  survival_variable <- survival_time_variable(survival_fit_all)
  data_predict_all <- align_ordinal_levels(data_predict_all, long_fit_all)
  data_predict_all <- drop_missing_longitudinal(data_predict_all, long_fit_all,
                                                outcome_variables(survival_fit_all), survival_variable_all)
  data_predict_all <- subset_at_risk(data_predict_all, survival_variable, s)
  ids <- prediction_patient_ids(data_predict_all, long_fit_all)
  has_cr <- !is.null(survival_fit_all$glm_fit)
  use_copula <- !is.null(long_fit_all$biomarker_type) && any(long_fit_all$biomarker_type == "ordinal")
  yt <- if (use_copula) conditionalYTCopula else conditionalYT
  ydt <- if (use_copula) conditionalYDTCopula else conditionalYDT

  upper <- integration_upper_bound(data_predict_all, long_fit_all, survival_fit_all, s)
  grid <- prepare_infinity_grid(data_predict_all, long_fit_all, survival_fit_all, s, upper, n_inf)
  mids <- grid$predict.time.infinity
  log_f <- function(d, l) {
    if (has_cr) {
      ydt(d, long_fit_all, survival_fit_all, l_i = l, survival_variable, time_variable,
          survival_variable_all, survival_trans_function)
    } else {
      yt(d, long_fit_all, l_i = l, survival_variable, time_variable, survival_variable_all,
         survival_trans_function)
    }
  }
  f_inf <- log_f(data_predict_all, mids)
  D_inf <- if (has_cr) conditionalDT(data_predict_all, long_fit_all, survival_fit_all, l_i = mids) else
    list(matrix(1, length(mids), length(ids)))
  gap <- interval_gap_grid(data_predict_all, long_fit_all, survival_fit_all, s, time_variable, n_gap)
  f_gap <- D_gap <- NULL
  if (!is.null(gap)) {
    f_gap <- gap$eval(log_f)
    D_gap <- if (has_cr) {
      gap$eval(function(d, l) conditionalDT(d, long_fit_all, survival_fit_all, l_i = l), fill = 0)
    } else list(matrix(1, n_gap, length(ids)))
  }
  shift <- do.call(patient_log_shift, c(f_inf, f_gap))
  K <- length(f_inf)

  stats::setNames(lapply(seq_along(ids), function(j) {
    cell <- function(f, D, S) exp(f[, j] - shift[j]) * D[, j] * S[, j]
    mass_inf <- sapply(seq_len(K), function(k) cell(f_inf[[k]], D_inf[[k]], grid$S_T_all_infinity))
    edges <- grid$predict.time.infinity.1
    mass <- matrix(mass_inf, ncol = K)
    if (!is.null(gap) && gap$has_gap[j]) {
      mass_gap <- matrix(sapply(seq_len(K), function(k) cell(f_gap[[k]], D_gap[[k]], gap$S)), ncol = K)
      edges <- c(seq(gap$last_visit[j], s, length.out = n_gap + 1), edges[-1])
      mass <- rbind(mass_gap, mass)
    }
    mass[!is.finite(mass)] <- 0
    list(edges = edges, mass = mass / sum(mass))
  }), ids)
}

#' Probability of a case and of a control from the posterior of T
#'
#' @description Model-based case/control status of one subject in the
#' window \eqn{(s, e]}, given what was observed: the event in \eqn{(L, R]}
#' with type \code{cause}, or no event by \eqn{L} (\code{R = Inf}). With
#' \eqn{F_k} the posterior cumulative mass of event type \eqn{k}
#' (\code{event_time_posterior()}, linear within cells), the probability
#' of being a case of type \code{k} is the posterior mass of
#' \eqn{(\max(L, s), \min(R, e)]} for the observed type over that of
#' \eqn{(L, R]}, and of being a control (event-free at \eqn{e}) the mass
#' beyond \eqn{e} over that of \eqn{(L, R]}. Yang et al. (2026), Biometrical
#' Journal 68:e70108.
#' @param post One element of \code{event_time_posterior()}.
#' @param L,R,cause The subject's observed outcome.
#' @param s,e The window.
#' @param k The event type evaluated.
#' @return \code{c(case = , control = )}.
#' @keywords internal
window_status_probs <- function(post, L, R, cause, s, e, k) {
  cum <- rbind(0, apply(post$mass, 2, cumsum))
  F_at <- function(x, types) {
    x <- min(max(x, post$edges[1]), post$edges[length(post$edges)])
    sum(vapply(types, function(t) stats::approx(post$edges, cum[, t], xout = x, ties = "ordered")$y,
               numeric(1)))
  }
  all_types <- seq_len(ncol(post$mass))
  lo <- max(L, post$edges[1])
  if (is.finite(R)) {
    d <- if (is.na(cause)) 1L else cause
    Z <- F_at(R, d) - F_at(lo, d)
    if (Z <= 0) {
      # no posterior mass on the observed interval: split it uniformly
      len <- R - lo
      case <- if (d == k && len > 0) max(0, min(R, e) - max(lo, s)) / len else 0
      control <- if (len > 0) max(0, R - max(lo, e)) / len else as.numeric(R > e)
      return(c(case = case, control = control))
    }
    case <- if (d == k) max(0, F_at(min(R, e), d) - F_at(max(lo, s), d)) / Z else 0
    control <- max(0, F_at(R, d) - F_at(max(lo, e), d)) / Z
  } else {
    Z <- 1 - F_at(lo, all_types)
    if (Z <= 0) return(c(case = 0, control = 1))
    case <- max(0, F_at(e, k) - F_at(max(lo, s), k)) / Z
    control <- max(0, 1 - F_at(max(lo, e), all_types)) / Z
  }
  c(case = min(1, case), control = min(1, control))
}

#' IPCW AUC and Brier score for interval-censored outcomes at one landmark
#'
#' @description Uses only subjects whose status in \eqn{(s, e]} is certain:
#' cases of type \code{k} (last negative visit \eqn{L \ge s}, event detected
#' at \eqn{R \le e}), the other event type detected the same way (a non-case
#' in the Brier score), and controls (\eqn{L \ge e}). They are weighted by the
#' inverse probability of still being followed up -- at \eqn{R} for an
#' event, at \eqn{e} for a control -- from the Kaplan--Meier estimate of the
#' end of follow-up among the subjects at risk at \eqn{s}, as in Yang et al.
#' (2026), Biometrical Journal 68:e70108. Subjects whose interval straddles
#' \eqn{s} or \eqn{e} get weight 0.
#' @param risk Predicted risks of event type \code{k}.
#' @param L,R,followup,cause Observed outcome (see
#'   \code{performance_outcome_interval()}).
#' @param k,s,horizon Event type, landmark and window length.
#' @return As \code{performance_metrics()}, plus \code{n_uncertain}.
#' @keywords internal
performance_metrics_interval_ipcw <- function(risk, L, R, followup, cause, k, s, horizon) {
  end <- s + horizon
  detected_in_window <- is.finite(R) & L >= s & R <= end
  is_case <- detected_in_window & (is.na(cause) | cause == k)
  is_control <- L >= end

  cens <- survfit(Surv(followup, as.numeric(!is.finite(R))) ~ 1)
  G <- stats::stepfun(cens$time, c(1, cens$surv))
  G_minus <- function(t) G(t - sqrt(.Machine$double.eps))
  w <- numeric(length(risk))
  w[detected_in_window] <- 1 / G_minus(R[detected_in_window])
  w[is_control] <- 1 / G(end)
  w[!is.finite(w)] <- 0

  brier <- sum(w * (as.numeric(is_case) - risk)^2) / length(risk)
  auc <- NA_real_
  if (any(is_case) && any(is_control)) {
    rc <- risk[is_case]
    wc <- w[is_case]
    rn <- risk[is_control]
    concordant <- vapply(rc, function(r) sum(rn < r) + 0.5 * sum(rn == r), numeric(1))
    auc <- sum(wc * concordant) / (sum(wc) * length(rn))
  }
  list(auc = auc, brier = brier, n_at_risk = length(risk), n_cases = sum(is_case),
       n_uncertain = sum(!detected_in_window & !is_control))
}

#' Model-based AUC and Brier score for interval-censored outcomes
#'
#' @description Every subject at risk counts, as a case with probability
#' \code{p_case} and a control with probability \code{p_control}
#' (\code{window_status_probs()}): the expected squared error for the Brier
#' score, and case/control pairs of different subjects weighted by
#' \eqn{p_i q_j} for the AUC (Yang et al., 2026). These weights come from the
#' model being evaluated, so the measures lean optimistic when that model is
#' misspecified.
#' @param risk Predicted risks of the event type evaluated.
#' @param p_case,p_control Case and control probabilities per subject.
#' @return As \code{performance_metrics()}, with \code{n_cases} the expected
#'   number of cases.
#' @keywords internal
performance_metrics_interval_model <- function(risk, p_case, p_control) {
  brier <- mean(p_case * (1 - risk)^2 + (1 - p_case) * risk^2)
  wins <- outer(risk, risk, function(a, b) (a > b) + 0.5 * (a == b))
  pair <- outer(p_case, p_control)
  diag(pair) <- 0
  auc <- if (sum(pair) > 0) sum(pair * wins) / sum(pair) else NA_real_
  list(auc = auc, brier = brier, n_at_risk = length(risk), n_cases = sum(p_case))
}

#' Nonparametric estimate of the event-time distribution from interval data
#'
#' @description Self-consistency (EM) estimate in the spirit of Turnbull
#' (1976) for interval-censored event times \eqn{(L, R]}, left-truncated at
#' \eqn{V} (each subject is only in the sample because \eqn{T > V}), with an
#' optional known event type per detected event (the competing-risks
#' estimate of Hudgens, Satten & Longini, 2001, Biometrics 57:74). Mass sits
#' on the cells between consecutive distinct time points (all \eqn{L},
#' \eqn{R}, \eqn{V} and \code{extra}); within an observed interval it is not
#' identified, and the EM's choice there is what the estimate reports.
#' @param L,R,V Interval and truncation time per subject (\code{R = Inf}
#'   when no event was detected).
#' @param cause Event type per subject (1..K, \code{NA} without an event),
#'   or \code{NULL}.
#' @param extra Additional time points to put cell boundaries at.
#' @param tol,max_iter EM convergence settings.
#' @return A list with \code{edges} and \code{mass} (cells x event types,
#'   the last cell running to \code{Inf}).
#' @keywords internal
npmle_interval <- function(L, R, V, cause = NULL, extra = NULL, tol = 1e-9, max_iter = 5000) {
  K <- if (is.null(cause)) 1L else max(2L, max(cause, na.rm = TRUE))
  cause_i <- if (is.null(cause)) ifelse(is.finite(R), 1L, NA_integer_) else cause
  pts <- sort(unique(c(0, L, R[is.finite(R)], V, extra)))
  edges <- c(pts, Inf)
  J <- length(edges) - 1
  lo <- edges[-(J + 1)]
  hi <- edges[-1]
  n <- length(L)
  # cell j is inside (L_i, R_i] / beyond V_i
  inside <- outer(L, lo, "<=") & outer(R, hi, ">=")
  beyond <- outer(V, lo, "<=")
  # which event types each subject's observation allows (all when no event)
  allowed <- matrix(1, n, K)
  known <- !is.na(cause_i)
  allowed[known, ] <- 0
  allowed[cbind(which(known), cause_i[known])] <- 1
  A <- inside * 1
  not_beyond <- (!beyond) * 1
  mass <- matrix(1 / (J * K), J, K)
  for (it in seq_len(max_iter)) {
    # E-step: each subject's share of every (cell, type) inside its
    # observation, plus the truncated "ghost" subjects with T <= V_i
    den <- rowSums((A %*% mass) * allowed)
    share <- allowed / ifelse(den > 0, den, Inf)
    kept <- c((1 - not_beyond) %*% rowSums(mass))
    ghosts <- c(crossprod(not_beyond, 1 / ifelse(kept > 0, kept, Inf)))
    expected <- mass * (crossprod(A, share) + ghosts)
    new <- expected / sum(expected)
    converged <- max(abs(new - mass)) < tol
    mass <- new
    if (converged) break
  }
  list(edges = edges, mass = mass)
}

#' Observed risk per risk group for interval-censored outcomes
#'
#' @description The interval-censored counterpart of
#' \code{calibration_groups()}: within each group of predicted risk, the
#' observed risk of an event of type \code{k} in \eqn{(s, e]} for a subject
#' known to be event-free at their last visit \eqn{V} is
#' \eqn{(F_k(e) - F_k(s)) / (1 - F(V))} with \eqn{F} from
#' \code{npmle_interval()}, averaged over the group's subjects -- the same
#' quantity \code{predictRisk()} predicts. No confidence interval is given.
#' @param risk Predicted risks of event type \code{k}.
#' @param L,R,cause Observed outcome (see \code{performance_outcome_interval()}).
#' @param V Each subject's last visit up to the landmark.
#' @param k,s,end Event type, landmark and window end.
#' @param n_groups Number of groups.
#' @return As \code{calibration_groups()}.
#' @keywords internal
calibration_groups_interval <- function(risk, L, R, cause, V, k, s, end, n_groups) {
  n <- length(risk)
  n_groups <- min(n_groups, n)
  group <- ceiling(rank(risk, ties.method = "first") * n_groups / n)
  competing <- any(!is.na(cause))
  do.call(rbind, lapply(seq_len(n_groups), function(g) {
    in_g <- group == g
    fit <- npmle_interval(L[in_g], R[in_g], V[in_g], if (competing) cause[in_g], extra = c(s, end))
    cum <- rbind(0, apply(fit$mass, 2, cumsum))
    F_at <- function(x, types) {
      j <- findInterval(x, fit$edges)
      sum(cum[j, types])
    }
    all_types <- seq_len(ncol(fit$mass))
    kk <- if (competing) k else 1L
    num <- F_at(end, kk) - F_at(s, kk)
    obs <- mean(vapply(V[in_g], function(v) {
      den <- 1 - F_at(v, all_types)
      if (den > 0) num / den else NA_real_
    }, numeric(1)), na.rm = TRUE)
    data.frame(group = g, predicted = mean(risk[in_g]), observed = obs,
               lower = NA_real_, upper = NA_real_, n = sum(in_g),
               n_cases = sum(is.finite(R[in_g]) & L[in_g] >= s & R[in_g] <= end &
                               (!competing | cause[in_g] == k)))
  }))
}
