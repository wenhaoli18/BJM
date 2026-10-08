#' Fit the backward joint model to interval-censored event times
#'
#' @description \strong{Experimental.} With interval censoring nobody's
#' event time \eqn{T} is known exactly -- only \eqn{L < T \le R} -- so the
#' longitudinal sub-model \eqn{f(Y \mid T)}, which has \eqn{T} as a
#' covariate, cannot be fit on the subjects with an observed event the way
#' it is for right-censored data. \code{fitIntervalBJM()} fills in \eqn{T}
#' by stochastic EM, then by multiple imputation:
#' \enumerate{
#'   \item fit the marginal survival model to the intervals
#'   (\code{survivalSub()} with \code{Surv(L, R, type = "interval2")}); it
#'   does not depend on \eqn{T}, so it is fit once;
#'   \item start every interval-censored \eqn{T} at its interval midpoint
#'   and fit \code{longitudinalSub()};
#'   \item draw each subject's \eqn{T} from
#'   \eqn{f(T \mid L < T \le R, Y, X) \propto f(Y \mid T, X) f(T \mid X)}
#'   on \eqn{(L, R]} (\code{imputeEventTime()}), refit
#'   \code{longitudinalSub()}, and repeat for \code{n_burnin} iterations;
#'   \item keep iterating, saving the next \code{n_imputations} draws and
#'   their \code{longitudinalSub()} fits.
#' }
#' Step 3 uses the backward model's own factorization: \eqn{f(Y \mid T)}
#' on a grid of \eqn{T} is exactly what dynamic prediction already
#' computes (\code{conditionalYT()}), so each draw costs about one
#' \code{predictRisk()} call per subject.
#'
#' By default only subjects whose event was observed (finite \eqn{R}) enter
#' the longitudinal fit, mirroring the complete-case fit used for
#' right-censored data, which is valid when censoring is independent of
#' \eqn{(T, Y)} given the covariates. \code{include_right_censored = TRUE}
#' also imputes \eqn{T \in (L, \infty)} for right-censored subjects and
#' uses them; their draws then rely on the survival model's extrapolated
#' tail.
#'
#' The saved draws come from successive iterations of one chain with the
#' model parameters at their point estimates, not from their posterior, so
#' the between-imputation variance in \code{pooled} can be somewhat
#' understated.
#'
#' @param data_survival One row per subject: the interval endpoints, the
#'   survival covariates and the subject id.
#' @param data_fit_all As in \code{\link{longitudinalSub}}: one data frame
#'   per biomarker (or one shared data frame), long format, for every
#'   subject. The \code{event_time} column and the columns named in
#'   \code{survival_variable_all} are (re)filled from the imputed \eqn{T}.
#' @param form_marginal_surv \code{Surv(L, R, type = "interval2") ~
#'   covariates}.
#' @param event_time Name of the event-time column used in
#'   \code{long_sub_fixed}.
#' @param long_sub_fixed,long_sub_random As in \code{\link{longitudinalSub}}.
#' @param time_variable Name of the visit-time column.
#' @param survival_variable_all,survival_trans_function Transformed event-time
#'   columns, as in \code{\link{predictRisk}} (\code{NULL} if none).
#' @param form_conditional_cr Optional competing-risks (event-type) model, as
#'   in \code{\link{survivalSub}}, e.g. \code{type ~ Tev + x}. The event type
#'   (0/1, in \code{data_survival} and in every \code{data_fit_all} data
#'   frame, alongside the covariates of this formula) is taken as known once
#'   the event is detected; each subject's \eqn{T} is then drawn from
#'   \eqn{f(Y \mid T, D) P(D \mid T) f(T)}, and this model, which uses
#'   \eqn{T}, is refit on every draw.
#' @param n_burnin Stochastic-EM iterations before draws are saved.
#' @param n_imputations Number of saved draws (and longitudinal fits).
#' @param n_grid Grid points per subject interval.
#' @param baseline,df Baseline hazard of the survival model (see
#'   \code{\link{survivalSub}}).
#' @param include_right_censored See Description.
#' @param seed Optional integer seed (the global RNG state is restored).
#'
#' @return An object of class \code{"fitIntervalBJM"}, a list with
#'   \code{survival_fit_all} (with competing risks, the event-type model fit
#'   to the first imputation; \code{survival_fit_all_list} has one per
#'   imputation); \code{long_fit_all_list} (one
#'   \code{longitudinalSub()} fit per imputation); \code{long_fit_all} (the
#'   first of them, for use with \code{predictRisk()} and friends);
#'   \code{pooled} (\code{\link{poolLongitudinalSub}} of the fits, when
#'   \code{n_imputations > 1}); \code{imputed_T} (a subjects x imputations
#'   matrix); and \code{trace} (each iteration's fixed effects, to check
#'   the chain has settled).
#' @export
fitIntervalBJM <- function(data_survival, data_fit_all, form_marginal_surv, event_time,
                           long_sub_fixed, long_sub_random, time_variable,
                           survival_variable_all = NULL, survival_trans_function = NULL,
                           form_conditional_cr = NULL,
                           n_burnin = 10, n_imputations = 5, n_grid = 40,
                           baseline = c("spline", "piecewise"), df = NULL,
                           include_right_censored = FALSE, seed = NULL) {
  if (!is.null(seed)) local_r_seed(seed)
  if (inherits(long_sub_fixed, "formula")) long_sub_fixed <- list(long_sub_fixed)
  if (inherits(long_sub_random, "formula")) long_sub_random <- list(long_sub_random)
  M <- length(long_sub_fixed)
  if (is.data.frame(data_fit_all)) data_fit_all <- rep(list(data_fit_all), M)
  id <- as.character(nlme::splitFormula(long_sub_random[[1]], "|")[[2]])[2]

  has_cr <- length(form_conditional_cr) != 0
  if (has_cr && include_right_censored) {
    stop("include_right_censored = TRUE is not supported with competing risks (the event type of a right-censored subject is unknown).",
         call. = FALSE)
  }
  event_type <- if (has_cr) all.vars(form_conditional_cr[[2]])[1] else NULL
  bounds <- subject_intervals(data_survival, form_marginal_surv, id, event_type)
  # the event-type model needs an event time: start from interval midpoints
  with_event_time <- function(T_values) {
    ds <- data_survival
    ds[[event_time]] <- unname(T_values[as.character(ds[[id]])])
    ds
  }
  mid <- stats::setNames(ifelse(is.finite(bounds$R), (bounds$L + bounds$R) / 2, NA_real_), bounds$id)
  survival_fit_all <- survivalSub(with_event_time(mid), form_marginal_surv, form_conditional_cr,
                                  event_time = event_time, baseline = baseline, df = df)
  if (!is_interval_censored(survival_fit_all)) {
    stop("`form_marginal_surv` must have an interval-censored outcome, Surv(L, R, type = \"interval2\").",
         call. = FALSE)
  }
  refit_event_type <- function(T_values) {
    if (!has_cr) return(survival_fit_all)
    out <- survival_fit_all
    out$glm_fit <- fit_event_type_interval(with_event_time(T_values), form_marginal_surv,
                                           form_conditional_cr, event_time)
    out
  }
  use <- is.finite(bounds$R) | include_right_censored
  bounds <- bounds[use, , drop = FALSE]
  data_fit_all <- lapply(data_fit_all, function(d) d[as.character(d[[id]]) %in% bounds$id, , drop = FALSE])

  T_cur <- stats::setNames(ifelse(is.finite(bounds$R), (bounds$L + bounds$R) / 2, NA_real_), bounds$id)
  if (anyNA(T_cur)) {
    # right-censored subjects start one median residual life past L
    T_cur[is.na(T_cur)] <- bounds$L[is.na(T_cur)] +
      stats::median(bounds$R[is.finite(bounds$R)] - bounds$L[is.finite(bounds$R)])
  }

  # burn-in fits are thrown away, so their Sigma_fit EM convergence
  # warnings are muffled; the saved fits keep theirs
  fit_long <- function(T_values, quiet) {
    filled <- fill_event_time(data_fit_all, id, T_values, event_time,
                              survival_variable_all, survival_trans_function)
    withCallingHandlers(longitudinalSub(filled, long_sub_fixed, long_sub_random),
                        warning = function(w) {
                          if (quiet && grepl("Sigma_fit", conditionMessage(w))) invokeRestart("muffleWarning")
                        })
  }
  # the starting fit can fail like any later one (nlme::lme() trouble with
  # a random-effect variance near zero): retry from event times drawn
  # uniformly within each interval
  long_fit_all <- tryCatch(fit_long(T_cur, quiet = TRUE), error = function(e) e)
  for (attempt in seq_len(5)) {
    if (!inherits(long_fit_all, "error")) break
    T_try <- T_cur
    finite <- is.finite(bounds$R) & bounds$R > bounds$L
    T_try[finite] <- stats::runif(sum(finite), bounds$L[finite], bounds$R[finite])
    long_fit_all <- tryCatch(fit_long(T_try, quiet = TRUE), error = function(e) e)
    if (!inherits(long_fit_all, "error")) T_cur <- T_try
  }
  if (inherits(long_fit_all, "error")) {
    stop(sprintf(paste0(
      "longitudinalSub() failed on the starting event times and on 5 random restarts: %s ",
      "This usually means a random-effect variance is close to zero; a simpler long_sub_random ",
      "(e.g. ~ 1 | id) may help."), conditionMessage(long_fit_all)), call. = FALSE)
  }

  n_iter <- n_burnin + n_imputations
  trace <- vector("list", n_iter)
  imputed_T <- matrix(NA_real_, nrow(bounds), n_imputations, dimnames = list(bounds$id, NULL))
  long_fit_all_list <- vector("list", n_imputations)
  survival_fit_all_list <- vector("list", n_imputations)
  for (it in seq_len(n_iter)) {
    # an unlucky draw occasionally makes nlme::lme() fail (e.g. "Singularity
    # in backsolve"); draw again from the same fit rather than abort the chain
    for (attempt in seq_len(5)) {
      T_new <- imputeEventTime(data_fit_all, long_fit_all, survival_fit_all, bounds, time_variable,
                               survival_variable_all, survival_trans_function, n_grid = n_grid)
      new_fit <- tryCatch(fit_long(T_new, quiet = it <= n_burnin), error = function(e) e)
      if (!inherits(new_fit, "error")) break
    }
    if (inherits(new_fit, "error")) {
      stop(sprintf("longitudinalSub() failed on 5 successive imputations at iteration %d: %s",
                   it, conditionMessage(new_fit)), call. = FALSE)
    }
    T_cur <- T_new
    long_fit_all <- new_fit
    # the event-type model uses T too: refit it on the same draw
    survival_fit_all <- refit_event_type(T_cur)
    trace[[it]] <- c(unlist(lapply(long_fit_all$lfit, nlme::fixef)),
                     if (has_cr) stats::setNames(stats::coef(survival_fit_all$glm_fit),
                                                 paste0("cr.", names(stats::coef(survival_fit_all$glm_fit)))))
    if (it > n_burnin) {
      imputed_T[, it - n_burnin] <- T_cur[bounds$id]
      long_fit_all_list[[it - n_burnin]] <- long_fit_all
      survival_fit_all_list[[it - n_burnin]] <- survival_fit_all
    }
  }

  out <- list(survival_fit_all = survival_fit_all_list[[1]], long_fit_all_list = long_fit_all_list,
              survival_fit_all_list = survival_fit_all_list,
              long_fit_all = long_fit_all_list[[1]],
              pooled = if (n_imputations > 1) poolLongitudinalSub(long_fit_all_list),
              imputed_T = imputed_T, trace = do.call(rbind, trace))
  class(out) <- "fitIntervalBJM"
  out
}

#' Draw interval-censored event times given the biomarker history
#'
#' @description \strong{Experimental.} One draw of each subject's event time
#' from \eqn{f(T \mid L < T \le R, Y, X) \propto f(Y \mid T, X) f(T \mid
#' X)}: \eqn{f(Y \mid T)} from the longitudinal sub-model on a grid of
#' \code{n_grid} intervals tiling \eqn{(L, R]} (\code{conditionalYT()}),
#' times the survival model's probability of each grid interval
#' (\code{marginalT()}); an interval is drawn with those weights and
#' \eqn{T} uniformly within it. A right-censored subject (\eqn{R = \infty})
#' uses the same equal-probability grid from \eqn{L} out to the survival
#' model's tail as \code{predictRisk()}. A subject with \eqn{L = R} keeps
#' \eqn{T = L}. Used by \code{\link{fitIntervalBJM}}.
#'
#' @param data_fit_all Longitudinal data, one data frame per biomarker.
#' @param long_fit_all Output of \code{longitudinalSub()}.
#' @param survival_fit_all Interval-censored output of \code{survivalSub()}.
#' @param bounds Data frame with columns \code{id}, \code{L}, \code{R}
#'   (and, with competing risks, the observed event type \code{D}).
#' @param time_variable,survival_variable_all,survival_trans_function As in
#'   \code{predictRisk()}.
#' @param n_grid Grid intervals per subject.
#' @param n_draws Independent draws per subject (the weights are computed
#'   once and reused).
#' @return A named numeric vector of draws, one per row of \code{bounds};
#'   with \code{n_draws > 1}, a matrix with one column per draw.
#' @keywords internal
imputeEventTime <- function(data_fit_all, long_fit_all, survival_fit_all, bounds, time_variable,
                            survival_variable_all = NULL, survival_trans_function = NULL,
                            n_grid = 40, n_draws = 1) {
  id <- as.character(nlme::splitFormula(long_fit_all$long_sub_random[[1]], "|")[[2]])[2]
  survival_variable <- survival_time_variable(survival_fit_all)
  use_copula <- !is.null(long_fit_all$biomarker_type) && any(long_fit_all$biomarker_type == "ordinal")
  conditionalYT_fun <- if (use_copula) conditionalYTCopula else conditionalYT
  conditionalYDT_fun <- if (use_copula) conditionalYDTCopula else conditionalYDT
  has_cr <- !is.null(survival_fit_all$glm_fit)
  if (has_cr && (is.null(bounds$D) || anyNA(bounds$D[is.finite(bounds$R)]))) {
    stop("With competing risks, `bounds` needs the event type D of every subject with an observed event.",
         call. = FALSE)
  }

  data_fit_all <- align_ordinal_levels(data_fit_all, long_fit_all)
  data_fit_all <- suppressWarnings(drop_missing_longitudinal(data_fit_all, long_fit_all, survival_variable,
                                                             survival_variable_all))
  out <- matrix(NA_real_, nrow(bounds), n_draws, dimnames = list(bounds$id, NULL))
  for (j in seq_len(nrow(bounds))) {
    L <- bounds$L[j]
    R <- bounds$R[j]
    if (is.finite(R) && R <= L) {
      out[j, ] <- L
      next
    }
    data_j <- lapply(data_fit_all, function(d) d[as.character(d[[id]]) == bounds$id[j], , drop = FALSE])
    has_rows <- vapply(data_j, nrow, integer(1)) > 0
    if (!any(has_rows)) {
      # nothing usable at all (not even covariates): uniform on the interval
      out[j, ] <- if (is.finite(R)) stats::runif(n_draws, L, R) else L
      next
    }
    # a biomarker without usable rows: f(Y | T) is left out, so T is drawn
    # from the survival model alone (covariates read from another biomarker)
    if (!all(has_rows)) data_j <- rep(data_j[which(has_rows)[1]], length(data_j))
    if (is.finite(R)) {
      edges <- seq(L, R, length.out = n_grid + 1)
      mass <- marginalT(data_j, long_fit_all, survival_fit_all, l_i = edges, upper_bound = R)[, 1]
    } else {
      upper <- integration_upper_bound(data_j, long_fit_all, survival_fit_all, L, min_upper = L)
      grid <- prepare_infinity_grid(data_j, long_fit_all, survival_fit_all, L, upper, n_grid - 1)
      edges <- grid$predict.time.infinity.1
      mass <- grid$S_T_all_infinity[, 1]
    }
    mids <- (edges[-1] + edges[-length(edges)]) / 2
    if (has_cr) {
      # the event type is known: weight by P(D = d | T) and use f(Y | T, D = d)
      k <- bounds$D[j] + 1
      mass <- mass * conditionalDT(data_j, long_fit_all, survival_fit_all, l_i = mids)[[k]][, 1]
      log_f <- if (all(has_rows)) {
        conditionalYDT_fun(data_j, long_fit_all, survival_fit_all, l_i = mids, survival_variable,
                           time_variable, survival_variable_all, survival_trans_function)[[k]][, 1]
      } else rep(0, length(mids))
    } else {
      log_f <- if (all(has_rows)) {
        conditionalYT_fun(data_j, long_fit_all, l_i = mids, survival_variable, time_variable,
                          survival_variable_all, survival_trans_function)[[1]][, 1]
      } else rep(0, length(mids))
    }
    w <- exp(log_f - max(log_f[is.finite(log_f)])) * pmax(mass, 0)
    w[!is.finite(w)] <- 0
    k <- if (sum(w) > 0) {
      sample.int(length(w), n_draws, replace = TRUE, prob = w)
    } else {
      sample.int(length(w), n_draws, replace = TRUE)
    }
    out[j, ] <- stats::runif(n_draws, edges[k], edges[k + 1])
  }
  if (n_draws == 1) stats::setNames(out[, 1], bounds$id) else out
}

#' Interval endpoints per subject
#' @keywords internal
subject_intervals <- function(data_survival, form_marginal_surv, id, event_type = NULL) {
  assert_vars_in_data(c(id, event_type), data_survival, "the id and event-type variables", "data_survival")
  y <- eval(form_marginal_surv[[2]], data_survival, environment(form_marginal_surv))
  b <- interval_bounds(y)
  keep <- !is.na(b$L)
  out <- data.frame(id = as.character(data_survival[[id]])[keep], L = b$L[keep], R = b$R[keep],
                    stringsAsFactors = FALSE)
  if (!is.null(event_type)) out$D <- data_survival[[event_type]][keep]
  out
}

#' Write event times (and their transforms) into longitudinal data
#' @keywords internal
fill_event_time <- function(data_fit_all, id, T_values, event_time,
                            survival_variable_all = NULL, survival_trans_function = NULL) {
  lapply(data_fit_all, function(d) {
    d[[event_time]] <- unname(T_values[as.character(d[[id]])])
    for (k in seq_along(survival_variable_all)) {
      d[[survival_variable_all[[k]]]] <- vapply(d[[event_time]], survival_trans_function[[k]], numeric(1))
    }
    d
  })
}
