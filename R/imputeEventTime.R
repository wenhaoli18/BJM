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
#'   \item keep iterating, saving every \code{n_thin}-th draw and its
#'   \code{longitudinalSub()} fit until \code{n_imputations} are saved.
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
#' The saved draws come from one chain, so successive draws are correlated;
#' saving only every \code{n_thin}-th iteration makes them closer to
#' independent. They are drawn with the model parameters at their point
#' estimates, not from their posterior, so the between-imputation variance
#' in \code{pooled} can still be somewhat understated.
#'
#' For prediction, \code{long_fit_all} combines the imputations: its fixed
#' effects are the Rubin's-rules pooled estimates and its residual and
#' random-effects variances the averages over the imputations (see
#' \code{\link{poolLongitudinalSub}}), so a prediction does not depend on
#' which single imputation happened to come first. With competing risks the
#' event-type model of \code{survival_fit_all} likewise has its
#' coefficients averaged over the imputations.
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
#' @param n_thin Iterations between saved draws (1 saves every iteration
#'   after the burn-in). Larger values make the saved draws less correlated
#'   at the cost of \code{n_burnin + n_imputations * n_thin} iterations in
#'   all.
#' @param n_grid Grid points per subject interval.
#' @param baseline,df Baseline hazard of the survival model (see
#'   \code{\link{survivalSub}}).
#' @param include_right_censored See Description.
#' @param seed Optional integer seed (the global RNG state is restored).
#'
#' @return An object of class \code{"fitIntervalBJM"}, a list with
#'   \code{survival_fit_all} (with competing risks, its event-type model has
#'   the coefficients averaged over the imputations, while the rest of that
#'   \code{glm} object, e.g. its standard errors, is the first imputation's
#'   fit; \code{survival_fit_all_list} has one fit per imputation);
#'   \code{long_fit_all_list} (one \code{longitudinalSub()} fit per
#'   imputation); \code{long_fit_all} (the imputations combined, see
#'   Description, for use with \code{predictRisk()} and friends);
#'   \code{pooled} (\code{\link{poolLongitudinalSub}} of the fits, when
#'   \code{n_imputations > 1}); \code{imputed_T} (a subjects x imputations
#'   matrix); and \code{trace} (each iteration's fixed effects, to check
#'   the chain has settled).
#' @examples
#' \donttest{
#' # simulated study: visits every 1-3 years, an event known only to lie
#' # between the last visit without it (L) and the first visit with it (R)
#' set.seed(1)
#' n <- 60
#' x <- rnorm(n)
#' event_time <- (-log(runif(n)) / (0.05 * exp(0.5 * x)))^(1 / 1.5)
#' surv <- data.frame(id = seq_len(n), x = x, L = NA_real_, R = NA_real_)
#' long <- list()
#' for (i in seq_len(n)) {
#'   visits <- c(0, cumsum(runif(30, 1, 3)))
#'   visits <- visits[visits < 12]
#'   if (event_time[i] > max(visits)) {
#'     surv$L[i] <- max(visits)
#'     seen <- visits
#'   } else {
#'     k <- findInterval(event_time[i], visits)
#'     surv$L[i] <- visits[k]
#'     surv$R[i] <- visits[k + 1]
#'     seen <- visits[visits <= event_time[i]]
#'   }
#'   long[[i]] <- data.frame(id = i, year = seen, x = x[i],
#'                           y = 1 + 0.3 * seen - 0.2 * event_time[i] + 0.4 * x[i] +
#'                             rnorm(1, 0, 0.5) + rnorm(length(seen), 0, 0.3))
#' }
#' long <- do.call(rbind, long)
#'
#' fit <- fitIntervalBJM(surv, long, Surv(L, R, type = "interval2") ~ x,
#'                       event_time = "Tev", long_sub_fixed = y ~ year + Tev + x,
#'                       long_sub_random = ~ year | id, time_variable = "year",
#'                       n_burnin = 2, n_imputations = 2, n_thin = 1, n_grid = 15,
#'                       seed = 1)
#' fit$pooled
#'
#' # the fit plugs into the prediction functions; Tev is unknown for a new patient
#' history <- long[long$id == 1 & long$year <= 2, ]
#' history$Tev <- NA_real_
#' predictRisk(history, fit$long_fit_all, fit$survival_fit_all,
#'             prediction_time = 2, horizon = 2, time_variable = "year",
#'             survival_variable_all = NULL, survival_trans_function = NULL,
#'             bandcount1 = 10, bandcount2 = 20)
#' }
#' @export
fitIntervalBJM <- function(data_survival, data_fit_all, form_marginal_surv, event_time,
                           long_sub_fixed, long_sub_random, time_variable,
                           survival_variable_all = NULL, survival_trans_function = NULL,
                           form_conditional_cr = NULL,
                           n_burnin = 10, n_imputations = 5, n_thin = 5, n_grid = 40,
                           baseline = c("spline", "piecewise"), df = NULL,
                           include_right_censored = FALSE, seed = NULL) {
  assert_positive_integer(n_thin, "n_thin")
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

  # after the burn-in, every n_thin-th iteration is saved
  n_iter <- n_burnin + n_imputations * n_thin
  saved_index <- function(it) {
    k <- (it - n_burnin) / n_thin
    if (it > n_burnin && k == round(k)) k else NA_integer_
  }
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
      new_fit <- tryCatch(fit_long(T_new, quiet = is.na(saved_index(it))), error = function(e) e)
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
    k <- saved_index(it)
    if (!is.na(k)) {
      imputed_T[, k] <- T_cur[bounds$id]
      long_fit_all_list[[k]] <- long_fit_all
      survival_fit_all_list[[k]] <- survival_fit_all
    }
  }

  # predictions use the imputations combined, not whichever came first
  pooled <- if (n_imputations > 1) poolLongitudinalSub(long_fit_all_list)
  out <- list(survival_fit_all = pool_event_type_fit(survival_fit_all_list),
              long_fit_all_list = long_fit_all_list,
              survival_fit_all_list = survival_fit_all_list,
              long_fit_all = if (n_imputations > 1) pooled$long_fit_all else long_fit_all_list[[1]],
              pooled = pooled,
              imputed_T = imputed_T, trace = do.call(rbind, trace))
  class(out) <- "fitIntervalBJM"
  out
}

#' Average the event-type model over the imputations
#'
#' @description With competing risks \code{fitIntervalBJM()} refits the
#' event-type \code{glm} on every imputation, since it uses the event time.
#' For prediction the coefficients are averaged over the imputations (the
#' Rubin's-rules point estimate); the rest of the returned \code{glm} object
#' is the first imputation's fit.
#'
#' @param survival_fit_all_list One \code{survivalSub.BJM} fit per imputation.
#' @return The first fit, with the averaged event-type coefficients.
#' @keywords internal
pool_event_type_fit <- function(survival_fit_all_list) {
  out <- survival_fit_all_list[[1]]
  if (is.null(out$glm_fit) || length(survival_fit_all_list) == 1) return(out)
  coefs <- vapply(survival_fit_all_list, function(f) stats::coef(f$glm_fit),
                  numeric(length(stats::coef(out$glm_fit))))
  out$glm_fit$coefficients <- rowMeans(matrix(coefs, nrow = length(stats::coef(out$glm_fit))))
  names(out$glm_fit$coefficients) <- names(stats::coef(survival_fit_all_list[[1]]$glm_fit))
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
