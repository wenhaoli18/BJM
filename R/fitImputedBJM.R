#' Fit the backward joint model with censored subjects in the longitudinal sub-model
#'
#' @description \strong{Experimental.} The usual right-censored fit estimates
#' the longitudinal sub-model \eqn{f(Y \mid T)} on the subjects with an
#' observed event only. That is valid when censoring is independent of
#' \eqn{(T, Y)}, but \eqn{f(Y \mid T)} is then learned only over the event
#' times seen among the events, which are the early ones when many subjects
#' are censored. Prediction needs \eqn{f(Y \mid T)} for every \eqn{T} after
#' the prediction time, so for a patient doing well it extrapolates the
#' event-time effect beyond the data. With a large censored fraction this
#' shifts the predictive distribution of \eqn{T} towards early events: risks
#' are overestimated, increasingly so at later landmarks, and in short
#' follow-up the ranking of patients can be lost too.
#'
#' \code{fitImputedBJM()} also uses the censored subjects, who carry the
#' information on long event times. It fills in their \eqn{T > C} by
#' stochastic EM, then by multiple imputation, as \code{\link{fitIntervalBJM}}
#' does for interval censoring:
#' \enumerate{
#'   \item fit the survival sub-model (\code{survivalSub()}) once, and the
#'   longitudinal sub-model on the subjects with an event (the usual fit);
#'   \item draw each censored subject's \eqn{T} from
#'   \eqn{f(T \mid T > C, Y, X) \propto f(Y \mid T, X) f(T \mid X)}
#'   (\code{imputeEventTime()}); with competing risks, draw the event type
#'   \eqn{D} jointly, from \eqn{f(Y \mid T, D) P(D \mid T) f(T)};
#'   \item refit \code{longitudinalSub()} on all subjects and repeat for
#'   \code{n_burnin} iterations;
#'   \item keep iterating, saving every \code{n_thin}-th fit until
#'   \code{n_imputations} are saved, and pool them for prediction (see
#'   \code{\link{poolLongitudinalSub}}).
#' }
#' The survival sub-model, including the event-type model, is the usual one
#' and is not refit. The draws of \eqn{T} beyond the last follow-up rely on
#' the survival model's extrapolated tail.
#'
#' Using the censored subjects makes the event-time effect matter where it
#' was never observed, so a term linear in \eqn{T} can drift over the long
#' event times now in the fit. A bounded effect, e.g.
#' \code{I(exp(-years / 5))} in place of \code{years} in every term of
#' \code{long_sub_fixed}, levels off for long \eqn{T}; choose its scale by
#' comparing \code{\link{observedLongLogLik}} across candidate fits. On
#' \code{pbc3} (landmarks 1-3 years, 3-year window), the usual fit
#' overestimates the cumulative incidence of death by 8-32\% in 5-fold
#' cross-validation and this fit with \code{exp(-years / 5)} by -1 to 7\%,
#' with the same AUC.
#'
#' @param data_survival One row per subject: event or censoring time, event
#'   indicator, the survival covariates and the subject id.
#' @param data_fit_all As in \code{\link{longitudinalSub}}: one data frame
#'   per biomarker (or one shared data frame), long format, for
#'   \strong{every} subject, with the event (or censoring) time in the
#'   column used in \code{long_sub_fixed}. The censored subjects' event-time
#'   column, event-type column and the columns named in
#'   \code{survival_variable_all} are (re)filled from the draws.
#' @param form_marginal_surv,form_conditional_cr As in
#'   \code{\link{survivalSub}}: \code{Surv(time, status) ~ covariates} and
#'   the optional event-type model. The event-time variable of
#'   \code{form_marginal_surv} must be the one used in \code{long_sub_fixed}.
#' @param long_sub_fixed,long_sub_random As in \code{\link{longitudinalSub}}.
#' @param time_variable Name of the visit-time column.
#' @param survival_variable_all,survival_trans_function Transformed event-time
#'   columns, as in \code{\link{predictRisk}} (\code{NULL} if none).
#' @param n_burnin Stochastic-EM iterations before fits are saved.
#' @param n_imputations Number of saved fits.
#' @param n_thin Iterations between saved fits.
#' @param n_grid Grid points per subject for the draws.
#' @param baseline,df Baseline hazard of the survival model (see
#'   \code{\link{survivalSub}}).
#' @param seed Optional integer seed (the global RNG state is restored).
#'
#' @return An object of class \code{"fitImputedBJM"}, a list with
#'   \code{survival_fit_all}; \code{long_fit_all} (the imputations combined,
#'   for \code{predictRisk()} and friends); \code{long_fit_all_list} (one fit
#'   per imputation); \code{pooled} (\code{\link{poolLongitudinalSub}} of the
#'   fits, when \code{n_imputations > 1}); \code{imputed_T} (and, with
#'   competing risks, \code{imputed_D}), censored subjects x imputations;
#'   \code{trace} (each iteration's fixed effects, to check the chain has
#'   settled); and \code{settings}. \code{print()} shows a short summary.
#' @seealso \code{\link{observedLongLogLik}} to compare fits;
#'   \code{\link{fitIntervalBJM}} for interval-censored event times.
#' @examples
#' \donttest{
#' data(pbc3)
#' pbc3_surv <- pbc3[!duplicated(pbc3$id), ]
#' fit <- fitImputedBJM(pbc3_surv, pbc3, Surv(years, status3) ~ age + sex,
#'                      long_sub_fixed = serBilir ~ year + age + sex +
#'                        I(exp(-years / 5)) + I(year * exp(-years / 5)),
#'                      long_sub_random = ~ year | id, time_variable = "year",
#'                      n_burnin = 2, n_imputations = 2, n_thin = 1, seed = 1)
#' fit
#'
#' # compare with the usual fit on the subjects with an event
#' usual <- longitudinalSub(pbc3[pbc3$status3 == 1, ],
#'                          list(serBilir ~ year + age + sex + years + years * year),
#'                          list(~ year | id))
#' observedLongLogLik(fit, pbc3_surv, pbc3, time_variable = "year")
#' observedLongLogLik(list(long_fit_all = usual, survival_fit_all = fit$survival_fit_all),
#'                    pbc3_surv, pbc3, time_variable = "year")
#'
#' history <- pbc3[pbc3$id == 2 & pbc3$year <= 3, ]
#' predictRisk(history, fit$long_fit_all, fit$survival_fit_all,
#'             prediction_time = 3, horizon = 3, time_variable = "year",
#'             survival_variable_all = NULL, survival_trans_function = NULL,
#'             bandcount1 = 10, bandcount2 = 10)
#' }
#' @export
fitImputedBJM <- function(data_survival, data_fit_all, form_marginal_surv,
                          long_sub_fixed, long_sub_random, time_variable,
                          form_conditional_cr = NULL,
                          survival_variable_all = NULL, survival_trans_function = NULL,
                          n_burnin = 5, n_imputations = 5, n_thin = 2, n_grid = 40,
                          baseline = c("spline", "piecewise"), df = NULL, seed = NULL) {
  assert_positive_integer(n_thin, "n_thin")
  assert_positive_integer(n_imputations, "n_imputations")
  if (!is.null(seed)) local_r_seed(seed)
  if (inherits(long_sub_fixed, "formula")) long_sub_fixed <- list(long_sub_fixed)
  if (inherits(long_sub_random, "formula")) long_sub_random <- list(long_sub_random)
  M <- length(long_sub_fixed)
  if (is.data.frame(data_fit_all)) data_fit_all <- rep(list(data_fit_all), M)
  id <- as.character(nlme::splitFormula(long_sub_random[[1]], "|")[[2]])[2]

  y <- eval(form_marginal_surv[[2]], data_survival, environment(form_marginal_surv))
  if (!inherits(y, "Surv") || attr(y, "type") != "right") {
    stop("`form_marginal_surv` must be Surv(time, status) ~ ... (right-censored); for interval censoring use fitIntervalBJM().",
         call. = FALSE)
  }
  survival_fit_all <- survivalSub(data_survival, form_marginal_surv, form_conditional_cr,
                                  baseline = baseline, df = df)
  event_time <- survival_time_variable(survival_fit_all)
  has_cr <- length(form_conditional_cr) != 0
  event_type <- if (has_cr) all.vars(form_conditional_cr[[2]])[1] else NULL
  ids <- as.character(data_survival[[id]])
  status <- y[, "status"]
  keep <- !is.na(status) & !is.na(ids)
  event_ids <- ids[keep & status == 1]
  censored <- keep & status == 0
  if (!any(censored)) stop("No censored subjects: the usual fit already uses everyone.", call. = FALSE)
  bounds <- data.frame(id = ids[censored], L = y[censored, "time"], R = Inf, stringsAsFactors = FALSE)
  if (has_cr) bounds$D <- NA_real_

  data_event <- lapply(data_fit_all, function(d) d[as.character(d[[id]]) %in% event_ids, , drop = FALSE])
  data_cens <- lapply(data_fit_all, function(d) d[as.character(d[[id]]) %in% bounds$id, , drop = FALSE])
  # the event type is set inside the conditional densities; a placeholder
  # keeps the censored rows from being dropped as incomplete before the
  # first draw
  if (has_cr) data_cens <- lapply(data_cens, function(d) { d[[event_type]] <- 0; d })

  fill_cens <- function(T_values, D_values) {
    filled <- fill_event_time(data_cens, id, T_values, event_time,
                              survival_variable_all, survival_trans_function)
    if (has_cr) filled <- lapply(filled, function(d) {
      d[[event_type]] <- unname(D_values[as.character(d[[id]])])
      d
    })
    Map(rbind, data_event, filled)
  }
  # burn-in fits are thrown away, so their Sigma_fit EM convergence
  # warnings are muffled; the saved fits keep theirs
  fit_long <- function(data, quiet) {
    withCallingHandlers(longitudinalSub(data, long_sub_fixed, long_sub_random),
                        warning = function(w) {
                          if (quiet && grepl("Sigma_fit", conditionMessage(w))) invokeRestart("muffleWarning")
                        })
  }
  long_fit_all <- fit_long(data_event, quiet = TRUE)

  n_iter <- n_burnin + n_imputations * n_thin
  saved_index <- function(it) {
    k <- (it - n_burnin) / n_thin
    if (it > n_burnin && k == round(k)) k else NA_integer_
  }
  trace <- vector("list", n_iter)
  imputed_T <- matrix(NA_real_, nrow(bounds), n_imputations, dimnames = list(bounds$id, NULL))
  imputed_D <- if (has_cr) imputed_T
  long_fit_all_list <- vector("list", n_imputations)
  for (it in seq_len(n_iter)) {
    # an unlucky draw occasionally makes nlme::lme() fail (e.g. "Singularity
    # in backsolve"); draw again from the same fit rather than abort the chain
    for (attempt in seq_len(5)) {
      T_new <- imputeEventTime(data_cens, long_fit_all, survival_fit_all, bounds, time_variable,
                               survival_variable_all, survival_trans_function, n_grid = n_grid)
      D_new <- attr(T_new, "event_type")
      new_fit <- tryCatch(fit_long(fill_cens(T_new, D_new), quiet = is.na(saved_index(it))),
                          error = function(e) e)
      if (!inherits(new_fit, "error")) break
    }
    if (inherits(new_fit, "error")) {
      stop(sprintf("longitudinalSub() failed on 5 successive imputations at iteration %d: %s",
                   it, conditionMessage(new_fit)), call. = FALSE)
    }
    long_fit_all <- new_fit
    trace[[it]] <- unlist(lapply(long_fit_all$lfit, nlme::fixef))
    k <- saved_index(it)
    if (!is.na(k)) {
      imputed_T[, k] <- T_new[bounds$id]
      if (has_cr) imputed_D[, k] <- D_new[bounds$id]
      long_fit_all_list[[k]] <- long_fit_all
    }
  }

  pooled <- if (n_imputations > 1) poolLongitudinalSub(long_fit_all_list)
  combined <- if (n_imputations > 1) pooled$long_fit_all else long_fit_all_list[[1]]
  combined$imputed_censored <- TRUE
  out <- list(survival_fit_all = survival_fit_all, long_fit_all = combined,
              long_fit_all_list = long_fit_all_list, pooled = pooled,
              imputed_T = imputed_T, imputed_D = imputed_D, trace = do.call(rbind, trace),
              settings = list(n_burnin = n_burnin, n_imputations = n_imputations, n_thin = n_thin,
                              n_events = length(event_ids), n_censored = nrow(bounds)))
  class(out) <- "fitImputedBJM"
  out
}

#' Observed-data log-likelihood of the longitudinal sub-model
#'
#' @description \strong{Experimental.} Compares longitudinal sub-models
#' fitted to the same right-censored data, e.g. different forms of the
#' event-time effect in \code{long_sub_fixed}, or a fit on the subjects with
#' an event against \code{\link{fitImputedBJM}}. A subject with an event
#' contributes \eqn{\log f(Y \mid T, D)} at the observed \eqn{T} (and event
#' type); a censored subject contributes
#' \eqn{\log \int_C^\infty \sum_d f(Y \mid t, d) P(d \mid t) \, dF(t) -
#' \log S(C)}, so the censored subjects count too. The factors
#' \eqn{f(T)} and \eqn{P(D \mid T)} of the events do not involve the
#' longitudinal sub-model and are left out, so only differences between fits
#' sharing \code{survival_fit_all} are meaningful. Higher is better; the fits
#' compared should have the same number of parameters, or the difference be
#' judged with that in mind.
#'
#' @param fit A \code{"fitImputedBJM"} object (its log-likelihood is averaged
#'   over the imputations' fits), or a list with \code{long_fit_all} and
#'   \code{survival_fit_all}.
#' @param data_survival,data_fit_all The data the fit was made on, as in
#'   \code{\link{fitImputedBJM}}: every subject, with the observed event or
#'   censoring time.
#' @param time_variable,survival_variable_all,survival_trans_function As in
#'   \code{\link{predictRisk}}.
#' @param n_grid Grid points per censored subject.
#' @return The log-likelihood (a number).
#' @export
observedLongLogLik <- function(fit, data_survival, data_fit_all, time_variable,
                               survival_variable_all = NULL, survival_trans_function = NULL,
                               n_grid = 40) {
  survival_fit_all <- fit$survival_fit_all
  fits <- if (inherits(fit, "fitImputedBJM")) fit$long_fit_all_list else list(fit$long_fit_all)
  assert_class(survival_fit_all, "survivalSub.BJM", "fit$survival_fit_all", "survivalSub")
  if (is_interval_censored(survival_fit_all)) {
    stop("observedLongLogLik() supports right-censored survival fits only.", call. = FALSE)
  }
  form <- survival_fit_all$form_marginal_surv
  y <- eval(form[[2]], data_survival, environment(form))
  event_time <- survival_time_variable(survival_fit_all)
  has_cr <- !is.null(survival_fit_all$glm_fit)
  event_type <- if (has_cr) all.vars(survival_fit_all$form_conditional_cr[[2]])[1]
  mean(vapply(fits, function(long_fit_all) {
    id <- as.character(nlme::splitFormula(long_fit_all$long_sub_random[[1]], "|")[[2]])[2]
    M <- length(long_fit_all$lfit)
    data_all <- if (is.data.frame(data_fit_all)) rep(list(data_fit_all), M) else data_fit_all
    ids <- as.character(data_survival[[id]])
    if (has_cr) data_all <- lapply(data_all, function(d) {
      d[[event_type]][is.na(d[[event_type]])] <- 0
      d
    })
    data_all <- align_ordinal_levels(data_all, long_fit_all)
    data_all <- suppressWarnings(drop_missing_longitudinal(data_all, long_fit_all, event_time,
                                                           survival_variable_all))
    ll <- 0
    for (j in seq_along(ids)) {
      if (is.na(y[j, "status"])) next
      data_j <- lapply(data_all, function(d) d[as.character(d[[id]]) == ids[j], , drop = FALSE])
      if (!all(vapply(data_j, nrow, integer(1)) > 0)) next
      if (y[j, "status"] == 1) {
        log_f <- if (has_cr) {
          conditionalYDT(data_j, long_fit_all, survival_fit_all, l_i = y[j, "time"], event_time,
                         time_variable, survival_variable_all,
                         survival_trans_function)[[data_survival[[event_type]][j] + 1]][1, 1]
        } else {
          conditionalYT(data_j, long_fit_all, l_i = y[j, "time"], event_time, time_variable,
                        survival_variable_all, survival_trans_function)[[1]][1, 1]
        }
        ll <- ll + log_f
      } else {
        tw <- subject_time_weights(data_j, TRUE, long_fit_all, survival_fit_all, y[j, "time"], Inf,
                                   NA, time_variable, survival_variable_all, survival_trans_function,
                                   n_grid)
        mx <- max(tw$log_w[is.finite(tw$log_w)])
        ll <- ll + mx + log(sum(exp(tw$log_w - mx))) - log(sum(pmax(tw$mass, 0)))
      }
    }
    ll
  }, numeric(1)))
}

#' Print a fitImputedBJM fit
#'
#' @description A short summary of a \code{\link{fitImputedBJM}} fit: the
#' imputation settings and the longitudinal sub-model's fixed effects
#' (pooled over the imputations, with the fraction of missing information).
#' @param x A \code{"fitImputedBJM"} object.
#' @param digits Number of significant digits.
#' @param ... Not used.
#' @return \code{x}, invisibly.
#' @keywords internal
#' @export
print.fitImputedBJM <- function(x, digits = 4, ...) {
  st <- x$settings
  cat("\nBackward joint model with imputed censored event times (fitImputedBJM)\n\n")
  cat(sprintf("  %d subjects with an event, %d censored with T imputed; %d imputation(s) after %d burn-in iterations, keeping %s\n",
              st$n_events, st$n_censored, st$n_imputations, st$n_burnin,
              if (st$n_thin == 1) "every iteration" else sprintf("one in every %d iterations", st$n_thin)))
  if (!is.null(x$pooled)) {
    cat("\nLongitudinal sub-model, pooled over the imputations (Rubin's rules):\n")
    for (k in seq_along(x$pooled$pooled)) {
      tab <- x$pooled$pooled[[k]][, c("estimate", "std_error", "fmi")]
      cat(sprintf("\n  [%d] %s\n", k, names(x$pooled$pooled)[k]))
      print(signif(as.matrix(tab), digits))
    }
  } else {
    cat("\nLongitudinal sub-model (one imputation), fixed effects:\n")
    for (k in seq_along(x$long_fit_all$lfit)) print(signif(nlme::fixef(x$long_fit_all$lfit[[k]]), digits))
  }
  cat("\nUse x$long_fit_all and x$survival_fit_all for prediction; x$trace to check the chain has settled.\n")
  invisible(x)
}

bjm_state <- new.env(parent = emptyenv())

#' Warn once when the longitudinal sub-model ignores many censored subjects
#'
#' @description Called by \code{predictRisk()}: if the longitudinal fit has
#' no more subjects than the survival fit has events while at least 30\% of
#' the survival subjects are censored, it was most likely fit on the
#' subjects with an event only, and risks can be overestimated (see
#' \code{\link{fitImputedBJM}}). Warns once per session; set
#' \code{options(BJM.check_complete_case = FALSE)} to silence it.
#' @keywords internal
check_complete_case_fit <- function(long_fit_all, survival_fit_all) {
  if (!isTRUE(getOption("BJM.check_complete_case", TRUE)) || isTRUE(bjm_state$warned_complete_case) ||
      isTRUE(long_fit_all$imputed_censored) || is_interval_censored(survival_fit_all) ||
      is.null(survival_fit_all$coxph_fit)) return(invisible())
  n <- survival_fit_all$coxph_fit$n
  n_event <- survival_fit_all$coxph_fit$nevent
  n_long <- unname(long_fit_all$lfit[[1]]$dims$ngrps[1])
  if (is.null(n) || is.null(n_event) || is.null(n_long)) return(invisible())
  if (n_long <= n_event && (n - n_event) / n >= 0.3) {
    bjm_state$warned_complete_case <- TRUE
    warning(sprintf(paste0(
      "The longitudinal sub-model has %d subjects, no more than the %d with an event, while %.0f%% of the ",
      "%d subjects in the survival sub-model are censored. Fitted on the events only, f(Y | T) extrapolates ",
      "the event-time effect to the long event times of the censored, which can overestimate risk. ",
      "fitImputedBJM() uses the censored subjects too. (Shown once per session; ",
      "options(BJM.check_complete_case = FALSE) turns it off.)"),
      n_long, n_event, 100 * (n - n_event) / n, n), call. = FALSE)
  }
  invisible()
}
