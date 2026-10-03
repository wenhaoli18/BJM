#' Simulate future trajectories from the backward joint model
#'
#' @description Draws \code{n_sim} complete futures for each patient from the
#' fitted backward joint model, conditional on the patient's longitudinal
#' history up to \code{prediction_time} and on being event-free at
#' \code{prediction_time}: an event time (and, with competing risks, an
#' event type), and the values of every biomarker at the requested
#' \code{times}. Where \code{\link{predictRisk}} and
#' \code{\link{predictLongitudinal}} summarize the predictive distribution,
#' this returns samples from it, so it can be used to generate synthetic
#' patients or to propagate the prediction into any other quantity (e.g. a
#' patient's "digital twin" under the fitted model).
#'
#' Each draw follows the model's factorization
#' \eqn{f(T, D)\, f(Y \mid T, D)}:
#' \enumerate{
#'   \item \eqn{(T, D)} is drawn from its posterior given the history,
#'   proportional to \eqn{f(Y_{hist} \mid T, D)\, f(D \mid T)\, f(T)}
#'   restricted to \eqn{T > } \code{prediction_time}, evaluated on the same
#'   integration grid as \code{\link{predictRisk}}'s denominator;
#'   \item the random effects \eqn{b} are drawn from their Gaussian
#'   posterior given the history and the drawn \eqn{(T, D)};
#'   \item the biomarkers at \code{times} are drawn from the longitudinal
#'   sub-model given \eqn{(T, D, b)}, adding residual error.
#' }
#'
#' A patient with no non-missing biomarker measurement up to
#' \code{prediction_time} (e.g. a single baseline row with the biomarkers set
#' to \code{NA}) is simulated from the model's prior given their baseline
#' covariates, i.e. as a new synthetic patient.
#'
#' @details
#' The event time is drawn on the grid of \code{bandcount2} intervals that
#' \code{\link{predictRisk}} integrates over (each holding about the same
#' share of the conditional survival probability) and reported as the
#' midpoint of its interval, so its resolution is about
#' \code{1 / bandcount2} of the conditional distribution.
#'
#' The biomarkers at \code{times} are evaluated with the covariates of the
#' patient's last row in \code{data_predict_all} (with \code{time_variable}
#' set to each element of \code{times}), so a time-varying covariate is
#' carried forward at its last value.
#'
#' The backward model describes the biomarkers before the event. With
#' \code{truncate = TRUE} (the default) a value at a time at or after the
#' drawn event time is \code{NA}; with \code{truncate = FALSE} it is still
#' drawn from the sub-model, which is the quantity whose density
#' \code{\link{predictLongitudinal}} returns.
#'
#' Only fits with continuous biomarkers are supported; a fit with an ordinal
#' biomarker (see \code{\link{longitudinalSubCopula}}) gives an error.
#'
#' @inheritParams predictRisk
#' @param times Numeric vector of times (on the scale of
#'   \code{time_variable}) at which to simulate the biomarkers, each at least
#'   \code{prediction_time}.
#' @param n_sim Number of futures to draw per patient.
#' @param bandcount2 The number of intervals the event-time distribution
#'   after \code{prediction_time} is discretized into (see Details and
#'   \code{\link{predictRisk}}).
#' @param truncate If \code{TRUE}, biomarker values at or after the drawn
#'   event time are \code{NA} (see Details).
#' @param seed Optional integer seed, for reproducible draws. The caller's
#'   random number stream is left unchanged.
#'
#' @return A data frame of class \code{"simulateTrajectory.BJM"}, with one
#' row per patient, draw and element of \code{times}, and columns: the
#' patient id (named as in \code{long_sub_random}); \code{sim}, the draw
#' number; \code{time_variable}; \code{event_time}, the drawn event time;
#' with competing risks, the drawn event type (named as the response of
#' \code{form_conditional_cr}, with the same 0/1 coding); and one column per
#' biomarker, named by its response variable. Attributes
#' \code{"prediction_time"} and \code{"times"} record the call's values.
#'
#' @examples
#' \donttest{
#' data(pbc3)
#'
#' survival_fit_all = survivalSub(pbc3[!duplicated(pbc3$id), ],
#'                                Surv(years, status3) ~ age + sex, NULL)
#'
#' long_sub_fixed = list(
#'   "long1" = serBilir ~ year + age + sex + (years) + (years) * year,
#'   "long2" = albumin ~ year + age + sex + (years) + (years) * year)
#' long_sub_random = list("long1" = ~ year | id, "long2" = ~ year | id)
#' data_fit_all = list(pbc3[pbc3$status3 == 1, ], pbc3[pbc3$status3 == 1, ])
#' long_fit_all = longitudinalSub(data_fit_all, long_sub_fixed, long_sub_random)
#'
#' trans = survivalTrans(c(1, 3, 5, 7))
#'
#' # Patient 2's history up to year 3: 200 possible futures over years 3-8
#' history = pbc3[pbc3$id == 2 & pbc3$year <= 3, ]
#' sims = simulateTrajectory(list(history, history), long_fit_all, survival_fit_all,
#'                           prediction_time = 3, times = seq(3, 8, by = 0.5),
#'                           time_variable = "year", trans$survival_variable_all,
#'                           trans$survival_trans_function, n_sim = 200, seed = 1)
#' head(sims)
#'
#' # Probability of an event within 2 years, from the draws
#' mean(sims$event_time[!duplicated(sims$sim)] <= 5)
#' }
#'
#' @export
simulateTrajectory <- function(data_predict_all, long_fit_all, survival_fit_all,
                               prediction_time, times, time_variable,
                               survival_variable_all, survival_trans_function,
                               n_sim = 100, bandcount2 = 100, truncate = TRUE, seed = NULL) {

  assert_class(long_fit_all, "longitudinalSub.BJM", "long_fit_all", "longitudinalSub")
  assert_class(survival_fit_all, "survivalSub.BJM", "survival_fit_all", "survivalSub")
  if (!is.null(long_fit_all$biomarker_type) && any(long_fit_all$biomarker_type == "ordinal")) {
    stop("simulateTrajectory() supports only fits with continuous biomarkers; this fit has an ordinal biomarker.",
         call. = FALSE)
  }
  M <- length(long_fit_all$lfit)
  assert_data_list(data_predict_all, "data_predict_all", M, allow_bare_df = TRUE)
  if (!is.list(data_predict_all) || is.data.frame(data_predict_all)) {
    data_predict_all <- rep(list(data_predict_all), each = M)
  }
  assert_scalar_numeric(prediction_time, "prediction_time")
  if (!is.numeric(times) || length(times) == 0 || any(!is.finite(times))) {
    stop("`times` must be a non-empty numeric vector of finite times.", call. = FALSE)
  }
  if (any(times < prediction_time)) {
    stop(sprintf("Every element of `times` must be at least prediction_time = %g.", prediction_time),
         call. = FALSE)
  }
  times <- sort(unique(times))
  assert_string(time_variable, "time_variable")
  assert_positive_integer(n_sim, "n_sim")
  assert_positive_integer(bandcount2, "bandcount2")
  if (!is.logical(truncate) || length(truncate) != 1 || is.na(truncate)) {
    stop("`truncate` must be TRUE or FALSE.", call. = FALSE)
  }
  assert_survival_trans(survival_variable_all, survival_trans_function, probe_value = prediction_time)
  if (!is.null(seed)) local_r_seed(seed)

  data_predict_all <- drop_after_prediction_time(data_predict_all, time_variable, prediction_time)

  lfit <- long_fit_all$lfit
  Sigma <- long_fit_all$Sigma_fit
  id <- as.character(nlme::splitFormula(long_fit_all$long_sub_random[[1]], "|")[[2]])[2]
  survival_variable <- as.character(formula(survival_fit_all$coxph_fit)[[2]])[2]
  has_cr <- length(survival_fit_all$form_conditional_cr) != 0
  event_type_variable <- if (has_cr) all.vars(survival_fit_all$form_conditional_cr[[2]]) else NULL

  ### the outcome columns are simulated, so a new patient need not have
  ### them; conditionalDT() needs them to exist
  for (i in seq_len(M)) {
    for (v in c(survival_variable, event_type_variable)) {
      if (!v %in% names(data_predict_all[[i]])) data_predict_all[[i]][[v]] <- NA_real_
    }
  }
  data_predict_all <- subset_at_risk(data_predict_all, survival_variable, prediction_time)
  patient_ids <- prediction_patient_ids(data_predict_all, long_fit_all)

  ### event-time grid and its prior weights, shared with predictRisk()
  upper_bound <- integration_upper_bound(data_predict_all, long_fit_all, survival_fit_all,
                                         prediction_time, min_upper = max(times))
  infinity_grid <- prepare_infinity_grid(data_predict_all, long_fit_all, survival_fit_all,
                                          prediction_time, upper_bound, bandcount2)
  l_grid <- infinity_grid$predict.time.infinity
  log_prior <- list(log(pmax(infinity_grid$S_T_all_infinity, 0)))
  d_values <- NA
  if (has_cr) {
    D_T <- conditionalDT(data_predict_all, long_fit_all, survival_fit_all, l_i = l_grid)
    log_prior <- list(log_prior[[1]] + log(D_T[[1]]), log_prior[[1]] + log(D_T[[2]]))
    d_values <- c(0, 1)
  }

  sigma2 <- vapply(lfit, function(f) f$sigma^2, numeric(1))
  response_names <- vapply(lfit, function(f) as.character(formula(f)[[2]]), character(1))
  n_times <- length(times)

  out <- vector("list", length(patient_ids))
  for (p in seq_along(patient_ids)) {
    rows_p <- lapply(data_predict_all, function(d) d[as.character(d[[id]]) == patient_ids[p], , drop = FALSE])
    history <- lapply(seq_len(M), function(i) complete_history_rows(rows_p[[i]], i, long_fit_all, survival_variable,
                                                                     c(event_type_variable, unlist(survival_variable_all))))
    template <- rows_p[[1]][nrow(rows_p[[1]]), , drop = FALSE]
    future <- lapply(seq_len(M), function(i) {
      d <- template[rep(1, n_times), , drop = FALSE]
      d[[time_variable]] <- times
      d
    })

    y_hist <- unlist(lapply(seq_len(M), function(i) history[[i]][[response_names[i]]]))
    Z_hist <- stacked_random_design(history, long_fit_all$long_sub_random)
    Z_fut <- stacked_random_design(future, long_fit_all$long_sub_random)
    r_hist <- rep(sigma2, vapply(history, nrow, integer(1)))
    r_fut <- rep(sigma2, each = n_times)

    ### posterior of b given the history is Gaussian, with a covariance that
    ### does not depend on (T, D): b | . ~ N(K (y - mu(T, D)), Sigma - K Z Sigma)
    if (length(y_hist) > 0) {
      V_y <- Z_hist %*% Sigma %*% t(Z_hist) + diag(r_hist, length(r_hist))
      V_y_chol <- chol(V_y)
      K <- Sigma %*% t(Z_hist) %*% chol2inv(V_y_chol)
      V_b <- Sigma - K %*% Z_hist %*% Sigma
    } else {
      K <- NULL
      V_b <- Sigma
    }
    V_b_sqrt <- symmetric_sqrt(V_b)

    ### log weight of each (T, D) cell and the means it implies
    cells <- expand.grid(k = seq_along(l_grid), d = seq_along(d_values))
    log_w <- numeric(nrow(cells))
    mu_hist_cells <- mu_fut_cells <- vector("list", length(d_values))
    for (j in seq_along(d_values)) {
      mu_hist_cells[[j]] <- stacked_mean_on_grid(history, long_fit_all, l_grid, d_values[j], survival_variable,
                                                 event_type_variable, survival_variable_all, survival_trans_function)
      mu_fut_cells[[j]] <- stacked_mean_on_grid(future, long_fit_all, l_grid, d_values[j], survival_variable,
                                                event_type_variable, survival_variable_all, survival_trans_function)
      log_lik <- if (length(y_hist) > 0) gaussian_logdens_chol(V_y_chol, y_hist - mu_hist_cells[[j]]) else 0
      log_w[cells$d == j] <- log_lik + log_prior[[j]][, p]
    }
    log_w[!is.finite(log_w)] <- -Inf
    if (all(log_w == -Inf)) {
      warning(sprintf("Patient %s has no event time with positive posterior weight; skipped.", patient_ids[p]),
              call. = FALSE)
      next
    }
    w <- exp(log_w - max(log_w))
    drawn <- cells[sample.int(nrow(cells), n_sim, replace = TRUE, prob = w), , drop = FALSE]

    Y <- matrix(NA_real_, M * n_times, n_sim)
    for (s in seq_len(n_sim)) {
      k <- drawn$k[s]
      j <- drawn$d[s]
      b_mean <- if (is.null(K)) rep(0, ncol(Sigma)) else K %*% (y_hist - mu_hist_cells[[j]][, k])
      b <- b_mean + V_b_sqrt %*% stats::rnorm(ncol(Sigma))
      Y[, s] <- mu_fut_cells[[j]][, k] + Z_fut %*% b + stats::rnorm(length(r_fut), sd = sqrt(r_fut))
    }

    event_time <- l_grid[drawn$k]
    res <- data.frame(rep(patient_ids[p], n_sim * n_times),
                      rep(seq_len(n_sim), each = n_times),
                      rep(times, n_sim),
                      rep(event_time, each = n_times))
    names(res) <- c(id, "sim", time_variable, "event_time")
    if (has_cr) res[[event_type_variable]] <- rep(d_values[drawn$d], each = n_times)
    after_event <- res[[time_variable]] >= res$event_time
    for (i in seq_len(M)) {
      values <- c(Y[(i - 1) * n_times + seq_len(n_times), , drop = FALSE])
      if (truncate) values[after_event] <- NA
      res[[response_names[i]]] <- values
    }
    out[[p]] <- res
  }

  out <- do.call(rbind, out)
  rownames(out) <- NULL
  attr(out, "prediction_time") <- prediction_time
  attr(out, "times") <- times
  class(out) <- c("simulateTrajectory.BJM", "data.frame")
  out
}

#' A patient's history rows usable for conditioning
#'
#' @description The rows of one biomarker's data with the biomarker and
#' every covariate of its sub-model observed (the outcome columns, which are
#' set at each candidate event time, excepted). May have no rows.
#' @keywords internal
complete_history_rows <- function(d, i, long_fit_all, survival_variable, other_outcomes) {
  vars <- model_vars(long_fit_all$long_sub_fixed[[i]], long_fit_all$long_sub_random[[i]])
  vars <- intersect(setdiff(vars, c(survival_variable, other_outcomes)), names(d))
  d[rowSums(is.na(d[vars])) == 0, , drop = FALSE]
}

#' Random-effects design stacked across biomarkers
#'
#' @description Like \code{random_effects_design()}, but a biomarker may
#' contribute no rows (its block is then empty and its random effects
#' simply do not enter).
#' @keywords internal
stacked_random_design <- function(rows, long_sub_random) {
  blocks <- lapply(seq_along(rows), function(i) {
    ffk <- nlme::splitFormula(long_sub_random[[i]], "|")[[1]]
    q <- ncol(model.matrix(ffk, model.frame(ffk, rows[[i]][rep(1, max(nrow(rows[[i]]), 1)), , drop = FALSE],
                                            na.action = stats::na.pass)))
    if (nrow(rows[[i]]) == 0) return(matrix(0, 0, q))
    model.matrix(ffk, model.frame(ffk, rows[[i]], na.action = stats::na.pass))
  })
  q <- vapply(blocks, ncol, integer(1))
  n <- vapply(blocks, nrow, integer(1))
  Z <- matrix(0, sum(n), sum(q))
  row_off <- c(0, cumsum(n))
  col_off <- c(0, cumsum(q))
  for (i in seq_along(blocks)) {
    if (n[i] > 0) Z[row_off[i] + seq_len(n[i]), col_off[i] + seq_len(q[i])] <- blocks[[i]]
  }
  Z
}

#' Fixed-effects means at every candidate event time
#'
#' @description For each biomarker's rows, \eqn{X(l)\beta} with the
#' survival-time columns (and the event type, with competing risks) set to
#' each element of \code{l_grid}, stacked across biomarkers.
#' @return A matrix with one row per stacked observation and one column per
#'   element of \code{l_grid}.
#' @keywords internal
stacked_mean_on_grid <- function(rows, long_fit_all, l_grid, d, survival_variable, event_type_variable,
                                 survival_variable_all, survival_trans_function) {
  lfit <- long_fit_all$lfit
  blocks <- lapply(seq_along(rows), function(i) {
    data_i <- rows[[i]]
    if (nrow(data_i) == 0) return(matrix(0, 0, length(l_grid)))
    if (!is.null(event_type_variable)) data_i[[event_type_variable]] <- d
    response <- as.character(formula(lfit[[i]])[[2]])
    data_i[[response]][is.na(data_i[[response]])] <- 0
    data_i <- set_survival_columns(data_i, survival_variable, l_grid[1], survival_variable_all,
                                   survival_trans_function)
    terms_i <- lfit[[i]]$terms
    xlev_i <- if (!is.null(long_fit_all$xlevels)) long_fit_all$xlevels[[i]] else NULL
    mf <- model.frame(terms_i, data_i, xlev = xlev_i, na.action = stats::na.pass)
    in_place <- time_columns_bare(terms_i, c(survival_variable, unlist(survival_variable_all)))
    beta <- nlme::fixef(lfit[[i]])
    vapply(l_grid, function(l) {
      mf <- survival_model_frame_at(mf, data_i, in_place, terms_i, xlev_i, survival_variable, l,
                                    survival_variable_all, survival_trans_function)
      c(model.matrix(terms_i, mf, contrasts.arg = lfit[[i]]$contrasts) %*% beta)
    }, numeric(nrow(data_i)))
  })
  blocks <- lapply(blocks, function(b) matrix(b, ncol = length(l_grid)))
  do.call(rbind, blocks)
}

#' Gaussian log density from a Cholesky factor
#'
#' @description Zero-mean multivariate normal log density of each column of
#' \code{E}, with covariance \eqn{R^T R} given by its upper Cholesky factor
#' \code{R}.
#' @keywords internal
gaussian_logdens_chol <- function(R, E) {
  E <- as.matrix(E)
  quad <- colSums(backsolve(R, E, transpose = TRUE)^2)
  -0.5 * (nrow(R) * log(2 * pi) + 2 * sum(log(diag(R))) + quad)
}

#' Square root of a covariance matrix
#'
#' @description \eqn{S} with \eqn{S S^T = V}, from the eigendecomposition so
#' a positive semi-definite \eqn{V} (e.g. a posterior covariance that has
#' lost definiteness to rounding) works too.
#' @keywords internal
symmetric_sqrt <- function(V) {
  eig <- eigen((V + t(V)) / 2, symmetric = TRUE)
  eig$vectors %*% diag(sqrt(pmax(eig$values, 0)), ncol(V))
}
