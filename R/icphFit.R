#' Proportional hazards model for interval-censored event times
#'
#' @description Internal fitter behind \code{survivalSub()} when the
#' outcome is interval-censored (\code{Surv(L, R, type = "interval2")}).
#' \code{\link[survival]{coxph}} cannot fit such data, so the marginal model
#' is a proportional hazards model with a piecewise-constant baseline
#' hazard, \eqn{h(t \mid x) = \lambda_k e^{x\beta}} on the \eqn{k}-th of
#' \code{n_pieces} intervals, fit by maximum likelihood. The cut points are
#' quantiles of the finite, positive interval endpoints; the last piece
#' extends to infinity. Each subject contributes
#' \eqn{S(L) - S(R)} for \eqn{L < T \le R}, \eqn{S(L)} when right-censored
#' (\eqn{R = \infty}), \eqn{1 - S(R)} when left-censored, and
#' \eqn{h(T) S(T)} for an exactly observed \eqn{T}.
#'
#' The baseline cumulative hazard is then tabulated on a fine grid in the
#' same \code{hazard}/\code{time} layout as \code{\link[survival]{basehaz}},
#' so the prediction code can use it exactly as it uses the Breslow
#' estimate of a right-censored fit (see \code{survival_accessors}).
#'
#' @param formula \code{Surv(L, R, type = "interval2") ~ covariates}.
#'   \code{strata()} terms are not supported.
#' @param data One row per subject.
#' @param event_time Name of the event-time column the longitudinal
#'   sub-model uses.
#' @param n_pieces Number of pieces of the baseline hazard.
#' @return An object of class \code{"icph.BJM"}.
#' @keywords internal
icphFit <- function(formula, data, event_time, n_pieces = 6) {
  if (length(survival::untangle.specials(stats::terms(formula, specials = "strata"), "strata")$vars) > 0) {
    stop("strata() terms are not supported for an interval-censored survival sub-model.", call. = FALSE)
  }
  mf <- model.frame(formula, data, na.action = stats::na.omit)
  y <- model.response(mf)
  terms_rhs <- stats::delete.response(stats::terms(mf))
  X <- model.matrix(terms_rhs, mf)
  contrasts <- attr(X, "contrasts")
  X <- X[, colnames(X) != "(Intercept)", drop = FALSE]

  bounds <- interval_bounds(y)
  L <- bounds$L
  R <- bounds$R
  endpoints <- c(L[L > 0], R[is.finite(R)])
  if (length(unique(endpoints)) < n_pieces) {
    stop(sprintf("Too few distinct interval endpoints (%d) for n_pieces = %d.",
                 length(unique(endpoints)), n_pieces), call. = FALSE)
  }
  cuts <- c(0, unique(stats::quantile(endpoints, probs = seq_len(n_pieces - 1) / n_pieces,
                                      names = FALSE)), Inf)
  n_pieces <- length(cuts) - 1

  exact <- is.finite(R) & L == R
  p <- ncol(X)
  negloglik <- function(theta) {
    lambda <- exp(theta[seq_len(n_pieces)])
    eta <- if (p > 0) c(X %*% theta[n_pieces + seq_len(p)]) else rep(0, nrow(X))
    -sum(icph_loglik_terms(L, R, exact, eta, lambda, cuts))
  }
  # crude constant-hazard start: events per unit of (left-endpoint) follow-up
  n_events <- sum(is.finite(R))
  rate0 <- max(n_events, 1) / sum(pmax(ifelse(is.finite(R), (L + R) / 2, L), 1e-8))
  start <- c(rep(log(rate0), n_pieces), rep(0, p))
  opt <- stats::optim(start, negloglik, method = "BFGS", hessian = TRUE,
                      control = list(maxit = 1000, reltol = 1e-12))
  if (opt$convergence != 0) {
    warning("The interval-censored survival sub-model did not converge.", call. = FALSE)
  }

  coef_names <- colnames(X)
  vc <- tryCatch(solve(opt$hessian), error = function(e) {
    matrix(NA_real_, length(start), length(start))
  })
  beta <- stats::setNames(opt$par[n_pieces + seq_len(p)], coef_names)
  vcov_beta <- vc[n_pieces + seq_len(p), n_pieces + seq_len(p), drop = FALSE]
  dimnames(vcov_beta) <- list(coef_names, coef_names)
  lambda <- exp(opt$par[seq_len(n_pieces)])

  max_followup <- max(endpoints)
  grid <- sort(unique(c(seq(0, max_followup, length.out = 2001), cuts[is.finite(cuts) & cuts <= max_followup])))
  cum_basehaz <- data.frame(hazard = piecewise_cumhaz(grid, lambda, cuts), time = grid)

  out <- list(coefficients = beta, var = vcov_beta, lambda = lambda, cuts = cuts,
              loglik = -opt$value, n = nrow(X), n_exact = sum(exact),
              n_interval = sum(is.finite(R) & !exact),
              n_right = sum(!is.finite(R)),
              terms = terms_rhs, xlevels = stats::.getXlevels(stats::terms(mf), mf),
              contrasts = contrasts, cum_basehaz = cum_basehaz,
              max_followup = max_followup, event_time = event_time,
              L = L, R = R, formula = formula)
  class(out) <- "icph.BJM"
  out
}

#' Interval endpoints from an interval-type \code{Surv} object
#'
#' @description \eqn{L < T \le R} for every subject: \eqn{R = \infty} when
#' right-censored, \eqn{L = 0} when left-censored, \eqn{L = R} when \eqn{T}
#' is observed exactly.
#' @param y A \code{Surv} object of type \code{"interval"}.
#' @return A list with numeric vectors \code{L} and \code{R}.
#' @keywords internal
interval_bounds <- function(y) {
  status <- y[, "status"]
  t1 <- y[, "time1"]
  t2 <- y[, "time2"]
  L <- ifelse(status == 2, 0, t1)
  R <- ifelse(status == 0, Inf, ifelse(status == 1 | status == 2, t1, t2))
  list(L = unname(L), R = unname(R))
}

#' Cumulative hazard of a piecewise-constant hazard
#' @param t Times.
#' @param lambda Hazard on each piece.
#' @param cuts Piece boundaries, starting at 0 and ending at \code{Inf}.
#' @keywords internal
piecewise_cumhaz <- function(t, lambda, cuts) {
  widths <- diff(cuts)
  exposure <- pmin(pmax(outer(t, cuts[-length(cuts)], "-"), 0), matrix(widths, length(t), length(widths), byrow = TRUE))
  exposure[!is.finite(t), ] <- Inf
  c(exposure %*% lambda)
}

#' Per-subject log-likelihood of the interval-censored PH model
#' @keywords internal
icph_loglik_terms <- function(L, R, exact, eta, lambda, cuts) {
  e <- exp(eta)
  HL <- piecewise_cumhaz(L, lambda, cuts) * e
  out <- -HL
  fin <- is.finite(R) & !exact
  if (any(fin)) {
    dH <- (piecewise_cumhaz(R[fin], lambda, cuts) * e[fin]) - HL[fin]
    out[fin] <- out[fin] + log(-expm1(-pmax(dH, 1e-300)))
  }
  if (any(exact)) {
    piece <- findInterval(R[exact], cuts, left.open = TRUE, rightmost.closed = TRUE)
    piece <- pmax(piece, 1)
    out[exact] <- out[exact] + log(lambda[piece]) + eta[exact]
  }
  out
}

#' Linear predictor of an interval-censored PH fit
#' @param fit An \code{icph.BJM} object.
#' @param newdata One row per patient.
#' @keywords internal
icph_lp <- function(fit, newdata) {
  mf <- model.frame(fit$terms, newdata, na.action = stats::na.pass, xlev = fit$xlevels)
  X <- model.matrix(fit$terms, mf, contrasts.arg = fit$contrasts)
  X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  if (length(fit$coefficients) == 0) return(rep(0, nrow(X)))
  c(X[, names(fit$coefficients), drop = FALSE] %*% fit$coefficients)
}
