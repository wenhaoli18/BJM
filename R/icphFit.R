#' Proportional hazards model for interval-censored event times
#'
#' @description Internal fitter behind \code{survivalSub()} when the
#' outcome is interval-censored (\code{Surv(L, R, type = "interval2")}).
#' \code{\link[survival]{coxph}} cannot fit such data, so the marginal model
#' is a proportional hazards model \eqn{H(t \mid x) = H_0(t) e^{x\beta}}
#' with a parametric baseline, fit by maximum likelihood. Each subject
#' contributes \eqn{S(L) - S(R)} for \eqn{L < T \le R}, \eqn{S(L)} when
#' right-censored (\eqn{R = \infty}), \eqn{1 - S(R)} when left-censored, and
#' \eqn{h(T) S(T)} for an exactly observed \eqn{T}. Two baselines:
#' \describe{
#'   \item{\code{"spline"} (default)}{Royston--Parmar: \eqn{\log H_0(t)} is a
#'   natural cubic spline in \eqn{\log t} with \code{df} degrees of freedom
#'   (\code{df - 1} internal knots at quantiles of the log interval
#'   endpoints, boundary knots at their extremes). \code{df = 1} is a
#'   Weibull model.}
#'   \item{\code{"piecewise"}}{a piecewise-constant hazard on \code{df}
#'   pieces cut at quantiles of the interval endpoints.}
#' }
#' The spline is the default because interval-censored data say little
#' about the shape of the hazard \emph{within} a visit interval, and that
#' shape is exactly what \code{imputeEventTime()} draws \eqn{T} from: the
#' piecewise-constant hazard is flat there, which in simulations biased the
#' longitudinal sub-model fit by \code{fitIntervalBJM()} when the true
#' hazard was increasing, while the smooth spline carries the trend across
#' neighbouring intervals into each one.
#'
#' The baseline cumulative hazard is then tabulated on a fine grid in the
#' same \code{hazard}/\code{time} layout as \code{\link[survival]{basehaz}},
#' so the prediction code can use it exactly as it uses the Breslow
#' estimate of a right-censored fit (see \code{survival_accessors}).
#'
#' @param formula \code{Surv(L, R, type = "interval2") ~ covariates}. A
#'   \code{strata()} term gives each stratum its own baseline (spline knots
#'   or cut points from that stratum's interval endpoints), with the
#'   regression coefficients shared, as in a stratified Cox model.
#' @param data One row per subject.
#' @param event_time Name of the event-time column the longitudinal
#'   sub-model uses.
#' @param baseline \code{"spline"} or \code{"piecewise"}.
#' @param df Spline degrees of freedom, or number of pieces.
#' @param lower_tail Spline only: how \eqn{\log H_0} is extrapolated before
#'   the earliest interval endpoint (the lower boundary knot), where it never
#'   enters the likelihood. \code{"spline"} (default) continues the spline's
#'   own linear tail; \code{"weibull"} uses the slope of the nested Weibull
#'   fit (\code{df = 1}) instead. In simulations where a third to a half of
#'   the events fell before the first visit, \code{"weibull"} made
#'   \code{fitIntervalBJM()}'s estimates less variable but biased them when
#'   the true early hazard was not Weibull-shaped (Gompertz), so it is kept
#'   only for sensitivity analysis: the shape there is not identified by
#'   interval-censored data, and every choice is an assumption.
#' @return An object of class \code{"icph.BJM"}.
#' @keywords internal
icphFit <- function(formula, data, event_time, baseline = c("spline", "piecewise"), df = 3,
                    lower_tail = c("spline", "weibull")) {
  baseline <- match.arg(baseline)
  lower_tail <- match.arg(lower_tail)
  tt <- stats::terms(formula, specials = "strata", data = data)
  mf <- model.frame(tt, data, na.action = stats::na.omit)
  y <- model.response(mf)
  strata_info <- survival::untangle.specials(tt, "strata")
  strata_var <- if (length(strata_info$vars) > 0) strata_info$vars else NULL
  # covariates without the strata() term
  labels <- attr(stats::delete.response(tt), "term.labels")
  x_labels <- setdiff(labels, strata_var)
  terms_rhs <- stats::delete.response(stats::terms(
    stats::reformulate(if (length(x_labels) > 0) x_labels else "1", env = environment(formula))))
  X <- model.matrix(terms_rhs, mf)
  contrasts <- attr(X, "contrasts")
  X <- X[, colnames(X) != "(Intercept)", drop = FALSE]

  bounds <- interval_bounds(y)
  L <- bounds$L
  R <- bounds$R
  stratum <- if (is.null(strata_var)) rep("all", length(L)) else as.character(mf[[strata_var]])
  strata_levels <- unique(stratum)
  if (!is.null(strata_var)) strata_levels <- sort(strata_levels)
  rows <- lapply(strata_levels, function(s) which(stratum == s))

  # one baseline per stratum, on that stratum's own interval endpoints
  endpoints_of <- function(r) c(L[r][L[r] > 0], R[r][is.finite(R[r])])
  make_base <- function(endpoints, tail_slope = NULL) {
    if (baseline == "spline") rp_baseline(endpoints, df, tail_slope) else piecewise_baseline(endpoints, df)
  }
  bases <- lapply(seq_along(strata_levels), function(h) {
    endpoints <- endpoints_of(rows[[h]])
    if (length(unique(endpoints)) < df + 1) {
      stop(sprintf("Too few distinct interval endpoints (%d) for df = %d%s.",
                   length(unique(endpoints)), df,
                   if (is.null(strata_var)) "" else sprintf(" in stratum %s", strata_levels[h])),
           call. = FALSE)
    }
    make_base(endpoints)
  })
  k <- bases[[1]]$n_par
  H <- length(bases)
  gamma_index <- function(h) (h - 1) * k + seq_len(k)

  exact <- is.finite(R) & L == R
  fin <- is.finite(R) & !exact
  p <- ncol(X)
  negloglik <- function(theta) {
    e <- exp(if (p > 0) c(X %*% theta[H * k + seq_len(p)]) else rep(0, nrow(X)))
    total <- 0
    for (h in seq_len(H)) {
      g <- theta[gamma_index(h)]
      base <- bases[[h]]
      if (!base$monotone(g)) return(1e100)
      r <- rows[[h]]
      HL <- base$cumhaz(L[r], g) * e[r]
      ll <- -HL
      f <- fin[r]
      if (any(f)) {
        dH <- base$cumhaz(R[r][f], g) * e[r][f] - HL[f]
        ll[f] <- ll[f] + log(-expm1(-pmax(dH, 1e-300)))
      }
      x <- exact[r]
      if (any(x)) ll[x] <- ll[x] + base$loghaz(R[r][x], g) + log(e[r][x])
      total <- total - sum(ll)
    }
    if (is.finite(total)) total else 1e100
  }
  # constant-hazard start: events per unit of (midpoint) follow-up
  start <- c(unlist(lapply(seq_len(H), function(h) {
    r <- rows[[h]]
    rate0 <- max(sum(is.finite(R[r])), 1) /
      sum(pmax(ifelse(is.finite(R[r]), (L[r] + R[r]) / 2, L[r]), 1e-8))
    bases[[h]]$start(rate0)
  })), rep(0, p))
  if (baseline == "spline" && df > 1) {
    # start from the Weibull fit (df = 1, nested in this one: same intercept
    # and scaled log-t column, zero weight on the extra spline terms). From
    # the constant-hazard start, BFGS sometimes wandered into a non-monotone
    # region and stopped at a nonsensical baseline.
    weibull <- icphFit(formula, data, event_time, "spline", 1)
    start <- c(unlist(lapply(weibull$gammas, function(g) c(g, rep(0, df - 1)))), weibull$coefficients)
    if (lower_tail == "weibull") {
      bases <- lapply(seq_len(H), function(h) {
        make_base(endpoints_of(rows[[h]]),
                  tail_slope = weibull$gammas[[h]][2] / weibull$bases[[h]]$scale[1])
      })
    }
  }
  opt <- stats::optim(start, negloglik, method = "BFGS", hessian = TRUE,
                      control = list(maxit = 2000, reltol = 1e-12))
  if (opt$convergence != 0) {
    warning("The interval-censored survival sub-model did not converge.", call. = FALSE)
  }

  coef_names <- colnames(X)
  vc <- tryCatch(solve(opt$hessian), error = function(e) {
    matrix(NA_real_, length(start), length(start))
  })
  beta <- stats::setNames(opt$par[H * k + seq_len(p)], coef_names)
  vcov_beta <- vc[H * k + seq_len(p), H * k + seq_len(p), drop = FALSE]
  dimnames(vcov_beta) <- list(coef_names, coef_names)
  gammas <- stats::setNames(lapply(seq_len(H), function(h) opt$par[gamma_index(h)]), strata_levels)

  max_followup <- max(endpoints_of(seq_along(L)))
  cum_basehaz <- do.call(rbind, lapply(seq_len(H), function(h) {
    upper <- max(endpoints_of(rows[[h]]))
    knots <- bases[[h]]$knot_times
    grid <- sort(unique(c(seq(0, upper, length.out = 2001), knots[knots <= upper])))
    H0 <- bases[[h]]$cumhaz(grid, gammas[[h]])
    if (is.unsorted(H0)) {
      warning("The fitted baseline cumulative hazard is not monotone; it was made monotone with cummax().",
              call. = FALSE)
      H0 <- cummax(H0)
    }
    out <- data.frame(hazard = H0, time = grid)
    if (!is.null(strata_var)) out$strata <- strata_levels[h]
    out
  }))

  out <- list(coefficients = beta, var = vcov_beta, baseline = baseline, df = df,
              gamma = gammas[[1]], base = bases[[1]], gammas = gammas, bases = bases,
              strata_var = strata_var, strata_levels = if (!is.null(strata_var)) strata_levels,
              loglik = -opt$value, n = nrow(X), n_exact = sum(exact),
              n_interval = sum(fin), n_right = sum(!is.finite(R)),
              n_before_first_visit = sum(fin & L == 0),
              terms = terms_rhs, xlevels = stats::.getXlevels(terms_rhs, mf),
              contrasts = contrasts, cum_basehaz = cum_basehaz,
              max_followup = max_followup, event_time = event_time,
              L = L, R = R, formula = formula)
  if (baseline == "piecewise") {
    out$lambda <- exp(gammas[[1]])
    out$cuts <- bases[[1]]$cuts
  }
  class(out) <- "icph.BJM"
  out
}

#' Stratum of each patient for a stratified interval-censored fit
#' @param fit An \code{icph.BJM} object.
#' @param newdata One row per patient.
#' @return Character labels in the fit's \code{strata()} format, or
#'   \code{NULL} when the fit is not stratified.
#' @keywords internal
icph_strata <- function(fit, newdata) {
  if (is.null(fit$strata_var)) return(NULL)
  env <- list2env(list(strata = survival::strata), parent = environment(fit$formula))
  as.character(eval(str2lang(fit$strata_var), newdata, env))
}

#' Baselines for \code{icphFit()}
#'
#' @description Each returns a list with \code{n_par}, \code{cumhaz(t, g)}
#' (baseline cumulative hazard, 0 at \code{t <= 0}, \code{Inf} at
#' \code{Inf}), \code{loghaz(t, g)} (log baseline hazard, for exactly
#' observed times), \code{monotone(g)} (whether \code{cumhaz} is
#' increasing), \code{start(rate0)} (parameters of a constant hazard
#' \code{rate0}) and \code{knot_times}.
#' @param endpoints Finite, positive interval endpoints.
#' @param df Spline degrees of freedom / number of pieces.
#' @name icph_baselines
#' @keywords internal
NULL

#' @describeIn icph_baselines Royston--Parmar spline for \eqn{\log H_0} in
#'   \eqn{\log t}. Parameters: intercept, then one per basis column; the
#'   basis columns are scaled by their standard deviation over
#'   \code{endpoints} so the optimizer sees comparable scales.
rp_baseline <- function(endpoints, df, tail_slope = NULL) {
  x <- log(endpoints)
  knots <- c(min(x), if (df > 1) stats::quantile(x, seq_len(df - 1) / df, names = FALSE), max(x))
  scale <- apply(rp_basis(x, knots), 2, stats::sd)
  scale[!is.finite(scale) | scale == 0] <- 1
  s_fun <- function(t, g) {
    out <- rep(-Inf, length(t))
    pos <- t > 0 & is.finite(t)
    B <- sweep(rp_basis(log(t[pos]), knots), 2, scale, "/")
    out[pos] <- g[1] + c(B %*% g[-1])
    if (!is.null(tail_slope)) {
      below <- pos & log(pmax(t, 1e-300)) < knots[1]
      if (any(below)) {
        s_kmin <- g[1] + c(sweep(rp_basis(knots[1], knots), 2, scale, "/") %*% g[-1])
        out[below] <- s_kmin + tail_slope * (log(t[below]) - knots[1])
      }
    }
    out[is.infinite(t) & t > 0] <- Inf
    out
  }
  list(
    n_par = df + 1,
    knots = knots, scale = scale,
    cumhaz = function(t, g) exp(s_fun(t, g)),
    loghaz = function(t, g) {
      dB <- sweep(rp_basis(log(t), knots, deriv = TRUE), 2, scale, "/")
      slope <- c(dB %*% g[-1])
      if (!is.null(tail_slope)) slope[log(t) < knots[1]] <- tail_slope
      ifelse(slope > 0, s_fun(t, g) + log(pmax(slope, 1e-300)) - log(t), -Inf)
    },
    # the spline is linear in log t outside the boundary knots, so H0 is
    # increasing everywhere iff its slope is positive between them. The
    # likelihood alone does not enforce this before the first visit, where
    # there are no interval endpoints: unconstrained fits occasionally had
    # H0 rising towards t = 0.
    monotone = function(g) {
      xs <- seq(knots[1], knots[length(knots)], length.out = 101)
      all(sweep(rp_basis(xs, knots, deriv = TRUE), 2, scale, "/") %*% g[-1] > 0)
    },
    # H0(t) = rate0 * t: log H0 = log(rate0) + log t, i.e. a unit slope on log t
    start = function(rate0) c(log(rate0), scale[1], rep(0, df - 1)),
    knot_times = exp(knots)
  )
}

#' @describeIn icph_baselines Piecewise-constant hazard; parameters are the
#'   log hazards of the pieces.
piecewise_baseline <- function(endpoints, df) {
  cuts <- c(0, unique(stats::quantile(endpoints, probs = seq_len(df - 1) / df, names = FALSE)), Inf)
  n <- length(cuts) - 1
  list(
    n_par = n, cuts = cuts,
    cumhaz = function(t, g) piecewise_cumhaz(t, exp(g), cuts),
    loghaz = function(t, g) {
      piece <- pmax(findInterval(t, cuts, left.open = TRUE, rightmost.closed = TRUE), 1)
      g[piece]
    },
    monotone = function(g) TRUE,
    start = function(rate0) rep(log(rate0), n),
    knot_times = cuts[is.finite(cuts) & cuts > 0]
  )
}

#' Restricted cubic spline basis (Royston--Parmar)
#'
#' @description Columns \eqn{x} and, for each internal knot \eqn{k_j},
#' \eqn{(x - k_j)_+^3 - \lambda_j (x - k_{min})_+^3 - (1 - \lambda_j)
#' (x - k_{max})_+^3}, \eqn{\lambda_j = (k_{max} - k_j) / (k_{max} -
#' k_{min})}: a cubic spline that is linear beyond the boundary knots. With
#' \code{deriv = TRUE}, the derivative of each column in \eqn{x}.
#' @param x Values (log times).
#' @param knots Boundary and internal knots, sorted, boundary first and last.
#' @param deriv Return the derivative instead.
#' @keywords internal
rp_basis <- function(x, knots, deriv = FALSE) {
  kmin <- knots[1]
  kmax <- knots[length(knots)]
  internal <- knots[-c(1, length(knots))]
  pos3 <- function(u) if (deriv) 3 * pmax(u, 0)^2 else pmax(u, 0)^3
  cols <- lapply(internal, function(k) {
    lam <- (kmax - k) / (kmax - kmin)
    pos3(x - k) - lam * pos3(x - kmin) - (1 - lam) * pos3(x - kmax)
  })
  cbind(if (deriv) rep(1, length(x)) else x, do.call(cbind, cols))
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
