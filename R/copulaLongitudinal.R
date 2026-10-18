#' Resolve the per-biomarker type used by \code{longitudinalSub()}
#'
#' @description If \code{biomarker_type} is supplied explicitly, it always
#' wins (after validation). Otherwise each biomarker's type is auto-detected
#' from its own response column in \code{data_fit_all}: a \code{factor} (or
#' \code{ordered factor}) response is treated as \code{"ordinal"}; anything
#' else is treated as \code{"continuous"}. This mirrors exactly the priority
#' rule requested for \code{\link{longitudinalSub}}: explicit
#' \code{biomarker_type} overrides auto-detection, and auto-detection is only
#' consulted when \code{biomarker_type} is \code{NULL}.
#' @keywords internal
resolve_biomarker_type <- function(biomarker_type, data_fit_all, long_sub_fixed, M) {
  if (!is.null(biomarker_type)) {
    assert_biomarker_type(biomarker_type, M)
    return(biomarker_type)
  }
  vapply(seq_len(M), function(m) {
    resp_name <- all.vars(long_sub_fixed[[m]])[1]
    resp <- data_fit_all[[m]][[resp_name]]
    if (is.factor(resp)) "ordinal" else "continuous"
  }, character(1))
}

#' Convert an nlme-style random-effects formula to an lme4-style bar term
#'
#' @description \code{longitudinalSub()}'s \code{long_sub_random} argument
#' uses \code{nlme::lme()}'s \code{~ terms | group} formula convention.
#' \code{ordinal::clmm()} (used to fit the marginal model for a categorical/
#' ordinal biomarker -- see \code{fit_marginal_ordinal()}) instead expects
#' random effects written as an \code{lme4}-style \code{(terms | group)} bar
#' term embedded directly in the model formula. This helper translates one
#' into the other so that the exact same \code{long_sub_random} argument
#' drives both \code{nlme::lme()} (continuous biomarkers) and
#' \code{ordinal::clmm()} (ordinal biomarkers), with no second random-effects
#' argument for the user to keep in sync.
#' @return A length-1 character string, e.g. \code{"(year | id)"}.
#' @keywords internal
nlme_random_to_lme4_bars <- function(random_formula) {
  parts <- nlme::splitFormula(random_formula, "|")
  terms_str <- deparse(parts[[1]][[2]])
  group_str <- as.character(parts[[2]])[2]
  sprintf("(%s | %s)", terms_str, group_str)
}

#' Fit the marginal model for a single ordinal/categorical biomarker
#'
#' @description Fits a cumulative link mixed model (probit link) via
#' \code{ordinal::clmm()}, used as the marginal model for a categorical
#' biomarker in the Gaussian-copula extension of \code{longitudinalSub()}
#' (see \code{longitudinalSubCopula()}). The probit link is used
#' specifically (rather than the more common logit) because it integrates
#' directly with the latent-Gaussian-score representation the copula
#' extension relies on: the underlying continuous latent variable implied by
#' a probit cumulative link model is, by construction, Gaussian, which is
#' exactly the representation needed to fold a categorical biomarker into
#' the same random-effects covariance structure as the continuous
#' biomarkers.
#'
#' This follows the classical two-stage IFM (Inference Functions for
#' Margins; Joe 2005) strategy already used throughout \code{BJM}: fit each
#' biomarker's marginal model separately first (here, exactly as
#' \code{nlme::lme()} does for continuous biomarkers), then re-estimate the
#' joint random-effects covariance across all biomarkers afterwards (see
#' \code{longitudinalSubVarCopula()}).
#' @return A fitted \code{ordinal::clmm} object.
#' @keywords internal
fit_marginal_ordinal <- function(fixed_formula, random_formula, data) {
  assert_package_installed("ordinal", "Fitting an ordinal/categorical biomarker with `longitudinalSub()`")

  resp_name <- all.vars(fixed_formula)[1]
  if (is.factor(data[[resp_name]]) && !is.ordered(data[[resp_name]])) {
    # A binary (2-level) factor has only one cumulative-link threshold, so
    # swapping which level is "low" vs "high" only flips the sign of every
    # coefficient/threshold -- log-likelihood and model fit are unaffected.
    # The arbitrary-order problem this warning exists for is specific to 3+
    # level factors, where the *relative* order of middle categories cannot
    # be fixed by any sign flip; skip the warning (but still coerce to
    # ordered, since ordinal::clmm() requires it) for the 2-level case.
    if (nlevels(data[[resp_name]]) > 2) {
      warning(sprintf(paste0(
        "Biomarker response '%s' is an unordered factor; treating its factor level ",
        "order (%s) as the ordinal category order for the Gaussian-copula fit. Use ",
        "an ordered factor, or pass `biomarker_type` explicitly, to make this order ",
        "intentional."
      ), resp_name, paste(levels(data[[resp_name]]), collapse = " < ")), call. = FALSE)
    }
    data[[resp_name]] <- factor(data[[resp_name]], levels = levels(data[[resp_name]]), ordered = TRUE)
  }

  bar_term <- nlme_random_to_lme4_bars(random_formula)
  full_formula <- stats::as.formula(
    paste(paste(deparse(fixed_formula), collapse = " "), "+", bar_term)
  )

  ordinal::clmm(full_formula, data = data, link = "probit", Hess = TRUE)
}

#' Joint random-effects covariance with ordinal biomarkers (Gaussian copula)
#'
#' @description Re-estimates the shared random-effects covariance \code{D}
#' across all biomarkers by maximum likelihood, holding each biomarker's fixed
#' effects, residual variance and (for an ordinal biomarker) thresholds at
#' their separate-fit values, as \code{longitudinalSubVar()} does for
#' continuous ones. An ordinal observation is the category of a latent score
#' \eqn{Z = \eta + \epsilon}, \eqn{\epsilon \sim N(0, 1)} (probit link,
#' Albert & Chib), observed only through the thresholds that bracket it, so
#' the posterior of the random effects given the data is not Gaussian and the
#' EM's \eqn{E[b b^T \mid y]} has no closed form. It is computed by
#' \itemize{
#'   \item \code{copula_D_aghq()}, an EM with adaptive Gauss-Hermite
#'   quadrature over the ordinal biomarkers' random effects (the continuous
#'   biomarkers' random effects are integrated out exactly), when those
#'   number at most \code{max_quadrature_dim}; deterministic, with the usual
#'   convergence check;
#'   \item otherwise \code{copula_D_saem()}, a stochastic-approximation EM
#'   (SAEM) that draws the latent scores, whose cost does not grow with that
#'   dimension.
#' }
#' Both target the same maximum-likelihood estimate, but with the default
#' number of iterations SAEM stopped short of it: for two ordinal biomarkers
#' with a random intercept and slope each (four ordinal random effects), the
#' random-effects variances were off by about 5\% in simulations, which
#' quadrature (\eqn{5^4 = 625} nodes per subject) removes; hence the default
#' \code{max_quadrature_dim = 4}. The earlier version
#' plugged in the mean of each latent score given its category, which ignored
#' the score's own variance and underestimated the ordinal biomarkers'
#' random-effects variances (by about 40\% for a random intercept and 75\%
#' for a random slope in simulations).
#' @param thetaLong List with the stacked fixed effects \code{beta}, the
#'   starting \code{D} and the residual variances \code{sigma2} (1 for an
#'   ordinal biomarker).
#' @param l Per-subject design bookkeeping built by
#'   \code{longitudinalSubCopula()}; \code{l$yik[[m]]} holds an ordinal
#'   biomarker's category codes.
#' @param biomarker_type,thresholds As stored in the fit.
#' @param tol.em,max.iter Convergence tolerance on \code{D} and iteration cap
#'   of the quadrature EM.
#' @param max_quadrature_dim Largest number of ordinal random effects handled
#'   by quadrature.
#' @param verbose Print the change in \code{D} at every iteration.
#' @return A list with the estimate \code{D}.
#' @keywords internal
longitudinalSubVarCopula <- function(thetaLong, l, biomarker_type, thresholds,
                                     tol.em = 1e-04, max.iter = 500, max_quadrature_dim = 4,
                                     verbose = FALSE) {
  setup <- copula_estep_setup(thetaLong, l, biomarker_type, thresholds)
  ### an EM started from a singular D stays on that boundary; the separate
  ### fits can give one (e.g. clmm() estimating a perfectly correlated random
  ### intercept and slope from few visits), so start from its diagonal then
  D_start <- thetaLong$D
  ev <- eigen(D_start, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) < 1e-8 * max(ev)) D_start[] <- diag(diag(D_start), nrow(D_start))
  if (length(setup$ro) <= max_quadrature_dim) {
    copula_D_aghq(D_start, setup, tol.em = tol.em, max.iter = max.iter, verbose = verbose)
  } else {
    copula_D_saem(D_start, setup, verbose = verbose)
  }
}

#' Per-subject pieces of the copula E-step that do not depend on D
#'
#' @description For every subject: the ordinal rows' fixed-effect linear
#' predictor, the thresholds bracketing each observed category and the
#' random-effects design of those rows; and the continuous rows' residuals
#' from the fixed effects, their random-effects design and residual
#' precisions. Also the positions \code{ro}/\code{rc} of the ordinal and
#' continuous biomarkers' random effects in \code{D}.
#' @inheritParams longitudinalSubVarCopula
#' @return A list with \code{subjects} (one entry per subject), \code{ro}
#'   and \code{rc}.
#' @keywords internal
copula_estep_setup <- function(thetaLong, l, biomarker_type, thresholds) {
  ord <- which(biomarker_type == "ordinal")
  r_off <- c(0, cumsum(l$r))
  ro <- unlist(lapply(ord, function(m) r_off[m] + seq_len(l$r[m])))
  rc <- setdiff(seq_len(sum(l$r)), ro)
  ids <- names(l$yi)
  subjects <- lapply(ids, function(i) {
    n_im <- l$nik[[i]]
    off <- c(0, cumsum(n_im))
    rows <- lapply(seq_along(n_im), function(m) off[m] + seq_len(n_im[m]))
    o <- unlist(rows[ord])
    cc <- setdiff(seq_len(sum(n_im)), o)
    lo <- hi <- numeric(0)
    for (m in ord) {
      alpha <- c(-Inf, thresholds[[m]], Inf)
      y <- l$yik[[m]][[i]]
      lo <- c(lo, alpha[y])
      hi <- c(hi, alpha[y + 1])
    }
    mu <- as.vector(l$Xi[[i]] %*% thetaLong$beta)
    Z <- l$Zi[[i]]
    list(mu_o = mu[o], lo = lo, hi = hi, Zo = Z[o, ro, drop = FALSE],
         Zc = Z[cc, rc, drop = FALSE], resid_c = l$yi[[i]][cc] - mu[cc],
         Rinv_c = 1 / rep(thetaLong$sigma2, n_im)[cc])
  })
  names(subjects) <- ids
  list(subjects = subjects, ro = ro, rc = rc)
}

#' Log probability of an ordinal category and its derivatives
#'
#' @description \eqn{\log P(lo < \eta + \epsilon \le hi)} for
#' \eqn{\epsilon \sim N(0, 1)}, with its first and second derivatives in
#' \eqn{\eta}; evaluated on the side of the mean where the interval lies, so
#' a category far in a tail keeps its precision.
#' @param eta Linear predictors; \code{lo}, \code{hi} the category bounds.
#' @param lo,hi Lower and upper thresholds (\code{-Inf}/\code{Inf} at the ends).
#' @return A list with \code{logP}, \code{g} (first derivative) and \code{h}
#'   (second derivative, \eqn{\le 0}).
#' @keywords internal
probit_box_derivs <- function(eta, lo, hi) {
  u <- lo - eta
  v <- hi - eta
  upper_tail <- (u + v) / 2 > 0
  P <- ifelse(upper_tail, stats::pnorm(-u) - stats::pnorm(-v), stats::pnorm(v) - stats::pnorm(u))
  P <- pmax(P, 1e-300)
  fu <- stats::dnorm(u)
  fv <- stats::dnorm(v)
  ufu <- ifelse(is.finite(u), u * fu, 0)
  vfv <- ifelse(is.finite(v), v * fv, 0)
  g <- (fu - fv) / P
  list(logP = log(P), g = g, h = (ufu - vfv) / P - g^2)
}

#' Gauss-Hermite quadrature rule for the standard normal
#'
#' @description Nodes and weights (summing to 1) of the \code{k}-point rule
#' for \eqn{E f(X)}, \eqn{X \sim N(0, 1)}, by the Golub-Welsch eigenvalue
#' method.
#' @param k Number of nodes.
#' @return A list with \code{x} (nodes) and \code{w} (weights).
#' @keywords internal
gauss_hermite_normal <- function(k) {
  J <- matrix(0, k, k)
  if (k > 1) {
    off <- sqrt(seq_len(k - 1))
    J[cbind(seq_len(k - 1), 2:k)] <- off
    J[cbind(2:k, seq_len(k - 1))] <- off
  }
  e <- eigen(J, symmetric = TRUE)
  list(x = e$values, w = e$vectors[1, ]^2)
}

#' Copula EM with adaptive Gauss-Hermite quadrature
#'
#' @description Given the ordinal random effects \eqn{b_o}, the continuous
#' ones \eqn{b_c} have a Gaussian posterior with mean \eqn{a + B b_o} and
#' covariance \eqn{A_c}, and the continuous data contribute a Gaussian factor
#' in \eqn{b_o}; only \eqn{b_o} is integrated numerically. For each subject
#' the posterior mode of \eqn{b_o} (Newton, warm-started from the previous
#' iteration) and the curvature there place a product Gauss-Hermite grid,
#' which gives \eqn{E[b b^T \mid y]} block by block; their average over
#' subjects is the next \eqn{D}. Uses 7 nodes per dimension up to two ordinal
#' random effects and 5 beyond. The EM is accelerated by SQUAREM (Varadhan and
#' Roland, 2008), falling back to plain EM steps where an extrapolated
#' \eqn{D} is not positive definite; convergence is a change in \eqn{D},
#' scaled by \eqn{\sqrt{D_{jj} D_{kk}}} so that it does not depend on the
#' biomarkers' units, below \code{tol.em} over one accelerated cycle (three
#' EM updates).
#' @param D Starting covariance.
#' @param setup Output of \code{copula_estep_setup()}.
#' @inheritParams longitudinalSubVarCopula
#' @return A list with \code{D}.
#' @keywords internal
copula_D_aghq <- function(D, setup, tol.em = 1e-04, max.iter = 500, verbose = FALSE) {
  ro <- setup$ro
  rc <- setup$rc
  q <- length(ro)
  gh <- gauss_hermite_normal(if (q <= 2) 7 else 5)
  grid <- as.matrix(expand.grid(rep(list(gh$x), q)))
  log_w_grid <- rowSums(log(as.matrix(expand.grid(rep(list(gh$w), q))))) + 0.5 * rowSums(grid^2)
  n_nodes <- nrow(grid)
  mode <- lapply(setup$subjects, function(s) rep(0, q))

  ### one EM update of D (also moves the posterior modes it warm-starts from)
  em_update <- function(D) {
    Doo_inv <- pd_inverse(D[ro, ro, drop = FALSE])
    if (length(rc) > 0) {
      C <- D[rc, ro, drop = FALSE] %*% Doo_inv           # E[b_c | b_o] = C b_o
      Vc_inv <- pd_inverse(D[rc, rc, drop = FALSE] - C %*% D[ro, rc, drop = FALSE])
    }
    S <- matrix(0, nrow(D), ncol(D))
    for (i in names(setup$subjects)) {
      s <- setup$subjects[[i]]
      # the Gaussian part of log p(b_o, y_c): -0.5 b_o' P b_o + h' b_o
      if (length(rc) > 0) {
        ZRZ <- crossprod(s$Zc, s$Rinv_c * s$Zc)
        Ac <- solve(ZRZ + Vc_inv)
        a <- as.vector(Ac %*% crossprod(s$Zc, s$Rinv_c * s$resid_c))
        B <- C - Ac %*% ZRZ %*% C                        # E[b_c | y_c, b_o] = a + B b_o
        RZ <- s$Rinv_c * s$Zc
        W_inv <- diag(s$Rinv_c, length(s$Rinv_c)) - RZ %*% Ac %*% t(RZ)   # Var(y_c | b_o)^-1
        ZC <- s$Zc %*% C
        P <- Doo_inv + t(ZC) %*% W_inv %*% ZC
        h <- as.vector(t(ZC) %*% W_inv %*% s$resid_c)
      } else {
        P <- Doo_inv
        h <- rep(0, q)
      }
      # posterior mode of b_o and the curvature there
      bo <- mode[[i]]
      for (step_i in seq_len(50)) {
        dd <- probit_box_derivs(s$mu_o + as.vector(s$Zo %*% bo), s$lo, s$hi)
        H <- P - crossprod(s$Zo, dd$h * s$Zo)
        step <- solve(H, as.vector(crossprod(s$Zo, dd$g)) - as.vector(P %*% bo) + h)
        bo <- bo + step
        if (max(abs(step)) < 1e-8) break
      }
      dd <- probit_box_derivs(s$mu_o + as.vector(s$Zo %*% bo), s$lo, s$hi)
      H <- P - crossprod(s$Zo, dd$h * s$Zo)
      mode[[i]] <<- bo
      nodes <- sweep(grid %*% chol(solve(H)), 2, bo, "+")   # rows b_o + R'x, R'R = H^-1
      eta <- rep(s$mu_o, each = n_nodes) + as.vector(nodes %*% t(s$Zo))
      log_lik <- rowSums(matrix(probit_box_derivs(eta, rep(s$lo, each = n_nodes),
                                                  rep(s$hi, each = n_nodes))$logP, n_nodes))
      lw <- log_lik - 0.5 * rowSums((nodes %*% P) * nodes) + as.vector(nodes %*% h) + log_w_grid
      w <- exp(lw - max(lw))
      w <- w / sum(w)
      S[ro, ro] <- S[ro, ro] + crossprod(nodes, w * nodes)
      if (length(rc) > 0) {
        mc <- sweep(nodes %*% t(B), 2, a, "+")
        S[rc, ro] <- S[rc, ro] + crossprod(mc, w * nodes)
        S[rc, rc] <- S[rc, rc] + Ac + crossprod(mc, w * mc)
      }
    }
    if (length(rc) > 0) S[ro, rc] <- t(S[rc, ro])
    D_new <- S / length(setup$subjects)
    dimnames(D_new) <- dimnames(D)
    D_new
  }
  is_pd <- function(V) {
    ev <- eigen((V + t(V)) / 2, symmetric = TRUE, only.values = TRUE)$values
    min(ev) > 1e-10 * max(ev)
  }

  ### EM accelerated by SQUAREM (Varadhan and Roland, 2008): two EM updates,
  ### an extrapolation along them, then an EM update from there; plain EM
  ### where the extrapolation leaves the positive-definite matrices
  n_updates <- 0
  repeat {
    D1 <- em_update(D)
    D2 <- em_update(D1)
    r <- D1 - D
    v <- D2 - 2 * D1 + D
    alpha <- if (sum(v^2) > 0) min(-1, -sqrt(sum(r^2) / sum(v^2))) else -1
    D_extra <- D - 2 * alpha * r + alpha^2 * v
    if (!is_pd(D_extra)) D_extra <- D2
    D_new <- em_update(D_extra)
    n_updates <- n_updates + 3
    ### on the correlation scale, so the check does not depend on the
    ### biomarkers' units (an absolute change never got below the tolerance
    ### for a biomarker with variances in the thousands)
    scale <- sqrt(pmax(diag(D_new), .Machine$double.eps))
    change <- max(abs(D_new - D) / outer(scale, scale))
    D <- D_new
    if (verbose) message(sprintf("EM updates %d, max|D change| = %.6f", n_updates, change))
    if (change < tol.em) break
    if (n_updates >= max.iter) {
      warning(sprintf(paste0(
        "The EM estimate of the joint random-effects covariance (Sigma_fit) did not converge ",
        "within %d iterations (last max change %.3g, tolerance %g); returning the last iterate."),
        max.iter, change, tol.em), call. = FALSE)
      break
    }
  }
  list(D = D)
}

#' Copula EM by stochastic approximation (SAEM)
#'
#' @description Fallback of \code{longitudinalSubVarCopula()} for many ordinal
#' random effects. Each iteration draws every ordinal latent score from its
#' truncated normal given the subject's current random effects (one
#' dimension at a time, as they are independent given \eqn{b}), computes
#' \eqn{E[b b^T \mid y, z]} exactly, draws new random effects from their
#' Gaussian posterior, and updates \eqn{D} by stochastic approximation: step
#' 1 for \code{n_burnin} iterations, then \eqn{1/k}, averaging the next
#' \code{n_average}. The draws use a fixed seed, with the caller's random
#' number state restored, so a fit is reproducible.
#' @param D Starting covariance.
#' @param setup Output of \code{copula_estep_setup()}.
#' @param n_burnin,n_average Iterations with step 1 and with step \eqn{1/k}.
#' @inheritParams longitudinalSubVarCopula
#' @return A list with \code{D}.
#' @keywords internal
copula_D_saem <- function(D, setup, n_burnin = 100, n_average = 200, verbose = FALSE) {
  local_r_seed(1L)
  ro <- setup$ro
  rc <- setup$rc
  q <- nrow(D)
  b <- lapply(setup$subjects, function(s) rep(0, q))
  # each subject's random-effects design and residual precision, continuous
  # rows first, then the ordinal (latent) ones
  design <- lapply(setup$subjects, function(s) {
    Z <- matrix(0, length(s$resid_c) + length(s$mu_o), q)
    Z[seq_along(s$resid_c), rc] <- s$Zc
    Z[length(s$resid_c) + seq_along(s$mu_o), ro] <- s$Zo
    Rinv <- c(s$Rinv_c, rep(1, length(s$mu_o)))
    list(Z = Z, Rinv = Rinv, ZRZ = crossprod(Z, Rinv * Z))
  })
  Q <- D
  for (k in seq_len(n_burnin + n_average)) {
    Dinv <- pd_inverse(D)
    S <- matrix(0, q, q)
    for (i in names(setup$subjects)) {
      s <- setup$subjects[[i]]
      g <- design[[i]]
      z <- rtnorm_vec(s$mu_o + as.vector(s$Zo %*% b[[i]][ro]), s$lo, s$hi)
      # all random effects given the continuous data and the drawn latent scores
      A <- solve(g$ZRZ + Dinv)
      m <- as.vector(A %*% crossprod(g$Z, g$Rinv * c(s$resid_c, z - s$mu_o)))
      S <- S + A + tcrossprod(m)
      b[[i]] <- m + as.vector(t(chol(A)) %*% stats::rnorm(q))
    }
    gamma <- if (k <= n_burnin) 1 else 1 / (k - n_burnin)
    Q <- Q + gamma * (S / length(setup$subjects) - Q)
    D <- Q
    if (verbose) message(sprintf("SAEM iter %d", k))
  }
  list(D = D)
}

#' Inverse of a covariance matrix that may be close to singular
#'
#' @description By the eigendecomposition, with eigenvalues floored at
#' \code{1e-10} times the largest, so a covariance estimate approaching the
#' boundary (e.g. a random-effect variance tending to 0) still gives a finite,
#' positive-definite precision.
#' @param V Symmetric matrix.
#' @return Its (regularized) inverse.
#' @keywords internal
pd_inverse <- function(V) {
  e <- eigen((V + t(V)) / 2, symmetric = TRUE)
  values <- pmax(e$values, 1e-10 * max(e$values, 0), .Machine$double.xmin)
  e$vectors %*% (t(e$vectors) / values)
}

#' Vectorized draws from truncated standard-variance normals
#'
#' @description \eqn{N(mean, 1)} truncated to \eqn{(lower, upper]}, by
#' inversion on the side of the mean where the interval lies.
#' @param mean,lower,upper Vectors of equal length.
#' @return A numeric vector of draws.
#' @keywords internal
rtnorm_vec <- function(mean, lower, upper) {
  a <- lower - mean
  b <- upper - mean
  flip <- (a + b) / 2 > 0
  a2 <- ifelse(flip, -b, a)
  b2 <- ifelse(flip, -a, b)
  pa <- stats::pnorm(a2)
  pb <- stats::pnorm(b2)
  u <- pa + stats::runif(length(mean)) * pmax(pb - pa, 0)
  e <- pmin(pmax(stats::qnorm(u), a2), b2)
  mean + ifelse(flip, -e, e)
}

#' Fit a multivariate longitudinal sub-model with mixed continuous/ordinal
#' biomarkers (Gaussian-copula extension)
#'
#' @description Internal counterpart to \code{longitudinalSubGaussian()},
#' used by the \code{\link{longitudinalSub}} dispatcher whenever at least one
#' biomarker is categorical/ordinal. Each biomarker is still fit on its own
#' first, exactly following \code{longitudinalSubGaussian()}'s IFM strategy:
#' continuous biomarkers via \code{nlme::lme()} (byte-identical code path to
#' \code{longitudinalSubGaussian()}), ordinal biomarkers via
#' \code{fit_marginal_ordinal()} (\code{ordinal::clmm()} with a probit link).
#' The resulting per-biomarker fixed effects, thresholds, and random-effects
#' covariance blocks seed the same \code{D}/\code{beta}/\code{sigma2}
#' bookkeeping \code{longitudinalSubGaussian()} builds, which is then handed
#' to \code{longitudinalSubVarCopula()} (rather than
#' \code{longitudinalSubVar()}) to jointly re-estimate \code{D} across all
#' biomarkers, including the ordinal ones' latent Gaussian scores.
#' @keywords internal
longitudinalSubCopula <- function(data_fit_all, long_sub_fixed, long_sub_random, biomarker_type, M) {
  id <- as.character(nlme::splitFormula(long_sub_random[[1]], "|")[[2]])[2]

  lfit       <- list()
  mf.fixed   <- list()
  xlevels    <- list()
  yik        <- list()
  Xik        <- list()
  nk         <- vector(length = M)
  Xik.list   <- list()
  nik.list   <- list()
  Zik        <- list()
  Zik.list   <- list()
  thresholds <- vector("list", M)

  ### different biomarkers have different missing samples, we have to get the intersect of them
  unique_num <- list()
  for (m in 1:M) {
    data.fit.one <- data_fit_all[[m]]
    data.fit.one <- data.fit.one[!as.logical(rowSums(data.frame(is.na(data.fit.one[model_vars(long_sub_fixed[[m]], long_sub_random[[m]])])))), ]
    unique_num[[m]] <- unique(unlist(data.fit.one[id]))
  }
  if (M != 1) {
    all_biomarker_num <- unlist(unique_num[[1]])
    for (vec in unique_num[-1]) {
      all_biomarker_num <- intersect(all_biomarker_num, unlist(vec))
    }
  } else {
    all_biomarker_num <- unique_num[[1]]
  }

  for (m in 1:M) {
    data.fit.one <- data_fit_all[[m]]

    if (biomarker_type[m] == "continuous") {
      lfit[[m]] <- lme_with_retry(nlme::lme(fixed = long_sub_fixed[[m]], random = long_sub_random[[m]],
                                            data = data.fit.one, method = "ML",
                                            control = nlme::lmeControl(opt = "optim"), na.action = na.omit),
                                  as.character(long_sub_fixed[[m]][[2]]))
      lfit[[m]]$call$fixed <- eval(lfit[[m]]$call$fixed)
      ### see longitudinalSubGaussian(): levels from the full fitting data,
      ### and the EM design built from lme()'s own terms (training basis)
      xlevels[[m]]  <- training_xlevels(lfit[[m]]$terms, data.fit.one, long_sub_fixed[[m]])

      data.fit.one <- data.fit.one[!as.logical(rowSums(data.frame(is.na(data.fit.one[model_vars(long_sub_fixed[[m]], long_sub_random[[m]])])))), ]
      data.fit.one <- data.fit.one[unlist(data.fit.one[id]) %in% unlist(all_biomarker_num), ]

      mf.fixed[[m]] <- model.frame(lfit[[m]]$terms, data.fit.one[, all.vars(long_sub_fixed[[m]])],
                                   xlev = xlevels[[m]])
      yik[[m]] <- by(model.response(mf.fixed[[m]], "numeric"), factor(data.fit.one[[id]]), as.vector)
      Xik[[m]] <- data.frame("id2" = factor(data.fit.one[[id]]),
                             model.matrix(lfit[[m]]$terms, mf.fixed[[m]],
                                          contrasts.arg = lfit[[m]]$contrasts))
    } else {
      lfit[[m]] <- fit_marginal_ordinal(long_sub_fixed[[m]], long_sub_random[[m]], data.fit.one)

      data.fit.one <- data.fit.one[!as.logical(rowSums(data.frame(is.na(data.fit.one[model_vars(long_sub_fixed[[m]], long_sub_random[[m]])])))), ]
      data.fit.one <- data.fit.one[unlist(data.fit.one[id]) %in% unlist(all_biomarker_num), ]

      mf.fixed[[m]] <- model.frame(stats::terms(long_sub_fixed[[m]]), data.fit.one[, all.vars(long_sub_fixed[[m]])])
      xlevels[[m]]  <- .getXlevels(stats::terms(long_sub_fixed[[m]]), mf.fixed[[m]])
      ### category codes (1..K), not the latent score: longitudinalSubVarCopula()
      ### integrates over the latent score within each category's thresholds.
      ### numbered by the categories clmm() actually fit (lfit$y.levels, to
      ### which the thresholds belong): clmm() drops a level nobody has, and
      ### the factor's own codes then no longer matched the thresholds (an
      ### unused lowest or middle level shifted every category above it)
      y_codes <- match(as.character(model.response(mf.fixed[[m]])), lfit[[m]]$y.levels)
      yik[[m]] <- by(y_codes, factor(data.fit.one[[id]]), as.vector)

      beta_names_m <- names(lfit[[m]]$beta)
      full_mm <- model.matrix(long_sub_fixed[[m]], data.fit.one)
      Xik[[m]] <- data.frame("id2" = factor(data.fit.one[[id]]), full_mm[, beta_names_m, drop = FALSE])
      thresholds[[m]] <- lfit[[m]]$alpha
    }

    nk[m] <- nrow(Xik[[m]])
    Xik.list[[m]] <- by(Xik[[m]], Xik[[m]]$id2, function(u) as.matrix(u[, -1]))
    nik.list[[m]] <- by(Xik[[m]], Xik[[m]]$id2, nrow)

    ffk <- nlme::splitFormula(long_sub_random[[m]], "|")[[1]]
    Zik[[m]] <- data.frame("id2" = factor(data.fit.one[[id]]), model.matrix(ffk, data.fit.one))
    Zik.list[[m]] <- by(Zik[[m]], Zik[[m]]$id2, function(u) as.matrix(u[, -1]))
  }

  # Flatten lists to length = n
  yi <- sapply(names(yik[[1]]), function(i) {
    unlist(lapply(yik, "[[", i))
  }, USE.NAMES = TRUE, simplify = FALSE)

  Xi <- sapply(names(Xik.list[[1]]), function(i) {
    as.matrix(Matrix::bdiag(lapply(Xik.list, "[[", i)))
  }, USE.NAMES = TRUE, simplify = FALSE)

  Zi <- sapply(names(Zik.list[[1]]), function(i) {
    as.matrix(Matrix::bdiag(lapply(Zik.list, "[[", i)))
  }, USE.NAMES = TRUE, simplify = FALSE)

  Zit <- lapply(Zi, t)

  XtXi <- lapply(Xi, crossprod)
  XtX.inv <- solve(Reduce("+", XtXi))

  Xtyi <- mapply(function(x, y) crossprod(x, y), x = Xi, y = yi, SIMPLIFY = FALSE)
  XtZi <- mapply(function(x, z) crossprod(x, z), x = Xi, z = Zi, SIMPLIFY = FALSE)

  nik <- sapply(names(nik.list[[1]]), function(i) {
    unlist(lapply(nik.list, "[[", i))
  }, USE.NAMES = TRUE, simplify = FALSE)

  p <- sapply(1:M, function(i) ncol(Xik[[i]]) - 1)
  r <- sapply(1:M, function(i) ncol(Zik[[i]]) - 1)

  l <- list(yi = yi, Xi = Xi, Zi = Zi, Zit = Zit, nik = nik, yik = yik,
            Xik.list = Xik.list, Zik.list = Zik.list, XtX.inv = XtX.inv,
            Xtyi = Xtyi, XtZi = XtZi, p = p, r = r, M = M, n = n_subjects(yi), nk = nk)

  get_vc <- function(m) {
    if (biomarker_type[m] == "continuous") nlme::getVarCov(lfit[[m]]) else nlme::VarCorr(lfit[[m]])[[1]]
  }

  D <- Matrix::bdiag(lapply(1:M, function(m) {
    vc <- get_vc(m)
    matrix(vc, dim(vc))
  }))
  D <- as.matrix(D)
  D.names <- c()
  for (m in 1:M) {
    D.names.k <- paste0(rownames(get_vc(m)), "_", m)
    D.names <- c(D.names, D.names.k)
  }
  rownames(D) <- colnames(D) <- D.names

  beta.1 <- do.call("c", lapply(1:M, function(m) {
    if (biomarker_type[m] == "continuous") fixef(lfit[[m]]) else lfit[[m]]$beta
  }))
  names(beta.1) <- paste0(names(beta.1), "_", rep(1:M, p))

  ### sigma2 is fixed at 1 for ordinal biomarkers (the identification
  ### constraint for a probit latent variable) and is never updated by
  ### longitudinalSubVarCopula(), exactly mirroring how longitudinalSubVar()
  ### never updates sigma2 for continuous biomarkers either.
  sigma2 <- sapply(1:M, function(m) if (biomarker_type[m] == "continuous") lfit[[m]]$sigma^2 else 1)

  ### see longitudinalSubGaussian(): identity starting value if the
  ### per-biomarker D is singular, with a warning rather than a bare message
  D_start <- tryCatch({ solve(D); D }, error = function(e) {
    warning(sprintf(paste0(
      "The random-effects covariance from the separate per-biomarker fits is singular (%s); ",
      "starting the joint EM from the identity matrix instead. Sigma_fit may be unreliable: ",
      "consider simplifying the random-effects structure."), conditionMessage(e)), call. = FALSE)
    diag(1, dim(D)[1])
  })
  out <- longitudinalSubVarCopula(thetaLong = list("beta" = beta.1, "D" = D_start, "sigma2" = sigma2),
                                   l = l, biomarker_type = biomarker_type, thresholds = thresholds,
                                   tol.em = 1e-04, verbose = FALSE)

  Sigma_fit <- out$D

  long_fit_all <- list(lfit = lfit, Sigma_fit = Sigma_fit,
                        long_sub_fixed = long_sub_fixed, long_sub_random = long_sub_random,
                        xlevels = xlevels, biomarker_type = biomarker_type, thresholds = thresholds)

  class(long_fit_all) <- "longitudinalSub.BJM"
  return(long_fit_all)
}
