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

#' Impute the latent Gaussian score for an ordinal biomarker (inner E-step)
#'
#' @description The Gaussian-copula extension represents each ordinal
#' biomarker observation as a censored/interval-observed draw from an
#' underlying continuous latent Gaussian variable (the classical
#' Albert & Chib probit data-augmentation representation): the category
#' actually observed only tells us the latent score fell in the interval
#' implied by the fitted thresholds, not its exact value. Before each outer
#' E-step/M-step update of the shared random-effects covariance matrix
#' \code{D} (see \code{longitudinalSubVarCopula()}), this "inner" E-step
#' plugs in, for every ordinal observation, the mean of that latent score's
#' truncated-normal conditional distribution given its observed category and
#' the current fixed-effect and (from the previous outer E-step) predicted
#' random-effect contribution.
#'
#' This is a deterministic moment-matching plug-in (not a full Bayesian/MCEM
#' draw), which is a documented v1 simplification: it is the mean-field
#' approximation classically used to bootstrap probit/ordinal EM algorithms,
#' but does not propagate the extra uncertainty a stochastic draw would.
#' @param m Index of the ordinal biomarker being imputed.
#' @param l The per-subject design-matrix bookkeeping list built by
#'   \code{longitudinalSubCopula()} (same shape as the \code{l} passed to
#'   \code{longitudinalSubVar()}); \code{l$yik[[m]]} holds each subject's
#'   observed category codes (\code{1..K}) for biomarker \code{m}.
#' @param beta Current stacked fixed-effect vector (all biomarkers).
#' @param alpha_m Numeric vector of \code{K - 1} cumulative-link thresholds
#'   for biomarker \code{m} (\code{lfit[[m]]$alpha}).
#' @param Eb_prev Named list of each subject's posterior mean random-effects
#'   vector from the previous outer E-step, or \code{NULL} on the first
#'   iteration (random-effect contribution is then treated as 0).
#' @param beta_idx_m Integer indices of \code{beta} belonging to biomarker
#'   \code{m}.
#' @param r_idx_m Integer indices of the stacked random-effects vector
#'   belonging to biomarker \code{m}.
#' @return A named list (one element per subject id, matching
#'   \code{l$yik[[m]]}) of imputed latent-score vectors.
#' @keywords internal
impute_latent_ordinal <- function(m, l, beta, alpha_m, Eb_prev, beta_idx_m, r_idx_m) {
  ids <- names(l$yik[[m]])
  full_alpha <- c(-Inf, alpha_m, Inf)
  beta_m <- beta[beta_idx_m]

  out <- vector("list", length(ids))
  names(out) <- ids
  for (i in ids) {
    X_i_m <- as.matrix(l$Xik.list[[m]][[i]])
    eta <- as.vector(X_i_m %*% beta_m)
    if (!is.null(Eb_prev) && length(r_idx_m) > 0) {
      Z_i_m <- as.matrix(l$Zik.list[[m]][[i]])
      b_i_m <- Eb_prev[[i]][r_idx_m]
      eta <- eta + as.vector(Z_i_m %*% b_i_m)
    }

    y_cat <- l$yik[[m]][[i]]
    lower <- full_alpha[y_cat]
    upper <- full_alpha[y_cat + 1]

    ### truncated-standard-normal moment-matching plug-in for the latent
    ### residual epsilon ~ N(0, 1) given a < epsilon < b (see description
    ### above); pnorm()/dnorm() correctly return 0/0/1 at +-Inf so the
    ### open-ended bottom/top categories need no special-casing beyond the
    ### denominator floor guarding against extreme-tail numerical underflow.
    a <- lower - eta
    b <- upper - eta
    denom <- pmax(stats::pnorm(b) - stats::pnorm(a), 1e-8)
    eps_hat <- (stats::dnorm(a) - stats::dnorm(b)) / denom

    out[[i]] <- eta + eps_hat
  }
  out
}

#' EM re-estimation of the shared random-effects covariance (Gaussian-copula
#' extension)
#'
#' @description Extends \code{longitudinalSubVar()}'s EM algorithm with one
#' additional "inner" E-step (\code{impute_latent_ordinal()}), run before the
#' existing "outer" E-step (posterior mean/second-moment of the random
#' effects) at every iteration -- an ECM (Expectation-Conditional-Maximization)
#' scheme. The outer E-step and M-step formulas are otherwise completely
#' unchanged from \code{longitudinalSubVar()}: only the stacked outcome
#' vector \code{yi} fed into them differs, since ordinal biomarkers'
#' observed category codes are replaced by the inner E-step's imputed latent
#' scores before every outer-step update. As in \code{longitudinalSubVar()},
#' \code{sigma2} is never updated (for ordinal biomarkers it stays fixed at
#' 1, the identification constraint for a probit latent variable); only
#' \code{D} is re-estimated.
#' @keywords internal
longitudinalSubVarCopula <- function(thetaLong, l, biomarker_type, thresholds,
                                      tol.em = 1e-04, max.iter = 500, verbose = FALSE) {
  beta   <- thetaLong$beta
  D      <- thetaLong$D
  sigma2 <- thetaLong$sigma2
  n      <- l$n
  p      <- l$p
  r      <- l$r

  ord_idx <- which(biomarker_type == "ordinal")

  beta_bounds <- cumsum(c(0, p))
  r_bounds    <- cumsum(c(0, r))
  beta_idx <- lapply(seq_along(p), function(m) seq_len(p[m]) + beta_bounds[m])
  r_idx    <- lapply(seq_along(r), function(m) seq_len(r[m]) + r_bounds[m])

  yik_current <- l$yik
  Eb_prev <- NULL
  D_current <- D

  for (iter in seq_len(max.iter)) {
    if (length(ord_idx) > 0) {
      for (m in ord_idx) {
        yik_current[[m]] <- impute_latent_ordinal(
          m, l, beta, thresholds[[m]], Eb_prev, beta_idx[[m]], r_idx[[m]]
        )
      }
    }

    yi <- sapply(names(yik_current[[1]]), function(i) {
      unlist(lapply(yik_current, "[[", i))
    }, USE.NAMES = TRUE, simplify = FALSE)

    Dinv <- solve(D_current)
    Sigmai.inv <- lapply(l$nik, function(i) diag(x = rep(1 / sigma2, i), ncol = sum(i)))
    Ai <- mapply(function(zt, s, z) solve((zt %*% s %*% z) + Dinv),
                 z = l$Zi, zt = l$Zit, s = Sigmai.inv, SIMPLIFY = FALSE)
    Eb <- mapply(function(a, z, s, y, X) as.vector(a %*% (z %*% s %*% (y - X %*% beta))),
                 a = Ai, z = l$Zit, s = Sigmai.inv, y = yi, X = l$Xi, SIMPLIFY = FALSE)
    EbbT <- mapply(function(v, e) v + tcrossprod(e), v = Ai, e = Eb, SIMPLIFY = FALSE)
    D.new <- Reduce("+", EbbT) / n

    diff <- max(abs(D.new - D_current))
    D_current <- D.new
    Eb_prev <- Eb

    if (verbose) message(sprintf("ECM iter %d, max|D change| = %.6f", iter, diff))
    if (diff < tol.em) break
    if (iter == max.iter) {
      warning(sprintf(paste0(
        "The ECM estimate of the joint random-effects covariance (Sigma_fit) did not converge ",
        "within %d iterations (last max change %.3g, tolerance %g); returning the last iterate."),
        max.iter, diff, tol.em), call. = FALSE)
    }
  }

  list(D = D_current, Eb = Eb_prev, yik = yik_current)
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
    data.fit.one <- data.fit.one[!as.logical(rowSums(data.frame(is.na(data.fit.one[all.vars(long_sub_fixed[[m]])])))), ]
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
      lfit[[m]] <- nlme::lme(fixed = long_sub_fixed[[m]], random = long_sub_random[[m]],
                             data = data.fit.one, method = "ML",
                             control = nlme::lmeControl(opt = "optim"), na.action = na.omit)
      lfit[[m]]$call$fixed <- eval(lfit[[m]]$call$fixed)
      ### see longitudinalSubGaussian(): levels from the full fitting data,
      ### and the EM design built from lme()'s own terms (training basis)
      xlevels[[m]]  <- training_xlevels(lfit[[m]]$terms, data.fit.one, long_sub_fixed[[m]])

      data.fit.one <- data.fit.one[!as.logical(rowSums(data.frame(is.na(data.fit.one[all.vars(formula(lfit[[m]]))])))), ]
      data.fit.one <- data.fit.one[unlist(data.fit.one[id]) %in% unlist(all_biomarker_num), ]

      mf.fixed[[m]] <- model.frame(lfit[[m]]$terms, data.fit.one[, all.vars(long_sub_fixed[[m]])],
                                   xlev = xlevels[[m]])
      yik[[m]] <- by(model.response(mf.fixed[[m]], "numeric"), factor(data.fit.one[[id]]), as.vector)
      Xik[[m]] <- data.frame("id2" = factor(data.fit.one[[id]]),
                             model.matrix(lfit[[m]]$terms, mf.fixed[[m]],
                                          contrasts.arg = lfit[[m]]$contrasts))
    } else {
      lfit[[m]] <- fit_marginal_ordinal(long_sub_fixed[[m]], long_sub_random[[m]], data.fit.one)

      data.fit.one <- data.fit.one[!as.logical(rowSums(data.frame(is.na(data.fit.one[all.vars(long_sub_fixed[[m]])])))), ]
      data.fit.one <- data.fit.one[unlist(data.fit.one[id]) %in% unlist(all_biomarker_num), ]

      mf.fixed[[m]] <- model.frame(stats::terms(long_sub_fixed[[m]]), data.fit.one[, all.vars(long_sub_fixed[[m]])])
      xlevels[[m]]  <- .getXlevels(stats::terms(long_sub_fixed[[m]]), mf.fixed[[m]])
      ### category codes (1..K), not the latent score -- longitudinalSubVarCopula()'s
      ### inner E-step turns these into imputed latent Gaussian scores every iteration.
      yik[[m]] <- by(as.numeric(model.response(mf.fixed[[m]])), factor(data.fit.one[[id]]), as.vector)

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
    Zik.list[[m]] <- by(Zik[[m]], c(Zik[[m]]$id2), function(u) as.matrix(u[, -1]))
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
