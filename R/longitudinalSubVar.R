#' Variance-covariance matrix
#' Reference: package "lmm" and package "joineRML" function \code{mvlme}.
#'
#' The EM is accelerated by SQUAREM (Varadhan and Roland, 2008), and stops
#' once the relative change in \code{D} over an accelerated cycle (three EM
#' updates) is at most \code{tol.em}, or after \code{max.iter} EM updates
#' with a warning.
#' @keywords internal
longitudinalSubVar <- function(thetaLong, l, tol.em, verbose, max.iter = 1000) {

  # Multivariate longitudinal data
  n <- l$n    # number of subjects
  beta <- thetaLong$beta
  sigma2 <- thetaLong$sigma2

  ### beta and sigma2 stay fixed, so each subject's Z' Sigma^-1 Z and
  ### Z' Sigma^-1 (y - X beta) are computed once
  ZSZ <- mapply(function(z, zt, nik_i) zt %*% (rep(1 / sigma2, nik_i) * z),
                z = l$Zi, zt = l$Zit, nik_i = l$nik, SIMPLIFY = FALSE)
  ZSr <- mapply(function(zt, y, X, nik_i) as.vector(zt %*% (rep(1 / sigma2, nik_i) * (y - X %*% beta))),
                zt = l$Zit, y = l$yi, X = l$Xi, nik_i = l$nik, SIMPLIFY = FALSE)

  ### one EM update: D_new = mean over subjects of E[b b^T | y]
  em_update <- function(D) {
    ### D can become singular (e.g. a random-effect variance collapsing to
    ### 0); fall back to the identity for this iteration's precision, but say
    ### so -- this used to be reported only as a bare message().
    Dinv <- tryCatch(solve(D), error = function(e) {
      warning(sprintf(paste0(
        "The joint random-effects covariance became singular during EM (%s); using the ",
        "identity matrix in its place for this iteration. Sigma_fit may be unreliable: ",
        "consider simplifying the random-effects structure."), conditionMessage(e)), call. = FALSE)
      diag(1, dim(D)[1])
    })
    EbbT <- mapply(function(zsz, zsr) {
      A <- solve(zsz + Dinv)          # Var(b | y)
      e <- as.vector(A %*% zsr)       # E[b | y]
      A + tcrossprod(e)
    }, zsz = ZSZ, zsr = ZSr, SIMPLIFY = FALSE)
    D_new <- Reduce("+", EbbT) / n
    dimnames(D_new) <- dimnames(D)
    D_new
  }
  is_pd <- function(V) {
    ev <- eigen((V + t(V)) / 2, symmetric = TRUE, only.values = TRUE)$values
    min(ev) > 1e-10 * max(ev)
  }

  ### EM accelerated by SQUAREM (Varadhan and Roland, 2008): two EM updates,
  ### an extrapolation along them, then an EM update from there (plain EM
  ### where the extrapolation is not positive definite). Plain EM needed
  ### hundreds of iterations, and hit the cap of 1000 on ordinary data.
  D <- thetaLong$D
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

    # Relative parameter change over the cycle
    delta <- max(abs(D - D_new) / (abs(D) + 1e-03))
    D <- D_new
    if (verbose) print(D)
    if (delta <= tol.em) break
    if (n_updates >= max.iter) {
      warning(sprintf(paste0(
        "The EM estimate of the joint random-effects covariance (Sigma_fit) did not converge ",
        "within %d iterations (last max relative change %.3g, tolerance %g); returning the ",
        "last iterate."), max.iter, delta, tol.em), call. = FALSE)
      break
    }
  }

  list("D" = D)
}
