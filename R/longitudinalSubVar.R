#' Variance-covariance matrix
#' Reference: package "lmm" and package "joineRML" function \code{mvlme}.
#'
#' The EM loop stops after \code{max.iter} iterations even if the relative
#' change in \code{D} is still above \code{tol.em} (previously it had no
#' limit, so a fit that never met the tolerance never returned), and warns
#' when that happens.
#' @keywords internal
longitudinalSubVar <- function(thetaLong, l, tol.em, verbose, max.iter = 1000) {
  
  # Multivariate longitudinal data
  yi <- l$yi
  Xi <- l$Xi
  XtX.inv <- l$XtX.inv
  Xtyi <- l$Xtyi
  XtZi <- l$XtZi
  Zi <- l$Zi
  Zit <- l$Zit
  nik <- l$nik
  yik <- l$yik
  Xik.list <- l$Xik.list
  Zik.list <- l$Zik.list
  n <- l$n    # number of subjects
  p <- l$p    # vector of fixed effect dims
  r <- l$r    # vector of random effect dims
  M <- l$M    # number of longitudinal markers
  nk <- l$nk  # vector of number of observations per outcome
  
  delta <- 1
  beta <- thetaLong$beta
  sigma2 <- thetaLong$sigma2
  iter <- 0
  while (delta > tol.em) {
    iter <- iter + 1
    if (iter > max.iter) {
      warning(sprintf(paste0(
        "The EM estimate of the joint random-effects covariance (Sigma_fit) did not converge ",
        "within %d iterations (last max relative change %.3g, tolerance %g); returning the ",
        "last iterate."), max.iter, delta, tol.em), call. = FALSE)
      break
    }
    
    # Input parameter estimates
    D <- thetaLong$D
    
    ### E step
    
    # Inverse-Sigma_i (error precision matrix; diagonal matrix)
    Sigmai.inv <- lapply(nik, function(i) {
      diag(x = rep(1 / sigma2, i), ncol = sum(i))
    })
    
    # MVN covariance matrix for [b | y]
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
    Ai <- mapply(FUN = function(zt, s, z) {
      solve((zt %*% s %*% z) + Dinv)
    },
    z = Zi, zt = Zit, s = Sigmai.inv,
    SIMPLIFY = FALSE)
    
    # MVN mean vector for [y | b]
    Eb <- mapply(function(a, z, s, y, X) {
      as.vector(a %*% (z %*% s %*% (y - X %*% beta)))
    },
    a = Ai, z = Zit, s = Sigmai.inv, y = yi, X = Xi,
    SIMPLIFY = FALSE)
    
    # E[bb^T]
    EbbT <- mapply(function(v, e) {
      v + tcrossprod(e)
    },
    v = Ai, e = Eb,
    SIMPLIFY = FALSE)

    # M-step
    # D
    D.new <- Reduce("+", EbbT) / n
    rownames(D.new) <- colnames(D.new) <- rownames(D)
    
    #-----------------------------------------------------
    
    thetaLong.new <- list("D" = D.new)
    
    # Relative parameter change
    delta <- sapply(c("D"), function(i) {
      abs(thetaLong[[i]] - thetaLong.new[[i]]) / (abs(thetaLong[[i]]) + 1e-03)
    })
    delta <- max(unlist(delta))
    
    thetaLong <- thetaLong.new
    
    if (verbose) {
      print(thetaLong.new)
    }
    
  }
  
  return(thetaLong.new)
  
}
