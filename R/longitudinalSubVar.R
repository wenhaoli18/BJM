#' Variance-covariance matrix
#' Reference: package "lmm" and package "joineRML" function \code{mvlme}.
#' @keywords internal
longitudinalSubVar <- function(thetaLong, l, tol.em, verbose) {
  
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
  while (delta > tol.em) {
    
    # Input parameter estimates
    D <- thetaLong$D
    
    ### E step
    
    # Inverse-Sigma_i (error precision matrix; diagonal matrix)
    Sigmai.inv <- lapply(nik, function(i) {
      diag(x = rep(1 / sigma2, i), ncol = sum(i))
    })
    
    # MVN covariance matrix for [b | y]
    #if('try-error' %in% class(solve(D))) {
    #  Dinv = diag(1, 6, 6)
    #} else{
    #  Dinv <- solve(D)
    #}
    result.solve.D <- tryCatch({
      solve(D)
      # code that may produce an error
    }, error = function(e) {
      # code to handle the error, such as printing a message
      message("An error occurred: ", conditionMessage(e))
      NULL  # return NULL to indicate that an error occurred
    })
    
    if (!is.null(result.solve.D)) {
      # code to execute if no error occurs
      # use the 'result' variable here, if necessary
      # ...
      Dinv <- solve(D)
    } else {
      # code to execute if an error occurs
      # ...
      Dinv <- diag(1, dim(D)[1])
    }
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
