#' Factor one patient's joint longitudinal covariance
#'
#' @description The stacked (across every biomarker) observation vector of
#' one patient has covariance
#' \code{Sigma_all = A \%*\% D \%*\% t(A) + diag(r_diag)}: \code{A} is the
#' \code{N x r} random-effects design (\code{random_effects_design()}),
#' \code{D} the \code{r x r} random-effects covariance (\code{Sigma_fit})
#' and \code{r_diag} each observation's residual variance. This factors it
#' once so that \code{cov_logdens()} can evaluate the Gaussian log-density
#' at any number of residual vectors.
#'
#' With \code{method = "woodbury"}, \code{Sigma_all} is never formed. Writing
#' \code{D = G \%*\% t(G)} and \code{B = A \%*\% G}, the matrix determinant
#' lemma and the Woodbury identity give
#' \deqn{\log|\Sigma| = \sum \log r_j + \log|K|,\quad
#'   e'\Sigma^{-1}e = e'R^{-1}e - u'K^{-1}u,}
#' with \code{K = I + t(B) \%*\% R^{-1} \%*\% B} (\code{r x r}) and
#' \code{u = t(B) \%*\% R^{-1} \%*\% e}. Factoring costs
#' \code{O(N r^2 + r^3)} instead of \code{O(N^3)}, each density
#' \code{O(N r)} instead of \code{O(N^2)}, and memory \code{O(N r)} instead
#' of \code{O(N^2)}. \code{D^{-1}} is never needed: \code{G} is
#' \code{t(chol(D))}, or, when \code{D} is only positive semi-definite (a
#' random-effect variance collapsing to 0 during EM), its eigenvectors
#' scaled by the square roots of the non-negative eigenvalues, which still
#' reproduces \code{D} exactly.
#'
#' With \code{method = "dense"}, \code{Sigma_all} is formed and Cholesky
#' factored, as \code{mvtnorm::dmvnorm()} does. \code{"auto"} picks
#' \code{"woodbury"} only when every residual variance is positive and
#' \code{N >= 3 r}: in timings the two break even around \code{N = 2.5 r}
#' (the Woodbury steps carry larger constants), and dense is faster below. Both give the same density up to
#' floating-point rounding.
#'
#' @param A Random-effects design matrix, \code{N x r}.
#' @param D Random-effects covariance matrix, \code{r x r}.
#' @param r_diag Residual variances, length \code{N}.
#' @param method \code{"auto"}, \code{"woodbury"} or \code{"dense"}.
#' @return An object for \code{cov_logdens()}.
#' @keywords internal
cov_factor <- function(A, D, r_diag, method = c("auto", "woodbury", "dense")) {
  method <- match.arg(method)
  N <- nrow(A)
  if (method == "auto") {
    method <- if (N >= 3 * ncol(A) && all(r_diag > 0)) "woodbury" else "dense"
  }

  if (method == "dense") {
    Sigma_all <- A %*% D %*% t(A) + diag(r_diag, N)
    ### a non-positive-definite Sigma_all is handled in cov_logdens() the
    ### way mvtnorm::dmvnorm() handles it
    dec <- tryCatch(base::chol(Sigma_all), error = function(e) NULL)
    logdet <- if (is.null(dec)) NA_real_ else 2 * sum(log(diag(dec)))
    return(list(method = "dense", N = N, dec = dec, logdet = logdet))
  }

  G <- tryCatch(t(base::chol(D)), error = function(e) {
    eig <- eigen(D, symmetric = TRUE)
    eig$vectors %*% diag(sqrt(pmax(eig$values, 0)), ncol(D))
  })
  B <- A %*% G
  r_inv <- 1 / r_diag
  K <- crossprod(B, B * r_inv)
  diag(K) <- diag(K) + 1
  K_chol <- base::chol(K)
  list(method = "woodbury", N = N, B = B, r_inv = r_inv, K_chol = K_chol,
       logdet = sum(log(r_diag)) + 2 * sum(log(diag(K_chol))))
}

#' Gaussian log-density from a factored covariance
#'
#' @description Evaluates the \code{N}-variate normal log-density with
#' covariance factored by \code{cov_factor()} at each column of a residual
#' matrix (observation minus mean).
#' @param fac The result of \code{cov_factor()}.
#' @param E Residual matrix, \code{N x m} (or a length-\code{N} vector).
#' @return A length-\code{m} vector of log-densities. If the covariance is
#' not positive definite (\code{"dense"} only), \code{-Inf}, or \code{Inf}
#' for a zero residual, as \code{mvtnorm::dmvnorm()} returns.
#' @keywords internal
cov_logdens <- function(fac, E) {
  E <- as.matrix(E)
  if (fac$method == "dense") {
    if (is.null(fac$dec)) {
      return(ifelse(colSums(E != 0) == 0, Inf, -Inf))
    }
    quad <- colSums(backsolve(fac$dec, E, transpose = TRUE)^2)
  } else {
    RE <- E * fac$r_inv
    u <- backsolve(fac$K_chol, crossprod(fac$B, RE), transpose = TRUE)
    quad <- colSums(E * RE) - colSums(u^2)
  }
  -0.5 * (fac$N * log(2 * pi) + fac$logdet + quad)
}

#' Gaussian log-densities of many points under many means
#'
#' @description For each mean vector in \code{means}, the log-density of
#' every row of \code{x}, with covariance factored by \code{cov_factor()}.
#' Used by \code{conditionalYTBio()}/\code{conditionalYDTBio()}, where the
#' rows of \code{x} are the candidate biomarker values and the means are
#' the survival-time grid points; the covariance does not depend on either,
#' so it is factored once per patient.
#' @param fac The result of \code{cov_factor()}.
#' @param x Matrix with one row per point, \code{N} columns.
#' @param means List of length-\code{N} mean vectors.
#' @return A list the same length as \code{means}; element \code{k} holds
#' the log-densities of the rows of \code{x} under \code{means[[k]]}.
#' @keywords internal
cov_logdens_means <- function(fac, x, means) {
  if (is.vector(x)) x <- matrix(x, ncol = length(x))
  if (ncol(x) != fac$N) stop("x and the covariance have non-conforming size")
  tx <- t(x)
  lapply(means, function(mean) {
    mean <- c(mean)
    if (length(mean) != fac$N) stop("x and mean have non-conforming size")
    out <- cov_logdens(fac, tx - mean)
    names(out) <- rownames(x)
    out
  })
}
