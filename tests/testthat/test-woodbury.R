# cov_factor()/cov_logdens() evaluate the Gaussian log-density of one
# patient's stacked observations with covariance A D A' + diag(r_diag),
# either through the Woodbury identity (never forming the N x N matrix) or
# densely. Both must agree with mvtnorm::dmvnorm() on the explicit matrix.

random_cov_parts <- function(N, r, seed, singular = FALSE) {
  set.seed(seed)
  A <- matrix(rnorm(N * r), N, r)
  L <- matrix(rnorm(r * r), r, r)
  D <- tcrossprod(L) / r
  # a random-effect variance collapsed to 0, as can happen during EM: D is
  # then only positive semi-definite and chol(D) fails
  if (singular) D[1, ] <- D[, 1] <- 0
  # different residual variances, as with several biomarkers
  r_diag <- runif(N, 0.2, 2)
  list(A = A, D = D, r_diag = r_diag,
       Sigma = A %*% D %*% t(A) + diag(r_diag, N))
}

expect_matches_dmvnorm <- function(parts, method) {
  N <- nrow(parts$A)
  E <- matrix(rnorm(N * 4), N, 4)
  fac <- cov_factor(parts$A, parts$D, parts$r_diag, method = method)
  expect_equal(fac$method, method)
  ref <- apply(E, 2, function(e) {
    mvtnorm::dmvnorm(e, mean = rep(0, N), sigma = parts$Sigma, log = TRUE)
  })
  expect_equal(cov_logdens(fac, E), ref, tolerance = 1e-10)
}

test_that("Woodbury and dense log-densities match mvtnorm::dmvnorm", {
  for (method in c("woodbury", "dense")) {
    expect_matches_dmvnorm(random_cov_parts(N = 40, r = 6, seed = 1), method)   # N > r
    expect_matches_dmvnorm(random_cov_parts(N = 5, r = 12, seed = 2), method)   # N < r
    expect_matches_dmvnorm(random_cov_parts(N = 1, r = 2, seed = 3), method)    # one observation
  }
})

test_that("a singular random-effects covariance is handled without inverting it", {
  parts <- random_cov_parts(N = 30, r = 4, seed = 4, singular = TRUE)
  expect_error(chol(parts$D))
  expect_matches_dmvnorm(parts, "woodbury")
})

test_that("auto picks Woodbury only when it is cheaper", {
  wide <- random_cov_parts(N = 30, r = 4, seed = 5)
  narrow <- random_cov_parts(N = 10, r = 4, seed = 5)  # N < 3 r
  expect_equal(cov_factor(wide$A, wide$D, wide$r_diag)$method, "woodbury")
  expect_equal(cov_factor(narrow$A, narrow$D, narrow$r_diag)$method, "dense")
  # the dense path does not need a positive residual variance
  expect_equal(cov_factor(wide$A, wide$D, c(0, wide$r_diag[-1]))$method, "dense")
})

test_that("a non-positive-definite covariance gives dmvnorm's -Inf/Inf", {
  A <- matrix(0, 2, 1)
  fac <- cov_factor(A, matrix(1), c(0, 0), method = "dense")
  expect_equal(cov_logdens(fac, cbind(c(1, 0), c(0, 0))), c(-Inf, Inf))
})

test_that("cov_logdens_means evaluates every row under every mean", {
  parts <- random_cov_parts(N = 20, r = 4, seed = 6)
  x <- matrix(rnorm(3 * 20), 3, 20)
  means <- list(rnorm(20), rnorm(20))
  fac <- cov_factor(parts$A, parts$D, parts$r_diag)
  out <- cov_logdens_means(fac, x, means)
  ref <- lapply(means, function(m) mvtnorm::dmvnorm(x, m, parts$Sigma, log = TRUE))
  expect_equal(out, ref, tolerance = 1e-10)
})
