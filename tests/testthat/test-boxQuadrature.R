# The ordinal box probability of a patient's observed categories is computed
# by adaptive Gauss-Hermite quadrature when one ordinal random effect is
# involved (see mixed_density_copula_prep()), and by mvtnorm::pmvnorm()
# otherwise.

reference_box <- function(lo, up, mean, sigma) {
  log(as.numeric(mvtnorm::pmvnorm(lo, up, mean, sigma = sigma,
                                  algorithm = mvtnorm::GenzBretz(maxpts = 5e5, abseps = 1e-14,
                                                                 releps = 1e-7))))
}

test_that("box_prob_quadrature matches a high-precision pmvnorm, also far in the tails", {
  set.seed(1)
  thr <- c(-Inf, 0, 1.2, 2.5, Inf)
  for (spread in c(0, 1.8)) {
    n_o <- 6
    l <- stats::rnorm(n_o, 0, 1.2)
    S <- tcrossprod(l) + diag(n_o)
    cats <- sample(1:4, n_o, TRUE)
    lo <- thr[cats]
    up <- thr[cats + 1]
    N <- 20
    m <- matrix((pmax(lo, -1) + pmin(up, 3.5)) / 2 + spread * stats::rnorm(n_o * N), n_o, N)
    quad <- box_prob_quadrature(m, matrix(lo, n_o, N), matrix(up, n_o, N), latent_factor(S))
    ref <- vapply(seq_len(N), function(k) reference_box(lo, up, m[, k], S), numeric(1))
    expect_lt(max(abs(exp(quad - ref) - 1)), 1e-5)
  }
})

test_that("latent_factor finds a single random effect and refuses the rest", {
  l <- c(0.5, 1, -0.3)
  expect_equal(abs(latent_factor(tcrossprod(l) + diag(3))), abs(l), tolerance = 1e-10)
  expect_equal(latent_factor(diag(3)), rep(0, 3))
  # two random effects, or a latent residual variance other than 1: pmvnorm
  Z <- cbind(1, c(0, 1, 2))
  expect_null(latent_factor(Z %*% diag(c(1, 0.2)) %*% t(Z) + diag(3)))
  expect_null(latent_factor(0.5 * diag(3)))
})

test_that("the mixed density by quadrature agrees with pmvnorm, and is reproducible", {
  set.seed(2)
  # two continuous and three ordinal rows sharing one random effect (+ a
  # continuous-only one)
  A <- cbind(c(1, 1, 0.8, 1.1, 0.9), c(0.5, 1, 0, 0, 0))
  Sigma_all <- A %*% diag(c(1.5, 0.3)) %*% t(A) + diag(c(0.4, 0.4, 1, 1, 1))
  cc <- 1:2
  oo <- 3:5
  lo <- c(NA, NA, -Inf, 0.5, 1.5)
  up <- c(NA, NA, 0.5, 1.5, Inf)
  y <- c(0.3, 1.1, NA, NA, NA)
  mu <- matrix(stats::rnorm(5 * 30, 0.6, 0.7), 5, 30)
  quad <- mixed_density_copula_prep(Sigma_all, cc, oo)(mu, y, lo, up)
  mc <- mixed_density_copula_prep(Sigma_all, cc, oo, box = "pmvnorm")(mu, y, lo, up)
  expect_equal(quad, mc, tolerance = 1e-3)
  expect_identical(mixed_density_copula_prep(Sigma_all, cc, oo)(mu, y, lo, up), quad)
  # one column at a time gives the same as all at once
  expect_equal(mixed_density_copula_prep(Sigma_all, cc, oo)(mu[, 7], y, lo, up), quad[7])
})

test_that("predictions with one ordinal random effect no longer depend on the random seed", {
  skip_if_not_installed("ordinal")
  data(pbc3, envir = environment())
  pbc3$edema <- factor(pbc3$edema, labels = c("none", "controlled", "persistent"), ordered = TRUE)
  survival_fit <- survivalSub(pbc3[!duplicated(pbc3$id), ], Surv(years, status3) ~ age + sex, NULL)
  long_fit <- suppressWarnings(longitudinalSub(
    pbc3[pbc3$status3 == 1, ],
    list(serBilir ~ year + age + years, edema ~ year + age + years),
    list(~ year | id, ~ 1 | id)))
  patient <- pbc3[pbc3$id %in% c(2, 5) & pbc3$year <= 3, ]
  patient$years <- NA
  pr <- function(seed) {
    set.seed(seed)
    predictRisk(list(patient, patient), long_fit, survival_fit, prediction_time = 3, horizon = 2,
                time_variable = "year", list(), list(), bandcount1 = 10, bandcount2 = 20)
  }
  expect_identical(pr(1), pr(2))
})

test_that("interpolate_log_box follows a smooth curve and is -Inf where the probability is 0", {
  x <- seq(-3, 3, length.out = 15)
  f <- function(x) stats::pnorm(x, log.p = TRUE) + stats::pnorm(-0.5 * x, log.p = TRUE)
  x_out <- seq(-3, 3, length.out = 200)
  expect_lt(max(abs(interpolate_log_box(x, f(x), x_out) - f(x_out))), 1e-3)
  y <- f(x)
  y[1:3] <- -Inf
  out <- interpolate_log_box(x, y, x_out)
  expect_true(all(out[x_out < x[4]] == -Inf))
  expect_true(all(is.finite(out[x_out >= x[4]])))
  expect_equal(interpolate_log_box(x, rep(-Inf, 15), x_out), rep(-Inf, 200))
})

test_that("a continuous biomarker's prediction with interpolated box probabilities matches computing every value", {
  skip_if_not_installed("ordinal")
  skip_on_cran()
  # the ordinal biomarker has a random intercept and slope: Monte Carlo box
  # probabilities, interpolated over the candidate values of serBilir
  fx <- setup_copula_predict_fixture()
  survival_fit_all <- survivalSub(fx$pbc3[!duplicated(fx$pbc3$id), ], Surv(years, status3) ~ age + sex, NULL)
  raw <- fx$pbc3[fx$pbc3$id %in% c(2, 5) & fx$pbc3$year <= 3, ]
  data_predict_all <- list(fx$to_cat(raw), fx$to_cat(raw))
  pred <- function() {
    predictLongitudinal(1, data_predict_all, fx$long_fit_all, survival_fit_all, prediction_time = 3,
                        horizon = 1, time_variable = "year", list(), list(), bandcount2 = 10,
                        bandcount3 = 60)
  }
  interpolated <- pred()
  every_value <- copula_logdens_grid
  local_mocked_bindings(copula_logdens_grid = function(...) every_value(..., n_exact = Inf))
  full <- pred()
  normalize <- function(d) sweep(d, 2, colSums(d), "/")
  expect_equal(normalize(interpolated$Y_density), normalize(full$Y_density), tolerance = 2e-3)
  expect_equal(interpolated$Y_predict, full$Y_predict, tolerance = 2e-3)
})
