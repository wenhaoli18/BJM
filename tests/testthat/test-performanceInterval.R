# Performance measures for interval-censored fits (experimental).

test_that("npmle_interval reduces to Kaplan-Meier and Aalen-Johansen for (nearly) exact times", {
  set.seed(1)
  n <- 300
  T <- stats::rexp(n, 0.2)
  C <- stats::runif(n, 0, 15)
  V <- stats::runif(n, 0, 2)
  keep <- T > V & C > V
  T <- T[keep]
  C <- C[keep]
  V <- V[keep]
  time <- pmin(T, C)
  status <- as.numeric(T <= C)
  L <- ifelse(status == 1, time - 1e-6, time)
  R <- ifelse(status == 1, time, Inf)
  tt <- c(3, 6, 9)

  # left truncation at V
  fit <- npmle_interval(L, R, V, extra = tt)
  F_at <- function(x) c(0, cumsum(fit$mass[, 1]))[findInterval(x, fit$edges)]
  km <- survival::survfit(Surv(V, time, status) ~ 1)
  expect_equal(1 - (F_at(tt) - F_at(min(V))) / (1 - F_at(min(V))),
               summary(km, times = tt)$surv, tolerance = 1e-6)

  # competing risks, no truncation
  cause <- ifelse(status == 1, sample(1:2, length(status), TRUE, prob = c(0.3, 0.7)), NA)
  fit2 <- npmle_interval(L, R, rep(0, length(L)), cause, extra = tt)
  cum <- apply(fit2$mass, 2, cumsum)
  F2 <- function(x, k) c(0, cum[, k])[findInterval(x, fit2$edges)]
  aj <- summary(survival::survfit(Surv(time, factor(ifelse(status == 0, 0, cause), 0:2)) ~ 1), times = tt)
  expect_equal(sapply(tt, F2, k = 1), unname(aj$pstate[, 2]), tolerance = 1e-6)
  expect_equal(sapply(tt, F2, k = 2), unname(aj$pstate[, 3]), tolerance = 1e-6)
})

test_that("window_status_probs splits an interval by the posterior mass", {
  # uniform posterior on (0, 10], one event type
  post <- list(edges = 0:10, mass = matrix(0.1, 10, 1))
  # event in (2, 6], window (4, 5]: (4, 5] is a quarter of the interval (a
  # case), (5, 6] another quarter (event-free at the window's end, a control)
  p <- window_status_probs(post, L = 2, R = 6, cause = NA, s = 4, e = 5, k = 1)
  expect_equal(unname(p), c(0.25, 0.25))
  # certain case and certain control
  expect_equal(unname(window_status_probs(post, 4, 5, NA, 3, 6, 1)), c(1, 0))
  expect_equal(unname(window_status_probs(post, 7, Inf, NA, 3, 6, 1)), c(0, 1))
  # no event by L = 5 (R = Inf), window (4, 7]: T uniform on (5, 10]
  expect_equal(unname(window_status_probs(post, 5, Inf, NA, 4, 7, 1)), c(0.4, 0.6))
  # with two event types, a detected event of the other type is never a case
  post2 <- list(edges = 0:10, mass = matrix(0.05, 10, 2))
  expect_equal(unname(window_status_probs(post2, 2, 6, 2L, 4, 5, 1)), c(0, 0.25))
})

test_that("the IPCW and model-based measures use only what they should", {
  risk <- c(0.9, 0.2, 0.6, 0.1)
  # case detected in (2, 4]; control seen at 5; straddling case; censored early
  L <- c(2.5, 5, 1, 3)
  R <- c(3.5, Inf, 3, Inf)
  m <- performance_metrics_interval_ipcw(risk, L, R, ifelse(is.finite(R), R, L), rep(NA, 4), 1, 2, 2)
  expect_equal(m$n_cases, 1)
  expect_equal(m$n_uncertain, 2)
  expect_equal(m$auc, 1)

  mb <- performance_metrics_interval_model(risk, p_case = c(1, 0, 1, 0), p_control = c(0, 1, 0, 1))
  expect_equal(mb$brier, mean(c(0.1, 0.2, 0.4, 0.1)^2))
  expect_equal(mb$auc, 1)
})
