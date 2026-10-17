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

test_that("performancePlot and calibrationPlot evaluate an interval-censored fit", {
  skip_on_cran()
  set.seed(8)
  n <- 150
  x <- stats::rnorm(n)
  T <- (-log(stats::runif(n)) / (0.05 * exp(0.5 * x)))^(1 / 1.5)
  surv <- data.frame(id = seq_len(n), x = x, L = NA_real_, R = NA_real_)
  long <- list()
  for (i in seq_len(n)) {
    v <- c(0, cumsum(stats::runif(30, 0.5, 1.5)))
    v <- v[v < 10]
    if (T[i] > max(v)) {
      surv$L[i] <- max(v)
      obs <- v
    } else {
      k <- findInterval(T[i], v)
      surv$L[i] <- v[k]
      surv$R[i] <- v[k + 1]
      obs <- v[v <= T[i]]
    }
    b <- stats::rnorm(2, 0, c(0.5, 0.1))
    long[[i]] <- data.frame(id = i, year = obs, x = x[i],
                            y = 1 + 0.3 * obs - 0.2 * T[i] + 0.4 * x[i] + b[1] + b[2] * obs +
                              stats::rnorm(length(obs), 0, 0.3))
  }
  long <- do.call(rbind, long)
  fit <- suppressWarnings(fitIntervalBJM(surv, long, Surv(L, R, type = "interval2") ~ x, "Tev",
                                         y ~ year + Tev + x, ~ year | id, "year",
                                         n_burnin = 2, n_imputations = 1, seed = 1))
  ev <- merge(long, surv[, c("id", "L", "R")], by = "id")

  p_model <- performancePlot(ev, fit$long_fit_all, fit$survival_fit_all, c(2, 4), 2, "year", NULL, NULL,
                             bandcount1 = 15, bandcount2 = 25)
  p_ipcw <- performancePlot(ev, fit$long_fit_all, fit$survival_fit_all, c(2, 4), 2, "year", NULL, NULL,
                            bandcount1 = 15, bandcount2 = 25, interval_method = "ipcw")
  for (p in list(p_model, p_ipcw)) {
    expect_s3_class(p, "ggplot")
    expect_equal(nrow(p$data), 4)
    expect_true(all(p$data$value >= 0 & p$data$value <= 1, na.rm = TRUE))
  }
  # the model-based version counts every subject at risk (expected cases),
  # IPCW only those whose status in the window is certain
  expect_true(all(p_model$data$n_cases >= p_ipcw$data$n_cases))
  # subjects at risk at s: followed up past s with no event detected by s
  followup <- ifelse(is.na(surv$R), surv$L, surv$R)
  expect_equal(p_model$data$n_at_risk[p_model$data$landmark == 2][1], sum(followup > 2))

  cal <- calibrationPlot(ev, fit$long_fit_all, fit$survival_fit_all, c(2, 4), 2, "year", NULL, NULL,
                         n_groups = 3, bandcount1 = 15, bandcount2 = 25)
  expect_s3_class(cal, "ggplot")
  expect_equal(nrow(cal$data), 6)
  expect_true(all(cal$data$observed >= 0 & cal$data$observed <= 1))
  # the highest predicted-risk group has a higher observed risk than the lowest
  for (s in c(2, 4)) {
    g <- cal$data[cal$data$landmark == s, ]
    expect_lt(g$observed[g$group == 1], g$observed[g$group == 3])
  }
})
