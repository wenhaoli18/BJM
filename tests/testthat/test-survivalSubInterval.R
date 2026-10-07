# Interval-censored survival sub-model (experimental).

# Weibull PH event times observed only between irregular visits.
simulate_interval_data <- function(n, seed = 1, gompertz = FALSE, first_visit = 0.5) {
  set.seed(seed)
  x1 <- stats::rnorm(n)
  x2 <- stats::rbinom(n, 1, 0.5)
  lp <- 0.5 * x1 - 0.7 * x2
  T <- if (gompertz) {
    log(1 + stats::rexp(n) / (0.08 * exp(lp))) / 0.25
  } else {
    (-log(stats::runif(n)) / (0.1 * exp(lp)))^(1 / 1.5)
  }
  C <- stats::runif(n, 2, 12)
  L <- R <- numeric(n)
  for (i in seq_len(n)) {
    v <- cumsum(c(stats::runif(1, first_visit, first_visit + 1), stats::runif(19, 0.5, 1.5)))
    v <- c(0, v[v < C[i]])
    if (T[i] > max(v)) {
      L[i] <- max(v)
      R[i] <- NA
    } else {
      k <- findInterval(T[i], v)
      L[i] <- v[k]
      R[i] <- v[k + 1]
    }
  }
  data.frame(id = seq_len(n), x1 = x1, x2 = x2, L = L, R = R)
}

test_that("an interval2 outcome is fit by the interval-censored PH model", {
  d <- simulate_interval_data(800)
  fit <- survivalSub(d, Surv(L, R, type = "interval2") ~ x1 + x2, NULL, event_time = "T")

  expect_s3_class(fit, "survivalSub.BJM")
  expect_null(fit$coxph_fit)
  expect_s3_class(fit$ic_fit, "icph.BJM")
  expect_true(is_interval_censored(fit))
  expect_identical(survival_time_variable(fit), "T")
  expect_equal(unname(fit$ic_fit$coefficients), c(0.5, -0.7), tolerance = 0.25)
  expect_equal(fit$ic_fit$n_right + fit$ic_fit$n_interval + fit$ic_fit$n_exact, 800)

  lp <- survival_lp(fit, d[1:3, ])
  expect_equal(lp, c(as.matrix(d[1:3, c("x1", "x2")]) %*% fit$ic_fit$coefficients))
  bh <- survival_cum_basehaz(fit)
  expect_named(bh, c("hazard", "time"))
  expect_false(is.unsorted(bh$hazard))
  expect_null(survival_patient_strata(fit, d))
})

test_that("the piecewise-constant hazard model matches its closed form on exact data", {
  # one piece and no covariates: exponential MLE = events / total time
  set.seed(3)
  t <- stats::rexp(300, 0.4)
  d <- data.frame(L = t, R = t)
  fit <- icphFit(Surv(L, R, type = "interval2") ~ 1, d, "T", baseline = "piecewise", df = 1)
  expect_equal(fit$lambda, 300 / sum(t), tolerance = 1e-5)
})

test_that("the spline baseline with df = 1 is a Weibull model and fits a Weibull exactly", {
  # exact Weibull times: log H0 = log(0.1) + 1.5 log t, so the MLE of the
  # spline's intercept and (unscaled) slope estimate those
  set.seed(4)
  t <- (-log(stats::runif(2000)) / 0.1)^(1 / 1.5)
  d <- data.frame(L = t, R = t)
  fit <- icphFit(Surv(L, R, type = "interval2") ~ 1, d, "T", baseline = "spline", df = 1)
  expect_equal(fit$gamma[1], log(0.1), tolerance = 0.05)
  expect_equal(fit$gamma[2] / fit$base$scale[1], 1.5, tolerance = 0.05)
  # and its log hazard matches the derivative of its cumulative hazard
  h <- exp(fit$base$loghaz(c(1, 3), fit$gamma))
  eps <- 1e-6
  dH <- (fit$base$cumhaz(c(1, 3) + eps, fit$gamma) - fit$base$cumhaz(c(1, 3) - eps, fit$gamma)) / (2 * eps)
  expect_equal(h, dH, tolerance = 1e-6)
})

test_that("the spline baseline recovers an increasing hazard that the piecewise one flattens", {
  d <- simulate_interval_data(800, seed = 7)
  spline <- survivalSub(d, Surv(L, R, type = "interval2") ~ x1 + x2, NULL, event_time = "T")
  piecewise <- survivalSub(d, Surv(L, R, type = "interval2") ~ x1 + x2, NULL, event_time = "T",
                           baseline = "piecewise")
  expect_identical(spline$ic_fit$baseline, "spline")
  # the baseline choice barely moves the regression coefficients
  expect_equal(spline$ic_fit$coefficients, piecewise$ic_fit$coefficients, tolerance = 0.02)
  tt <- c(1, 3, 6)
  true_H0 <- 0.1 * tt^1.5
  err <- function(f) abs(cumulative_baseline_at(f$ic_fit$cum_basehaz, tt) - true_H0)
  expect_lt(err(spline)[1], err(piecewise)[1])
  expect_equal(cumulative_baseline_at(spline$ic_fit$cum_basehaz, tt), true_H0, tolerance = 0.25)
})

test_that("interval-censored fits reject unsupported options", {
  d <- simulate_interval_data(200)
  d$type <- stats::rbinom(200, 1, 0.5)
  expect_error(survivalSub(d, Surv(L, R, type = "interval2") ~ x1, type ~ x1, event_time = "T"),
               "Competing risks")
  expect_error(survivalSub(d, Surv(L, R, type = "interval2") ~ x1, NULL), "event_time")
  expect_error(survivalSub(d, Surv(L, R, type = "interval2") ~ x1, NULL, event_time = "T", df = 0), "df")
  expect_error(survivalSub(d, Surv(L, R, type = "interval2") ~ x1 + strata(x2), NULL, event_time = "T"),
               "strata")
})

test_that("print, summary and plot work for an interval-censored fit", {
  d <- simulate_interval_data(200)
  fit <- survivalSub(d, Surv(L, R, type = "interval2") ~ x1 + x2, NULL, event_time = "T")
  expect_output(print(fit), "interval-censored")
  before <- sum(!is.na(d$R) & d$L == 0)
  expect_equal(fit$ic_fit$n_before_first_visit, before)
  expect_output(print(fit), sprintf("Events before first visit : %d \\(%.0f%% of events\\)",
                                    before, 100 * before / sum(!is.na(d$R))))
  s <- utils::capture.output(out <- summary(fit))
  expect_true(any(grepl("Spline knots", s)))
  pw <- survivalSub(d, Surv(L, R, type = "interval2") ~ x1 + x2, NULL, event_time = "T", baseline = "piecewise")
  expect_true(any(grepl("Baseline hazard by piece", utils::capture.output(summary(pw)))))
  expect_equal(rownames(out$ic_coefficients), c("x1", "x2"))
  expect_s3_class(plot(fit), "ggplot")
  expect_s3_class(plot(fit, which = "basehaz"), "ggplot")
})

test_that("fitIntervalBJM imputes event times inside each interval and fits f(Y | T)", {
  skip_on_cran()
  set.seed(5)
  n <- 80
  x <- stats::rnorm(n)
  T <- (-log(stats::runif(n)) / (0.05 * exp(0.5 * x)))^(1 / 1.5)
  surv <- data.frame(id = seq_len(n), x = x, L = NA_real_, R = NA_real_)
  long <- list()
  for (i in seq_len(n)) {
    v <- c(0, cumsum(stats::runif(30, 1, 3)))
    v <- v[v < 12]
    if (T[i] > max(v)) {
      surv$L[i] <- max(v)
      obs <- v
    } else {
      k <- findInterval(T[i], v)
      surv$L[i] <- v[k]
      surv$R[i] <- v[k + 1]
      obs <- v[v <= T[i]]
    }
    y <- 1 + 0.3 * obs - 0.2 * T[i] + 0.4 * x[i] + stats::rnorm(1, 0, 0.5) + stats::rnorm(length(obs), 0, 0.3)
    long[[i]] <- data.frame(id = i, year = obs, x = x[i], y = y)
  }
  long <- do.call(rbind, long)

  # 80 subjects is too few for the Sigma_fit EM to meet its tolerance
  fit <- suppressWarnings(fitIntervalBJM(surv, long, Surv(L, R, type = "interval2") ~ x, "Tev",
                        y ~ year + Tev + x, ~ year | id, "year",
                        survival_variable_all = list("Tev_1"),
                        survival_trans_function = list(function(t) abs(t - 1)),
                        n_burnin = 2, n_imputations = 2, n_grid = 15, seed = 1))

  expect_s3_class(fit, "fitIntervalBJM")
  observed <- surv[!is.na(surv$R), ]
  expect_equal(rownames(fit$imputed_T), as.character(observed$id))
  inside <- fit$imputed_T > observed$L & fit$imputed_T <= observed$R
  expect_true(all(inside))
  expect_length(fit$long_fit_all_list, 2)
  expect_s3_class(fit$pooled, "poolLongitudinalSub.BJM")
  expect_equal(nrow(fit$trace), 4)
  expect_true("Tev" %in% names(nlme::fixef(fit$long_fit_all$lfit[[1]])))

  # the fit plugs into the existing prediction functions
  hist <- long[long$id == observed$id[1] & long$year <= 2, ]
  hist$Tev <- NA_real_
  risk <- predictRisk(list(hist), fit$long_fit_all, fit$survival_fit_all,
                      prediction_time = 2, horizon = 2, time_variable = "year",
                      survival_variable_all = list("Tev_1"),
                      survival_trans_function = list(function(t) abs(t - 1)),
                      bandcount1 = 10, bandcount2 = 20)
  expect_true(risk$risk_prob_1 > 0 && risk$risk_prob_1 < 1)

  # several draws per subject, all inside the subject's interval
  bounds <- subject_intervals(surv, Surv(L, R, type = "interval2") ~ x, "id")
  bounds <- bounds[is.finite(bounds$R), ]
  long_filled <- fill_event_time(list(long[long$id %in% bounds$id, ]), "id", fit$imputed_T[, 1], "Tev",
                                 list("Tev_1"), list(function(t) abs(t - 1)))
  draws <- imputeEventTime(long_filled, fit$long_fit_all, fit$survival_fit_all, bounds, "year",
                           list("Tev_1"), list(function(t) abs(t - 1)), n_grid = 15, n_draws = 4)
  expect_equal(dim(draws), c(nrow(bounds), 4))
  expect_true(all(draws > bounds$L & draws <= bounds$R))

  # the seed makes the chain reproducible
  fit2 <- suppressWarnings(fitIntervalBJM(surv, long, Surv(L, R, type = "interval2") ~ x, "Tev",
                         y ~ year + Tev + x, ~ year | id, "year",
                         survival_variable_all = list("Tev_1"),
                         survival_trans_function = list(function(t) abs(t - 1)),
                         n_burnin = 2, n_imputations = 2, n_grid = 15, seed = 1))
  expect_identical(fit2$imputed_T, fit$imputed_T)
})

test_that("the spline baseline stays increasing before the first visit", {
  # With no interval endpoint before the first visit, the likelihood does
  # not stop log H0 from sloping downwards below the lower boundary knot,
  # which unconstrained fits occasionally did (H0 rising towards t = 0).
  # seeds/sizes on which the unconstrained fit did go non-monotone
  for (case in list(c(10, 200), c(72, 200), c(76, 300))) {
    d <- simulate_interval_data(case[2], seed = case[1], gompertz = TRUE, first_visit = 2)
    fit <- expect_no_warning(survivalSub(d, Surv(L, R, type = "interval2") ~ x1 + x2, NULL,
                                         event_time = "T"))
    H0 <- fit$ic_fit$base$cumhaz(c(0.1, 0.5, 1, 2, 4, 8), fit$ic_fit$gamma)
    expect_false(is.unsorted(H0))
    expect_true(fit$ic_fit$base$monotone(fit$ic_fit$gamma))
  }
  base <- rp_baseline(c(1, 2, 3, 5, 8, 13), df = 3)
  expect_false(base$monotone(c(0, -1, 0, 0)))
})

test_that("lower_tail = 'weibull' extrapolates before the first endpoint with the Weibull slope", {
  # Below the lower boundary knot (the earliest interval endpoint) the
  # spline is pure extrapolation that never enters the likelihood, so the
  # two lower_tail choices must give the same fit everywhere else.
  d <- simulate_interval_data(300, seed = 10, gompertz = TRUE, first_visit = 2)
  f <- Surv(L, R, type = "interval2") ~ x1 + x2
  anchored <- icphFit(f, d, "T", lower_tail = "weibull")
  own <- icphFit(f, d, "T", lower_tail = "spline")
  weibull <- icphFit(f, d, "T", df = 1)

  expect_equal(anchored$loglik, own$loglik)
  expect_equal(anchored$coefficients, own$coefficients)
  kmin <- exp(anchored$base$knots[1])
  above <- c(kmin, 3, 6)
  expect_equal(anchored$base$cumhaz(above, anchored$gamma), own$base$cumhaz(above, own$gamma))

  below <- c(0.5, 1)
  slope <- diff(log(anchored$base$cumhaz(below, anchored$gamma))) / diff(log(below))
  expect_equal(slope, weibull$gamma[2] / weibull$base$scale[1])
  h <- exp(anchored$base$loghaz(below, anchored$gamma))
  eps <- 1e-6
  dH <- (anchored$base$cumhaz(below + eps, anchored$gamma) -
           anchored$base$cumhaz(below - eps, anchored$gamma)) / (2 * eps)
  expect_equal(h, dH, tolerance = 1e-6)
})
