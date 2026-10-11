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

test_that("interval-censored fits reject invalid input", {
  d <- simulate_interval_data(200)
  d$type <- stats::rbinom(200, 1, 0.5)
  # the event-type model may use T only if T has been filled in
  expect_error(survivalSub(d, Surv(L, R, type = "interval2") ~ x1, type ~ T + x1, event_time = "T"),
               "interval censoring leaves unknown")
  d$type[which(!is.na(d$R))[1]] <- 2
  expect_error(survivalSub(d, Surv(L, R, type = "interval2") ~ x1, type ~ x1, event_time = "T"),
               "0 or 1")
  expect_error(survivalSub(d, Surv(L, R, type = "interval2") ~ x1, NULL), "event_time")
  expect_error(survivalSub(d, Surv(L, R, type = "interval2") ~ x1, NULL, event_time = "T", df = 0), "df")
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

test_that("interval-censored predictions condition on being event-free at the last visit", {
  skip_on_cran()
  set.seed(5)
  n <- 150
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
    long[[i]] <- data.frame(id = i, year = obs, x = x[i],
                            y = 1 + 0.3 * obs - 0.2 * T[i] + 0.4 * x[i] + stats::rnorm(1, 0, 0.5) +
                              stats::rnorm(length(obs), 0, 0.3))
  }
  long <- do.call(rbind, long)
  fit <- suppressWarnings(fitIntervalBJM(surv, long, Surv(L, R, type = "interval2") ~ x, "Tev",
                                         y ~ year + Tev + x, ~ year | id, "year",
                                         n_burnin = 2, n_imputations = 1, seed = 1))
  ids <- unique(long$id[long$year > 0 & long$year < 4])[1:4]
  hist <- long[long$id %in% ids & long$year <= 4, ]
  hist$Tev <- NA_real_
  last_visit <- tapply(hist$year, hist$id, max)
  expect_true(all(last_visit < 4))
  pr <- function(d, s, h, b) {
    predictRisk(list(d), fit$long_fit_all, fit$survival_fit_all, s, h, "year", NULL, NULL,
                bandcount1 = b, bandcount2 = b)
  }

  # the event may already have happened between the last visit V and s = 4,
  # so P(V < T <= 6 | T > V) splits into "already, undetected" + "in (4, 6]"
  at_4 <- pr(hist, 4, 2, 40)
  expect_true(all(at_4$prob_undetected_1 > 0 & at_4$prob_undetected_1 < 1))
  expect_output(print(at_4), "Already, undetected")
  for (i in as.character(ids)) {
    from_V <- pr(hist[hist$id == i, ], last_visit[[i]], 6 - last_visit[[i]], 200)
    expect_equal(unname(from_V$risk_prob_1), unname(at_4$prob_undetected_1[i] + at_4$risk_prob_1[i]),
                 tolerance = 5e-3)
  }

  # predicting at the last visit itself leaves nothing undetected
  one <- hist[hist$id == ids[1], ]
  expect_equal(unname(pr(one, last_visit[[1]], 2, 20)$prob_undetected_1), 0)

  # the biomarker's predictive density still integrates to one
  bio <- dynamicPredictionBio(1, list(hist), fit$long_fit_all, fit$survival_fit_all, 4, 1, "year",
                              NULL, NULL, bandcount2 = 20, bandcount3 = 200)
  expect_equal(unname(colSums(bio$Y_density) * diff(bio$Y_all)[1]), rep(1, length(ids)),
               tolerance = 0.02)

  # simulateTrajectory() only supports predicting at the last visit
  expect_error(simulateTrajectory(list(one), fit$long_fit_all, fit$survival_fit_all, 4, 5, "year",
                                  NULL, NULL, n_sim = 5, bandcount2 = 20),
               "last visit")
  expect_no_error(simulateTrajectory(list(one), fit$long_fit_all, fit$survival_fit_all,
                                     last_visit[[1]], last_visit[[1]] + 1, "year",
                                     NULL, NULL, n_sim = 5, bandcount2 = 20, seed = 1))
})

test_that("strata() gives each stratum its own baseline with shared coefficients", {
  # two strata with different Weibull baselines, a common effect of x1
  set.seed(11)
  n <- 1200
  x1 <- stats::rnorm(n)
  g <- stats::rbinom(n, 1, 0.5)
  T <- ifelse(g == 1,
              (-log(stats::runif(n)) / (0.1 * exp(0.5 * x1)))^(1 / 1.5),
              (-log(stats::runif(n)) / (0.2 * exp(0.5 * x1)))^(1 / 0.8))
  L <- R <- numeric(n)
  for (i in seq_len(n)) {
    v <- c(0, cumsum(stats::runif(20, 0.5, 1.5)))
    v <- v[v < stats::runif(1, 2, 12)]
    if (T[i] > max(v)) {
      L[i] <- max(v)
      R[i] <- NA
    } else {
      k <- findInterval(T[i], v)
      L[i] <- v[k]
      R[i] <- v[k + 1]
    }
  }
  d <- data.frame(id = seq_len(n), x1 = x1, g = g, L = L, R = R)
  fit <- survivalSub(d, Surv(L, R, type = "interval2") ~ x1 + strata(g), NULL, event_time = "T")
  ic <- fit$ic_fit

  expect_named(ic$coefficients, "x1")
  expect_equal(unname(ic$coefficients), 0.5, tolerance = 0.2)
  expect_equal(ic$strata_levels, c("g=0", "g=1"))
  bh <- survival_cum_basehaz(fit)
  expect_named(bh, c("hazard", "time", "strata"))
  expect_setequal(unique(bh$strata), c("g=0", "g=1"))
  tt <- c(1, 3, 6)
  H0 <- function(s) cumulative_baseline_at(bh[bh$strata == s, c("hazard", "time")], tt)
  expect_equal(H0("g=1"), 0.1 * tt^1.5, tolerance = 0.3)
  expect_equal(H0("g=0"), 0.2 * tt^0.8, tolerance = 0.3)

  # prediction code reads each patient's own stratum
  expect_equal(survival_patient_strata(fit, d[1:4, ]), paste0("g=", d$g[1:4]))
  m <- marginalT(list(d[1:4, ]), list(long_sub_random = list(~ 1 | id)), fit, l_i = c(1, 2, 3))
  for (j in 1:4) {
    s <- paste0("g=", d$g[j])
    S <- exp(-cumulative_baseline_at(bh[bh$strata == s, c("hazard", "time")], c(1, 2, 3)) *
               exp(ic$coefficients * d$x1[j]))
    expect_equal(m[, j], -diff(S))
  }

  # a stratum not seen in fitting is an error, as for a stratified Cox model
  expect_error(marginalT(list(transform(d[1, ], g = 5)), list(long_sub_random = list(~ 1 | id)), fit,
                         l_i = c(1, 2)), "Stratum")
  expect_output(print(fit), "Strata            : g=0, g=1")
  expect_true(any(grepl("Spline knots \\(time scale\\) \\[g=1\\]", utils::capture.output(summary(fit)))))
})

# Interval-censored competing risks: all-cause Weibull PH T, event type
# D | T ~ logistic(-1 + 0.3 T + 0.5 x) known once detected, and
# y | T, D = 1 + 0.3 t - 0.2 T + 0.5 D + 0.4 x + random intercept and slope
# + noise.
simulate_interval_cr <- function(n, seed, gap = c(1, 3)) {
  set.seed(seed)
  x <- stats::rnorm(n)
  T <- (-log(stats::runif(n)) / (0.05 * exp(0.5 * x)))^(1 / 1.5)
  D <- stats::rbinom(n, 1, stats::plogis(-1 + 0.3 * T + 0.5 * x))
  surv <- data.frame(id = seq_len(n), x = x, L = NA_real_, R = NA_real_, type = NA_real_, trueT = T)
  long <- list()
  for (i in seq_len(n)) {
    v <- c(0, cumsum(stats::runif(30, gap[1], gap[2])))
    v <- v[v < 12]
    if (T[i] > max(v)) {
      surv$L[i] <- max(v)
      obs <- v
    } else {
      k <- findInterval(T[i], v)
      surv$L[i] <- v[k]
      surv$R[i] <- v[k + 1]
      surv$type[i] <- D[i]
      obs <- v[v <= T[i]]
    }
    b <- stats::rnorm(2, 0, c(0.5, 0.1))
    long[[i]] <- data.frame(id = i, year = obs, x = x[i], type = surv$type[i],
                            y = 1 + 0.3 * obs - 0.2 * T[i] + 0.5 * D[i] + 0.4 * x[i] +
                              b[1] + b[2] * obs + stats::rnorm(length(obs), 0, 0.3))
  }
  list(surv = surv, long = do.call(rbind, long))
}

test_that("fitIntervalBJM handles competing risks with the event type known at detection", {
  skip_on_cran()
  d <- simulate_interval_cr(250, seed = 21)
  fit <- suppressWarnings(fitIntervalBJM(d$surv, d$long, Surv(L, R, type = "interval2") ~ x, "Tev",
                                         y ~ year + Tev + type + x, ~ year | id, "year",
                                         form_conditional_cr = type ~ Tev + x,
                                         n_burnin = 3, n_imputations = 2, seed = 1))
  observed <- d$surv[!is.na(d$surv$R), ]
  expect_true(all(fit$imputed_T > observed$L & fit$imputed_T <= observed$R))
  expect_length(fit$survival_fit_all_list, 2)
  expect_s3_class(fit$survival_fit_all$glm_fit, "glm")
  expect_equal(nrow(fit$survival_fit_all$glm_fit$data), nrow(observed))
  expect_true(all(c("cr.(Intercept)", "cr.Tev", "cr.x") %in% colnames(fit$trace)))
  # rough recovery (the formal check is a simulation study)
  expect_equal(unname(stats::coef(fit$survival_fit_all$glm_fit)["Tev"]), 0.3, tolerance = 0.5)
  expect_equal(unname(nlme::fixef(fit$long_fit_all$lfit[[1]])[c("Tev", "type")]), c(-0.2, 0.5),
               tolerance = 0.3)

  # predictions: per event type, P(V < T <= 6, D = d | T > V) predicted at the
  # last visit V splits into "already, undetected" + "in (4, 6]" at s = 4
  ids <- unique(d$long$id[d$long$year > 0 & d$long$year < 4])[1:3]
  hist <- d$long[d$long$id %in% ids & d$long$year <= 4, ]
  hist$Tev <- NA_real_
  hist$type <- NA_real_
  last_visit <- tapply(hist$year, hist$id, max)
  pr <- function(h, s, hz, b) {
    predictRisk(list(h), fit$long_fit_all, fit$survival_fit_all, s, hz, "year", NULL, NULL,
                bandcount1 = b, bandcount2 = b)
  }
  at_4 <- pr(hist, 4, 2, 40)
  expect_false(is.null(at_4$prob_undetected_2))
  expect_output(print(at_4), "Cause 2 undetected")
  for (i in as.character(ids)) {
    from_V <- pr(hist[hist$id == i, ], last_visit[[i]], 6 - last_visit[[i]], 200)
    expect_equal(unname(from_V$risk_prob_1), unname(at_4$prob_undetected_1[i] + at_4$risk_prob_1[i]),
                 tolerance = 5e-3)
    expect_equal(unname(from_V$risk_prob_2), unname(at_4$prob_undetected_2[i] + at_4$risk_prob_2[i]),
                 tolerance = 5e-3)
  }

  bio <- dynamicPredictionBio(1, list(hist), fit$long_fit_all, fit$survival_fit_all, 4, 1, "year",
                              NULL, NULL, bandcount2 = 20, bandcount3 = 200)
  expect_equal(unname(colSums(bio$Y_density) * diff(bio$Y_all)[1]), rep(1, length(ids)),
               tolerance = 0.02)

  expect_error(fitIntervalBJM(d$surv, d$long, Surv(L, R, type = "interval2") ~ x, "Tev",
                              y ~ year + Tev + type + x, ~ year | id, "year",
                              form_conditional_cr = type ~ Tev + x, include_right_censored = TRUE),
               "not supported with competing risks")
})
