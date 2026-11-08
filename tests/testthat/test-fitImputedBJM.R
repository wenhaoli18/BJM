# Right-censored study where the biomarker depends on the event time; follow-up
# ends at 6, so about half the subjects are censored.
sim_right_censored <- function(n = 80, seed = 5) {
  set.seed(seed)
  x <- stats::rnorm(n)
  T <- (-log(stats::runif(n)) / (0.05 * exp(0.5 * x)))^(1 / 1.5)
  C <- pmin(stats::runif(n, 3, 6), 6)
  surv <- data.frame(id = seq_len(n), x = x, Tobs = pmin(T, C), status = as.numeric(T <= C))
  long <- do.call(rbind, lapply(seq_len(n), function(i) {
    v <- seq(0, surv$Tobs[i], by = 1)
    data.frame(id = i, year = v, x = x[i], Tobs = surv$Tobs[i], status = surv$status[i],
               y = 1 + 0.3 * v + 2 * exp(-T[i] / 3) + 0.4 * x[i] + stats::rnorm(1, 0, 0.5) +
                 stats::rnorm(length(v), 0, 0.3))
  }))
  list(surv = surv, long = long)
}

test_that("fitImputedBJM imputes censored event times beyond censoring and uses everyone", {
  skip_on_cran()
  d <- sim_right_censored()
  # 80 subjects is too few for the Sigma_fit EM to meet its tolerance
  fit <- suppressWarnings(fitImputedBJM(d$surv, d$long, Surv(Tobs, status) ~ x,
                                        y ~ year + I(exp(-Tobs / 3)) + x, ~ year | id, "year",
                                        n_burnin = 2, n_imputations = 2, n_thin = 1, n_grid = 15,
                                        seed = 1))
  expect_s3_class(fit, "fitImputedBJM")
  cens <- d$surv[d$surv$status == 0, ]
  expect_equal(rownames(fit$imputed_T), as.character(cens$id))
  expect_true(all(fit$imputed_T > cens$Tobs))
  expect_null(fit$imputed_D)
  expect_equal(unname(fit$long_fit_all_list[[1]]$lfit[[1]]$dims$ngrps[1]), nrow(d$surv))
  expect_equal(nrow(fit$trace), 4)
  expect_true(isTRUE(fit$long_fit_all$imputed_censored))
  printed <- capture.output(print(fit))
  expect_true(any(grepl(sprintf("%d censored with T imputed", nrow(cens)), printed)))

  # the fit plugs into predictRisk() without the complete-case warning
  old_opt <- options(BJM.check_complete_case = TRUE)
  on.exit(options(old_opt), add = TRUE)
  history <- d$long[d$long$id == cens$id[1] & d$long$year <= 2, ]
  history$Tobs <- NA_real_
  expect_no_warning(r <- predictRisk(history, fit$long_fit_all, fit$survival_fit_all, 2, 2, "year",
                                     NULL, NULL, bandcount1 = 10, bandcount2 = 10))
  expect_true(r$risk_prob_1 > 0 && r$risk_prob_1 < 1)

  ll <- observedLongLogLik(fit, d$surv, d$long, "year", n_grid = 15)
  expect_true(is.finite(ll))
})

test_that("fitImputedBJM draws the event type of censored subjects under competing risks", {
  skip_on_cran()
  data(pbc3, envir = environment())
  ds <- pbc3[!duplicated(pbc3$id), ]
  fit <- suppressWarnings(fitImputedBJM(
    ds, pbc3, Surv(years, status3) ~ age + sex,
    serBilir ~ year + age + I(exp(-years / 5)) + status4 + status4:I(exp(-years / 5)),
    ~ year | id, "year", form_conditional_cr = status4 ~ years + age + sex,
    n_burnin = 1, n_imputations = 2, n_thin = 1, n_grid = 10, seed = 2))
  cens <- ds[ds$status3 == 0, ]
  expect_true(all(fit$imputed_T > cens$years))
  expect_true(all(fit$imputed_D %in% c(0, 1)))
  expect_equal(dim(fit$imputed_D), c(nrow(cens), 2))
  ll <- observedLongLogLik(fit, ds, pbc3, "year", n_grid = 10)
  expect_true(is.finite(ll))
})

test_that("fitImputedBJM rejects interval-censored outcomes and data without censoring", {
  d <- sim_right_censored(n = 40)
  expect_error(fitImputedBJM(d$surv[d$surv$status == 1, ], d$long[d$long$status == 1, ],
                             Surv(Tobs, status) ~ x, y ~ year + Tobs, ~ year | id, "year"),
               "No censored subjects")
  d$surv$R <- ifelse(d$surv$status == 1, d$surv$Tobs, NA)
  expect_error(fitImputedBJM(d$surv, d$long, Surv(Tobs, R, type = "interval2") ~ x,
                             y ~ year + Tobs, ~ year | id, "year"),
               "fitIntervalBJM")
})

test_that("predictRisk warns once when the longitudinal fit ignores many censored subjects", {
  fx <- setup_dp_fixture()
  old_opt <- options(BJM.check_complete_case = TRUE)
  old <- isTRUE(BJM:::bjm_state$warned_complete_case)
  on.exit({
    options(old_opt)
    assign("warned_complete_case", old, envir = BJM:::bjm_state)
  }, add = TRUE)
  assign("warned_complete_case", FALSE, envir = BJM:::bjm_state)
  call_predict <- function() {
    predictRisk(fx$data_predict_all, fx$long_fit_all, fx$survival_fit_all, 5, 2, "year",
                fx$survival_variable_all, fx$survival_trans_function,
                bandcount1 = 10, bandcount2 = 10)
  }
  expect_warning(call_predict(), "fitImputedBJM")
  expect_no_warning(call_predict())
})
