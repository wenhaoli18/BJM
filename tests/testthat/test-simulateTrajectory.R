fx <- setup_dp_fixture()

simulate_fx <- function(...) {
  simulateTrajectory(fx$data_predict_all, fx$long_fit_all, fx$survival_fit_all,
                     time_variable = "year", survival_variable_all = fx$survival_variable_all,
                     survival_trans_function = fx$survival_trans_function, ...)
}

test_that("simulateTrajectory() returns one row per patient, draw and time", {
  sims <- simulate_fx(prediction_time = 5, times = c(5, 6, 7), n_sim = 20, seed = 1)
  expect_s3_class(sims, "simulateTrajectory.BJM")
  expect_equal(nrow(sims), 20 * 3)
  expect_named(sims, c("id", "sim", "year", "event_time", "event", "status4", "serBilir", "albumin"))
  expect_true(all(sims$event_time > 5))
  expect_true(all(sims$event %in% c(0, 1)))
  expect_true(all(sims$status4[sims$event == 1] %in% c(0, 1)))
  expect_true(all(is.na(sims$status4[sims$event == 0])))
  ### values at or after the drawn event time are NA, the rest are not
  after <- sims$year >= sims$event_time
  expect_true(all(is.na(sims$serBilir[after])))
  expect_false(anyNA(sims$serBilir[!after]))
})

test_that("simulateTrajectory() is reproducible with a seed and leaves the RNG alone", {
  set.seed(42)
  before <- .Random.seed
  a <- simulate_fx(prediction_time = 5, times = 6, n_sim = 10, seed = 3)
  expect_identical(.Random.seed, before)
  b <- simulate_fx(prediction_time = 5, times = 6, n_sim = 10, seed = 3)
  expect_identical(a, b)
})

test_that("simulated event times agree with predictRisk()", {
  sims <- simulate_fx(prediction_time = 5, times = 7, n_sim = 4000, bandcount2 = 80,
                      max_event_time = Inf, seed = 1)
  risk <- predictRisk(fx$data_predict_all, fx$long_fit_all, fx$survival_fit_all,
                      prediction_time = 5, horizon = 2, time_variable = "year",
                      fx$survival_variable_all, fx$survival_trans_function,
                      bandcount1 = 40, bandcount2 = 80)
  expect_equal(mean(sims$event_time <= 7 & sims$status4 == 0), unname(risk$risk_prob_1), tolerance = 0.02)
  expect_equal(mean(sims$event_time <= 7 & sims$status4 == 1), unname(risk$risk_prob_2), tolerance = 0.01)
})

test_that("untruncated simulated biomarkers agree with predictLongitudinal()", {
  sims <- simulate_fx(prediction_time = 5, times = 7, n_sim = 4000, bandcount2 = 80,
                      max_event_time = Inf, truncate = FALSE, seed = 1)
  pl <- predictLongitudinal(bio_i = 1, fx$data_predict_all, fx$long_fit_all, fx$survival_fit_all,
                            prediction_time = 5, horizon = 2, time_variable = "year",
                            fx$survival_variable_all, fx$survival_trans_function,
                            bandcount2 = 80, bandcount3 = 200)
  w <- pl$Y_density[, 1] / sum(pl$Y_density[, 1])
  m <- sum(w * pl$Y_all)
  s <- sqrt(sum(w * (pl$Y_all - m)^2))
  expect_equal(mean(sims$serBilir), m, tolerance = 0.05 / abs(m))
  expect_equal(sd(sims$serBilir), s, tolerance = 0.1)
})

test_that("a patient with no biomarker history is simulated from the prior", {
  newp <- fx$data_predict_all[[1]][1, ]
  newp$serBilir <- NA
  newp$albumin <- NA
  newp$year <- 0
  newp$years <- NA
  sims <- simulateTrajectory(newp, fx$long_fit_all, fx$survival_fit_all, prediction_time = 0,
                             times = 0:2, time_variable = "year",
                             survival_variable_all = fx$survival_variable_all,
                             survival_trans_function = fx$survival_trans_function,
                             n_sim = 50, seed = 1)
  expect_equal(nrow(sims), 50 * 3)
  expect_true(all(sims$event_time > 0))
  expect_false(anyNA(sims$serBilir[sims$year < sims$event_time]))
})

test_that("simulateTrajectory() rejects times before prediction_time", {
  expect_error(simulate_fx(prediction_time = 5, times = c(4, 6)), "at least prediction_time")
})

test_that("event times past max_event_time are reported as event-free through it", {
  last_follow_up <- max(fx$survival_fit_all$coxph_fit$y[, 1])
  sims <- simulate_fx(prediction_time = 5, times = c(6, 20), n_sim = 500, seed = 1)
  expect_equal(attr(sims, "max_event_time"), last_follow_up)
  expect_true(all(sims$event_time <= last_follow_up))
  censored <- sims$event == 0
  expect_true(any(censored))
  expect_true(all(sims$event_time[censored] == last_follow_up))
  ### nothing is simulated past the end of follow-up
  expect_true(all(is.na(sims$serBilir[sims$year == 20])))

  expect_error(simulate_fx(prediction_time = 5, times = 6, max_event_time = 4), "later than prediction_time")
})

test_that("without history, the share event-free through max_event_time is the Cox survival", {
  new_patient <- data.frame(id = "A", year = 0, age = 50, sex = 1, serBilir = NA_real_, albumin = NA_real_)
  sims <- simulateTrajectory(new_patient, fx$long_fit_all, fx$survival_fit_all, prediction_time = 0,
                             times = 0, time_variable = "year",
                             survival_variable_all = fx$survival_variable_all,
                             survival_trans_function = fx$survival_trans_function,
                             n_sim = 4000, seed = 1)
  cox_surv <- summary(survival::survfit(fx$survival_fit_all$coxph_fit, newdata = new_patient),
                      times = attr(sims, "max_event_time"), extend = TRUE)$surv
  expect_equal(mean(sims$event == 0), cox_surv, tolerance = 0.03 / cox_surv)
})

test_that("synthetic patients drawn from the prior keep biomarkers near the observed range", {
  data(pbc3, envir = environment())
  new_patients <- data.frame(id = c("A", "B"), year = 0, age = c(45, 60), sex = c(0, 1),
                             serBilir = NA_real_, albumin = NA_real_)
  sims <- simulateTrajectory(new_patients, fx$long_fit_all, fx$survival_fit_all, prediction_time = 0,
                             times = 0, time_variable = "year",
                             survival_variable_all = fx$survival_variable_all,
                             survival_trans_function = fx$survival_trans_function,
                             n_sim = 1000, seed = 1)
  ### extrapolating the sub-model past follow-up put 5% of these below -8
  observed <- range(pbc3$serBilir[pbc3$year == 0])
  simulated <- stats::quantile(sims$serBilir, c(0.01, 0.99), na.rm = TRUE)
  expect_gt(simulated[[1]], observed[1] - 2)
  expect_lt(simulated[[2]], observed[2] + 2)
})

test_that("rtmvnorm_box()'s Gibbs fallback matches its rejection draws", {
  set.seed(1)
  S <- matrix(c(1, 0.6, 0.3, 0.6, 1, 0.5, 0.3, 0.5, 1), 3)
  lower <- c(0.5, -Inf, 1)
  upper <- c(Inf, 0, 2)
  rejection <- rtmvnorm_box(20000, c(0, 0, 0), S, lower, upper)
  gibbs <- rtmvnorm_box(20000, c(0, 0, 0), S, lower, upper, max_proposals = 0)
  expect_true(all(gibbs > lower & gibbs <= upper))
  expect_equal(rowMeans(gibbs), rowMeans(rejection), tolerance = 0.02)
  expect_equal(apply(gibbs, 1, sd), apply(rejection, 1, sd), tolerance = 0.03)
})

setup_ordinal_sim_fixture <- function() {
  data(pbc3, envir = environment())
  d <- pbc3[pbc3$status3 == 1, ]
  breaks <- stats::quantile(d$albumin, c(0, 1 / 3, 2 / 3, 1), na.rm = TRUE)
  to_cat <- function(dd) {
    dd$albumin_cat <- cut(dd$albumin, breaks = breaks, include.lowest = TRUE,
                          labels = c("low", "mid", "high"), ordered_result = TRUE)
    dd
  }
  d <- to_cat(d)
  long_fit_all <- longitudinalSub(list(d, d),
                                  list(serBilir ~ year + age + sex + years, albumin_cat ~ year + age + sex + years),
                                  list(~ year | id, ~ year | id))
  survival_fit_all <- survivalSub(pbc3[!duplicated(pbc3$id), ], Surv(years, status3) ~ age + sex, NULL)
  history <- to_cat(pbc3[pbc3$id == 2 & pbc3$year <= 3, ])
  list(long_fit_all = long_fit_all, survival_fit_all = survival_fit_all,
       data_predict_all = list(history, history))
}

### fit once, on first use, for every ordinal test below
ordinal_sim_fixture <- local({
  fixture <- NULL
  function() {
    if (is.null(fixture)) fixture <<- setup_ordinal_sim_fixture()
    fixture
  }
})

test_that("an ordinal biomarker is simulated as an ordered factor matching predictLongitudinal()", {
  testthat::skip_if_not_installed("ordinal")
  ox <- ordinal_sim_fixture()
  sims <- simulateTrajectory(ox$data_predict_all, ox$long_fit_all, ox$survival_fit_all,
                             prediction_time = 3, times = 6, time_variable = "year",
                             survival_variable_all = list(), survival_trans_function = list(),
                             n_sim = 4000, bandcount2 = 40, max_event_time = Inf, truncate = FALSE, seed = 1)
  expect_true(is.ordered(sims$albumin_cat))
  expect_equal(levels(sims$albumin_cat), c("low", "mid", "high"))

  pl <- predictLongitudinal(bio_i = NULL, ox$data_predict_all, ox$long_fit_all, ox$survival_fit_all,
                            prediction_time = 3, horizon = 3, time_variable = "year",
                            survival_variable_all = list(), survival_trans_function = list(),
                            bandcount2 = 40, bandcount3 = 200)
  probs <- pl[[2]]$Y_density[, 1] / sum(pl[[2]]$Y_density[, 1])
  expect_equal(as.numeric(prop.table(table(sims$albumin_cat))), probs, tolerance = 0.03)
  w <- pl[[1]]$Y_density[, 1] / sum(pl[[1]]$Y_density[, 1])
  expect_equal(mean(sims$serBilir), sum(w * pl[[1]]$Y_all), tolerance = 0.05)
})

test_that("with an ordinal history, simulated event times agree with predictRisk()", {
  testthat::skip_if_not_installed("ordinal")
  ox <- ordinal_sim_fixture()
  sims <- simulateTrajectory(ox$data_predict_all, ox$long_fit_all, ox$survival_fit_all,
                             prediction_time = 3, times = 6, time_variable = "year",
                             survival_variable_all = list(), survival_trans_function = list(),
                             n_sim = 20000, bandcount2 = 160, max_event_time = Inf, seed = 1)
  risk <- predictRisk(ox$data_predict_all, ox$long_fit_all, ox$survival_fit_all,
                      prediction_time = 3, horizon = 3, time_variable = "year",
                      survival_variable_all = list(), survival_trans_function = list(),
                      bandcount1 = 40, bandcount2 = 160)
  expect_equal(mean(sims$event_time <= 6), unname(risk$risk_prob_1), tolerance = 0.01 / risk$risk_prob_1)
})

test_that("a synthetic patient can be drawn from a fit with an ordinal biomarker", {
  testthat::skip_if_not_installed("ordinal")
  ox <- ordinal_sim_fixture()
  new_patient <- data.frame(id = "A", year = 0, age = 50, sex = 1, serBilir = NA_real_, albumin_cat = NA)
  sims <- simulateTrajectory(new_patient, ox$long_fit_all, ox$survival_fit_all, prediction_time = 0,
                             times = 0:3, time_variable = "year", survival_variable_all = list(),
                             survival_trans_function = list(), n_sim = 50, seed = 1)
  before_event <- sims$year < sims$event_time
  expect_false(anyNA(sims$albumin_cat[before_event]))
  expect_true(all(is.na(sims$albumin_cat[!before_event])))
})

test_that("event times are spread within their interval rather than rounded to its midpoint", {
  sims <- simulate_fx(prediction_time = 5, times = 6, n_sim = 4000, bandcount2 = 20,
                      max_event_time = Inf, seed = 1)
  ### a coarse grid used to bias this by the mass of the interval straddling 7
  risk <- predictRisk(fx$data_predict_all, fx$long_fit_all, fx$survival_fit_all,
                      prediction_time = 5, horizon = 2, time_variable = "year",
                      fx$survival_variable_all, fx$survival_trans_function,
                      bandcount1 = 40, bandcount2 = 20)
  expect_equal(mean(sims$event_time <= 7), unname(risk$risk_prob_1 + risk$risk_prob_2), tolerance = 0.02)
  expect_gt(length(unique(sims$event_time)), 1000)
})
