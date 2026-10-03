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
  expect_named(sims, c("id", "sim", "year", "event_time", "status4", "serBilir", "albumin"))
  expect_true(all(sims$event_time > 5))
  expect_true(all(sims$status4 %in% c(0, 1)))
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
  sims <- simulate_fx(prediction_time = 5, times = 7, n_sim = 4000, bandcount2 = 80, seed = 1)
  risk <- predictRisk(fx$data_predict_all, fx$long_fit_all, fx$survival_fit_all,
                      prediction_time = 5, horizon = 2, time_variable = "year",
                      fx$survival_variable_all, fx$survival_trans_function,
                      bandcount1 = 40, bandcount2 = 80)
  expect_equal(mean(sims$event_time <= 7 & sims$status4 == 0), unname(risk$risk_prob_1), tolerance = 0.02)
  expect_equal(mean(sims$event_time <= 7 & sims$status4 == 1), unname(risk$risk_prob_2), tolerance = 0.01)
})

test_that("untruncated simulated biomarkers agree with predictLongitudinal()", {
  sims <- simulate_fx(prediction_time = 5, times = 7, n_sim = 4000, bandcount2 = 80,
                      truncate = FALSE, seed = 1)
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
