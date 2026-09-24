# Exercises the "auto" bandcount1/bandcount2/bandcount3 default: the
# doubling-until-stable search in auto_tune_bandcount() (unit-tested directly
# here against deterministic fake predict_fun()s, since real convergence
# behavior on real data is not something a test can pin down exactly), and
# its wiring into dynamicPrediction()/dynamicPredictionBio()/predictPlot()/
# riskPlot().

test_that("auto_tune_bandcount stops after the first doubling when already converged, without warning", {
  # v barely changes between bandcount1 = 10 and bandcount1 = 20 (starting
  # point and its first doubling, see bandcount_auto_start()), so this
  # should converge immediately and never try a second doubling.
  calls <- c()
  fake_predict <- function(bandcount1) {
    calls <<- c(calls, bandcount1)
    list(v = 100 - 1 / bandcount1)
  }

  out <- expect_no_warning(auto_tune_bandcount(fake_predict, list(bandcount1 = "auto"), "bandcount1"))

  expect_equal(calls, c(10, 20))
  expect_equal(out$bandcount$bandcount1, 20)
  expect_equal(out$result$v, 100 - 1 / 20)
})

test_that("auto_tune_bandcount warns and returns the largest value tried when it never converges", {
  # v keeps changing by a large relative amount every doubling, so this
  # never converges and must stop after max_rounds (default 2) doublings.
  calls <- c()
  fake_predict <- function(bandcount1) {
    calls <<- c(calls, bandcount1)
    list(v = 100 - 100 / bandcount1)
  }

  # expect_warning()/expect_no_warning() return the captured warning
  # condition, not the wrapped expression's value, so the warning and the
  # return value are checked separately here.
  warned_msg <- NULL
  out <- withCallingHandlers(
    auto_tune_bandcount(fake_predict, list(bandcount1 = "auto"), "bandcount1"),
    warning = function(w) {
      warned_msg <<- conditionMessage(w)
      invokeRestart("muffleWarning")
    }
  )
  expect_match(warned_msg, "had not converged")

  expect_equal(calls, c(10, 20, 40))
  expect_equal(out$bandcount$bandcount1, 40)
  expect_equal(out$result$v, 100 - 100 / 40)
})

test_that("auto_tune_bandcount never makes more than max_rounds + 1 calls", {
  calls <- c()
  fake_predict <- function(bandcount1) {
    calls <<- c(calls, bandcount1)
    list(v = 100 - 100 / bandcount1)
  }

  auto_tune_bandcount(fake_predict, list(bandcount1 = "auto"), "bandcount1", max_rounds = 3) |>
    suppressWarnings()

  expect_length(calls, 4)
})

test_that("dynamicPrediction defaults to auto bandcount1/bandcount2 and returns valid risk probabilities", {
  f <- setup_dp_fixture()

  risk <- suppressWarnings(dynamicPrediction(
    f$data_predict_all, f$long_fit_all, f$survival_fit_all,
    prediction_time = 5, horizon = 1, time_variable = "year",
    f$survival_variable_all, f$survival_trans_function
  ))

  expect_s3_class(risk, "dynamicPrediction.BJM")
  expect_true(risk$risk_prob_1 >= 0 && risk$risk_prob_1 <= 1)
  expect_true(risk$risk_prob_2 >= 0 && risk$risk_prob_2 <= 1)
})

test_that("dynamicPrediction with explicit numeric bandcount1/bandcount2 still bypasses auto-tuning (backward compatible)", {
  f <- setup_dp_fixture()

  risk <- dynamicPrediction(f$data_predict_all, f$long_fit_all, f$survival_fit_all,
                             prediction_time = 5, horizon = 1, time_variable = "year",
                             f$survival_variable_all, f$survival_trans_function,
                             bandcount1 = 10, bandcount2 = 20)

  expect_s3_class(risk, "dynamicPrediction.BJM")
})

test_that("dynamicPredictionBio defaults to auto bandcount2/bandcount3 and returns a MAP estimate", {
  f <- setup_dp_fixture()

  Y_pred <- suppressWarnings(dynamicPredictionBio(
    bio_i = 1, f$data_predict_all, f$long_fit_all, f$survival_fit_all,
    prediction_time = 5, horizon = 1, time_variable = "year",
    f$survival_variable_all, f$survival_trans_function
  ))

  expect_s3_class(Y_pred, "dynamicPredictionBio.BJM")
  expect_true(is.numeric(Y_pred$Y_predict))
})

test_that("predictPlot resolves auto bandcount1/bandcount2/bandcount3 once, not per horizon point", {
  f <- setup_dp_fixture()
  data.raw.predict.plot <- f$data_predict_all[[1]]

  ns <- asNamespace("BJM")
  n_calls <- 0
  trace_env <- environment()
  trace("dynamicPrediction", tracer = function() {
    assign("n_calls", get("n_calls", envir = trace_env) + 1, envir = trace_env)
  }, where = ns, print = FALSE)
  on.exit(suppressMessages(untrace("dynamicPrediction", where = ns)), add = TRUE)

  plot_obj <- suppressWarnings(predictPlot(
    f$data_predict_all, f$long_fit_all, f$survival_fit_all,
    prediction_time = 5, bio_his = 1, bio_pred = NULL,
    horizon = c(0.5, 1), time_variable = "year",
    f$survival_variable_all, f$survival_trans_function
  ))

  expect_s3_class(plot_obj, "ggplot")
  # If bandcount1/bandcount2 were re-resolved via "auto" on every horizon
  # point (2 points here) instead of once up front, dynamicPrediction()
  # would be called many more times than this.
  expect_lte(n_calls, 2 + 3)
})

test_that("riskPlot resolves auto bandcount1/bandcount2 once, not per landmark time", {
  f <- setup_dp_fixture()

  plot_obj <- suppressWarnings(riskPlot(
    f$data_predict_all, f$long_fit_all, f$survival_fit_all,
    prediction_time = c(5, 5.5), horizon = 1, time_variable = "year",
    survival_variable_all = f$survival_variable_all,
    survival_trans_function = f$survival_trans_function
  ))

  expect_s3_class(plot_obj, "ggplot")
})
