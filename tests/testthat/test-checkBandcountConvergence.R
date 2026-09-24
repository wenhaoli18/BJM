# Exercises checkBandcountConvergence(), the lightweight "compare base vs.
# scaled bandcount" diagnostic added as an alternative to purely manual
# bandcount1/bandcount2/bandcount3 tuning.

test_that("checkBandcountConvergence rejects a predict_fun that is not dynamicPrediction/dynamicPredictionBio", {
  f <- setup_dp_fixture()

  expect_error(
    checkBandcountConvergence(predictPlot, f$data_predict_all, f$long_fit_all, f$survival_fit_all,
                               prediction_time = 5, horizon = 1, time_variable = "year",
                               f$survival_variable_all, f$survival_trans_function,
                               bandcount_args = list(bandcount1 = 10)),
    "must be .dynamicPrediction. or .dynamicPredictionBio."
  )
})

test_that("checkBandcountConvergence rejects an empty/unnamed bandcount_args", {
  f <- setup_dp_fixture()

  expect_error(
    checkBandcountConvergence(dynamicPrediction, f$data_predict_all, f$long_fit_all, f$survival_fit_all,
                               prediction_time = 5, horizon = 1, time_variable = "year",
                               f$survival_variable_all, f$survival_trans_function,
                               bandcount_args = list()),
    "non-empty named list"
  )
})

test_that("checkBandcountConvergence rejects a bandcount_args name that isn't a bandcount argument of predict_fun", {
  f <- setup_dp_fixture()

  expect_error(
    checkBandcountConvergence(dynamicPrediction, f$data_predict_all, f$long_fit_all, f$survival_fit_all,
                               prediction_time = 5, horizon = 1, time_variable = "year",
                               f$survival_variable_all, f$survival_trans_function,
                               bandcount_args = list(bandcount3 = 50)),
    "not bandcount argument"
  )
})

test_that("checkBandcountConvergence rejects multiplier <= 1", {
  f <- setup_dp_fixture()

  expect_error(
    checkBandcountConvergence(dynamicPrediction, f$data_predict_all, f$long_fit_all, f$survival_fit_all,
                               prediction_time = 5, horizon = 1, time_variable = "year",
                               f$survival_variable_all, f$survival_trans_function,
                               bandcount_args = list(bandcount1 = 10, bandcount2 = 20), multiplier = 1),
    "`multiplier` must be greater than 1"
  )
})

test_that("checkBandcountConvergence runs dynamicPrediction at base and scaled bandcount1/bandcount2 and compares risk_prob_1/2", {
  f <- setup_dp_fixture()

  check <- checkBandcountConvergence(
    dynamicPrediction, f$data_predict_all, f$long_fit_all, f$survival_fit_all,
    prediction_time = 5, horizon = 1, time_variable = "year",
    f$survival_variable_all, f$survival_trans_function,
    bandcount_args = list(bandcount1 = 10, bandcount2 = 20)
  )

  expect_s3_class(check, "checkBandcountConvergence.BJM")
  expect_equal(check$base_bandcount, list(bandcount1 = 10, bandcount2 = 20))
  expect_equal(check$scaled_bandcount, list(bandcount1 = 20, bandcount2 = 40))
  expect_true(all(c("risk_prob_1", "risk_prob_2") %in% names(check$by_field)))
  expect_true(is.numeric(check$max_rel_diff))
  expect_true(is.logical(check$converged))
  expect_output(print(check), "Bandcount convergence check")
})

test_that("checkBandcountConvergence can isolate a single bandcount argument, holding the other fixed via ...", {
  f <- setup_dp_fixture()

  check <- checkBandcountConvergence(
    dynamicPrediction, f$data_predict_all, f$long_fit_all, f$survival_fit_all,
    prediction_time = 5, horizon = 1, time_variable = "year",
    f$survival_variable_all, f$survival_trans_function,
    bandcount1 = 10,
    bandcount_args = list(bandcount2 = 20)
  )

  expect_equal(check$base_bandcount, list(bandcount2 = 20))
  expect_equal(check$scaled_bandcount, list(bandcount2 = 40))
})

test_that("checkBandcountConvergence works with dynamicPredictionBio, comparing Y_predict but not the bandcount3-sized Y_density/Y_all grid", {
  f <- setup_dp_fixture()

  check <- checkBandcountConvergence(
    dynamicPredictionBio, bio_i = 1, f$data_predict_all, f$long_fit_all, f$survival_fit_all,
    prediction_time = 5, horizon = 1, time_variable = "year",
    f$survival_variable_all, f$survival_trans_function,
    bandcount2 = 20,
    bandcount_args = list(bandcount3 = 50)
  )

  expect_equal(check$base_bandcount, list(bandcount3 = 50))
  expect_equal(check$scaled_bandcount, list(bandcount3 = 100))
  expect_true("Y_predict" %in% names(check$by_field))
  expect_false("Y_all" %in% names(check$by_field))
  expect_false("Y_density" %in% names(check$by_field))
})
