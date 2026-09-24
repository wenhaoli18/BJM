# setup_dp_fixture() is defined in helper-fixtures.R (auto-sourced by testthat
# and shared across all test-*.R files).

test_that("dynamicPrediction returns risk probabilities in [0, 1] with competing risks", {
  f <- setup_dp_fixture()

  risk <- dynamicPrediction(f$data_predict_all, f$long_fit_all, f$survival_fit_all,
                             prediction_time = 5, horizon = 1, time_variable = "year",
                             f$survival_variable_all, f$survival_trans_function,
                             bandcount1 = 10, bandcount2 = 20)

  expect_s3_class(risk, "dynamicPrediction.BJM")
  expect_length(risk, 2)
  expect_true(risk[[1]] >= 0 && risk[[1]] <= 1)
  expect_true(risk[[2]] >= 0 && risk[[2]] <= 1)
  expect_true((risk[[1]] + risk[[2]]) <= 1)
})

test_that("dynamicPrediction risk with a zero-length horizon is exactly 0", {
  f <- setup_dp_fixture()

  risk <- dynamicPrediction(f$data_predict_all, f$long_fit_all, f$survival_fit_all,
                             prediction_time = 5, horizon = 0, time_variable = "year",
                             f$survival_variable_all, f$survival_trans_function,
                             bandcount1 = 10, bandcount2 = 20)

  expect_equal(risk[[1]], 0)
  expect_null(risk[[2]])
})

test_that("dynamicPredictionBio returns a MAP estimate and a density grid", {
  f <- setup_dp_fixture()

  Y_pred <- dynamicPredictionBio(bio_i = 1, f$data_predict_all, f$long_fit_all, f$survival_fit_all,
                                  prediction_time = 5, horizon = 1, time_variable = "year",
                                  f$survival_variable_all, f$survival_trans_function,
                                  bandcount2 = 20, bandcount3 = 50)

  expect_s3_class(Y_pred, "dynamicPredictionBio.BJM")
  expect_length(Y_pred, 3)
  expect_true(is.numeric(Y_pred[[1]]))
  expect_equal(nrow(Y_pred[[2]]), length(unlist(Y_pred[[3]])))
  expect_true(all(Y_pred[[2]] >= 0))
})

test_that("print and summary methods run without error", {
  f <- setup_dp_fixture()
  risk <- dynamicPrediction(f$data_predict_all, f$long_fit_all, f$survival_fit_all,
                             prediction_time = 5, horizon = 1, time_variable = "year",
                             f$survival_variable_all, f$survival_trans_function,
                             bandcount1 = 10, bandcount2 = 20)

  expect_output(print(risk), "Dynamic Prediction")
  expect_output(summary(risk), "Summary Statistics")
})
