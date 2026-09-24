# Exercises the friendly-error-message input validation added on top of the
# core pipeline functions. Each test checks that a common user mistake now
# fails fast with a specific, actionable message instead of a cryptic error
# from deep inside model-fitting/indexing code.

test_that("survivalSub rejects a non-data.frame data_survival_fitting", {
  expect_error(
    survivalSub(list(a = 1), Surv(years, status3) ~ age, NULL),
    "`data_survival_fitting` must be a data.frame"
  )
})

test_that("survivalSub rejects a form_marginal_surv referencing a missing variable", {
  data(pbc3, envir = environment())
  data_survival_fitting <- pbc3[!duplicated(pbc3$id), ]

  expect_error(
    survivalSub(data_survival_fitting, Surv(years, status3) ~ age + not_a_column, NULL),
    "not_a_column"
  )
})

test_that("survivalSub rejects a non-formula form_conditional_cr", {
  data(pbc3, envir = environment())
  data_survival_fitting <- pbc3[!duplicated(pbc3$id), ]

  expect_error(
    survivalSub(data_survival_fitting, Surv(years, status3) ~ age, "status4 ~ age"),
    "`form_conditional_cr` must be a formula"
  )
})

test_that("longitudinalSub rejects mismatched long_sub_fixed/long_sub_random lengths", {
  data(pbc3, envir = environment())
  long_sub_fixed <- list(serBilir ~ year + age, albumin ~ year + age)
  long_sub_random <- list(~ year | id)

  expect_error(
    longitudinalSub(list(pbc3, pbc3), long_sub_fixed, long_sub_random),
    "long_sub_fixed.*long_sub_random"
  )
})

test_that("longitudinalSub rejects a data_fit_all list of the wrong length", {
  data(pbc3, envir = environment())
  long_sub_fixed <- list(serBilir ~ year + age, albumin ~ year + age)
  long_sub_random <- list(~ year | id, ~ year | id)

  expect_error(
    longitudinalSub(list(pbc3), long_sub_fixed, long_sub_random),
    "`data_fit_all` has 1 element"
  )
})

test_that("longitudinalSub rejects a fixed-effects formula referencing a missing variable", {
  data(pbc3, envir = environment())

  expect_error(
    longitudinalSub(pbc3, serBilir ~ year + not_a_column, ~ year | id),
    "not_a_column"
  )
})

# Shared fixture reused from test-dynamicPrediction.R (sourced into the same
# environment by testthat).

test_that("dynamicPrediction rejects a long_fit_all/survival_fit_all swap", {
  f <- setup_dp_fixture()

  expect_error(
    dynamicPrediction(f$data_predict_all, f$survival_fit_all, f$long_fit_all,
                       prediction_time = 5, horizon = 1, time_variable = "year",
                       f$survival_variable_all, f$survival_trans_function,
                       bandcount1 = 10, bandcount2 = 20),
    "`long_fit_all` must be the output of longitudinalSub"
  )
})

test_that("dynamicPrediction accepts a bare data.frame for data_predict_all, reusing it for every biomarker", {
  f <- setup_dp_fixture()

  risk_bare <- dynamicPrediction(f$data_predict_all[[1]], f$long_fit_all, f$survival_fit_all,
                                  prediction_time = 5, horizon = 1, time_variable = "year",
                                  f$survival_variable_all, f$survival_trans_function,
                                  bandcount1 = 10, bandcount2 = 20)
  risk_list <- dynamicPrediction(f$data_predict_all, f$long_fit_all, f$survival_fit_all,
                                  prediction_time = 5, horizon = 1, time_variable = "year",
                                  f$survival_variable_all, f$survival_trans_function,
                                  bandcount1 = 10, bandcount2 = 20)

  expect_equal(risk_bare$risk_prob_1, risk_list$risk_prob_1)
  expect_equal(risk_bare$risk_prob_2, risk_list$risk_prob_2)
})

test_that("dynamicPrediction rejects a data_predict_all list of the wrong length", {
  f <- setup_dp_fixture()

  expect_error(
    dynamicPrediction(f$data_predict_all[1], f$long_fit_all, f$survival_fit_all,
                       prediction_time = 5, horizon = 1, time_variable = "year",
                       f$survival_variable_all, f$survival_trans_function,
                       bandcount1 = 10, bandcount2 = 20),
    "`data_predict_all` has 1 element"
  )
})

test_that("dynamicPrediction rejects a missing time_variable column", {
  f <- setup_dp_fixture()

  expect_error(
    dynamicPrediction(f$data_predict_all, f$long_fit_all, f$survival_fit_all,
                       prediction_time = 5, horizon = 1, time_variable = "not_a_column",
                       f$survival_variable_all, f$survival_trans_function,
                       bandcount1 = 10, bandcount2 = 20),
    "not_a_column"
  )
})

test_that("dynamicPrediction rejects a non-positive bandcount1", {
  f <- setup_dp_fixture()

  expect_error(
    dynamicPrediction(f$data_predict_all, f$long_fit_all, f$survival_fit_all,
                       prediction_time = 5, horizon = 1, time_variable = "year",
                       f$survival_variable_all, f$survival_trans_function,
                       bandcount1 = 0, bandcount2 = 20),
    "`bandcount1` must be a single positive number, or the string .auto."
  )
})

test_that("dynamicPrediction rejects mismatched survival_variable_all/survival_trans_function", {
  f <- setup_dp_fixture()

  expect_error(
    dynamicPrediction(f$data_predict_all, f$long_fit_all, f$survival_fit_all,
                       prediction_time = 5, horizon = 1, time_variable = "year",
                       f$survival_variable_all, f$survival_trans_function[1:2],
                       bandcount1 = 10, bandcount2 = 20),
    "survival_variable_all.*survival_trans_function"
  )
})

test_that("assert_survival_trans skips the probe when probe_value = NULL, for backward compatibility", {
  expect_no_error(
    assert_survival_trans(list("Tyears1"), list(function(x) stop("boom")), probe_value = NULL)
  )
})

test_that("assert_survival_trans accepts a well-formed transform when probed", {
  expect_no_error(
    assert_survival_trans(list("Tyears1"), list(function(x) abs(x - 1)), probe_value = 5)
  )
})

test_that("dynamicPrediction rejects a survival_trans_function that throws an error, instead of failing deep inside the per-patient prediction grid", {
  f <- setup_dp_fixture()
  bad_trans <- f$survival_trans_function
  bad_trans[[2]] <- function(x) stop("boom")

  expect_error(
    dynamicPrediction(f$data_predict_all, f$long_fit_all, f$survival_fit_all,
                       prediction_time = 5, horizon = 1, time_variable = "year",
                       f$survival_variable_all, bad_trans,
                       bandcount1 = 10, bandcount2 = 20),
    "survival_trans_function\\[\\[2\\]\\].*failed"
  )
})

test_that("dynamicPrediction rejects a survival_trans_function that returns a character value", {
  f <- setup_dp_fixture()
  bad_trans <- f$survival_trans_function
  bad_trans[[1]] <- function(x) "not a number"

  expect_error(
    dynamicPrediction(f$data_predict_all, f$long_fit_all, f$survival_fit_all,
                       prediction_time = 5, horizon = 1, time_variable = "year",
                       f$survival_variable_all, bad_trans,
                       bandcount1 = 10, bandcount2 = 20),
    "survival_trans_function\\[\\[1\\]\\].*single, finite numeric value"
  )
})

test_that("dynamicPrediction rejects a survival_trans_function that returns a length > 1 vector", {
  f <- setup_dp_fixture()
  bad_trans <- f$survival_trans_function
  bad_trans[[3]] <- function(x) c(x, x + 1)

  expect_error(
    dynamicPrediction(f$data_predict_all, f$long_fit_all, f$survival_fit_all,
                       prediction_time = 5, horizon = 1, time_variable = "year",
                       f$survival_variable_all, bad_trans,
                       bandcount1 = 10, bandcount2 = 20),
    "survival_trans_function\\[\\[3\\]\\].*single, finite numeric value"
  )
})

test_that("dynamicPrediction rejects a survival_trans_function that returns a non-finite value", {
  f <- setup_dp_fixture()
  bad_trans <- f$survival_trans_function
  # log() of a negative number is NaN, not an error -- would otherwise
  # silently propagate into model.matrix() at prediction time
  bad_trans[[4]] <- function(x) log(x - 10)

  expect_error(
    suppressWarnings(
      dynamicPrediction(f$data_predict_all, f$long_fit_all, f$survival_fit_all,
                         prediction_time = 5, horizon = 1, time_variable = "year",
                         f$survival_variable_all, bad_trans,
                         bandcount1 = 10, bandcount2 = 20)
    ),
    "survival_trans_function\\[\\[4\\]\\].*single, finite numeric value"
  )
})

test_that("dynamicPredictionBio rejects an out-of-range bio_i", {
  f <- setup_dp_fixture()

  expect_error(
    dynamicPredictionBio(bio_i = 5, f$data_predict_all, f$long_fit_all, f$survival_fit_all,
                          prediction_time = 5, horizon = 1, time_variable = "year",
                          f$survival_variable_all, f$survival_trans_function,
                          bandcount2 = 20, bandcount3 = 50),
    "`bio_i` must be between 1 and 2"
  )
})

test_that("dynamicPredictionBio rejects a non-integer bio_i", {
  f <- setup_dp_fixture()

  expect_error(
    dynamicPredictionBio(bio_i = 1.5, f$data_predict_all, f$long_fit_all, f$survival_fit_all,
                          prediction_time = 5, horizon = 1, time_variable = "year",
                          f$survival_variable_all, f$survival_trans_function,
                          bandcount2 = 20, bandcount3 = 50),
    "`bio_i` must be a single integer"
  )
})
