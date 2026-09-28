# predictLongitudinal() is the single exported entry point for
# biomarker-value prediction: it dispatches to the internal
# dynamicPredictionBio() (single bio_i) or dynamicPredictionBioAll() (bio_i
# NULL or length > 1). These tests confirm the dispatch itself -- that
# predictLongitudinal() produces exactly the same result as calling the
# internal engine directly -- not the prediction math itself, which is
# already covered by test-predictRisk.R/test-dynamicPredictionBioAll.R.

test_that("predictLongitudinal dispatches to dynamicPredictionBio for a single bio_i", {
  f <- setup_dp_fixture()

  set.seed(1)
  direct <- dynamicPredictionBio(bio_i = 1, f$data_predict_all, f$long_fit_all, f$survival_fit_all,
                                  prediction_time = 5, horizon = 1, time_variable = "year",
                                  f$survival_variable_all, f$survival_trans_function,
                                  bandcount2 = 10, bandcount3 = 50)

  set.seed(1)
  via_wrapper <- predictLongitudinal(bio_i = 1, f$data_predict_all, f$long_fit_all, f$survival_fit_all,
                                      prediction_time = 5, horizon = 1, time_variable = "year",
                                      f$survival_variable_all, f$survival_trans_function,
                                      bandcount2 = 10, bandcount3 = 50)

  expect_s3_class(via_wrapper, "dynamicPredictionBio.BJM")
  expect_equal(via_wrapper, direct)
})

test_that("predictLongitudinal dispatches to dynamicPredictionBioAll when bio_i names more than one biomarker", {
  f <- setup_dp_fixture()

  set.seed(1)
  direct <- dynamicPredictionBioAll(bio_i = c(1, 2), f$data_predict_all, f$long_fit_all, f$survival_fit_all,
                                     prediction_time = 5, horizon = 1, time_variable = "year",
                                     f$survival_variable_all, f$survival_trans_function,
                                     bandcount2 = 10, bandcount3 = 50)

  set.seed(1)
  via_wrapper <- predictLongitudinal(bio_i = c(1, 2), f$data_predict_all, f$long_fit_all, f$survival_fit_all,
                                      prediction_time = 5, horizon = 1, time_variable = "year",
                                      f$survival_variable_all, f$survival_trans_function,
                                      bandcount2 = 10, bandcount3 = 50)

  expect_s3_class(via_wrapper, "dynamicPredictionBioAll.BJM")
  expect_equal(via_wrapper, direct)
})

test_that("predictLongitudinal dispatches to dynamicPredictionBioAll when bio_i is left NULL (every biomarker)", {
  f <- setup_dp_fixture()

  set.seed(1)
  direct <- dynamicPredictionBioAll(bio_i = NULL, f$data_predict_all, f$long_fit_all, f$survival_fit_all,
                                     prediction_time = 5, horizon = 1, time_variable = "year",
                                     f$survival_variable_all, f$survival_trans_function,
                                     bandcount2 = 10, bandcount3 = 50)

  set.seed(1)
  via_wrapper <- predictLongitudinal(data_predict_all = f$data_predict_all, long_fit_all = f$long_fit_all,
                                      survival_fit_all = f$survival_fit_all,
                                      prediction_time = 5, horizon = 1, time_variable = "year",
                                      survival_variable_all = f$survival_variable_all,
                                      survival_trans_function = f$survival_trans_function,
                                      bandcount2 = 10, bandcount3 = 50)

  expect_s3_class(via_wrapper, "dynamicPredictionBioAll.BJM")
  expect_length(via_wrapper, 2)
  expect_equal(via_wrapper, direct)
})
