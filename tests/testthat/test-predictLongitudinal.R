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

test_that("the predicted mode does not snap to the bandcount3 grid", {
  f <- setup_dp_fixture()
  coarse <- predictLongitudinal(1, f$data_predict_all, f$long_fit_all, f$survival_fit_all,
                                5, 1, "year", f$survival_variable_all, f$survival_trans_function,
                                bandcount2 = 20, bandcount3 = 50)
  fine <- predictLongitudinal(1, f$data_predict_all, f$long_fit_all, f$survival_fit_all,
                              5, 1, "year", f$survival_variable_all, f$survival_trans_function,
                              bandcount2 = 20, bandcount3 = 2000)
  # a 50-point grid step here is about 0.26; the refined mode is far closer
  expect_equal(coarse$Y_predict, fine$Y_predict, tolerance = 1e-3)
})

test_that("density_mode recovers a sharp peak between grid points", {
  Y_all <- seq(-1, 1, length.out = 201)
  dens <- dnorm(Y_all, 0.03, 0.1)  # peak ~4
  expect_equal(density_mode(Y_all, dens), 0.03, tolerance = 1e-6)
})

test_that("auto bandcount tuning converges for predictLongitudinal on pbc3", {
  f <- setup_dp_fixture()
  expect_no_warning(
    predictLongitudinal(NULL, f$data_predict_all, f$long_fit_all, f$survival_fit_all,
                        5, 1, "year", f$survival_variable_all, f$survival_trans_function))
})

test_that("a single observed value does not collapse the candidate grid", {
  f <- setup_dp_fixture()
  # Only the baseline measurement: the observed range is 0, which used to
  # make Y_all a single point and return the observed value as the
  # prediction.
  baseline <- pbc3[pbc3$id == 2 & pbc3$year == 0, ]
  pred <- predictLongitudinal(1, list(baseline, baseline), f$long_fit_all, f$survival_fit_all,
                              0.5, 1, "year", f$survival_variable_all, f$survival_trans_function,
                              bandcount2 = 20, bandcount3 = 50)
  expect_gt(length(pred$Y_all), 40)
  expect_true(is.finite(pred$Y_predict))
  expect_false(isTRUE(all.equal(unname(pred$Y_predict), baseline$serBilir)))
  # the predictive density is negligible at the edges of the grid
  expect_lt(max(pred$Y_density[c(1, nrow(pred$Y_density)), ]), 1e-3 * max(pred$Y_density))

  fine <- predictLongitudinal(1, list(baseline, baseline), f$long_fit_all, f$survival_fit_all,
                              0.5, 1, "year", f$survival_variable_all, f$survival_trans_function,
                              bandcount2 = 20, bandcount3 = 1000)
  expect_equal(pred$Y_predict, fine$Y_predict, tolerance = 5e-3)
})
