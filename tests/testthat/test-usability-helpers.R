# Exercises the usability improvements added on top of the core pipeline:
# survivalTrans() (avoids hand-writing survival_variable_all/
# survival_trans_function), bare-data.frame auto-repeat for the various
# data_*_all arguments, and the cmtPlot() id_variable/condi_time2event fixes.

test_that("survivalTrans builds the same survival_variable_all/survival_trans_function pair as hand-writing it", {
  trans <- survivalTrans(c(1, 3, 5, 7))

  expect_equal(trans$survival_variable_all, list("Tyears1", "Tyears2", "Tyears3", "Tyears4"))
  expect_equal(trans$survival_trans_function[[1]](2), abs(2 - 1))
  expect_equal(trans$survival_trans_function[[3]](10), abs(10 - 5))
})

test_that("survivalTrans respects a custom prefix", {
  trans <- survivalTrans(c(2, 4), prefix = "basis")
  expect_equal(trans$survival_variable_all, list("basis1", "basis2"))
})

test_that("survivalTrans rejects empty or non-numeric cut_points", {
  expect_error(survivalTrans(numeric(0)), "`cut_points` must be a non-empty numeric vector")
  expect_error(survivalTrans(c("a", "b")), "`cut_points` must be a non-empty numeric vector")
  expect_error(survivalTrans(c(1, NA)), "`cut_points` must be a non-empty numeric vector")
})

test_that("dynamicPredictionBio accepts a bare data.frame for data_predict_all", {
  f <- setup_dp_fixture()

  bio_bare <- dynamicPredictionBio(bio_i = 1, f$data_predict_all[[1]], f$long_fit_all, f$survival_fit_all,
                                    prediction_time = 5, horizon = 1, time_variable = "year",
                                    f$survival_variable_all, f$survival_trans_function,
                                    bandcount2 = 20, bandcount3 = 50)
  bio_list <- dynamicPredictionBio(bio_i = 1, f$data_predict_all, f$long_fit_all, f$survival_fit_all,
                                    prediction_time = 5, horizon = 1, time_variable = "year",
                                    f$survival_variable_all, f$survival_trans_function,
                                    bandcount2 = 20, bandcount3 = 50)

  expect_equal(bio_bare$Y_predict, bio_list$Y_predict)
})

test_that("cmtPlot picks the midpoint of time_variable when condi_time2event is NULL", {
  data(pbc3, envir = environment())

  cp <- cmtPlot(data_plot_all = pbc3, condi_time2event = NULL,
                event_type_variable = NULL, event_type = NULL,
                bio_variable = "serBilir", time_variable = "year",
                survival_variable = "years", interval_time = 1 / 12)

  expect_s3_class(cp, "ggplot")
})

test_that("cmtPlot deduplicates by id_variable, not a literal 'id_variable' column", {
  data(pbc3, envir = environment())
  pbc3_renamed <- pbc3
  names(pbc3_renamed)[names(pbc3_renamed) == "id"] <- "patient_id"

  cp <- cmtPlot(data_plot_all = pbc3_renamed, condi_time2event = 5,
                event_type_variable = NULL, event_type = NULL,
                bio_variable = "serBilir", time_variable = "year",
                survival_variable = "years", interval_time = 1 / 12,
                id_variable = "patient_id")

  expect_s3_class(cp, "ggplot")
})

test_that("cmtPlot rejects an id_variable not present in data_plot_all", {
  data(pbc3, envir = environment())

  expect_error(
    cmtPlot(data_plot_all = pbc3, condi_time2event = 5,
            event_type_variable = NULL, event_type = NULL,
            bio_variable = "serBilir", time_variable = "year",
            survival_variable = "years", interval_time = 1 / 12,
            id_variable = "not_a_column"),
    "not_a_column"
  )
})
