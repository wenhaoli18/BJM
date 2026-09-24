# Shared fixture for dynamicPrediction / dynamicPredictionBio tests.
setup_dp_fixture <- function() {
  data(pbc3, envir = environment())

  data_survival_fitting <- pbc3[!duplicated(pbc3$id), ]
  survival_fit_all <- survivalSub(data_survival_fitting,
                                   Surv(years, status3) ~ age + sex,
                                   status4 ~ years + age + sex)

  long_sub_fixed <- list(
    "long1" = serBilir ~ year + age + sex + (years) + (years) * year,
    "long2" = albumin ~ year + age + sex + (years) + (years) * year
  )
  long_sub_random <- list("long1" = ~ year | id, "long2" = ~ year | id)
  data_fit_all <- list(pbc3[pbc3$status3 == 1, ], pbc3[pbc3$status3 == 1, ])
  long_fit_all <- longitudinalSub(data_fit_all, long_sub_fixed, long_sub_random)

  survival_variable_all <- list("Tyears1", "Tyears2", "Tyears3", "Tyears4")
  survival_trans_function <- list(
    fun1 = function(x) abs(x - 1),
    fun2 = function(x) abs(x - 3),
    fun3 = function(x) abs(x - 5),
    fun4 = function(x) abs(x - 7)
  )

  data.raw.predict <- pbc3[pbc3$id == 2, ]
  data_predict_all <- list(data.raw.predict, data.raw.predict)

  list(survival_fit_all = survival_fit_all, long_fit_all = long_fit_all,
       survival_variable_all = survival_variable_all,
       survival_trans_function = survival_trans_function,
       data_predict_all = data_predict_all)
}

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
