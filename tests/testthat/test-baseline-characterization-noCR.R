# Characterization (golden-master) test for the no-competing-risk code path.
#
# The main baseline (test-baseline-characterization.R) always fits a
# competing-risk model (form_conditional_cr non-NULL), so it only ever
# exercises conditionalYDT/conditionalYDTBio. conditionalYT/conditionalYTBio
# (the form_conditional_cr = NULL branch) had no coverage at all before this
# test, which would have made an internal-helper extraction across that pair
# unverifiable. Captures numeric output from the package before that
# extraction (see testdata/baseline_noCR.rds).

test_that("no-competing-risk pipeline output matches pre-refactor baseline", {
  skip_on_cran()

  baseline <- readRDS(test_path("testdata", "baseline_noCR.rds"))

  data(pbc3, envir = environment())

  data_survival_fitting <- pbc3[!duplicated(pbc3$id), ]
  survival_fit_all <- survivalSub(data_survival_fitting,
                                   Surv(years, status3) ~ age + sex,
                                   NULL)

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

  risk_pred <- dynamicPrediction(data_predict_all, long_fit_all, survival_fit_all,
                                  prediction_time = 5, horizon = 1, time_variable = "year",
                                  survival_variable_all, survival_trans_function,
                                  bandcount1 = 10, bandcount2 = 20)
  expect_equal(unname(unlist(risk_pred)), baseline$risk_pred, tolerance = 1e-6)

  Y_pred <- dynamicPredictionBio(bio_i = 1, data_predict_all, long_fit_all, survival_fit_all,
                                  prediction_time = 5, horizon = 1, time_variable = "year",
                                  survival_variable_all, survival_trans_function,
                                  bandcount2 = 20, bandcount3 = 50)
  expect_equal(Y_pred[[1]], baseline$Y_predict_mode, tolerance = 1e-6)
  expect_equal(c(sum = sum(Y_pred[[2]]), mean = mean(Y_pred[[2]]), max = max(Y_pred[[2]])),
               baseline$Y_density_summary, tolerance = 1e-6)
  expect_equal(range(Y_pred[[3]]), baseline$Y_all_range, tolerance = 1e-6)
})
