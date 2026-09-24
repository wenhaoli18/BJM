# Characterization (golden-master) test.
#
# Captures the full survivalSub -> longitudinalSub -> dynamicPrediction ->
# dynamicPredictionBio -> riskPlot/predictPlot/cmtPlot pipeline on the
# package's own pbc3 data, and compares it against numeric output captured
# from the package before the 0.2.0 API refactor (see testdata/baseline.rds).
#
# This exists to catch numeric drift introduced while extracting shared
# helpers out of the conditionalYT/conditionalYDT/dynamicPrediction family
# of functions. It intentionally accesses fitted-model internals positionally
# via [[ ]] so that it keeps working across the named-list API change.

test_that("full pipeline output matches pre-refactor baseline", {
  skip_on_cran()

  baseline <- readRDS(test_path("testdata", "baseline.rds"))

  data(pbc3, envir = environment())

  data_survival_fitting <- pbc3[!duplicated(pbc3$id), ]
  form_marginal_surv  <- Surv(years, status3) ~ age + sex
  form_conditional_cr <- status4 ~ years + age + sex
  survival_fit_all <- survivalSub(data_survival_fitting, form_marginal_surv, form_conditional_cr)

  expect_equal(coef(survival_fit_all[[1]]), baseline$survival_fit_all_coef, tolerance = 1e-6)
  expect_equal(coef(survival_fit_all[[3]]), baseline$survival_fit_all_glm_coef, tolerance = 1e-6)

  long_sub_fixed <- list(
    "long1" = serBilir ~ year + age + sex +  (years) + (years) * year,
    "long2" = prothrombin ~ year + age + sex + (years) + (years) * year,
    "long3" = albumin ~ year + age + age * year + sex + (years) + (years) * year
  )
  long_sub_random <- list(
    "long1" =  ~ year| id,
    "long2" =  ~ year| id,
    "long3" =  ~ year| id
  )
  survival_variable_all <- list("Tyears1", "Tyears2", "Tyears3", "Tyears4")
  survival_trans_function <- list(
    fun1 = function(x){abs(x - 1)},
    fun2 = function(x){abs(x - 3)},
    fun3 = function(x){abs(x - 5)},
    fun4 = function(x){abs(x - 7)}
  )

  data_fit_all <- list()
  for (i in seq_len(length(long_sub_fixed))) {
    data_fit_all[[i]] <- pbc3[pbc3$status3 == 1, ]
  }

  long_fit_all <- longitudinalSub(data_fit_all, long_sub_fixed, long_sub_random)

  expect_equal(lapply(long_fit_all[[1]], nlme::fixef), baseline$long_fit_all_beta, tolerance = 1e-6)
  expect_equal(sapply(long_fit_all[[1]], function(u) u$sigma), baseline$long_fit_all_sigma, tolerance = 1e-6)
  expect_equal(long_fit_all[[2]], baseline$long_fit_all_D, tolerance = 1e-6)

  i_PID <- 2
  data.raw.predict <- pbc3[pbc3$id == i_PID, ]
  data_predict_all <- list(data.raw.predict, data.raw.predict, data.raw.predict)

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

  rp <- riskPlot(data_predict_all, long_fit_all, survival_fit_all,
                 prediction_time = 5, bio_i = 1,
                 horizon = 1, time_variable = "year",
                 survival_variable_all, survival_trans_function,
                 bandcount1 = 10, bandcount2 = 20)
  expect_equal(sapply(rp$layers, function(l) class(l$geom)[1]), baseline$riskPlot_layer_classes)

  pp <- predictPlot(data_predict_all, long_fit_all, survival_fit_all,
                     prediction_time = 5, horizon = seq(0.5, 1, 0.5), time_variable = "year",
                     survival_variable_all, survival_trans_function,
                     bandcount1 = 10, bandcount2 = 10, bandcount3 = 50,
                     bio_his = 1, bio_pred = 1, density = 1)
  expect_equal(sapply(pp$layers, function(l) class(l$geom)[1]), baseline$predictPlot_layer_classes)

  cp <- cmtPlot(data_plot_all = pbc3, condi_time2event = 5,
                event_type_variable = NULL, event_type = NULL,
                bio_variable = "serBilir", time_variable = "year",
                survival_variable = "years",
                interval_time = 1/12)
  expect_equal(sapply(cp$layers, function(l) class(l$geom)[1]), baseline$cmtPlot_layer_classes)
})
