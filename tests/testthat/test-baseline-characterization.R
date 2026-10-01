# Characterization (golden-master) test.
#
# Captures the full survivalSub -> longitudinalSub -> predictRisk ->
# dynamicPredictionBio -> riskPlot/predictPlot/cmtPlot pipeline on the
# package's own pbc3 data, and compares it against numeric output captured
# from the package before the 0.2.0 API refactor (see testdata/baseline.rds).
#
# This exists to catch numeric drift introduced while extracting shared
# helpers out of the conditionalYT/conditionalYDT/predictRisk family
# of functions. It intentionally accesses fitted-model internals positionally
# via [[ ]] so that it keeps working across the named-list API change.
#
# The predictRisk()/dynamicPredictionBio() entries (risk_pred, Y_predict_mode,
# Y_density_summary, Y_all_range) were regenerated in 0.3.0, when those
# functions began dropping measurements after prediction_time: the original
# baseline had been computed from patient 2's full history. The regenerated
# values equal the pre-change code's output on the truncated history.
#
# risk_pred and Y_density_summary were regenerated again later in 0.3.0, when
# the denominator's integration upper limit stopped being twice the predicted
# patient's own (future) survival time and its grid became equal-mass and
# started at prediction_time; at this test's coarse bandcounts that moves
# them by about 2%. Both versions converge to the same value as the grid is
# refined. risk_pred moved again (about -15%) when the numerator's grid was
# fixed to tile exactly (prediction_time, prediction_time + horizon]: the old
# grid integrated over horizon + 1 interval, so at bandcount1 = 10 the stored
# value was about 18% above the fine-grid answer; the new one is within 0.5%.
# Y_predict_mode moved (1.144 -> about 1.21) when the predicted mode stopped
# snapping to the bandcount3 grid: at this test's bandcount3 = 50 (grid step
# 0.26) the snapped mode was 5% off, while the refined one matches a
# 5000-point grid to 4 decimals.
# All three then moved in the 6th significant digit when the baseline hazard
# past the last follow-up time stopped being tabulated on a 0.005 grid and
# was evaluated directly.
# risk_pred then rose about 6% (and Y_predict_mode/Y_density_summary moved
# in the 4th significant digit) when the Breslow cumulative hazard began to
# be read as a step function instead of at the nearest tabulated time: the
# nearest time to prediction_time = 5 is an event at 5.002, whose jump was
# counted before the window instead of inside it. The two readings converge
# to different values, so this is a correction, not grid noise.

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
  # predictRisk()/dynamicPredictionBio() condition on the history up to
  # prediction_time only; riskPlot()/predictPlot() below take the full
  # history and truncate per landmark time themselves.
  history_predict_all <- lapply(data_predict_all, function(d) d[d$year <= 5, ])

  risk_pred <- predictRisk(history_predict_all, long_fit_all, survival_fit_all,
                                  prediction_time = 5, horizon = 1, time_variable = "year",
                                  survival_variable_all, survival_trans_function,
                                  bandcount1 = 10, bandcount2 = 20)
  expect_equal(unname(unlist(risk_pred)), baseline$risk_pred, tolerance = 1e-6)

  Y_pred <- dynamicPredictionBio(bio_i = 1, history_predict_all, long_fit_all, survival_fit_all,
                                  prediction_time = 5, horizon = 1, time_variable = "year",
                                  survival_variable_all, survival_trans_function,
                                  bandcount2 = 20, bandcount3 = 50)
  expect_equal(unname(Y_pred[[1]]), baseline$Y_predict_mode, tolerance = 1e-6)
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
