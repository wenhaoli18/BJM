# riskPlot()/predictPlot() objects have to survive ggplot_build(): several
# failures only showed when the plot was drawn, not when it was built.

render_ok <- function(p) {
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
}

test_that("riskPlot() draws a patient whose event time is unknown", {
  skip_on_cran()
  f <- setup_dp_fixture()
  # a new patient: survival time (and event type) not known yet. The line
  # marking the event time used to get NA and fail to draw.
  d <- lapply(f$data_predict_all, function(x) { x$years <- NA; x$status4 <- NA; x })
  p <- riskPlot(d, f$long_fit_all, f$survival_fit_all, prediction_time = c(3, 4, 5), bio_i = 1,
                horizon = 1, time_variable = "year", f$survival_variable_all,
                f$survival_trans_function, bandcount1 = 10, bandcount2 = 20)
  render_ok(p)
})

test_that("riskPlot() skips landmark times after the patient's event", {
  skip_on_cran()
  f <- setup_dp_fixture()
  # patient 1 died at 1.1 years; landmarks after that have no prediction
  d <- rep(list(pbc3[pbc3$id == 1, ]), 2)
  p <- riskPlot(d, f$long_fit_all, f$survival_fit_all, prediction_time = c(0.5, 1, 2), bio_i = 1,
                horizon = 0.5, time_variable = "year", f$survival_variable_all,
                f$survival_trans_function, bandcount1 = 10, bandcount2 = 20)
  render_ok(p)
})

test_that("plots keep a positive risk axis when every biomarker value is <= 0", {
  skip_on_cran()
  f <- setup_dp_fixture()
  # serBilir is log-scale; shift this patient's history below 0. The risk
  # axis scale used to be 2 * max(biomarker), here negative.
  d <- lapply(f$data_predict_all, function(x) { x$serBilir <- x$serBilir - 5; x$years <- NA; x })
  p <- predictPlot(d, f$long_fit_all, f$survival_fit_all, prediction_time = 5, horizon = c(0, 1),
                   time_variable = "year", f$survival_variable_all, f$survival_trans_function,
                   bandcount1 = 10, bandcount2 = 20, bandcount3 = 50)
  render_ok(p)
  k <- which(vapply(p$layers, function(l) any(grepl("probType1", deparse(l$mapping$y))), logical(1)))[1]
  scaled <- ggplot2::ggplot_build(p)$data[[k]]$y
  expect_true(all(scaled >= 0))
  expect_true(scaled[2] > 0)

  p <- riskPlot(d, f$long_fit_all, f$survival_fit_all, prediction_time = c(3, 5), bio_i = 1,
                horizon = 1, time_variable = "year", f$survival_variable_all,
                f$survival_trans_function, bandcount1 = 10, bandcount2 = 20)
  render_ok(p)
})

test_that("predictPlot() and riskPlot() plot an ordinal biomarker's history", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  fx <- setup_copula_predict_fixture()
  survival_fit_all <- survivalSub(pbc3[!duplicated(pbc3$id), ], Surv(years, status3) ~ age + sex, NULL)
  d <- fx$to_cat(pbc3[pbc3$id == 2 & pbc3$year <= 5, ])
  d$years <- NA
  d <- list(d, d)
  set.seed(1)
  p <- predictPlot(d, fx$long_fit_all, survival_fit_all, prediction_time = 5, horizon = c(0, 1),
                   time_variable = "year", NULL, NULL, bandcount1 = 10, bandcount2 = 20,
                   bio_his = 2, bio_pred = 2)
  render_ok(p)
  history <- Filter(function(l) "longitudinal" %in% names(l$data), p$layers)[[1]]$data$longitudinal
  expect_true(all(history %in% 1:3))

  p <- riskPlot(d, fx$long_fit_all, survival_fit_all, prediction_time = c(4, 5), bio_i = 2,
                horizon = 1, time_variable = "year", NULL, NULL, bandcount1 = 10, bandcount2 = 20)
  render_ok(p)
})
