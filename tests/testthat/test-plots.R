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

test_that("spaghettiPlot() draws baseline- and event-aligned trajectories", {
  p <- spaghettiPlot(pbc3, "serBilir", "year", event_type_variable = "status2")
  render_ok(p)
  expect_equal(nrow(p$data), sum(!is.na(pbc3$serBilir)))

  p <- spaghettiPlot(pbc3, "serBilir", "year", survival_variable = "years",
                     event_type_variable = "status2", align = "event", smooth = FALSE)
  render_ok(p)
  expect_true(all(p$data$x >= 0))

  set.seed(1)
  p <- spaghettiPlot(pbc3, "albumin", "year", n_subjects = 10)
  render_ok(p)
  expect_equal(length(unique(p$data$id)), 10)
})

test_that("spaghettiPlot() validates its inputs", {
  expect_error(spaghettiPlot(pbc3, "serBilir", "year", align = "event"), "survival_variable")
  expect_error(spaghettiPlot(pbc3, "serBilir", "year", event_type_variable = "ascites"),
               "constant within each subject")
  expect_error(spaghettiPlot(pbc3, "not_a_column", "year"), "not_a_column")
})

test_that("cifPlot() matches 1 - Kaplan-Meier for a single event type", {
  p <- cifPlot(pbc3, "years", "status2", conf_int = FALSE)
  render_ok(p)
  d <- pbc3[!duplicated(pbc3$id), ]
  km <- survival::survfit(survival::Surv(years, status2) ~ 1, data = d)
  expect_equal(p$data$cif[-1], 1 - km$surv)
})

test_that("cifPlot() handles competing risks coded with 0 or NA as censoring", {
  p5 <- cifPlot(pbc3, "years", "status5", event_labels = c("1" = "Death", "2" = "Transplant"))
  p4 <- cifPlot(pbc3, "years", "status4", censor_value = NA,
                event_labels = c("0" = "Death", "1" = "Transplant"))
  render_ok(p5)
  render_ok(p4)
  expect_equal(p5$data$cif, p4$data$cif)
  expect_equal(levels(p5$data$event), c("Death", "Transplant"))
  # the two incidences plus survival sum to at most 1 at every time
  last <- tapply(p5$data$cif, p5$data$event, function(v) v[length(v)])
  expect_true(sum(last) <= 1)
  # the confidence band is defined from time 0, including before the first transplant
  band <- ggplot2::ggplot_build(p5)$data[[1]]
  expect_false(anyNA(band$ymin))
  expect_equal(min(band$x), 0)
})

test_that("cifPlot() stratifies by group and validates its inputs", {
  p <- cifPlot(pbc3, "years", "status2", group_variable = "drug")
  render_ok(p)
  expect_setequal(levels(p$data$group), levels(pbc3$drug))

  expect_warning(cifPlot(pbc3, "years", "status4"), "censor_value = NA")
  expect_error(cifPlot(pbc3, "years", "status5", event_labels = c("1" = "Death")), "no label")
  expect_error(cifPlot(pbc3, "years", "status2", censor_value = 1:2), "single value")
  expect_error(cifPlot(pbc3, "years", "ascites"), "constant within each subject")
})

test_that("plot.longitudinalSub.BJM() draws every diagnostic", {
  skip_on_cran()
  f <- setup_dp_fixture()
  for (w in c("residuals", "qq", "ranef", "corr")) render_ok(plot(f$long_fit_all, which = w))
  corr <- plot(f$long_fit_all, which = "corr")$data
  expect_equal(nrow(corr), 16)
  expect_true(all(abs(corr$value) <= 1 + 1e-8))
  expect_true("serBilir: (Intercept)" %in% levels(corr$col))
  ranef <- plot(f$long_fit_all, which = "ranef")$data
  expect_setequal(levels(ranef$panel), c("serBilir: (Intercept)", "serBilir: year",
                                         "albumin: (Intercept)", "albumin: year"))
})

test_that("plot.longitudinalSub.BJM() leaves ordinal biomarkers out of residual plots", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  fx <- setup_copula_predict_fixture()
  expect_message(p <- plot(fx$long_fit_all), "albumin_cat")
  render_ok(p)
  expect_equal(levels(p$data$biomarker), "serBilir")
  render_ok(plot(fx$long_fit_all, which = "ranef"))
  render_ok(plot(fx$long_fit_all, which = "corr"))
})

test_that("plot.survivalSub.BJM() draws hazard/odds ratios and the baseline hazard", {
  d <- pbc3[!duplicated(pbc3$id), ]
  s_cr <- survivalSub(d, Surv(years, status3) ~ age + sex, status4 ~ years + age + sex)
  p <- plot(s_cr)
  render_ok(p)
  expect_equal(nlevels(p$data$model), 2)
  hr <- p$data[p$data$model == "Survival model: hazard ratio", ]
  expect_equal(unname(hr$estimate), unname(exp(coef(s_cr$coxph_fit))))

  s <- survivalSub(d, Surv(years, status3) ~ age + sex, NULL)
  expect_equal(nlevels(plot(s)$data$model), 1)
  p <- plot(s, which = "basehaz")
  render_ok(p)
  expect_equal(p$data$hazard[1], 0)
  render_ok(plot(survivalSub(d, Surv(years, status3) ~ age + strata(sex), NULL), which = "basehaz"))
})

test_that("plot.predictLongitudinal.BJM() draws densities and selects patients", {
  skip_on_cran()
  f <- setup_dp_fixture()
  d <- pbc3[pbc3$id %in% c(2, 4, 5) & pbc3$year <= 3, ]
  d$years <- NA
  pred <- predictLongitudinal(list(d, d), f$long_fit_all, f$survival_fit_all, prediction_time = 3,
                              horizon = 1, time_variable = "year", f$survival_variable_all,
                              f$survival_trans_function, bio_i = 1, bandcount2 = 20, bandcount3 = 40)
  p <- plot(pred)
  render_ok(p)
  expect_equal(levels(p$data$subject), names(pred$Y_predict))
  expect_equal(levels(plot(pred, subject = "4")$data$subject), "4")
  expect_equal(levels(plot(pred, subject = 2:3)$data$subject), names(pred$Y_predict)[2:3])
  expect_error(plot(pred, subject = "999"), "999")
  expect_error(plot(pred, subject = 10), "between 1 and 3")
})

test_that("plot.predictLongitudinal.BJM() draws category probabilities for an ordinal biomarker", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  fx <- setup_copula_predict_fixture()
  s <- survivalSub(pbc3[!duplicated(pbc3$id), ], Surv(years, status3) ~ age + sex, NULL)
  d <- fx$to_cat(pbc3[pbc3$id %in% c(2, 4) & pbc3$year <= 5, ])
  d$years <- NA
  set.seed(1)
  pred <- predictLongitudinal(list(d, d), fx$long_fit_all, s, prediction_time = 5, horizon = 1,
                              time_variable = "year", NULL, NULL, bio_i = 2, bandcount2 = 20)
  p <- plot(pred)
  render_ok(p)
  expect_equal(levels(p$data$category), c("low", "mid", "high"))
})

test_that("performance_metrics() reduces to Mann-Whitney AUC and plain Brier without censoring", {
  set.seed(1)
  n <- 60
  time <- runif(n, 1, 5)
  status <- rep(1, n)
  risk <- runif(n)
  m <- performance_metrics(risk, time, status, rep(NA, n), 1, s = 0.5, horizon = 2)
  case <- time <= 2.5
  expect_equal(m$auc, mean(outer(risk[case], risk[!case], ">")))
  expect_equal(m$brier, mean((case - risk)^2))
  expect_equal(m$n_cases, sum(case))
})

test_that("performance_metrics() matches timeROC and riskRegression under censoring", {
  skip_if_not_installed("timeROC")
  skip_if_not_installed("riskRegression")
  skip_if_not_installed("prodlim")
  set.seed(42)
  n <- 300
  x <- rnorm(n)
  t1 <- rexp(n, 0.3 * exp(0.8 * x)); t2 <- rexp(n, 0.15); cc <- runif(n, 0, 6)
  time <- pmin(t1, t2, cc)
  event <- ifelse(cc < pmin(t1, t2), 0, ifelse(t1 < t2, 1, 2))
  risk <- plogis(0.8 * x - 0.5)
  d <- data.frame(time = time, event = event, status = as.numeric(event > 0))
  Hist <- prodlim::Hist  # Score() needs Hist() visible from the formula

  m <- performance_metrics(risk, time, d$status, rep(NA, n), 1, 0, 2)
  expect_equal(m$auc, unname(timeROC::timeROC(time, d$status, risk, cause = 1, times = 2)$AUC[2]))
  sc <- riskRegression::Score(list(m = matrix(risk, ncol = 1)), Hist(time, status) ~ 1,
                              data = d, times = 2, metrics = "brier", null.model = FALSE)
  expect_equal(m$brier, sc$Brier$score$Brier)

  cause <- ifelse(event == 0, NA, event)
  m <- performance_metrics(risk, time, d$status, cause, 1, 0, 2)
  expect_equal(m$auc, unname(timeROC::timeROC(time, event, risk, cause = 1, times = 2)$AUC_1[2]))
  sc <- riskRegression::Score(list(m = matrix(risk, ncol = 1)), Hist(time, event) ~ 1,
                              data = d, times = 2, cause = 1, metrics = "brier", null.model = FALSE)
  expect_equal(m$brier, sc$Brier$score$Brier)
})

test_that("performancePlot() evaluates each landmark and event type", {
  skip_on_cran()
  f <- setup_dp_fixture()
  ids <- unique(pbc3$id)[1:80]
  d <- pbc3[pbc3$id %in% ids, ]
  p <- suppressWarnings(performancePlot(d, f$long_fit_all, f$survival_fit_all, prediction_time = c(2, 4),
                                        horizon = 2, time_variable = "year", f$survival_variable_all,
                                        f$survival_trans_function, bandcount1 = 10, bandcount2 = 40))
  render_ok(p)
  perf <- p$data
  expect_equal(nrow(perf), 2 * 2 * 2)
  expect_equal(nlevels(perf$cause), 2)
  # the subjects predicted at each landmark are exactly those event-free then
  first <- d[!duplicated(d$id), ]
  expect_equal(perf$n_at_risk[perf$landmark == 2][1], sum(first$years > 2))
  expect_true(all(perf$value[perf$measure == "AUC"] >= 0 & perf$value[perf$measure == "AUC"] <= 1,
                  na.rm = TRUE))

  expect_error(performancePlot(d, f$long_fit_all, f$survival_fit_all, prediction_time = 2, horizon = 0,
                               time_variable = "year", f$survival_variable_all, f$survival_trans_function),
               "positive")
  d$years <- NULL
  expect_error(performancePlot(d, f$long_fit_all, f$survival_fit_all, prediction_time = 2, horizon = 2,
                               time_variable = "year", f$survival_variable_all, f$survival_trans_function),
               "years")
})

test_that("a subject without a prediction is left out of performancePlot(), not NA for everyone", {
  skip_on_cran()
  f <- setup_dp_fixture()
  ids <- unique(pbc3$id)[1:60]
  d <- pbc3[pbc3$id %in% ids, ]
  d$albumin[d$id == ids[3]] <- NA  # never has albumin: predictRisk() gives NA
  expect_message(
    p <- suppressWarnings(performancePlot(d, f$long_fit_all, f$survival_fit_all, prediction_time = 2,
                                          horizon = 2, time_variable = "year", f$survival_variable_all,
                                          f$survival_trans_function, bandcount1 = 10, bandcount2 = 20)),
    "1 at-risk subject\\(s\\) without a prediction")
  # Brier is defined whenever anyone is evaluated (AUC can still be NA for a
  # cause with no cases in the window)
  expect_false(anyNA(p$data$value[p$data$measure == "Brier score"]))
  first <- d[!duplicated(d$id), ]
  expect_equal(p$data$n_at_risk[1], sum(first$years > 2) - 1)
})

test_that("calibration_groups() gives 1 - Kaplan-Meier and Aalen-Johansen per risk group", {
  set.seed(3)
  n <- 200
  x <- rnorm(n)
  t1 <- rexp(n, 0.3 * exp(0.8 * x)); t2 <- rexp(n, 0.15); cc <- runif(n, 0, 6)
  time <- pmin(t1, t2, cc)
  event <- ifelse(cc < pmin(t1, t2), 0, ifelse(t1 < t2, 1, 2))
  status <- as.numeric(event > 0)
  risk <- plogis(0.8 * x - 1)
  grp <- ceiling(rank(risk, ties.method = "first") * 4 / n)

  cg <- calibration_groups(risk, time, status, rep(NA, n), 1, end = 2, n_groups = 4)
  km <- vapply(1:4, function(g) {
    1 - summary(survival::survfit(survival::Surv(time[grp == g], status[grp == g]) ~ 1), times = 2)$surv
  }, numeric(1))
  expect_equal(cg$observed, km)
  expect_equal(cg$predicted, as.vector(tapply(risk, grp, mean)))
  expect_equal(sum(cg$n), n)
  expect_true(all(cg$lower <= cg$observed & cg$observed <= cg$upper))

  cause <- ifelse(event == 0, NA, event)
  cg <- calibration_groups(risk, time, status, cause, 2, end = 2, n_groups = 4)
  aj <- vapply(1:4, function(g) {
    f <- survival::survfit(survival::Surv(time[grp == g], factor(event[grp == g], 0:2)) ~ 1)
    unname(summary(f, times = 2)$pstate[1, match("2", f$states)])
  }, numeric(1))
  expect_equal(cg$observed, aj)
})

test_that("calibrationPlot() draws one panel per landmark and event type", {
  skip_on_cran()
  f <- setup_dp_fixture()
  ids <- unique(pbc3$id)[1:80]
  d <- pbc3[pbc3$id %in% ids, ]
  p <- suppressWarnings(calibrationPlot(d, f$long_fit_all, f$survival_fit_all, prediction_time = c(2, 4),
                                        horizon = 2, time_variable = "year", f$survival_variable_all,
                                        f$survival_trans_function, n_groups = 4,
                                        bandcount1 = 10, bandcount2 = 40))
  render_ok(p)
  cal <- p$data
  expect_equal(nrow(cal), 2 * 2 * 4)
  expect_equal(nlevels(cal$cause), 2)
  first <- d[!duplicated(d$id), ]
  expect_equal(sum(cal$n[cal$landmark == 2 & cal$cause == levels(cal$cause)[1]]), sum(first$years > 2))
  expect_true(all(cal$predicted >= 0 & cal$predicted <= 1))
  expect_error(calibrationPlot(d, f$long_fit_all, f$survival_fit_all, prediction_time = 2, horizon = 2,
                               time_variable = "year", f$survival_variable_all,
                               f$survival_trans_function, n_groups = 0), "n_groups")
})

test_that("the risk axis is labelled only within [0, 1]", {
  # a biomarker below 0 used to show negative "probabilities" on the risk axis
  d <- data.frame(time = 1:3, y = c(-0.5, 0.2, 1))
  p <- ggplot2::ggplot(d, ggplot2::aes(time, y)) + ggplot2::geom_point() + risk_axis_scale(2)
  labels <- as.numeric(ggplot2::ggplot_build(p)$layout$panel_params[[1]]$y.sec$get_labels())
  expect_true(all(labels >= 0 & labels <= 1))
  # 1 is kept although pretty() can return it as 1 + 2e-10
  d$y <- c(0, 1, 2)
  p <- ggplot2::ggplot(d, ggplot2::aes(time, y)) + ggplot2::geom_point() + risk_axis_scale(2)
  expect_true(1 %in% as.numeric(ggplot2::ggplot_build(p)$layout$panel_params[[1]]$y.sec$get_labels()))
  # an ordinal biomarker's axis shows its categories
  p <- ggplot2::ggplot(d, ggplot2::aes(time, y)) + ggplot2::geom_point() + risk_axis_scale(2, c("low", "mid", "high"))
  expect_equal(ggplot2::ggplot_build(p)$layout$panel_params[[1]]$y$get_labels(), c("low", "mid", "high"))
})

test_that("cmtPlot() draws without messages, and the baseline plot names its reference", {
  data(pbc3, envir = environment())
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_message(print(cmtPlot(pbc3, 5, NULL, NULL, "serBilir", "year", "years", 1 / 12,
                                       status_variable = "status3")))
  sf <- survivalSub(pbc3[!duplicated(pbc3$id), ], Surv(years, status3) ~ age + sex, NULL)
  expect_match(plot(sf, which = "basehaz")$labels$y, "all covariates 0")
})

test_that("predictPlot() draws an ordinal biomarker's category probabilities", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  data(pbc3, envir = environment())
  q <- stats::quantile(pbc3$albumin[pbc3$status3 == 1], c(0, 1 / 3, 2 / 3, 1), na.rm = TRUE)
  cut3 <- function(x) cut(x, q, include.lowest = TRUE, labels = c("low", "mid", "high"), ordered_result = TRUE)
  d1 <- pbc3[pbc3$status3 == 1, ]
  d1$alb <- cut3(d1$albumin)
  fit <- longitudinalSub(d1, list(serBilir ~ year + age + sex + years, alb ~ year + age + sex + years),
                         list(~ year | id, ~ year | id))
  sf <- survivalSub(pbc3[!duplicated(pbc3$id), ], Surv(years, status3) ~ age + sex, NULL)
  p2 <- pbc3[pbc3$id == 2, ]
  p2$alb <- cut3(p2$albumin)
  p <- predictPlot(p2, fit, sf, prediction_time = 5, bio_his = 2, bio_pred = 2, horizon = c(1, 2),
                   time_variable = "year", list(), list(), bandcount1 = 10, bandcount2 = 20)
  probs <- Filter(function(l) "prob" %in% names(l$data), p$layers)
  expect_length(probs, 1)
  sums <- tapply(probs[[1]]$data$prob, probs[[1]]$data$time, sum)
  expect_equal(as.vector(sums), c(1, 1), tolerance = 1e-8)
  # no quantile bands for a handful of categories
  expect_false(any(vapply(p$layers, function(l) inherits(l$geom, "GeomRibbon"), logical(1))))
  render_ok(p)
})

test_that("plot.predictLongitudinal.BJM() leaves out the empty tails of the value grid", {
  skip_on_cran()
  f <- setup_dp_fixture()
  pl <- predictLongitudinal(1, f$data_predict_all, f$long_fit_all, f$survival_fit_all, 5, 1, "year",
                            f$survival_variable_all, f$survival_trans_function, 20, 200)
  shown <- range(plot(pl)$data$value)
  expect_gt(shown[1], min(pl$Y_all))
  expect_lt(shown[2], max(pl$Y_all))
})

test_that("riskPlot()'s end-of-follow-up line tells an event from a censoring", {
  data(pbc3, envir = environment())
  first <- pbc3[!duplicated(pbc3$id), ]
  sf <- survivalSub(first, Surv(years, status3) ~ age + sex, NULL)
  sfc <- survivalSub(first, Surv(years, status3) ~ age + sex, status4 ~ years + age + sex)
  line_of <- function(id, fit, type = NULL) {
    layers <- event_time_vline(list(pbc3[pbc3$id == id, ]), fit, type)
    if (length(layers) == 0) return(NULL)
    c(label = layers[[1]]$data$end, colour = layers[[1]]$aes_params$colour,
      linetype = unname(layers[[2]]$palette(1)))
  }
  # patient 2 is censored: it used to get the same red line as an event
  expect_equal(line_of(2, sf), c(label = "Censored", colour = "grey40", linetype = "dashed"))
  expect_equal(line_of(4, sf), c(label = "Event", colour = "red", linetype = "solid"))
  expect_equal(line_of(4, sfc, "status4")[["label"]], "Event type1")
  expect_equal(line_of(5, sfc, "status4")[["label"]], "Event type2")
  # without the status (a new patient) nothing is drawn
  newp <- pbc3[pbc3$id == 4, ]
  newp$status3 <- NULL
  expect_length(event_time_vline(list(newp), sf), 0)
})

test_that("cmtPlot() uses patients with an event, on the data's own time scale", {
  data(pbc3, envir = environment())
  cp <- cmtPlot(pbc3, 5, NULL, NULL, "serBilir", "year", "years", 1 / 12, status_variable = "status3")
  pts <- cp$layers[[1]]$data
  # the means of the 17 patients with an event in (4, 5] (23 with censored
  # ones), per month (t - 1/12, t], plotted at the month's end in years
  ev <- pbc3[pbc3$years > 4 & pbc3$years <= 5 & pbc3$status3 == 1, ]
  expect_equal(length(unique(ev$id)), 17)
  ends <- seq(min(ev$year), 5, by = 1 / 12)
  bin <- findInterval(ev$year, c(-Inf, ends), left.open = TRUE)
  expected <- tapply(ev$serBilir, bin, mean)
  expect_equal(pts$time, ends[as.integer(names(expected))])
  expect_equal(pts$probEvent, as.vector(expected))
  # months without a measurement are left out, not copied from the previous one
  expect_lt(nrow(pts), length(ends))
  expect_lte(max(pts$time), 5)
  expect_equal(cp$labels$x, "year")
  # censored patients enter only, with a warning, when no status is given
  expect_warning(cmtPlot(pbc3, 5, NULL, NULL, "serBilir", "year", "years", 1 / 12), "status_variable")
})

test_that("cmtPlot() gives the same trajectory whatever the time unit", {
  data(pbc3, envir = environment())
  in_months <- pbc3
  in_months$year <- in_months$year * 12
  in_months$years <- in_months$years * 12
  by_year <- cmtPlot(pbc3, 5, NULL, NULL, "serBilir", "year", "years", 1 / 12, status_variable = "status3")
  by_month <- cmtPlot(in_months, 60, NULL, NULL, "serBilir", "year", "years", 1,
                      status_variable = "status3", event_window = 12)
  expect_equal(by_month$layers[[1]]$data$probEvent, by_year$layers[[1]]$data$probEvent)
  expect_equal(by_month$layers[[1]]$data$time, 12 * by_year$layers[[1]]$data$time)
  # competing risks: one trajectory per event type, from those types only
  cr <- cmtPlot(pbc3, 5, "status4", c(0, 1), "albumin", "year", "years", 1 / 4)
  expect_setequal(unique(unlist(lapply(cr$layers, function(l) l$data$line))), c("Event1", "Event2"))
})
