# setup_dp_fixture() is defined in helper-fixtures.R (auto-sourced by testthat
# and shared across all test-*.R files).

test_that("predictRisk returns risk probabilities in [0, 1] with competing risks", {
  f <- setup_dp_fixture()

  risk <- predictRisk(f$data_predict_all, f$long_fit_all, f$survival_fit_all,
                             prediction_time = 5, horizon = 1, time_variable = "year",
                             f$survival_variable_all, f$survival_trans_function,
                             bandcount1 = 10, bandcount2 = 20)

  expect_s3_class(risk, "predictRisk.BJM")
  expect_length(risk, 2)
  expect_true(risk[[1]] >= 0 && risk[[1]] <= 1)
  expect_true(risk[[2]] >= 0 && risk[[2]] <= 1)
  expect_true((risk[[1]] + risk[[2]]) <= 1)
})

test_that("predictRisk risk with a zero-length horizon is exactly 0", {
  f <- setup_dp_fixture()

  risk <- predictRisk(f$data_predict_all, f$long_fit_all, f$survival_fit_all,
                             prediction_time = 5, horizon = 0, time_variable = "year",
                             f$survival_variable_all, f$survival_trans_function,
                             bandcount1 = 10, bandcount2 = 20)

  # one zero per at-risk patient, for each cause (this fixture has competing risks)
  expect_equal(risk[[1]], c("2" = 0))
  expect_equal(risk[[2]], c("2" = 0))
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
  risk <- predictRisk(f$data_predict_all, f$long_fit_all, f$survival_fit_all,
                             prediction_time = 5, horizon = 1, time_variable = "year",
                             f$survival_variable_all, f$survival_trans_function,
                             bandcount1 = 10, bandcount2 = 20)

  expect_output(print(risk), "Dynamic Prediction")
  expect_output(summary(risk), "Summary Statistics")
})

test_that("measurements after prediction_time are dropped with a warning", {
  f <- setup_dp_fixture()
  # The fixture holds patient 2's history up to year 5; add the later visits.
  history <- f$data_predict_all
  full <- rep(list(pbc3[pbc3$id == 2, ]), 2)
  n_after <- sum(full[[1]]$year > 5)
  expect_gt(n_after, 0)

  expect_warning(
    risk_full <- predictRisk(full, f$long_fit_all, f$survival_fit_all,
                        prediction_time = 5, horizon = 1, time_variable = "year",
                        f$survival_variable_all, f$survival_trans_function,
                        bandcount1 = 10, bandcount2 = 20),
    sprintf("%d row\\(s\\) in \\[\\[1\\]\\], %d row\\(s\\) in \\[\\[2\\]\\]", n_after, n_after))
  expect_no_warning(
    truncated <- predictRisk(history, f$long_fit_all, f$survival_fit_all,
                             prediction_time = 5, horizon = 1, time_variable = "year",
                             f$survival_variable_all, f$survival_trans_function,
                             bandcount1 = 10, bandcount2 = 20))
  expect_equal(risk_full, truncated)

  expect_warning(
    bio_full <- predictLongitudinal(bio_i = 1, full, f$long_fit_all, f$survival_fit_all,
                                    prediction_time = 5, horizon = 1, time_variable = "year",
                                    f$survival_variable_all, f$survival_trans_function,
                                    bandcount2 = 20, bandcount3 = 50),
    "after prediction_time")
  bio_truncated <- predictLongitudinal(bio_i = 1, history, f$long_fit_all, f$survival_fit_all,
                                       prediction_time = 5, horizon = 1, time_variable = "year",
                                       f$survival_variable_all, f$survival_trans_function,
                                       bandcount2 = 20, bandcount3 = 50)
  expect_equal(bio_full, bio_truncated)

  # bio_i of length > 1 goes through dynamicPredictionBioAll(); one warning per call.
  expect_warning(
    all_full <- predictLongitudinal(bio_i = NULL, full, f$long_fit_all, f$survival_fit_all,
                                    prediction_time = 5, horizon = 1, time_variable = "year",
                                    f$survival_variable_all, f$survival_trans_function,
                                    bandcount2 = 20, bandcount3 = 50),
    "after prediction_time")
  all_truncated <- predictLongitudinal(bio_i = NULL, history, f$long_fit_all, f$survival_fit_all,
                                       prediction_time = 5, horizon = 1, time_variable = "year",
                                       f$survival_variable_all, f$survival_trans_function,
                                       bandcount2 = 20, bandcount3 = 50)
  expect_equal(all_full, all_truncated)
})

test_that("auto bandcount tuning warns about dropped measurements only once", {
  f <- setup_dp_fixture()
  full <- rep(list(pbc3[pbc3$id == 2, ]), 2)
  warnings_seen <- character()
  withCallingHandlers(
    predictRisk(full, f$long_fit_all, f$survival_fit_all,
                prediction_time = 5, horizon = 1, time_variable = "year",
                f$survival_variable_all, f$survival_trans_function),
    warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  expect_equal(sum(grepl("after prediction_time", warnings_seen)), 1)
})

predict_fixture_risk <- function(f, data, prediction_time = 5) {
  predictRisk(data, f$long_fit_all, f$survival_fit_all,
              prediction_time = prediction_time, horizon = 1, time_variable = "year",
              f$survival_variable_all, f$survival_trans_function,
              bandcount1 = 10, bandcount2 = 20)
}

expect_valid_risk <- function(out) {
  probs <- c(out$risk_prob_1, out$risk_prob_2)
  expect_true(all(is.finite(probs)))
  expect_true(all(probs >= 0 & probs <= 1))
}

test_that("the prediction does not depend on the patient's own (future) event time", {
  # The integration upper limit used to be twice the predicted patient's
  # recorded survival time: information not available at prediction_time,
  # which changed the risk, and for a time just after prediction_time
  # failed with "wrong sign in 'by'" (or inflated the risk several-fold).
  f <- setup_dp_fixture()
  base <- predict_fixture_risk(f, f$data_predict_all)
  expect_valid_risk(base)
  for (t in c(5.2, 8, 30, NA)) {
    d <- lapply(f$data_predict_all, function(x) { x$years <- t; x })
    expect_equal(predict_fixture_risk(f, d), base)
  }
})

test_that("a patient whose event time and status are not yet known can be predicted", {
  f <- setup_dp_fixture()
  base <- predict_fixture_risk(f, f$data_predict_all)

  no_status <- lapply(f$data_predict_all, function(d) { d$status3 <- NA; d$status4 <- NA; d })
  expect_equal(predict_fixture_risk(f, no_status), base)

  unknown <- lapply(no_status, function(d) { d$years <- NA; d })
  expect_equal(predict_fixture_risk(f, unknown), base)

  expect_no_error(
    predictLongitudinal(bio_i = 1, unknown, f$long_fit_all, f$survival_fit_all,
                        prediction_time = 5, horizon = 1, time_variable = "year",
                        f$survival_variable_all, f$survival_trans_function,
                        bandcount2 = 20, bandcount3 = 50))
})

test_that("a missing biomarker value only drops that measurement, not the prediction", {
  f <- setup_dp_fixture()
  with_na <- f$data_predict_all
  with_na[[2]]$albumin[1] <- NA
  without_row <- f$data_predict_all
  without_row[[2]] <- without_row[[2]][-1, ]

  out <- predict_fixture_risk(f, with_na)
  expect_valid_risk(out)
  expect_equal(out, predict_fixture_risk(f, without_row))
})

test_that("a patient with no observation of some biomarker is reported", {
  f <- setup_dp_fixture()
  # alongside a predictable patient: a warning names the unpredictable one
  other <- pbc3[pbc3$id == 6 & pbc3$year <= 5, ]
  d <- lapply(f$data_predict_all, function(x) rbind(x, other))
  d[[2]]$albumin[d[[2]]$id == 2] <- NA
  expect_warning(out <- predict_fixture_risk(f, d), "Patient\\(s\\) 2 have no non-missing measurement")
  expect_true(is.na(out$risk_prob_1[["2"]]))
  expect_true(is.finite(out$risk_prob_1[["6"]]))
  # on its own: an error, not a failure deeper down
  d <- f$data_predict_all
  d[[2]]$albumin <- NA
  expect_error(suppressWarnings(predict_fixture_risk(f, d)), "No patient in data_predict_all has a non-missing")
})

test_that("no patient at risk at prediction_time is an error", {
  f <- setup_dp_fixture()
  d <- lapply(f$data_predict_all, function(x) { x$years <- 1; x })
  expect_error(predict_fixture_risk(f, d), "No patient in data_predict_all is at risk")
  expect_error(predictLongitudinal(1, d, f$long_fit_all, f$survival_fit_all, prediction_time = 5,
                                   horizon = 1, time_variable = "year", f$survival_variable_all,
                                   f$survival_trans_function, bandcount2 = 20, bandcount3 = 50),
               "No patient in data_predict_all is at risk")
})

test_that("the event type may enter long_sub_fixed and is unknown for a new patient", {
  f <- setup_dp_fixture()
  long_fit_all <- longitudinalSub(
    list(pbc3[pbc3$status3 == 1, ], pbc3[pbc3$status3 == 1, ]),
    list(serBilir ~ year + age + sex + years + years:year + status4,
         albumin ~ year + age + sex + years + years:year + status4),
    list(~ year | id, ~ year | id))
  d <- lapply(f$data_predict_all, function(x) { x$status4 <- NA; x$years <- NA; x })
  out <- predictRisk(d, long_fit_all, f$survival_fit_all, prediction_time = 5, horizon = 2,
                     time_variable = "year", f$survival_variable_all, f$survival_trans_function,
                     bandcount1 = 10, bandcount2 = 20)
  expect_valid_risk(out)
  expect_equal(names(out$risk_prob_1), "2")
})

pbc3_risk <- function(data, long_fit_all, survival_fit_all) {
  predictRisk(data, long_fit_all, survival_fit_all, prediction_time = 3, horizon = 3,
              time_variable = "year", list("Tyears1"), list(function(x) abs(x - 1)),
              bandcount1 = 10, bandcount2 = 20)
}

pbc3_long_fit <- function(long_sub_random = list(~ year | id, ~ year | id), data = pbc3) {
  longitudinalSub(data[data$status3 == 1, ],
                  list(serBilir ~ year + age + sex + years + years * year,
                       albumin ~ year + age + sex + years + years * year),
                  long_sub_random)
}

test_that("predictions do not depend on the type of the patient id", {
  surv <- survivalSub(pbc3[!duplicated(pbc3$id), ], Surv(years, status3) ~ age + sex, NULL)
  ref_fit <- pbc3_long_fit()
  ref <- pbc3_risk(pbc3[pbc3$id == 2 & pbc3$year <= 3, ], ref_fit, surv)
  ref_bio <- predictLongitudinal(1, pbc3[pbc3$id == 2 & pbc3$year <= 3, ], ref_fit, surv,
                                 3, 3, "year", list("Tyears1"), list(function(x) abs(x - 1)), 20, 50)

  # ids whose factor codes differ from their labels, as well as numeric ids
  for (make_id in list(function(x) as.numeric(as.character(x)),
                       function(x) factor(as.numeric(as.character(x)) + 1000))) {
    d <- pbc3
    d$id <- make_id(d$id)
    fit <- pbc3_long_fit(data = d)
    patient <- d[as.character(d$id) == as.character(make_id(factor(2))) & d$year <= 3, ]
    expect_equal(pbc3_risk(patient, fit, surv), ref, ignore_attr = TRUE)
    bio <- predictLongitudinal(1, patient, fit, surv, 3, 3, "year",
                               list("Tyears1"), list(function(x) abs(x - 1)), 20, 50)
    expect_equal(bio$Y_predict, ref_bio$Y_predict, ignore_attr = TRUE)
  }
})

test_that("biomarkers with different random-effects structures can be predicted", {
  surv <- survivalSub(pbc3[!duplicated(pbc3$id), ], Surv(years, status3) ~ age + sex, NULL)
  fit <- pbc3_long_fit(list(~ 1 | id, ~ year | id))
  patient <- pbc3[pbc3$id == 2 & pbc3$year <= 3, ]
  expect_valid_risk(pbc3_risk(patient, fit, surv))
  bio <- predictLongitudinal(1, patient, fit, surv, 3, 3, "year",
                             list("Tyears1"), list(function(x) abs(x - 1)), 20, 50)
  expect_true(all(is.finite(bio$Y_predict)))
})

test_that("an event-type model without the survival time gives valid predictions", {
  surv <- survivalSub(pbc3[!duplicated(pbc3$id), ], Surv(years, status3) ~ age + sex,
                      status4 ~ age + sex)
  expect_valid_risk(pbc3_risk(pbc3[pbc3$id == 2 & pbc3$year <= 3, ], pbc3_long_fit(), surv))
})

test_that("a stratified Cox model uses each patient's own stratum", {
  surv <- survivalSub(pbc3[!duplicated(pbc3$id), ], Surv(years, status3) ~ age + strata(sex), NULL)
  fit <- pbc3_long_fit()
  first_rows <- pbc3[!duplicated(pbc3$id) & pbc3$years > 3, ]
  ids <- c(first_rows$id[first_rows$sex == 0][1], first_rows$id[first_rows$sex == 1][1])
  both <- pbc3[pbc3$id %in% ids & pbc3$year <= 3, ]

  out <- pbc3_risk(both, fit, surv)
  expect_valid_risk(out)
  # predicting the two patients together equals predicting each one alone
  # (each uses its own stratum's baseline hazard; upper_bound is shared, so
  # compare the marginal survival densities rather than the risks)
  grid <- seq(3, 10, 1)
  data_both <- list(both, both)
  joint <- marginalT(data_both, fit, surv, grid, 20)
  patients <- unique(as.character(both$id))  # marginalT's column order
  for (k in seq_along(patients)) {
    one <- pbc3[pbc3$id == patients[k] & pbc3$year <= 3, ]
    expect_equal(joint[, k], c(marginalT(list(one, one), fit, surv, grid, 20)))
  }
  # and the two strata really do differ
  expect_false(isTRUE(all.equal(joint[, 1], joint[, 2])))
})

test_that("the denominator grid starts at prediction_time and holds equal mass per interval", {
  f <- setup_dp_fixture()
  d <- f$data_predict_all
  upper <- integration_upper_bound(d, f$long_fit_all, f$survival_fit_all, prediction_time = 5)
  grid <- prepare_infinity_grid(d, f$long_fit_all, f$survival_fit_all, 5, upper, bandcount2 = 20)
  edges <- grid$predict.time.infinity.1
  expect_equal(range(edges), c(5, upper))
  expect_true(all(diff(edges) > 0))
  # one patient here, so every interval carries (almost) the same share of
  # that patient's survival mass beyond prediction_time
  mass <- grid$S_T_all_infinity[, 1]
  expect_lt(max(mass) / min(mass), 1.5)
})

test_that("the numerator grid covers exactly the prediction window", {
  # It used to integrate over horizon + one interval, inflating the risk by
  # roughly 1 / bandcount1 (about 10-18% at bandcount1 = 10).
  f <- setup_dp_fixture()
  coarse <- predictRisk(f$data_predict_all, f$long_fit_all, f$survival_fit_all, 5, 1, "year",
                        f$survival_variable_all, f$survival_trans_function, 10, 40)
  fine <- predictRisk(f$data_predict_all, f$long_fit_all, f$survival_fit_all, 5, 1, "year",
                      f$survival_variable_all, f$survival_trans_function, 160, 40)
  expect_equal(unlist(coarse[1:2]), unlist(fine[1:2]), tolerance = 0.02)
})

test_that("auto bandcount tuning converges on pbc3 without a warning", {
  f <- setup_dp_fixture()
  expect_no_warning(
    auto <- predictRisk(f$data_predict_all, f$long_fit_all, f$survival_fit_all, 5, 1, "year",
                        f$survival_variable_all, f$survival_trans_function))
  fine <- predictRisk(f$data_predict_all, f$long_fit_all, f$survival_fit_all, 5, 1, "year",
                      f$survival_variable_all, f$survival_trans_function, 320, 640)
  expect_equal(unlist(auto[1:2]), unlist(fine[1:2]), tolerance = 0.02)
})

test_that("predictions are named by patient id, for the at-risk patients only", {
  surv <- survivalSub(pbc3[!duplicated(pbc3$id), ], Surv(years, status3) ~ age + sex,
                      status4 ~ years + age + sex)
  fit <- pbc3_long_fit()
  # patient 10's event (years = 0.14) is before prediction_time = 3
  d <- pbc3[pbc3$id %in% c(4, 10, 2) & pbc3$year <= 3, ]
  out <- pbc3_risk(d, fit, surv)
  expect_named(out$risk_prob_1, c("2", "4"))
  expect_named(out$risk_prob_2, c("2", "4"))
  # the integration grid is shared by the patients predicted together, so a
  # patient's value differs slightly between calls at this coarse grid
  expect_equal(unname(out$risk_prob_1["2"]),
               unname(pbc3_risk(d[d$id == 2, ], fit, surv)$risk_prob_1), tolerance = 0.03)

  bio <- predictLongitudinal(1, d, fit, surv, 3, 3, "year", list("Tyears1"),
                             list(function(x) abs(x - 1)), 20, 50)
  expect_named(bio$Y_predict, c("2", "4"))
  expect_equal(colnames(bio$Y_density), c("2", "4"))
})

test_that("horizon <= 0 gives a zero risk for every at-risk patient", {
  surv <- survivalSub(pbc3[!duplicated(pbc3$id), ], Surv(years, status3) ~ age + sex,
                      status4 ~ years + age + sex)
  d <- pbc3[pbc3$id %in% c(2, 4) & pbc3$year <= 3, ]
  out <- predictRisk(d, pbc3_long_fit(), surv, prediction_time = 3, horizon = 0,
                     time_variable = "year", list(), list(), 10, 20)
  expect_equal(out$risk_prob_1, c("2" = 0, "4" = 0))
  expect_equal(out$risk_prob_2, c("2" = 0, "4" = 0))
})

test_that("predictions do not depend on the biomarkers' units (no underflow/overflow)", {
  # The conditional densities scale with det(2 * pi * Sigma), which grows with
  # the biomarkers' units and the number of observations. A "+ 1e-20" in the
  # no-competing-risk denominator then swamped it (risk 0.35 -> 1e-26 after
  # just multiplying the biomarkers by 100), and det() overflowed to Inf
  # (risk 0, or NaN with competing risks) for larger scales.
  d1 <- pbc3[!duplicated(pbc3$id), ]
  surv <- survivalSub(d1, Surv(years, status3) ~ age + sex, NULL)
  surv_cr <- survivalSub(d1, Surv(years, status3) ~ age + sex, status4 ~ years + age + sex)
  predict_in_units <- function(k) {
    d <- pbc3
    d$serBilir <- d$serBilir * k
    d$albumin <- d$albumin * k
    fit <- longitudinalSub(d[d$status3 == 1, ],
                           list(serBilir ~ year + age + sex + years, albumin ~ year + age + sex + years),
                           list(~ year | id, ~ year | id))
    patient <- d[d$id == 15 & d$year <= 8, ]  # 10 visits
    list(risk = predictRisk(patient, fit, surv, 8, 2, "year", list(), list(), 10, 20)$risk_prob_1,
         risk_cr = unlist(predictRisk(patient, fit, surv_cr, 8, 2, "year", list(), list(), 10, 20)[1:2]),
         y = predictLongitudinal(1, patient, fit, surv, 8, 2, "year", list(), list(), 20, 200)$Y_predict / k)
  }
  ref <- predict_in_units(1)
  expect_gt(ref$risk, 0.1)
  for (k in c(100, 1e8)) {
    scaled <- predict_in_units(k)
    expect_equal(scaled$risk, ref$risk, tolerance = 1e-3)
    expect_equal(scaled$risk_cr, ref$risk_cr, tolerance = 1e-3)
    expect_equal(scaled$y, ref$y, tolerance = 1e-3)
  }
})

test_that("predictRisk rejects a negative horizon", {
  f <- setup_dp_fixture()
  expect_error(predictRisk(f$data_predict_all, f$long_fit_all, f$survival_fit_all,
                           prediction_time = 5, horizon = -1, time_variable = "year",
                           f$survival_variable_all, f$survival_trans_function,
                           bandcount1 = 10, bandcount2 = 20),
               "must be >= 0")
})

test_that("the baseline cumulative hazard is read as a right-continuous step function", {
  table <- data.frame(hazard = c(0.1, 0.3, 0.6), time = c(1, 2, 3))
  # before the first time, at a jump, just before and just after a jump
  expect_equal(cumulative_baseline_at(table, c(0.5, 1, 1.99, 2, 2.01)),
               c(0, 0.1, 0.1, 0.3, 0.3))
  # order of the table does not matter
  expect_equal(cumulative_baseline_at(table[3:1, ], 1.99), 0.1)
  # past the last time: continues from the last value with the slope of the
  # least-squares line through the table
  slope <- unname(coef(lm(hazard ~ time, table))[2])
  expect_equal(cumulative_baseline_at(table, 4), 0.6 + slope * 1)
})

test_that("the extrapolated baseline cumulative hazard is continuous and non-decreasing", {
  # convex: the least-squares line lies below the last value at the last time
  convex <- data.frame(hazard = c(0.01, 0.04, 0.09, 0.16, 0.5), time = 1:5)
  # concave: it lies above
  concave <- data.frame(hazard = c(0.4, 0.6, 0.7, 0.75, 0.78), time = 1:5)
  for (table in list(convex, concave)) {
    H <- cumulative_baseline_at(table, c(5, 5 + 1e-8, 5.5, 8))
    expect_equal(H[2], H[1], tolerance = 1e-6)
    expect_true(all(diff(H) >= 0))
  }
  # a decreasing least-squares slope cannot occur for a cumulative hazard,
  # but a single-row table has no slope at all: held constant
  expect_equal(cumulative_baseline_at(data.frame(hazard = 0.2, time = 1), c(1, 3)), c(0.2, 0.2))
  # tail_time() inverts the same extrapolation
  t_hit <- tail_time(convex, lp = 0, s = 2, tail_prob = 0.01)
  H_s <- cumulative_baseline_at(convex, 2)
  expect_equal(exp(-(cumulative_baseline_at(convex, t_hit) - H_s)), 0.01, tolerance = 1e-8)
})
