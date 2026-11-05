test_that("baseline transport to the source data reproduces the source fit", {
  fx <- setup_dp_fixture()
  data(pbc3, envir = environment())
  first <- pbc3[!duplicated(pbc3$id), ]
  src <- fx$survival_fit_all
  tr <- transportBJM(src, "baseline", data_target = first)
  expect_s3_class(tr, "transportBJM")
  expect_s3_class(tr, "survivalSub.BJM")
  bh <- survival::basehaz(src$coxph_fit, centered = FALSE)
  expect_equal(tr$transport$cum_basehaz$hazard, bh$hazard, tolerance = 1e-10)
  expect_equal(tr$transport$cum_basehaz$time, bh$time)
  args <- list(fx$data_predict_all, fx$long_fit_all, prediction_time = 5, horizon = 1,
               time_variable = "year", survival_variable_all = fx$survival_variable_all,
               survival_trans_function = fx$survival_trans_function, bandcount1 = 10, bandcount2 = 20)
  r_src <- do.call(predictRisk, c(args[1:2], list(src), args[-(1:2)]))
  r_tr <- do.call(predictRisk, c(args[1:2], list(tr), args[-(1:2)]))
  expect_equal(r_tr$risk_prob_1, r_src$risk_prob_1, tolerance = 1e-10)
  expect_equal(r_tr$risk_prob_2, r_src$risk_prob_2, tolerance = 1e-10)
})

test_that("recalibration on the source data has slope 1", {
  data(pbc3, envir = environment())
  first <- pbc3[!duplicated(pbc3$id), ]
  src <- survivalSub(first, Surv(years, status3) ~ age + sex, NULL)
  tr <- transportBJM(src, "recalibrate", data_target = first)
  expect_equal(tr$transport$lp_slope, 1, tolerance = 1e-4)
})

test_that("a higher target baseline hazard raises the predicted risk", {
  data(pbc3, envir = environment())
  first <- pbc3[!duplicated(pbc3$id), ]
  src <- survivalSub(first, Surv(years, status3) ~ age + sex, NULL)
  target <- first
  target$years <- target$years / 2   # events twice as fast
  tr <- transportBJM(src, "baseline", data_target = target)
  expect_gt(survival_cum_basehaz(tr)$hazard[nrow(tr$transport$cum_basehaz)] /
              max(survival_cum_basehaz(src)$hazard[survival_cum_basehaz(src)$time <= max(target$years)]), 1)
  expect_equal(survival_max_followup(tr), max(target$years))
})

test_that("unlabeled transport runs, with and without competing risks", {
  fx <- setup_dp_fixture()
  data(pbc3, envir = environment())
  hist <- pbc3[pbc3$id %in% unique(pbc3$id)[1:60], ]
  hist$years <- NA
  tr_cr <- transportBJM(fx$survival_fit_all, "unlabeled", data_predict_all = hist,
                        long_fit_all = fx$long_fit_all, time_variable = "year",
                        survival_variable_all = fx$survival_variable_all,
                        survival_trans_function = fx$survival_trans_function, n_grid = 60)
  expect_true(is.finite(tr_cr$transport$alpha))
  expect_true(is.finite(tr_cr$transport$alpha_se))
  expect_equal(tr_cr$transport$cum_basehaz$hazard,
               survival::basehaz(fx$survival_fit_all$coxph_fit, centered = FALSE)$hazard *
                 exp(tr_cr$transport$alpha))
  ### the same histories as source reference: the shift is exactly 0
  tr_ref <- transportBJM(fx$survival_fit_all, "unlabeled", data_predict_all = hist,
                         long_fit_all = fx$long_fit_all, time_variable = "year",
                         survival_variable_all = fx$survival_variable_all,
                         survival_trans_function = fx$survival_trans_function, n_grid = 60,
                         data_source = hist)
  expect_equal(tr_ref$transport$alpha, 0)
  expect_output(print(tr_cr), "unlabeled")
})

test_that("transportBJM() rejects unsupported input", {
  data(pbc3, envir = environment())
  first <- pbc3[!duplicated(pbc3$id), ]
  src <- survivalSub(first, Surv(years, status3) ~ age + sex, NULL)
  expect_error(transportBJM(src, "baseline"), "data_target")
  expect_error(transportBJM(src, "unlabeled"), "data_predict_all")
  strat <- survivalSub(first, Surv(years, status3) ~ age + strata(sex), NULL)
  expect_error(transportBJM(strat, "baseline", data_target = first), "stratified")
  noev <- first; noev$status3 <- 0
  expect_error(transportBJM(src, "baseline", data_target = noev), "no events")
})
