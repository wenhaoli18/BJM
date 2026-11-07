# backwardTimeBasis(): a piecewise-linear basis |u - c| in the backward time
# u = T - s (event time minus visit time), added to long_sub_fixed as
# I(abs(years - year - c)) terms.
#
# The prediction functions evaluate the fixed-effects mean at every
# integration grid point of the event time by rebuilding the design from the
# formula on rows whose event-time column is set to the grid value (see
# fixed_mean_grid()). A backward-time term therefore needs no
# survival_variable_all/survival_trans_function entry. These tests check
# that the terms are evaluated row by row (with each row's own visit time)
# at each grid point, against means computed by hand from the coefficients.

test_that("backwardTimeBasis appends one |T - s - c| term per cut point", {
  f <- backwardTimeBasis(serBilir ~ year + age + sex + years, c(1, 3))
  expect_s3_class(f, "formula")
  expect_equal(attr(terms(f), "term.labels"),
               c("year", "age", "sex", "years",
                 "I(abs(years - year - 1))", "I(abs(years - year - 3))"))
  expect_identical(as.character(f[[2]]), "serBilir")
})

test_that("backwardTimeBasis keeps list names, the formula environment and custom column names", {
  env <- new.env()
  f1 <- serBilir ~ year + years
  environment(f1) <- env
  out <- backwardTimeBasis(list(a = f1, b = albumin ~ fuyr + T_ev), 0.5,
                           survival_variable = "T_ev", time_variable = "fuyr")[c("a", "b")]
  expect_named(out, c("a", "b"))
  expect_identical(environment(out$a), env)
  expect_true("I(abs(T_ev - fuyr - 0.5))" %in% attr(terms(out$b), "term.labels"))
  # the default column names are used for the first formula as well
  expect_true("I(abs(years - year - 0.5))" %in%
                attr(terms(backwardTimeBasis(list(a = f1), 0.5)$a), "term.labels"))
})

test_that("backwardTimeBasis rejects malformed input", {
  f <- serBilir ~ year + years
  expect_error(backwardTimeBasis(f, numeric(0)), "`cut_points` must be a non-empty numeric vector")
  expect_error(backwardTimeBasis(f, c(1, NA)), "`cut_points` must be a non-empty numeric vector")
  expect_error(backwardTimeBasis(f, c(1, Inf)), "`cut_points` must be a non-empty numeric vector")
  expect_error(backwardTimeBasis(f, c(1, 1)), "must not contain duplicates")
  expect_error(backwardTimeBasis("serBilir ~ year", 1), "must be a formula or a non-empty list")
  expect_error(backwardTimeBasis(list(), 1), "must be a formula or a non-empty list")
  expect_error(backwardTimeBasis(~ year + years, 1), "must be a two-sided formula")
  expect_error(backwardTimeBasis(f, 1, survival_variable = "year"), "two different columns")
  expect_error(backwardTimeBasis(backwardTimeBasis(f, 1), c(1, 2)), "already contains")
})

## a two-biomarker fit on pbc3 with a backward-time basis
build_backward_fit <- function() {
  data(pbc3, envir = environment())
  long_sub_fixed <- backwardTimeBasis(
    list("long1" = serBilir ~ year + age + sex + years,
         "long2" = albumin ~ year + age + sex + years),
    cut_points = c(1, 3)
  )
  long_sub_random <- list("long1" = ~ year | id, "long2" = ~ year | id)
  data_fit_all <- list(pbc3[pbc3$status3 == 1, ], pbc3[pbc3$status3 == 1, ])
  long_fit_all <- longitudinalSub(data_fit_all, long_sub_fixed, long_sub_random)
  data_survival_fitting <- pbc3[!duplicated(pbc3$id), ]
  survival_fit_all <- survivalSub(data_survival_fitting, Surv(years, status3) ~ age + sex, NULL)
  list(pbc3 = pbc3, long_fit_all = long_fit_all, survival_fit_all = survival_fit_all)
}

## the fixed-effects mean written out by hand from the coefficients, with the
## backward-time columns computed from each row's own visit time
manual_mean <- function(lfit_i, d, l) {
  b <- lfit_i$coefficients$fixed
  as.numeric(b[["(Intercept)"]] + b[["year"]] * d$year + b[["age"]] * d$age + b[["sex"]] * d$sex +
               b[["years"]] * l +
               b[["I(abs(years - year - 1))"]] * abs(l - d$year - 1) +
               b[["I(abs(years - year - 3))"]] * abs(l - d$year - 3))
}

test_that("fixed_mean_grid evaluates the backward-time basis row by row at every grid point", {
  fx <- build_backward_fit()
  lfit1 <- fx$long_fit_all$lfit[[1]]
  d <- fx$pbc3[fx$pbc3$id %in% c(2, 5), ]
  l_i <- c(4.5, 6, 9.25)
  # max_rows = 5 forces several chunks of stacked rows
  for (max_rows in c(2e5, 5)) {
    mu <- fixed_mean_grid(d, lfit1, NULL, l_i, "years", list(), list(), max_rows = max_rows)
    for (k in seq_along(l_i)) {
      expect_equal(mu[, k], manual_mean(lfit1, d, l_i[k]), tolerance = 1e-10)
    }
  }
  # the basis genuinely varies within a patient (different visit times)
  expect_gt(length(unique(round(abs(l_i[1] - d$year[d$id == 2] - 1), 8))), 1)
})

test_that("conditionalYT with a backward-time basis matches an independent dmvnorm ground truth", {
  fx <- build_backward_fit()
  long_fit_all <- fx$long_fit_all
  num_i <- 2
  d <- fx$pbc3[fx$pbc3$id == num_i & fx$pbc3$year <= 5, ]
  data_predict_all <- list(d, d)
  l_i <- c(5.5, 7, 10)

  out <- conditionalYT(data_predict_all, long_fit_all, l_i, survival_variable = "years",
                       time_variable = "year", survival_variable_all = list(),
                       survival_trans_function = list())

  patient <- select_patient_longitudinal_data(data_predict_all, "id", num_i, 2, "year")
  A <- random_effects_design(patient$data_num_i_list, long_fit_all$long_sub_random, long_fit_all$Sigma_fit)
  r_diag <- unlist(lapply(1:2, function(i) rep(long_fit_all$lfit[[i]]$sigma^2,
                                                nrow(patient$data_num_i_list[[i]]))))
  Sigma_all <- A %*% long_fit_all$Sigma_fit %*% t(A) + diag(r_diag, length(r_diag))
  y <- c(patient$data_num_i_list[[1]]$serBilir, patient$data_num_i_list[[2]]$albumin)
  for (k in seq_along(l_i)) {
    mu <- c(manual_mean(long_fit_all$lfit[[1]], patient$data_num_i_list[[1]], l_i[k]),
            manual_mean(long_fit_all$lfit[[2]], patient$data_num_i_list[[2]], l_i[k]))
    expect_equal(out[[1]][k, 1], mvtnorm::dmvnorm(y, mu, Sigma_all, log = TRUE), tolerance = 1e-8)
  }
})

test_that("predictRisk, predictLongitudinal and simulateTrajectory run with a backward-time basis", {
  fx <- build_backward_fit()
  d <- fx$pbc3[fx$pbc3$id %in% c(2, 3, 4) & fx$pbc3$year <= 3, ]
  data_predict_all <- list(d, d)

  risk <- predictRisk(data_predict_all, fx$long_fit_all, fx$survival_fit_all,
                      prediction_time = 3, horizon = 2, time_variable = "year",
                      survival_variable_all = list(), survival_trans_function = list(),
                      bandcount1 = 10, bandcount2 = 20)
  expect_s3_class(risk, "predictRisk.BJM")
  expect_true(all(is.finite(risk$risk_prob_1)))
  expect_true(all(risk$risk_prob_1 >= 0 & risk$risk_prob_1 <= 1))

  # the backward-time terms change the prediction (they are not silently dropped)
  data_fit_all <- list(fx$pbc3[fx$pbc3$status3 == 1, ], fx$pbc3[fx$pbc3$status3 == 1, ])
  long_plain <- longitudinalSub(data_fit_all,
                                list("long1" = serBilir ~ year + age + sex + years,
                                     "long2" = albumin ~ year + age + sex + years),
                                list("long1" = ~ year | id, "long2" = ~ year | id))
  risk_plain <- predictRisk(data_predict_all, long_plain, fx$survival_fit_all,
                            prediction_time = 3, horizon = 2, time_variable = "year",
                            survival_variable_all = list(), survival_trans_function = list(),
                            bandcount1 = 10, bandcount2 = 20)
  expect_false(isTRUE(all.equal(unname(risk$risk_prob_1), unname(risk_plain$risk_prob_1))))

  bio <- predictLongitudinal(bio_i = 1, data_predict_all, fx$long_fit_all, fx$survival_fit_all,
                             prediction_time = 3, horizon = 1, time_variable = "year",
                             survival_variable_all = list(), survival_trans_function = list(),
                             bandcount2 = 10, bandcount3 = 50)
  expect_s3_class(bio, "predictLongitudinal.BJM")

  sims <- simulateTrajectory(data_predict_all, fx$long_fit_all, fx$survival_fit_all,
                             prediction_time = 3, times = c(3.5, 4), n_sim = 5, seed = 1,
                             time_variable = "year", survival_variable_all = list(),
                             survival_trans_function = list())
  expect_s3_class(sims, "simulateTrajectory.BJM")
  expect_true(all(is.finite(sims$serBilir[!is.na(sims$serBilir)])))
})
