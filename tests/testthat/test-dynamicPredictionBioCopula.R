# Prediction-side tests for the Gaussian-copula (mixed continuous/ordinal)
# extension applied to *biomarker-value* prediction:
# conditionalYTBioCopula()/conditionalYDTBioCopula() and their dispatch from
# dynamicPredictionBio(). See test-predictRiskCopula.R for the
# event-risk-prediction counterpart (predictRisk()); this file covers
# predicting a future biomarker *value* from a fitted mixed model instead.
#
# Both a continuous and an ordinal bio_i target are supported: for an
# ordinal target, conditionalYTBioCopula()/conditionalYDTBioCopula() bracket
# the candidate category between its cumulative-link thresholds via the same
# generic biomarker_type-driven machinery used for every other ordinal row
# in the joint density (see conditionalDesignCopula.R), and
# compute_bio_marker_step() (R/dynamicPredictionBio.R) exposes the result as
# integer category codes rather than a numeric grid -- see the test below.
#
# Same independent-reconstruction methodology and tolerance rationale as
# test-predictRiskCopula.R (mvtnorm::pmvnorm()'s default GenzBretz
# algorithm is Monte Carlo, so both the package code and the reference are
# individually stochastic).

setup_copula_predict_fixture <- function() {
  data(pbc3, envir = environment())
  d <- pbc3[pbc3$status3 == 1, ]
  breaks <- stats::quantile(d$albumin, c(0, 1 / 3, 2 / 3, 1), na.rm = TRUE)
  to_cat <- function(dd) {
    dd$albumin_cat <- cut(dd$albumin, breaks = breaks, include.lowest = TRUE,
                           labels = c("low", "mid", "high"), ordered_result = TRUE)
    dd
  }
  d <- to_cat(d)

  long_sub_fixed <- list("m1" = serBilir ~ year + age + sex + years,
                          "m2" = albumin_cat ~ year + age + sex + years)
  long_sub_random <- list("m1" = ~ year | id, "m2" = ~ year | id)
  long_fit_all <- longitudinalSub(list(d, d), long_sub_fixed, long_sub_random)

  list(pbc3 = pbc3, breaks = breaks, to_cat = to_cat,
       long_fit_all = long_fit_all, long_sub_fixed = long_sub_fixed)
}

skip_if_ordinal <- function() testthat::skip_if_not_installed("ordinal")

# Independently recomputes conditionalYTBioCopula()'s target quantity for
# one patient/candidate biomarker value: the same mixed
# Gaussian-density-times-copula-box-probability construction as
# test-predictRiskCopula.R's reference_mixed_density(), but with the
# candidate value substituted for bio_i's added row instead of using the
# truly observed value there.
reference_bio_density <- function(long_fit_all, long_sub_fixed, num_i, data_predict_all,
                                   bio_i, l_i_value, time_new, Y_candidate) {
  lfit <- long_fit_all$lfit
  Sigma <- long_fit_all$Sigma_fit
  sigma.longitudinal <- c(lfit[[1]]$sigma, 1)

  patient_data <- select_patient_longitudinal_data(data_predict_all, "id", num_i, 2, "year")
  data_it <- patient_data$data_num_i_list
  for (i in 1:2) data_it[[i]]$years <- l_i_value

  bio_i_name <- as.character(formula(long_sub_fixed[[bio_i]])[[2]])
  data_it[[bio_i]] <- rbind(data_it[[bio_i]], data_it[[bio_i]][nrow(data_it[[bio_i]]), ])
  data_it[[bio_i]]$year[nrow(data_it[[bio_i]])] <- time_new
  data_it[[bio_i]][[bio_i_name]][nrow(data_it[[bio_i]])] <- Y_candidate
  rep_i <- patient_data$rep_num_i_list
  rep_i[[bio_i]] <- rep(1, nrow(data_it[[bio_i]]))

  X1 <- model.matrix(terms(lfit[[1]]), model.frame(terms(lfit[[1]]), data_it[[1]]))
  mu1 <- as.numeric(X1 %*% lfit[[1]]$coefficients$fixed)
  y1 <- data_it[[1]]$serBilir

  beta_names_2 <- names(lfit[[2]]$beta)
  X2 <- model.matrix(long_sub_fixed[[2]], data_it[[2]])[, beta_names_2, drop = FALSE]
  mu2 <- as.numeric(X2 %*% lfit[[2]]$beta)
  y2_cat <- as.numeric(data_it[[2]]$albumin_cat)
  full_alpha <- c(-Inf, lfit[[2]]$alpha, Inf)
  lower2 <- full_alpha[y2_cat]
  upper2 <- full_alpha[y2_cat + 1]

  A1 <- rbind(rep_i[[1]], data_it[[1]]$year)
  A2 <- rbind(rep_i[[2]], data_it[[2]]$year)
  n1 <- ncol(A1); n2 <- ncol(A2)
  A <- matrix(0, n1 + n2, 4)
  A[1:n1, 1:2] <- t(A1)
  A[(n1 + 1):(n1 + n2), 3:4] <- t(A2)
  sigvec <- c(rep(sigma.longitudinal[1]^2, n1), rep(sigma.longitudinal[2]^2, n2))
  Sigma_all <- A %*% Sigma %*% t(A) + diag(sigvec, length(sigvec))

  cc_idx <- 1:n1
  oo_idx <- (n1 + 1):(n1 + n2)
  Sigma_cc <- Sigma_all[cc_idx, cc_idx, drop = FALSE]
  Sigma_oo <- Sigma_all[oo_idx, oo_idx, drop = FALSE]
  Sigma_oc <- Sigma_all[oo_idx, cc_idx, drop = FALSE]
  Sigma_co <- Sigma_all[cc_idx, oo_idx, drop = FALSE]

  dens_c <- mvtnorm::dmvnorm(y1, mean = mu1, sigma = Sigma_cc)
  cond_mean_o <- mu2 + as.numeric(Sigma_oc %*% solve(Sigma_cc) %*% (y1 - mu1))
  cond_cov_o <- Sigma_oo - Sigma_oc %*% solve(Sigma_cc) %*% Sigma_co
  prob_o <- as.numeric(mvtnorm::pmvnorm(lower = lower2, upper = upper2, mean = cond_mean_o, sigma = cond_cov_o))

  dens_c * prob_o
}

test_that("conditionalYTBioCopula matches an independent from-scratch reconstruction", {
  skip_if_ordinal()
  fx <- setup_copula_predict_fixture()

  num_i <- 2
  raw <- fx$pbc3[fx$pbc3$id == num_i & fx$pbc3$year <= 3, ]
  data_predict_all <- list(fx$to_cat(raw), fx$to_cat(raw))
  bio_i <- 1
  l_i <- 6
  time_new <- 6
  Y_all <- c(1.5, 2.0, 2.5)

  set.seed(20260927)
  out <- conditionalYTBioCopula(Y_all, time_new, bio_i, data_predict_all, fx$long_fit_all, l_i,
                                 survival_variable = "years", time_variable = "year",
                                 survival_variable_all = list(), survival_trans_function = list())

  set.seed(20260927)
  for (Y_i in seq_along(Y_all)) {
    ref <- reference_bio_density(fx$long_fit_all, fx$long_sub_fixed, num_i, data_predict_all,
                                  bio_i, l_i, time_new, Y_all[Y_i])
    expect_equal(out[[1]][[Y_i]][1, 1], ref, tolerance = 1e-2)
  }
})

test_that("conditionalYDTBioCopula matches an independent from-scratch reconstruction (w0/w1)", {
  skip_if_ordinal()
  fx <- setup_copula_predict_fixture()

  data_survival_fitting <- fx$pbc3[!duplicated(fx$pbc3$id), ]
  survival_fit_all <- survivalSub(data_survival_fitting, Surv(years, status3) ~ age + sex,
                                   status4 ~ years + age + sex)

  num_i <- 2
  raw <- fx$pbc3[fx$pbc3$id == num_i & fx$pbc3$year <= 3, ]
  data_predict_all <- list(fx$to_cat(raw), fx$to_cat(raw))
  bio_i <- 1
  l_i <- 6
  time_new <- 6
  Y_all <- c(1.5, 2.0, 2.5)

  set.seed(20260927)
  out <- conditionalYDTBioCopula(Y_all, time_new, bio_i, data_predict_all, fx$long_fit_all,
                                  survival_fit_all, l_i, survival_variable = "years",
                                  time_variable = "year", survival_variable_all = list(),
                                  survival_trans_function = list())

  # event_type_variable ("status4") does not appear in long_sub_fixed here,
  # so it has no effect on the mean/design -- the w0/w1 branches should
  # therefore agree with each other and with the (event-type-agnostic)
  # reference reconstruction used above.
  set.seed(20260927)
  for (Y_i in seq_along(Y_all)) {
    ref <- reference_bio_density(fx$long_fit_all, fx$long_sub_fixed, num_i, data_predict_all,
                                  bio_i, l_i, time_new, Y_all[Y_i])
    expect_equal(out[[1]][[Y_i]][1, 1], ref, tolerance = 1e-2)
    expect_equal(out[[2]][[Y_i]][1, 1], ref, tolerance = 1e-2)
  }
})

test_that("dynamicPredictionBio dispatches to the copula path for a mixed fit and returns a sane prediction", {
  skip_if_ordinal()
  fx <- setup_copula_predict_fixture()

  data_survival_fitting <- fx$pbc3[!duplicated(fx$pbc3$id), ]
  survival_fit_all <- survivalSub(data_survival_fitting, Surv(years, status3) ~ age + sex, NULL)

  raw <- fx$pbc3[fx$pbc3$id == 2 & fx$pbc3$year <= 3, ]
  data_predict_all <- list(fx$to_cat(raw), fx$to_cat(raw))

  ypred <- dynamicPredictionBio(bio_i = 1, data_predict_all, fx$long_fit_all, survival_fit_all,
                                 prediction_time = 3, horizon = 3, time_variable = "year",
                                 survival_variable_all = list(), survival_trans_function = list(),
                                 bandcount2 = 10, bandcount3 = 15)

  expect_s3_class(ypred, "dynamicPredictionBio.BJM")
  expect_true(is.finite(ypred$Y_predict))
  expect_true(all(is.finite(ypred$Y_density)))
})

test_that("dynamicPredictionBio dispatches to the copula path under competing risk and returns a sane prediction", {
  skip_if_ordinal()
  fx <- setup_copula_predict_fixture()

  data_survival_fitting <- fx$pbc3[!duplicated(fx$pbc3$id), ]
  survival_fit_all <- survivalSub(data_survival_fitting, Surv(years, status3) ~ age + sex,
                                   status4 ~ years + age + sex)

  raw <- fx$pbc3[fx$pbc3$id == 2 & fx$pbc3$year <= 3, ]
  data_predict_all <- list(fx$to_cat(raw), fx$to_cat(raw))

  ypred <- dynamicPredictionBio(bio_i = 1, data_predict_all, fx$long_fit_all, survival_fit_all,
                                 prediction_time = 3, horizon = 3, time_variable = "year",
                                 survival_variable_all = list(), survival_trans_function = list(),
                                 bandcount2 = 10, bandcount3 = 15)

  expect_s3_class(ypred, "dynamicPredictionBio.BJM")
  expect_true(is.finite(ypred$Y_predict))
  expect_true(all(is.finite(ypred$Y_density)))
})

test_that("dynamicPredictionBio on an ordinal bio_i predicts its future category", {
  skip_if_ordinal()
  fx <- setup_copula_predict_fixture()

  data_survival_fitting <- fx$pbc3[!duplicated(fx$pbc3$id), ]
  survival_fit_all <- survivalSub(data_survival_fitting, Surv(years, status3) ~ age + sex, NULL)

  raw <- fx$pbc3[fx$pbc3$id == 2 & fx$pbc3$year <= 3, ]
  data_predict_all <- list(fx$to_cat(raw), fx$to_cat(raw))

  # bandcount3 = "auto" is passed through unused for an ordinal bio_i (its
  # candidate grid is fixed at its category count -- see
  # compute_bio_marker_step()), so this also exercises that bandcount3 need
  # not be a valid numeric doubling target in that case.
  ypred <- dynamicPredictionBio(bio_i = 2, data_predict_all, fx$long_fit_all, survival_fit_all,
                                 prediction_time = 3, horizon = 3, time_variable = "year",
                                 survival_variable_all = list(), survival_trans_function = list(),
                                 bandcount2 = 10, bandcount3 = "auto")

  expect_s3_class(ypred, "dynamicPredictionBio.BJM")
  # Y_all/Y_predict are integer category codes (1:3 for low/mid/high), not a
  # numeric grid.
  expect_equal(as.integer(ypred$Y_all), 1:3)
  expect_equal(attr(ypred$Y_all, "category_labels"), c("low", "mid", "high"))
  expect_true(ypred$Y_predict %in% 1:3)
  expect_true(all(is.finite(ypred$Y_density)))
  # Y_density's rows are per-category probabilities, so they sum to ~1 for
  # this single patient.
  expect_equal(sum(ypred$Y_density[, 1]), 1, tolerance = 0.1)
})

test_that("an all-continuous fit's dynamicPredictionBio is unaffected by the copula dispatch (hard requirement)", {
  data(pbc3, envir = environment())
  data_fit_all <- pbc3[pbc3$status3 == 1, ]
  long_sub_fixed <- list("m1" = serBilir ~ year + age + sex + years, "m2" = albumin ~ year + age + sex + years)
  long_sub_random <- list("m1" = ~ year | id, "m2" = ~ year | id)
  long_fit_all <- longitudinalSub(list(data_fit_all, data_fit_all), long_sub_fixed, long_sub_random)
  expect_null(long_fit_all$biomarker_type) # confirms use_copula would be FALSE

  data_survival_fitting <- pbc3[!duplicated(pbc3$id), ]
  survival_fit_all <- survivalSub(data_survival_fitting, Surv(years, status3) ~ age + sex, NULL)
  data_predict_all <- list(pbc3[pbc3$id == 2 & pbc3$year <= 3, ], pbc3[pbc3$id == 2 & pbc3$year <= 3, ])

  ypred <- dynamicPredictionBio(bio_i = 1, data_predict_all, long_fit_all, survival_fit_all,
                                 prediction_time = 3, horizon = 3, time_variable = "year",
                                 survival_variable_all = list(), survival_trans_function = list(),
                                 bandcount2 = 10, bandcount3 = 15)
  expect_s3_class(ypred, "dynamicPredictionBio.BJM")
  expect_true(is.finite(ypred$Y_predict))
})
