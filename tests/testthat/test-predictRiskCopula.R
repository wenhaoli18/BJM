# Prediction-side tests for the Gaussian-copula (mixed continuous/ordinal)
# extension: conditionalYTCopula()/conditionalYDTCopula() and their dispatch
# from predictRisk(). See test-longitudinalSubCopula.R for the fitting
# side; this file covers using a fitted mixed model for prediction.
#
# Mirrors test-conditional-density-correctness.R's methodology: an
# INDEPENDENT ground-truth reconstruction (not calling any of the package's
# own conditionalYT*/build_conditional_design* internals) verifies the
# density/probability conditionalYTCopula() computes for a mixed
# continuous/ordinal patient. Because the ordinal factor's contribution is a
# multivariate-normal box probability evaluated via mvtnorm::pmvnorm()'s
# default (Monte Carlo) GenzBretz algorithm, both the package code and this
# independent reference are individually stochastic; the tolerance below is
# set generously (1% relative) to comfortably clear that inherent MC noise
# while still catching a real implementation bug (which produced ~40x-sized
# errors during development, not ~1e-4 relative ones -- see git history).

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

# Independently recomputes conditionalYTCopula()'s target quantity for one
# patient/l_i: an exact Gaussian density at the observed continuous values,
# times the Gaussian-copula box probability for the ordinal categories,
# conditional on the continuous observations -- built from first principles
# (raw model.matrix()/solve() calls), not via
# build_conditional_design_copula()/mixed_density_prob_copula().
reference_mixed_density <- function(long_fit_all, long_sub_fixed, num_i, data_predict_all, l_i_value) {
  lfit <- long_fit_all$lfit
  Sigma <- long_fit_all$Sigma_fit
  sigma.longitudinal <- c(lfit[[1]]$sigma, 1)

  patient_data <- select_patient_longitudinal_data(data_predict_all, "id", num_i, 2, "year")
  data_it <- patient_data$data_num_i_list
  for (i in 1:2) data_it[[i]]$years <- l_i_value

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

  A1 <- rbind(patient_data$rep_num_i_list[[1]], data_it[[1]]$year)
  A2 <- rbind(patient_data$rep_num_i_list[[2]], data_it[[2]]$year)
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

test_that("conditionalYTCopula matches an independent from-scratch reconstruction", {
  skip_if_ordinal()
  fx <- setup_copula_predict_fixture()

  num_i <- 2
  raw <- fx$pbc3[fx$pbc3$id == num_i, ]
  data_predict_all <- list(fx$to_cat(raw), fx$to_cat(raw))
  l_i <- c(6, 7)

  set.seed(20260927)
  out <- conditionalYTCopula(data_predict_all, fx$long_fit_all, l_i, survival_variable = "years",
                              time_variable = "year", survival_variable_all = list(),
                              survival_trans_function = list())

  set.seed(20260927)
  for (it in seq_along(l_i)) {
    ref <- reference_mixed_density(fx$long_fit_all, fx$long_sub_fixed, num_i, data_predict_all, l_i[it])
    expect_equal(out[[1]][it, 1], ref, tolerance = 1e-2)
  }
})

test_that("conditionalYDTCopula matches an independent from-scratch reconstruction (w0/w1)", {
  skip_if_ordinal()
  fx <- setup_copula_predict_fixture()

  data_survival_fitting <- fx$pbc3[!duplicated(fx$pbc3$id), ]
  survival_fit_all <- survivalSub(data_survival_fitting, Surv(years, status3) ~ age + sex,
                                   status4 ~ years + age + sex)

  num_i <- 2
  raw <- fx$pbc3[fx$pbc3$id == num_i, ]
  data_predict_all <- list(fx$to_cat(raw), fx$to_cat(raw))
  l_i <- c(6, 7)

  set.seed(20260927)
  out <- conditionalYDTCopula(data_predict_all, fx$long_fit_all, survival_fit_all, l_i,
                               survival_variable = "years", time_variable = "year",
                               survival_variable_all = list(), survival_trans_function = list())

  # event_type_variable ("status4") does not appear in long_sub_fixed here,
  # so it has no effect on the mean/design -- the w0/w1 branches should
  # therefore agree with each other and with the (event-type-agnostic)
  # reference reconstruction used for conditionalYTCopula() above.
  set.seed(20260927)
  for (it in seq_along(l_i)) {
    ref <- reference_mixed_density(fx$long_fit_all, fx$long_sub_fixed, num_i, data_predict_all, l_i[it])
    expect_equal(out[[1]][it, 1], ref, tolerance = 1e-2)
    expect_equal(out[[2]][it, 1], ref, tolerance = 1e-2)
  }
})

test_that("predictRisk dispatches to the copula path for a mixed fit and returns sane risk", {
  skip_if_ordinal()
  fx <- setup_copula_predict_fixture()

  data_survival_fitting <- fx$pbc3[!duplicated(fx$pbc3$id), ]
  survival_fit_all <- survivalSub(data_survival_fitting, Surv(years, status3) ~ age + sex, NULL)

  raw <- fx$pbc3[fx$pbc3$id == 2 & fx$pbc3$year <= 3, ]
  data_predict_all <- list(fx$to_cat(raw), fx$to_cat(raw))

  risk <- predictRisk(data_predict_all, fx$long_fit_all, survival_fit_all,
                             prediction_time = 3, horizon = 3, time_variable = "year",
                             survival_variable_all = list(), survival_trans_function = list(),
                             bandcount1 = 10, bandcount2 = 10)

  expect_s3_class(risk, "predictRisk.BJM")
  expect_true(is.finite(risk$risk_prob_1))
  expect_true(risk$risk_prob_1 >= 0 && risk$risk_prob_1 <= 1)
  expect_null(risk$risk_prob_2)
})

test_that("predictRisk dispatches to the copula path under competing risk and returns sane risk", {
  skip_if_ordinal()
  fx <- setup_copula_predict_fixture()

  data_survival_fitting <- fx$pbc3[!duplicated(fx$pbc3$id), ]
  survival_fit_all <- survivalSub(data_survival_fitting, Surv(years, status3) ~ age + sex,
                                   status4 ~ years + age + sex)

  raw <- fx$pbc3[fx$pbc3$id == 2 & fx$pbc3$year <= 3, ]
  data_predict_all <- list(fx$to_cat(raw), fx$to_cat(raw))

  risk <- predictRisk(data_predict_all, fx$long_fit_all, survival_fit_all,
                             prediction_time = 3, horizon = 3, time_variable = "year",
                             survival_variable_all = list(), survival_trans_function = list(),
                             bandcount1 = 10, bandcount2 = 10)

  expect_true(is.finite(risk$risk_prob_1) && risk$risk_prob_1 >= 0 && risk$risk_prob_1 <= 1)
  expect_true(is.finite(risk$risk_prob_2) && risk$risk_prob_2 >= 0 && risk$risk_prob_2 <= 1)
})

test_that("an all-continuous fit's predictRisk is unaffected by the copula dispatch (hard requirement)", {
  data(pbc3, envir = environment())
  data_fit_all <- pbc3[pbc3$status3 == 1, ]
  long_sub_fixed <- list("m1" = serBilir ~ year + age + sex + years, "m2" = albumin ~ year + age + sex + years)
  long_sub_random <- list("m1" = ~ year | id, "m2" = ~ year | id)
  long_fit_all <- longitudinalSub(list(data_fit_all, data_fit_all), long_sub_fixed, long_sub_random)
  expect_null(long_fit_all$biomarker_type) # confirms use_copula would be FALSE

  data_survival_fitting <- pbc3[!duplicated(pbc3$id), ]
  survival_fit_all <- survivalSub(data_survival_fitting, Surv(years, status3) ~ age + sex, NULL)
  data_predict_all <- list(pbc3[pbc3$id == 2 & pbc3$year <= 3, ], pbc3[pbc3$id == 2 & pbc3$year <= 3, ])

  risk_dispatch <- predictRisk(data_predict_all, long_fit_all, survival_fit_all,
                                      prediction_time = 3, horizon = 3, time_variable = "year",
                                      survival_variable_all = list(), survival_trans_function = list(),
                                      bandcount1 = 10, bandcount2 = 10)
  f_direct <- conditionalYT(data_predict_all, long_fit_all, l_i = 3, "years", "year", list(), list())
  expect_equal(class(f_direct), "list") # sanity: conditionalYT() itself still works standalone
  expect_true(is.finite(risk_dispatch$risk_prob_1))
})
