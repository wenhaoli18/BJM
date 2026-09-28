# Independent numeric ground-truth check for the joint Y|T (and Y|D,T)
# conditional density computed by conditionalYT()/conditionalYDT().
#
# Neither function previously had a test that verified its output against an
# INDEPENDENT reference -- the only numeric checks were golden-master
# snapshots of the package's own prior output (test-baseline-characterization
# (-noCR).R), which cannot catch a bug that was already present when the
# snapshot was captured.
#
# This matters because conditionalYT()/conditionalYDT() used to compute the
# joint (across all M jointly-fit biomarkers) Gaussian quadratic form
# (Y - mu)' Sigma_all_solve (Y - mu) via a "trace trick": reducing it to
# sum(diag(...)) of an M x M matrix built from parameter_matrix/
# longitudinal_all_matrix sandwiched through Sigma_all_solve. That trace only
# sums the *diagonal* (i == j) blocks of the reduced matrix, silently
# dropping every cross-marker (i != j) contribution of Sigma_all_solve --
# i.e. it implicitly treated the biomarkers as conditionally independent
# given the random effects, even though Sigma_fit (the fitted joint
# random-effects covariance) is generally NOT block-diagonal across markers.
# Since Sigma_all is built from Sigma_fit precisely to capture that
# cross-marker correlation, the bug silently threw away real information
# whenever two or more jointly-fit biomarkers had correlated random effects
# -- the ordinary case for a joint model, not an edge case.
#
# These tests independently recompute the same density via
# mvtnorm::dmvnorm() on the full stacked mean/covariance (rather than via
# conditionalYT.R/conditionalYDT.R's own internal matrix algebra), for a
# patient/fit where the fitted Sigma_fit has genuinely nonzero cross-marker
# blocks -- confirming the fix, and guarding against this class of bug
# recurring.

build_two_marker_fit <- function() {
  data(pbc3, envir = environment())

  long_sub_fixed <- list(
    "long1" = serBilir ~ year + age + sex + years,
    "long2" = albumin ~ year + age + sex + years
  )
  long_sub_random <- list("long1" = ~ year | id, "long2" = ~ year | id)
  data_fit_all <- list(pbc3[pbc3$status3 == 1, ], pbc3[pbc3$status3 == 1, ])
  long_fit_all <- longitudinalSub(data_fit_all, long_sub_fixed, long_sub_random)

  list(pbc3 = pbc3, long_fit_all = long_fit_all)
}

# Independently recomputes the stacked-vector Gaussian density that
# conditionalYT()/conditionalYDT() are supposed to compute for one patient
# at one l_i grid point, evaluated at a given value of any extra covariates
# (e.g. event_type_variable for the competing-risk path) already baked into
# `data_num_i_list`.
dmvnorm_reference <- function(long_fit_all, data_num_i_list, Sigma_all, l_i_value, survival_variable) {
  lfit <- long_fit_all$lfit
  n_longitudinal <- length(lfit)

  mu_parts <- list()
  y_parts <- list()
  for (i in seq_len(n_longitudinal)) {
    d_i <- data_num_i_list[[i]]
    d_i[[survival_variable]] <- l_i_value
    terms_model <- lfit[[i]]$terms
    xlev_i <- if (!is.null(long_fit_all$xlevels)) long_fit_all$xlevels[[i]] else NULL
    mf_i <- model.frame(terms_model, d_i, xlev = xlev_i)
    X_i <- model.matrix(terms_model, mf_i, contrasts.arg = lfit[[i]]$contrasts)
    mu_parts[[i]] <- as.numeric(X_i %*% lfit[[i]]$coefficients$fixed)
    y_parts[[i]] <- as.numeric(d_i[[as.character(formula(lfit[[i]]))[2]]])
  }
  mvtnorm::dmvnorm(unlist(y_parts), mean = unlist(mu_parts), sigma = Sigma_all)
}

test_that("conditionalYT matches an independent mvtnorm::dmvnorm ground truth for correlated markers", {
  fx <- build_two_marker_fit()
  long_fit_all <- fx$long_fit_all

  # sanity: the fitted joint random-effects covariance actually has nonzero
  # cross-marker blocks for this fixture -- otherwise this test wouldn't
  # exercise the bug at all (any quadratic-form implementation would agree
  # when Sigma_fit is block-diagonal, per the toy check in the PR).
  expect_true(any(abs(long_fit_all$Sigma_fit[1:2, 3:4]) > 1e-6))

  num_i <- 2
  data_predict_all <- list(fx$pbc3[fx$pbc3$id == num_i, ], fx$pbc3[fx$pbc3$id == num_i, ])
  l_i <- c(6, 7)

  out <- conditionalYT(data_predict_all, long_fit_all, l_i, survival_variable = "years",
                        time_variable = "year", survival_variable_all = list(),
                        survival_trans_function = list())

  n_longitudinal <- length(long_fit_all$lfit)
  patient_data <- select_patient_longitudinal_data(data_predict_all, "id", num_i, n_longitudinal, "year")
  design <- build_conditional_design(patient_data$rep_num_i_list, patient_data$data_num_i_list,
                                      long_fit_all$lfit, long_fit_all$Sigma_fit,
                                      sapply(long_fit_all$lfit, function(u) u$sigma),
                                      "year", n_longitudinal)

  for (it in seq_along(l_i)) {
    ref <- dmvnorm_reference(long_fit_all, patient_data$data_num_i_list, design$Sigma_all,
                              l_i[it], "years")
    # conditionalYT() returns an unnamed list (return(f_Y_T_D = list(...))
    # discards the argument name since return() takes a single value), so
    # the density matrix is out[[1]], not out$f_Y_T_D.
    expect_equal(out[[1]][it, 1], ref, tolerance = 1e-8)
  }
})

test_that("conditionalYDT matches an independent mvtnorm::dmvnorm ground truth for correlated markers", {
  fx <- build_two_marker_fit()
  long_fit_all <- fx$long_fit_all
  pbc3 <- fx$pbc3

  data_survival_fitting <- pbc3[!duplicated(pbc3$id), ]
  survival_fit_all <- survivalSub(data_survival_fitting,
                                   Surv(years, status3) ~ age + sex,
                                   status4 ~ years + age + sex)

  num_i <- 2
  data_predict_all <- list(pbc3[pbc3$id == num_i, ], pbc3[pbc3$id == num_i, ])
  l_i <- c(6, 7)

  out <- conditionalYDT(data_predict_all, long_fit_all, survival_fit_all, l_i,
                         survival_variable = "years", time_variable = "year",
                         survival_variable_all = list(), survival_trans_function = list())

  n_longitudinal <- length(long_fit_all$lfit)
  event_type_variable <- as.character(formula(survival_fit_all$form_conditional_cr)[[2]])
  patient_data <- select_patient_longitudinal_data(data_predict_all, "id", num_i, n_longitudinal, "year")
  design <- build_conditional_design(patient_data$rep_num_i_list, patient_data$data_num_i_list,
                                      long_fit_all$lfit, long_fit_all$Sigma_fit,
                                      sapply(long_fit_all$lfit, function(u) u$sigma),
                                      "year", n_longitudinal)

  for (it in seq_along(l_i)) {
    data_num_i_list_1 <- lapply(patient_data$data_num_i_list, function(d) {
      d[[event_type_variable]] <- 1
      d
    })
    data_num_i_list_0 <- lapply(patient_data$data_num_i_list, function(d) {
      d[[event_type_variable]] <- 0
      d
    })
    ref1 <- dmvnorm_reference(long_fit_all, data_num_i_list_1, design$Sigma_all, l_i[it], "years")
    ref0 <- dmvnorm_reference(long_fit_all, data_num_i_list_0, design$Sigma_all, l_i[it], "years")

    # conditionalYDT() returns an unnamed list(f_Y_T_D_w0, f_Y_T_D_w1) (see
    # note above on return()'s argument name being discarded).
    expect_equal(out[[2]][it, 1], ref1, tolerance = 1e-8)
    expect_equal(out[[1]][it, 1], ref0, tolerance = 1e-8)
  }
})
