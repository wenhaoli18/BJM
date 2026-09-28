# Exercises poolLongitudinalSub(), the Rubin's-rules combination step for
# longitudinalSub() fits across the multiple completed datasets returned by
# imputeLongitudinal(..., impute = "multiple"). Input-validation and
# numeric-correctness tests below only call longitudinalSub() on small pbc3
# subsets (or hand-built fake fit objects) -- no torch is needed for any of
# that. Only the integration test at the bottom needs torch, to actually
# produce multiple imputed completions end to end.

# A handful of cheap, valid longitudinalSub.BJM fits (real nlme::lme() fits
# on small, different pbc3 subject subsets, so their estimates/SEs/dfcom
# genuinely differ -- like independent imputed completions would), reused
# by the validation tests below, which just need *some* valid fit objects
# to build malformed lists/mismatches around.
setup_pool_fixture <- function() {
  data(pbc3, envir = environment())
  d <- pbc3[pbc3$status3 == 1, ]
  ids <- unique(d$id)

  long_sub_fixed_2m <- list(
    "long1" = serBilir ~ year + age + sex,
    "long2" = albumin ~ year + age + sex
  )
  long_sub_random_2m <- list("long1" = ~ year | id, "long2" = ~ year | id)

  fit1 <- longitudinalSub(d[d$id %in% ids[1:80], ], long_sub_fixed_2m, long_sub_random_2m)
  fit2 <- longitudinalSub(d[d$id %in% ids[1:75], ], long_sub_fixed_2m, long_sub_random_2m)
  fit3 <- longitudinalSub(d[d$id %in% ids[1:85], ], long_sub_fixed_2m, long_sub_random_2m)

  list(fit1 = fit1, fit2 = fit2, fit3 = fit3, d = d, ids = ids,
       long_sub_fixed_2m = long_sub_fixed_2m, long_sub_random_2m = long_sub_random_2m)
}

# Builds a fake object with exactly the shape poolLongitudinalSub() reads
# from a longitudinalSub.BJM fit -- fixef.lme(object) is just
# `object$coefficients$fixed` (confirmed against the installed nlme), so a
# plain list tagged class "lme" round-trips through nlme::fixef() without
# needing a real model fit -- letting the numeric-correctness test below
# specify exact, known estimates/SEs/dfcom instead of whatever a real fit
# happens to produce.
make_fake_lme <- function(fixed, se, dfcom) {
  structure(
    list(coefficients = list(fixed = fixed),
         varFix = diag(se^2, nrow = length(se)),
         fixDF = list(terms = stats::setNames(rep(dfcom, length(se)), names(fixed)))),
    class = "lme"
  )
}

make_fake_longitudinal_fit <- function(fixed, se, dfcom, long_sub_fixed) {
  structure(
    list(lfit = list(make_fake_lme(fixed, se, dfcom)),
         Sigma_fit = NULL, long_sub_fixed = long_sub_fixed,
         long_sub_random = NULL, xlevels = NULL),
    class = "longitudinalSub.BJM"
  )
}

test_that("poolLongitudinalSub rejects a non-list long_fit_all_list", {
  expect_error(poolLongitudinalSub(5), "must be a list of")
})

test_that("poolLongitudinalSub rejects a long_fit_all_list with fewer than 2 fits", {
  f <- setup_pool_fixture()
  expect_error(poolLongitudinalSub(list(f$fit1)), "at least 2")
})

test_that("poolLongitudinalSub rejects a list element that is not a longitudinalSub.BJM object", {
  f <- setup_pool_fixture()
  expect_error(poolLongitudinalSub(list(f$fit1, 5)), "must be the output of longitudinalSub")
})

test_that("poolLongitudinalSub rejects fits with different numbers of longitudinal outcomes", {
  f <- setup_pool_fixture()
  fit_1m <- longitudinalSub(f$d[f$d$id %in% f$ids[1:80], ],
                             list("long1" = serBilir ~ year + age + sex),
                             list("long1" = ~ year | id))

  expect_error(poolLongitudinalSub(list(f$fit1, fit_1m)), "same set of biomarkers")
})

test_that("poolLongitudinalSub rejects fits with different fixed-effect coefficient names", {
  f <- setup_pool_fixture()
  fit_a <- longitudinalSub(f$d[f$d$id %in% f$ids[1:80], ],
                            list("long1" = serBilir ~ year + age + sex),
                            list("long1" = ~ year | id))
  fit_b <- longitudinalSub(f$d[f$d$id %in% f$ids[1:80], ],
                            list("long1" = serBilir ~ year + age + sex + albumin),
                            list("long1" = ~ year | id))

  expect_error(poolLongitudinalSub(list(fit_a, fit_b)), "different fixed-effect coefficient names")
})

test_that("poolLongitudinalSub reproduces textbook Rubin's rules / Barnard-Rubin df exactly", {
  fixed <- c("(Intercept)" = 1.0)
  long_sub_fixed <- list("long1" = y ~ 1)

  estimates <- c(1.0, 1.2, 0.8)
  se <- c(0.1, 0.1, 0.1)
  dfcom <- 100

  fits <- lapply(estimates, function(est) {
    make_fake_longitudinal_fit(stats::setNames(est, "(Intercept)"), se[1], dfcom, long_sub_fixed)
  })

  pooled <- poolLongitudinalSub(fits)
  expect_s3_class(pooled, "poolLongitudinalSub.BJM")
  expect_equal(pooled$m, 3)
  tab <- pooled$pooled[["long1"]]
  expect_equal(rownames(tab), "(Intercept)")

  # hand-computed reference values (Rubin 1987 Ch. 3; Barnard & Rubin 1999),
  # written out independently of poolLongitudinalSub()'s own implementation
  m <- 3
  qbar <- mean(estimates)
  ubar <- mean(se^2)
  b <- stats::var(estimates)
  t_var <- ubar + (1 + 1 / m) * b
  se_expected <- sqrt(t_var)
  riv_expected <- (1 + 1 / m) * b / ubar
  lambda <- (1 + 1 / m) * b / t_var
  df_old <- (m - 1) / lambda^2
  df_obs <- (dfcom + 1) / (dfcom + 3) * dfcom * (1 - lambda)
  df_expected <- df_old * df_obs / (df_old + df_obs)
  fmi_expected <- (riv_expected + 2 / (dfcom + 3)) / (riv_expected + 1)
  stat_expected <- qbar / se_expected
  p_expected <- 2 * stats::pt(-abs(stat_expected), df = df_expected)

  expect_equal(tab$estimate, qbar)
  expect_equal(tab$std_error, se_expected)
  expect_equal(tab$riv, riv_expected)
  expect_equal(tab$df, df_expected)
  expect_equal(tab$fmi, fmi_expected)
  expect_equal(tab$statistic, stat_expected)
  expect_equal(tab$p_value, p_expected)
})

test_that("poolLongitudinalSub with identical per-completion estimates gives zero between-imputation variance", {
  # b = 0 exercises the lambda-clamping guard in barnard_rubin_df() (a
  # coefficient perfectly agreed on by every completion still has to return
  # a finite, sensible df, not divide by zero)
  fixed <- c("(Intercept)" = 2.0)
  long_sub_fixed <- list("long1" = y ~ 1)
  fits <- lapply(1:4, function(i) {
    make_fake_longitudinal_fit(fixed, 0.2, dfcom = 50, long_sub_fixed)
  })

  pooled <- poolLongitudinalSub(fits)
  tab <- pooled$pooled[["long1"]]
  expect_equal(tab$estimate, 2.0)
  expect_equal(tab$std_error, 0.2)
  expect_true(is.finite(tab$df))
  expect_gt(tab$df, 0)
  expect_equal(tab$riv, 0)
})

test_that("poolLongitudinalSub warns and uses the smallest dfcom when completions disagree", {
  fixed <- c("(Intercept)" = 1.0)
  long_sub_fixed <- list("long1" = y ~ 1)
  fit_a <- make_fake_longitudinal_fit(fixed, 0.1, dfcom = 100, long_sub_fixed)
  fit_b <- make_fake_longitudinal_fit(stats::setNames(1.1, "(Intercept)"), 0.1, dfcom = 40, long_sub_fixed)

  expect_warning(pooled <- poolLongitudinalSub(list(fit_a, fit_b)), "differ across")
  # dfcom feeds df_obs upward with dfcom, so using the smaller (40, not
  # 100) value should show up as a smaller resulting pooled df than if the
  # larger value had been used instead
  fixed_100 <- make_fake_longitudinal_fit(fixed, 0.1, dfcom = 100, long_sub_fixed)
  fixed_100b <- make_fake_longitudinal_fit(stats::setNames(1.1, "(Intercept)"), 0.1, dfcom = 100, long_sub_fixed)
  pooled_100 <- suppressWarnings(poolLongitudinalSub(list(fixed_100, fixed_100b)))
  expect_lt(pooled$pooled[["long1"]]$df, pooled_100$pooled[["long1"]]$df)
})

test_that("print.poolLongitudinalSub.BJM runs without error and reports FMI", {
  f <- setup_pool_fixture()
  # f$fit1/fit2/fit3 are deliberately fit on different-sized subject subsets
  # (to give genuinely different per-completion estimates/SEs), which also
  # means their dfcom values legitimately differ -- expected here, and
  # already covered by its own test above, so it is suppressed rather than
  # re-asserted in this print-focused test.
  pooled <- suppressWarnings(poolLongitudinalSub(list(f$fit1, f$fit2, f$fit3)))
  expect_output(print(pooled), "Rubin's-Rules-Pooled")
  expect_output(print(pooled), "Fraction of missing information")
})

test_that("poolLongitudinalSub(imputeLongitudinal(impute = 'multiple') completions) pools end to end (requires torch)", {
  skip_if_not_installed("torch")
  f <- setup_impute_fixture()

  imputed <- imputeLongitudinal(f$data_fit_all, f$long_sub_fixed, f$long_sub_random,
                                 f$time_variable, n_imputations = 3, latent_dim = 4,
                                 hidden_units = c(16, 8), epochs = 15,
                                 importance_samples = 5, impute = "multiple", seed = 5)

  long_fit_all_list <- lapply(imputed$data_fit_all_list, function(d) {
    longitudinalSub(d, f$long_sub_fixed, f$long_sub_random)
  })

  pooled <- poolLongitudinalSub(long_fit_all_list)
  expect_s3_class(pooled, "poolLongitudinalSub.BJM")
  expect_equal(pooled$m, 3)
  expect_named(pooled$pooled, names(f$long_sub_fixed))

  for (tab in pooled$pooled) {
    expect_true(all(is.finite(tab$estimate)))
    expect_true(all(tab$std_error > 0))
    expect_true(all(tab$df > 0))
    expect_true(all(tab$fmi >= 0))
  }

  expect_output(print(pooled), "poolLongitudinalSub")
})
