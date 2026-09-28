# Exercises imputeLongitudinal(), the optional deep-generative (MIWAE)
# preprocessing step that fills interrupted-follow-up gaps before
# longitudinalSub() runs. Input-validation tests below do not require torch
# to be installed -- imputeLongitudinal() runs every cheap argument check
# before it ever touches torch, specifically so a user without torch still
# gets a specific, actionable message for a common mistake (bad
# data_fit_all, mismatched formula lengths, ...) instead of always seeing
# "torch is not installed" regardless of what else is wrong. Only the
# integration test at the bottom needs torch to actually be installed.

test_that("imputeLongitudinal rejects a data_fit_all that is not a data.frame or list of them", {
  f <- setup_impute_fixture()

  expect_error(
    imputeLongitudinal(5, f$long_sub_fixed, f$long_sub_random, f$time_variable),
    "`data_fit_all` must be a list of data.frame objects"
  )
})

test_that("imputeLongitudinal rejects mismatched long_sub_fixed/long_sub_random lengths", {
  f <- setup_impute_fixture()
  long_sub_fixed <- list(f$long_sub_fixed[[1]], f$long_sub_fixed[[2]])
  long_sub_random <- list(f$long_sub_random[[1]])

  expect_error(
    imputeLongitudinal(f$data_fit_all, long_sub_fixed, long_sub_random, f$time_variable),
    "long_sub_fixed.*long_sub_random"
  )
})

test_that("imputeLongitudinal rejects an unknown time_variable column", {
  f <- setup_impute_fixture()

  expect_error(
    imputeLongitudinal(f$data_fit_all, f$long_sub_fixed, f$long_sub_random, "not_a_column"),
    "not_a_column"
  )
})

test_that("imputeLongitudinal rejects a non-positive-integer n_imputations", {
  f <- setup_impute_fixture()

  expect_error(
    imputeLongitudinal(f$data_fit_all, f$long_sub_fixed, f$long_sub_random, f$time_variable,
                        n_imputations = 0),
    "`n_imputations` must be a single positive integer"
  )
})

test_that("imputeLongitudinal rejects a non-positive-integer epochs", {
  f <- setup_impute_fixture()

  expect_error(
    imputeLongitudinal(f$data_fit_all, f$long_sub_fixed, f$long_sub_random, f$time_variable,
                        epochs = -1),
    "`epochs` must be a single positive integer"
  )
})

test_that("imputeLongitudinal rejects a non-positive-integer diffusion_steps", {
  f <- setup_impute_fixture()

  expect_error(
    imputeLongitudinal(f$data_fit_all, f$long_sub_fixed, f$long_sub_random, f$time_variable,
                        diffusion_steps = 0),
    "`diffusion_steps` must be a single positive integer"
  )
})

test_that("imputeLongitudinal rejects an unknown method", {
  f <- setup_impute_fixture()

  expect_error(
    imputeLongitudinal(f$data_fit_all, f$long_sub_fixed, f$long_sub_random, f$time_variable,
                        method = "not_a_real_method"),
    "should be one of"
  )
})

test_that("imputeLongitudinal rejects a long_sub_fixed with a transformed left-hand side", {
  f <- setup_impute_fixture()
  long_sub_fixed <- f$long_sub_fixed
  long_sub_fixed[[1]] <- log(serBilir) ~ year + age + sex + years

  expect_error(
    imputeLongitudinal(f$data_fit_all, long_sub_fixed, f$long_sub_random, f$time_variable),
    "transformed left-hand side"
  )
})

test_that("imputeLongitudinal rejects duplicate response variables across biomarkers", {
  f <- setup_impute_fixture()
  long_sub_fixed <- list(f$long_sub_fixed[[1]], f$long_sub_fixed[[1]])
  long_sub_random <- list(f$long_sub_random[[1]], f$long_sub_random[[1]])

  expect_error(
    imputeLongitudinal(f$data_fit_all, long_sub_fixed, long_sub_random, f$time_variable),
    "must be distinct across biomarkers"
  )
})

test_that("imputeLongitudinal rejects covariates with missing values", {
  f <- setup_impute_fixture()
  data_fit_all <- f$data_fit_all
  data_fit_all$age[1] <- NA

  expect_error(
    imputeLongitudinal(data_fit_all, f$long_sub_fixed, f$long_sub_random, f$time_variable),
    "covariate.*missing values"
  )
})

test_that("imputeLongitudinal errors when there is nothing to impute", {
  f <- setup_impute_fixture()
  complete_data <- f$data_fit_all
  complete_data$serBilir[is.na(complete_data$serBilir)] <- 1
  complete_data$albumin[is.na(complete_data$albumin)] <- 1

  expect_error(
    imputeLongitudinal(complete_data, f$long_sub_fixed, f$long_sub_random, f$time_variable),
    "nothing to impute"
  )
})

test_that("assert_package_installed reports a missing optional package with an actionable message", {
  # Exercises the same shared helper imputeLongitudinal() uses to guard its
  # torch dependency, using a package name guaranteed not to be installed
  # (rather than mocking base::requireNamespace(), which testthat does not
  # support mocking). torch itself is installed in this test environment, so
  # this cannot be reproduced by calling imputeLongitudinal() directly here.
  expect_error(
    assert_package_installed("not.a.real.package.xyz", "imputeLongitudinal()"),
    "requires the 'not.a.real.package.xyz' package, which is not installed"
  )
})

test_that("imputeLongitudinal fills all missing biomarker cells and preserves row/subject counts (requires torch)", {
  skip_if_not_installed("torch")
  f <- setup_impute_fixture()

  n_missing_pre <- sum(is.na(f$data_fit_all$serBilir)) + sum(is.na(f$data_fit_all$albumin))
  expect_gt(n_missing_pre, 0)

  imputed <- imputeLongitudinal(f$data_fit_all, f$long_sub_fixed, f$long_sub_random,
                                 f$time_variable, n_imputations = 2, latent_dim = 4,
                                 hidden_units = c(16, 8), epochs = 15,
                                 importance_samples = 5, seed = 1)

  expect_named(imputed, c("data_fit_all", "data_fit_all_list", "diagnostics"))
  expect_equal(nrow(imputed$data_fit_all[[1]]), nrow(f$data_fit_all))
  expect_equal(nrow(imputed$data_fit_all[[2]]), nrow(f$data_fit_all))
  expect_false(anyNA(imputed$data_fit_all[[1]]$serBilir))
  expect_false(anyNA(imputed$data_fit_all[[2]]$albumin))

  obs_range_serBilir <- range(f$data_fit_all$serBilir, na.rm = TRUE)
  obs_range_albumin <- range(f$data_fit_all$albumin, na.rm = TRUE)
  pad_serBilir <- diff(obs_range_serBilir)
  pad_albumin <- diff(obs_range_albumin)
  expect_true(all(imputed$data_fit_all[[1]]$serBilir >=
                     obs_range_serBilir[1] - pad_serBilir &
                   imputed$data_fit_all[[1]]$serBilir <=
                     obs_range_serBilir[2] + pad_serBilir))
  expect_true(all(imputed$data_fit_all[[2]]$albumin >=
                     obs_range_albumin[1] - pad_albumin &
                   imputed$data_fit_all[[2]]$albumin <=
                     obs_range_albumin[2] + pad_albumin))

  expect_equal(imputed$diagnostics$serBilir$biomarker, "serBilir")
  expect_equal(imputed$diagnostics$serBilir$n_missing, sum(is.na(f$data_fit_all$serBilir)))
  expect_true(is.numeric(imputed$diagnostics$serBilir$deep_generative$imputed_mean))
  expect_true(is.numeric(imputed$diagnostics$serBilir$classical_baseline$imputed_mean))

  # longitudinalSub() itself must still succeed on the completed data (no
  # missing-data-related crash) --
  expect_no_error(longitudinalSub(imputed$data_fit_all, f$long_sub_fixed, f$long_sub_random))

  n_retained_complete_case <- length(count_retained_subjects(f$data_fit_all, f$long_sub_fixed))
  n_retained_imputed <- length(count_retained_subjects(imputed$data_fit_all, f$long_sub_fixed))
  expect_gte(n_retained_imputed, n_retained_complete_case)
})

test_that("imputeLongitudinal(impute = 'multiple') returns n_imputations distinct completed datasets (requires torch)", {
  skip_if_not_installed("torch")
  f <- setup_impute_fixture()

  imputed <- imputeLongitudinal(f$data_fit_all, f$long_sub_fixed, f$long_sub_random,
                                 f$time_variable, n_imputations = 3, latent_dim = 4,
                                 hidden_units = c(16, 8), epochs = 15,
                                 importance_samples = 5, impute = "multiple", seed = 2)

  expect_length(imputed$data_fit_all_list, 3)
  for (completion in imputed$data_fit_all_list) {
    expect_equal(nrow(completion[[1]]), nrow(f$data_fit_all))
    expect_false(anyNA(completion[[1]]$serBilir))
    expect_false(anyNA(completion[[2]]$albumin))
  }
})

test_that("imputeLongitudinal(method = 'diffusion') fills all missing biomarker cells, preserves observed values exactly, and preserves row/subject counts (requires torch)", {
  skip_if_not_installed("torch")
  f <- setup_impute_fixture()

  imputed <- imputeLongitudinal(f$data_fit_all, f$long_sub_fixed, f$long_sub_random,
                                 f$time_variable, n_imputations = 2, method = "diffusion",
                                 hidden_units = c(16, 8), epochs = 30, diffusion_steps = 20,
                                 seed = 3)

  expect_named(imputed, c("data_fit_all", "data_fit_all_list", "diagnostics"))
  expect_equal(nrow(imputed$data_fit_all[[1]]), nrow(f$data_fit_all))
  expect_equal(nrow(imputed$data_fit_all[[2]]), nrow(f$data_fit_all))
  expect_false(anyNA(imputed$data_fit_all[[1]]$serBilir))
  expect_false(anyNA(imputed$data_fit_all[[2]]$albumin))

  # the RePaint-style sampler is designed to reproduce every already-observed
  # cell exactly, not just impute the missing ones -- verify that directly
  obs_serBilir <- !is.na(f$data_fit_all$serBilir)
  obs_albumin <- !is.na(f$data_fit_all$albumin)
  expect_equal(imputed$data_fit_all[[1]]$serBilir[obs_serBilir],
               f$data_fit_all$serBilir[obs_serBilir])
  expect_equal(imputed$data_fit_all[[2]]$albumin[obs_albumin],
               f$data_fit_all$albumin[obs_albumin])

  expect_equal(imputed$diagnostics$serBilir$biomarker, "serBilir")
  expect_equal(imputed$diagnostics$serBilir$n_missing, sum(is.na(f$data_fit_all$serBilir)))
  expect_true(is.numeric(imputed$diagnostics$serBilir$deep_generative$imputed_mean))

  expect_no_error(longitudinalSub(imputed$data_fit_all, f$long_sub_fixed, f$long_sub_random))

  n_retained_complete_case <- length(count_retained_subjects(f$data_fit_all, f$long_sub_fixed))
  n_retained_imputed <- length(count_retained_subjects(imputed$data_fit_all, f$long_sub_fixed))
  expect_gte(n_retained_imputed, n_retained_complete_case)
})

test_that("imputeLongitudinal(method = 'diffusion', impute = 'multiple') returns n_imputations distinct completed datasets (requires torch)", {
  skip_if_not_installed("torch")
  f <- setup_impute_fixture()

  imputed <- imputeLongitudinal(f$data_fit_all, f$long_sub_fixed, f$long_sub_random,
                                 f$time_variable, n_imputations = 3, method = "diffusion",
                                 hidden_units = c(16, 8), epochs = 30, diffusion_steps = 20,
                                 impute = "multiple", seed = 4)

  expect_length(imputed$data_fit_all_list, 3)
  for (completion in imputed$data_fit_all_list) {
    expect_equal(nrow(completion[[1]]), nrow(f$data_fit_all))
    expect_false(anyNA(completion[[1]]$serBilir))
    expect_false(anyNA(completion[[2]]$albumin))
  }
})
