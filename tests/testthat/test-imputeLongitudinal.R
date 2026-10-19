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
  skip_on_cran()
  skip_if_no_torch()
  f <- setup_impute_fixture()

  n_missing_pre <- sum(is.na(f$data_fit_all$serBilir)) + sum(is.na(f$data_fit_all$albumin))
  expect_gt(n_missing_pre, 0)

  imputed <- imputeLongitudinal(f$data_fit_all, f$long_sub_fixed, f$long_sub_random,
                                 f$time_variable, n_imputations = 2, latent_dim = 4,
                                 hidden_units = c(16, 8), epochs = 15,
                                 importance_samples = 5, method = "miwae", seed = 1)

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
  skip_if_no_torch()
  f <- setup_impute_fixture()

  imputed <- imputeLongitudinal(f$data_fit_all, f$long_sub_fixed, f$long_sub_random,
                                 f$time_variable, n_imputations = 3, latent_dim = 4,
                                 hidden_units = c(16, 8), epochs = 15,
                                 importance_samples = 5, impute = "multiple", method = "miwae",
                                 seed = 2)

  expect_length(imputed$data_fit_all_list, 3)
  for (completion in imputed$data_fit_all_list) {
    expect_equal(nrow(completion[[1]]), nrow(f$data_fit_all))
    expect_false(anyNA(completion[[1]]$serBilir))
    expect_false(anyNA(completion[[2]]$albumin))
  }
})

test_that("imputeLongitudinal(method = 'diffusion') fills all missing biomarker cells, preserves observed values exactly, and preserves row/subject counts (requires torch)", {
  skip_if_no_torch()
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
  skip_if_no_torch()
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

test_that("imputeLongitudinal(seed =) is reproducible and leaves the caller's RNG stream alone (requires torch)", {
  skip_if_no_torch()
  f <- setup_impute_fixture()
  impute <- function() {
    imputeLongitudinal(f$data_fit_all, f$long_sub_fixed, f$long_sub_random,
                       f$time_variable, n_imputations = 2, latent_dim = 4,
                       hidden_units = c(16, 8), epochs = 10,
                       importance_samples = 5, impute = "multiple", method = "miwae", seed = 3)
  }
  # MIWAE's importance resampling draws with R's RNG; seeding torch alone
  # used to leave these completions different from run to run.
  set.seed(10)
  first <- impute()
  after_first <- runif(1)
  set.seed(20)
  second <- impute()
  expect_equal(second$data_fit_all_list, first$data_fit_all_list)

  set.seed(10)
  expect_equal(runif(1), after_first)
})

test_that("the imputation model uses the visit time when long_sub_fixed does", {
  # the visit time is also the key that aligns rows across biomarkers; it used
  # to be dropped from the model along with the id
  expect_equal(imputation_covariates(list(c("year", "age", "years"), c("year", "sex")), "id"),
               c("year", "age", "years", "sex"))
})

test_that("imputed values follow a biomarker's trend over visit time (requires torch)", {
  skip_on_cran()
  skip_if_no_torch()
  set.seed(1)
  d <- expand.grid(year = 0:5, id = 1:40)
  d$age <- stats::rnorm(40)[d$id]
  d$y <- 2 * d$year + 0.5 * d$age + stats::rnorm(nrow(d), sd = 0.3)
  missing <- sample(nrow(d), 40)
  d$y[missing] <- NA
  out <- imputeLongitudinal(d, y ~ year + age, ~ 1 | id, time_variable = "year", epochs = 300,
                            method = "miwae", seed = 1)
  filled <- out$data_fit_all[[1]]
  filled <- filled[order(filled$id, filled$year), ]
  ref <- d[order(d$id, d$year), ]
  m <- is.na(ref$y)
  # without the visit time: correlation -0.24 with it, RMSE 4.1 against the truth
  expect_gt(stats::cor(filled$y[m], ref$year[m]), 0.9)
  expect_lt(sqrt(mean((filled$y[m] - (2 * ref$year[m] + 0.5 * ref$age[m]))^2)), 1)
})

test_that("imputeLongitudinal rejects a non-numeric (e.g. ordinal) biomarker", {
  f <- setup_impute_fixture()
  d <- f$data_fit_all
  d$albumin <- cut(d$albumin, 3, ordered_result = TRUE)
  # it used to be coerced to numeric, turning the whole column into NaN
  expect_error(imputeLongitudinal(d, f$long_sub_fixed, f$long_sub_random, time_variable = "year"),
               "`albumin` is not numeric")
})

test_that("imputeLongitudinal rejects repeated (id, time) rows", {
  f <- setup_impute_fixture()
  d <- rbind(f$data_fit_all, f$data_fit_all[1, ])
  expect_error(imputeLongitudinal(d, f$long_sub_fixed, f$long_sub_random, time_variable = "year"),
               "repeating an \\(id, year\\) pair")
})

test_that("the diffusion backend imputes from the row's other biomarkers (requires torch)", {
  skip_on_cran()
  skip_if_no_torch()
  set.seed(2)
  n <- 150
  d <- expand.grid(year = 0:5, id = seq_len(n))
  d$age <- stats::rnorm(n)[d$id]
  d$y1 <- 1 + 0.5 * d$year + 0.4 * d$age + stats::rnorm(n, 0, 0.7)[d$id] + stats::rnorm(nrow(d), 0, 0.4)
  d$y2 <- 0.8 * d$y1 + 0.3 * d$age + stats::rnorm(nrow(d), 0, 0.3)
  truth <- d$y2
  m <- stats::runif(nrow(d)) < stats::plogis(-1.5 + (d$y1 - mean(d$y1)))   # MAR given y1
  d$y2[m] <- NA
  out <- imputeLongitudinal(d, list(y1 ~ year + age, y2 ~ year + age), list(~ 1 | id, ~ 1 | id),
                            time_variable = "year", method = "diffusion", seed = 1)
  filled <- out$data_fit_all[[2]]
  filled <- filled[order(filled$id, filled$year), ]
  ref <- d[order(d$id, d$year), ]
  # residual SD given y1 is 0.3; the earlier sampler, which barely used y1,
  # gave an RMSE of about 0.55
  expect_lt(sqrt(mean((filled$y2[is.na(ref$y2)] - truth[order(d$id, d$year)][is.na(ref$y2)])^2)), 0.42)
})

test_that("method = 'mixed' (the default) fills every missing cell from the mixed model, without torch", {
  f <- setup_impute_fixture()
  imputed <- imputeLongitudinal(f$data_fit_all, f$long_sub_fixed, f$long_sub_random, f$time_variable,
                                n_imputations = 3, impute = "multiple", seed = 1)
  expect_length(imputed$data_fit_all_list, 3)
  observed <- !is.na(f$data_fit_all$serBilir)
  for (completion in imputed$data_fit_all_list) {
    expect_false(anyNA(completion[[1]]$serBilir))
    expect_false(anyNA(completion[[2]]$albumin))
    key <- match(paste(f$data_fit_all$id, f$data_fit_all$year), paste(completion[[1]]$id, completion[[1]]$year))
    expect_equal(completion[[1]]$serBilir[key][observed], f$data_fit_all$serBilir[observed])
  }
  # draws differ between completions, and seed makes them reproducible
  expect_false(isTRUE(all.equal(imputed$data_fit_all_list[[1]], imputed$data_fit_all_list[[2]])))
  again <- imputeLongitudinal(f$data_fit_all, f$long_sub_fixed, f$long_sub_random, f$time_variable,
                              n_imputations = 3, impute = "multiple", seed = 1)
  expect_equal(again$data_fit_all_list, imputed$data_fit_all_list)
})

test_that("mixed-model imputation follows each subject's own level", {
  # subjects differ in level (random intercept SD 2); a deleted value should
  # come back near the subject's own other values, not the population mean
  set.seed(4)
  d <- expand.grid(year = 0:5, id = 1:60)
  level <- stats::rnorm(60, 0, 2)
  d$y <- level[d$id] + 0.5 * d$year + stats::rnorm(nrow(d), sd = 0.3)
  truth <- d$y
  missing <- sample(nrow(d), 60)
  d$y[missing] <- NA
  out <- imputeLongitudinal(d, y ~ year, ~ 1 | id, time_variable = "year", seed = 1)
  filled <- out$data_fit_all[[1]]
  key <- match(paste(d$id, d$year)[missing], paste(filled$id, filled$year))
  error <- filled$y[key] - truth[missing]
  # the subject's level explains most of the variance; the error is close to
  # the residual SD (0.3), far below the spread of the levels (2)
  expect_lt(sqrt(mean(error^2)), 0.6)
})

test_that("history_features summarises a subject's other visits, never the visit itself", {
  aligned <- data.frame(id = c(1, 1, 1, 2), year = c(0, 1, 3, 0))
  resp <- cbind(y = c(10, NA, 14, 5))
  mask <- !is.na(resp)
  h <- history_features(aligned, "id", "year", resp, mask)
  expect_equal(unname(h[1, ]), c(14, log1p(1), 14, log1p(3)))     # visit 1: only visit 3 is "other"
  expect_equal(unname(h[2, ]), c(12, log1p(2), 10, log1p(1)))     # nearest other: year 0
  # subject 2 has no other visit: the overall mean, a count of 0, the largest distance
  expect_equal(unname(h[4, 1:2]), c(mean(c(10, 14, 5)), 0))
  expect_error(imputeLongitudinal(aligned, y ~ year, ~ 1 | id, "year", history = NA), "history")
})
