# Shared fixture: pbc3 with an ordered-factor version of albumin (terciles),
# used as the "ordinal" biomarker throughout this file. Discretizing a real
# continuous biomarker (rather than simulating one from scratch) keeps the
# mixed continuous/ordinal model grounded in realistic covariate structure.
setup_copula_fixture <- function() {
  data(pbc3, envir = environment())
  d <- pbc3[pbc3$status3 == 1, ]
  d$albumin_cat <- cut(d$albumin,
                        breaks = stats::quantile(d$albumin, c(0, 1 / 3, 2 / 3, 1), na.rm = TRUE),
                        include.lowest = TRUE, labels = c("low", "mid", "high"),
                        ordered_result = TRUE)

  long_sub_fixed <- list("m1" = serBilir ~ year + age + sex, "m2" = albumin_cat ~ year + age + sex)
  long_sub_random <- list("m1" = ~ year | id, "m2" = ~ year | id)

  list(data_fit_all = d, long_sub_fixed = long_sub_fixed, long_sub_random = long_sub_random)
}

test_that("assert_biomarker_type validates length, type, and missingness", {
  expect_error(assert_biomarker_type(c("continuous"), 2), "length 2")
  expect_error(assert_biomarker_type(c(1, 2), 2), "length 2")
  expect_error(assert_biomarker_type(c("continuous", NA), 2), "length 2")
  expect_error(assert_biomarker_type(c("continuous", "nominal"), 2), "must only contain")
  expect_true(assert_biomarker_type(c("continuous", "ordinal"), 2))
})

test_that("resolve_biomarker_type auto-detects from factor/ordered response columns", {
  fx <- setup_copula_fixture()
  data_fit_all_norm <- list(fx$data_fit_all, fx$data_fit_all)

  detected <- resolve_biomarker_type(NULL, data_fit_all_norm, fx$long_sub_fixed, 2)
  expect_equal(detected, c("continuous", "ordinal"))
})

test_that("resolve_biomarker_type gives explicit biomarker_type priority over auto-detection", {
  fx <- setup_copula_fixture()
  data_fit_all_norm <- list(fx$data_fit_all, fx$data_fit_all)

  # data would auto-detect as c("continuous", "ordinal"); explicit says otherwise
  explicit <- resolve_biomarker_type(c("continuous", "continuous"), data_fit_all_norm, fx$long_sub_fixed, 2)
  expect_equal(explicit, c("continuous", "continuous"))
})

skip_if_ordinal <- function() testthat::skip_if_not_installed("ordinal")

test_that("all-continuous input is completely unaffected by biomarker_type (hard requirement)", {
  data(pbc3, envir = environment())
  data.fit <- pbc3[pbc3$status3 == 1, ]

  fit_dispatch <- longitudinalSub(data.fit, serBilir ~ year + age + sex, ~ year | id)
  fit_direct   <- longitudinalSubGaussian(data.fit, serBilir ~ year + age + sex, ~ year | id)
  expect_identical(fit_dispatch, fit_direct)
  # confirms the dispatcher doesn't tack on any new fields for the all-continuous path
  expect_equal(names(fit_dispatch), c("lfit", "Sigma_fit", "long_sub_fixed", "long_sub_random", "xlevels"))

  long_sub_fixed <- list("long1" = serBilir ~ year + age + sex, "long2" = albumin ~ year + age + sex)
  long_sub_random <- list("long1" = ~ year | id, "long2" = ~ year | id)
  fit_multi_dispatch <- longitudinalSub(list(data.fit, data.fit), long_sub_fixed, long_sub_random)
  fit_multi_direct   <- longitudinalSubGaussian(list(data.fit, data.fit), long_sub_fixed, long_sub_random)
  expect_identical(fit_multi_dispatch, fit_multi_direct)

  # explicit biomarker_type = all "continuous" must also route to the untouched path
  fit_explicit_cont <- longitudinalSub(data.fit, serBilir ~ year + age + sex, ~ year | id,
                                        biomarker_type = "continuous")
  expect_identical(fit_explicit_cont, fit_direct)
})

test_that("a mixed continuous/ordinal fit runs, auto-detects, and returns a sane joint model", {
  skip_if_ordinal()
  fx <- setup_copula_fixture()

  fit <- longitudinalSub(fx$data_fit_all, fx$long_sub_fixed, fx$long_sub_random)

  expect_s3_class(fit, "longitudinalSub.BJM")
  expect_equal(fit$biomarker_type, c("continuous", "ordinal"))
  expect_s3_class(fit$lfit[[1]], "lme")
  expect_s3_class(fit$lfit[[2]], "clmm")

  # random intercept + slope for each of the 2 markers -> 4x4 joint D
  expect_equal(dim(fit$Sigma_fit), c(4, 4))
  expect_true(isSymmetric(unname(fit$Sigma_fit)))
  expect_true(all(is.finite(fit$Sigma_fit)))
  # a valid covariance matrix must be positive semi-definite
  ev <- eigen(fit$Sigma_fit, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(ev > -1e-8))

  expect_length(fit$thresholds[[2]], 2) # 3 ordinal categories -> 2 thresholds
  expect_true(is.null(fit$thresholds[[1]])) # no thresholds for the continuous marker
})

test_that("explicit biomarker_type = auto-detected result agrees, and can override it", {
  skip_on_cran()
  skip_if_ordinal()
  fx <- setup_copula_fixture()

  fit_auto <- longitudinalSub(fx$data_fit_all, fx$long_sub_fixed, fx$long_sub_random)
  fit_explicit <- longitudinalSub(fx$data_fit_all, fx$long_sub_fixed, fx$long_sub_random,
                                   biomarker_type = c("continuous", "ordinal"))

  expect_equal(fit_auto$lfit[[1]]$coefficients$fixed, fit_explicit$lfit[[1]]$coefficients$fixed)
  expect_equal(fit_auto$lfit[[2]]$beta, fit_explicit$lfit[[2]]$beta)
  expect_equal(fit_auto$Sigma_fit, fit_explicit$Sigma_fit)

  # forcing the ordered-factor biomarker through the continuous (lme) path
  # is nonsensical and errors -- demonstrating the explicit override actually
  # took priority over auto-detection (which would have used "ordinal" here)
  expect_error(
    longitudinalSub(fx$data_fit_all, fx$long_sub_fixed, fx$long_sub_random,
                     biomarker_type = c("continuous", "continuous"))
  )
})

test_that("print and summary handle a mixed continuous/ordinal fit without erroring", {
  skip_on_cran()
  skip_if_ordinal()
  fx <- setup_copula_fixture()
  fit <- longitudinalSub(fx$data_fit_all, fx$long_sub_fixed, fx$long_sub_random)

  expect_output(print(fit), "continuous")
  expect_output(print(fit), "ordinal")
  expect_output(print(fit), "Cumulative thresholds")
  expect_output(summary(fit), "Correlations")
})

test_that("biomarker_type = 'ordinal' accepts a numeric response, ordered by value", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  f <- setup_copula_fixture()
  d <- f$data_fit_all
  d$albumin_score <- as.integer(d$albumin_cat) - 1L # 0/1/2
  as_factor <- suppressWarnings(longitudinalSub(d, f$long_sub_fixed, f$long_sub_random))
  as_number <- suppressWarnings(longitudinalSub(
    d, list(f$long_sub_fixed[[1]], albumin_score ~ year + age + sex), f$long_sub_random,
    biomarker_type = c("continuous", "ordinal")))

  expect_equal(as_number$lfit[[2]]$y.levels, c("0", "1", "2"))
  expect_equal(unname(as_number$thresholds[[2]]), unname(as_factor$thresholds[[2]]))
  expect_equal(as_number$lfit[[2]]$beta, as_factor$lfit[[2]]$beta)
  expect_equal(unname(as_number$Sigma_fit), unname(as_factor$Sigma_fit))
})

test_that("an unused category of an ordinal biomarker does not shift the others", {
  skip_on_cran()
  skip_if_ordinal()
  fx <- setup_copula_fixture()
  ref <- longitudinalSub(fx$data_fit_all, fx$long_sub_fixed, fx$long_sub_random)
  # clmm() drops a category nobody has; the EM used the factor's own codes,
  # which then no longer matched the thresholds (an error for an unused
  # lowest or middle level)
  for (lev in list(c("none", "low", "mid", "high"), c("low", "lowmid", "mid", "high"))) {
    d <- fx$data_fit_all
    d$albumin_cat <- factor(as.character(d$albumin_cat), levels = lev, ordered = TRUE)
    fit <- longitudinalSub(d, fx$long_sub_fixed, fx$long_sub_random)
    expect_equal(fit$thresholds[[2]], ref$thresholds[[2]])
    expect_equal(fit$Sigma_fit, ref$Sigma_fit)
  }
})

# a continuous marker and an ordinal one (or two ordinal ones), each with a
# random intercept and slope, correlated across markers
simulate_copula_markers <- function(seed, n = 400, visits = 0:5, second_ordinal = FALSE) {
  D <- matrix(c(1.0, 0.1, 0.5, 0.05,
                0.1, 0.2, 0.05, 0.08,
                0.5, 0.05, 1.0, 0.1,
                0.05, 0.08, 0.1, 0.2), 4)
  set.seed(seed)
  b <- matrix(stats::rnorm(4 * n), n) %*% chol(D)
  d <- expand.grid(year = visits, id = seq_len(n))
  d$x <- stats::rnorm(n)[d$id]
  lin1 <- 1 + 0.3 * d$year + 0.5 * d$x + b[d$id, 1] + b[d$id, 2] * d$year
  d$y <- if (second_ordinal) {
    cut(lin1 - 1 + stats::rnorm(nrow(d)), c(-Inf, 0, 1.5, Inf), labels = c("a", "b", "c"), ordered_result = TRUE)
  } else {
    lin1 + stats::rnorm(nrow(d), sd = 0.5)
  }
  z <- 0.2 * d$year + 0.5 * d$x + b[d$id, 3] + b[d$id, 4] * d$year + stats::rnorm(nrow(d))
  d$o <- cut(z, c(-Inf, -0.5, 1, Inf), labels = c("a", "b", "c"), ordered_result = TRUE)
  d
}

test_that("Sigma_fit recovers an ordinal biomarker's random-effects variances (quadrature EM)", {
  skip_on_cran()
  skip_if_ordinal()
  d <- simulate_copula_markers(1)
  fit <- longitudinalSub(d, list(y ~ year + x, o ~ year + x), list(~ year | id, ~ year | id))
  # truth 1, 0.2 (ordinal intercept, slope) and 0.5 (cross-marker intercept
  # covariance); plugging in the latent scores' means gave about 0.63 and 0.05
  expect_equal(fit$Sigma_fit[3, 3], 1, tolerance = 0.2)
  expect_equal(fit$Sigma_fit[4, 4], 0.2, tolerance = 0.3)
  expect_equal(fit$Sigma_fit[1, 3], 0.5, tolerance = 0.2)
  expect_equal(diag(fit$Sigma_fit)[1:2], c(1, 0.2), tolerance = 0.15, ignore_attr = TRUE)
})

test_that("many ordinal random effects fall back to a reproducible SAEM", {
  skip_on_cran()
  skip_if_ordinal()
  d <- simulate_copula_markers(1, second_ordinal = TRUE)
  set.seed(42)
  before <- .Random.seed
  fit <- longitudinalSub(d, list(y ~ year + x, o ~ year + x), list(~ year | id, ~ year | id))
  expect_identical(.Random.seed, before)   # the caller's random numbers are left alone
  again <- longitudinalSub(d, list(y ~ year + x, o ~ year + x), list(~ year | id, ~ year | id))
  expect_identical(fit$Sigma_fit, again$Sigma_fit)
  expect_equal(diag(fit$Sigma_fit), c(1, 0.2, 1, 0.2), tolerance = 0.3, ignore_attr = TRUE)
})

test_that("the quadrature and SAEM estimates of Sigma_fit agree", {
  skip_on_cran()
  skip_if_ordinal()
  d <- simulate_copula_markers(2, n = 300)
  fit <- longitudinalSub(d, list(y ~ year + x, o ~ year + x), list(~ year | id, ~ year | id))
  saem <- local({
    ns <- asNamespace("BJM")
    estep <- get("longitudinalSubVarCopula", ns)
    testthat::local_mocked_bindings(
      longitudinalSubVarCopula = function(...) estep(..., max_quadrature_dim = 0),
      .package = "BJM")
    longitudinalSub(d, list(y ~ year + x, o ~ year + x), list(~ year | id, ~ year | id))
  })
  expect_equal(fit$Sigma_fit, saem$Sigma_fit, tolerance = 0.15)
})

test_that("the copula EM does not depend on a biomarker's units", {
  skip_on_cran()
  skip_if_ordinal()
  d <- simulate_copula_markers(3, n = 200)
  fit <- longitudinalSub(d, list(y ~ year + x, o ~ year + x), list(~ year | id, ~ year | id))
  d$y <- d$y * 1000
  # an absolute convergence check never got below the tolerance at this scale
  expect_no_warning(scaled <- longitudinalSub(d, list(y ~ year + x, o ~ year + x), list(~ year | id, ~ year | id)))
  units <- c(1000, 1000, 1, 1)
  expect_equal(scaled$Sigma_fit / outer(units, units), fit$Sigma_fit, tolerance = 1e-3)
})
