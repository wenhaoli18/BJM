test_that("longitudinalSub fits a univariate model when given a single formula", {
  data(pbc3, envir = environment())
  data.fit <- pbc3[pbc3$status3 == 1, ]

  fit <- longitudinalSub(data.fit,
                          serBilir ~ year + age + sex,
                          ~ year | id)

  expect_s3_class(fit, "longitudinalSub.BJM")
  expect_length(fit, 5)
  expect_length(fit[[1]], 1)
  expect_s3_class(fit[[1]][[1]], "lme")
  # D is the 2x2 random-effects covariance matrix (intercept + slope)
  expect_equal(dim(fit[[2]]), c(2, 2))
})

test_that("longitudinalSub fits a multivariate model and returns a joint D matrix", {
  data(pbc3, envir = environment())
  data.fit <- pbc3[pbc3$status3 == 1, ]

  long_sub_fixed <- list(
    "long1" = serBilir ~ year + age + sex,
    "long2" = albumin ~ year + age + sex
  )
  long_sub_random <- list("long1" = ~ year | id, "long2" = ~ year | id)

  fit <- longitudinalSub(list(data.fit, data.fit), long_sub_fixed, long_sub_random)

  expect_s3_class(fit, "longitudinalSub.BJM")
  expect_length(fit[[1]], 2)
  # random intercept + slope for each of the 2 markers -> 4x4 joint D
  expect_equal(dim(fit[[2]]), c(4, 4))
  expect_true(isSymmetric(unname(fit[[2]])))
})

test_that("print and summary methods run without error", {
  data(pbc3, envir = environment())
  data.fit <- pbc3[pbc3$status3 == 1, ]
  fit <- longitudinalSub(data.fit, serBilir ~ year + age + sex, ~ year | id)

  expect_output(print(fit), "Longitudinal Sub-model")
  expect_output(summary(fit), "Correlations")
})
