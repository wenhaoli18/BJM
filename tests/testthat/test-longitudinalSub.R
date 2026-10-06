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

test_that("longitudinalSub accepts numeric and character patient ids", {
  data(pbc3, envir = environment())
  data.fit <- pbc3[pbc3$status3 == 1, ]
  long_sub_fixed <- list(serBilir ~ year + age + sex, albumin ~ year + age + sex)
  long_sub_random <- list(~ year | id, ~ year | id)
  ref <- longitudinalSub(data.fit, long_sub_fixed, long_sub_random)

  for (id_type in c("numeric", "character")) {
    d <- data.fit
    d$id <- as(as.character(d$id), id_type)
    fit <- longitudinalSub(d, long_sub_fixed, long_sub_random)
    expect_equal(fit$Sigma_fit, ref$Sigma_fit, ignore_attr = TRUE)
  }
})

test_that("biomarkers may have different random-effects structures", {
  data(pbc3, envir = environment())
  data.fit <- pbc3[pbc3$status3 == 1, ]
  fit <- longitudinalSub(data.fit,
                         list(serBilir ~ year + age + sex, albumin ~ year + age + sex),
                         list(~ 1 | id, ~ year | id))
  # random intercept for marker 1 + intercept and slope for marker 2
  expect_equal(dim(fit$Sigma_fit), c(3, 3))
})

test_that("Sigma_fit averages over the subjects actually retained in the fit", {
  # 50 subjects never have albumin, so they are excluded from the joint fit.
  # At EM convergence Sigma_fit must equal the mean, over the *retained*
  # subjects, of each subject's posterior E[b b^T] -- computed here by hand
  # from the returned fit. It used to divide by all 169 subjects in the
  # first biomarker's data instead, underestimating Sigma_fit several-fold.
  data(pbc3, envir = environment())
  d <- pbc3[pbc3$status3 == 1, ]
  set.seed(1)
  never <- sample(unique(as.character(d$id)), 50)
  d$albumin[d$id %in% never] <- NA
  long_sub_fixed <- list(serBilir ~ year + age + sex, albumin ~ year + age + sex)
  fit <- longitudinalSub(d, long_sub_fixed, list(~ year | id, ~ year | id))

  D <- fit$Sigma_fit
  beta <- lapply(fit$lfit, nlme::fixef)
  sigma2 <- vapply(fit$lfit, function(u) u$sigma^2, numeric(1))
  retained <- setdiff(unique(as.character(d$id)), never)
  expect_length(retained, 119)

  EbbT <- lapply(retained, function(i) {
    rows <- lapply(long_sub_fixed, function(f) {
      di <- d[d$id == i, ]
      di[stats::complete.cases(di[all.vars(f)]), ]
    })
    y <- unlist(Map(function(f, di) di[[all.vars(f)[1]]], long_sub_fixed, rows))
    X <- as.matrix(Matrix::bdiag(Map(function(f, di) model.matrix(f, di), long_sub_fixed, rows)))
    Z <- as.matrix(Matrix::bdiag(lapply(rows, function(di) cbind(1, di$year))))
    S_inv <- diag(rep(1 / sigma2, vapply(rows, nrow, numeric(1))))
    A <- solve(t(Z) %*% S_inv %*% Z + solve(D))
    Eb <- A %*% t(Z) %*% S_inv %*% (y - X %*% unlist(beta))
    A + Eb %*% t(Eb)
  })
  expect_equal(unname(D), unname(Reduce(`+`, EbbT) / length(retained)), tolerance = 1e-3)
})

test_that("the joint EM stops at max.iter, and a singular covariance is reported", {
  data(pbc3, envir = environment())
  d <- pbc3[pbc3$status3 == 1, ]
  seen <- new.env()
  orig <- BJM:::longitudinalSubVar
  testthat::local_mocked_bindings(longitudinalSubVar = function(thetaLong, l, ...) {
    seen$theta <- thetaLong
    seen$l <- l
    orig(thetaLong, l, ...)
  })
  longitudinalSub(d, list(serBilir ~ year, albumin ~ year), list(~ year | id, ~ year | id))

  expect_warning(orig(seen$theta, seen$l, tol.em = 1e-12, verbose = FALSE, max.iter = 2),
                 "did not converge within 2 iterations")
  singular <- seen$theta
  singular$D[] <- 0
  msgs <- character()
  withCallingHandlers(orig(singular, seen$l, tol.em = 1e-4, verbose = FALSE, max.iter = 3),
                      warning = function(w) {
                        msgs <<- c(msgs, conditionMessage(w))
                        invokeRestart("muffleWarning")
                      })
  expect_true(any(grepl("became singular", msgs)))
})

test_that("a missing value in a variable used only in long_sub_random drops that row", {
  data(pbc3, envir = environment())
  d <- pbc3[pbc3$status3 == 1, ]
  d$visit_time <- d$year
  d$visit_time[c(5, 50, 300)] <- NA
  # this used to fail with "arguments imply differing number of rows"
  with_na <- longitudinalSub(d, serBilir ~ year + age + years, ~ visit_time | id)
  dropped <- longitudinalSub(d[!is.na(d$visit_time), ], serBilir ~ year + age + years, ~ visit_time | id)
  expect_equal(with_na$Sigma_fit, dropped$Sigma_fit)
  expect_equal(nlme::fixef(with_na$lfit[[1]]), nlme::fixef(dropped$lfit[[1]]))
})
