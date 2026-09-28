# Tests for dynamicPredictionBioAll(): the batch counterpart to
# dynamicPredictionBio() that computes the bio_i-independent shared step
# (compute_bio_shared_step()) once and reuses it across every requested
# biomarker, instead of recomputing it once per biomarker the way calling
# dynamicPredictionBio() in a loop would (see R/dynamicPredictionBioAll.R).
#
# Three things are checked: (1) it returns exactly what looping
# dynamicPredictionBio() per biomarker would, so the refactor is behavior-
# preserving; (2) the shared step really is computed only once, not once per
# biomarker; (3) the bio_i = NULL default / explicit-ordinal-target error
# handling for a mixed continuous/ordinal (Gaussian-copula) fit.

setup_copula_bioall_fixture <- function() {
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

test_that("dynamicPredictionBioAll matches per-biomarker dynamicPredictionBio calls (all-continuous, competing risk)", {
  fx <- setup_dp_fixture()

  individual <- lapply(1:2, function(b) {
    dynamicPredictionBio(bio_i = b, fx$data_predict_all, fx$long_fit_all, fx$survival_fit_all,
                          prediction_time = 3, horizon = 3, time_variable = "year",
                          fx$survival_variable_all, fx$survival_trans_function,
                          bandcount2 = 10, bandcount3 = 15)
  })

  batch <- dynamicPredictionBioAll(bio_i = c(1, 2), fx$data_predict_all, fx$long_fit_all, fx$survival_fit_all,
                                    prediction_time = 3, horizon = 3, time_variable = "year",
                                    fx$survival_variable_all, fx$survival_trans_function,
                                    bandcount2 = 10, bandcount3 = 15)

  expect_s3_class(batch, "dynamicPredictionBioAll.BJM")
  expect_length(batch, 2)

  bio_names <- vapply(fx$long_fit_all$long_sub_fixed, function(f) as.character(formula(f)[[2]]), character(1))
  for (b in 1:2) {
    expect_s3_class(batch[[bio_names[b]]], "dynamicPredictionBio.BJM")
    expect_equal(batch[[bio_names[b]]]$Y_predict, individual[[b]]$Y_predict)
    expect_equal(batch[[bio_names[b]]]$Y_density, individual[[b]]$Y_density)
    expect_equal(batch[[bio_names[b]]]$Y_all, individual[[b]]$Y_all)
  }

  expect_equal(attr(batch, "bandcount2"), 10)
  expect_equal(unname(attr(batch, "bandcount3")), c(15, 15))
})

test_that("dynamicPredictionBioAll matches per-biomarker dynamicPredictionBio calls (all-continuous, no competing risk)", {
  fx <- setup_dp_fixture()
  data(pbc3, envir = environment())
  data_survival_fitting <- pbc3[!duplicated(pbc3$id), ]
  survival_fit_all <- survivalSub(data_survival_fitting, Surv(years, status3) ~ age + sex, NULL)

  individual <- lapply(1:2, function(b) {
    dynamicPredictionBio(bio_i = b, fx$data_predict_all, fx$long_fit_all, survival_fit_all,
                          prediction_time = 3, horizon = 3, time_variable = "year",
                          fx$survival_variable_all, fx$survival_trans_function,
                          bandcount2 = 10, bandcount3 = 15)
  })

  batch <- dynamicPredictionBioAll(bio_i = NULL, fx$data_predict_all, fx$long_fit_all, survival_fit_all,
                                    prediction_time = 3, horizon = 3, time_variable = "year",
                                    fx$survival_variable_all, fx$survival_trans_function,
                                    bandcount2 = 10, bandcount3 = 15)

  bio_names <- vapply(fx$long_fit_all$long_sub_fixed, function(f) as.character(formula(f)[[2]]), character(1))
  expect_setequal(names(batch), bio_names)
  for (b in 1:2) {
    expect_equal(batch[[bio_names[b]]]$Y_predict, individual[[b]]$Y_predict)
    expect_equal(batch[[bio_names[b]]]$Y_density, individual[[b]]$Y_density)
  }
})

test_that("dynamicPredictionBioAll computes the shared step only once across multiple biomarkers", {
  fx <- setup_dp_fixture()

  orig <- compute_bio_shared_step
  call_count <- 0
  testthat::local_mocked_bindings(
    compute_bio_shared_step = function(...) {
      call_count <<- call_count + 1
      orig(...)
    }
  )

  dynamicPredictionBioAll(bio_i = c(1, 2), fx$data_predict_all, fx$long_fit_all, fx$survival_fit_all,
                           prediction_time = 3, horizon = 3, time_variable = "year",
                           fx$survival_variable_all, fx$survival_trans_function,
                           bandcount2 = 10, bandcount3 = 15)

  expect_equal(call_count, 1)
})

test_that("dynamicPredictionBioAll with explicit numeric bandcounts resolves bandcount3 per marker under 'auto'", {
  fx <- setup_dp_fixture()

  orig <- compute_bio_shared_step
  call_count <- 0
  testthat::local_mocked_bindings(
    compute_bio_shared_step = function(...) {
      call_count <<- call_count + 1
      orig(...)
    }
  )

  # suppressWarnings(): bandcount3's doubling search may not fully converge
  # within max_rounds for this small fixture -- irrelevant to what this test
  # checks (the shared step's call count, and that per-marker bandcount3
  # values were resolved at all), and already covered as a warn-not-error
  # path by test-autoBandcount.R.
  batch <- suppressWarnings(
    dynamicPredictionBioAll(bio_i = c(1, 2), fx$data_predict_all, fx$long_fit_all, fx$survival_fit_all,
                             prediction_time = 3, horizon = 3, time_variable = "year",
                             fx$survival_variable_all, fx$survival_trans_function,
                             bandcount2 = 10, bandcount3 = "auto")
  )

  # bandcount2 is fixed (not "auto"), so the shared step is still built once,
  # even though bandcount3 is auto-tuned independently per marker.
  expect_equal(call_count, 1)
  expect_length(attr(batch, "bandcount3"), 2)
})

test_that("dynamicPredictionBioAll bio_i = NULL default includes ordinal biomarkers, for a mixed fit", {
  skip_if_not_installed("ordinal")
  fx <- setup_copula_bioall_fixture()

  data_survival_fitting <- fx$pbc3[!duplicated(fx$pbc3$id), ]
  survival_fit_all <- survivalSub(data_survival_fitting, Surv(years, status3) ~ age + sex, NULL)
  raw <- fx$pbc3[fx$pbc3$id == 2 & fx$pbc3$year <= 3, ]
  data_predict_all <- list(fx$to_cat(raw), fx$to_cat(raw))

  batch <- dynamicPredictionBioAll(bio_i = NULL, data_predict_all, fx$long_fit_all, survival_fit_all,
                                    prediction_time = 3, horizon = 3, time_variable = "year",
                                    survival_variable_all = list(), survival_trans_function = list(),
                                    bandcount2 = 10, bandcount3 = 15)

  expect_length(batch, 2)
  expect_setequal(names(batch), c("serBilir", "albumin_cat"))
  expect_s3_class(batch[["albumin_cat"]], "dynamicPredictionBio.BJM")
  expect_true(all(batch[["albumin_cat"]]$Y_predict %in% seq_len(3)))
  # bandcount3 is not applicable to an ordinal marker's fixed category grid.
  expect_true(is.na(unname(attr(batch, "bandcount3")["2"])))
})

test_that("dynamicPredictionBioAll on an explicit bio_i vector including an ordinal marker predicts it instead of erroring", {
  skip_if_not_installed("ordinal")
  fx <- setup_copula_bioall_fixture()

  data_survival_fitting <- fx$pbc3[!duplicated(fx$pbc3$id), ]
  survival_fit_all <- survivalSub(data_survival_fitting, Surv(years, status3) ~ age + sex, NULL)
  raw <- fx$pbc3[fx$pbc3$id == 2 & fx$pbc3$year <= 3, ]
  data_predict_all <- list(fx$to_cat(raw), fx$to_cat(raw))

  batch <- dynamicPredictionBioAll(bio_i = c(1, 2), data_predict_all, fx$long_fit_all, survival_fit_all,
                                    prediction_time = 3, horizon = 3, time_variable = "year",
                                    survival_variable_all = list(), survival_trans_function = list(),
                                    bandcount2 = 10, bandcount3 = 15)

  expect_length(batch, 2)
  expect_s3_class(batch[["albumin_cat"]], "dynamicPredictionBio.BJM")
  expect_true(all(batch[["albumin_cat"]]$Y_predict %in% seq_len(3)))
  expect_equal(attr(batch[["albumin_cat"]]$Y_all, "category_labels"), c("low", "mid", "high"))
})

test_that("dynamicPredictionBioAll for the continuous marker of a mixed fit matches dynamicPredictionBio", {
  skip_if_not_installed("ordinal")
  fx <- setup_copula_bioall_fixture()

  data_survival_fitting <- fx$pbc3[!duplicated(fx$pbc3$id), ]
  survival_fit_all <- survivalSub(data_survival_fitting, Surv(years, status3) ~ age + sex, NULL)
  raw <- fx$pbc3[fx$pbc3$id == 2 & fx$pbc3$year <= 3, ]
  data_predict_all <- list(fx$to_cat(raw), fx$to_cat(raw))

  set.seed(20260927)
  single <- dynamicPredictionBio(bio_i = 1, data_predict_all, fx$long_fit_all, survival_fit_all,
                                  prediction_time = 3, horizon = 3, time_variable = "year",
                                  survival_variable_all = list(), survival_trans_function = list(),
                                  bandcount2 = 10, bandcount3 = 15)

  set.seed(20260927)
  batch <- dynamicPredictionBioAll(bio_i = 1, data_predict_all, fx$long_fit_all, survival_fit_all,
                                    prediction_time = 3, horizon = 3, time_variable = "year",
                                    survival_variable_all = list(), survival_trans_function = list(),
                                    bandcount2 = 10, bandcount3 = 15)

  expect_equal(batch[["serBilir"]]$Y_predict, single$Y_predict)
  expect_equal(batch[["serBilir"]]$Y_density, single$Y_density)
})

test_that("dynamicPredictionBioAll for the ordinal marker of a mixed fit matches dynamicPredictionBio", {
  skip_if_not_installed("ordinal")
  fx <- setup_copula_bioall_fixture()

  data_survival_fitting <- fx$pbc3[!duplicated(fx$pbc3$id), ]
  survival_fit_all <- survivalSub(data_survival_fitting, Surv(years, status3) ~ age + sex, NULL)
  raw <- fx$pbc3[fx$pbc3$id == 2 & fx$pbc3$year <= 3, ]
  data_predict_all <- list(fx$to_cat(raw), fx$to_cat(raw))

  # bio_i = 2 alone (not c(1, 2)) so the sequence of mvtnorm::pmvnorm()
  # Monte-Carlo draws lines up identically between the single call and the
  # batch call's one-marker loop -- both do shared step, then marker 2,
  # with nothing else in between to consume the RNG differently.
  set.seed(20260927)
  single <- dynamicPredictionBio(bio_i = 2, data_predict_all, fx$long_fit_all, survival_fit_all,
                                  prediction_time = 3, horizon = 3, time_variable = "year",
                                  survival_variable_all = list(), survival_trans_function = list(),
                                  bandcount2 = 10, bandcount3 = 15)

  set.seed(20260927)
  batch <- dynamicPredictionBioAll(bio_i = 2, data_predict_all, fx$long_fit_all, survival_fit_all,
                                    prediction_time = 3, horizon = 3, time_variable = "year",
                                    survival_variable_all = list(), survival_trans_function = list(),
                                    bandcount2 = 10, bandcount3 = 15)

  expect_equal(batch[["albumin_cat"]]$Y_predict, single$Y_predict)
  expect_equal(batch[["albumin_cat"]]$Y_density, single$Y_density)
  expect_equal(batch[["albumin_cat"]]$Y_all, single$Y_all)
  expect_true(is.na(unname(attr(batch, "bandcount3")["2"])))
})
