# See helper-rightcensor-reference.R: right-censored outputs must be
# unchanged by the interval-censoring refactor.
test_that("right-censored outputs match the pre-interval-censoring reference", {
  skip_on_cran()
  baseline <- readRDS(test_path("testdata", "baseline_rightcensor.rds"))
  current <- right_censor_reference_outputs()
  # 1e-8 rather than exact: the reference was saved on one platform, and
  # BLAS/compiler differences elsewhere move the last few digits
  expect_equal(current, baseline, tolerance = 1e-8)
})
