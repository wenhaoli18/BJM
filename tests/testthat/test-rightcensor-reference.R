# See helper-rightcensor-reference.R: right-censored outputs must be
# unchanged by the interval-censoring refactor.
test_that("right-censored outputs match the pre-interval-censoring reference", {
  skip_on_cran()
  baseline <- readRDS(test_path("testdata", "baseline_rightcensor.rds"))
  current <- right_censor_reference_outputs()
  expect_equal(current, baseline, tolerance = 1e-12)
})
