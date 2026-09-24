test_that("survivalSub fits a marginal Cox model without competing risks", {
  data(pbc3, envir = environment())
  data_survival_fitting <- pbc3[!duplicated(pbc3$id), ]

  fit <- survivalSub(data_survival_fitting, Surv(years, status3) ~ age + sex, NULL)

  expect_s3_class(fit, "survivalSub.BJM")
  expect_length(fit, 4)
  expect_s3_class(fit[[1]], "coxph")
  expect_null(fit[[3]])
  expect_null(fit[[4]])
  expect_equal(names(coef(fit[[1]])), c("age", "sex"))
})

test_that("survivalSub fits the competing-risks logistic sub-model when requested", {
  data(pbc3, envir = environment())
  data_survival_fitting <- pbc3[!duplicated(pbc3$id), ]

  fit <- survivalSub(data_survival_fitting,
                      Surv(years, status3) ~ age + sex,
                      status4 ~ years + age + sex)

  expect_s3_class(fit, "survivalSub.BJM")
  expect_s3_class(fit[[3]], "glm")
  expect_identical(family(fit[[3]])$family, "binomial")
  # the GLM should only be fit on subjects who experienced the primary event
  expect_equal(nrow(fit[[3]]$data), sum(data_survival_fitting$status3 != 0))
})

test_that("print and summary methods run without error", {
  data(pbc3, envir = environment())
  data_survival_fitting <- pbc3[!duplicated(pbc3$id), ]
  fit <- survivalSub(data_survival_fitting, Surv(years, status3) ~ age + sex,
                      status4 ~ years + age + sex)

  expect_output(print(fit), "Cox PH")
  expect_output(summary(fit), "Baseline cumulative hazard")
})
