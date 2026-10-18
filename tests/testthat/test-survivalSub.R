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

test_that("the fitted coxph model does not depend on the caller's data variable being named data_survival_fitting", {
  data(pbc3, envir = environment())
  # Deliberately not named `data_survival_fitting` (survivalSub()'s own
  # parameter name), to catch model.frame.coxph() re-evaluating the fit's
  # captured call in the formula's environment, where this object would not
  # exist under any other name.
  surv_train_data <- pbc3[!duplicated(pbc3$id), ]

  fit <- survivalSub(surv_train_data, Surv(years, status3) ~ age + sex, NULL)

  expect_no_error(survival::basehaz(fit$coxph_fit, centered = FALSE))
})

test_that("print and summary methods run without error", {
  data(pbc3, envir = environment())
  data_survival_fitting <- pbc3[!duplicated(pbc3$id), ]
  fit <- survivalSub(data_survival_fitting, Surv(years, status3) ~ age + sex,
                      status4 ~ years + age + sex)

  expect_output(print(fit), "Cox PH")
  expect_output(summary(fit), "Baseline cumulative hazard")
})

test_that("the survival time may only enter the event-type model as a main effect", {
  data(pbc3, envir = environment())
  d <- pbc3[!duplicated(pbc3$id), ]
  expect_error(survivalSub(d, Surv(years, status3) ~ age + sex, status4 ~ log(years) + age),
               "`log\\(years\\)`")
  expect_error(survivalSub(d, Surv(years, status3) ~ age + sex, status4 ~ years * age),
               "`years:age`")
  expect_no_error(survivalSub(d, Surv(years, status3) ~ age + sex, status4 ~ years + age))
  expect_no_error(survivalSub(d, Surv(years, status3) ~ age + sex, status4 ~ age + sex))
})

test_that("the event rows come from the Surv() status, whatever its form", {
  data(pbc3, envir = environment())
  d <- pbc3[!duplicated(pbc3$id), ]
  ref <- survivalSub(d, Surv(years, status3) ~ age + sex, status4 ~ years + age + sex)

  # an expression as the status (previously "undefined columns selected")
  expr <- survivalSub(d, Surv(years, status3 == 1) ~ age + sex, status4 ~ years + age + sex)
  expect_equal(coef(expr$glm_fit), coef(ref$glm_fit))

  # 1 = censored / 2 = event coding (previously every subject entered the glm)
  d$status12 <- d$status3 + 1
  coded <- survivalSub(d, Surv(years, status12) ~ age + sex, status4 ~ years + age + sex)
  expect_equal(coef(coded$glm_fit), coef(ref$glm_fit))
})

test_that("print() counts cause 1 as the event type predictRisk() calls risk_prob_1", {
  data(pbc3, envir = environment())
  first <- pbc3[!duplicated(pbc3$id), ]
  fit <- survivalSub(first, Surv(years, status3) ~ age + sex, status4 ~ years + age + sex)
  printed <- capture.output(print(fit))
  # status4 = 0 is death (140 patients), the type risk_prob_1 refers to
  expect_true(any(grepl("Cause-1 events \\(status4 = 0\\) : 140", printed)))
  expect_true(any(grepl("Cause-2 events \\(status4 = 1\\) : 29", printed)))
  first$cause <- factor(ifelse(first$status4 == 1, "transplant", "death"))
  fit_factor <- survivalSub(first, Surv(years, status3) ~ age + sex, cause ~ years + age + sex)
  printed <- capture.output(print(fit_factor))
  expect_true(any(grepl("Cause-1 events \\(cause = death\\) : 140", printed)))
  expect_true(any(grepl("Cause-2 events \\(cause = transplant\\) : 29", printed)))
})
