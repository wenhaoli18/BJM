test_that("survivalSubExternal borrows the external tail under a generation term", {
  data(pbc3, envir = environment())
  surv <- pbc3[!duplicated(pbc3$id), ]
  id <- as.numeric(as.character(surv$id))
  study <- surv[id %% 2 == 0, ]
  study$status3[study$years > 4] <- 0
  study$years <- pmin(study$years, 4)
  external <- surv[id %% 2 == 1, ]
  fit <- survivalSubExternal(study, external, Surv(years, status3) ~ age + sex)
  expect_s3_class(fit, "survivalSub.BJM")
  expect_true("generation" %in% names(stats::coef(fit$coxph_fit)))
  expect_equal(fit$external$max_followup_internal, 4)
  expect_gt(fit$external$max_followup_external, 10)
  expect_equal(BJM:::survival_max_followup(fit), max(external$years))
  expect_named(fit$external$ph_test, c("chisq", "df", "p"))
  # prediction data without the generation column are taken as the study
  nd <- study[1:3, ]
  expect_equal(BJM:::survival_lp(fit, nd),
               BJM:::survival_lp(fit, transform(nd, generation = 1)))
  expect_error(survivalSubExternal(transform(study, generation = 1), external,
                                   Surv(years, status3) ~ age), "already a column")
})
