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

test_that("tail_hr fixes the hazard ratio between the cohorts", {
  data(pbc3, envir = environment())
  surv <- pbc3[!duplicated(pbc3$id), ]
  id <- as.numeric(as.character(surv$id))
  study <- surv[id %% 2 == 0, ]
  study$status3[study$years > 4] <- 0
  study$years <- pmin(study$years, 4)
  external <- surv[id %% 2 == 1, ]
  free <- survivalSubExternal(study, external, Surv(years, status3) ~ age + sex)
  expect_false(free$external$tail_hr_fixed)
  expect_equal(free$external$hazard_ratio, free$external$hazard_ratio_estimated)

  # tail_hr = 1: the generations pooled, as one Cox model without the generation term
  one <- survivalSubExternal(study, external, Surv(years, status3) ~ age + sex, tail_hr = 1,
                             splice = FALSE)
  pooled <- survival::coxph(Surv(years, status3) ~ age + sex,
                            rbind(study[c("years", "status3", "age", "sex")],
                                  external[c("years", "status3", "age", "sex")]), ties = "breslow")
  expect_true(one$external$tail_hr_fixed)
  expect_equal(one$external$hazard_ratio, 1)
  expect_equal(one$external$hazard_ratio_estimated, free$external$hazard_ratio_estimated)
  expect_equal(stats::coef(one$coxph_fit)[c("age", "sex")], stats::coef(pooled), tolerance = 1e-6)
  expect_equal(BJM:::survival_cum_basehaz(one)$hazard,
               survival::basehaz(pooled, centered = FALSE)$hazard, tolerance = 1e-6)

  # another value: the study's linear predictor is shifted by log(tail_hr)
  half <- survivalSubExternal(study, external, Surv(years, status3) ~ age + sex, tail_hr = 0.5)
  expect_equal(unname(exp(stats::coef(half$coxph_fit)[["generation"]])), 0.5)
  nd <- study[1:3, ]
  b <- stats::coef(half$coxph_fit)
  expect_equal(unname(BJM:::survival_lp(half, nd)),
               b[["age"]] * nd$age + b[["sex"]] * nd$sex + log(0.5))

  expect_error(survivalSubExternal(study, external, Surv(years, status3) ~ age, tail_hr = -1),
               "single positive number")
  expect_error(survivalSubExternal(study, external, Surv(years, status3) ~ age, tail_hr = c(1, 2)),
               "single positive number")
})

test_that("splice uses the study's own baseline within its follow-up", {
  data(pbc3, envir = environment())
  surv <- pbc3[!duplicated(pbc3$id), ]
  id <- as.numeric(as.character(surv$id))
  study <- surv[id %% 2 == 0, ]
  study$status3[study$years > 4] <- 0
  study$years <- pmin(study$years, 4)
  external <- surv[id %% 2 == 1, ]

  # no covariates: the study's cumulative hazard within follow-up is its Nelson-Aalen estimate
  fit <- survivalSubExternal(study, external, Surv(years, status3) ~ 1)
  expect_true(fit$external$splice)
  bh <- BJM:::survival_cum_basehaz(fit)
  lp <- unique(BJM:::survival_lp(fit, study))
  na <- survival::survfit(Surv(years, status3) ~ 1, study, stype = 2, ctype = 1)
  within <- bh$time <= 4
  expect_equal(exp(lp) * bh$hazard[within],
               stats::approx(na$time, na$cumhaz, xout = bh$time[within], method = "constant", rule = 2)$y)
  # beyond follow-up: the joint fit's increments, continuous at tau
  joint <- survival::basehaz(fit$coxph_fit, centered = FALSE)
  h_tau <- max(joint$hazard[joint$time <= 4])
  beyond <- bh[!within, ]
  expect_equal(beyond$hazard - max(bh$hazard[within]),
               joint$hazard[match(beyond$time, joint$time)] - h_tau)
  expect_equal(max(beyond$time), max(external$years))

  # splice = FALSE keeps the joint fit's baseline throughout
  stacked <- survivalSubExternal(study, external, Surv(years, status3) ~ age + sex, splice = FALSE)
  expect_null(stacked$external$cum_basehaz)
  expect_equal(BJM:::survival_cum_basehaz(stacked), survival::basehaz(stacked$coxph_fit, centered = FALSE))
  # with covariates and a fixed hazard ratio
  sp <- survivalSubExternal(study, external, Surv(years, status3) ~ age + sex, tail_hr = 1)
  expect_true(all(diff(BJM:::survival_cum_basehaz(sp)$hazard) >= 0))
  expect_error(survivalSubExternal(study, external, Surv(years, status3) ~ age, splice = NA), "TRUE or FALSE")
  expect_error(survivalSubExternal(study, external, Surv(years, status3) ~ age + strata(sex)), "stratified")
})
