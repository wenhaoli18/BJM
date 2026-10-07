# Right-censored reference outputs guarding the interval-censoring refactor.
#
# The interval-censoring extension routes every read of the survival
# sub-model (linear predictor, baseline cumulative hazard, outcome variable
# names) through accessor functions instead of touching `coxph_fit`
# directly. For right-censored data those accessors must reproduce the
# pre-refactor numbers exactly, so this computes a broad set of outputs --
# with and without competing risks, a stratified Cox model, several
# patients at once -- and test-rightcensor-reference.R compares them with
# testdata/baseline_rightcensor.rds, generated from master before the
# refactor by running, from the package root:
#
#   devtools::load_all()
#   source("tests/testthat/helper-rightcensor-reference.R")
#   saveRDS(right_censor_reference_outputs(),
#           "tests/testthat/testdata/baseline_rightcensor.rds")
right_censor_reference_outputs <- function() {
  data(pbc3, envir = environment())
  surv_data <- pbc3[!duplicated(pbc3$id), ]

  fits <- list(
    noCR = survivalSub(surv_data, Surv(years, status3) ~ age + sex, NULL),
    CR = survivalSub(surv_data, Surv(years, status3) ~ age + sex, status4 ~ years + age + sex),
    strata = survivalSub(surv_data, Surv(years, status3) ~ age + strata(sex), NULL)
  )

  long_sub_fixed <- list(
    "long1" = serBilir ~ year + age + sex + (years) + (years) * year,
    "long2" = albumin ~ year + age + sex + (years) + (years) * year
  )
  long_sub_random <- list("long1" = ~ year | id, "long2" = ~ year | id)
  data_fit_all <- list(pbc3[pbc3$status3 == 1, ], pbc3[pbc3$status3 == 1, ])
  long_fit_all <- longitudinalSub(data_fit_all, long_sub_fixed, long_sub_random)

  survival_variable_all <- list("Tyears1", "Tyears2", "Tyears3", "Tyears4")
  survival_trans_function <- list(
    fun1 = function(x) abs(x - 1),
    fun2 = function(x) abs(x - 3),
    fun3 = function(x) abs(x - 5),
    fun4 = function(x) abs(x - 7)
  )

  # several patients at risk at year 5, history up to year 5
  ids <- head(surv_data$id[surv_data$years > 6], 5)
  history <- pbc3[pbc3$id %in% ids & pbc3$year <= 5, ]
  data_predict_all <- list(history, history)
  one_patient <- lapply(data_predict_all, function(d) d[d$id == ids[1], ])

  out <- list()
  for (nm in names(fits)) {
    fit <- fits[[nm]]
    at_risk <- subset_at_risk(data_predict_all, "years", 5)
    upper <- integration_upper_bound(at_risk, long_fit_all, fit, 5, min_upper = 6)
    grid <- prepare_infinity_grid(at_risk, long_fit_all, fit, 5, upper, 15)
    out[[nm]] <- list(
      upper_bound = upper,
      grid_edges = grid$predict.time.infinity.1,
      S_T_infinity = grid$S_T_all_infinity,
      marginalT = marginalT(at_risk, long_fit_all, fit, l_i = seq(5, 6, by = 0.1), upper),
      risk = unclass(predictRisk(data_predict_all, long_fit_all, fit,
                                 prediction_time = 5, horizon = 1, time_variable = "year",
                                 survival_variable_all, survival_trans_function,
                                 bandcount1 = 10, bandcount2 = 20)),
      bio = dynamicPredictionBio(1, one_patient, long_fit_all, fit,
                                 prediction_time = 5, horizon = 1, time_variable = "year",
                                 survival_variable_all, survival_trans_function,
                                 bandcount2 = 20, bandcount3 = 30),
      traj = simulateTrajectory(one_patient, long_fit_all, fit,
                                prediction_time = 5, times = c(5.5, 6, 7), time_variable = "year",
                                survival_variable_all, survival_trans_function,
                                n_sim = 20, bandcount2 = 30, seed = 1),
      printed = utils::capture.output(print(fit), summary(fit))
    )
  }
  out
}
