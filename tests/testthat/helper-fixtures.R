# Shared fixture for dynamicPrediction / dynamicPredictionBio tests.
# testthat auto-sources helper-*.R files before running any test-*.R file,
# making this available across all test files (unlike a top-level function
# defined inside a test-*.R file, which is not reliably shared).
setup_dp_fixture <- function() {
  data(pbc3, envir = environment())

  data_survival_fitting <- pbc3[!duplicated(pbc3$id), ]
  survival_fit_all <- survivalSub(data_survival_fitting,
                                   Surv(years, status3) ~ age + sex,
                                   status4 ~ years + age + sex)

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

  data.raw.predict <- pbc3[pbc3$id == 2, ]
  data_predict_all <- list(data.raw.predict, data.raw.predict)

  list(survival_fit_all = survival_fit_all, long_fit_all = long_fit_all,
       survival_variable_all = survival_variable_all,
       survival_trans_function = survival_trans_function,
       data_predict_all = data_predict_all)
}
