# Shared fixture for predictRisk / dynamicPredictionBio tests.
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

# Shared fixture for imputeLongitudinal() tests: a small, artificially
# MCAR-masked slice of pbc3 (a handful of subjects, so the MIWAE trains
# quickly under test). serBilir/albumin cells are set to NA at random,
# independent of every other variable, which satisfies MAR.
setup_impute_fixture <- function(seed = 1) {
  data(pbc3, envir = environment())

  data_fit_all <- pbc3[pbc3$status3 == 1, ]
  keep_ids <- unique(data_fit_all$id)[1:60]
  data_fit_all <- data_fit_all[data_fit_all$id %in% keep_ids, ]

  long_sub_fixed <- list(
    "long1" = serBilir ~ year + age + sex + years,
    "long2" = albumin ~ year + age + sex + years
  )
  long_sub_random <- list("long1" = ~ year | id, "long2" = ~ year | id)

  set.seed(seed)
  masked <- data_fit_all
  n <- nrow(masked)
  masked$serBilir[sample.int(n, size = max(1, floor(0.12 * n)))] <- NA
  masked$albumin[sample.int(n, size = max(1, floor(0.12 * n)))] <- NA

  list(data_fit_all = masked, long_sub_fixed = long_sub_fixed,
       long_sub_random = long_sub_random, time_variable = "year")
}

# longitudinalSub() does not return the number of subjects it actually
# retained (its `yi` field is an internal intermediate used only to compute
# `Sigma_fit`, never part of the returned object) -- so subject retention
# has to be recomputed directly, by replicating the same complete-case +
# cross-biomarker-intersection rule longitudinalSub() itself documents (see
# NEWS.md, "BJM 0.3.0"): a subject is retained only if it has at least one
# non-missing observation of every biomarker's own fixed-effect variables.
count_retained_subjects <- function(data_fit_all, long_sub_fixed, id = "id") {
  if (is.data.frame(data_fit_all)) {
    data_fit_all <- rep(list(data_fit_all), length(long_sub_fixed))
  }
  ids_by_marker <- lapply(seq_along(long_sub_fixed), function(m) {
    d <- data_fit_all[[m]]
    complete <- d[rowSums(is.na(d[all.vars(long_sub_fixed[[m]])])) == 0, ]
    unique(complete[[id]])
  })
  Reduce(intersect, ids_by_marker)
}
