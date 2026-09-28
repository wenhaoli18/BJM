#' Predict multiple future biomarker values from fitted sub-models
#'
#' @description
#' \strong{Internal multi-biomarker engine} behind
#' \code{\link{predictLongitudinal}} -- call \code{predictLongitudinal()}
#' directly instead (it dispatches here automatically when \code{bio_i}
#' names more than one biomarker, or is left \code{NULL}).
#'
#' Batch counterpart to \code{\link{dynamicPredictionBio}}:
#' predicts every requested biomarker in \code{bio_i} from the same fitted
#' longitudinal/survival sub-models and the same \code{data_predict_all}, but
#' computes the bio_i-\strong{independent} part of the pipeline (restricting
#' to at-risk patients, the survival-side integration grid, and the
#' denominator conditional density -- see \code{compute_bio_shared_step()})
#' only \strong{once} and reuses it across every biomarker, instead of
#' recomputing it once per biomarker the way calling
#' \code{dynamicPredictionBio()} in a loop would. This matters most under
#' the Gaussian-copula path (see \code{\link{longitudinalSubCopula}}), where
#' that denominator involves \code{mvtnorm::pmvnorm()} Monte-Carlo
#' evaluations that are otherwise the dominant cost of a multi-biomarker
#' prediction run.
#'
#' \code{bandcount2} (controlling the shared survival-integration grid) and
#' \code{bandcount3} (controlling each biomarker's own candidate-value grid)
#' are auto-tuned separately when left at their default \code{"auto"}: if
#' \code{bandcount2 = "auto"}, it is resolved \strong{once}, using a single
#' representative biomarker (the first one in \code{bio_i}), via the same
#' doubling-until-stable check \code{\link{dynamicPredictionBio}} uses; the
#' shared step is then built once at that resolved value. If
#' \code{bandcount3 = "auto"}, it is then resolved independently for
#' \strong{every} biomarker in \code{bio_i} (their candidate-value grids
#' need not converge at the same resolution), reusing the once-computed
#' shared step for every doubling round rather than rebuilding it -- except
#' for any \strong{ordinal} biomarker, for which \code{bandcount3} tuning is
#' always skipped (its candidate grid is fixed at its category count; see
#' \code{\link{dynamicPredictionBio}}) and \code{"bandcount3"} is recorded
#' as \code{NA} for that biomarker. See Details in
#' \code{\link{dynamicPredictionBio}} for the general auto-bandcount
#' rationale.
#'
#' @param data_predict_all See \code{\link{dynamicPredictionBio}}.
#' @param long_fit_all See \code{\link{dynamicPredictionBio}}.
#' @param survival_fit_all See \code{\link{dynamicPredictionBio}}.
#' @param prediction_time See \code{\link{dynamicPredictionBio}}.
#' @param horizon See \code{\link{dynamicPredictionBio}}.
#' @param time_variable See \code{\link{dynamicPredictionBio}}.
#' @param survival_variable_all See \code{\link{dynamicPredictionBio}}.
#' @param survival_trans_function See \code{\link{dynamicPredictionBio}}.
#' @param bandcount2 See \code{\link{dynamicPredictionBio}}.
#' @param bandcount3 See \code{\link{dynamicPredictionBio}}.
#'
#' @param bio_i Integer vector of biomarkers to predict. May include
#' continuous and/or ordinal biomarkers (see
#' \code{\link{dynamicPredictionBio}}). Defaults to \code{NULL}, meaning
#' every biomarker in \code{long_fit_all}.
#' @return A named list of \code{"dynamicPredictionBio.BJM"} objects (see
#' \code{\link{dynamicPredictionBio}}), one per requested biomarker, named
#' by that biomarker's response-variable name; with attributes
#' \code{"bandcount2"} (the single resolved/used \code{bandcount2}) and
#' \code{"bandcount3"} (a named numeric vector of the resolved/used
#' \code{bandcount3} for each biomarker, or \code{NA} for an ordinal
#' biomarker). Classed \code{"dynamicPredictionBioAll.BJM"}.
#'
#' @examples
#'
#' \donttest{
#' data(pbc3)
#'
#' data_survival_fitting =  pbc3[!duplicated(pbc3$id), ]
#'
#' form_marginal_surv = Surv(years, status3) ~ age + sex
#' form_conditional_cr = NULL
#'
#' survival_fit_all = survivalSub(data_survival_fitting, form_marginal_surv,
#'                                form_conditional_cr)
#'
#' long_sub_fixed = list(
#'   "long1" = serBilir ~ year + age + sex +  (years) + (years) * year,
#'   "long2" = prothrombin ~ year + age + sex + (years) + (years) * year,
#'   "long3" = albumin ~ year + age + age * year + sex + (years) + (years) * year,
#'   "long4" = alkaline ~ year + age + sex + (years) + (years) * year,
#'   "long5" = SGOT ~ year + age + sex + (years) + (years) * year,
#'   "long6" = platelets ~ year + age + sex + (years)  + (years) * year)
#'
#' long_sub_random =list(
#'   "long1" =  ~ year| id,
#'   "long2" =  ~ year| id,
#'   "long3" =  ~ year| id,
#'   "long4" =  ~ year| id,
#'   "long5" =  ~ year| id,
#'   "long6" =  ~ year| id)
#'
#' survival_variable_all = list(
#'   "Tyears1",  "Tyears2", "Tyears3", "Tyears4"
#' )
#'
#' survival_trans_function = list(
#'   fun1 = function(x){abs(x - 1)},
#'   fun2 = function(x){abs(x - 3)},
#'   fun3 = function(x){abs(x - 5)},
#'   fun4 = function(x){abs(x - 7)}
#' )
#'
#' # Complete case analysis
#' data_fit_all = list()
#' for(i in seq_len(length(long_sub_fixed))){
#'   data_fit_all[[i]] = pbc3[pbc3$status3 == 1, ]
#' }
#'
#' # fitting longitudinal submodel
#' long_fit_all = longitudinalSub(data_fit_all, long_sub_fixed, long_sub_random)
#'
#' i_PID = 2
#' data.raw.predict.1 = pbc3[pbc3$id == i_PID, ]
#'
#' data_predict_all = list()
#' for(i in seq_len(length(long_sub_fixed))){
#'   data_predict_all[[i]] = data.raw.predict.1[data.raw.predict.1$year <= 3,]
#' }
#'
#' Y_predict_all = dynamicPredictionBioAll(bio_i = NULL, data_predict_all, long_fit_all,
#'                                         survival_fit_all, prediction_time = 3,
#'                                         horizon = 3, time_variable = "year",
#'                                         survival_variable_all, survival_trans_function,
#'                                         bandcount2 = 40, bandcount3 = 400)
#'
#' }
#'
#' @keywords internal
dynamicPredictionBioAll <- function(bio_i = NULL, data_predict_all, long_fit_all, survival_fit_all,
                                     prediction_time, horizon, time_variable,
                                     survival_variable_all, survival_trans_function,
                                     bandcount2 = "auto", bandcount3 = "auto") {

  assert_class(long_fit_all, "longitudinalSub.BJM", "long_fit_all", "longitudinalSub")
  assert_class(survival_fit_all, "survivalSub.BJM", "survival_fit_all", "survivalSub")

  n_long <- length(long_fit_all$lfit)
  if (is.null(bio_i)) {
    # Both continuous and ordinal biomarkers are predictable (see
    # dynamicPredictionBio()), so the default is simply every biomarker.
    bio_i <- seq_len(n_long)
  } else {
    for (b in bio_i) {
      assert_index(b, n_long, "bio_i", "longitudinal outcomes in long_fit_all")
    }
  }

  assert_data_list(data_predict_all, "data_predict_all", n_long, allow_bare_df = TRUE)
  if (!is.list(data_predict_all) || is.data.frame(data_predict_all)) {
    data_predict_all <- rep(list(data_predict_all), each = n_long)
  }
  assert_scalar_numeric(prediction_time, "prediction_time")
  assert_scalar_numeric(horizon, "horizon")
  assert_string(time_variable, "time_variable")
  assert_bandcount(bandcount2, "bandcount2")
  assert_bandcount(bandcount3, "bandcount3")
  assert_survival_trans(survival_variable_all, survival_trans_function, probe_value = prediction_time)

  coxph_fit = survival_fit_all$coxph_fit
  survival_variable = as.character(formula(coxph_fit)[[2]])[2]
  for (i in seq_along(data_predict_all)) {
    assert_vars_in_data(time_variable, data_predict_all[[i]],
                         "time_variable", sprintf("data_predict_all[[%d]]", i))
    assert_vars_in_data(survival_variable, data_predict_all[[i]],
                         "the survival-time variable used to fit survival_fit_all",
                         sprintf("data_predict_all[[%d]]", i))
  }

  # bandcount2 is shared across every biomarker: if "auto", resolve it once
  # using only the first requested biomarker -- exactly what a single
  # dynamicPredictionBio() call would do -- via the existing, already-tested
  # auto_tune_bandcount() machinery, then reuse the resolved value to build
  # the shared step once below, for every biomarker.
  resolved_bandcount2 <- bandcount2
  if (identical(bandcount2, "auto")) {
    probe_bandcount3 <- if (identical(bandcount3, "auto")) unname(bandcount_auto_start["bandcount3"]) else bandcount3
    probe_args <- list(bio_i = bio_i[1], data_predict_all = data_predict_all, long_fit_all = long_fit_all,
                        survival_fit_all = survival_fit_all, prediction_time = prediction_time,
                        horizon = horizon, time_variable = time_variable,
                        survival_variable_all = survival_variable_all,
                        survival_trans_function = survival_trans_function,
                        bandcount2 = "auto", bandcount3 = probe_bandcount3)
    resolved_bandcount2 <- auto_tune_bandcount(dynamicPredictionBio, probe_args, "bandcount2")$bandcount$bandcount2
  }

  shared <- compute_bio_shared_step(data_predict_all, long_fit_all, survival_fit_all,
                                     prediction_time, time_variable, survival_variable_all,
                                     survival_trans_function, resolved_bandcount2)

  results <- list()
  bandcount3_used <- c()
  for (b in bio_i) {
    b_is_ordinal <- !is.null(long_fit_all$biomarker_type) && long_fit_all$biomarker_type[b] == "ordinal"
    if (b_is_ordinal) {
      # An ordinal marker's candidate grid is fixed at its category count,
      # not controlled by bandcount3 at all (see compute_bio_marker_step()),
      # so bandcount3 auto-tuning is skipped for it regardless of the
      # bandcount3 argument's value.
      marker_result <- compute_bio_marker_step(shared, b, long_fit_all, survival_fit_all,
                                                 prediction_time, horizon, time_variable,
                                                 survival_variable_all, survival_trans_function,
                                                 bandcount3)
      bandcount3_used[as.character(b)] <- NA_real_
    } else if (identical(bandcount3, "auto")) {
      tuned <- auto_tune_marker_bandcount3(shared, b, long_fit_all, survival_fit_all,
                                            prediction_time, horizon, time_variable,
                                            survival_variable_all, survival_trans_function)
      marker_result <- tuned$result
      bandcount3_used[as.character(b)] <- tuned$bandcount3
    } else {
      marker_result <- compute_bio_marker_step(shared, b, long_fit_all, survival_fit_all,
                                                 prediction_time, horizon, time_variable,
                                                 survival_variable_all, survival_trans_function,
                                                 bandcount3)
      bandcount3_used[as.character(b)] <- bandcount3
    }
    out <- list(Y_predict = marker_result$Y_predict, Y_density = marker_result$Y_density,
                Y_all = marker_result$Y_all)
    class(out) <- "dynamicPredictionBio.BJM"
    bio_name <- as.character(formula(long_fit_all$long_sub_fixed[[b]])[[2]])
    results[[bio_name]] <- out
  }

  attr(results, "bandcount2") <- resolved_bandcount2
  attr(results, "bandcount3") <- bandcount3_used
  class(results) <- "dynamicPredictionBioAll.BJM"
  results
}
